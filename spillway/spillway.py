import os
import numpy as np
from pandas import read_csv
from matplotlib import pyplot as plt
from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.interpolate import interp1d
from matplotlib.patches import Polygon
from warnings import warn

class Spillway(object):

	"""
	A prismatic spillway of some shape that connects two reservoirs of water with some volume and head. 
	At the moment, it is rectangular, so that the only parameter that matters is the width of the spillway. 
	All elevations are made in reference to a common datum, which is the "bedrock," a flat surface 
	underneath the bottom of the spillway. There is also a rim, which is the maximum elevation water can rise to. 
	The spillway has an unerodible bed with some slope. We track the water elevation at the entrance,
	the exit, and the profile along the spillway.
	"""

	def __init__(self):
		self.S = None		# the spillway slope
		self.b = None		# the spillway width
		self.Zc = None		# section factor for critical flow
		self.K = None		# conveyance factor for uniform flow
		self.L = None		# the spillway length
		self.H = None		# the height of the spillway rim
		self.Q = None		# the discharge through the spillway
		self.x = None		# the x coordinate along the spillway, 0 is the inlet, increasing toward the exit. Negative is the upstream reservoir.
		self.eta = None		# the surface of the bed (for both reservoirs and spillway)
		self.ya = None		# the upstream basin reservoir level.
		self.yb = None		# the downstream basin reservoir level.
		self.za = None		# the upstream basin floor level.
		self.zb = None		# the downstream basin floor level.
		self.t = None 		# the time vector through history (s)
		self.ya_t = None 	# the history of the upstream basin reservoir level.
		self.yb_t = None 	# the history of the downstream basin reservoir level.
		self.make_submerged_weir_correction()

	def set_defaults(self):
		"""
		Set a few default constants and dimensions, to ease development. 
		"""
		self.L = 2.44					# meters
		self.La = 4.5					# meters the length of the upstream basin
		self.Lb = 11 - self.La - self.L	# meters the flume is 11 m long,
		self.za = 0.1					# meters, this value provides for just-supercritical flow if yn = 2.5 cm
		self.zb = 0						# meters
		self.delH = self.za - self.zb	# meters
		self.H = 0.381					# meters
		self.n = 0.015					# metric units for n
		self.Cd = 0.35					# from engineering literature, units uncertain
		self.B = 0.1524					# meters
		self.g = 9.8 					# meters / s^2
		self.set_b(0.0324)				# meters
		self.set_S()				 	# dimensionless, based on delH and L
		self.set_x(x = np.arange(-self.La, self.Lb + self.L, 0.01)) # meters 
		self.set_eta(self.za)			# the profile in m along x

	def set_b(self, b):
		"""
		This is to set the width of the spillway
		"""
		self.b = b

	def set_y_values(self, ya, yb):
		"""
		This is to set the upstream and downstream basin levels.
		"""
		if ya + self.za > self.H:
			raise ValueError("upstream water elevation is above the lip of the flume")
		elif yb + self.zb > self.H:
			raise ValueError("downstream water elevation is above the lip of the flume")
		else:
			self.ya = ya
			self.yb = yb

	def set_za(self, za):
		self.za = za
		self.set_S()
		self.set_eta(self.za)			# the profile in m along x

	def set_S(self): 
		"""
		importantly assumes that zb is zero and that the difference in the elevation of the floor of each reservoir is relative to downstream.
		"""
		self.S = self.za / self.L

	def set_x(self, x = None, low = None, high = None, step = None):
		"""
		Initialize the x domain

		Either supply an array or specify the beginning and ending value as well as the step.
		"""
		if x is not None:
			self.x = x
		elif (x is None) and ((low or high or step) is not None): 
			self.x = np.arange(low, high, step)
		else:
			raise ValueError("x must be either a 1D array or values to make such an array.")

	def set_eta(self, za, delH = None, L = None, S = None):
		"""
		Initialize the bed surface.

		In this model, the bed in the reservoir is flush with the spillway at the entrance and exit.
		"""
		if self.H is None:
			pass
		elif za >= self.H:
			raise ValueError("upstream elevation is above or at the maximum value ")
		if self.x is None:
			raise ValueError("specify an x domain first.")
		if self.L is None:
			self.L = L
		if self.S is None:
			self.S = S
		self.za = za
		self.eta = np.full(self.x.shape, 0, dtype = float)
		self.eta[self.x < 0] = 0
		inSpillway = (self.x <= self.L) & (self.x >= 0)
		self.eta[inSpillway] = self.za - self.S * self.x[inSpillway]
		lowest_eta = np.min(self.eta[inSpillway])
		if lowest_eta < 0:
			raise ValueError("choice of slope takes downstream reservoir below zero-datum")
		else:
			self.eta[self.x > self.L] = lowest_eta

	def find_normal_discharge(self, y):
		"""
		Find the uniform flow discharge of a rectangular channel with given slope, width, and roughness provided a flow depth.
		Can be sub- or super-critical.
		"""
		aa = (self.b * y) / (2 * y + self.b)
		bb = aa ** (2/3)
		Q = (np.sqrt(self.S) * self.b * y * bb) / self.n
		return Q

	def find_critical_discharge(self, y):
		"""
		Find the critical flow discharge of a rectangular channel with given width provided a flow depth.
		"""
		Q = y**(3/2) * self.b * np.sqrt(self.g)
		return Q 

	def gradually_varied_flow_diffEQ(self, x = 0, y = 0, Q = 0):
		"""
		The equation for gradually-varied flow along a rectangular channel of known width, b.
		Expresses the gradient of the water surface in terms of the current elevation and the
		critical and normal flow depths for that channel.
		returns the spatial derivative of the water depth.
		requires x variable for use with standard python ODE solvers
		"""
		yc = self.find_critical_flow_depth(Q)
		yn = self.find_normal_flow_depth(Q)

		AA = (yn / y)**(10/3)
		BB = ((2*y + self.b)/(2*yn + self.b))**(4/3)
		CC = (yc / y)**3
		dydx = self.S * (1 - AA * BB)/(1 - CC)
		return dydx

	def find_froude(self, Q, y):
		"""
		Find the Froude number of a given discharge and flow depth,
		that is flowing through a prismatic rectangular channel of known width and slope.
		"""
		uu = Q / (self.b * y)
		fr = uu / np.sqrt(self.g * y)
		return fr

	def find_normal_flow_depth(self, Q):
		"""
		For a rectangular channel, the uniform-flow formula cannot be re-expressed as an explicit formula for y
		so here, we solve for it using a root-finding algorithm.
		At the moment, we use a bracket to find the root that only goes from a flow depth of zero to the 
		maximum depth of the flume. This is a sensical bound for this sepecific setup. 
		An error here likely indicates that H is higher than the lip of the flume.
		"""
		objective = lambda y, Q: self.find_normal_discharge(y) - Q
		sol = root_scalar(objective, args = (Q), bracket = [0, self.H])
		return sol.root

	def find_critical_flow_depth(self, Q):
		"""
		Using the section factor for critical flow, and rearranging the flow equation to find the critical depth.
		"""
		y = (Q / (self.b * np.sqrt(self.g))) ** (2/3)
		return y

	def find_section_factor(self, y):
		"""
		This assumes that critical flow is established at some point along the spillway. 
		This calculation is agnostic to the point at which this occurs. 
		Thus, the `y` value specified here is the depth at the critical section, either the top of the spillway or the bottom of the spillway. 
		"""
		self.Zc = self.b * y * np.sqrt(y)
		return self.Zc

	def critical_discharge_given_slope(self, Qmax = 1e-2, Qmin = 1e-4):
		"""
		Takes a window of discharge conditions, and knowing the current geometry of the channel, 
		finds the discharge of a critical flow
		"""
		def objective(q):
			yn = self.find_normal_flow_depth(q)
			Fr = self.find_froude(q, yn)
			return Fr - 1
		sol = root_scalar(objective, bracket = [Qmin, Qmax])
		self.Qc = sol.root
		return self.Qc

	def check_upstream_supercritical(self):
		"""
		"""
		# first is to find the critical depth at the outlet if it is just draining
		yc = self.ya * 2/3
		# next, compute that discharge assuming you have critical flow. 
		qc = self.find_critical_discharge(yc)
		# next, find normal depth for that discharge
		yn = self.find_normal_flow_depth(qc)
		crit_test = yn < yc # flow is supercritical if yn < yc
		return crit_test

	def falling_rate_upstream(self):
		"""
		returns the relationship between the level of the basin and the discharge
		"""
		basin_area = self.B * self.La
		dydt = -self.Q / basin_area
		return dydt

	def rising_rate_downstream(self):
		"""
		returns the relationship between the level of the basin and the discharge
		"""
		basin_area = self.B * self.Lb
		dydt = self.Q / basin_area
		return dydt

	def make_subcritical_discharge_interpolator(self, Qvals, eval_points = 20):
		"""
		A function that takes the geometry of the spillway system, and computes
		an interpolater for the discharge of the spillway, given the upstream 
		and downstream pool elevations.
		"""
		self.yb_lookup = np.array([], dtype = float)
		self.ya_lookup = np.array([], dtype = float)
		self.Q_lookup = np.array([], dtype = float)
		for q in Qvals:
			yb, ya, QQ = self.compute_subcritical_discharge_curve(q, eval_points = eval_points)
			self.yb_lookup = np.concatenate((self.yb_lookup, yb))
			self.ya_lookup = np.concatenate((self.ya_lookup, ya))
			self.Q_lookup = np.concatenate((self.Q_lookup, QQ))
		lookup_obj = bisplrep(self.yb_lookup, self.ya_lookup, self.Q_lookup)
		def lookup_function(yb, ya):
			qval = bisplev(yb, ya, lookup_obj)
			return qval
		self.interpolator = lookup_function
		return lookup_function

	def compute_subcritical_discharge_curve(self, Q, ymin = 1e-3, ymax = None, eval_points = 20):
		"""
		compute a q-constant curve across the space of ya and yb
		Assuming that the flow in the canal is subcritical
		"""
		# we will limit the max downstream depth to 5 cm below the flume lip
		if ymax is None:
			ymax = self.H - 0.05

		yb = np.linspace(ymin, ymax, eval_points)		
		ya = np.full(yb.shape, 0.0)
		yc = self.find_critical_flow_depth(Q)
		yn = self.find_normal_flow_depth(Q)

		if yn <= yc:
			ya = np.array([3/2 * yc for y in yb])
			y_allowed = yb < (ya + self.S * self.L)
			ya = ya[y_allowed]
			yb = yb[y_allowed]
		else:
			for i, y in enumerate(yb):
				if y <= yc:
					ya[i] = self.subcritical_upstream_depth(Q, yc*1.001)
				else:
					ya[i] = self.subcritical_upstream_depth(Q, y)

		return yb, ya, np.full(yb.shape, Q)

	def compute_crit_depth_line(self, Qmax, eval_points = 20):
		"""
		"""
		# we use the critical discharge for the slope of this channel as the starting 
		# place, so that the line pertains only to subcritical flows where the 
		# downstream drawdown goes into freefall.
		Qvals = np.linspace(self.Qc * 1.05, Qmax, eval_points)
		
		y_norm = np.array([self.find_normal_flow_depth(q) for q in Qvals])

		# y_norm = self.find_normal_flow_depth(Qvals)

		yb_crit = self.find_critical_flow_depth(Qvals)
		ya_crit = np.full(yb_crit.shape, 0.0)

		for q, (i, y), yn in zip(Qvals, enumerate(yb_crit), y_norm):
			if y >= yn:
				ya_crit[i] = 3/2 * y
			else:
				ya_crit[i] = self.subcritical_upstream_depth(q, y * 1.001)

		self.ya_crit = ya_crit
		self.yb_crit = yb_crit
		self.Q_crit = Qvals

		return ya_crit, yb_crit, Qvals

	def compute_subcritical_flow_profile(self, Q, yb):
		"""
		Use the gradually varied flow equation to find the flow profile 
		discharge and downstream depth are known.
		"""
		ys = solve_ivp(self.gradually_varied_flow_diffEQ, [self.L, 0], [yb], args = [Q], max_step = self.L/100)
		return np.flip(ys.t), np.flip(ys.y[0]) # return the list, but flip since we computed backwards

	def compute_supercritical_flow_profile(self, Q):
		"""
		Use the gradually varied flow equation to find the flow profile for an S2 curve
		"""
		ys = solve_ivp(self.gradually_varied_flow_diffEQ, [0, self.L], [yc * 0.99], args = [Q], max_step = self.L/100)
		return ys.t, ys.y[0]

	def subcritical_upstream_depth(self, Q, yb):
		"""
		Use the gradually varied flow equation to find the upstream depth when the 
		discharge and downstream depth are known.
		"""
		ys = solve_ivp(self.gradually_varied_flow_diffEQ, [self.L, 0], [yb], args = [Q], max_step = self.L/100)
		return ys.y[0][-1] # return the last element in the list, because we computed backwards

	def find_subcritical_discharge(self, Qmin, Qmax):
		"""
		Takes a window of discharge conditions. By knowing the current geometry of the channel and
		the elevation of downstream reservoir, finds the discharge that matches the upstream reservoir
		"""
		def objective(q):
			ya = self.subcritical_upstream_depth(q, self.yb)
			return self.ya - ya
		sol = root_scalar(objective, bracket = [Qmin, Qmax])
		self.Q = sol.root
		return self.Q

	def make_submerged_weir_correction(self):
		"""
		"""
		this_file_location = os.path.dirname(os.path.realpath(__file__))
		ES207_data_file = os.path.join(this_file_location, 'data', 'ES-207.csv')
		self.ES207_data = read_csv(ES207_data_file)
		the_cols = self.ES207_data.columns
		exx = self.ES207_data.loc[:, the_cols[0]].to_numpy()
		why = self.ES207_data.loc[:, the_cols[1]].to_numpy()
		self.submerge_interp = interp1d(exx, why, kind = 'cubic')

	def compute_discharge_arbitrary(self):
		"""
		"""
		# Should I make a helper function for this?
		if self.yb + self.zb >= self.ya + self.za:
			self.Q = 0
		elif self.check_upstream_supercritical():
			
			yc = 2/3 * self.ya

			# here we select a condition that in order for flow to be supercritical and run freely,
			# critical flow has to be achieved. At the moment, this criterion simply specifies that 
			# there has to be a downstream gradient

			############################################################################
			# BUT: here is the place to experiment with this criterion.
			############################################################################
			
			if self.za > self.yb + self.zb:
				self.Q = self.find_critical_discharge(yc)
			else:
				Qfree = self.find_critical_discharge(yc)
				# print((self.ya + self.za, self.yb, (self.yb - self.za) / self.ya))
				self.Q = Qfree * self.submerge_interp((self.yb - self.za) / self.ya)
		else:
			qc = self.find_critical_discharge(self.yb)
			qn = self.find_normal_discharge(self.ya)
			yc = self.find_critical_flow_depth(qn)

			if self.yb > self.ya:
				self.Q = self.find_subcritical_discharge(0, 0.999 * qc)
			elif (self.yb < self.ya) and (self.yb >= yc):
				self.Q = self.find_subcritical_discharge(0, 0.999 * qc)
			elif (self.yb < yc):
				self.Q = qn
			else:
				print((self.ya, self.yb, self.Q))
				raise ValueError('something has produced an nan')
		return self.Q

	def spillway_ode(self, t, y):
		"""
		"""
		
		# take ODE solver vector and update object
		self.yb = y[0]
		self.ya = y[1]
		
		# compute the discharge
		Q = self.compute_discharge_arbitrary()
		# self.find_broadcrested_weir_discharge() ## this is a simpler alternative. not as physically-based

		# compute rise and fall of downstream pools
		dybdt = self.rising_rate_downstream()
		dyadt = self.falling_rate_upstream()
		return [dybdt, dyadt]
	
	def evolve_system(self, t = 30, time_step = np.inf):
		"""
		Use the diffeq system to construct a history of the rising and falling waters for 30 s (default)
		"""
		
		ys = solve_ivp(self.spillway_ode, [0, t], [self.yb, self.ya], max_step = time_step)
		
		self.t = ys.t
		self.y_t = ys.y
		self.ya_t = ys.y[1]
		self.yb_t = ys.y[0]

		return ys.t, ys.y

########## Plotting Functions ##########

	def plot_one_Q_constant_curve(self, yb_vals, ya_vals, Q = None):
		"""
		"""
		ax = self.qConstAx

		ax.plot(yb_vals, ya_vals, '-', linewidth = 3, color = '#1f77b4')

		self.qConstAx = ax

		return ax

	def plot_computed_Q_constant_curves(self):
		"""
		need to compute the interpolator before this will work
		"""
		for q in np.unique(self.Q_lookup):
			thisQ = self.Q_lookup == q
			self.plot_one_Q_constant_curve(self.yb_lookup[thisQ], self.ya_lookup[thisQ])


	def make_base_Q_curve_plot(self, add_annotation = False):
		"""
		make a base plate on which to plot q constant curves
		"""
		ydummy = np.arange(0, self.H, 0.01)

		fig, ax = plt.subplots(figsize = (5,5))
		ax.set_xlim(0, self.H)
		ax.set_ylim(0, self.H)
		ax.plot(ydummy, ydummy, '0.5')
		ax.plot(ydummy + self.S * self.L, ydummy, '0.5')
		ax.plot(2/3*ydummy, ydummy, '0.8')

		if add_annotation:
			ax.annotate(
				"Normal Depth Line", 
				(self.H / 2 * 1.05, self.H / 2 * 0.95), 
				rotation = 45, color = '0.5', ha = 'center', va = 'center'
			)
			head_diff = self.S * self.L
			limit_anno_center = ((self.H + head_diff) / 2 * 1.05, (self.H - head_diff) / 2 * 0.95)
			ax.annotate(
				"Limit Depth Line", 
				limit_anno_center, 
				rotation = 45, color = '0.5', ha = 'center', va = 'center'
			)

			crit_anno_center = ((self.H*1/3) * 0.95, (self.H) / 2 * 1.05)
			ax.annotate(
				"Approximate Critical Depth Line", 
				crit_anno_center, 
				rotation = 56.3, color = '0.8', ha = 'center', va = 'center'
			)

		ax.set_xlabel('Downstream depth (m)')
		ax.set_ylabel('Upstream depth (m)')

		self.qConstFig = fig
		self.qConstAx = ax

		return fig, ax

	def plot_crit_depth_line(self, add_annotation = False):
		"""
		"""
		if not hasattr(self, 'ya_crit'):
			raise RuntimeError("critical-depth line has not been computed.\
			run `compute_crit_depth_line` with a maximum discharge")
		
		ax = self.qConstAx

		ax.plot(self.yb_crit, self.ya_crit, 'r-.')

		if add_annotation:
			crit_anno_center = ((self.H*1/3) * 0.95, (self.H) / 2 * 1.05)
			ax.annotate(
				"Critical Depth Line", 
				crit_anno_center, 
				rotation = 56.3, color = 'r', ha = 'center', va = 'center'
			)

		self.qConstAx = ax

		return ax

	def plot_flow_profile(self, Q, x, y):
		"""
		"""
		yc = self.find_critical_flow_depth(Q)
		yn = self.find_normal_flow_depth(Q)

		fig, ax = plt.subplots(figsize = (15,3))
		xs = np.linspace(0, self.L, 100)
		bed = np.flip(xs * self.S)

		yprofile = np.interp(xs, x, y)
		ax.plot(xs, bed, 'k')
		ax.plot(xs, bed + yc, 'r--')
		ax.plot(xs, bed + yn, 'k--')
		ax.plot(xs, bed + yprofile)
		ax.plot(xs[-1], bed[-1] + yprofile[-1], 'ko')
		ax.plot(xs[0], bed[0] + yprofile[0], 'ro')

		ax.set_xlabel('Downstream distance (m)')
		ax.set_ylabel('Elevation (m)')

		return fig, ax

	def show_current_condition(self):
		"""
		This function just makes a cartoon plot of the apparatus given the current conditions that are set.
		"""
		fig, ax = plt.subplots(figsize = (15, 4))
		ax.plot(self.x, self.eta, 'k')
		ax.set_ylim([-0.02, 0.4])
		for s in ax.spines:
			ax.spines[s].set_visible(False)
			
		# show the spillway filled in as a triangle-shaped weir. 

		vertsX = [
			0, 0, self.L, self.L + self.Lb,
			self.L + self.Lb, 0

		]
		vertsY = [
			0, self.za, 
			self.zb, self.zb, self.zb, self.zb
		]

		verts = np.stack((vertsX, vertsY))
		verts.transpose().shape
		pp = Polygon(verts.transpose(), facecolor = '0.8')
		ax.add_patch(pp)

		# show the initial water level filling the upper pool. 

		vertsAX = [
			-self.La, 0, 
			0, -self.La
		]
		vertsAY = [
			self.ya + self.za, self.ya + self.za, 
			0, 0
		]

		vertsA = np.stack((vertsAX, vertsAY))
		vertsA.transpose().shape
		ppA = Polygon(vertsA.transpose(), facecolor = '#61b8e0')
		ax.add_patch(ppA)

		# show the initial water level filling the lower pool. 

		if self.yb + self.zb < self.za:
			bintersect = ((self.za) - (self.zb + self.yb)) / self.S
			vertsBX = [
				bintersect, self.L + self.Lb, 
				self.L + self.Lb, self.L, bintersect
			]
			vertsBY = [
				self.yb + self.zb, self.yb + self.zb, 
				self.zb, self.zb, self.yb + self.zb
			]
		elif self.yb + self.zb >= self.za:
			vertsBX = [
				0, 0, self.L + self.Lb, 
				self.L + self.Lb, self.L, 0
			]
			vertsBY = [
				self.za, self.yb + self.zb, self.yb + self.zb, 
				self.zb, self.zb, self.za
			]


		vertsB = np.stack((vertsBX, vertsBY))
		vertsB.transpose().shape
		ppB = Polygon(vertsB.transpose(), facecolor = '#61b8e0')
		ax.add_patch(ppB)

		# now plot the apparatus lines on top.

		ax.plot((0, 0), (self.za, self.za + self.ya + 0.05), 'k-')
		ax.plot((-self.La, self.Lb + self.L), (0,0), 'k-')
		ax.plot((-self.La, -self.La), (0, self.H), 'k-')
		ax.plot((self.Lb + self.L, self.Lb + self.L), (0, self.H), 'k-')

		ax.set_xlabel('Downstream distance (m)')
		ax.set_ylabel('Elevation (m)')

		return fig, ax
