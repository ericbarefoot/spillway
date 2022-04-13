import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spillway",
    version="0.0.9000",
    author="Eric A. Barefoot",
    author_email="eabarefoot@gmail.com",
    description="Discharge and water surface profiles on spillways",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ericbarefoot/spillway",
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Open-Channel Hydraulics",
        "Intended Audience :: Science/Research",
    ],
    keywords='channel spillway open-channel hydraulics',
    include_package_data=True,
    package_data={'': ['data/*.csv']},
    # project_urls={
    #     'Model page': 'https://csdms.colorado.edu/wiki/Model:GRLP',
    # },
)
