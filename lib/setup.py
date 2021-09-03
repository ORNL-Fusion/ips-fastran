from setuptools import setup, find_packages

setup(
   name='Namelist',
   version='0.1',
   packages=[''],
   package_dir={'':'namelist'}
)

setup(
   name="pyps",
   version='0.1',
   packages=[''],
   package_dir={'':'pyps'},
   package_data={'': ['_pyplasmastate.so']},
) 

setup(
   name="omfit",
   version='0.1',
   packages=["omfit", "omfit.classes"],
) 
