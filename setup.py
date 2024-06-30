from setuptools import setup

setup(
   name='orbital_mechs',
   version='0.1.1',
   description='Computational Spaceflight Dynamics',
   author='Nicolas Ehrlich',
   author_email='nicolas.ehrlich@yahoo.de',
   packages=['orbital_mechs'],
   package_data={'orbital_mechs': ['kernels/*']},
   install_requires=[
       'numpy', 
       'scipy', 
       'spiceypy', 
       'lamberthub'],
)
