from setuptools import setup

setup(
   name='orbital_mechs',
   version='0.1',
   description='Computational Spaceflight Dynamics',
   author='Nicolas Ehrlich',
   author_email='nicolas.ehrlich@yahoo.de',
   packages=['orbital_mechs'],  #same as name
   install_requires=['numpy'], #external packages as dependencies
)
