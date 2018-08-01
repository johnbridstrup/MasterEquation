import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='StateSim',
    author='John Bridstrup',
    email='john.bridstrup@gmail.com',
    version='0.1dev',
    packages=setuptools.find_packages(),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
)