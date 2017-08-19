import os.path
from setuptools import setup, find_packages


script_path = os.path.dirname(__file__)
with open(os.path.join(script_path, 'README.md')) as f:
    readme = f.read()

with open(os.path.join(script_path, 'LICENSE')) as f:
    license = f.read()

setup(
    name='diana',
    version='3.9',
    description='Package to predict computationally drug combinations',
    long_description=readme,
    author='Joaquim Aguirre-Plans',
    author_email='joaquim.aguirre@upf.edu',
    url='',
    license=license,
    packages=find_packages(exclude=('diana'))
)