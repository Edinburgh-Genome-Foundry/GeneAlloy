import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('genealloy/version.py').read()) # loads __version__

setup(name='genealloy',
      version=__version__,
      author='Peter Vegh',
    description='GeneAlloy helps designing overlapping sequences.',
    long_description=open('README.md').read(),
    license='MIT',
    keywords="biology",
    packages= find_packages(exclude='docs'),
    install_requires= ['biopython'])
