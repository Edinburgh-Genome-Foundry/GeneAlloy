import ez_setup

from setuptools import setup, find_packages

ez_setup.use_setuptools()

version = {}
with open("genealloy/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="genealloy",
    version=version["__version__"],
    author="Peter Vegh",
    description="GeneAlloy helps designing overlapping sequences.",
    long_description=open("README.md").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/GeneAlloy",
    keywords="biology",
    packages=find_packages(exclude="docs"),
    install_requires=["biopython"],
)
