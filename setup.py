import setuptools
import os

# Load package info from __version__.py
about = dict()
here = os.path.abspath(os.path.dirname(__file__))
pkg_path = os.path.join(here, "lymphgenerator")
with open(os.path.join(pkg_path, "__version__.py")) as f:
    exec(f.read(), about)

setuptools.setup(
    name = about["__title__"],
    version = about["__version__"],
    description = about["__description__"],
    url = about["__url__"],
    author = about["__author__"],
    package_dir = {"lymphgenerator": pkg_path},
    install_requires = ['polars', 'SigProfilerAssignment', 'seaborn'],
    author_email = about["__author_email__"],
    license = about["__license__"],
    packages = setuptools.find_packages(),
    zip_safe = False
)
