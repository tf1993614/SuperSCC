import sys

if sys.version_info < (3, 6):
    sys.exit("SuperSCC requires Python >= 3.6")

from setuptools import setup, find_packages
from pathlib import Path

version = {}
with open("SuperSCC/_version.py") as fp:
    exec(fp.read(), version)

setup(
    name="SuperSCC",
    version=version["__version__"],
    author="Feng Tang",
    author_email="fengtang614@qq.com",
    license="BSD",
    description="SuperSCC: Super single cell clustering",
    long_description=Path("README.md").read_text("utf-8"),
    long_description_content_type="text/markdown",
    url="https://github.com/tf1993614/SuperSCC",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
       "dcor==0.6", 
       "scanpy==1.10.1", 
       "pandas==2.2.2", 
       "numpy==1.26.4", 
       "scipy==1.13.0", 
       "scikit-learn==1.2.2", 
       "igraph==0.11.4", 
       "leidenalg==0.10.2", 
       "plotly==5.22.0", 
       "rpy2==3.5.17", 
       "langchain",
       "langchain-community"
    ]
)
