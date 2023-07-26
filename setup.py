from setuptools import dist, setup, Extension

# bootstrap numpy; can we workaround this? 
dist.Distribution().fetch_build_eggs(["numpy>=1.14.5"])

# should be fine now
import numpy as np


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="qsp",
    version="0.1",
    author="Caoqi Fan",
    author_email="fcq1116@gmail.com",
    description="Quantitative Systems Pharmacology",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Ephraim-usc/qsp.git",
    packages=["qsp"],
    python_requires=">=3",
    install_requires=[
        "numpy>=1.14.5",
        "pandas",
    ],
    ext_modules=[],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
