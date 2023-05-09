import os
from distutils.core import Extension, setup


setup(
    name="fcdimen",
    version="0.1.0",
    description="Analyzing dimensionality of materials structure with force constants matrix",
    author="Mohammad Bagheri",
    author_email="Mohammad.Bagheri@oulu.fi",
    url="https://github.com/FCDimen/",
    install_requires=["numpy","phonopy", "matplotlib", "ase", "networkx"],
    packages=[
        "fcdimen",
        "fcdimen.functions",
    ],
    scripts=[
        "scripts/fcdimen",
    ],
    classifiers=[
          "Development Status :: 1 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: Creative Commons Attribution-ShareAlike 4.0 International License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11",
          "Topic :: Scientific/Engineering :: Materials",
          "Topic :: Scientific/Engineering :: Physics",
          ],
)
