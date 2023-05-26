import os
from distutils.core import Extension, setup


setup(
    name="fcdimen",
    version="0.1.0",
    description="Analyzing dimensionality of materials structure using force constants",
    author="Mohammad Bagheri, Ethan Berger, Hannu-Pekka Komsa",
    author_email="Mohammad.Bagheri@oulu.fi",
    url="https://github.com/FCDimen/",
    install_requires=["numpy >= 1.22.3","phonopy >= 2.8.1", "matplotlib >= 2.2.2", "ase >= 3.22.1", "networkx >= 2.7.1"],
    packages=[
        "fcdimen",
        "fcdimen.functions",
    ],
    scripts=[
        "scripts/fcdimen",
    ],
    classifiers=[
          "Development Status :: 1 - Beta",
          "License :: OSI Approved :: BSD License",
          "Intended Audience :: Science/Research",
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
