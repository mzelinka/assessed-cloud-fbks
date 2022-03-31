import glob
import os
import subprocess
import sys

from setuptools import find_packages, setup

packages = find_packages()
scripts = ["assessed_cloud_fbks/main.py"]
data_files = glob.glob("data/*.json")

setup(
    name="assessed_cloud_fbks",
    version="1",
    author="mzelinka",
    description="cloud feedbacks",
    long_description="",
    url="https://github.com/mzelinka/assessed-cloud-fbks",
    packages=packages,
    scripts=scripts,
    data_files=data_files
)
