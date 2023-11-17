from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="funstools",  # Replace with your own username
    version="0.6.2",
    author="Shinyoung Kim",
    author_email="radioshiny@gmail.com",
    description="FUNStools is python toolkit for FUNS project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/radioshiny/funstools",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.9",
    ],
    python_requires='>=2.7',
)
