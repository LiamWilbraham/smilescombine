import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="smilescombine",
    version="0.0.1",
    author="Liam Wilbraham",
    author_email="l.wilbraham@ucl.ac.uk",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LiamWilbraham/smilescombine",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ),
)
