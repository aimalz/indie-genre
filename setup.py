from setuptools import setup

setup(
    name="pzgen",
    version="0.1",
    author="Alex Malz",
    author_email="aimalz@nyu.edu",
    url = "https://github.com/aimalz/pzgen",
    packages=["pzgen"],
    description="Generic forward model for probability spaces of true and estimated scalars",
    long_description=open("README.md").read(),
    package_data={"": ["README.md", "LICENSE"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU GPL v3.0 License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        ],
    install_requires=["matplotlib", "numpy", "pomegranate"]
)
