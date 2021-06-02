from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='openmm-toolkit',
    version='0.1.0',
    author='Zhenyu Wei',
    author_email='zhenyuwei99@gmail.com',
    description='This repo contains samplers can be used together with openmm',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    url='https://github.com/zhenyuwei/openmm-toolkit',
    project_urls={
        "Source Code": "https://github.com/zhenyuwei/openmm-toolkit",
    },
    packages=find_packages(),
    package_data={
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=[
        'pytest',
        'pytest-xdist',
        'nvgpu'
    ],
    python_requires='>=3.7'
)
