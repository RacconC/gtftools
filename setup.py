from setuptools import setup, find_packages
with open('README.md','r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name="gtftools",
    version="0.9.0",
    author="Hong-Dong Li",
    author_email="hongdong@csu.edu.cn",
    description="gtftools provides a set of functions to compute or extract various features of gene models.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['gtftools = gtftools.gtftools:main']
    },
    package_data = {'gtftools': ['test_data/*.txt', 'test_data/*.gtf'],},
    zip_safe = False
)