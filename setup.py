from setuptools import setup, find_packages
from setuptools.extension import Extension


with open('README.md') as infile:
    long_description = infile.read()


setup(
    name='BioPandas',
    version='1.1.0',  # major.minor.maintenance
    description='Import genomic data to get a custom Pandas & Biopython hybrid class object with fancy shortcuts to make Machine Learning preprocessing easy!',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/tmsincomb/BioPandas',
    author='Troy Sincomb',
    author_email='troysincomb@gmail.com',
    license='MIT',
    keywords='bio biopandas',
    packages=find_packages('BioPandas'),
    # include_package_data=True,  # try this out: might be the reason packages didnt break since it wont run without this.
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 1 - ALPHA',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    #  TODO: add classifiers for machine learning/variant caller https://pypi.org/classifiers/
    install_requires=[
        'biopython',
        'pandas',
        'pysam',
    ],
    # entry_points={
    #     'console_scripts': [
    #         'biopandas=BioPandas.core:main',
    #     ],
    # },
)