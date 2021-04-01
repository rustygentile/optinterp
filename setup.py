from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='optinterp',
    version='0.0.5',
    description='Optimal Interpolation Nodes',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/rustygentile/optinterp',
    author='Rusty Gentile',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords=['interpolation', 'approximation', 'curve fitting'],
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.6, <4',
    install_requires=['scipy==1.5.3'],
    project_urls={
        'Source': 'https://github.com/rustygentile/optinterp',
    },
)
