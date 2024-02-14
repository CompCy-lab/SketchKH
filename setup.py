from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='sketchKH',
    version='0.1.1',
    description='Distribution-based sketching of single-cell samples',
    author='CompCy Lab',
    author_email='compcylab@gmail.com',
    url='https://github.com/CompCy-lab/SketchKH',
    license='MIT',
    python_requires='>=3.6',
    long_description=long_description,
    long_description_content_type = 'text/markdown',
    packages = find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
    ],
    keywords=[
        'cytometry',
        'single-cell bioinformatics',
        'clinical prediction',
        'compression',
        'computational biology',
    ],
    install_requires=[
        'anndata>=0.7.6',
        'numpy>=1.22.4',
        'scipy>=1.7.1',
        'tqdm',
    ],
    ext_modules=[],
)