""" setup module """

from setuptools import setup, find_packages
import os.path

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pymzid',
    version='0.3.2',
    description='pymzid reads in mzidentml files from proteomics search engine',

    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/Lau-Lab/pymzid/',

    author='Edward Lau',
    author_email='edward.lau@cuanschutz.edu',

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='scientific proteomics mass-spectrometry mzid mzidentml',  # Optional

    packages=find_packages(),

    python_requires='>=3.6, <4',

    install_requires=['pandas>=1,<2',
                      'tqdm>=4,<5'
                      ],  # external packages as dependencies

    entry_points={
        'console_scripts': [
            'pymzid=pymzid.__main__:main',
        ],
    },

    project_urls={
        'Source': 'https://github.com/Lau-Lab/pymzid',
        'Edward Lau Lab': 'https://www.laulab.net',
    },

    data_files=[('tests',
                 [os.path.join('tests', 'data', 'comet_percolator', 'percolator.target.mzid'),
                  os.path.join('tests', 'data', 'msgfplus', '20180216_BSA.mzid'),
                  ]),
                ],

)
