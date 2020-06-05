# Pymzid

Reads in mzid files from protein identification results in mass spectrometry/proteomics experiments

## Installation and Usage

Install Python 3.7+ and pip. See instructions on Python website for specific instructions for your operating system.

Pymzid can be installed from PyPI via pip. We recommend using a virtual environment.

    $ pip install pymzid

Launch as a standalone:

    $ python -m pymzid path/to/mzid -o path/to/out

Alternatively:
    
    $ pymzid
    
Use as a module:

    from pymzid.read_mzid import Mzid
    
    mzid = Mzid("path/to/mzid")
    mzid.read_psm()
    
gives a pandas object under mzid.psm_df

To test that the installation can load test data files in tests/data:

    $ pip install tox
    $ tox
   
To run the test Percolator data and print the output to home:

    $ python -m pymzid tests/data/comet_percolator/percolator.target.mzid -o ~  
    
### Dependencies

Pymzid is tested in Python 3.7 and 3.8 and uses the following packages:

```
pandas==1.0.4
tqdm==4.46.1
```


## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)

See also the list of [contributors](https://github.com/Lau-Lab/pymzid/graphs/contributors) who participated in this project.


## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/Lau-Lab/pymzid/blob/master/LICENSE.md) file for details

