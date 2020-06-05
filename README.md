# Pymzid

Reads in mzid files from protein identification results in mass spectrometry/proteomics experiments

## Usage

Launch as a standalone:

    $ python -m pymzid path/to/mzid -o path/to/out
    
Use as a module:

    from pymzid.read_mzid import Mzid
    
    mzid = Mzid("path/to/mzid")
    mzid.read_psm()
    
gives a pandas object under mzid.psm_df
   
To run the test Percolator data and print the output to home:

    $ python -m pymzid tests/data/comet_percolator/percolator.target.mzid -o ~  
    

## Contributing

Please contact us if you wish to contribute, and submit pull requests to us.


## Authors

* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)

See also the list of [contributors](https://github.com/Lau-Lab/pymzid/graphs/contributors) who participated in this project.


## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/Lau-Lab/pymzid/blob/master/LICENSE.md) file for details

