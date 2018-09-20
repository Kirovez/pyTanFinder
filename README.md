# pyTanFinder

pyTanFinder is a python tool to search tandem repeats in genomic or long-read data by running, parsing and processing [TandemRepeatFinder (Benson, 1999)](https://tandem.bu.edu/trf/trf.html) output

The program was developed in frame of the RNF grant (№17-14-01189)

Our website: []()
## Getting Started

Download pyTanFinder from [Github repository](https://github.com/Kirovez/pyTanFinder) and proceed with installation instruction

### Prerequisites

pyTanFinder requires:
* TRF executable file ([Can be donwload here](https://tandem.bu.edu/trf/trf.download.html)). 
For linux system it has to be made executable before run (e.g. run command `chmod +x trf409.linux64`)
* balstn and makeblastdb programs. The paths to these programs can be set via `-bp` and `-mp` flags, respectively
* python v3.6
python packages to be installed: biopython, networkx
(run command: `pip install biopython networkx`)
* for Ubuntu tkinter is required (run command `sudo apt-get install python3-tk`)

## Running pyTanFinder


## Authors

* **Ilya Kirov** 


## License

This project is licensed under the MIT License
