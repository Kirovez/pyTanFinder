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
* blastn and makeblastdb programs. The paths to these programs can be set via `-bp` and `-mp` flags, respectively
* python v3.6
python packages to be installed: biopython, networkx
(run command: `pip install biopython networkx`)
* for Ubuntu tkinter is required (run command `sudo apt-get install python3-tk`)

## Running pyTanFinder
usage: pyTanFinder.py fasta [-h] [-minM MINMONLENGTH] [-maxM MAXMONLENGTH]
                      [-minMN MINMONNUM] [-minA MINABUNDANCY] [-px PREFIX]
                      [--no_blast] [--no_html] [--no_runTRF] [-tp TRF_PATH]
                      [-bp BLAST_PATH] [-mp MAKE_BLAST]

positional arguments:

**fasta**  - fasta file where tandem repeats will be found

optional arguments:

**-h, --help** show this help message and exit

**-minM** --minMonLength minimum length for tandem repeat

  -maxM MAXMONLENGTH, --maxMonLength MAXMONLENGTH
                        maximum length for tandem repeat
  -minMN MINMONNUM, --minMonNum MINMONNUM
                        minimum number of repetitions
  -minA MINABUNDANCY, --minAbundancy MINABUNDANCY
                        minimum abundancy of a repeat in genome, bp
  -px PREFIX, --prefix PREFIX
                        prefix that will be used for file and sequence names
  --no_blast            do not run blast (useful when it already was run
                        before)
  --no_html             do not make html report
  --no_runTRF           do not run TRF (useful when it already was run before)
  -tp TRF_PATH, --trf_path TRF_PATH
                        path to TRF executable
  -bp BLAST_PATH, --blast_path BLAST_PATH
                        path to blastn executable
  -mp MAKE_BLAST, --make_blast MAKE_BLAST
                        path to makeblastdb executable


## Authors

* **Ilya Kirov** 


## License

This project is licensed under the MIT License
