# pyTanFinder


![pyTanFinder](https://github.com/Kirovez/pyTanFinder/blob/master/picture.png)


pyTanFinder is a python tool to search tandem repeats in genomic or long-read data by running, parsing and processing [TandemRepeatFinder (Benson, 1999)](https://tandem.bu.edu/trf/trf.html) output.

The pyTanFinder paper has been [published] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6302065/). Please cite:

Kirov, I., Gilyok, M., Knyazev, A., & Fesenko, I. (2018). Pilot satellitome analysis of the model plant, Physcomitrellapatens, revealed a transcribed and high-copy IGS related tandem repeat. Comparative Cytogenetics, 12(4), 493.

The program was developed in frame of the RNF grant (№17-14-01189). Our website: [Plantprotlab.com](http://plantprotlab.com/)


## Pipeline
pyTanFinder is user friendly command line tool to run TRF software and parse the results followed by clustering of similar tandem repeats. The output of this program is a fasta file of all tandem repeats and table containing unique TR sequences with the estimated abundancy in genome. In addition, pyTanFinder also genearates html report containing the histograms of distribution of TR monomer size and number of connections of each monomer into individual cluster. 

![pipeline](https://github.com/Kirovez/pyTanFinder/blob/master/pipeline.png)


## Getting Started

Download pyTanFinder from [Github repository](https://github.com/Kirovez/pyTanFinder) and proceed with installation instruction

### Prerequisites

pyTanFinder requires:
* TRF executable file ([Can be donwload here](https://tandem.bu.edu/trf/trf.download.html)). 
For linux system it has to be made executable before run (e.g. run command `chmod +x trf409.linux64`)
* blastn and makeblastdb programs. The paths to these programs can be set via `-bp` and `-mp` flags, respectively
* python v3.6
python packages to be installed: biopython, networkx
(run command: `pip install matplotlib biopython networkx`)
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

**-minM** - minimum length for tandem repeat

**-maxM** - maximum length for tandem repeat

**-minMN** - minimum number of repetitions

**-minA** - minimum abundancy of a repeat in genome, bp

**-px** - prefix that will be used for file and sequence names

**--no_blast** - do not run blast (useful when it already was run before)

**--no_html** - do not make html report

**--no_runTRF** - do not run TRF (useful when it already was run before)

**-tp** - path to TRF executable

**-bp** - path to blastn executable

**-mp** - path to makeblastdb executable


## Authors

**Ilya Kirov** 


## License

This project is licensed under the MIT License
