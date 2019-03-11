"""
1. This script allows to run TRF program (Bennetsen, 1999)
2. It parses the results and write fasta file of all tandem repeat
3. CLustering the results and select the most represented TR sequence
"""
from pathlib import Path
import os
from TRF_run_parse import TRF_run_parse
from TRF_merger import TRFmerger
from distutils.spawn import find_executable


class pyTanFinder():
    def __init__(self, fastaFile, minMonLength = 30, maxMonLength = 1000, minMonNum = 5, minAbundancy = 10000,
                 prefix='XX', do_blast=True, writeHTML = True, runTRF=True,
                 separator_in_sequence = "<*>",
                 trf_path='./trf409.linux64',
                 blast_path = r'blastn',
                 make_blast = r'makeblastdb'):
        self.file_start = fastaFile
        self.minMonLength = minMonLength
        self.maxMonLength = maxMonLength
        self.minMonNum = minMonNum
        self.minAbundancy = minAbundancy
        self.prefix = prefix
        self.runTRF = runTRF
        self.blast_path = blast_path
        self.do_blast = do_blast
        self.make_blast = make_blast
        self.writeHTML = writeHTML
        self.trf_path = trf_path
        self.separator_in_sequence = separator_in_sequence
        self.rootOut = "out_pyTanFinder"

        self.name_for_trf = '{0}.2.7.7.80.10.20.{1}.dat'.format(self.getFileName(self.file_start), self.maxMonLength)  # the numbers should be changed if other parametres will be used
        self.outFasta = '{0}/FILTERED_tandem_repeats_{1}.fasta'.format(self.rootOut, prefix)

        ####merger files#####
        self.merger_out_dir = '{}/out_merger'.format(self.rootOut)
        self.outputDir = "{0}/merger_output_{1}".format(self.merger_out_dir, self.prefix)
        self.clustering_dir = "{0}/Clustering".format(self.outputDir)
        self.imgs_folder = "{0}/Clustering/imgs".format(self.outputDir)
        self.__checkPars()
        self.main()
        
    def getFileName(self,file):
        if '/'in file:
            return file.split('/')[-1]
        else:
            return file
        
    def __checkPars(self):
        if not self.__checkFileExist(self.file_start):
            raise IOError("FASTA file was not found!")
        if not self.__checkFileExist(self.blast_path):
            raise IOError("BLAST executable was not found!")
        if not self.__checkFileExist(self.make_blast):
            raise IOError("makeblastdb executable was not found!")
        if not self.__checkFileExist(self.trf_path):
            raise IOError("TRF executable was not found!")
        if not self.__checkFolderExist(self.rootOut):
            os.makedirs(self.rootOut)
        if not self.__checkFolderExist(self.merger_out_dir):
            os.makedirs(self.merger_out_dir)
        if not self.__checkFolderExist(self.outputDir):
            os.mkdir(self.outputDir)
        if not self.__checkFolderExist(self.clustering_dir):
            os.mkdir(self.clustering_dir)
        if not self.__checkFolderExist(self.imgs_folder):
            os.mkdir(self.imgs_folder)


    def __checkFileExist(self, outputFilename):
        my_file = Path(outputFilename)
        if my_file.is_file():
            return True
        #Check whether `name` is on PATH.
        elif find_executable(outputFilename) is not None:
            return True
        return False

    def __checkFolderExist(self, folder):
        my_folder = Path(folder)
        if my_folder.is_dir():
            return True
        return False

    def main(self):
        # run TRF if it is required and parse the result file
        # output id fasta file with selected monomers
        trf_run = TRF_run_parse(self.file_start, self.outFasta,
                                self.runTRF, self.trf_path,
                                self.name_for_trf, self.minMonLength,
                                self.maxMonLength, self.minMonNum,
                                self.prefix,
                                self.separator_in_sequence)
        trf_run.run()


        trf_merger = TRFmerger(self.outFasta, self.prefix,
                                self.rootOut, self.blast_path,
                                self.make_blast, self.outputDir,
                               self.separator_in_sequence,
                               self.minAbundancy,
                               do_blast=self.do_blast,
                                writeHTML = self.writeHTML)
        trf_merger.merge()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='fasta file where tandem repeats will be found')
    parser.add_argument('-minM', '--minMonLength', help='minimum length for tandem repeat', default=20)
    parser.add_argument('-maxM', '--maxMonLength', help='maximum length for tandem repeat', default=2000)
    parser.add_argument('-minMN', '--minMonNum', help='minimum number of repetitions', default=5)
    parser.add_argument('-minA', '--minAbundancy', help='minimum abundancy of a repeat in genome, bp', default=10000)
    parser.add_argument('-px', '--prefix', help='prefix that will be used for file and sequence names',default="XXX")
    parser.add_argument('--no_blast', help='do not run blast (useful when it already was run before)', action='store_true')
    parser.add_argument('--no_html', help='do not make html report', action='store_true')
    parser.add_argument('--no_runTRF', help='do not run TRF (useful when it already was run before)', action='store_true')
    parser.add_argument('-tp', '--trf_path', help='path to TRF executable', default='./trf409.linux64')
    parser.add_argument('-bp', '--blast_path', help='path to blastn executable', default='blastn')
    parser.add_argument('-mp', '--make_blast', help='path to makeblastdb executable', default='makeblastdb')

    pars = parser.parse_args()

    pyTanFinder(pars.fasta, prefix=pars.prefix,
            do_blast=pars.no_blast==False,
                runTRF=pars.no_runTRF==False,
                writeHTML=pars.no_html == False,
                minMonLength = int(pars.minMonLength),
                maxMonLength = int(pars.maxMonLength),
                minMonNum = int(pars.minMonNum),
                minAbundancy = int(pars.minAbundancy),
                trf_path = pars.trf_path,
                blast_path = pars.blast_path,
                make_blast = pars.make_blast)
