import os
class TRF_run_parse():

    def __init__(self,file_start, outFasta, runTRF,trf_path,
                 name_for_trf, minMonLength, maxMonLength,
                 minMonNum, prefix,separator_in_sequence):
        self.file_start = file_start
        self.runTRF = runTRF
        self.trf_path = trf_path
        self.name_for_trf = name_for_trf
        self.outFasta = outFasta
        self.minMonLength, self.maxMonLength, self.minMonNum = minMonLength, maxMonLength, minMonNum
        self.prefix = prefix
        self.separator_in_sequence = separator_in_sequence

    def run(self):
        if self.runTRF:
            cmd = r'{0} {1} 2 7 7 80 10 20 {2} -f -h -m'.format(self.trf_path, self.file_start,self.maxMonLength)
            print(cmd)
            os.system(cmd)

        with open(self.name_for_trf) as fasta, open(self.outFasta, 'w') as out:
            tr_number = 0
            seq_id = []
            q = 0
            a = [str(line) for line in fasta]

            for lin in a:
                # identification of start of the sequence
                # name of the sequence
                if lin[0] == 'S':
                    seq_id = lin.split(' ')[1].rstrip()

                # identification of the sequence
                ## if line starts with integer then it corresponds to TR line
                check = lin.split(' ')[0]
                if check.isdigit() == True:  # m != -1:
                    split_lin = lin.split(' ')
                    tr_number += 1  # number of the sequence
                    seq_seq = split_lin[13]  # lin[m:]  # sequence itself from start to the end
                    s = str(tr_number)  # number of the sequence
                    Rn = split_lin[3]  # time of monomer repeats

                    if self.maxMonLength >= len(seq_seq) >= self.minMonLength and \
                                    float(Rn) >= self.minMonNum:  # length and number of repeat filter

                        q += 1  # showed number of sequences
                        # print('>' + self.prefix + seq_id + self.separator_in_sequence +
                        #       s + self.separator_in_sequence +
                        #       str(Rn) + self.separator_in_sequence + str(
                        #     len(seq_seq)) + '\n' + seq_seq)
                        out.write('>' + self.prefix + seq_id + self.separator_in_sequence + s + self.separator_in_sequence + str(Rn) + self.separator_in_sequence + str(
                            len(seq_seq)) + '\n' + seq_seq + '\n')

        print('Total number of sequences found by TRF: ', str(tr_number))
        print('Number of sequences did not meet the filtration criterion: ' + str(tr_number - q))
        print('Number of sequences in the output file {}: '.format(self.outFasta) + str(q))

