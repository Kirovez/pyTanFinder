"""
This script takes fasta file from TRFparserCMD
and clusters similar TRs and choose one fasta sequence per cluster
return fasta of the most representable Rs and table
"""

from Bio import SeqIO
from collections import defaultdict
import os
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from matplotlib import pyplot as plt
from networkx import *
from HTML_writer import HTMLreport

class Graph_selected():
    def __init__(self):
        self.graphs = []

class TRFmerger():
    def __init__(self, trffasta, prefix,
                 outputDir, blast_path,
                 make_blast, merger_out_dir,
                    separator_in_sequence,
                 minAbundancy,
                 do_blast=True,
                 evalue = 1e-5,
                 coverage = 80,
                 similarity = 80,
                 writeHTML = False):
        self.trffasta = trffasta
        self.bp = blast_path
        self.bl = do_blast
        self.mb = make_blast
        self.parsingBl = False
        self.prefix = prefix
        self.writeHTML = writeHTML
        self.rootout = outputDir
        self.similarity = similarity
        self.evalue = evalue
        self.coverage = coverage
        self.minAbundancy = minAbundancy
        self.seq_no_hits = []
        self.merger_out_dir = merger_out_dir
        self.separator_in_sequence = separator_in_sequence
        self.blast_out = '{0}/{1}_blast.out.xml'.format(self.merger_out_dir, self.prefix)
        self.blast_mask_out = '{0}/parsed_BLAST_results_{1}'.format(self.merger_out_dir, self.prefix)
        self.query_in_fasta  = '{0}/query_in_fasta_{1}'.format(self.merger_out_dir, self.prefix)
        self.no_hits = '{0}/NO_hits_{1}'.format(self.merger_out_dir, self.prefix)
        #self.clans = '{0}/CLANS_{1}'.format(self.merger_out_dir, self.prefix)
        self.abundancy = "{0}/RepeatAbundance_{1}.txt".format(self.merger_out_dir, self.prefix)
        self.results = '{0}/fasta_{1}_after_clustering.fasta'.format(self.merger_out_dir, self.prefix)
        self.table_out = '{0}/after_clusterin_summary_{1}.txt'.format(self.merger_out_dir, self.prefix)
        self.html_file = "{0}/Clustering/Clustering_BLAST_data_{1}.html".format(self.merger_out_dir, self.prefix)
        self.HTMLrep = HTMLreport("{0}/Clustering/imgs".format(self.merger_out_dir),
                                  self.html_file)


    def getTandemLEngthInContig(self,id):
        """
        :param id: tandem repeat id e.g AeTaMCGU01017908.1_276934_6.0_75
        :return:multiply length and number of repeats
        """
        sp = id.split(self.separator_in_sequence)
        return int(sp[-1]) * float(sp[-2])

    def BLAST(self):
        if self.bl:
            ################make blast db#############################
            cmd = r'{0} -in {1} -dbtype nucl'.format(self.mb, self.trffasta)
            os.system(cmd)
            ################run blast#############################
            print("BLAST is running...")
            blast = NcbiblastnCommandline(self.bp, query=self.trffasta, db=self.trffasta, out=self.blast_out, outfmt=5, word_size = 11, evalue = 0.001)
            stdout, stderr = blast()
            print("BLAST finished. \n Parsing ......")

    def getCoverage(self, seq_len, intervals):
        if not intervals:
            return 0
        intervals.sort(key=lambda interval: interval[0])
        merged = [intervals[0]]
        for current in intervals:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)

        covered_bp = 0
        for merged_intervals in merged:
            covered_bp += abs(merged_intervals[0] - merged_intervals[1])

        return covered_bp *100/seq_len


    def progress(self, count, total, status=''):
        bar_len = 100
        filled_len = int(round(bar_len * count / float(total)))
        percents = round(100.0 * count / float(total), 1)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        print('[%s] %s%s ...%s\r' % (bar, percents, '%', status), end='')

    def parseBLAST(self):
        total_bases = defaultdict(dict)
        no_hits_ind = 0
        query_indexed_db = SeqIO.index(self.trffasta, 'fasta')
        print('Number of query sequences: ' + str(len(query_indexed_db)))
        query_list = []
        list_sbjct = []
        g=0
        cal = 0
        query_count = 0
        number_of_quries_in_fasta_file =0
        list_of_printed=[]

        print("BLAST filtering starts.....")
        with open (self.blast_mask_out, 'w') as BLAST_MASK, \
                open(self.query_in_fasta,'w') as query_in_fasta, \
                open(self.no_hits, 'w') as no_hits, \
                open(self.abundancy, "w") as abundance:

            file2 = open(self.blast_out)
            s1 = SearchIO.parse(file2, 'blast-xml')
            count_hits = 0
            seq_num = len(query_indexed_db)
            for recor in s1:
                total_bases[recor.id][recor.id] = self.getTandemLEngthInContig(recor.id)
                number_of_hits=0
                query_count+=1
                sequence = str(query_indexed_db[recor.id].seq)

                self.progress(query_count, seq_num)


                for HSP in recor:
                    indicator = 0
                    hs = HSP.hsps
                    cov_intervals = []
                    for u in hs:
                        if u.evalue<=self.evalue and u.ident_num*100/len(u.query) >= self.similarity:
                            cov_intervals.append([u.query_range[0],u.query_range[1]])
                            #sequence = sequence.replace(sequence[u.query_range[0]:u.query_range[1]], 'n' * len(u.query))
                            indicator += 1
                            g+=1
                            count_hits+=1
                            cal+=1
                    #coverage = sequence.count('n')*100/len(sequence)
                    #print(coverage, self.getCoverage(len(sequence), cov_intervals))
                    #Coverage calculation based on masking results
                    if self.getCoverage(len(sequence), cov_intervals)>=self.coverage and indicator >= 1:
                        number_of_hits+=1
                        BLAST_MASK.write('{0}{1}{2}\n'.format(recor.id,self.separator_in_sequence,str(HSP.id)))
                        # for_clans = recor.id + ' ' + HSP.id + ':'+str(e_value) + '\n'
                        # clans.write(for_clans)
                        if HSP.id != recor.id and HSP.id not in list_sbjct:
                                list_sbjct.append(HSP.id)
                        if recor.id not in list_of_printed:
                            number_of_quries_in_fasta_file +=1
                            SeqIO.write(query_indexed_db[recor.id],query_in_fasta,'fasta')
                            list_of_printed.append(recor.id)
                            query_list.append(recor.id)

                        ### calculate coverage of hit and add to query coverage
                        if recor.id != HSP.id:
                            total_bases[recor.id][HSP.id] = 0.0
                            total_bases[recor.id][HSP.id] = self.getTandemLEngthInContig(HSP.id)
                if number_of_hits==0:
                    no_hits_ind+=1

            ### OUTPUT TABLE WRITING ###
            abundance.write("TR id \t TR total abundancy \t Number of TR similar \n")

            # NO hits file writing
            for seq in SeqIO.parse(self.trffasta,'fasta'):
                if seq.id not in query_list:
                    SeqIO.write(seq,no_hits,'fasta')
                    abundance.write(seq.id + "\t" + str(self.getTandemLEngthInContig(seq.id)) + "\t" + '0' + "\n")
                    self.seq_no_hits.append(seq)

            for tandem in total_bases:
                abund = round(sum([total_bases[tandem][i] for i in total_bases[tandem]]))
                number_merged_tandems = len(total_bases[tandem])
                abundance.write(tandem + "\t" + str(abund) + "\t" + str(number_merged_tandems) + "\n")


        print('Hit pairs in HIT file: ' + str(cal))
        print('Number of query sequences in fasta file: '+str(number_of_quries_in_fasta_file))
        print('Number of queries with no hits: ' + str(no_hits_ind))
        print('Number of different sbjct sequences ' + str(len(list_sbjct)))

    def getGraph(self):
        print("Graph generation....")
        edge = []
        with open(self.blast_mask_out) as filtered_blast:
            for lines in filtered_blast:
                lines = lines.rstrip()
                edge.append((lines.split(self.separator_in_sequence + self.prefix)[0],
                             self.prefix + lines.split(self.separator_in_sequence + self.prefix)[1]))
        G = nx.Graph()
        G.add_edges_from(edge)
        return G

    def merger(self):
        cnt=-1
        monomer_length = []
        repeats = []
        in_fasta = 0
        query_indexed_db1  = SeqIO.index(self.trffasta,'fasta')
        list_cluster_len = []

        G = self.getGraph() # it opens mask file and generates Graph

        with open (self.results, 'w') as res ,\
        open(self.table_out, "w") as tabOut, \
                open (self.html_file, "w") as htmlfile:
            tabOut.write("Cluster\tID\tSequence\tMonomer length\tAccumulated abundance\tNumber of connections\n")
            connected_components = nx.connected_components(G)
            print("Number of clusters: ", number_connected_components(G))
            # iterate through clusters
            for num, clusters in enumerate(connected_components):
                with open("{0}\Clustering\CL{1}.fasta".format(self.merger_out_dir, num), "w") as outfile:
                    list_cluster_len.append(len(clusters))
                    cnt+=1

                    accumul_abundace = 0.0
                    monmer_length = []
                    degree_list = []

                    # selection of the node in the cluster with max degree
                    node_list = []

                    # sort sequences in cluster by the length from small to big
                    clusters = sorted(clusters, key=lambda x: int(x.split(self.separator_in_sequence)[-1]))

                    for n in clusters:
                        accumul_abundace += self.getTandemLEngthInContig(n)
                        degree_list.append(degree(G)[n])
                        node_list.append(n)
                        monmer_length.append(float(n.split(self.separator_in_sequence)[3]))
                    # write sequence node with max degree
                    max_degree = max(degree_list)

                    if accumul_abundace >= self.minAbundancy:
                        for node_num, i1 in enumerate(clusters):
                            if degree(G)[i1] == max_degree:

                                if self.writeHTML:
                                    cluster_name = "CL{}".format(num)
                                    print("Generation of html for",cluster_name)
                                    sequence_name = i1
                                    self.HTMLrep.fun_write_HTML(cluster_name,
                                           sequence_name,
                                           node_num,
                                           monmer_length,
                                           max_degree,
                                           len(node_list),
                                           degree_list,
                                            accumul_abundace)

                                sp = i1.split(self.separator_in_sequence)
                                moomer_len = float(sp[3].split("##")[0])
                                monomer_length.append(moomer_len)
                                repeats.append(int(float(sp[2])))

                                seqTR = query_indexed_db1[i1]
                                seqTR.description += " accumulated_abundance: {} ".format(accumul_abundace) + " connections: {}".format(degree(G)[i1])

                                SeqIO.write(seqTR,res,'fasta')

                                tabOut.write("CL{}".format(num) + "\t" +
                                             seqTR.id +
                                             "\t" + str(seqTR.seq) +
                                             "\t" + str(moomer_len).replace(".",",") \
                                             + "\t" + str(accumul_abundace).replace(".",",") + "\t" \
                                             + str(degree(G)[i1]) + "\t" + ",".join(node_list) + "\n")

                                SeqIO.write(seqTR,outfile,"fasta")
                                for seq in node_list:
                                    if seq != seqTR.id:
                                        SeqIO.write(query_indexed_db1[seq],outfile,"fasta")
                                in_fasta+=1
                                break

            for no_hits_seq in self.seq_no_hits:
                tabOut.write("CL-" + "\t" +
                             no_hits_seq.id +
                             "\t" + str(no_hits_seq.seq) +
                             "\t" + str(len(no_hits_seq.seq)) \
                             + "\t" + str(self.getTandemLEngthInContig(no_hits_seq.id)) + "\t" \
                             + '0' + "\t" + '-' + "\n")
                SeqIO.write(no_hits_seq, res, 'fasta')


        print('in fasta file: ' + str(in_fasta))

        plt.plot(monomer_length,repeats,'ro')
        plt.ylabel('Number of repeats in contig')
        plt.xlabel('Monomer length')
        plt.savefig(self.merger_out_dir + '/length_vs_repeats.png',format='png')

        print('Minimum monomer length: ', min(monomer_length))
        print('Maximum monomer length: ',max(monomer_length))

        self.HTMLrep.writeHeaderFooter(header=False)

    def merge(self):
        #self.BLAST()
        #self.parseBLAST()
        self.merger()
