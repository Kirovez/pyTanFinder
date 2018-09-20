
from matplotlib import pyplot as plt

class HTMLreport():

    def __init__(self, imgs_dir, out_html_file ):
        self.out_html_file = open(out_html_file, "w")
        self.imgs_dir = imgs_dir
        self.writeHeaderFooter(header=True)

    def monomerLength(self,id):
        return float(id.split("_")[3].split("##")[0])

    def writeHeaderFooter(self, header=True):
        htmlfile = self.out_html_file
        if header:
            htmlfile.write("<html>\n")
            htmlfile.write("<head>\n")
            htmlfile.write(
                '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta.2/css/bootstrap.min.css" integrity="sha384-PsH8R72JQ3SOdhVi3uxftmaW6Vc51MKb0q5P2rRUpPvrszuE4W1povHYgTpBfshb" crossorigin="anonymous">\n')
            htmlfile.write(
                '<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>\n')
            htmlfile.write(
                '<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta.2/js/bootstrap.min.js" integrity="sha384-alpBpkh1PFOepccYVYDB4do5UnbKysX5WZXm3XxPqe5iKTfUKjNkCk9SaVuEZflJ" crossorigin="anonymous"></script>\n')
            htmlfile.write("</head>\n")
            htmlfile.write("<body>\n")

            htmlfile.write("<div class='container'>\n")
        else:
            htmlfile.write("</div>\n")
            htmlfile.write("</body>\n")
            htmlfile.write("</html>\n")

    def fun_write_HTML(self, cluster_name,
                       sequence_name,
                       node_num,
                       monmer_length,
                       max_degree,
                       total_seq_in_cluster,
                       degree_list,
                       abundancy):
        # print(cluster_name,
        #                sequence_name,
        #                node_num,
        #                monmer_length,
        #                max_degree,
        #                total_seq_in_cluster,
        #                degree_list)
        htmlfile = self.out_html_file

        ### monomer length figure ###
        plt.hist(monmer_length)
        plt.title('Monomer length distribution in cluster')
        mon_length_figure = self.imgs_dir + '/Monomer_length_in_cluster{0}.png'.format(cluster_name)
        plt.ylabel('Number of sequences')
        plt.xlabel('Monomer length')
        plt.savefig(mon_length_figure, format='png')
        plt.close()
        #########################

        ### degree length figure ###
        plt.hist(degree_list)
        plt.title('Distributions of number of intra sequence connections in the cluster')
        degree_figure = self.imgs_dir + '/Degree_in_cluster{0}.png'.format(cluster_name)
        plt.ylabel('Number of sequences')
        plt.xlabel('Number of connections')
        plt.savefig(degree_figure, format='png')
        plt.close()
        #########################

        htmlfile.write("<div class='row bg-secondary'>\n")
        htmlfile.write('<div class="col-md-9"> <h1> Cluster {0} ({1}bp) </h1> </div>\n'.format(cluster_name, int(abundancy)))
        htmlfile.write("</div>\n")

        htmlfile.write("<div class='row bg-light'>\n")
        htmlfile.write('<div> <p>Selected sequence <b>{0}</b> :</p> </div>'
                       '<ul>'
                       '<li> Connected with {1} ({2}%) other sequences in the cluster </li>\n'.
                       format(sequence_name, max_degree - 1, round(((max_degree - 1) * 100) / total_seq_in_cluster, 1)))
        htmlfile.write('<li> Monomer_length: {} </li>\n'.format(monmer_length[node_num]))
        htmlfile.write('<li> Copy number: {} </li>\n'.format(int(abundancy/monmer_length[node_num])))
        htmlfile.write('</ul>')
        htmlfile.write("<div class='col-md-6'>\n")
        htmlfile.write('<img src = "' + mon_length_figure.replace(self.imgs_dir, "imgs") + '" alt = "monomer length distribution in cluster">\n')
        htmlfile.write("</div>\n")


        htmlfile.write("<div class='col-md-6'>\n")
        htmlfile.write('<img src = "' + degree_figure.replace(self.imgs_dir, "imgs")  + '" alt ="number of connections distribution in cluster">\n')
        htmlfile.write("</div>\n")
        htmlfile.write("</div>\n")
        htmlfile.write("<div class='row'>\n</div>\n")
