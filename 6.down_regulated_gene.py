from optparse import OptionParser
import os
import subprocess
import datetime
from Bio import SeqIO
import pandas as pd

def Options():
    parser = OptionParser(usage="RMT_searching.py [-h, -v] -i[--input_file = ] -L[--low_threshold = ] -H[--high_threshold = ] -o[--output_file = ]\
", version="%prog 1.3")
    parser.add_option("-i", "--input_file", action="store", dest="inputFile",
                      help="The DEG_result containing differential_genes.csv.",
                      default=False)
    parser.add_option("-s", "--seq", action="store", dest="seq",
                      help="the sequence for transcript_index;",
                      default=False)
    # ///
    (options, args) = parser.parse_args()
    input = options.inputFile
    seq = options.seq
    return (input,seq)

###################################################################################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
###################################################################################################
def seq_extract(seq,input):
    print(">>>>>>>>>> Extracting differential gene sequences...................")
    outputfile = "differential_gene_sequence_down.fasta"
    data = pd.read_csv(os.path.abspath(input) + "/" + "differential_genes.csv")
    df1= data.set_index(['GeneID'])['change'].to_dict()
    df2 = {k: v for k, v in df1.items() if v =='Down'}
    #df2 = {k: v for k, v in df1.items() if v == 'Stable'}
    for key in df2.keys():
        with open(input +'/differential_genes_id_down.txt', 'a') as f:
            f.write(key)
            f.write('\r\n')
    path = os.path.abspath(input + "/" + "differential_genes_id_down.txt")
    cmd2 = "seqkit grep -f "+ path +' ' +seq + " -o " + input + "/"+ outputfile
    run_command(cmd2)
    path2 = os.path.abspath(input + "/" + "differential_gene_sequence_down.fasta")
    path3 = os.path.abspath(input + "/" + "differential_gene_sequence_down.pep")
    cmd3 = "seqkit translate " + path2 + " > " + path3
    run_command(cmd3)
    print(">>>>>>>>>> Differential down-regulated gene sequences extraction finished...................")
####################################################################################################################
def main():
    input,seq = Options()
    seq_extract(seq,input)

if __name__ == "__main__":
    main()