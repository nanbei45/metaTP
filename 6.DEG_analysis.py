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
                      help="The correlation coefficients matrix; CSV-based format is necessary.",
                      default=False)
    parser.add_option("-g", "--group_informaton", action="store", dest="group_informaton",
                      help="The group information; csv-based file format is necessary [required].",
                      default=False)
    parser.add_option("-n", "--group_name", action="store", dest="group_name",
                      help="Experimental group/control group;which is the same as that in the groupfile;such as Rhizosphere/Bulk[required].",
                      default=False)
    parser.add_option("-p", "--cut_off_pvalue", action="store", type = "float", dest="cut_off_pvalue",
                      help="the value of cut_off_pvalue",
                      default=0.001)
    parser.add_option("-f", "--cut_off_logFC", action="store", dest="cut_off_logFC",type="float",
                      help="the value of cut_off_logFC;",
                      default=1)
    parser.add_option("-s", "--seq", action="store", dest="seq",
                      help="the sequence for transcript_index;",
                      default=False)
    parser.add_option("-o", "--output_folder", action="store", dest="outputFile",
                      help="The output result folder",
                      default=False)
    # ///
    (options, args) = parser.parse_args()
    input = options.inputFile
    group = options.group_informaton
    group_name = options.group_name
    seq = options.seq
    cut_off_pvalue = options.cut_off_pvalue
    cut_off_logFC = options.cut_off_logFC
    output = options.outputFile
    return (input, group, group_name,seq,cut_off_pvalue,cut_off_logFC,output)

###################################################################################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
###################################################################################################
def DEG(input,group,group_name,cut_off_pvalue,cut_off_logFC,output):
    ##Rscript DEG_analysis.r -i bacteria_asv_occ15_table.csv -g bacteria_group4.csv -n Rhizosphere/Bulk -o DEG_result
    ##Rscript DEG_analysis.r -i ./test_sra_data/transcripts_quant/transcript_abundance_quantification_table.csv -g ./test_sra_data/transcripts_quant/sample_group.csv -n rhizosphere/bulk -o DEG_result
    cmd1 = "Rscript DEG_analysis.r " + "-i " + input+ " -g " + group + " -n " +group_name + " -o " +output +" -p " + str(cut_off_pvalue) + " -f " + str(cut_off_logFC)
    run_command(cmd1)
##################################################################################################################
def seq_extract(seq,output):
    print(">>>>>>>>>> Extracting differential gene sequences...................")
    outputfile = "differential_gene_sequence.fasta"
    data = pd.read_csv(os.path.abspath(output) + "/" + "differential_genes.csv")
    df1= data.set_index(['GeneID'])['change'].to_dict()
    df2 = {k: v for k, v in df1.items() if v =='Up' or v == 'Down'}
    #df2 = {k: v for k, v in df1.items() if v == 'Stable'}
    for key in df2.keys():
        with open(output +'/differential_genes_id.txt', 'a') as f:
            f.write(key)
            f.write('\r\n')
    path1 = os.path.abspath(output + "/" + "differential_genes_id.txt")
    cmd2 = "seqkit grep -f "+ path1 +' ' +seq + " -o " + output + "/"+ outputfile
    run_command(cmd2)
    path2 = os.path.abspath(output + "/" + "differential_gene_sequence.fasta")
    path3 = os.path.abspath(output + "/" + "differential_gene_sequence.pep")
    cmd3 = "seqkit translate " + path2 + " > " + path3
    run_command(cmd3)
    print(">>>>>>>>>> Differential gene sequences extraction finished...................")
####################################################################################################################
def main():
    input, group, group_name,seq,cut_off_pvalue,cut_off_logFC,output = Options()
    DEG(input,group,group_name,cut_off_pvalue,cut_off_logFC,output)
    seq_extract(seq, output)

if __name__ == "__main__":
    main()