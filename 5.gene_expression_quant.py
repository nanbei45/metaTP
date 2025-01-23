#!/usr/bin/python
__author__ = "Wang,Yansu"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Wang,yansu"
__email__ = "Wangys_c@hotmail.com"

import subprocess
import datetime
import os
import shutil
import argparse
import re
import pandas as pd
import csv
################################################################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help ='the fastq folder',required=True,type = str)
    parser.add_argument("-index", "--index", help="the index folder", required=True,type = str)
    parser.add_argument("-p", "--threads", help="thread",type =int,default=24)
    parser.add_argument("-o", "--output",dest="output",
                      help="the output folder",
                      default=False)
    args = parser.parse_args()
    input = args.input
    index = args.index
    threads = args.threads
    output = args.output
    return (input,index,threads,output)
###################################################################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
############################################################################################################################
def create_Folder(filename):
    file_path = os.getcwd() + '/' + filename
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
        os.mkdir(file_path)
    else:
        os.mkdir(file_path)
    return (file_path)
############################################################################################################################
def salmon_quant(input,index,threads,output):
    create_Folder(output + "/transcripts_quant")
    files = os.listdir(input)
    for file in files:
        print(file.split(".",1)[1])
        if file.split(".")[1] == "fastq.gz":
            file2 = file.split(".")[0]
            cmd1 = "salmon quant -i " + index + " -l SF -r " + input+"/"+file2 + "_rmrRNA.fastq.gz" + " -p " + str(threads) + " -o " + output + "/transcripts_quant/" + file2+ "_quant_result"
            run_command(cmd1)
        if file.split(".",1)[1] == "fastq.1.gz":
            file2 = file.split("_")[0]
            cmd2 = "salmon quant -i " + index + " -l IU -1 " + input+ "/"+ file2 + "_rmrRNA.fastq.1.gz" + " -2 " + input+ "/"+ file2 +"_rmrRNA.fastq.2.gz" +" -p " + str(threads) + " -o " + output + "/transcripts_quant/" + file2+ "_quant_result"
            run_command(cmd2)
        else:
            pass
##############################################################################################################################
def sort_key(s):
    #sort_strings_with_embedded_numbers
    re_digits = re.compile(r'(\d+)')
    pieces = re_digits.split(s)  # 切成数字与非数字
    pieces[1::2] = map(int, pieces[1::2])  # 将数字部分转成整数
    return pieces
def merge_tables (output):
    files = os.listdir(output + "/transcripts_quant")
    for file in files:
       if file.endswith("_quant_result") is True:
           sample_name = file.split("_", 1)[0]
           dict1 = {}
           dict1["GeneID"] = sample_name
           print(dict1)
           for file1 in os.listdir(output+"/transcripts_quant" + "/"+file):
               if file1.endswith(".sf") ==True:
                   file2 = open(output+"/transcripts_quant"+"/"+file+"/"+file1, 'r')
                   lines = file2.readlines()[2:]
                   for line in lines:
                       a = line.split()
                       if float(a[3]) >0:
                           dict1[a[0]] = a[3]
                   file2.close()
               output_file_name = output+"/transcripts_quant"+"/"+sample_name + "_gene_matrix.csv"
               with open (output_file_name,"w",newline ='')as f1:
                   writer =csv.writer(f1)
                   for row in dict1.items():
                       writer.writerow(row)
    gene_matrix_filelist = []
    for file in os.listdir(output+"/transcripts_quant"):
        if file.split("_",1)[1] == "gene_matrix.csv":
            gene_matrix_filelist.append(file)
    counts = [pd.read_csv(output+"/transcripts_quant"+"/"+f, index_col=0, header=0, skiprows=0)for f in gene_matrix_filelist]
    matrix = pd.concat(counts, axis=1).fillna(0)
    #排序
    print(matrix.index.to_series().str.rsplit("_").str[-1].str.extract(pat=r'(\w+)|m.'))
    matrix2 = matrix.sort_index().reindex(
        matrix.index.to_series().str.rsplit("_").str[-1].str.extract(pat=r'(\w+)|m.')[0].astype(float).sort_values().index
    )
    matrix2.to_csv(output+"/transcripts_quant"+"/"+"transcript_abundance_quantification_table.csv", sep=",")

#########################################################################################################################

def main():
    input,index,threads,output = parse_args()
    salmon_quant(input,index,threads,output)
    merge_tables (output)

if __name__ == "__main__":
    main()