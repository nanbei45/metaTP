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
from optparse import OptionParser

##########################################
def MakeOption():
    # make option
    parser = OptionParser(usage="%prog [-h] [-v] -i[--input=]-o[--output]",
                          version="%prog 1.2")
    parser.add_option("-i", "--input", action="store", dest="input",
                      help="the SRR_Acc_List.txt",
                      default=False)
    parser.add_option("-o", "--output", action="store", dest="output",
                      help="the fastq folder ",
                      default=False)
    (options, args) = parser.parse_args()

    # extract option from command line
    input = options.input
    output = options.output
    return (input,output)
#############################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
####################################################################
def create_Folder(filename):
    file_path = os.getcwd() + '/' + filename
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
        #os.removedirs(file_path)
        os.mkdir(file_path)
    else:
        os.mkdir(file_path)
    return (file_path)
#########################################################################3
def prefetch_sra(inputfile,outfile):
    create_Folder(outfile)
    with open(inputfile, "r", encoding='utf-8') as f:  # 打开文本
        for sra_acc in f.readlines():
            sra_acc = sra_acc.strip('\n')
            cmd1 = "prefetch " + sra_acc + ' --max-size 200G --output-directory ' + outfile
            run_command(cmd1)
    return()
###################################################################################
def sra2fastq(inputfile,outfile):
    create_Folder('./' + outfile + '/fastq')
    with open(inputfile, "r", encoding='utf-8') as f:  # 打开文本
        for sra_acc in f.readlines():
            sra_acc = sra_acc.strip('\n')
            os.system("parallel-fastq-dump --sra-id " + sra_acc + " --threads 24 --outdir " + './' + outfile + '/fastq' + " --split-3 --gzip")
###################################################################################
def main():
    input,output = MakeOption()
    prefetch_sra(input,output)
    sra2fastq(input,output)
if __name__ == "__main__":
    main()



