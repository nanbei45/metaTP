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
                      help="the fastq folder",
                      default=False)
    parser.add_option("-o", "--output", action="store", dest="output",
                      help="the QC folder",
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
        os.mkdir(file_path)
    else:
        os.mkdir(file_path)
    return (file_path)
#################################################################
def QC_test(inputfolder,outfolder):
    create_Folder(outfolder)
    files = os.listdir(inputfolder)
    for file in files:
        if file.split(".")[1] == "fastq":
            if file.split(".",2)[1] == ".fastq.gz":
                file1 = file.split(".")[0]
                cmd2 = "fastqc -t 24 -o " + outfolder + ' --nogroup ' + inputfolder + "/" + file1 + ".fastq.gz"
                run_command(cmd2)
            else:
                file1 = file.split("_")[0]
                cmd1 = "fastqc -t 24 -o " + outfolder + ' --nogroup ' + inputfolder +"/"+file1 +"_1.fastq.gz" \
                       +" " +inputfolder +"/"+file1 +"_2.fastq.gz"
                run_command(cmd1)

    return()
##############################################################
def main():
    input,output = MakeOption()
    QC_test(input,output)
if __name__ == "__main__":
    main()