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
################################################################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help ='the fastq folder',required=True,type = str)
    parser.add_argument("-k", "--kmer", help="kmer", type=int, default=31)
    parser.add_argument("-p", "--threads", help="thread",type =int,default=24)
    parser.add_argument("-o", "--output",dest="output",
                      help="the output folder",
                      default=False)
    args = parser.parse_args()
    input = args.input
    kmer = args.kmer
    threads = args.threads
    output = args.output
    return (input,kmer,threads,output)
###################################################################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
######################################################################################################
def create_Folder(filename):
    file_path = os.getcwd() + '/' + filename
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
        os.mkdir(file_path)
    else:
        os.mkdir(file_path)
    return (file_path)
##########################################################################################################
def salmon_index_create(input,kmer,threads,output):
    create_Folder(output +"salmon_index")
    cmd3 = "salmon index -t " + input + " -i " + output + "/transcripts_index" + " -k " + str(kmer) + " -p " + str(threads)
    run_command(cmd3)
#############################################################################################################
def main():
    input,kmer,threads,output = parse_args()
    salmon_index_create(input,kmer,threads,output)
if __name__ == "__main__":
    main()