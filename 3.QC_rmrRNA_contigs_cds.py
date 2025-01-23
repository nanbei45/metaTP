#!/usr/bin/python
__author__ = "Wang,Yansu"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Wang,yansu"
__email__ = "Wangys_c@hotmail.com"

import subprocess
import datetime
import os
import argparse
import re
from Bio import SeqIO
from os import system
##########################################
# def MakeOption():
#     # make option
#     parser = OptionParser(usage="%prog [-h] [-v] -i[--input=]-threads[--threads=]-phred[--phred=]-ILLUMINACLIP[--ILLUMINACLIP=]-LEADING[--LEADING=]-TRAILING[--TRAILING=]-SLIDINGWINDOW[--SLIDINGWINDOW=]-MINLEN[--MINLEN=]-o[--output]",
#                           version="%prog 1.2")
#     parser.add_option("-i", "--input", action="store", dest="input",
#                       help="the fastq folder",
#                       default=False)
#     parser.add_option("-m", "--m", action="store", dest="method",
#                       help="Paried-ended sequencing (PE) or Single-ended sequencing(SE)",
#                       default=False)
#     parser.add_option("--threads", "--threads", action="store", dest="threads",
#                       help="thread",
#                       default=16)
#     parser.add_option("-phred", "--phred", action="store", dest="phred",
#                       help="-phred33 or -phred64",
#                       default="-phred33")
#     parser.add_option("-ILLUMINACLIP", "--ILLUMINACLIP", action="store", dest="ILLUMINACLIP",
#                       help="paths for adapters",
#                       default="./trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10")
#     parser.add_option("-LEADING", "--LEADING", action="store", dest="LEADING",
#                       help="Remove front-end low-quality sequences or N bases,",
#                       default=5)
#     parser.add_option("-TRAILING", "--TRAILING", action="store", dest="TRAILING",
#                       help="Remove tail low-quality sequences or N bases,",
#                       default=5)
#     parser.add_option("-SLIDINGWINDOW", "--SLIDINGWINDOW", action="store", dest="SLIDINGWINDOW",
#                       help="SLIDINGWINDOW:4:15, scan sequences using a 4-base width sliding window and trim when the average mass per base is below 15",
#                       default="4:15")
#     parser.add_option("-MINLEN", "--MINLEN", action="store", dest="MINLEN",
#                       help="MINLEN:25, Remove the sequences that less than 25 bases",
#                       default=25)
#     parser.add_option("-o", "--output", action="store", dest="output",
#                       help="the output folder",
#                       default=False)
#     (options, args) = parser.parse_args()
#
#     input = options.input
#     method = options.method
#     threads = options.threads
#     phred = options.phred
#     ILLUMINACLIP = options.ILLUMINACLIP
#     LEADING = options.LEADING
#     TRAILING = options.TRAILING
#     SLIDINGWINDOW = options.SLIDINGWINDOW
#     MINLEN = options.MINLEN
#     output = options.output
#     return (input,threads,method,phred,ILLUMINACLIP,LEADING,TRAILING,SLIDINGWINDOW,MINLEN,output)

##################################################################################3
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help ='the fastq folder',required=True,type = str)
    parser.add_argument("-threads", "--threads", help="thread",type =int,default=24)
    parser.add_argument("-phred", "--phred", action="store_true", dest="phred",help="-phred33 or -phred64",default="-phred33")
    parser.add_argument("-ILLUMINACLIP", "--ILLUMINACLIP", action="store",help="paths for adapters",
                                            default="trimmomatic_adapters/TruSeq3")
    parser.add_argument("-LEADING", "--LEADING", dest="LEADING",
                      help="Remove front-end low-quality sequences or N bases,",type =int,default=5)
    parser.add_argument("-TRAILING", "--TRAILING", dest="TRAILING",
                      help="Remove tail low-quality sequences or N bases,",type =int,
                      default=5)
    parser.add_argument("-SLIDINGWINDOW", "--SLIDINGWINDOW", dest="SLIDINGWINDOW",
                      help="SLIDINGWINDOW:4:15, scan sequences using a 4-base width sliding window and trim when the average mass per base is below 15",
                      default="4:15")
    parser.add_argument("-MINLEN", "--MINLEN", dest="MINLEN",
                      help="MINLEN:25, Remove the sequences that less than 25 bases",type =int,
                      default=25)
    parser.add_argument("-MINIPRO", "--MINIPRO", dest="MINIPRO",
                        help="minimum protein length (default: 100)", type=int,
                        default=100)
    parser.add_argument("-o", "--output",dest="output",
                      help="the output folder",
                      default=False)
    args = parser.parse_args()
    input = args.input
    threads = args.threads
    phred = args.phred
    ILLUMINACLIP = args.ILLUMINACLIP
    LEADING = args.LEADING
    TRAILING = args.TRAILING
    SLIDINGWINDOW = args.SLIDINGWINDOW
    MINLEN = args.MINLEN
    MINIPRO = args.MINIPRO
    output = args.output
    return (input,threads,phred,ILLUMINACLIP,LEADING,TRAILING,SLIDINGWINDOW,MINLEN,MINIPRO,output)
#############################################################################################
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
#############################################################################################
def create_Folder(filename):
    file_path = os.getcwd() + '/' + filename
    if os.path.exists(file_path):
        print("The folder already exists.")
        print("exit")
        exit(1)
    else:
        os.mkdir(file_path)
    return (file_path)
################################################################################################
def QC_contol(inputfolder,outfolder,threads,phred,ILLUMINACLIP,LEADING,TRAILING,SLIDINGWINDOW,MINLEN):
    create_Folder(outfolder+"/QC_control")
    files = os.listdir(inputfolder)
    for file in files:
        if file.split(".")[1] == "fastq" or "fq":
            file1 = file
            if re.search ("_",str(file1)) == None:
                file2 = file.split(".")[0]
                cmd1 = "trimmomatic SE -threads " + str(threads) +" "+ str(phred) + " " + inputfolder +"/" + file2+ ".fastq.gz" +" " +outfolder+"/QC_control/"+file2+"_clean.fastq.gz "\
                   + "ILLUMINACLIP:"+ str(ILLUMINACLIP) +"-SE.fa:2:30:10"+ " LEADING:"+str(LEADING) + " TRAILING:" + str(TRAILING) + " SLIDINGWINDOW:" + str(SLIDINGWINDOW) +" MINLEN:" +str(MINLEN)
                run_command(cmd1)
            if re.search("_1",str(file1)) is not None:
                file2 = file.split("_")[0]
                cmd2 = "trimmomatic PE -threads " + str(threads) + " " + str(
                    phred) + " " + inputfolder + "/" + file2 + "_1.fastq.gz" + " " + inputfolder + "/" + file2 + "_2.fastq.gz" \
                       + " " + outfolder +"/"+"QC_control" + "/" + file2 + "_clean_forward_paired.fastq.gz " \
                       + " " + outfolder +"/"+ "QC_control" + "/" + file2 + "_clean_forward_unpaired.fastq.gz" \
                       + " " + outfolder +"/"+ "QC_control" + "/" + file2 + "_clean_reverse_paired.fastq.gz" \
                       + " " + outfolder +"/"+ "QC_control" + "/" + file2 + "_clean_reverse_unpaired.fastq.gz " \
                       + "ILLUMINACLIP:" + str(ILLUMINACLIP)+"-PE-2.fa:2:30:10" + " LEADING:" + str(LEADING) + " TRAILING:" + str(
                    TRAILING) + " SLIDINGWINDOW:" + str(SLIDINGWINDOW) + " MINLEN:" + str(MINLEN)
                run_command(cmd2)
        else:
            pass
#############################################################################################
def rmrRNA(outfolder):
    create_Folder(outfolder + "/rmrRNA")
    files = os.listdir(outfolder+"/QC_control")
    for file in files:
        if file.split("_")[1] == "clean.fastq.gz" :
            file2 = file.split("_")[0]
            cmd3 = "bowtie2 -x " + "./rRNA_databases/all_rRNA --un-gz " + outfolder+ "/rmrRNA/" + file2 +"_rmrRNA.fastq.gz" + " -U " + outfolder+"/QC_control/"+ file2 + "_clean.fastq.gz"+ " -p 16 -S "+ outfolder+"/rmrRNA/"+ file2 +"_rRNA.sam"
            run_command(cmd3)
            cmd_31 = "rm -rf " +outfolder +"/rmrRNA/"+file2 + "_rRNA.sam"
            run_command(cmd_31)
        if file.split("_",1)[1] == "clean_forward_paired.fastq.gz":
            file2 = file.split("_", 1)[0]
            cmd4 = "bowtie2 -x " + "./rRNA_databases/all_rRNA -1 " + outfolder+"/QC_control/" + file2 + "_clean_forward_paired.fastq.gz"+ " -2 " + outfolder+"/QC_control/"+ file2 + "_clean_reverse_paired.fastq.gz" + \
                   " --un-conc-gz " +  outfolder + "/rmrRNA/"+ file2 + "_rmrRNA.fastq.gz" +  " -p 16 -S " + outfolder + "/rmrRNA/" +file2 + "_rRNA.sam"
            run_command(cmd4)
            cmd_41 = "rm -rf " +outfolder +"/rmrRNA/"+file2 + "_rRNA.sam"
            run_command(cmd_41)
###########################################################################################
def to_contigs(outfolder):
    create_Folder(outfolder + "/megahit")
    path = os.path.abspath(outfolder + "/megahit")
    files = os.listdir(outfolder + "/rmrRNA")
    os.chdir(outfolder + "/rmrRNA")
    for file in files:
        if file.split("_")[1] == "rmrRNA.fastq.gz":
            file2 = file.split("_")[0]
            cmd5 = "megahit -r "+file2 + "_rmrRNA.fastq.gz -o " + file2 + "_megahit_out"
            run_command(cmd5)
            os.replace(file2 + "_megahit_out", path+"/"+file2 + "_megahit_out")
            #shutil.move(file2 + "_megahit_out", path)
        elif file.split("_")[1] == "rmrRNA.fastq.1.gz":
            file2 = file.split("_", 1)[0]
            cmd5 = "megahit -1 "+ file2 + "_rmrRNA.fastq.1.gz -2 " +file2 + "_rmrRNA.fastq.2.gz -o " + file2 + "_megahit_out"
            run_command(cmd5)
            os.replace(file2 + "_megahit_out", path+"/"+file2 + "_megahit_out")
            #shutil.move(file2 +"_megahit_out", path)
    os.chdir(os.path.split(os.path.realpath(__file__))[0])

#######################################################################################################################
def change_id(inputfasta,outputfasta,sample_id,outfolder):
    outputfasta = open(outfolder + "/megahit/"+outputfasta, "w")
    for each in SeqIO.parse(inputfasta, "fasta"):
        id = each.description
        seq = each.seq
        outputfasta.write(">%s\n%s\n" % (sample_id + "_" + id, seq))
    outputfasta.close()
    return outputfasta

def contigs_change_id(outfolder,MINIPRO):
    os.chdir(outfolder + "/rmrRNA")
    os.chdir(os.path.pardir)
    print(os.getcwd())
    os.chdir(os.path.pardir)
    print(os.getcwd())
    folders = os.listdir(outfolder + "/megahit")
    for folder in folders:
        if folder.split("_",1)[1] == "megahit_out":
            id = folder.split("_")[0]
            files = os.listdir(outfolder + "/megahit/" + folder)
            for file in files:
                if file == "final.contigs.fa":
                    change_id(outfolder + "/megahit/"+folder+"/"+file,id+"_"+file,id,outfolder)
                else:
                    pass
    print(">>>>>>>>>>>>>>>>>>>>>>>extracting coding sequence>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print(os.getcwd())
    os.chdir(outfolder + "/megahit")
    print(os.getcwd())
    folders = os.listdir("./")
    for file in folders:
        if file.split("_")[1] == "final.contigs.fa":
            cmd6 = "TransDecoder.LongOrfs -t " + file + " -m " + str(MINIPRO)
            run_command(cmd6)

##########################################################################################################
def change_id2(inputfasta,outputfasta):
    outputfasta = open(outputfasta, "w")
    for each in SeqIO.parse(inputfasta, "fasta"):
        id = each.description.split(" ")[0]
        seq = each.seq
        outputfasta.write(">%s\n%s\n" % (id, seq))
    outputfasta.close()
    return outputfasta

def rmdup(outfolder):
    print(">>>>>>>>>>>>>>>>>>>>>>>removing repetitive sequences>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    os.chdir(os.path.pardir)
    print(os.getcwd())
    os.chdir(os.path.pardir)
    print(os.getcwd())
    folders = os.listdir(outfolder + "/megahit")
    for folder in folders:
        if os.path.splitext(folder)[-1] == ".transdecoder_dir":
            id = folder.split("_")[0]
            files = os.listdir(outfolder + "/megahit/" + folder)
            for file in files:
                if file == "longest_orfs.cds":
                    new_file = id+"_"+file
                    os.rename(outfolder + "/megahit/" + folder+"/"+file,outfolder + "/megahit"+"/"+new_file)
                else:
                    pass
    os.chdir(outfolder + "/megahit")
    system("cat *longest_orfs.cds > all_longest_orfs_cds.fasta")
    system("seqkit rmdup -s all_longest_orfs_cds.fasta -o all_longest_orfs_cds_rmdup.fasta")
    inputfasta =  "all_longest_orfs_cds_rmdup.fasta"
    outputfasta = "all_longest_orfs_cds_rmdup_id.fasta"
    change_id2(inputfasta,outputfasta)

#######################################################################################################
def main():
    input, threads, phred, ILLUMINACLIP, LEADING, TRAILING, SLIDINGWINDOW, MINLEN, MINIPRO,output = parse_args()
    QC_contol(input,output,threads,phred,ILLUMINACLIP,LEADING,TRAILING,SLIDINGWINDOW,MINLEN)
    rmrRNA(output)
    to_contigs(output)
    contigs_change_id(output,MINIPRO)
    rmdup(output)

if __name__ == "__main__":
    main()
