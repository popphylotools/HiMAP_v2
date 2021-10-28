#!/usr/bin/env python
# coding=utf-8
""""""
#align mafft
#raxml model
#raxml-ng

import os,argparse,sys
import logging
import configparser
import multiprocessing as mp
import itertools
import subprocess
from Bio import SeqIO


#parser=MyParser()
parser = argparse.ArgumentParser()
parser.add_argument('--configPath', help='configPath', default='./config.ini')
args = parser.parse_args()
# load parameters from config file
config = configparser.ConfigParser()
config.read(args.configPath)

input_dir = config["Paths"]["final_exon_dir"]
output_dir = config["Paths"]["exon_phylo_dir"]
num_cores = mp.cpu_count()

if input_dir[-1] != "/":
		input_dir += "/"
if output_dir[-1] != "/":
		output_dir += "/"

"""
Functions
"""

#Multitask function
def multitask_mafft_raxml_ng(input_list,output_list,clusterID_list,number_processors=mp.cpu_count()):
    # create pool
    with mp.Pool(min(len(input_list), number_processors)) as p:
        p.starmap(mafft_raxml_ng_wraper, zip(input_list,output_list,clusterID_list))

#Create an input using the first word as header
def fasta_one_word_header(infile,outfile):
	
	with open(outfile,"w") as out:
		for record in SeqIO.parse(infile, "fasta"):
			id = record.id.split()[0]
			seq = str(record.seq)
			out.write(">" + id + "\n" +
				      seq + "\n")
	return outfile

#Create a directory
def create_dir(dir_path):
	if not os.path.exists(dir_path):
		os.makedirs(dir_path)
	return dir_path

#Align using mafft and activating the flag localpair option
def align_mafft(infile, outfile):
	output_dir = os.path.dirname(outfile)
	clusterID = os.path.basename(output_dir)

	if not os.path.exists(outfile):
		#Run the analysis
		with open(outfile,"wb") as out_handle:
			cmd = ["mafft", "--localpair","--maxiterate","1000","--thread","1","--quiet", infile]
			p = subprocess.Popen(cmd, stdout=out_handle)
			out = p.wait()
		#logging.info(clusterID + ": mafft alignment has finished the analysis.")
		return outfile
	
	else:
		#logging.info(clusterID + ": mafft alignment had previously been run.")
		print(infile)
		print(outfile)
		return outfile


#Get the best fit model
def modeltest_ng(infile,selector_criterion="BIC"):
	#Test the evolutionary model 0.
	#Creating the directory structure
	#output_cluster_dir = create_dir(output_cluster_dir)

	output_dir = os.path.dirname(infile)
	clusterID = os.path.basename(output_dir)
	output_prefix = os.path.splitext(infile)[0]
	outlog = output_prefix + ".log"

	if not os.path.exists(outlog):
		#Run the analysis
		cmd = ["modeltest-ng", "-i",infile,"-p","1","T","raxml","-o",output_prefix]
		p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		out = p.wait()
		#logging.info(clusterID + ": mafft alignment has finished the analysis.")
	
	with open(outlog,"r") as out:
		lines = list(out.readlines())[-6:-3]
		for line in lines:
			list_lines = line.strip().split()
			if list_lines[0] == selector_criterion:
				model = list_lines[1]
				break
	
	#Cleaning steps
	#
	#
	return model

#phylogeny inference using raxml-ng
def raxml_ng(infile,model,bootstrap="200"):
	output_dir = os.path.dirname(infile)
	clusterID = os.path.basename(output_dir)
	outfile = os.path.dirname(infile) + ".raxml.support"

	if not os.path.exists(outfile):
		#Run the analysis
		cmd = ["raxml-ng","seed","1","--all","--threads","1","--bs-trees",bootstrap,"--msa",infile,"--model",model]
		p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		out = p.wait()
		#logging.info(clusterID + ": mafft alignment has finished the analysis.")
		return True
	
	else:
		#logging.info(clusterID + ": mafft alignment had previously been run.")
		return True

#Create the lists with the input, output files and clusterIDs.
def prepare_multitask(input_dir,output_dir):
	input_list = []
	output_list = []
	clusterID_list = []
	for filename in os.listdir(input_dir):
		clusterID = os.path.splitext(filename)[0]
		extension = os.path.splitext(filename)[1]
		if extension == ".fasta":
			input_fasta = input_dir + filename
			output_cluster_dir = output_dir + clusterID + "/"
			outfile = output_cluster_dir + filename
			input_list.append(input_fasta)
			output_list.append(outfile)
			clusterID_list.append(clusterID)
	return input_list,output_list,clusterID_list

#Join all functions
def mafft_raxml_ng_wraper(input_fasta,outfile,clusterID):
	#Make the directory
	output_cluster_dir = os.path.dirname(outfile)
	create_dir(output_cluster_dir)
	input_fasta_one_word_header = fasta_one_word_header(input_fasta,outfile)
	#Align the sequences
	outfile_aln = output_cluster_dir + "/" + clusterID + "_aln.fasta"
	aln_file = align_mafft(input_fasta_one_word_header,outfile_aln)
	#Get the best fit model
	model = modeltest_ng(aln_file) # Change the default criterion "BIC"
	#Get the exon phylogeny
	raxml_ng(aln_file,model) # Change the default bootstrap replicates "200"
	
	return True


"""
Main
"""
input_list,output_list,clusterID_list = prepare_multitask(input_dir,output_dir)
multitask_mafft_raxml_ng(input_list,output_list,clusterID_list,num_cores) 
print("Step_05 is complete")
