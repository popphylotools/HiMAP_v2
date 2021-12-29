#!/usr/bin/env python
import os
import configparser
import subprocess
config = configparser.ConfigParser()
from context import project_dir
config.read(os.path.join(project_dir, 'config.ini'))

"""
FUNCTIONS
"""

#Concatenate the exon trees in a new file
def get_trees(exon_phylo_dir,speciestree_dir):
	path_trees = os.path.join(speciestree_dir,"exon_trees_support.tree") 
	
	with open(path_trees, 'w') as output_newick:
		for exon_ID in os.listdir(exon_phylo_dir):
			tree_path = os.path.join(exon_phylo_dir,exon_ID,exon_ID + "_aln.fasta.raxml.support")
			if os.path.isfile(tree_path):
				with open(tree_path, 'r') as input_newick:
					output_newick.write(input_newick.readline())

	return path_trees

def run_astral(speciestree_dir,astralpath,path_trees,max_mem):
	output_speciestree = os.path.join(speciestree_dir,"speciestree.tre")
	log_err_file = os.path.join(speciestree_dir,"speciestree.log")
	java_max_mem =  "-Xmx" + str(max_mem) + "M"

	#Get the species_tree
	with open(log_err_file, "wb") as err_log:
		cmd = ["java",java_max_mem,"-jar",astralpath,"-i",path_trees,"-o",output_speciestree] 
		p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err_log)
		out = p.communicate()
	
	#Annotate the species tree (quartet support)
	output_speciestree_qs = os.path.join(speciestree_dir,"speciestree_quartetsupport.tre")
	log_err_file_qs = os.path.join(speciestree_dir,"speciestree_quartetsupport.log")

	with open(log_err_file_qs, "wb") as err_log:
		cmd = ["java",java_max_mem,"-jar",astralpath,"-q",output_speciestree,"-t","1","-i",path_trees,"-o",output_speciestree_qs] 
		p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err_log)
		out = p.communicate()

	return

"""
MAIN
"""

#Get the paths from the config file
speciestree_dir = config["Paths"]["speciestree_dir"]
exon_phylo_dir = config["Paths"]["exon_phylo_dir"]
astralpath = config["External_program"]["astral_path"]
max_mem = config["External_program"]["max_mem"]

#Concatenate exon trees
path_trees = get_trees(exon_phylo_dir,speciestree_dir)

#Run ASTRAL
run_astral(speciestree_dir,astralpath,path_trees,max_mem)
