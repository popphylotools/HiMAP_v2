# HiMAP_v2

## Oksana's notes
How to run HiMAP v2:
1. Download HiMAP2 package (GitHub repo eventually)
2. unzip HiMAP2_master.zip
3. rename HiMAP2_master to HiMAP2
4. cd HiMAP2
5. create conda environment using bash script for linux or mac os
6. activate environment
   conda activate himap2
7. Set up directory system by running 
   ./bin/00_initilize_working_dir.py
8. Upload fasta files to "./data/fasta" directory
   Upload gff files to "./data/gff" directory
   Names of taxa should begin with the Annotation source name (e.g. NCBI, GeMoMa); these annotation sources should also be specified in the configuration file
   Sample input data files are provided in the "input" directory
8. Create gff data bases, extract CDS, select longest isoform, and prepare amino acid sequences for Orthofinder
   ./bin/01_sequence_extraction.py
9. Find orthogroups by running Orthofinder (input - amino acid sequences in the "./data/pep_fasta" directory; output - Orthogroups.tsv file in "./data/orthofinder" directory as well as all other output data file from Orthofinder run (./data/orthofinder/Orthofinder).
   ./bin/02_orthofinder.py
10. Select orthogroups to keep based on the specifications in the config file (single copy orthologs, number of copies for core, supplementary, and outgroup)(input - "Orthogroups.tsv" file in "./data/orthofinder"; output - "./data/orthofinder/keep_orthos.csv")
    ./bin/03_ortho_selection.py
11. Create and align padded exons, raw exons, and final filtered exons using settings specified in the configuration file (output - padded exons are in "./data/core_alignment"; raw exons are in "./data/supplementary_alignment"; filtered exons - "./data/ortho_cds_component"):
    ./bin/04_alignments_and_filtering.py
12. Final output from this part of the pipeline is in "./data/ortho_cds_component".


## Julian To-do notes
remove hyphy2 files
remove tapir/hyphy2 components of env building scripts
		ARE THERE OTHER MODULES IN THE MAIN ENV THAT ARE ONLY USED BY TAPIR? ugg, that might be apain to figure out...
    Also remove parts of environment files related to tapir
remove stuff for HiMAP_2.0_modulefile, since we're not using that now??

Create a "run_himap" script with logging (e.g. "step01 is complete"")

Rename the data dirs (ortho_cds_component)

multithreading for mafft steps more than what is done? Use gnu parallel?

Edit all .sh scripts to deal with envs overwriting each other's priority?

## Julian notes:
Conda envs for running on LCC cluster:
`#tweaked all .sh scripts to deal with envs overwriting each other
module load ccs/conda/python-3.7.3
conda env create -f environment.yml --prefix ./manual_envs3
conda install --prefix ./manual_envs3 -c bioconda orthofinder=2.5.4`

