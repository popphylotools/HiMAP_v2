# HiMap_2.0

## Oksana's directions
1. Download HiMAP2 package (GitHub repo eventually)
2. unzip HiMAP2_master.zip
3. rename HiMAP2_master to HiMAP2
4. cd HiMAP2
5. create conda environment using bash script for linux or mac os
	 We don't really need this anymore if we add orthofinder (with version) to the environment.yml file, as that'll load mafft
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

## Julian to-do list:
remove hyphy2 files
remove tapir/hyphy2 components of env building scripts
		ARE THERE OTHER MODULES IN THE MAIN ENV THAT ARE ONLY USED BY TAPIR? ugg, that might be apain to figure out...
		Figure out env issues 
remove stuff for HiMAP_2.0_modulefile, since we're not using that now??

Create a "run_himap" script with logging (e.g. "step01 is complete"")

Rename the data dirs (ortho_cds_component)

multithreading for mafft-like stuff



## Julian's notes
conda env for LCC:
```
#tweaked all .sh scripts to deal with envs overwriting each other
module load ccs/conda/python-3.7.3
conda env create -f environment.yml --prefix ./manual_envs3
conda install --prefix ./manual_envs3 -c bioconda orthofinder=2.5.4
```


Carlos' env:
```
name: HiMAP2
channels:
  - anaconda
  - defaults
  - conda-forge
  - bioconda
dependencies:
  - argcomplete=1.12.3=pyhd3eb1b0_0
  - argh=0.26.2=py38_0
  - backcall=0.2.0=py_0
  - blas=1.0=mkl
  - ca-certificates=2020.10.14=0
  - certifi=2020.6.20=py38_0
  - configparser=5.0.1=py_0
  - decorator=4.4.2=py_0
  - gffutils=0.10.1=pyh864c0ab_1
  - importlib-metadata=4.8.1=py38h06a4308_0
  - importlib_metadata=4.8.1=hd3eb1b0_0
  - intel-openmp=2021.3.0=h06a4308_3350
  - ipython=7.18.1=py38h5ca1d4c_0
  - ipython_genutils=0.2.0=py38_0
  - jedi=0.18.0=py38h06a4308_1
  - ld_impl_linux-64=2.33.1=h53a641e_7
  - libedit=3.1.20191231=h14c3975_1
  - libffi=3.3=he6710b0_2
  - libgcc-ng=9.1.0=hdf63c60_0
  - libstdcxx-ng=9.1.0=hdf63c60_0
  - mkl=2021.3.0=h06a4308_520
  - mkl-service=2.4.0=py38h7f8727e_0
  - mkl_fft=1.3.0=py38h42c9631_2
  - mkl_random=1.2.2=py38h51133e4_0
  - ncurses=6.2=he6710b0_1
  - numpy=1.21.2=py38h20f2e39_0
  - numpy-base=1.21.2=py38h79a1101_0
  - openssl=1.1.1h=h7b6447c_0
  - pandas=1.1.3=py38he6710b0_0
  - parso=0.8.0=py_0
  - pexpect=4.8.0=py38_0
  - pickleshare=0.7.5=py38_1000
  - pip=20.2.4=py38_0
  - prompt-toolkit=3.0.8=py_0
  - ptyprocess=0.6.0=py38_0
  - pyfaidx=0.6.2=pyh5e36f6f_0
  - pygments=2.7.1=py_0
  - python=3.8.5=h7579374_1
  - python-dateutil=2.8.1=py_0
  - pytz=2020.1=py_0
  - readline=8.0=h7b6447c_0
  - setuptools=50.3.0=py38hb0f4dca_1
  - simplejson=3.17.3=py38h7f8727e_2
  - six=1.16.0=pyhd3eb1b0_0
  - sqlite=3.33.0=h62c20be_0
  - tk=8.6.10=hbc83047_0
  - traitlets=5.0.5=py_0
  - wcwidth=0.2.5=py_0
  - wheel=0.35.1=py_0
  - xz=5.2.5=h7b6447c_0
  - zipp=3.6.0=pyhd3eb1b0_0
  - zlib=1.2.11=h7b6447c_3
  ```
