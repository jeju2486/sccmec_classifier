# sccmec_classifier
goal: classifying sccmec type automatically

code I used for sccmec_typing. It uses minimap2 as the main alignment tool which is faster but less accurate than BLAST. 

*warning*: I did not test if this runnable in other computer yet.
*warning*: This is only currently runnable in Linux.

(##todo-list: make it as conda package)
(##todo-list: add subtyping result)
(##todo-list: add mec class D)

# Installation
First you need to clone the github 

```ruby
git clone https://github.com/jeju2486/sccmec_classifier.git
```

Assuming that you have Anaconda3 for environmental setup. Run the code below to build the environment.

```ruby
python setup.py
```
Now you are ready run the code. Here is the example way to run the code. you should change the each parameter to your file/directory address.

```ruby
bash run_script.sh -q gene_db.fasta -r "$reference_dir" -o "$output_dir" -s "$save_temp"
```

parameters are
* *-q* : query_sequence file. Here the code uses gene_db.fasta as the default. But you can change it to your sequences.
* *-r* : directory address to assemblies or reference genomes.
* *-o* : output directory
* *-s* : if you want to save the intermediate temporary files. "TRUE" as default.

# output
* *`./sam`* includes sam files
* *`./bed`* includes bed files and genome_size files
* *`./temp`* it includes the temporary files, only if you set the save_temp as TRUE

# Introduction
According to the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/), there are two factores that decides the sccmec complex type: mec complex class and ccr genes

`run_distinguish_type.sh` code distinguishes the mec complex type by the order order of IS sequences and mecA/RI/I.

## Second approach : direct approach
From [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/), the downloaded mec classes sequences from reference genomes are saved in 'mec_class.fasta'

`run_distinguish_type_direct.sh` code distinguishes the mec cmoplex type by this direct approach.

## ccr type
ccr types are decided by the best match with reference genes. the match percentages are present in the result file. I gave 1 point for match 0.5 for soft-cliping and -2 for in/del

