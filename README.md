# sccmec_classifier
goal: classifying sccmec type automatically

code I used for sccmec_typing. There are two methods

*warning: this code is currently only runable under ARC server (#todo-list: make it runable generally)

# First approach : gene approach
According to the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/). There are three type of mec class: mec complex A, B and C.

`run_distinguish_type.sh` code distinguishes the mec complex type by the order order of IS sequences and mecA/RI/I.

# Second approach : direct approach
From [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/), the downloaded mec classes sequences from reference genomes are saved in 'mec_class.fasta'

`run_distinguish_type_direct.sh` code distinguishes the mec cmoplex type by this direct approach.

## ccr type
ccr types are decided by the best match with reference genes. the match percentages are present in the result file. I gave 1 point for match 0.5 for soft-cliping and -2 for in/del

# output
* *`./sam`* includes sam files
* *`./bed`* includes bed files and genome_size files
* `ccr_info.tsv` has the information of ccr genes and their matching percentage
* `mec_info.tsv` has the information of mec complex
* `combined_results.tsv` name-based merged file of `ccr_info.tsv` and `mec_info.tsv`.