# sccmec_classifier
goal: classifying sccmec type automatically

code I used for sccmec_typing. There are two methods

# First approach : gene approach
According to the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/). There are three type of mec class: mec complex A, B and C.

'run_distinguish_type.sh' code distinguishes the mec complex type by the order order of IS sequences and mecA/RI/I.

# Second approach : direct approach
From [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/), the downloaded mec classes sequences from reference genomes are saved in 'mec_class.fasta'

'run_distinguish_type_direct.sh' code distinguishes the mec cmoplex type by this direct approach.

## ccr type
ccr types are decided by the best match with reference genes. the match percentages are present in the result file.