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
* **`./sam`** includes sam files
* **`./bed`** includes bed files and genome_size files
* **`./temp`** it includes the temporary files, only if you set the save_temp as TRUE
* **`./search_stat.tsv`** tab-delimited stat file. It inclues how much length it matched and which contigs/chromosome the gene was found
* **`./best_result.tsv`** tab-delimited result file. It inclues the best mec gene, mec complex class, and ccr combination with the percentage of match.

# Introduction
According to the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/), there are two factores that decides the sccmec complex type: mec complex class and ccr genes. 

## Step1. preparing gene database(DB)
Here the genes are collected from [paper, Table2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772726/) and saved as multi-fasta format in gene_db.fasta file.

| sccmec type   | Representative Strain <br> (Genebank Accession) |
| ------------- |:---------------------:|
|I              |[NCTC10442(AB033763)](https://www.ncbi.nlm.nih.gov/nuccore/AB033763)  |
|II             |[N315 (D86934)](https://www.ncbi.nlm.nih.gov/nuccore/D86934)  |
|III            |[85/2082 (AB037671)](https://www.ncbi.nlm.nih.gov/nuccore/AB037671)  |
|IV             |[CA05 (AB063172)](https://www.ncbi.nlm.nih.gov/nuccore/AB063172)  |
|V              |[WIS (AB121219)](https://www.ncbi.nlm.nih.gov/nuccore/AB121219)  |
|VI             |[HDE288 (AF411935)](https://www.ncbi.nlm.nih.gov/nuccore/AF411935)  |
|VII            |[P5747/2002(AB373032)](https://www.ncbi.nlm.nih.gov/nuccore/AB373032)  |
|VIII           |[C10682 (FJ390057)](https://www.ncbi.nlm.nih.gov/nuccore/FJ390057)  |
|IX             |[JCSC6943 (AB505628)](https://www.ncbi.nlm.nih.gov/nuccore/B505628)  |
|X              |[JCSC6945 (AB505630)](https://www.ncbi.nlm.nih.gov/nuccore/AB505630)  |
|XI             |[LGA251(FR821779, WGS)](https://www.ncbi.nlm.nih.gov/nuccore/FR821779)  |
|XII            |[BA01611 (KR187111)](https://www.ncbi.nlm.nih.gov/nuccore/KR187111)  |
|XIII           |[55-99-44 (MG674089)](https://www.ncbi.nlm.nih.gov/nuccore/MG674089) |
|XIV            |[SC792 (LC440647)](https://www.ncbi.nlm.nih.gov/nuccore/LC440647)  |

## Step2. Run the multiple sequence alignment (minimap2)
By using minimap2, which is faster than BLAST but less accurate, it finds the matching contigs and length. Results will be saved in `./sam` folder. 

## Step3. Filter the results
There are two potential cases of incorect matching result.
1. If the assembly file are made of many contigs and if the matching is ended at the first of sequence or end of the sequence, alignment cannot find the proper matching length.
2. Same genes are detected in similar position multiple times

To fix these errors, first, I filtered out the duplicated results. After that, program checks if there is any sequences ended at the start of end of the contigs. If there is, they cut out the matched part of the gene and re-run the code. It runs at maximum three times. The final result will be saved in `./sam/{reference_name}_combined.sam`

## Step4. Create the stat file.
By reading the final `./sam/{reference_name}_combined.sam`. it will generate the stat file. Stat file inclues how much length it matches out of total length and which contigs it is located. All results will be saved in tab delimited format in file name `serach_stat.tsv`.

## Step5. Find the best match
By reading the `serach_stat.tsv`, it will choose the best mec gene, mec class complex, ccr complex gene. If there is many best matches within 10% differences, it prioritise the one shares the same contigs with mec gene. If the best ccr gene match is ccrA or ccrB it will combine the two result into one and print out the average matching length. All results will be saved in tab delimited format in file name `best_result.tsv`.

## Step6. check the error
This is to-do list. I did it by manually so far. Prepare a cup of Coffee.


