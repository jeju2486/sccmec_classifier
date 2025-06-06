# SCCmecClassifier

Automated **SCCmec** typing for *Staphylococcus* assemblies

<!-- badges are welcome here -->

---

## Why SCCmecClassifier?

* **Fast.** Uses `minimap2` for rapid alignments.
* **Accurate.** Post‑processing filters correct most edge‑cases caused by split contigs and duplicated hits.
* **Reproducible.** Single command, deterministic outputs.
* **Extensible.** Drop‑in databases make adding new mec/ccr variants trivial.

> **Note** The workflow is currently tested on Linux. macOS should work with minimal tweaks; Windows support is not planned.

---

## Installation

```bash
# 1. clone
$ git clone https://github.com/jeju2486/sccmec_classifier.git
$ cd sccmec_classifier

# 2. create & activate environment
$ conda env create -f environment.yml
$ conda activate sccmec_classifier

# 3. (optional) legacy install
$ python setup.py install
```

---

## Quick‑start

```bash
bash run_script.sh \
  -q path/to/gene_db.fasta \
  -r path/to/assemblies/ \
  -o results/ \
  -s TRUE
```

| Flag | Required | Description                            | Default         |
| ---- | -------- | -------------------------------------- | --------------- |
| `-q` | Yes      | Query FASTA with mec/ccr genes         | `gene_db.fasta` |
| `-r` | Yes      | Directory containing genome assemblies | —               |
| `-o` | Yes      | Output directory                       | `./results`     |
| `-s` | No       | Keep intermediate files                | `TRUE`          |

---

## Output layout

```
results/
├── sam/               # raw & merged SAM alignments
├── bed/               # BED + genome size files
├── temp/              # intermediates (if -s TRUE)
├── search_stat.tsv    # per‑gene alignment stats
└── best_result.tsv    # inferred SCCmec type(s)
```

Both TSV files are ready for downstream parsing in R, pandas, etc.

---

## Gene database

Curated from the supplemental data of Iorio *et al.* 2022; stored in `gene_db.fasta` and reproduced here for convenience.

| SCCmec type | Representative strain | GenBank accession |
| ----------- | --------------------- | ----------------- |
| I           | NCTC 10442            | AB033763          |
| II          | N315                  | D86934            |
| III         | 85/2082               | AB037671          |
| IV          | CA05                  | AB063172          |
| V           | WIS                   | AB121219          |
| VI          | HDE288                | AF411935          |
| VII         | P5747/2002            | AB373032          |
| VIII        | C10682                | FJ390057          |
| IX          | JCSC6943              | AB505628          |
| X           | JCSC6945              | AB505630          |
| XI          | LGA251                | FR821779          |
| XII         | BA01611               | KR187111          |
| XIII        | 55‑99‑44              | MG674089          |
| XIV         | SC792                 | LC440647          |

---

## Workflow at a glance

1. **Align genes** — `minimap2` maps mec/ccr genes to each assembly (`sam/`).
2. **Filter & trim** — remove duplicates, fix contig‑edge truncations (≤3 iterations).
3. **Merge SAMs** — produce `<sample>_combined.sam`.
4. **Summarise coverage** — generate `search_stat.tsv` (percent length & contig).
5. **Call SCCmec** — choose best mec complex + ccr combo; merge `ccrA/B` if on same contig; break ties by proximity.
6. **Quality flags** — low‑coverage or ambiguous calls are flagged for review (more automation coming).

---

## Roadmap

* [ ] Publish on **conda‑forge**
* [ ] Add sub‑type support (e.g. IVa, IVb…)
* [ ] Incorporate mec class **D**
* [ ] CI tests across platforms
* [ ] Optional GUI wrapper

Contributions and pull requests very welcome!

---

## Key References

* **Database source:** Iorio, M. L., *et al.* (2022) *Front. Microbiol.* 13:826456 – comprehensive catalog of SCCmec elements.
* Li, H. (2018) *Minimap2: pairwise alignment for nucleotide sequences.* *Bioinformatics* 34, 3094‑3100 – alignment engine powering SCCmecClassifier.

---

## License

MIT – see `LICENSE` for full terms.
