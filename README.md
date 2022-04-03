# The bestrhodopsin analysis pipeline

This repository includes a [Snakemake](https://snakemake.readthedocs.io/) workflow used for bioinformatic analysis for the paper "Rhodopsin-bestrophin fusion proteins from unicellular algae form gigantic pentameric ion channels". Please refer to the paper for further details. Additional data are available from the paper and [from this repository](https://doi.org/10.5281/zenodo.5119843).

The entry point is `workflow/Snakemake`, so the workflow can be launched by simply calling `snakemake -c[number_of_cores]` in the root directory. Most of the dependencies are taken care of with `conda` (see `workflow/envs` for details) -- use `snakemake`'s `--use-conda`. However, a number of third-party programs have to be installed manually (the used versions and the executables expected to be in the `$PATH`):

* [Geneconv](https://www.math.wustl.edu/~sawyer/geneconv/index.html#progfiles) v1.81a - `geneconv`
* [Root Digger](https://github.com/computations/root_digger) v1.7.0-7-gccbe87e - `rd`
* [busco](https://gitlab.com/ezlab/busco/-/releases) v5.beta.1 (for other versions conda can be used) - `busco`
* [TtreeShrink](https://github.com/uym2/TreeShrink) v1.3.6 - `run_treeshrink.py`
* [ASTRAL](https://github.com/smirarab/ASTRAL) v5.7.4 - specify the path to `astral.version.jar` in `config/config.yaml` (also notice that `java` is not installed via conda)
* [ERaBLE](http://www.atgc-montpellier.fr/erable/usersguide.php) v1.0 - `erable`
* [SequenceBouncer](https://github.com/corydunnlab/SequenceBouncer) v1.18 - `SequenceBouncer.py`
* [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-4.1) v4.1 - `signalp` (versions >4 are not supported)
* [TargetP](https://services.healthtech.dtu.dk/service.php?TargetP-2.0) v2.0 - `targetp`
* [ASAFind](https://rocaplab.ocean.washington.edu/tools/asafind/) v1.1.7 - `ASAFind.py`
* [InterProScan](https://www.ebi.ac.uk/interpro/download/) v5.48-83.0 (available from conda, but we keep a global installation) - `interproscan.sh`

## The workflow

The workflow consists of four independent components:

* Codon-based analyses:
- codon and protein phylogenies of the domains
- different rooting strategies
- recombination analyses

* Global phylogenies for:
- rhodopsins and
- bestrophins

* Species phylogenies for:
- chlorophytes,
- dinoflagellates and
- haptophytes

* Structural alignment for:
- rhodopsins and
- bestrophins

## Input files

The input files are provided in the `input/` folder. In particular:

* inputs for the codon-based analysis are codon sequences for bestrhodopsin domains: `input/codons/domains/bestrophins.fasta` and `input/codons/domains/bestrophins.fasta`. Included in the fasta files are also outgroup sequences trimmed to homologous regions.

* inputs for the species phylogenies are expected in the `input/species/` directory, but the complete fasta files with all protein sequences for the assemblies are not included in the repository. Instead the pipeline can be started from the extracted orthogroups that can be downloaded from [the data repository](https://doi.org/10.5281/zenodo.5119843).

* inputs for the global bestrophin phylogeny were downloaded from Uniref and Pfam and are provided as `input/global_phylogeny/bestrophins/uniprot_aln.fasta` (Pfam bestrophin alignment), `input/global_phylogeny/bestrophins/uniref50.fasta` (Uniref50 sequences of bestrophins), `input/global_phylogeny/bestrophins/uniref50.txt` (list of bestrophins from Uniref50 in tabular-format).

* inputs for the global rhodopsin phylogeny are also available in `input/global_phylogeny/rhodopsins`. The file `rhodopsin_selected_8TMs_from_uniref50.fasta` is derived from uniref50 sequences matched to the microbial rhodopsin Pfam profile.

## Output files

The output files are collected in the `output/` folder:

* `output/codons/rooted_trees.svg` includes rooted domain trees:
- of the rhodopsin domains based on the codon alignment and the same topology with branch lengths optimized based on the protein alignment;
- of the bestrophin domains - protein tree including two outgroup sequences (A0A6T5JQZ6 and A0A6T5CHT1), codon tree with TGD-B, codon tree w/o TGD-B and the same tree topology with branch lengths optimized based on the protein alignment. The last three trees are rooted with Root Digger with the numbers next to the root indicating Likelihood Weight Ratio of the root placement.

* `output/codons/rhodopsin_recombination.svg` covers recombination analyses, from top to bottom:
- GENECONV analysis
- pairwise nucleotide identity profiles for selected sequences averaged over sliding windows of 15 bp
- GARD analysis

* `output/global_phylogeny/` includes global phylogeny trees (newick and annotated pdf)

* `output/species/*.svg` are visualizations of the species trees
* `output/species/species_chronos.nwk` chronos-scaled species trees in newick format
