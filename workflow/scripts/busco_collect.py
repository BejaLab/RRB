
from glob import glob
from os import path, mkdir
from collections import defaultdict
from Bio import SeqIO

buscos = snakemake.input
output = str(snakemake.output)

min_species = snakemake.params['min_species']

files = defaultdict(list)
names = {}
for busco in buscos:
	single_copy_files = glob(path.join(busco, '*.faa'))
	name = path.normpath(busco).split(path.sep)[-4] # {path}/{name}/run_{lineage}/busco_sequences/single_copy_busco_sequences/
	for fname in single_copy_files:
		orthogroup, ext = path.splitext(path.basename(fname))
		files[orthogroup].append(fname)
		names[fname] = name

assert names, "No single_copy_busco_sequences found"

mkdir(output)
for orthogroup, fnames in files.items():
	if len(fnames) >= min_species:
		dir_name = path.join(output, orthogroup)
		file_name = path.join(dir_name, 'sequences.faa')
		mkdir(dir_name)
		with open(file_name, 'w') as fh:
			for fname in fnames:
				fasta = SeqIO.parse(fname, 'fasta')
				record = next(fasta)
				record.id = names[fname]
				SeqIO.write(record, fh, 'fasta')
