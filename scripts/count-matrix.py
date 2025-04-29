import pandas as pd

def get_column():
    if str(snakemake.params.strand) == "False":
        return 1 #non stranded protocol
    elif str(snakemake.params.strand) == "True":
        return 2 #3rd column
    elif str(snakemake.params.strand) == "Reverse":
        return 3 #4th column, usually for Illumina truseq
    else:
        raise ValueError(("'strandedness' should be False, True or Reverse, instead has the value {}").format(repr(snakemake.params.strand)))


counts = [pd.read_table(f, index_col=0, usecols=[0, get_column()], 
          header=None, skiprows=4) 
          for f, snakemake.params.strand in zip(snakemake.input, [snakemake.params.strand]*len(snakemake.input))]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
