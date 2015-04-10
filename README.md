# Clone-based MSA Builder

Build a multiple sequence alignment (MSA) from a database of clones and a given
query sequence.

## Usage

Copy `config.template.json` to `config.json` and modify clone names, clone
paths, query sequences, and repeat sequences for your data set.

Run the MSA builder by running the following commands.

```bash
. config.sh
snakemake
```

This will produce a multiple sequence alignment per species and a pairwise
identity plot across each alignment.
