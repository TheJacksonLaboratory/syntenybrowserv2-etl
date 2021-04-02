### Synteny Browser Database Creation Directory
This directory contains all of the scripts required to load a fully-functional database for running
the Synteny Browser application.

A full set of source files for the Mus musculus, Homo sapiens and Rattus norvegicus species is available for download at [this link](http://www.informatics.jax.org/downloads/SyntenyBrowser/). The downloaded files can be added to a separate directory, such as './data', for example.

A more in-depth breakdown of the data files and the data each contains can be found [here](http://syntenybrowser.jax.org/docs/data-prep) and directions for loading a database can be found here [here](../README.md).

### File Purpose Summaries
Below are brief descriptions of the scripts used to load a database.

* `import_cytogenetic_loc.py` - parses cytogenetic bands data from GFF3 formatted source and loads the records into a 'cytogenetic_band' table
* `import_features.py` - parses feature data from GFF3 formatted source and loads the records into a 'feature' table. The 'feature' table contains features that are not genes, such as QTLs, mRNA, etc.
* `import_genes_exons.py` - parses gene and exon data from GFF3 formatted source and loads the records into 'gene' and 'exon' tables
* `import_homologs.py` - parses homologs data from a custom, TAB separated source and loads the records into a 'homolog' table
* `import_ontology.py` - parses ontology terms and mappings data from OBO and GAF formatted sources and loads the records into a 'on_terms', 'on_pairs', and 'gene_ontology_map' tables
* `import_synteny_blocks.py` - parses syntenic blocks data from custom, TAB separated source and loads the records into a 'syntenic_block' table
* `import_variants.py` - parses SNP variants data from VCF formatted source and loads the records into a 'snp_variant' table.
* `flex_open.py` - contains a utility function that assists in opening .gz and non-.gz compressed files