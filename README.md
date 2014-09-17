MINA
====

Integrative network analysis for gene-gene interactions associated with clinical outcome across multiple genomic profiles

## A brief instruction how to use MINA

### Fetching source code from git repository

<pre>
git clone https://github.com/hhjeong/MINA
</pre>

### Compliation of source codes

<pre>
cd MINA
make
</pre>

### Running MINA with toy example

<pre>
cd bin
./MINA
</pre>

### Setting parameters
* geneinfo
 * A file path for symbols of genes for every profiles.
* profiles
 * File paths for expressions/alterations of genomic profiles.
* clinical
 * A file path for outcomes for each patients.
* maxPerm
 * Number of iterations for permutation tests.
* alpha
 * alpha threshold of single profile / intersection network / union network.

### Example of param.txt
<pre>
geneinfo: ../sample/symbol.txt
profiles: ../sample/CNA.txt ../sample/mRNA.txt ../sample/METH.txt
clinical: ../sample/clinical.txt
maxPerm: 30
alpha: 0.0 0.7 1.0
</pre>

