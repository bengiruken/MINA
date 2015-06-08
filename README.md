MINA
====

Integrative network analysis for gene-gene interactions associated with clinical outcome across multiple genomic profiles

## A brief instruction how to use MINA

### Fetching source code from git repository

<pre>
git clone https://github.com/hhjeong/MINA
</pre>

### Compilation of source codes

<pre>
cd MINA
make
</pre>

### Running MINA with toy example

In windows
```
test.bat
```

In linux
```
./test.sh
```

### Detail of the command line parameter

* `-s`: name of gene symbol file
* `-ip`: name of a profile file
* `-io`: name of a clincial outcome file
* `-o`: path of the output
* `dist` or `network`: type of output
* `-perm`: number of iterations to get threshold values
* `-alpha`: stringent parameters to construct integrative networks
* `-dlo`, `-dhi` : boundary of distribution assessment

## Running examples

* Assessing mutual information distribution for each profile
```
./MINA -s symbol.txt -ip A.txt -ip B.txt -io outcome.txt -o output/ dist
```

* Constructing integrative network of the profiles
```
./MINA -s symbol.txt -ip A.txt -ip B.txt -io outcome.txt -maxperm 30 -alpha 0.8
```