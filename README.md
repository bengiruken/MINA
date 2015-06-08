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

### Details of Input
* symbol file - It consists of N lines with symbols of features
```
feature1
feature2
feature3
feature4
feature5
feature6
feature7
feature8
feature9
feature10
```

* outcome file - This file contains M binary outcomes 
 * i.e. 0 represents long-term survival patients, 1 represents short-term survival patients
```
1
0
1
1
0
0
```

* profile files - The type of files consist of N rows with M columns 
 * A value at i-th row of j-th column corresponds to measure of i-th feature for j-th patient
```
0.437280581	0.54533823	0.674094347	0.621944035	0.635085602	0.544135653
0.806529961	0.700986151	0.489553203	0.9280643	0.600577472	0.857410819
0.116087783	0.09777806	0.146282316	0.081376467	0.098742532	0.085923183
0.098671101	0.01743013	0.011476686	0.034055921	0.09660923	0.031223755
0.89914231	0.965998414	0.976762825	0.972373626	0.830364778	0.935723614
0.022942364	0.019296143	0.020128928	0.015178173	0.023613841	0.021235995
0.023414639	0.026617573	0.033094332	0.021436811	0.016342239	0.023211977
0.040137216	0.030910384	0.037995903	0.029218459	0.029634722	0.028493865
0.015701958	0.013569464	0.020047829	0.015802287	0.014216741	0.014052416
0.936975594	0.855736673	0.98083047	0.882792997	0.803765363	0.928080435
```

* Constraints
 * Number of rows of the profile files should be same as number of rows of the symbol file.
 * Number of columns of the profile files should be same as number of rows of the outcome file.


### Running MINA with toy example

In windows
```
test.bat
```

In linux
```
./test.sh
```

### Details of the command line parameter

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