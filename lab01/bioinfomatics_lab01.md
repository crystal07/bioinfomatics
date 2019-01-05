# Bioinfomatics Lab01
## Summary
**object** :  get motif from gibbs sampling
**author** : Sujeong Shim(2016025532)
**submit** : 2018.09.23
**environment** : mac os & python3.6

## Usage
	python3 lab01_gibbs.py [input file] [output file] [sampling count]
* This program has argument including input file name, output file name, sampling count.
* 'input file name', 'output file' name is neccesary, but 'sampling count' is not neccesary.
* But this program perform gibbs sampling as much as sampling count and choose motif that has best score in each sampling.
* If you don't input sampling count, this perform only 1 sampling.

## Input & Output File
#### input file
* __k__ : the length of motif
* __t__ : the number of string in a input file
* __m__ : max length of each string
* strings

#### output file
* motifs
* a profile with t motifs found
 * this profile matrix is sorted
 * order A, C, G, T


## Other
* This program is slow. Please, wait until program end normally
* This program generate motifs only consisting of upper case
* This program iterate 3000 in each sampling