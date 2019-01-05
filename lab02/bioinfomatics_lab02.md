# Bioinfomatics Lab02
## Summary
**object** :  get local, global alignment to measure sequence similarity
**author** : Sujeong Shim (2016025532)
**submit** : 2018.10.21
**environment** : mac os & python3.6

## Usage
	python3 lab02.py [sequence file 1] [sequence file 2] [score file] [output file]
* This program has argument including two sequence files name, score file name output file name.
* `sequence files`, `score file` name is neccesary, but `output file name` is not neccesary. If `output file name` is not given, it has default name `output.txt`.
* This program perform global alignment algorithm, local alignment algorithm at the same time.

## Input & Output File
#### input file
* Sequence File
  * first line : start with `>` and include sequence name
  * second line : DNA sequence consising of A, C, G and T
* Score File
  * Match, mismatch, gap score with format like below
    * match=`match score`
    * mismatch=`mismatch score`
    * gap=`gap score`

#### output file
* This program make two file start with each `global_`, `local_`. If file name is output.txt, this program generate global_output.txt, local_ouput.txt
* File Content
  * match, mismatch, gap count
  * score
  * result of alignment
    * If filename start with `global_`, it contains result of global alignment.
    * If filename start with `local_`, it containes result of local alignemnt.

## Other
* This program generate alignments only consisting of upper case. If you input lower, it will be changed to upper case