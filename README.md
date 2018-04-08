# Cw1_Malovichko

## Description
A simple CC tool for raw read trimming (for Illumina-Sanger coding) imitating basic Trimmomatic features.

## CC Syntax
**Overview**:
```bash
python Cw_1_Malovichko.py -i <input> -he[int] -c[int] -swl[int] -swq[int] -o [str]
```
**Arguments:**
+ **-i**: input file name/path.
+ **-he**: an integer defining number of nucleotides to be cut from the beginning (5'-end) of the read (an imitation of Trimmomatics HEADCROP). Default is 0
+ **-c**: an integer defining number of nucleotides to be cut from the ending (3'-end) of the read (an imitation of Trimmomatic's CROP). Default is 0.
+ **-swl**: an integer standing for size of sliding window. Default is 0.
+ **-swq**: an integer for sliding window baseline. If average quality within the window falls below the baseline, read is cut by the first nucleotide (not included). Default is 0.
+ **-o**: an optional argument for output file name/path. Default is 0. If not stated, returns formatted *.fastq* entries to stdout.
*Note: this function seems to be not optimized and writing file with from stdout with '>' operator goes faster than **-o** specification.* 

## Code structure
Code includes one *argparse* initializing chunk and four functions for input processing.

### HEADCROP
A function for HEADCROP imitation triggered by **-he** flag value. Requires two arguments: a SeqIO object (*read*) and number of nucleotides trimmed from the 5'-end (*drop_value*). Returns a formatted SeqIO object (*new_read*). Also, contains a two-line commentary showing an alternative (feature-by-feature) way for creating a SeqIO object.

### CROP
A CROP imitation triggered by **-c** flag value. Requires two arguments: a SeqIO object (*read*) and number of nucleotides trimmed from the 3'-end (*drop_value*). Returns a formatted SeqIO object (*new_read*). 

### SLIDINGWINDOW
A function for sliding window performance. Requires three arguments: a SeqIO object (*read*), size of the window provided by **-swl** argument (*length*) and sequence quality baseline provided by **-swq** argument )=(*baseline*). Performs CROP function once average quality falls below the baseline. Returns a formatted SeqIO object (*read*).

### trimming
Sequentially performs three aforementioned functions over every read, examines it whether its length is larger than 10 and either returns it to stdout or adds to list for further infile writing. Requires exactly one argument - an input file (*file*).

## Example
Two tests with to different read files were performed with the following command:
```bash
python Cw_1_Malovichko.py -i ~/Downloads/test_classwork{number}.fastq -he 3 -c 3 -swl 4 -swq 30 -o ~/Downloads/output_test{number}.fastq
```
An HTML FastQC report for all four files is stored at this repository.
