# demulti: phasing indexed reads

The script phases paired reads based on an index table into an output directory. It also outputs general stats (indices with 'N's, not in index table, mixed-up indices), and an index confusion matrix (relevant only if the two sides had an index from the table).

```
Usage: demulti.pl <input table> <index length> <input R1> <input R2> <output dir> [output prefix id 1, output prefix id2, ...]
Input:
 <input table>: input table with columns (phase,read1index,read2index)
 <index length>: length of index in nt
 <input R1>: R1 fastq
 <input R2>: R2 fastq
 <output dir>: Write output to this directory
 <output prefix id N>: Zero or more optional prefix ids
Output files (in OUTDIR):
 OUTDIR/PHASE_R[12].fastq: R1 and R2 of phased fastq files, with an optional PREFIX that is generated using one or more supplied prefix ids
 OUTDIR/stats.txt: general stats of phasing
 OUTDIR/matrix.txt: phase confusion matrix
Example 1: 
%> perl demulti.pl examples/index.txt 7 examples/R1.fastq examples/R2.fastq examples/out
Example 2: 
%> perl demulti.pl examples/index.txt 7 examples/R1.fastq examples/R2.fastq examples/out run7 sample
```
