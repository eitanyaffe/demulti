# demulti

```
Usage: demulti.pl <input table> <index length> <input R1> <input R2> <output dir>
Input:
 <input table>: input table with columns (phase,read1index,read2index)
 <index length>: length of index in nt
 <input R1>: R1 fastq
 <input R2>: R2 fastq
 <output dir>: Write output to this directory
Output files (in OUTDIR):
 OUTDIR/PHASE_R[12].fastq: R1 and R2 of phased fastq files
 OUTDIR/stats.txt: general stats of phasing
 OUTDIR/matrix.txt: phase confusion matrix
Example: 
%> perl demulti.pl examples/index.txt 7 examples/R1.fastq  examples/R2.fastq  examples/out
```
