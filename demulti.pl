#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
    print "Usage: $0 <input table> <index length> <input R1> <input R2> <output dir>\n";
    print "Input:\n";
    print " <input table>: input table with columns (phase,read1index,read2index)\n";
    print " <index length>: length of index in nt\n";
    print " <input R1>: R1 fastq\n";
    print " <input R2>: R2 fastq\n";
    print " <output dir>: Write output to this directory\n";
    print "Output files (in OUTDIR):\n";
    print " OUTDIR/PHASE_R[12].fastq: R1 and R2 of phased fastq files\n";
    print " OUTDIR/stats.txt: general stats of phasing\n";
    print " OUTDIR/matrix.txt: phase confusion matrix\n";
    print "Example: \n%> perl demulti.pl examples/index.txt 7 examples/R1.fastq  examples/R2.fastq  examples/out\n";
    exit 1;
}

my $ifn_table = $ARGV[0];
my $ilength = $ARGV[1];
my $ifn_R1 = $ARGV[2];
my $ifn_R2 = $ARGV[3];
my $odir = $ARGV[4];

(system('mkdir -p '.$odir) == 0) or die "unable to create directory $odir";

###############################################################################################
# read index table
###############################################################################################

print "reading index table: $ifn_table\n";
open(IN, $ifn_table) || die $ifn_table;
my $header = <IN>;
my %h = parse_header($header);

my %phases;
my %indices1;
my %indices2;

while (my $line = <IN>) {
    chomp $line;
    my @f = split("\t", $line);
    my $phase = $f[$h{phase}];
    my $i1 = $f[$h{read1index}];
    my $i2 = $f[$h{read2index}];
    ((length($i1) == $ilength) and (length($i2) == $ilength)) or die "expected index length: $ilength";
    $phases{$phase} = {};
    $phases{$phase}->{R1} = $i1;
    $phases{$phase}->{R2} = $i2;

    $indices1{$i1} = $phase;
    $indices2{$i2} = $phase;
    
    my $ofn1 = $odir."/".$phase."_R1.fastq";
    my $ofn2 = $odir."/".$phase."_R2.fastq";
    
    open(my $fh1, ">", $ofn1);
    open(my $fh2, ">", $ofn2);
    $phases{$phase}->{fh1} = $fh1;
    $phases{$phase}->{fh2} = $fh2;
}

###############################################################################################
# init confusion matrix
###############################################################################################

my %mat;
foreach my $index1 (keys %indices1) {
    my $phase1 = $indices1{$index1};
    $mat{$phase1} = {};
    foreach my $index2 (keys %indices2) {
	my $phase2 = $indices2{$index2};
	$mat{$phase1}->{$phase2} = 0;
    }
}

###############################################################################################
# travsere fastq
###############################################################################################

print "reading fastq R1: $ifn_R1\n";
print "reading fastq R2: $ifn_R2\n";

print "generating output in: $odir\n";

open(IN1, $ifn_R1) or die;
open(IN2, $ifn_R2) or die;

my $phase = -1;
my $iindex = 0;
my @reads1;
my @reads2;

my %stats;
$stats{total} = 0;
$stats{Ns} = 0;
$stats{no_index} = 0;
$stats{ok} = 0;
$stats{mixup} = 0;

while () {
    my $line1 = <IN1>;
    my $line2 = <IN2>;
    my $sindex = $iindex % 4;
    $iindex++;
    
    last if (!defined($line1) || !defined($line2));
    $reads1[$sindex] = $line1;
    $reads2[$sindex] = $line2;

    # determine phase
    if ($sindex == 1) {
	my $index1 = substr($line1, 0, $ilength);
	my $index2 = substr($line2, 0, $ilength);

	$stats{total}++;
	
	if ((index($index1, 'N') != -1) || (index($index2, 'N') != -1)) {
	    $stats{Ns}++;
	    $phase = -1;
	    next;
	}
	
	my $phase1 = defined($indices1{$index1}) ? $indices1{$index1} : -1;
	my $phase2 = defined($indices2{$index2}) ? $indices2{$index2} : -1;
	
	# non-index
	if ($phase1 == -1 or $phase2 == -1) {
	    $stats{no_index}++;
	    next;
	}

	$mat{$phase1}->{$phase2}++;
	
	if ($phase1 == $phase2) {
	    $stats{ok}++;
	    $phase = $phase1;
	} else {
	    $stats{mixup}++;
	    $phase = -1;
	}
	
    }

    # output sequence 
    if ($sindex == 3 && $phase != -1) {
	defined($phases{$phase}) or die "phase not defined";
	for (my $i=0; $i<4; $i++) {
	    my $fh1 = $phases{$phase}->{fh1};
	    my $fh2 = $phases{$phase}->{fh2};
	    print $fh1 $reads1[$i];
	    print $fh2 $reads2[$i];
	}
	$phase = -1;
    }
}

close(IN1);
close(IN2);

foreach my $phase (keys %phases) {
    my $fh1 = $phases{$phase}->{fh1};
    my $fh2 = $phases{$phase}->{fh2};
    close($fh1);
    close($fh2);
}

######################################################################################################
# output confusion
######################################################################################################

my $ofn_mat = $odir."/matrix.txt";

open(OUT, ">", $ofn_mat) or die;

# header
foreach my $phase (sort { $a <=> $b} keys %phases) {
    print OUT "\t", "i".$phase;
}
print OUT "\n";

foreach my $phase1 (sort { $a <=> $b} keys %phases) {
    print OUT "i".$phase1;
    foreach my $phase2 (sort { $a <=> $b} keys %phases) {
	print OUT "\t", $mat{$phase1}->{$phase2};
    }
    print OUT "\n";
}
close(OUT);

######################################################################################################
# stats
######################################################################################################

my $ofn_stats = $odir."/stats.txt";

open(OUT, ">", $ofn_stats) or die;
print OUT "type\tcount\n";
print OUT "Ns_found", "\t", $stats{Ns}, "\n";
print OUT "index_not_found", "\t", $stats{no_index}, "\n";
print OUT "index_mixup", "\t", $stats{mixup}, "\n";
print OUT "ok", "\t", $stats{ok}, "\n";
print OUT "TOTAL", "\t", $stats{total}, "\n";
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
