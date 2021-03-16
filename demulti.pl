#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <input table> <index length> <input R1> <input R2> <output prefix> <output confusion matrix> <output summary>\n";
	exit 1;
}

my $ifn_table = $ARGV[0];
my $ilength = $ARGV[1];
my $ifn_R1 = $ARGV[2];
my $ifn_R2 = $ARGV[3];
my $oprefix = $ARGV[4];
my $ofn_mat = $ARGV[5];

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
    
    my $ofn1 = $oprefix."_R1_".$phase;
    my $ofn2 = $oprefix."_R2_".$phase;
    
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

open(IN1, $ifn_R1) or die;
open(IN2, $ifn_R2) or die;

my $phase = -1;
my $iindex = 0;
my @reads1;
my @reads2;

while () {
    my $line1 = <IN1>;
    my $line2 = <IN2>;
    my $sindex = $iindex % 4;
    
    last if (!defined($line1) || !defined($line2));
    $reads1[$sindex] = $line1;
    $reads2[$sindex] = $line2;

    # determine phase
    if ($sindex == 1) {
	my $index1 = substr($line1, 0, $ilength);
	my $index2 = substr($line2, 0, $ilength);

	my $phase1 = defined($indices1{$index1}) ? $indices1{$index1} : -1;
	my $phase2 = defined($indices2{$index2}) ? $indices2{$index2} : -1;
	
	# skip
	if ($phase1 == -1 or $phase2 == -1) {
	    next;
	}
	
	if ($phase1 == $phase2) {
	    $phase = $phase1;
	} else {
	    $mat{$phase1}->{$phase2} = 0;
	}
	
    }

    # output sequence 
    if ($sindex == 3) {
	if ($phase != -1) {
	    defined($phases{$phase}) or die "phase not defined";
	    for (my $i=0; $i<4; $i++) {
		my $fh1 = $phases{$phase}->{fh1};
		my $fh2 = $phases{$phase}->{fh2};
		print $fh1 $reads1[$i];
		print $fh2 $reads2[$i];
	    }
	}
    }
    
    $iindex++;
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

open(OUT, $ofn_mat) or die;
print "generating confusion matrix: $ofn_mat\n";
print OUT "phase1\tphase2\tcount\n";
foreach my $index1 (keys %indices1) {
    my $phase1 = $indices1{$index1};
    $mat{$phase1} = {};
    foreach my $index2 (keys %indices2) {
	my $phase2 = $indices2{$index2};
	$mat{$phase1}->{$phase2} = 0;
	print OUT $phase1, "\t", $phase2, "\t", $mat{$phase1}->{$phase2}, "\n";
    }
}
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
