#!/usr/bin/perl
#
# Description: clip sequence from bed files; bed files should be filtered for length and/or E-value, and merged using bedtools merge to reduce the issues this script needs to address.  In other words, this scripts only clips sequence with proper "buffer" lengths before or after.
#
# Usage: 
# (*: mandatory input files)
#   $0 	    -i|input* 	<bed_file>
# 	    -s|sequence* 	<source_fasta>
#	    -m|mode	"both" (default), "5" or "3"
# 	    -b|before 	<int, nt or aa buffer before>
# 	    -a|after 	<int, nt or aa buffer after>
# 	    -e|extend	"partial" (default) or "full"
# 	    -p|print	<optional with file name.  If "-p", then print out the corresponding list seperated by tab>
# 	    > 		<output_fasta>
# 
# 
# YJ. 10/31/2014
# last modified 6/11/15

use strict;
use Getopt::Long;
use Bio::SeqIO;


## set variables

my ( $infile, $seqfile, $mode, $ext_b, $ext_a, $extend, $print );

GetOptions(
	"i|input=s" => \$infile,
	"s|sequence=s" => \$seqfile,
	"m|mode=s" => \$mode,
	"b|before=i" => \$ext_b, 
	"a|after=i" => \$ext_a,
	"e|extend=s" => \$extend,
	"p|print=s" => \$print,
);


if ( $infile eq "" || $seqfile eq "" ){
	die ("Please check the usage and give the correct input files!\n");
}

$ext_b ||= 0;
$ext_a ||= 0;
$mode ||= "both";
$extend ||= "partial";

## get sequence length in a hash
my %hseq_length;
my $in = Bio::SeqIO->new(-file => "$seqfile",
			 -format => 'Fasta');
						 
LOOP: while ( my $seq = $in->next_seq() ){
	my $seqID = $seq->display_id();
	$hseq_length{$seqID} = $seq->length();
}

#### test
#foreach my $i ( keys %hseq_length ){
#	print STDOUT "$i: $hseq_length{$i}\n";
#}


## parse bed file and switch the coordinates to 1-based instead of 0-based (bed)
my %clip_info;

open (BED, "<", $infile ) or die ("Cannot open input bed file: $!\n");
while ( my $line = <BED> ){
	chomp $line;
	my @fields = split "\t", $line;
#t	print STDOUT "$fields[0], $fields[1], $fields[2]\n";
	my $length = $hseq_length{$fields[0]};
	my ( $o_start, $o_end );

	# get start and end positions based on the -m option (1-based)
	if ( $mode eq "both" || $mode eq "5" ){
		$o_start = $fields[1] + 1; # to shift the starting positions to 1-based while end positions remain unchanged
	}elsif ( $mode eq "3" ){
		$o_start = $fields[2];
	}else {
		die "Incorrect -m parameter 1: both, 5 or 3\n";
	}

	if ( $mode eq "both" || $mode eq "3" ){
		$o_end = $fields[2];
	}elsif ( $mode eq "5" ){
		$o_end = $fields[1] + 1;
	}else {
		die "Incorrect -m parameter 2: both, 5 or 3\n";
	}

	my ( $n_start, $n_end );
	
#t	print STDOUT "Checking: $fields[0] length:  $hseq_length{$fields[0]} \n";
	unless ( $hseq_length{$fields[0]} ){
		print STDERR "Warning: sequence length for $fields[0] cannot be found; please check the input files!\n";
		next;
	}

	if ( $o_start  - $ext_b > 0 ){
		$n_start = $o_start  - $ext_b;
	}elsif ( $extend eq "partial" ){
		$n_start = 1;
	}else {
		print STDERR "Cannot extend 5' to the end: $fields[0] at $o_start: skipping...\n";
		next;
	}

	if ( $o_end  + $ext_a <= $length ){
		$n_end = $o_end + $ext_a - 1;
	}elsif ( $extend eq "partial" ){
		$n_end = $length - 1;
	}else {
		print STDERR "Cannot extend 3' to the end: $fields[0] at $o_end: skipping...\n";
		next;
	}

#### test
#t	print STDOUT "$fields[0],$n_start,$n_end\n";

	push @{ $clip_info{$fields[0]} }, [ $n_start, $n_end, $fields[1], $fields[2] ]; # added 5.8; modified 5.9...and 6.11
}

#### test
#foreach my $i ( keys %clip_info ){
#	foreach my $m ( @{$clip_info{$i}} ){
#		print STDOUT "$i: @$m\n";
#	}
#}


if ( $print ){
	open (TABLIST, ">", $print ) or die ("Cannot write to file: $!\n");
	print TABLIST "#bed_info\tclipped_seq\n"; # as a line of annotation, print what the two columns are in the file.
}

## read source fasta and clip the sequence.

my $seqin = Bio::SeqIO->new(-file => $seqfile,
			    -format => 'Fasta',
			);

while ( my $seq = $seqin->next_seq() ){
	my $source_seqID = $seq->display_id();
	if ( $clip_info{$source_seqID} ){
		foreach my $m ( @{ $clip_info{$source_seqID} } ){
			print STDOUT ">$source_seqID"."_$$m[0]"."_$$m[1]\n",
			$seq->subseq($$m[0], $$m[1]), "\n";
			if ( $print ){
				print TABLIST "$source_seqID"."_$$m[2]"."_$$m[3]","\t$source_seqID"."_$$m[0]"."_$$m[1]\n";
			}
		}
	}
}
