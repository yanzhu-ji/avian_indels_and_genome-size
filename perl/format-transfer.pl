#!/usr/bin/perl

# Description: 
#	tranfer BLAST (and HMMER) fmt6 output to bed format for better merging using BEDTools
# Usage:
# 	$0 (required) -inputformat <"fmt6", "gff3", "ucscCustomTable", "ucscCustomTableExon">  (see the end of script for specification).
#	(optional) -outputformat <"bed"> 

# file xxx.bed will be written automatically.

# YJ 8/4/14

use strict;
use Getopt::Long;

my ($inputType, $outputType);

GetOptions (
	'i|inputformat=s' => \$inputType,
	'o|outputformat=s' => \$outputType,
);

# set default values
$outputType ||= "bed";


my $infile = shift;

if ( $inputType eq "fmt6" && $outputType eq "bed" ){
	if ( $infile !~ /\.fmt6/){
		die "Please check the format of input file ", "'", $infile, "'", "\n";
	}
	fmt6_to_bed_print ( $infile );
}elsif ( $inputType eq "gff3" && $outputType eq "bed" ){
	if ( $infile !~ /\.gff3/ ){
		die "Please check the format of input file ", "'", $infile, "'", "\n";
	}
	gff3_to_bed_print ( $infile );
}elsif ( $inputType eq "ucscCustomTable" && $outputType eq "bed" ){
	if ( $infile !~ /\.txt/ ){
		die "Please check the format of input file ", "'", $infile, "'", "\n";
	}
	ucscCustomTable_to_bed_print ( $infile );
}elsif ( $inputType eq "ucscCustomTable" && $outputType eq "exon-bed") {
	if ( $infile !~ /.txt/ ){
		die "Please check the format of input file ", "'", $infile, "'", "\n";
	}
	ucscCustomTable_to_exon_print ( $infile );
}else {
	print STDERR "Please see the Usage!\n";
}

sub gff3_to_bed_print {
	my $infile = shift;
	my $outfile = $infile;
	$outfile =~ s/\.gff3/\.bed/;

	open (INFILE, "<", $infile ) or die ("Cannot open input file: $!\n");
	open (OUTFILE, ">", $outfile ) or die ("Cannot write to file: $!\n");
	my ( $strand, $starting_p, $ending_p );

	while ( my $line = <INFILE> ){
		next if ( $line !~ /\S/ );
		chomp $line;
		my @field = split " ", $line;
	
		$starting_p = $field[3] - 1;
		$ending_p = $field[4];
		$strand = $field[6];
		if ( $strand eq "?" ){
			$strand = ".";
		}
		my $name = join "_", ($field[1], $field[2], $field[8]);
		print OUTFILE "$field[0]\t$starting_p\t$ending_p\t$name\tna\t$strand\n";
		($starting_p, $ending_p, $strand, $name ) = "";
	}
	close INFILE;
	close OUTFILE;
	print STDOUT "Bed file $outfile has been written to the current directory.\n";
}

sub fmt6_to_bed_print {
	my $infile = shift;
	my $outfile = $infile;
	$outfile =~ s/\.fmt6/\.bed/;
	my ( $strand, $starting_p, $ending_p );

	open (INFILE, "<", $infile ) or die ("Cannot open input file: $!\n");
	open (OUTFILE, ">", $outfile ) or die ("Cannot write to file: $!\n");

#d	if ( $inputType eq "fmt6" && $outputType eq "bed" )

	while ( my $line = <INFILE> ){
		next if ( $line !~ /\S/ );
		chomp $line;
		my @field = split " ", $line;

		if ( $field[8] < $field[9] ){
			$starting_p = $field[8] - 1; # since starting positions of bed files are zero-based
			$ending_p = $field[9]; # but ending positions are one-based
			$strand = "+";
		}elsif ( $field[8] > $field[9] ){
			$starting_p = $field[9];
			$ending_p = $field[8] - 1;
			$strand = "-";
		}else {
			print STDERR "Weird: starting position equals to ending position!\n";
			print STDERR $line, "\n";
		}
		
		print OUTFILE "$field[1]\t$starting_p\t$ending_p\t$field[0]\tna\t$strand\n";
		($starting_p, $ending_p, $strand) = "";
	}
	close INFILE;
	close OUTFILE;
	print STDOUT "Bed file $outfile has been written to the current directory.\n";
}
#d	}elsif ($inputType eq "hmm-blah" ){
	# doing nothing temporarily
#d	}


sub ucscCustomTable_to_bed_print {
	my $infile = shift;
	my $outfile = $infile;
	$outfile =~ s/\.txt/\.bed/;

	open (INFILE, "<", $infile ) or die ("Cannot open input file: $!\n");
	open (OUTFILE, ">", $outfile ) or die ("Cannot write to file: $!\n");
	my ( $strand, $starting_p, $ending_p );

	while ( my $line = <INFILE> ){
		next if ( $line !~ /\S/ or $line =~ "#" );
		chomp $line;
		my @field = split " ", $line;
	
		$starting_p = $field[3];  
		$ending_p = $field[4];
		$strand = $field[2];
#		if ( $strand eq "?" ){
#			$strand = ".";
#		}
		my $chr = $field[1];
		my $id = $field[0];
		my $name2 = $field[11];
		my $score = $field[10];
		print OUTFILE "$chr\t$starting_p\t$ending_p\t$name2\t$score\t$strand\t$id\n";
		($chr, $starting_p, $ending_p, $name2, $strand, $score, $strand, $id ) = "";
	}
	close INFILE;
	close OUTFILE;
	print STDOUT "Bed file $outfile has been written to the current directory.\n";
}

sub ucscCustomTable_to_exon_print {
	my $infile = shift;
	my $outfile = $infile;
	$outfile =~ s/\.txt/_exon\.bed/;

	open (INFILE, "<", $infile ) or die ("Cannot open input file: $!\n");
	open (OUTFILE, ">", $outfile ) or die ("Cannot write to file: $!\n");

	while ( my $line = <INFILE> ){
		next if ( $line !~ /\S/ or $line =~ "#" );
		chomp $line;
		my @field = split " ", $line;

		
		my $chr = $field[1];

		my @exon_start = split ",", $field[8];
		my @exon_end = split ",", $field[9];

		if ( $#exon_start != $#exon_end ){
			print STDERR "Error: number of starts is not equal to number of ends: $line\n";
		}

		my $strand = $field[2];
		my $gene_id = $field[0];
		my $name2 = $field[11];
		my $score = $field[10];
		
		foreach my $i ( 0..$#exon_start ) {
			print OUTFILE "$chr\t$exon_start[$i]\t$exon_end[$i]\t$name2\t$score\t$strand\t$gene_id\n";
		}
#d		print OUTFILE "$chr\t$starting_p\t$ending_p\t$name2\t$score\t$strand\t$id\n";
		($chr, @exon_start, @exon_end, $name2, $strand, $score, $gene_id ) = "";
	}
	close INFILE;
	close OUTFILE;
	print STDOUT "Bed file $outfile has been written to the current directory.\n";
}


#Updates

# option ucscCustomTable is added 3/11/2016 to facilitate orthologous gene analysis near CR1 loci.
# 	ucscCustomTable is of the format "name, chrom, strand, txStart, txEnd, score, name2" 
#	as specified from UCSC Table Browser of RefGene (or xenoRefGene).

# option ucscCustomTableExon is added 3/22/2016 

# option ucscCustomTableExon is changed to "-i ucscCustomTable -o exon-bed" 4/4/2016


