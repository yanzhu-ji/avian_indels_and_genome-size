#!/usr/bin/perl
#
# Description: modify/update column values of a certain file based on a pre-defined reference file, i.e. all column values of the certain file will be modified to equivalent names.
#
# USAGE: $0 -i|infile	<file_to_be_modified,key_col,target_col>
# 	    -r|reference	<reference_table,reference_key_col,ref_value_col>  (e.g. "table",2,3)
#	    -p|print 	(if present, then print those without matches as well)
# 	    > 	<blah_out>
#
# Note: 1. the reference table is a tab-delimited file with multiple columns.
# 	2. content of key_col and ref_key_col should match; otherwise the unmatched rows will be deleted.
#	3. if headers (or other parts of the file to be modified) have sequence names, they will not be modified by the script (i.e. manual justification is needed in such cases).
# 
# YJ 11/7-8, 2014

use strict;
use Getopt::Long;

#### 0. read in options
my ( $in, $ref, $print );
GetOptions(
	"i|infile=s" => \$in,
	"r|reference=s" => \$ref,
	"p|print" => \$print,
);

my ( $infile, $in_key, $in_value) = split ",", $in;
my ( $reference, $ref_key, $ref_value ) = split ",", $ref;

# if ( $infile eq "" or $in_key eq "" or $ref_value eq "" or $reference eq "" or $ref_key eq "" or $ref_value eq "" ){
#	die ( "Abort: incorrect input.  Please check the usage!\n" );
#}


#### 1. read in reference table and create a hash with the column to be modified as keys
open (REF, "<", $reference ) or die ( "Cannot open $reference file: $!\n" );

my %info;
# check for uniqueness?
while ( my $line = <REF> ){
	chomp $line;
	if ( $line =~ /^#/ ){ next; }  # added 6/12/15 for the syntenic chromosomes txt file.
	my @field = split " ", $line; # edited 4/29: "\t" -> " ".  What the hell..
	my $index = $ref_key - 1;
	$info{$field[($ref_key - 1)]} = $field[($ref_value - 1 )];
}

close REF;


# foreach my $m (keys %info ){
#	print STDERR "$m\t$info{$m}\n";
# }

#### 2. read in infile and modify the selected column based on the hash

open (IN, "<", $infile ) or die ( "Cannot open $infile: $!\n" );

my $skip = 0;

print STDERR "Warning: equivalent name for the list below do(es) not exist.\n";
while ( my $line = <IN> ){
	chomp $line;
	if ( $line =~ /^#/ ){
		print STDOUT $line, "\n";
		next;
	}

	# if line is separated by one or more tabs, then cut fields by tabs as they may indicate empty columns that should not be combined.
	# else, cut fields by combined spaces 
	my @field;
	if ( $line =~ /\t{1,}/ ) {
		@field = split "\t", $line; 
	}else {
		@field = split " ", $line; 
	}
#	print STDERR "parsing:$line\n";

	if (  $info{$field[($in_key-1)]} ){
		$field[($in_value - 1 )] = $info{$field[($in_key-1)]};
		print STDOUT (join "\t", @field), "\n";
	}elsif ( $print ){
		$field[($in_value - 1 )] = "";
		print STDOUT (join "\t", @field), "\n";
		print STDERR "$field[($in_key-1)]\n";
		$skip ++;
	}else{
		print STDERR "$field[($in_key-1)]\n";
		$skip ++;
	}
}

close IN;

print STDERR "Total lines skipped: $skip.\n";

# updates
# "-p" option is added 3/28/2016