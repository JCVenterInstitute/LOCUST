#!/usr/bin/env perl
#Copy (C) 2016 The J. Craig Venter Institute (JCVI).  All rights reserved

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

use File::Slurp;
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

create_itol.pl -  A script to generate ITOL annotations based upon ST
 									derived by Locust

=head1 SYNOPSIS

  USAGE: download_mlst.pl --input_file <ST Output File from Locust>
                          --help

=head1 OPTIONS

B<--input_file, i>   : ST_all.out file from typer.pl run.

B<--type_size, t>			: Number of Alleles in Schema.

B<--attr_size, a>			: Number of Attributes in Schema.

B<--data_type, d>			:Dataset type to make.

B<--help, h>         : Display this help message.

=head1  DESCRIPTION

This program generates an annotation file for the Interactive Tree of Life.

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

my %opts;
GetOptions(\%opts,
	'help|h',
	'input_file|i=s',
	'type_size|t=i',
	'attr_size|a=i',
	'data_type|d=s',
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my @color_palette = ("rgba(102,194,165,1)", "rgba(252,141,98,1)", "rgba(141,160,203,1)", "rgba(231,138,195,1)", "rgba(166,216,84,1)", "rgba(255,217,47,1)", "rgba(229,196,148,1)", "rgba(179,179,179,1)", "rgba(228,26,28,1)", "rgba(55,126,184,1)", "rgba(77,175,74,1)", "rgba(152,78,163,1)", "rgba(255,127,0,1)", "rgba(255,255,51,1)", "rgba(166,86,40,1)", "rgba(247,129,191,1)", "rgba(153,153,153,1)", "rgba(141,211,199,1)", "rgba(255,255,179,1)", "rgba(190,186,218,1)", "rgba(251,128,114,1)", "rgba(128,177,211,1)", "rgba(253,180,98,1)", "rgba(179,222,105,1)", "rgba(252,205,229,1)", "rgba(217,217,217,1)", "rgba(188,128,189,1)", "rgba(204,235,197,1)", "rgba(255,237,111,1)", "rgba(137,197,218,1)", "rgba(218,87,36,1)", "rgba(116,217,68,1)", "rgba(206,80,202,1)", "rgba(63,73,33,1)", "rgba(192,113,124,1)", "rgba(203,213,136,1)", "rgba(95,127,199,1)", "rgba(103,55,112,1)", "rgba(211,217,62,1)", "rgba(56,51,62,1)", "rgba(80,133,120,1)", "rgba(215,193,177,1)", "rgba(104,144,48,1)", "rgba(173,111,59,1)", "rgba(205,155,205,1)", "rgba(209,66,133,1)", "rgba(109,222,136,1)", "rgba(101,41,38,1)", "rgba(127,220,192,1)", "rgba(200,66,72,1)", "rgba(133,105,213,1)", "rgba(94,115,143,1)", "rgba(209,163,61,1)", "rgba(138,124,100,1)", "rgba(89,152,97,1)", "rgba(102,194,165,0.5)", "rgba(252,141,98,0.5)", "rgba(141,160,203,0.5)", "rgba(231,138,195,0.5)", "rgba(166,216,84,0.5)", "rgba(255,217,47,0.5)", "rgba(229,196,148,0.5)", "rgba(179,179,179,0.5)", "rgba(228,26,28,0.5)", "rgba(55,126,184,0.5)", "rgba(77,175,74,0.5)", "rgba(152,78,163,0.5)", "rgba(255,127,0,0.5)", "rgba(255,255,51,0.5)", "rgba(166,86,40,0.5)", "rgba(247,129,191,0.5)", "rgba(153,153,153,0.5)", "rgba(141,211,199,0.5)", "rgba(255,255,179,0.5)", "rgba(190,186,218,0.5)", "rgba(251,128,114,0.5)", "rgba(128,177,211,0.5)", "rgba(253,180,98,0.5)", "rgba(179,222,105,0.5)", "rgba(252,205,229,0.5)", "rgba(217,217,217,0.5)", "rgba(188,128,189,0.5)", "rgba(204,235,197,0.5)", "rgba(255,237,111,0.5)", "rgba(137,197,218,0.5)", "rgba(218,87,36,0.5)", "rgba(116,217,68,0.5)", "rgba(206,80,202,0.5)", "rgba(63,73,33,0.5)", "rgba(192,113,124,0.5)", "rgba(203,213,136,0.5)", "rgba(95,127,199,0.5)", "rgba(103,55,112,0.5)", "rgba(211,217,62,0.5)", "rgba(56,51,62,0.5)", "rgba(80,133,120,0.5)", "rgba(215,193,177,0.5)", "rgba(104,144,48,0.5)", "rgba(173,111,59,0.5)", "rgba(205,155,205,0.5)", "rgba(209,66,133,0.5)", "rgba(109,222,136,0.5)", "rgba(101,41,38,0.5)", "rgba(127,220,192,0.5)", "rgba(200,66,72,0.5)", "rgba(133,105,213,0.5)", "rgba(94,115,143,0.5)", "rgba(209,163,61,0.5)", "rgba(138,124,100,0.5)", "rgba(89,152,97,0.5)", "rgba(102,194,165,0.75)", "rgba(252,141,98,0.75)", "rgba(141,160,203,0.75)", "rgba(231,138,195,0.75)", "rgba(166,216,84,0.75)", "rgba(255,217,47,0.75)", "rgba(229,196,148,0.75)", "rgba(179,179,179,0.75)", "rgba(228,26,28,0.75)", "rgba(55,126,184,0.75)", "rgba(77,175,74,0.75)", "rgba(152,78,163,0.75)", "rgba(255,127,0,0.75)", "rgba(255,255,51,0.75)", "rgba(166,86,40,0.75)", "rgba(247,129,191,0.75)", "rgba(153,153,153,0.75)", "rgba(141,211,199,0.75)", "rgba(255,255,179,0.75)", "rgba(190,186,218,0.75)", "rgba(251,128,114,0.75)", "rgba(128,177,211,0.75)", "rgba(253,180,98,0.75)", "rgba(179,222,105,0.75)", "rgba(252,205,229,0.75)", "rgba(217,217,217,0.75)", "rgba(188,128,189,0.75)", "rgba(204,235,197,0.75)", "rgba(255,237,111,0.75)", "rgba(137,197,218,0.75)", "rgba(218,87,36,0.75)", "rgba(116,217,68,0.75)", "rgba(206,80,202,0.75)", "rgba(63,73,33,0.75)", "rgba(192,113,124,0.75)", "rgba(203,213,136,0.75)", "rgba(95,127,199,0.75)", "rgba(103,55,112,0.75)", "rgba(211,217,62,0.75)", "rgba(56,51,62,0.75)", "rgba(80,133,120,0.75)", "rgba(215,193,177,0.75)", "rgba(104,144,48,0.75)", "rgba(173,111,59,0.75)", "rgba(205,155,205,0.75)", "rgba(209,66,133,0.75)", "rgba(109,222,136,0.75)", "rgba(101,41,38,0.75)", "rgba(127,220,192,0.75)", "rgba(200,66,72,0.75)", "rgba(133,105,213,0.75)", "rgba(94,115,143,0.75)", "rgba(209,163,61,0.75)", "rgba(138,124,100,0.75)", "rgba(89,152,97,0.75)", "rgba(102,194,165,0.25)", "rgba(252,141,98,0.25)", "rgba(141,160,203,0.25)", "rgba(231,138,195,0.25)", "rgba(166,216,84,0.25)", "rgba(255,217,47,0.25)", "rgba(229,196,148,0.25)", "rgba(179,179,179,0.25)", "rgba(228,26,28,0.25)", "rgba(55,126,184,0.25)", "rgba(77,175,74,0.25)", "rgba(152,78,163,0.25)", "rgba(255,127,0,0.25)", "rgba(255,255,51,0.25)", "rgba(166,86,40,0.25)", "rgba(247,129,191,0.25)", "rgba(153,153,153,0.25)", "rgba(141,211,199,0.25)", "rgba(255,255,179,0.25)", "rgba(190,186,218,0.25)", "rgba(251,128,114,0.25)", "rgba(128,177,211,0.25)", "rgba(253,180,98,0.25)", "rgba(179,222,105,0.25)", "rgba(252,205,229,0.25)", "rgba(217,217,217,0.25)", "rgba(188,128,189,0.25)", "rgba(204,235,197,0.25)", "rgba(255,237,111,0.25)", "rgba(137,197,218,0.25)", "rgba(218,87,36,0.25)", "rgba(116,217,68,0.25)", "rgba(206,80,202,0.25)", "rgba(63,73,33,0.25)", "rgba(192,113,124,0.25)", "rgba(203,213,136,0.25)", "rgba(95,127,199,0.25)", "rgba(103,55,112,0.25)", "rgba(211,217,62,0.25)", "rgba(56,51,62,0.25)", "rgba(80,133,120,0.25)", "rgba(215,193,177,0.25)", "rgba(104,144,48,0.25)", "rgba(173,111,59,0.25)", "rgba(205,155,205,0.25)", "rgba(209,66,133,0.25)", "rgba(109,222,136,0.25)", "rgba(101,41,38,0.25)", "rgba(127,220,192,0.25)", "rgba(200,66,72,0.25)", "rgba(133,105,213,0.25)", "rgba(94,115,143,0.25)", "rgba(209,163,61,0.25)", "rgba(138,124,100,0.25)", "rgba(89,152,97,0.25)");


my $input_file = $opts{input_file} || die "Couldn't open the input file.\n";

my $data_type = $opts{data_type};
my $num_of_alleles = $opts{type_size};
my $num_of_attributes = $opts{attr_size};

my @full_st_array;
open(my $fh, "<", $input_file) or die "Failed to open file: $!\n";
while(<$fh>) {
    chomp;
    push @full_st_array, $_;
}
close $fh;

my @columns_to_annotate = &find_columns_to_annotate($num_of_alleles, $num_of_attributes);
#Loop through columns we're going to annotate
foreach my $col (@columns_to_annotate){
	my ($out_hash, $out_file_name) = &create_annotation_hash($col, @full_st_array);
	my $colors_hash = create_color_hash(%$out_hash);
	my ($shape_str, $colors_str, $labels_str) = create_legend_info(%$colors_hash);
	my $out_file = "itol_" . $out_file_name . ".txt";
	write_color_itol_file(\%$out_hash, \%$colors_hash, $shape_str, $colors_str, $labels_str, $out_file);
}

sub find_columns_to_annotate{
	my @columns_to_use;
	#Push ST column
	push (@columns_to_use, 1);
	my ($num_of_alleles, $num_of_attributes) = @_;

	if ($num_of_attributes != 0){
		my $first_column = $num_of_alleles + 2;
		my $max_column = $num_of_alleles + 2 + $num_of_attributes;
		for (my $i = $first_column; $i < $max_column; $i++){
			push (@columns_to_use, $i);
		}
	}
	return @columns_to_use;
}

sub create_annotation_hash{
	my ($column, @row_values) = @_;
	my $header = 0;
	my (%out_hash, $out_file_prefix);
	foreach my $row (@row_values){
		my @split_row = split("\t", $row);
		if ($header == 0){
			$out_file_prefix = $split_row[$column];
			$header = 1;
		} else {
				my $sample = $split_row[0];
				if ( defined $split_row[$column]){
					$out_hash{$sample} = $split_row[$column];
					} else {
						$out_hash{$sample} = "UNKNOWN";
					}
			}
 		}
	return (\%out_hash, $out_file_prefix);
}

sub create_color_hash{
	my %st_hash = @_;
	my %color_hash;
	my $count = 0;
	foreach my $val (values %st_hash){
	if ($color_hash{$val}){
		next;
	}
	else{
		if ($count < scalar @color_palette){
			$color_hash{$val} = $color_palette[$count];
		}
		else {
			$color_hash{$val} = $color_palette[$count % scalar @color_palette];
		}
	}
	$count ++;
	}
	return \%color_hash;
	}

sub create_legend_info{
	my %color_hash = @_;
	my (@legend_shapes, @legend_colors, @legend_labels);
	while(my($k, $v) = each %color_hash){
		push(@legend_shapes, 1);
		push(@legend_colors, $v);
		push(@legend_labels, $k);
	}
 my $legend_shape_str = join("\t", @legend_shapes);
 my $legend_colors_str = join("\t", @legend_colors);
 my $legend_labels_str = join("\t", @legend_labels);

 return ($legend_shape_str, $legend_colors_str, $legend_labels_str);
}

sub write_color_itol_file{
	my ($st_hash, $color_hash, $legend_shape_str, $legend_colors_str, $legend_labels_str, $filename) = @_;
	my %st_hash = %$st_hash;
	my %color_hash = %$color_hash;
	unless ($legend_labels_str eq ""){
		open (my $oh, '>', $filename) or die "Couldn't open the outfile: $filename";
		print $oh "DATASET_COLORSTRIP\n";
		print $oh "SEPERATOR TAB\n";
		print $oh "DATASET_LABEL\tSequence Type\n";
		print $oh "LEGEND_TITLE\tSequence Type\n";
		print $oh "LEGEND_SHAPES\t$legend_shape_str\n";
		print $oh "LEGEND_COLORS\t$legend_colors_str\n";
		print $oh "LEGEND_LABELS\t$legend_labels_str\n";
		print $oh "COLOR\t#ff0000\n";
		print $oh "DATA\n";
		while (my ($key, $value) = each (%st_hash)){
			print $oh "$key\t$color_hash{$value}\t$value\n";
		}
		close $oh;
			if (scalar @color_palette <= keys %color_hash){
				print "Warning: There are less available colors than there are Sequence Types. Colors will be repeated.\n"
			}
		} else {
			open (my $oh, ">", $filename) or die "Couldn't open the outfile: $filename";
			print "No results were found to add to $filename.\n";
			print $oh "No results were found.";
			close $oh;
		}
}

sub write_text_itol_file{
	my ($st_hash, $filename) = @_;
	my %st_hash = %$st_hash;
	open (my $oh, ">", $filename) or die "Couldn't open the outfile: $filename";
	print $oh "DATASET_TEXT\n";
	print $oh "SEPERATOR TAB\n";
	print $oh "DATASET_LABEL\tSequence Type\n";
	print $oh "DATA\n";
	while (my ($key, $value) = each (%st_hash)) {
		print $oh "$key\t$value\t-1\t#000000\tnormal\t1\t0\n";
	}
}
