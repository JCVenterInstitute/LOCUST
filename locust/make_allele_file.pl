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

=head1 NAME

make_allele_file.pl 

=head1 SYNOPSIS

  USAGE: make_allele_file.pl --input_file <multifasta file>
                             --help

=head1 OPTIONS

B<--input_file, i>   : Multifasta file

B<--help, h>         : help documentation

=head1  DESCRIPTION

This program takes a multifasta and determines duplicates by headers and sequences. It produces
a non redundant multifasta file with unique headers which can then be used as input as an allele file into the LOCUST program.

=cut
    
use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use Path::Tiny;
use Bio::SeqIO;

my %opts;

GetOptions(\%opts, 'input_file|i=s',
	   'help|h') || die "Error getting options! $!";


pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

&check_params;

#my $fh = path($opts{input_file})->filehandle("<");
my $lfh = path("make_allele_file.log")->filehandle(">");

my ($current_seq, $current_header);
my ($seq_hsh);
my $allele_count;
my ($in_count,$nr_count) = (0,0);

#Loop through fasta file
#Collect header and sequence and store in hash

my $fh = Bio::SeqIO->new(-file => $opts{input_file},
		      -format => 'Fasta');

#Read sequences and store unique ones
#While keeping track of headers of all sequences
while(my $seq = $fh->next_seq){
    $in_count++;

    my $sequence = lc($seq->seq);
    my $header = $seq->display_id;

    if(exists $seq_hsh->{$sequence}){
	
	print $lfh "Duplicate sequence: $header\n";
	print $lfh "Matches:";

	foreach (keys %{$seq_hsh->{$sequence}}){
	    my @s = @{$seq_hsh->{$sequence}->{$_}};
	    print $lfh "$s[0]\n\n";
	    
	}
	
    }else{
	
	my($a,$n) = split(/_/,$header);
	push(@{$seq_hsh->{$sequence}->{$a}},$header);
	$allele_count->{$a} = 1;

    }
}

#Print new multi fasta files with non redundant sequences and
#new numbering
my $name = path($opts{input_file})->basename;
my @v = split(/\./,$name);
my $ofh = path("$v[0].nr")->filehandle(">");

#Go through each sequence to print
foreach my $seq (keys %$seq_hsh){
    
    my $seq_chunk;    
    my @allele_matches = keys %{$seq_hsh->{$seq}};    

    #Print warning of duplicate alleles w/ duplicate sequence
    my $m = "WARNING: " . join(",",@allele_matches) . " from have same sequence" ;
    print_log($m) if (scalar @allele_matches > 1);

    #remove trailing number if one
    foreach my $allele(@allele_matches){

	#Print new header line
	print $ofh ">$allele" . "_" . "$allele_count->{$allele}\n";
	
	#Create fasta format of sequence
	while(my $subset = substr($seq,0,60,"")){
	    $seq_chunk .= "$subset\n";
	}
	
	#Print sequence
	print $ofh "$seq_chunk";

	#increment allele count for new trailing number
	$allele_count->{$allele}++;
	
	#increment count of new entries being printed
	$nr_count++;
    } 
}

#Inform user of final and initial entry counts
my $m = "\n$opts{input_file}\nBefore cleanup: $in_count entries\n";
$m .= "After Cleanup: $nr_count entries";
print_log($m);

exit(0);

sub print_log{

    my $msg = shift;

    print $lfh "$msg\n";
    print STDOUT "$msg\n";
}
sub check_params{
    
    my $error;

    if($opts{input_file}){
	$error = "ERROR: File does not exist or is size zero\n" unless(-s $opts{input_file});
    }else{
	$error = "Usage: ./make_allele_file.pl --input_file <file> \n";
    }
    die($error) if $error;
}


