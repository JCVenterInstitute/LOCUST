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

use strict;
use Cwd;
#Keeps the only top hit for each query

use Data::Dumper;

my $inputfile = $ARGV[0];
my $multi_file = $ARGV[1];
my $multi_flag = $ARGV[2];

my @file_part = split ('/', $inputfile);
my ($genome, $ext) = split ('_hits.', $file_part[-1], 2);

$file_part[-1] = $genome . "_hits_top." . $ext;

my $outputfile = join ('/', (@file_part));

my %hits = ();
my %coords = ();
my %topbitscore = ();
my %toplenratio = ();
my %multi_copy = ();
my %bitscores = ();

open (INFILE, $inputfile) || die "Can't open $inputfile: $!";
open (OUTFILE, ">$outputfile") || die "Can't open $outputfile: $!";
open (my $MULTI_FILE_FH, ">$multi_file") if ($multi_file);

while (<INFILE>) {

	my $line = $_;
	
	next if ($line =~ /^#/);

	my @tokens = split("\t", $line);
	my $qid = _trim($tokens[0]);
	my $sid = _trim($tokens[1]);
	my $percent = _trim($tokens[2]);
	my $qstart = _trim($tokens[6]);
	my $qend = _trim($tokens[7]);
	my $sstart = _trim($tokens[8]);
	my $send = _trim($tokens[9]);
	my $bitscore = _trim($tokens[11]);
	my $qlen = _trim($tokens[12]);
	my $slen = _trim($tokens[13]);
	my $lenratio = ($qend - ($qstart - 1)) / $qlen;

	if ($sstart > $send) {
	    my $tmp = $send;
	    $send = $sstart;
	    $sstart = $tmp;
	}

	my($queryAllele,$scheme) = split(/\_/,$qid,2);

	$bitscores{$queryAllele}{$sid} = $bitscore;
	
	#Store hits above 90% id
	if($percent >= 90){
	    $multi_copy{$queryAllele}{$sid} = _trim($line);

	    if($multi_flag){
		$coords{$qid . ":" . $sid}{'start'} = $sstart;
		$coords{$qid . ":" . $sid}{'end'} = $send;
		$coords{$qid . ":" . $sid}{'sid'} = $sid;
	    }
	}

	if (($lenratio > $toplenratio{$queryAllele}) || 
	    ($lenratio == $toplenratio{$queryAllele})) {
	    
	    if($bitscore > $topbitscore{$queryAllele}){
		
		$toplenratio{$queryAllele} = $lenratio;
		$topbitscore{$queryAllele} = $bitscore;
		
		$coords{$qid . ":" . $sid}{'start'} = $sstart;
		$coords{$qid . ":" . $sid}{'end'} = $send;
		$coords{$qid . ":" . $sid}{'sid'} = $sid;

		#Store the top hit for this queryAllele
		$hits{$queryAllele} = $line;
	    }
	}

}

my @alleles = (keys %coords);
my $overlap_alleles;

for my $i (0 .. $#alleles) {
    my @values1 = split(":",$alleles[$i]);
    my $allele1 = $alleles[$i];
    my $start1 = $coords{$allele1}{'start'};
    my $end1 = $coords{$allele1}{'end'};
    my $sid1 = $coords{$allele1}{'sid'};
    my $length1 = ($end1 - $start1) + 1;

    for my $j (($i + 1) .. $#alleles) {
	my @values2 = split(":",$alleles[$j]);
	my $allele2 = $alleles[$j];
	my $start2 = $coords{$allele2}{'start'};
	my $end2 = $coords{$allele2}{'end'};
	my $sid2 = $coords{$allele2}{'sid'};
	my $length2 = ($end2 - $start2) + 1;
	my $overlap;

	if (($sid1 eq $sid2) && ($start2 < $end1) && ($end2 > $start1)) {#different alleles have overlapping matches in the genome
	    if ($end2 < $end1) {
		if ($start2 < $start1) {
		    $overlap = ($end2 - $start1) + 1;
		} else {
		    $overlap = ($end2 - $start2) + 1;
		}
	    } else {
		if ($start2 < $start1) {
		    $overlap = ($end1 - $start1) + 1;
		} else {
		    $overlap = ($end1 - $start2) + 1;
		}
	    }

	    if (($overlap >= 0.5 * $length1) || ($overlap >= 0.5 * $length2)) {
		my($a1,$scheme1) = split(/\_/,$allele1,2);
		my($a2,$scheme1) = split(/\_/,$allele2,2);

		if($a1 ne $a2){
		    print STDOUT "WARNING: alleles of two different genes ($values1[0] and $values2[0]) have matches that significantly overlap for genome $genome:\n";
		    print STDOUT $hits{$a1};
		    print STDOUT $hits{$a2};

		    #Determine higher bit score and only keep allele that is said to be true
		    my $a1bit = $coords{$allele1}{'bit'};
		    my $a2bit = $coords{$allele2}{'bit'};

		    if($a1bit > $a2bit){
			print STDOUT "INFO: Bit Score for $values1[0] greater than $values2[0]. Treating $values1[0] as true hit\n\n";
			$hits{$a2} = undef;
		    }else{
			print STDOUT "INFO: Bit Score for $values2[0] greater than $values1[0]. Treating $values2[0] as true hit\n\n";
			$hits{$a1} = undef;
		    }
		}else{

		    #Store overlaps of same allele/gene
		    my $a1bit = $coords{$a1}{'bit'};
		    my $a2bit = $coords{$a2}{'bit'};

		    #Remove overlapping multi copy alleles
		    if($a1bit > $a2bit){
			$multi_copy{$a2}{"$allele2:$sid2"} = undef if(exists $multi_copy{$a2}{"$allele2:$sid2"});
		    }else{
			$multi_copy{$a1}{"$allele1:$sid1"} = undef if(exists $multi_copy{$a1}{"$allele1:$sid1"});
		    }

		}
	    }
	}
    }
}

for my $key (keys %hits) {
    my $hit = $hits{$key};
    print OUTFILE $hit;
}

#Sort through same gene hits above 90%
#If they are not marked as an overlap
#Then print them as multi gene
my $print_g = 0;

foreach my $gene(keys %multi_copy){

    my $print_a = 0;

    #my $number = scalar keys $multi_copy{$gene};
    my $number = scalar keys $multi_copy{$gene};
    
    #Now loop and print information if there are multi copy
    if($number > 1){

	foreach my $id(keys $multi_copy{$gene}){
	    
	    my($a,$qid) = split(":",$id);
	    
	    if(defined $multi_copy{$gene}{$id}){
		
		my $hit = $multi_copy{$gene}{$id};
		
		#print header if necessary
		print $MULTI_FILE_FH "#Genome: $genome\n" unless $print_g;
		print $MULTI_FILE_FH "#Allele: $gene\n" unless $print_a;
		$print_a=1;
		$print_g=1;
		
		#print hit
		print $MULTI_FILE_FH "$hit\n";
		
		# IF FLAG PRINT TO HITS FILE
		# Don't print the hit that was the top hit as that was
		# already printed
		if($multi_flag){
		    print OUTFILE "$hit\n" unless (_trim($hits{$gene}) eq $hit);
		}
	    }
	}
    }
}

close (INFILE);
close (OUTFILE);

#subroutines

sub _trim{
    my $value = shift;

    $value =~ s/\s+$//;
    $value =~ s/^\s+//;

    return $value;

}
