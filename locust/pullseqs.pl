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

pullseqs.pl - generates a multifasta file of loci using a .1con file and a btab file as input.

=head1 USAGE

pullseqs.pl -fastafile <1con file> -blastfile <btab file> -scheme <scheme name>

=head1 OPTIONS

 -fastafile      .1con file for the genome (REQUIRED)
 -blastfile      btab input file with only top hits (REQUIRED)
 -seeds          Allele Seed Fasta File (REQUIRED)
 -max_mismatch   Maximum number of nucleotides that can be not matched on either end of the blast alignment
 -help           display this help message
 -debug          debug mode
 -verbose        lots and lots of output

=head1  DESCRIPTION

pullseqs.pl generates a multifasta file of the loci. The script uses cutFasta to carry out the heavy lifting. The output file is then analyzed by CLC to determine the sequence type (ST) for the strain being analyzed.

=head1  INPUT

User must supply three command line parameters. First, the ".1con" file for the genome. Second, a btab file identifying the candidate sequences in the genome. Third allele seed fasta file.

=head1  OUTPUT

Script outputs a multifasta file of the loci found in the genome.

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin);
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use IO::String;
use lib "$Bin";

############################ OPTIONS ###############################

my ($fastafile, $blastfile, $allele_file, $seeds, $max_mismatch, $help, $DEBUG, $VERBOSE);
my %opts = (
    "fastafile=s" => \$fastafile,
    "blastfile=s" => \$blastfile,
    "alleles=s" => \$allele_file,
    "help" => \$help,
    "debug" => \$DEBUG,
    "max_mismatch=s" => \$max_mismatch,
    "seeds=s" => \$seeds,
    "verbose" => \$VERBOSE,
    );

&GetOptions(%opts);

$max_mismatch = 5 unless($max_mismatch); #sets max number of nucleotides which can be unaligned on either side of the query allele blast match

###################### PROCESS ARGS AND OPTIONS #####################

# GRAB THE PROGRAM NAME
my $program = fileparse($0);

# CHECK FOR HELP REQUEST
if ($help) { &pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDOUT} ); }

# CHECK FOR DEBUG MODE
if ($DEBUG) { $DEBUG = 1; }

# CHECK FOR VERBOSE MODE
if ($DEBUG) {
    $VERBOSE = 1;
    $DEBUG = 1;
}

# VERIFY THAT REQUIRED 1con INPUT FILE SPECIFIED
if (!($fastafile)) { die "Did not supply required 1con file \"-fastafile <1con file>\".\n"; }

# VERIFY THAT REQUIRED btab INPUT FILE SPECIFIED
if (!($blastfile)) { die "Did not supply required btab file \"-blastfile <btab file>\".\n"; }

if (!($seeds)) { die "Did not provide seed fasta file \"--seeds <seed file>\".\n";}

######################### THE NUTS AND BOLTS #########################

my($blastname,$path,$suffix) = fileparse($blastfile);
my ($base, $ext) = split (/\./, $blastname);
my $outputfile = $path . "/" . $base . "_seqs.fa";
my $logfile = $path . "/" . $base . "_pullseqs.log";
my $pepfile = $path . "/" . $base . "_translated_seqs.fa";
if (-e $pepfile){
  unlink($pepfile);
}
open (INFILE, $blastfile) || die "Can't open $blastfile: $!";
open (OUTFILE, ">$outputfile") || die "Can't open $outputfile: $!";
open (LOGFILE, ">$logfile") || die "Can't open $logfile: $!";

while (<INFILE>) {

    chomp $_;
    my $line = $_;

    #print "Processing: $line\n";
    print LOGFILE "Processing: $line\n";

    if($line){
    	my @tokens = split(/\t/, $line);
    	my $querymatchedlength = ($tokens[7] - $tokens[6]) + 1;
    	my $genomematchedlength = ($tokens[8] <= $tokens[9] ? (($tokens[9] - $tokens[8]) + 1) : (($tokens[8] - $tokens[9]) + 1));
    	my $direction = ($tokens[8] <= $tokens[9] ? "forward" : "reverse");
    	my $perc_matched = $querymatchedlength / $tokens[12];
    	my $difference = $tokens[12] - $querymatchedlength;
    	my $short = 0; #true if we cannot retrieve a full length sequence

	#need to get query length form blast because we use the same name for multiple queries and so can't keep the value keyed by query name
	if ($difference < 0) {
	    die "Query match length ($querymatchedlength) from blast hit should not be greater than query length ($tokens[12]).\n";
	}

	#Only extend matches that have max_mismatch or less number of nucleotides that can be not matched on either end of the blast alignment
	if(($tokens[6] <= $max_mismatch) && (($tokens[12] - $tokens[7]) <= $max_mismatch)){

	    #print "Difference: $tokens[0] $tokens[12] - $querymatchedlength = $difference\n";
	    print LOGFILE "Difference: $tokens[0] $tokens[12] - $querymatchedlength = $difference\n";

	    if ($difference != 0) { # The matched region is not "full length", and should be extended

		#print "Original coords: $tokens[8] $tokens[9]\n";
		print LOGFILE "Original coords: $tokens[8] $tokens[9]\n";

		if ($direction eq "forward") {
		    $tokens[8] -= ($tokens[6]-1);
		    $tokens[9] += ($tokens[12] - $tokens[7]);
		    #make sure coordinates do not exceed the bounds of the subject sequence
		    if (($tokens[8] < 1) || ($tokens[9] > $tokens[13])) {
			print LOGFILE "Subject contig too short to pull allele: $tokens[6],$tokens[7];$tokens[12]:$tokens[8],$tokens[9];$tokens[13]\n";
			$short = 1;
		    }
		} elsif ($direction eq "reverse") {
		    $tokens[8] += ($tokens[6]-1);
		    $tokens[9] -= ($tokens[12] - $tokens[7]);
		    #make sure coordinates do not exceed the bounds of the subject sequence
		    if (($tokens[9] < 1) || ($tokens[8] > $tokens[13])) {
			print LOGFILE "Subject contig too short to pull allele: $tokens[6],$tokens[7];$tokens[12]:$tokens[8],$tokens[9];$tokens[13]\n";
			$short = 1;
		    }
		} else {
		    die "Ambiguous match direction: $direction";
		}
		if (!$short) {
		    #print "Adjusted coords: $tokens[8] $tokens[9]\n";
		    print LOGFILE "Adjusted coords: $tokens[8] $tokens[9]\n";
		}
	    }

	}else{

	    print LOGFILE "INFO: Did not extend blast hit. More than $max_mismatch nucleotides unaligned: $tokens[6],$tokens[7];$tokens[12]:$tokens[8],$tokens[9];$tokens[13]\n";
	    $short = 1;

	}

	# cutFasta command requires quotation marks on token[1] to properly parse the sequence identifier from the 'blastfile' *_hits_top.txt
	# e.g., ntkp04:1|cmr:2454  ==>  pipe characters break the script
	# e.g., gi|410113282|emb|CANR01000151.1|  ==> pipe characters break the script

	#Only pull sequences that have max_mismatch or less number of nucleotides that can be not matched on either end of the blast alignment
  if(!$short){
	    my $cutFasta = "$Bin/cutFasta -f $tokens[0] -x $tokens[8] -y $tokens[9] ";
	    $cutFasta .= "-s \"$tokens[1]\" $fastafile";

	    #print "cutFasta command: $cutFasta\n";
	    print LOGFILE "cutFasta command: $cutFasta\n";

	    my @seqlines = `$cutFasta`;
      my $header = $seqlines[0];
      my $fasta_sequence = join("", @seqlines[1 .. $#seqlines]);
      my $fasta_entry = $header . $fasta_sequence;

      open(my $stringfh, "<", \$fasta_entry) or die "Couldn't open sequence for reading: $!";
      my $seqio = Bio::SeqIO->new(
                                  -fh => $stringfh,
                                  -format => "fasta",
                                  );
      #Only 1 sequence
      my $sequence = $seqio->next_seq;
      my $translated_sequence = $sequence->translate(
                                                    -codontable_id => 11,
                                                    );
      my $nuc_length = $sequence->length;
      my $prot_length = $translated_sequence->length;

      my $pred_prot_seq = int($nuc_length / 3);
      if (-e $pepfile){
          open (PEPFILE, ">>", "$pepfile") || die "Can't open $pepfile: $!";
        } else {
          open (PEPFILE, ">", "$pepfile") || die "Can't open $pepfile: $!";
        }
        print PEPFILE $header;
        print PEPFILE $translated_sequence->seq . "\n";
        close(PEPFILE);
      if ($pred_prot_seq == $prot_length){
  	    for my $seq (@seqlines) {
  		      #print $seq;
  	    print OUTFILE $seq;
  	    print LOGFILE $seq;
  	    }
    } else {

        print OUTFILE "$header\n";
        print OUTFILE "PSEUDO\n";
        print LOGFILE "WARN: Printed PSEUDO as sequence because it had a premature stop when translated.\n";

    }
} else {
  if (($tokens[8] == $tokens[13]) && ($perc_matched < 1)){
      print OUTFILE ">$tokens[0]\n";
      print OUTFILE "5'PRTL\n";
      print LOGFILE "WARN: Printed 5'PRTL as sequence because one it was a non full length hit at the end of the contig.\n";
  } elsif (($tokens[9] == $tokens[13]) && ($perc_matched < 1)){
      print OUTFILE ">$tokens[0]\n";
      print OUTFILE "3'PRTL\n";
      print LOGFILE "WARN: Printed 3'PRTL as sequence because one it was a non full length hit at the end of the contig.\n";
    } else {
	    print OUTFILE ">$tokens[0]\n";
	    print OUTFILE "SHORT\n";
	    print LOGFILE "WARN: Printed SHORT as sequence because unaligned nucleotides were greater than $max_mismatch cutoff\n";
	}
}
	#print "\n\n";
	print LOGFILE "\n\n";
    }
}

close (INFILE);
close (OUTFILE);
close (LOGFILE);
