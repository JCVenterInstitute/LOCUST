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

# =============================== Pragmas ==================================
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use File::Slurp;
use IO::File;
use Cwd;
use Bio::SeqIO;
use Log::Log4perl qw(:easy);
use Data::Dumper;

Log::Log4perl->easy_init({level => $DEBUG, file => ">MLST_ST_finder.log"});

# =============================== Main =====================================

 MAIN:
{
    INFO("INITIALIZING SCRIPT ST_finder.pl");

    my $help;                               # Help message option flag.
    my $identifier = "Sample";              # identifier option.
    my $queryFastaFile;                     # query-fasta option.
    my $mlstAllelesFile;                    # MLST-alleles-fasta option.
    my $stSchemaFile;                       # ST-schema option.
    my $outputName = "./chosen_ST.txt";     # output option.
    my $new_alleles;                        #flag to not overwrite new allele fasta

    my %stMap;                              # hash of Sequence Type number to corresponding allele numbers.
    my %allelesFound;                       # hash of allele name to matching MLST allele.
    my @allelesOrdered;                     # array of allele names in order expected by ST-schema.
    my $query_sequences;                    # hash ref of query sequences

    # =========================== Acquire Parameters =======================
    INFO("Acquiring Parameters.");

    my $result = GetOptions("help"                  => \$help,
			    "identifier=s"          => \$identifier,
			    "query-fasta=s"         => \$queryFastaFile,
			    "MLST-alleles-fasta=s"  => \$mlstAllelesFile,
			    "ST-schema=s"           => \$stSchemaFile,
			    "new_alleles"           => \$new_alleles,
			    "output=s"              => \$outputName);

    # Print help message, if requested, and exit program.
    if ($help) {
	printHelpMsg();
	LOGEXIT("  Printed help message.");
    }

    # Log all parameters.

    if (defined $identifier) {
	INFO("  Using identitifer '$identifier'.");
    }

    checkFile("query-fasta",$queryFastaFile);
    checkFile("MLST-alleles-fasta",$mlstAllelesFile);
    checkFile("ST-schema",$stSchemaFile);

    my $outputDir = dirname($outputName);

    if (defined $outputName) {
	INFO("  Using output filename = '$outputName'.");
	$outputDir = dirname($outputName);
	LOGDIE("  Output directory '$outputDir' is not writable.") if (not (-w $outputDir));
    }

    # =========================== Fill stMap ===============================
    INFO("Filling in Sequence Type map.");

    # Parse through ST-schema file.
    my $stSchemaFile_fh = new IO::File "$stSchemaFile" or LOGDIE("  FATAL: Cannot open ST-schema file '$stSchemaFile'.");

    while (<$stSchemaFile_fh>) {
	chomp $_;
        my @data = split("\t",$_);
	my $ST = shift(@data);       # Remove the first element (Sequence Type) from the line and store it.

	# If this line is the header, store the line's array, which is the ordered list of allele names, in allelesOrdered.
	if ($ST eq "ST") {
	    @allelesOrdered = @data;
	    INFO("  Ordered Allele list: '@allelesOrdered'.")
	}
	# Otherwise, store the ST number in stMap as a value and the line (the ST's corresponding alleles) as a key.
	else {
	    $stMap{join(",",@data)} = $ST;
	}
    }

    $stSchemaFile_fh->close();

    # =========================== Acquire ST of each query =================
    INFO("Acquiring allele for each query and resulting Sequence Type");

    # Open filehandle for output file.
    my $output_fh = new IO::File ">$outputName" or LOGDIE("  FATAL: Cannot open output file '$outputName'.");
    print $output_fh "Sample\tST\t".join("\t",@allelesOrdered)."\n";

    # Open Bio::SeqIO object for queryFastaFile.
    my $inQueries = Bio::SeqIO->new(-file => "<$queryFastaFile", -format => "fasta",);

    #Open File for NEW Allele sequences
    my($file,$new_alleles_fh);

    unless($new_alleles){
	$file = $identifier . "_new.fa";
	$new_alleles_fh = new IO::File ">$outputDir/$file" or LOGDIE("  FATAL: Cannot open output file $file");
    }

    # Parse through the query.
    while (my $query = $inQueries->next_seq) {
	my $queryName = $query->primary_id;                      # Query name.
	my $querySeq = $query->seq;                              # Query sequence.

	#Split name to get queryAllele
	my($queryAllele,$scheme) = split(/\_/,$queryName,2);
	$allelesFound{$queryAllele} = "NEW";

       	$query_sequences->{$queryAllele}=$querySeq;
	INFO("  Examining query '$queryName'.");
  print $querySeq;
	# If the query's sequence begins with "SHORT", declare the queryAllele's hit as SHORT.
	if ($querySeq =~ /SHORT/) {
	    INFO("    Query '$queryName' has been truncated and is SHORT!");
	    $allelesFound{$queryAllele} = "SHORT";
	}
  # If the query's sequence begins with "SHORT", declare the queryAllele's hit as TRUNC.
  if ($querySeq =~ /5'PRTL/){
    INFO("    Query '$queryName''s five prime end is at the end of a contig and is TRUNCated!");
    $allelesFound{$queryAllele} = "5'PRTL";
  }
  if ($querySeq =~ /3'sPRTL/){
    INFO("    Query '$queryName''s three prime end is at the end of a contig and is TRUNCated!");
    $allelesFound{$queryAllele} = "3'PRTL";
  }
	# Otherwise, do the following.
	else {
	    # Open Bio::SeqIO object for mlstAllelesFile.
	    my $inMlstAlleles = Bio::SeqIO->new(-file => "<$mlstAllelesFile", -format => "fasta",);

	    # Parse through the MLST alleles.
	    while (my $mlstAllele = $inMlstAlleles->next_seq) {
		my $mlstAlleleName = $mlstAllele->primary_id;    # MLST Allele Name.
		my $mlstAlleleSeq = $mlstAllele->seq;            # MLST Allele Sequence.

		# If there is a match between the query and MLST allele sequences, do the following.
		if(lc($querySeq) eq lc($mlstAlleleSeq)){
		    # If the query allele is a member of the MLST allele (i.e. the query allele name forms the beginning of the MLST allele name), it is a proper match.
		    if (index($mlstAlleleName,$queryAllele) == 0) {
			INFO("    Query allele '$queryAllele' has found MLST allele '$mlstAlleleName'.");
			$allelesFound{$queryAllele} = ($mlstAlleleName =~ /(\d+)/g)[-1];
		    }
		    # Otherwise, something is wrong.
		    else {
			LOGDIE("    FATAL: Query allele '$queryAllele' has found mismatching MLST allele '$mlstAlleleName' in $identifier!");
		    }

		    last;
		}
	    }
	}

    }

    # stKey - comma-delimited string of allelesFound in order of allelesOrdered.
    my $stKey = "";

    foreach my $allele (@allelesOrdered) {

	#Print new allele sequences to fasta file
	if(exists $allelesFound{$allele}){
	    if($allelesFound{$allele} eq 'NEW'){
		print $new_alleles_fh ">$allele\n" unless $new_alleles;
		print $new_alleles_fh "$query_sequences->{$allele}\n" unless $new_alleles;
	    }
      if(($allelesFound{$allele} eq "5'PRTL") || ($allelesFound{$allele} eq "3'PRTL")){
		print $new_alleles_fh ">$allele\n" unless $new_alleles;
		print $new_alleles_fh "$query_sequences->{$allele}\n" unless $new_alleles;
	    }
	}else{
	    $allelesFound{$allele} = "MISSING";
	}

	$stKey .= "$allelesFound{$allele},";
    }
    chop $stKey;

    INFO("  Sequence Type key: '$stKey'.");

      # If the stKey exists in stMap, print both ST found and stKey to output file.
    if (exists $stMap{$stKey}) {
	INFO("  Found Sequence Type: '$stMap{$stKey}'.");
	print $output_fh "$identifier\t$stMap{$stKey}\t".join("\t",split(",",$stKey))."\n";
    }
    # Otherwise, print "UNKNOWN" and stKey to output file.
    # Note: this will happen if any of the alleles in stKey are SHORT or NEW.
    else {

	print $output_fh "$identifier\tUNKNOWN\t".join("\t",split(",",$stKey))."\n";
    }

    $output_fh->close();


}

# printHelpMsg - prints out help message when 'help' parameter used.
#	Input: None
#	Output: None
sub printHelpMsg {
	my $HELPTEXT = qq~

	MLST_ST_finder.pl [options]

	options:
		--identifier            Unique identifier for the sample (default "Sample").
		--query-fasta           FASTA file of query sequences.
		--MLST-alleles-fasta    FASTA file of MLST allele sequences.
		--ST-schema             Tab-delimited file of Sequence Type and corresponding alleles. Header: "ST<TAB><Allele Name><TAB><Allele Name>...".
		--output                The desired path and name for the output file (i.e. /usr/local/scratch/rsanka/output.txt). Default is "./chosen_ST.txt".
		--help                  Print this help message.
	~;

	print "$HELPTEXT\n";
}

# checkFile - checks if file exists and is readable.
#	Input:
#		- parameter: name of parameter of file.
#		- file: name of file.
#	Output: None
sub checkFile {
	my $parameter = shift;
	my $file = shift;

	if (defined $file) {
		INFO("  Using $parameter = '$file'.");
		LOGDIE("  FATAL: File '$file' for '$parameter' is unreadable.") if (not (-r $file));
	}
	else {
		LOGDIE("  FATAL: You must specify a file for '$parameter'!");
	}
}
