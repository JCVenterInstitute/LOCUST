#!/usr/bin/env perl
#Copy (C) 2016 The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Pratap Venepally

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

use warnings;
use strict;
#use open ':std', ':encoding(UTF-8)';
$|++;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use Cwd;
use Bio::SeqIO;
use Config::IniFiles;
use FindBin qw($Bin);
use lib "$Bin";

my %opts;

#Arguments/Options
GetOptions( \%opts, 'output|o=s',
	    'input_file|i=s',
	    'type|t=s',
	    'config|c=s',
	    'exclude_genome|e=s',
	    'seed_file|s=s') || die "Error getting options! $!";

my ($MUSCLE_CMD,$TRIM_CMD,$FASTTREE,$RAXML_CMD) = parse_config($opts{config});

unless($opts{seed_file} && $opts{type} && $opts{input_file} && $opts{output}){
        print "USAGE: ./tree_builder.pl --seed_file <seed.file> --input_file <genome.list> --type <raxml|fasttree> -output <output.dir>\n";
	exit;
}

if(!(lc($opts{type}) eq "raxml" or lc($opts{type}) eq "fasttree")){
    die("ERROR: the argument value '$opts{type}' for tree option is not valid. The valid entries are 'raxml' or 'fasttree'");
}
#OPEN LOG FILE
open(my $lfh, ">", "$opts{output}/tree_builder.log");

#make new directory for alleles
my $pwd = cwd();
my $OUTPUT = ($opts{output} =~ /^\//) ? $opts{output} : "$pwd/" . $opts{output};
print $OUTPUT . "\n";
mkdir ("$OUTPUT/alleles") unless (-d "$OUTPUT/alleles");
chdir ("$OUTPUT/alleles") or die "Cannot create or access directory '$OUTPUT/alleles': $!\n";

#grab genome and allele identifiers from input (genomes.txt) and seed (seed.fa) files
my @glines = read_file($opts{input_file});

#Store exclude genomes
my @ex_genomes = read_file($opts{exclude_genome}) if($opts{exclude_genome});
my $ex_hsh;

foreach (@ex_genomes){

    chomp;
    $ex_hsh->{$_} = 1;

}

my %genomes_seen = ();
my @genomes = ();
my $header = 1;

foreach my $gline (@glines){
    chomp $gline;

    my ($genome) = split(/\t/, $gline);

    if(exists $ex_hsh->{$genome}){
	print $lfh "Warning: Skipping $genome because it was included in the exclude genome list for having a MISSING or SHORT allele sequence\n";
    }else{
	#extract genome field

	if ($header) {#check if first line is a header
	    if ($gline =~ /^GENOME\tPATH/) {
		$header = 0;
		next;
	    }
	$header = 0;
	}

	$genome =~ s/\s+$//; #remove trailing spaces

	if(exists $genomes_seen{$genome}){
	    print $lfh "WARNING: a genome is allowed to occur only once in the genome list file - ignoring the second line.\n";
	    print $lfh "$gline\n$genomes_seen{$genome}\n";
	    next;
	}else{
	    $genomes_seen{$genome} = $gline;
	}

	if(-s "$OUTPUT/$genome/$genome" . "_hits.txt"){
	    push (@genomes, $genome);
	}else{
	    print $lfh "WARNING: Skipping $genome because there are no blast hits\n";
	}
    }
}

write_file ("genome_list.txt", map { "$_\n" } @genomes);

my @alines = read_file("$opts{seed_file}");
my %alleles = ();
#print "alleles: @alines\n";

foreach (@alines){
    chomp $_;
    #extract allele label
    if (/^>(\S+)/){
	my ($allele) = split(/\_/, $1, 2);
	$alleles{$allele} = 1;
	#print "alleles: %alleles\n";
    }
}
my @alleles_list = sort((keys %alleles));
write_file ("alleles_list.txt", map { "$_\n" } @alleles_list);

#combine alleles across genomes into individual multifasta files, adding identifiers to headers to keep alleles distinct.
foreach my $genome (@genomes){
    if(-s "$OUTPUT/$genome/$genome"."_hits_top_seqs.fa"){
	my $seqin = Bio::SeqIO->new (-format=>'fasta',-file=> "$OUTPUT/$genome/$genome"."_hits_top_seqs.fa");
	while (my $seqRec = $seqin->next_seq){
	    my ($id)  = split(/\_/, $seqRec->display_id(), 2);
	    $id = $genome."_".$id;
	    $seqRec -> display_id($id);
	    my $aout = Bio::SeqIO->new( -file => ">$id"."_allele.fasta", -format => 'Fasta');
	    $aout->write_seq($seqRec);
	}
    }
}

foreach my $allele (@alleles_list){
    my @afiles;
    my @afasta = glob ("*${allele}*allele.fasta");
    foreach my $afasta (@afasta){
	if(-s $afasta){
	    my $seqin = Bio::SeqIO->new (-format=>'fasta',-file=> $afasta);
	    while (my $seqRec = $seqin->next_seq){
		my $sequence  = $seqRec->seq();
		if ($sequence eq "" or length($sequence) < 10){
		    print $lfh "Allele sequence '$sequence' in '$afasta' is too short or not valid for an object of type ", ref($seqRec), "\n";
		    next;
		}
		else {
		    push(@afiles, $afasta);
		}
	    }
	}
    }
    my $all_afile = $allele."_inAllGenomes.txt";
    write_file($all_afile, map { "$_\n" } @afiles);
}

#Exclude allele.allGenomes.fa from multi-alignment if the allele is missing from any genome:
#count sequences present in the allele file & exclude the file if the number is less than the #of genomes
my @malign;

foreach my $allele (@alleles_list){
    #Add allele file prepend here?
    my $infile = $allele."_inAllGenomes.txt";
    my @allelefiles = read_file ($infile);

    if (scalar(@allelefiles) == scalar(@genomes)){
	push (@malign, $infile);
    }
    else {
	next;
    }
}

#generate multi-fasta for each allele in all genomes
foreach my $file (@malign){
    my ($alleleID) = $file =~ /^(\S+)_.*/;
    my @fasta = read_file ($file);
    foreach my $fasta (@fasta){
	#print "reading fasta file: $fasta\n";
	my $seqin = Bio::SeqIO->new (-format=>'fasta',-file=> $fasta);
	while (my $seqRec = $seqin->next_seq){
	    my $outfile = $alleleID."_AllGenomes.fa";
	    my $mfasta = Bio::SeqIO->new( -file => ">>$outfile", -format => 'Fasta');
	    $mfasta->write_seq($seqRec);
	}
    }
}

#multiple alignment of alleles from genomes: run muscle - iterate over multi-genome allele fasta files
#(slurp fasta from each for input to muscle)
my @malign_alleles = glob("*_AllGenomes.fa");
foreach my $afasta (@malign_alleles){
    my ($outprefix) = $afasta =~ /^(\w+)\.fa/;
    my $outaln = $outprefix.".afa";
    my $cmd = $MUSCLE_CMD . " -in $afasta -out $outaln";
    system("$cmd >/dev/null 2>&1") == 0 || die("ERROR: $cmd failed");
}

#trim multi-alignment fasta to remove gaps.
my @afa_files = glob ("*.afa");
foreach my $afa (@afa_files){
    my ($outprefix) = $afa =~ /^(\w+)\.afa/;
    my $trimaln = $outprefix."_trimmed_alignment.fasta";


    #Now set at .1 (1 - .1 = Remove positions with gaps in 90% or more of the sequences
    my $cmd = $TRIM_CMD . " -in $afa -out $trimaln -gt .1 -fasta";
    system($cmd) == 0 || die("ERROR: $cmd failed");
}

#split allele-specific trimmed alignment fasta files
foreach my $allele (@alleles_list){
    my @a_mfasta = glob ("*$allele*trimmed_alignment.fasta");
    foreach my $afasta (@a_mfasta){
	my $seqin = Bio::SeqIO->new (-format => 'fasta', -file => $afasta);
	while (my $seqRec = $seqin->next_seq){
	    my $id  = $seqRec->display_id();
	    my $afout = Bio::SeqIO->new( -file => ">$id"."_trimmed.fasta", -format => 'Fasta');
	    $afout->write_seq($seqRec);
	}
    }
}

#concatenate trimmed fasta for all alleles from each genome into single seequence
foreach my $genome (@genomes){
    foreach my $allele (@alleles_list){
	my $file = $genome."_".$allele."_trimmed.fasta";
	if (-e $file and -s $file){
	    my $cmd = "cat $file >>$genome"."_all_alleles.fasta";
	    system($cmd) == 0 || die("ERROR: $cmd failed");
	}
	else{
	    print $lfh "File '$file' does not exist\n:$!\n";
	    next;
	}
    }

    #remove individual fasta ID lines to join all records into one string
    if(-s $genome."_all_alleles.fasta"){
	my $mfasta = read_file ($genome."_all_alleles.fasta");
	($mfasta =~ s/>.*\n//g and $mfasta =~ s/\n//g) or die "Cannot remove newlines from '$mfasta':$!\n";
	my $newID = ">".$genome."\n";
	my $newfasta = $newID.$mfasta."\n";
	my $outfile = $genome."_joined_alleles.fasta";
	write_file ($outfile, $newfasta);
    }else{
	print $lfh "File $genome". "_all_alleles.fasta does not exist\n:$!\n";
    }
}

#concatenate joined alleles fasta from all genomes into single multifasta file
my @g_a_fasta = glob ("*joined_alleles.fasta");
foreach my $g_a_fasta (@g_a_fasta){
    my $ofile = "allGenomesJoinedAlleles.fasta";
    my $cmd = "cat $g_a_fasta >>$ofile";
    system($cmd) == 0 || die("ERROR: $cmd failed");
}

#Removing this step as it was said to be unnecessary
#multiple alignment of allGenomesJoinedAlleles.fasta
my $infasta = "allGenomesJoinedAlleles.fasta";
#my ($outprefix) = $infasta =~ /^(\w+)\.fasta/;
#my $outaln = $outprefix.".afa";
#my $cmd = $MUSCLE_CMD . " -in $infasta -out $outaln";
#system("$cmd >/dev/null 2>&1") == 0 || die("ERROR: $cmd failed");

if(-s $infasta){
    if ($opts{type} eq "raxml"){
	#generate ML tree
	my $cmd = $RAXML_CMD;
	$cmd = $cmd . " -p 1234 -f a -x 1234 -N 100 -m GTRGAMMA -n allGenomesJoinedAlleles.tree -s $infasta";
	print "cmd = <$cmd>\n";
	system("$cmd") == 0 || die("ERROR: Error running raxml");
    }

    if($opts{type} eq "fasttree"){
	#generate ML tree
	#my $cmd = $FASTTREE . " -nt -log allGenomesJoinedAlleles_fasttree.log allGenomesJoinedAlleles.afa >allGenomesJoinedAlleles_fasttree.tree";
	my $cmd = $FASTTREE . " -nt -log allGenomesJoinedAlleles_fasttree.log $infasta >allGenomesJoinedAlleles_fasttree.tree";
	system("$cmd") == 0 || die("ERROR: Error runing fasttree");
    }

    #move files and cleanup
    my $cmd = "mv allGenomesJoinedAlleles*fasta *tree* ../ && cd .. && rm -fr alleles";
    system($cmd) == 0 || die("ERROR: $cmd failed");
}else{

    print $lfh "File $infasta does not exist or is size zero. Tree building could not be run\n";
}
exit(0);

sub parse_config{

    my $cfg = Config::IniFiles->new(-file => "$opts{config}");
    my ($muscle,$trimal,$fasttree,$raxml);

    if($cfg->val('muscle','muscle')){
	$muscle = $cfg->val('muscle','muscle');
    }else{
	#Default uses JCVI installation
	$muscle = "/usr/local/bin/muscle";
    }

    if($cfg->val('fasttree','fasttree')){
	$fasttree = $cfg->val('fasttree','fasttree');
	if($fasttree =~ /^\/macOS/){
	    $fasttree = "$Bin/$fasttree";
	}
    }else{
	#Default uses JCVI installation
	$fasttree = "$Bin/FastTreeMP";
    }

    if($cfg->val('raxml','raxml')){
	$raxml = $cfg->val('raxml','raxml');
	if($raxml =~ /^\/macOS/){
	    $raxml = "$Bin/$raxml";
	}

    }else{
	#Default uses JCVI installation
	$raxml = "$Bin/raxmlHPC";
    }

    if($cfg->val('trimal','trimal')){
	$trimal = $cfg->val('trimal','trimal');
	if($trimal =~ /^\/macOS/){
	    $trimal = "$Bin/$trimal";
	}

    }else{
	#Default uses JCVI installation
	$trimal = "$Bin/trimal";
    }

    return ($muscle,$trimal,$fasttree,$raxml);
}
