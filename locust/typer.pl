#!/usr/bin/env perl

#Copy (C) 2016 The J. Craig Venter Institute (JCVI).  All rights reserved

#This program is free software: you can redistribute it and/or modify
#is under the terms of the GNU General Public License as published by
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
$|++;

=head1 NAME

typer.pl -  A Custom Sequence Locus Typer for Classifying Microbial Genotypic and Phenotypic Attributes.

=head1 SYNOPSIS

  USAGE: typer.pl --input_file <sequence file list>
                         --seed_file <seed sequences>
                         --scheme <scheme>
                         --alleles <allele sequencs>
                         --skip_blast
                         --append_schema
                         --novel_schema
                         --download_schema
                         --tree <type of tree>
                         --config <config file>
                         --org_search
                         --biosample_list
                         --accession_list
                         --genome_type <cg/wgs>
                         --max_contigs
                         --min_N50
                         --previous_output
                         --original_input_file
                         --retype
                         --hmm_model
                         --hmm_cutoff
                         --gb_list
                         --input_path
                         --multi_copy
                         --skip_itol
                         --help

=head1 OPTIONS

B<--input_file, i>   : Input of sequence files (<genome identifier><\t><path to seq file>)

B<--output, o>       : Directory to put output files in[Default: current working directory]

B<--previous_output,p> : Directory to existing result files from a previous typer run.

B<--alleles, a>      : File of allele sequences

B<--seed_file, s>    : File of seed sequences. By default the script creates a seed file from the allele file.

B<--scheme, m>       : scheme

B<--skip_blast>      : Skips blast if the files are already generated

B<--blast_length>    : Number of nucleotides on either end of the blast match that can not match the query allele and still be used [Default: 5]

B<--append_schema>  ; Option to extend schema with new alleles found in initial first run

B<--novel_schema>    : Creates a novel scheme

B<--download_schema>	: Download Schema and Alleles from PubMLST. Provide an organism to search for.

B<--tree, t>         : Type of tree building. Options: raxml,fasttree

B<--exclude_genomes> : Skips including genomes with short/missing alleles when building a tree.

B<--config, c>       : Config file to point to locally installed third party tools

B<--org_search>      : Term to be used to limit NCBI search based on the organism field in the assembly summary file.

B<--genome_type>     : Limit NCBI output to only complete genomes (cg) or incomplete genomes (wgs) based on 'assembly level' in the assembly summary file.

B<--biosample_list>  : Provide a list of biosample IDs to be used as input.

B<--accession_list>  : Provide a list of assembly accession IDs to be used as input.

B<--max_contigs>     : Omit genomes with more than this level of contigs.

B<--min_N50>         :  Minimum N50 for download.

B<--previous_output> : Used with increment. Points to previous LOCUST run output file area.

B<--original_input_file> : Used with increment. Previous LOCUST run input file.

B<--retype>         : Will re-run typer on previous genomes, for incremental mode

B<--hmm_model>      : Models of core hmm to limit results by

B<--hmm_cutoff>      : Percentage of HMMs genomes need to match [1-100]

B<--gb_list>         : List of gb file locations <genome><tab><gb location>

B<--input_path>      : Point to directory of fasta files

B<--multi_copy>      : Flag to pull multi copy genes

B<--skip_itol>			 : Skips making the itol annotation file

B<--help, h>         : Display this help message.

=head1  DESCRIPTION

This program runs the JCVI LOCUST pipeline. A Custom Sequence Locus Typer for Classifying Microbial Genotypic and Phenotypic Attributes.

=head1  INPUT

Instead of pointing to individual input files a user can pass in a path where the script will look for required fasta files. Required file names:

<genome>.fasta

=head1  PUBLICATION

Brinkac LM, Beck E, Inman J, Venepally P, Fouts DE, Sutton G. LOCUST: A Custom Sequence Locus Typer for Classifying Microbial Isolates. Bioinformatics (Oxford, England). 2017 Jan 27;

https://www.ncbi.nlm.nih.gov/pubmed/28130240

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Capture::Tiny qw{ capture capture_merged};
use Pod::Usage;
use File::Slurp;
use Data::Dumper;
use File::Basename;
use File::Path qw(mkpath remove_tree);
use File::Copy qw(move copy);
use File::Touch;
use Path::Tiny;
use Cwd;
use Bio::SeqIO;
use Config::IniFiles;
use FindBin qw($Bin);
use lib "$Bin";
use Clone 'clone';

use Data::Dumper;

my %opts;

#GLOBAL VARS
my $OUTPUT;
my $ORIG_DIR;
my $START_CWD = cwd;
my $CLEAN_FASTA = "perl $Bin/cleanFasta";
my @ORIG_GENOME_LIST;

#Arguments/Options
GetOptions( \%opts, 'input_file|i=s',
	    'seed_file|s=s',
	    'scheme|m=s',
	    'alleles|a=s',
	    'output|o=s',
	    'skip_blast',
	    'blast_length=s',
	    'append_schema',
	    'novel_schema',
            'download_schema=s',
	    'schema_alleles=i',
	    'config|c=s',
	    'tree|t=s',
	    'exclude_genomes',
	    'org_search=s',
	    'genome_type=s',
	    'biosample_list=s',
	    'accession_list=s',
	    'max_contigs=i',
	    'min_N50=i',
	    'previous_output|p=s',
	    'original_input_file=s',
	    'retype',
	    'hmm_model=s',
	    'hmm_cutoff=i',
	    'gb_list=s',
	    'input_path=s',
	    'multi_copy',
	    'skip_itol',
	    'help|h') || die "Error getting options! $!";


pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

#Check input params and set output directory
($OUTPUT,$ORIG_DIR) = &check_params;

#create type log file
my $log_dir  = "$OUTPUT/logs";
mkdir($log_dir) unless (-d $log_dir);

#Download Alleles and Schema from PubMLST
if ($opts{download_schema}){
	my $search_term = $opts{download_schema};
	my $num_of_schema_alleles = $opts{schema_alleles};
	my $cmd = "perl $Bin/download_mlst.pl";
	$cmd .= " -o '$search_term'";
	if ($num_of_schema_alleles){
		$cmd.= " -l $num_of_schema_alleles";
	}
	system($cmd) == 0 || die "\n";
	($opts{alleles}, $opts{scheme}) = get_downloaded_schema_and_alleles($OUTPUT, $START_CWD);
}

my $log_file = "$log_dir/typer.log";
#my $lfh = path($log_file)->filehandle(">");
open(my $lfh, ">", $log_file);

print $lfh "|Started typer.pl "  .  localtime . "\n";

#Open config file and assgin parameters
my($BLAST_CMD,$FORMATDB_EXEC) = parse_config($opts{config});

#Clean user input files to ensure proper formating
&clean_user_input;

#Store original genome list as an array
if($opts{original_input_file}){
    @ORIG_GENOME_LIST = read_file($opts{original_input_file});
}

map{$_ =~ s/\s+$//} @ORIG_GENOME_LIST;

my $combined_input_file;

#Save previous result files
&clean_current_output_dir if($opts{original_input_file});

#create seed file if not provided and requested and not already present
$opts{seed_file} = create_seed($opts{alleles}) unless ($opts{seed_file});
my $SEED_ALLELES = find_seed_alleles($opts{seed_file});

my ($sth,$nth);
my $new_allele_fa = "$OUTPUT/NOVEL_alleles.fa";

#Open files for final output
unless ($opts{novel_schema}) {

    #Open file handle to print all st_finder output to
    $sth = path("$OUTPUT/ST_all.out")->filehandle(">");
    $nth = path($new_allele_fa)->filehandle(">");

}

#Create/set input_file
my ($input_file,$gb_list);
if(($opts{org_search} ||
    $opts{biosample_list} ||
    $opts{accession_list}) &&
   !($opts{input_file})){

    #create input file from downloading from ncbi
    ($input_file,$gb_list) = download_from_ncbi($opts{org_search},
				     $opts{biosample_list},
				     $opts{accession_list});

}elsif($opts{input_path}){

    #if given directory create input file
    $input_file = create_input_file($opts{input_path});

}else{
    #default assignments/options
    $gb_list = $opts{gb_list};
    $input_file = $opts{input_file};
}

#If hmm parsing, run hmm search and modify input file to only
#include genomes that passed hmm search cutoff
my $hmm_passed;

if($opts{hmm_model}){
    $hmm_passed = run_hmm_search($opts{hmm_model}, $opts{hmm_cutoff},$gb_list);
    $input_file = create_hmm_input_file($hmm_passed,$input_file);
}

my $increment_combined_list = &combined_genome_lists($input_file,\@ORIG_GENOME_LIST);

#Run seq type with input file
my($st_files,$fa_files,$top_seqs_files,$multi_copy_logs);

if($opts{retype}){
    ($st_files,$fa_files,$top_seqs_files,$multi_copy_logs) = run_seq_type($increment_combined_list);
}else{
    ($st_files,$fa_files,$top_seqs_files,$multi_copy_logs) = run_seq_type($input_file);
}

#Cat all output files from each genome to into final output
unless($opts{novel_schema}){
    print $lfh "|Step: Combine all genome type files into one\n";
    #Cat all st_files into one

    &_cat($sth,$st_files);
    &_cat($nth,$fa_files);

    close $sth;
    close $nth;
}

#Cat all multi_copy logs
my $mcfh = path("$OUTPUT/multi_copy_hits.txt")->filehandle(">");
&_cat($mcfh,$multi_copy_logs);

#Run st finder on new alleles and add results to original
#output file
if($opts{append_schema}){
    print $lfh "|Step: Creating appended schema\n";
    my $append_outdir = "$OUTPUT/appended_schema";

    #Make append directory for out files
    mkdir($append_outdir);

    #Open combined alleles file
    my $combined_alleles = "$append_outdir/appended_alleles.fa";
    my $cfh = path($combined_alleles)->filehandle(">");

    my $nr_file = &clean_fasta_file($opts{alleles},$new_allele_fa,$append_outdir);

    #Combined new alleles with initial alleles
    my @aa_files;
    push(@aa_files,$nr_file);
    push(@aa_files, $opts{alleles});

    &_cat($cfh,\@aa_files); #make one large combined allele file

    #clean combined fasta file
    my $cmd = "$CLEAN_FASTA $combined_alleles";
    print $lfh "Running: $cmd\n";
    system($cmd) == 0 || die("ERROR: Problem running $cmd");

    unlink($nr_file);
    unlink($combined_alleles . "_orig");
    close $cfh;

    #Run seq type with new alleles
    print $lfh "|Step: Run seq typer on new alleles\n";
    my($nst_files,$nfa_files) = run_seq_type($input_file,$combined_alleles);

    #Cat all new ST files together
    my $new_st_file = "$append_outdir/tmp_allele_ST.out";
    my $nst = path($new_st_file)->filehandle(">");
    &_cat($nst,$nst_files);
    close $nst;

    #Add new ST types to $opts{mlst_scheme}
    print $lfh "|Step: Added new sequene types to user inputed scheme\n";
    &make_new_schema($new_st_file,$opts{scheme},$append_outdir);
    #unlink($new_st_file);

}elsif($opts{novel_schema}){

    print $lfh "|Step: Create novel files\n";
    &create_novel_files($fa_files);

}

#&remove_short_seq_stubs($top_seqs_files);

unless($opts{skip_itol}){
	my $st_file;
	if ($opts{append_schema}){
		$st_file = "$OUTPUT/appended_schema/append_allele_ST.out";
	}	elsif ($opts{novel_schema}){
		$st_file = "$OUTPUT/novel_schema/novel_ST_all.out";
	} else {
		$st_file = "$OUTPUT/ST_all.out";
	}
	&create_itol_file($st_file);
}

if($opts{tree}){
    my $tree_input_file;
    $opts{original_input_file} ? $tree_input_file = $increment_combined_list : $tree_input_file = $input_file;
    &create_tree($opts{tree},$tree_input_file);
    &strain_approximation();
}

&cleanup_files;

print $lfh "|Finished typer.pl "  .  localtime() . "\n";
#############################################
sub create_input_file{

    my $dir = shift;

    my $fh = path("$OUTPUT/directory.list")->filehandle(">");

    if(-d $opts{input_path}){

	my @files = glob("$opts{input_path}/*.fasta");

	foreach(@files){
	    my $file = $_;
	    my $name = path($file)->basename('.fasta');
	    print $fh "$name\t$file\n";
	}

    }else{
	die("ERROR: $opts{input_path} doesn't exist");
    }

    return("$OUTPUT/directory.list");
}
sub find_genome_files{
    my @list;

    if(-d $opts{input_path}){
	my @files = glob("$opts{input_path}/*.fasta");

	foreach(@files){
	    my($name,$path,$suffix) = fileparse($_);
	    my($pre,$post) = split(/\.fasta/,$name);
	    push(@list,"$pre\t$_");
	}
    }else{
	die("ERROR: $opts{input_path} doesn't exist");
    }

    return \@list;
}

sub get_downloaded_schema_and_alleles{
	my $allele_file, my $schema_file = "";
	my($output_dir) = @_;
	my $download_log_file = "$output_dir/logs/download.log";
	print "$download_log_file\n";
	open(my $dlf, "<", $download_log_file) || die "ERROR: Cannot open $download_log_file.\n";
	while(<$dlf>){
		chomp;
		my $line = $_;
		if ($line =~ /^MLST Profile/){
			my ($line_prefix, $downloaded_species) = split /: /, $line;
			chomp $downloaded_species;
			$allele_file = $downloaded_species . "_alleles.fa";
			$schema_file = $downloaded_species . "_schema.txt";
		}
	}
	close($dlf);
	return ($allele_file, $schema_file);
}

sub create_hmm_input_file{
    my($hmm_match,$input_file) = @_;

    my $hfh = path($hmm_match)->filehandle("<");
    my $ifh = path($input_file)->filehandle("<");
    my $nfh = path("$OUTPUT/hmm_parsed_input_file.txt")->filehandle(">");

    my $matched_hsh;

    while(<$hfh>){
	my ($genome,$perc) = split(/\t/,$_);
	$genome =~ s/\s+$//;

	$matched_hsh->{$genome} = $perc;
    }

    while(<$ifh>){
	my ($genome, $location) = split(/\t/,$_);
	$genome =~ s/\s+$//;

	print $nfh $_ if(exists $matched_hsh->{$genome});
    }

    return ("$OUTPUT/hmm_parsed_input_file.txt");

}
sub run_hmm_search{
    my ($model, $cutoff,$gb_list) = @_;

    #First run prep to create necessary pep files
    my $prep_exe = "$Bin/core_hmm_checker_prep.pl";
    my $core_exe = "$Bin/core_hmm_checker.pl";
    my $post_exe = "$Bin/core_hmm_checker_post.pl";

    #Run Prp with gb list
    my $cmd = "$prep_exe --gb_list $gb_list --output $OUTPUT";

    system($cmd) == 0 || die("ERROR: $cmd failed");

    #Run core check
    if(-s "$OUTPUT/pep/pep.list"){
	$cmd = "$core_exe --fasta_list $OUTPUT/pep/pep.list";
	$cmd .= " --cutoff $cutoff --output $OUTPUT";
	$cmd .= " --hmm $model";

	system($cmd) == 0 || die("ERROR: $cmd failed");
    }else{
	die("ERROR: Could not find pep.list file. Check that core_hmm_checker_prep.pl ran\n");
    }

    #Need new input file with only matched genomes
    if(-s "$OUTPUT/hmm/hmm_matched.txt"){
	return("$OUTPUT/hmm/hmm_matched.txt");
    }else{
	die("ERROR: Could not find hmm_matched. Check that core_hmm_checker.pl ran\n");
    }

}
sub find_seed_alleles{

    my $file = shift;
    my $seed_headers = `grep \">\" $file`;

    my @headers = split(">",$seed_headers);
    my $seed_alleles;

    foreach my $a (@headers){
	if($a){
	    $a =~ s/\s+$//;
	    my ($allele,$number) = split(/\_/,$a);
	    $seed_alleles->{$allele} = 1;
	}
    }

    my @seeds = keys %$seed_alleles;
    $seed_alleles = undef;

    return \@seeds;
}

sub combined_genome_lists{
    my ($current,$orig) = @_;
    my $combined_list = "$OUTPUT/increment_combined.list";

    my $fh = path($combined_list)->filehandle(">");
    my $in = path($current)->filehandle("<");

    foreach(<$in>){
	my $line  = $_;
	$line =~ s/\s+$//;

	print $fh "$line\n";
    }

    map{print $fh "$_\n"} @$orig;

    return $combined_list;
}

sub download_from_ncbi{
    print $lfh "|Step: Download sequence files from GenBank\n";

    my ($org_search,$biosample,$accession)  = @_;
    my $exe = "perl $Bin/ftp_download_ncbi_annotation.pl";

    my $base_cmd = "$exe --both --separate_downloads -o $OUTPUT";
    $base_cmd .= " --min_N50 $opts{min_n50}" if $opts{min_n50};
    $base_cmd .= " --max_contigs $opts{max_contigs}" if $opts{max_contigs};

    if($opts{genome_type}){
	if($opts{genome_type} eq 'cg'){
	    $base_cmd .= " --cg";
	}elsif($opts{genome_type} eq 'wgs'){
	    $base_cmd .= " --wgs";
	}
    }

    my (@fasta_files,@gb_files);

    #Find existing assembly file if it's there
    my $assembly_files = glob("$OUTPUT/*assembly_summary.txt");

    if($biosample){
	my $cmd = $base_cmd . " --output_prefix biosample --biosample_list $biosample";

	if($assembly_files){
	    #Use existing assembly file if it's there
	    $cmd .= " --assembly_summary_file $assembly_files";
	}else{
	    #Save assembly file to use later
	    $cmd .= " --preserve_assembly_summary";
	}

	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: Problem running $cmd");

	move("$OUTPUT/fasta.list","$OUTPUT/biosample.list");
	push(@fasta_files,"$OUTPUT/biosample.list");

	move("$OUTPUT/gb.list","$OUTPUT/biosample.gb");
	push(@gb_files, "$OUTPUT/biosample.gb");
    }

    if($accession){
	my $cmd = $base_cmd . " --output_prefix accession --accession_list $accession";

	if($assembly_files){
	    #Use existing assembly file if it's there
	    $cmd .= " --assembly_summary_file $assembly_files";
	}else{
	    #Save assembly file to use later
	    $cmd .= " --preserve_assembly_summary";
	}

	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: Problem running $cmd");

	move("$OUTPUT/fasta.list","$OUTPUT/accession.list");
	push(@fasta_files,"$OUTPUT/accession.list");

	move("$OUTPUT/gb.list", "$OUTPUT/accesion.gb");
	push(@gb_files, "$OUTPUT/accession.gb");
    }

    if($org_search){
	my $cmd = $base_cmd . " --output_prefix org_search --organism_search_term \"$org_search\"";

	if($assembly_files){
	    #Use existing assembly file if it's there
	    $cmd .= " --assembly_summary_file $assembly_files";
	}else{
	    #Save assembly file to use later
	    $cmd .= " --preserve_assembly_summary";
	}

	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: Problem running $cmd");

	move("$OUTPUT/fasta.list","$OUTPUT/org.list");
	push(@fasta_files,"$OUTPUT/org.list");

	move("$OUTPUT/gb.list", "$OUTPUT/org.gb");
	push(@gb_files, "$OUTPUT/org.gb");
    }

    #Cat all fasta list files together in one file
    my $output_file = "$OUTPUT/combined.list";
    my $fh = path($output_file)->filehandle(">");
    &_cat($fh,\@fasta_files);
    $fh = "";

    my $gb_list = "$OUTPUT/combined_gb.list";
    my $gfh = path($gb_list)->filehandle(">");
    &_cat($gfh, \@gb_files);


    if(-s $output_file){
	return ($output_file,$gb_list);
    }else{
	die( "ERROR: Problem download. No combined.list file created");
    }
}

sub clean_current_output_dir{
    #Sub used to clean up any exisiting files in order for incremental run to work
    print $lfh "|Step: Cleaning up output directory before incremental run\n";

    my $prev_dir = "$OUTPUT/previous_run";
    my $orig_dir = $ORIG_DIR;

    #make directory to store previous run files in
    unless(-d $prev_dir){
	print $lfh "|Step:Made $OUTPUT/previous_run directory\n";
	mkdir($prev_dir);
    }

    print $lfh "|Step: Copy existing out files to previous_run directory\n";

    #Move existing combined.list file to another name
    my $combined_list_file = "$prev_dir/original_combined.list";

    if(-s "$orig_dir/combined.list"){
	copy("$orig_dir/combined.list", $combined_list_file);
    }else{
	copy($opts{original_input_file}, "$prev_dir/original_input.list");
    }

    move("$orig_dir/ST_all.out", $prev_dir) if(-s "$orig_dir/ST_all.out");
    move("$orig_dir/NOVEL_alleles.fa", $prev_dir) if(-s "$orig_dir/NOVEL_alleles.fa");
    move("$orig_dir/novel_schema" , "$prev_dir/novel_schema") if (-d "$orig_dir/novel_schema");
    move("$orig_dir/append_schema" , "$prev_dir/append_schema") if (-d "$orig_dir/append_schema");

    my @trees = glob("$orig_dir/*tree*");
    foreach(@trees){
	move($_, "$prev_dir");
    }

    return($combined_list_file) if $combined_list_file;
}
sub clean_user_input{

    print $lfh "|Step: Cleaning user supply input\n";
    my $dos_exect = "$Bin/dos2unix";
    my @files = ($opts{alleles},$opts{biosample_list},$opts{seed},$opts{scheme},$opts{input},$opts{accession_list});

    foreach my $file(@files){

	if($file){
	    my $cmd = $dos_exect . " $file";
	    print $lfh "Running:$cmd\n";
	    system($cmd) == 0 || die( "Error: $cmd failed\n");
	}
    }

    my $cmd = $CLEAN_FASTA . " $opts{alleles}";
    system($cmd) == 0 || die();

    clean_identifiers($opts{input_file}) if ($opts{input_file});
}
sub clean_identifiers{
    my $file = shift;

    my @lines = path($file)->lines;
    my @identifiers;

    foreach(@lines){
	my @values = split(/\t/,$_);

	#Identifiers can only consist of alpha numeric and underscores
	unless($values[0] =~ /^([a-z]|[A-Z]|[0-9]|[_])*$/){
	    push(@identifiers,$values[0]);
	}
}


    if(scalar @identifiers > 0){
	print "ERROR: Identifiers can only contain alpha numeric characters and underscores. Please check the following in $file\n";
	print join("\n",@identifiers);
	print "\n";
	exit();
    }

}
sub parse_config{

    print $lfh "|Step: Parsing config file\n";
    my $cfg = Config::IniFiles->new( -file => "$opts{config}" );
    my ($blastn,$makedb, $cleanfasta);

    if($cfg->val('blast','blastn')){
	$blastn = $cfg->val('blast','blastn');
    }else{
	die("ERROR: $opts{config} did not have a value \"blastn\"\n");
    }

    if($cfg->val('blast','makeblastdb')){
	$makedb = $cfg->val('blast','makeblastdb');
    }else{
	die("ERROR: $opts{config} did not have a value \"makeblastdb\"\n");
    }

    return($blastn,$makedb);
}
sub create_tree{
    print $lfh "|Step: Creating $opts{tree} tree\n";

    my ($tree_type,$input_file) = @_;

    my $cmd = "perl $Bin/tree_builder.pl";
    $cmd .= " -i $input_file";
    $cmd .= " -s $opts{seed_file}";
    $cmd .= " -o $OUTPUT";
    $cmd .= " -c $opts{config}";
    $cmd .= " -t $opts{tree}";

    print $lfh "Running: $cmd\n";
    system($cmd);
}

sub strain_approximation{

	print $lfh "|Step: Generating st approximations.\n";

	my $cmd = "perl $Bin/strain_approximation.pl";
	$cmd .= " -i $opts{input_file}";
	$cmd .= " -s ST_all.out";
	system($cmd) == 0 || die("ERROR: $cmd failed");
}

sub median{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
sub create_seed{

    print $lfh "|Step: Create seed file from allele file\n";

    my $allele_file = shift;
    my $seed_file = "$OUTPUT/seeds_from_alleles.fa";
    my $output_blast = "$OUTPUT/alleles_vs_alleles.table";
    my $logfile = "$OUTPUT/temp_make_blast_db.log";

    if (-s $seed_file && $opts{skip_blast}){
	return $seed_file; #file already exists and --skip_blast means to not recreate it
    }

    if (!(-s $output_blast && $opts{skip_blast})){ #file already exists and --skip_blast means to not recreate it
	#Format fasta file into blast database
	my $cmd = $FORMATDB_EXEC . " -input_type fasta -in $allele_file -dbtype nucl -logfile $logfile";
	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: $cmd failed");
	#Need to customize blast tabular output to include the query length and subject length as the final columns - first twelve columns remain the same
	$cmd = $BLAST_CMD . " -task blastn -db $allele_file -query $allele_file -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -out $output_blast -qcov_hsp_perc 60 -max_target_seqs 2000 -num_threads 4";
	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: $cmd failed");
	unlink($logfile);
	unlink("$allele_file.nhr");
	unlink("$allele_file.nin");
	unlink("

file.nsq");
    }
    open(my $bfh, "<", $output_blast) || die "ERROR: Cannot open $output_blast.\n";
    my %redund_allele_cnt = ();# For each allele count the number of other alleles which make it redundant
    my %covered_allele = ();# For each allele mark whether it is in the set cover yet
    my %redund_graph = ();# The graph of redundancies between alleles
    my %bad_strand = ();# count for each allele of possibly being on the wrong strand
    my %bad_align = ();# count for each allele of possibly badly shifted alignment indicating a mistrimmed allele
    my @bad_allele = ();# array of alleles not to use in the seed because they appear to be on the wrong strand or poorly trimmed
    while(<$bfh>){
	chomp;
	my $line = $_;
	my @tokens = split("\t", $line);
	my $qid = _trim($tokens[0]);
	my $sid = _trim($tokens[1]);
	my $percent = $tokens[2];
	my $qstart = $tokens[6];
	my $qend = $tokens[7];
	my $sstart = $tokens[8];
	my $send = $tokens[9];
	my $qlen = $tokens[12];
	my $slen = $tokens[13];

	(my $q_allele, my $q_id) = split(/\_/, $qid, 2);
	(my $s_allele, my $s_id) = split(/\_/, $sid, 2);
	if (!defined $bad_strand{$q_allele}{$qid}){
	    $bad_strand{$q_allele}{$qid} = 0;
	}
	if (!defined $bad_align{$q_allele}{$qid}){
	    $bad_align{$q_allele}{$qid} = 0;
	}
	if (!defined $bad_strand{$s_allele}{$sid}){
	    $bad_strand{$s_allele}{$sid} = 0;
	}
	if (!defined $bad_align{$s_allele}{$sid}){
	    $bad_align{$s_allele}{$sid} = 0;
	}
	if (!defined $redund_allele_cnt{$qid}){
	    $redund_allele_cnt{$qid} = 1;
	    $covered_allele{$qid} = 0;
	}
	if (!defined $redund_allele_cnt{$sid}){
	    $redund_allele_cnt{$sid} = 1;
	    $covered_allele{$sid} = 0;
	}

	# Check if two alleles are redundant: >= 90% identity and full length match
	if ($qid eq $sid){
	    next;# Skip self matches
	}
	if (defined $redund_graph{$qid}{$sid}){
	    next;# We've already seen and counted this one from the opposite direction
	}
	if ($q_allele ne $s_allele){
	    if (($percent >= 80) && ((($qend - $qstart) >= 0.9 * ($qlen - 1)) || ((($send - $sstart) >= 0.9 * ($slen - 1))))){
		print STDERR "WARNING: two alleles for different genes are very similar: $qid and $sid.\n"
	    }
	    next;# skip matches between different genes
	}
	if ((($qend - ($qstart - 1)) < 0.9 * $qlen) && (abs($send - ($sstart - 1)) < 0.9 * $slen)){
	    next;# skip matches between partial matches
	}
	if ($sstart > $send) {
	    $bad_strand{$q_allele}{$qid}++;
	    #if ($q_allele eq "FUMC") {
		#print STDERR "Bad strand $qid:$sstart,$send,$slen\n";
	    #}
	    if (($percent >= 80) && ((($qstart - 1) != ($slen - $sstart)) || (($qlen - $qend) != ($send - 1)))){
		$bad_align{$q_allele}{$qid}++;
		#if ($q_allele eq "FUMC") {
		    #print STDERR "Bad alignment $qid,$qstart,$qend,$qlen:$sid,$sstart,$send,$slen\n";
		#}
	    }
	} elsif (($percent >= 80) && ((($qstart - 1) != ($sstart - 1)) || (($qlen - $qend) != ($slen - $send)))){
	    $bad_align{$q_allele}{$qid}++;
	    #if ($q_allele eq "FUMC") {
		#print STDERR "Bad alignment $qid,$qstart,$qend,$qlen:$sid,$sstart,$send,$slen\n";
	    #}
	} elsif (($percent >= 90) && (($qend - ($qstart - 1)) == $qlen) && (($send - ($sstart - 1)) == $slen)){
	    $redund_graph{$qid}{$sid} = $redund_graph{$sid}{$qid} = 1; # value doesn't matter - existence defines edge in graph - make symmetric for convenience
	    if (defined $redund_allele_cnt{$qid}){
		$redund_allele_cnt{$qid}++;
	    }else{
		$redund_allele_cnt{$qid} = 2;
		$covered_allele{$qid} = 0;
	    }
	    if (defined $redund_allele_cnt{$sid}){
		$redund_allele_cnt{$sid}++;
	    }else{
		$redund_allele_cnt{$sid} = 2;
		$covered_allele{$sid} = 0;
	    }
	    #if ($q_allele eq "FUMC") {
		#print STDERR "Redund: $qid($redund_allele_cnt{$qid}) $sid($redund_allele_cnt{$sid})\n";
	    #}
	}
    }
    close($bfh);
    #test for wrong strand and mistrimmed alleles
    my @med_array;# temporary storage to compute median
    my $index;
    foreach my $gene (keys %bad_strand){
	@med_array = ();# clear temporary storage to compute median
	$index = 0;
	#print STDERR "Bad strand counts for $gene:";
	foreach my $count (values $bad_strand{$gene}){
	    $med_array[$index++] = $count;
	    #print STDERR "\t$count";
	}
	my $bad_strand_median_cut = (2 * median(@med_array)) + 1;
	#print STDERR "\ncutoff: $bad_strand_median_cut\n";
	foreach my $allele (keys $bad_strand{$gene}){
	    #print STDERR "$allele: $bad_strand{$gene}{$allele}\n";
	    if (($bad_strand{$gene}{$allele} > $bad_strand_median_cut) && ($bad_strand{$gene}{$allele} > $redund_allele_cnt{$allele})) {
		push(@bad_allele,$allele);
		print STDERR "WARNING: $allele with bad strand count $bad_strand{$gene}{$allele} appears to be on the opposite strand from most of the alleles for the $gene gene. Not using as a seed.\n";
	    }
	}
    }
    foreach my $gene (keys %bad_align){
	@med_array = ();# clear temporary storage to compute median
	$index = 0;
	#print STDERR "Bad trim counts for $gene";
	foreach my $count (values $bad_align{$gene}){
	    $med_array[$index++] = $count;
	    #print STDERR "\t$count";
	}
	my $bad_align_median_cut = (2 * median(@med_array)) + 1;
	#print STDERR "\ncutoff: $bad_align_median_cut\n";
 	foreach my $allele (keys $bad_align{$gene}){
	    #print STDERR "$allele: $bad_align{$gene}{$allele}\n";
	    if (($bad_align{$gene}{$allele} > $bad_align_median_cut) && ($bad_align{$gene}{$allele} > $redund_allele_cnt{$allele})) {
		push(@bad_allele,$allele);
		print STDERR "WARNING: $allele with bad trim count $bad_align{$gene}{$allele} appears to be mistrimmed based on alignments to other alleles for the $gene gene. Not using as a seed.\n";
	    }
	}
    }
    while (scalar @bad_allele) { #iterate over removing bad alleles based on redundancies
	foreach my $allele (@bad_allele){# do not use bad alleles in the seed file
	    delete $redund_allele_cnt{$allele};
	    delete $redund_graph{$allele};
	}
	@bad_allele = (); #reset array
	foreach my $allele (keys %redund_allele_cnt){# check for consistency of bad alleles
	    if (defined $redund_graph{$allele}){
		foreach my $redund (keys $redund_graph{$allele}){
		    if (!defined $redund_allele_cnt{$redund}){
			print STDERR "WARNING: $redund was determined to be bad but was redundant with $allele which was determined to be good which is a conflict! Not using as a seed.\n";
			#remove this allele as well
			push(@bad_allele,$allele);
		    }
		}
	    }
	}
    }
    my %alleles_for_seed = ();# The set of alleles to keep for the seed file
    my @sorted_alleles = (sort {$redund_allele_cnt{$b} <=> $redund_allele_cnt{$a}} keys %redund_allele_cnt);# Sort from largest to smallest because we want to remove the least redundant alleles in favor of the more redundant which we believe to be more central while leaving a "set cover"
    while (scalar @sorted_alleles) { #iterate over choosing largest subset and removing redundant alleles
	my $allele = $sorted_alleles[0]; #largest subset
	#print STDERR "allele: $allele count: $redund_allele_cnt{$allele} $covered_allele{$allele}\n";
	if ($redund_allele_cnt{$allele} == 0){
	    last;#All remaining subsets have been completely covered
	}
	if (defined $redund_graph{$allele}) {
	    foreach my $redund (keys $redund_graph{$allele}){# Reduce allele counts for intersecting subsets
		if(!$covered_allele{$redund}) {
		    $covered_allele{$redund} = 1;# Mark as covered
		    if(!$covered_allele{$allele}) {
			$redund_allele_cnt{$redund} -= 2;
		    } else {
			$redund_allele_cnt{$redund}--;
		    }
		    #print STDERR "redund: $redund count: $redund_allele_cnt{$redund} $covered_allele{$redund}\n";
		    foreach my $covered (keys $redund_graph{$redund}){
			$redund_allele_cnt{$covered}--;
			#print STDERR "covered: $covered count: $redund_allele_cnt{$covered} $covered_allele{$covered}\n";
		    }
		} else {
		    if(!$covered_allele{$allele}) {
			$redund_allele_cnt{$redund}--;
			#print STDERR "redund2: $redund count: $redund_allele_cnt{$redund} $covered_allele{$redund}\n";
		    }
		}
	    }
	}
	$alleles_for_seed{$allele} = 1;# Add the current allele for the largest remaining subset to the set cover
	if(!$covered_allele{$allele}) {
	    $covered_allele{$allele} = 1;# Mark as covered
	    $redund_allele_cnt{$allele}--;
	}
	#print STDERR "allele2: $allele count: $redund_allele_cnt{$allele} $covered_allele{$allele}\n";
	#Remove the allele we added to the seeds file
	delete $redund_allele_cnt{$allele};
	@sorted_alleles = (sort {$redund_allele_cnt{$b} <=> $redund_allele_cnt{$a}} keys %redund_allele_cnt);# Sort from largest to smallest because we want to remove the least redundant alleles in favor of the more redundant which we believe to be more central while leaving a "set cover"
    }
    #foreach my $allele (keys %alleles_for_seed){
	#print "$allele:\n";
    #}
    open(my $afh, "<", $allele_file) || die "ERROR: Cannot open $allele_file.\n";
    open(my $sfh, ">", $seed_file) || die "ERROR: Cannot open $seed_file.\n";
    my $print_allele = 0;
    while(<$afh>){
	my $line = $_;
	if($line =~ /^>/){
	    my $allele = $line;
	    #print "$allele";
	    $allele =~ s/^>//;
	    $allele =~ s/\s.*//;
	    $allele =~ s/\s+//;
	    #print "$allele:";
	    if(defined $alleles_for_seed{$allele}){
		#print "1\n";
		$print_allele = 1;
	    }else{
		#print "0\n";
		$print_allele = 0;
	    }
	}
	if($print_allele){
	    print $sfh $line;
	}
    }
    close($afh);
    close($sfh);
    return $seed_file;
}
sub remove_short_seq_stubs{

    my $top_seqs_files = shift;

    foreach my $file (@$top_seqs_files){
	open(my $fh, "<", $file) || die "ERROR: Cannot open $file.\n";
	open(my $ofh, ">", "$file.cleaned") || die "ERROR: Cannot open $file.cleaned.\n";
	my $first;
	my $header;
	while(<$fh>){
	    my $line = $_;
	    if($line =~ /^>/){
		$header = $line;
		$first = 1;
		next;
	    }
	    if($line =~ /^SHORT$/ || $line =~ /^PSEUDO$/ || $line =~ /^5'PRTL$/ || $line =~ /^3'PRTL$/){
		next;
	    }
	    if($first){
		print $ofh $header;
		$first = 0;
	    }
	    print $ofh $line;
	}
	close $fh;
	close $ofh;
	unlink($file);
	link("$file.cleaned",$file);
	unlink("$file.cleaned");
    }

    return;
}
sub create_novel_files{

    my $fasta_files = shift;

    my $unique_sequences; #stores unique allele sequences
    my $genome_st; #stores what alleles are associated with what genome
    my $allele_number; #stores allele numbers to increment
    my @genome_names; #stores genome names

    #loop through each genomes top hit sequence file
    foreach my $file (@$fasta_files){

	#find the genome name to use in novel ST file
	my ($name,$path,$suffix) = fileparse($file);
	my $genome = $1 if ($name =~ /(.*)\_hits\_top\_seqs\.fa/);
	push(@genome_names,$genome);

	#open genome fasta file
	my $fh = path($file)->filehandle("<");
	my ($c_allele,$p_allele,$current_sequence) = ("","","");
	
	#loop through fasta file and store allele information
	#and sequence if it is unique. Determine correct allele
	#designation to assign

	my $unique_id = 0;

	while(<$fh>){

	    my $line = $_;
	    $line =~ s/\s+$//;
	    
	    if($line =~ /^>/){

		#Parse defline to determine allele name and reconfigure to correct format
		my ($h,$temp_allele) = split(/>/,$line);
		my @alleles = split(/\_/,$temp_allele, 2); #removes genome designation if appended to allele name

		$c_allele = $alleles[0];
		$c_allele =~ s/\s+$//;

		$unique_id++ if($current_sequence);

		$c_allele = $c_allele . ":" . $unique_id;
		$p_allele = $c_allele unless($p_allele);
		
		if($current_sequence){

		    #Remove all whitespace to make one string
		    $current_sequence =~ s/\s+//g;

		    my @values = split(":",$p_allele);
		    
		    if ($current_sequence  =~ /^SHORT/) {

			#Store allele SHORT for this genome
			$genome_st->{$genome}->{$values[0]}->{$p_allele} = "SHORT";
			
		    } elsif ($current_sequence =~ /^3'PRTL/){
			
			$genome_st->{$genome}->{$values[0]}->{$p_allele} = "3'PRTL";
			
		    } elsif ($current_sequence =~ /^5'PRTL/){
			
			$genome_st->{$genome}->{$values[0]}->{$p_allele} = "5'PRTL"; 
			
		    } elsif ($current_sequence =~ /^PSEUDO/){
			$genome_st->{$genome}->{$values[0]}->{$p_allele} = "PSEUDO";
			
		    } else {

			my $novel_id;

			if(exists $unique_sequences->{$current_sequence}->{$values[0]}){
			  
			    my($id,$orig_allele) = split("_",$unique_sequences->{$current_sequence}->{$values[0]});
			    $genome_st->{$genome}->{$values[0]}->{$orig_allele} = $id;
			    
			}else{
			
			    #Increment number based on original allele name, not unique value of allele and unique id
			    if(exists $allele_number->{$values[0]}){
				$allele_number->{$values[0]}++;
			    }else{
				$allele_number->{$values[0]} = 1;
			    }

			    #Store new seq and assign next incremented allele number
			    $novel_id = "NOVEL" . $allele_number->{$values[0]};
			    $unique_sequences->{$current_sequence}->{$values[0]} = $novel_id . "_$p_allele";

			    $genome_st->{$genome}->{$values[0]}->{$values[0]. ":" .$unique_id} = $novel_id;
			}

		    }

		    #reset current sequence and set current allele to previous allele
		    $current_sequence = "";
		    $p_allele = $c_allele;

		}
	    }else{

		#collect all the sequence lines in one string
		$current_sequence .= $line;

	    }

	}
	
	#Account for final sequence
	if($current_sequence){

	    $unique_id++;
	    
	    #Remove all whitespace to make one string
	    $current_sequence =~ s/\s+//g;

	    #split just in case multi copy
	    my @values = split(":",$p_allele);
	    
	    if ($current_sequence  =~ /^SHORT/) {

		#Store allele SHORT for this genome
		$genome_st->{$genome}->{$values[0]}->{$p_allele} = "SHORT";

	    }elsif ($current_sequence =~ /^3'PRTL/){

		$genome_st->{$genome}->{$values[0]}->{$p_allele} = "3'PRTL";
		
	    } elsif ($current_sequence =~ /^5'PRTL/){
		
		$genome_st->{$genome}->{$values[0]}->{$p_allele} = "5'PRTL";
		
	    } elsif ($current_sequence =~ /^PSEUDO/){
	
		$genome_st->{$genome}->{$values[0]}->{$p_allele} = "PSEUDO";
		
	    }else{

		my $novel_id;


		if(exists $unique_sequences->{$current_sequence}->{$values[0]}){

		    my($id,$orig_allele) = split("_",$unique_sequences->{$current_sequence}->{$values[0]});
		    $genome_st->{$genome}->{$values[0]}->{$orig_allele} = $id;
		
		}else{
		
		    #Split p_allele name on : which is used to help multi copy designation
		    my @values = split(":",$p_allele);

		    #Increment number based on original allele name, not unique value of allele and unique id
		    if(exists $allele_number->{$values[0]}){
			$allele_number->{$values[0]}++;
		    }else{
			$allele_number->{$values[0]} = 1;
		    }

		    #Store new seq and assign next incremented allele number
		    $novel_id = "NOVEL" . $allele_number->{$values[0]};
		    $unique_sequences->{$current_sequence}->{$values[0]} = $novel_id . "_$p_allele";

		    $genome_st->{$genome}->{$values[0]}->{$values[0]. ":" .$unique_id} = $novel_id;
		}
				
	    }
	}
    }
    #Modify genome_st to add _MC label where necessary
    if($opts{multi_copy}){
	foreach my $genome(keys %$genome_st){
	    
	    foreach my $allele(keys %{$genome_st->{$genome}}){
		
		my $count_alleles = scalar(keys %{$genome_st->{$genome}->{$allele}});
		my $labels;
		
		if($count_alleles > 1){ #means multi copy
		    
		    foreach my $unique(keys %{$genome_st->{$genome}->{$allele}}){
			
			my $label = $genome_st->{$genome}->{$allele}->{$unique};
			
			if(exists $labels->{$label}){
			    
			    $genome_st->{$genome}->{$allele}->{$unique} = $label . "_MC";
			    
			    my $orig_label = $labels->{$label};
			    $genome_st->{$genome}->{$allele}->{$orig_label} = $label . "_MC";
			    
			}else{
			    
			    $labels->{$label} = $unique;
			}
			
		    }
		    
		}				       
		
	    }
	
	}
    }
    
    #Sort allele and genome names to ensure proper ordering
    my @alleles = sort keys %$allele_number;
    @genome_names = sort @genome_names;

    #Create novel outdir
    my $novel_outdir = "$OUTPUT/novel_schema";
    mkdir($novel_outdir) unless (-d $novel_outdir);


    #Call novel print subs
    print_novel_schema($genome_st,\@genome_names,\@alleles,$novel_outdir);
    print_novel_fasta($unique_sequences,$novel_outdir);
}

sub print_novel_fasta{
    my ($sequences,$outdir) = @_;
    my $file = "$outdir/novel_alleles.fa";

    my $fh = path($file)->filehandle(">");

    foreach my $seq (keys %$sequences){
	
	foreach my $allele (keys %{$sequences->{$seq}}){
	    
	    if ($seq !~ /^SHORT/) {
		my @values = split(":",$allele);
		my ($novel_id,$orig_allele) = split("_",$sequences->{$seq}->{$allele});

		my $header = "$values[0]" . "_" . "$novel_id";
		print $fh ">$header\n";
		print $fh "$seq\n";
	    }
	}
    }

    #clean combined fasta file
    if(-s $file){
	my $cmd = $CLEAN_FASTA . " $file";
	system($cmd) == 0 || die("ERROR: Problem running $cmd");

	unlink($file . "_orig");
    }else{
	print $lfh "|WARNING: The novel_alleles.fa file is of size zero\n";
    }
}
sub print_novel_schema{

    my ($genome_st,$genome_names,$alleles,$outdir) = @_;
    my $ST;
    my $ST_number = 0;
    my $genome_alleles;
    my $genome_print;

    #Determine which alleles have multiple variants
    foreach my $genome(@$genome_names){
	
	my %locations;
	my $location_count = 0;
	my @multi_variants;
	my @base_string;
	
	foreach my $allele (sort @$SEED_ALLELES){
	    
	    $locations{$allele} = $location_count;
	    $location_count++;
	    
	    my @v;

	    if(exists $genome_st->{$genome}->{$allele}){
		@v = keys $genome_st->{$genome}->{$allele};
	    }else{
		$genome_st->{$genome}->{$allele}->{'1'} = "MISSING";
	    }

	    my $count = scalar(@v);

	    if($count > 1){ #Means multi copy for this allele
	
		push(@multi_variants,\@v);
		push(@base_string,"");

	    }else{

		if(@v){
      		    my $schema_num = $genome_st->{$genome}->{$allele}->{$v[0]};
		    push(@base_string,$schema_num);
		    
		}else{
		    push(@base_string,"MISSING");
		}
	    }
	}

	my @variants;
	push @variants, clone(\@base_string);
	
	foreach my $gene(@multi_variants){
	  
	    my($base_allele,$unique_num) = split(/:/,$gene->[0]);

	    my $values_to_add = scalar(@variants);
	    my $size = scalar(@$gene);
	    grow_variants( \@variants, $size );

	    my $allele_index = $locations{$base_allele};
	   
	    my $values_added = 0;
	    my $index = 0;  #increment this until it's time to move on.

	    for my $allele ( @$gene ) {
	    
		my $allele_schema = $genome_st->{$genome}->{$base_allele}->{$allele};

		while ( $values_added < $values_to_add ) {
		    # Fill the location for allele at this entry
		    $variants[$index][$allele_index] = $allele_schema;
		    
		    # then increment both index and values_added
		    $index++;
		    $values_added++;
		}

		$values_added = 0; # reset values added between alleles.
	    }
	}
	
	#Clean up @variants to remove duplicate type string
	my $unique_variants;

	foreach (@variants){
	    my $variant_value = $_;
	    my $variant_string = join("\t",@$variant_value);
	    
	    $unique_variants->{$variant_string} = 1;

	    unless($variant_string =~ /(MISSING|SHORT|PSEUDO|PRTL)/){
		$ST_number++;
		$ST->{$variant_string} = $ST_number;
	    }
	}

	$genome_print->{$genome} = $unique_variants;
    }

    #Parse genome st hash to come up with numbers assoicated
    #with allele combinations
    my $st_fh = path("$outdir/novel_ST_all.out")->filehandle(">");
    my $s_fh = path("$outdir/novel_schema.txt")->filehandle(">");

    #Print headers
    print $st_fh "Sample\tST\t" . join("\t", sort @$SEED_ALLELES) . "\n";
    print $s_fh "ST\t" . join("\t", sort @$SEED_ALLELES). "\n";
	   
    #Print ST
    foreach my $genome(keys %$genome_print){

	foreach my $v (keys $genome_print->{$genome}){
	    print $st_fh $genome . "\t";

	    if($v =~ /(MISSING|SHORT|PSEUDO|PRTL)/){
		print $st_fh "UNKNOWN\t";
	    }else{
		print $st_fh $ST->{$v} . "\t";
	    }
	   
	    print $st_fh $v;
	    print $st_fh "\n";

	}
    }
    
    #Print Schema
    foreach my $ST_key (sort {$ST->{$a} <=> $ST->{$b}} keys %$ST){
	print $s_fh $ST->{$ST_key} . "\t" . $ST_key . "\n";
    }
}
sub make_new_schema{
    my ($new_st_file, $orig_st, $outdir) = @_;

    my $max_st_num; #stores current max number
    my $schema; #stores all the new/orig ST types
    my $stAttributes; #store any columns that come after the schema

    my $new_file =  "$outdir/appended_scheme.txt";
    my $final_st_file = "$outdir/append_allele_ST.out";

    open(my $new_fh, ">", $new_file) || die "ERROR: Cannot open $new_file.\n";
    open(my $final_st_fh, ">", $final_st_file) || die "ERROR: Cannot open $final_st_file.\n";
    open(my $orig_fh, "<", $orig_st) || die "ERROR: Cannot open $orig_st.\n";

    #Determines the max ST Type Number found in
    #original file
    my ($ST_type_size,$ST_attr_size);
    my $header =1;
    my $unknown = "";

    while(my $line = <$orig_fh>){
	$line =~ s/\s+$//;

        my @data = split("\t",$line);
	my $ST = shift(@data);       # Remove the first element (Sequence Type) from the line and store it.

	# If this line is the header, store the line's array, which is the ordered list of allele names, in allelesOrdered.
	if ($header) {
	    $header = 0;
	    if ($ST eq "ST") {
		$ST_type_size = scalar @data;
		$ST_attr_size = 0;
	    } elsif ($ST =~ /^ST\((\d+)\)$/) {
		$ST_type_size = $1;
		if ($ST_type_size > (scalar @data)) {
		    die "ERROR: fewer header columns than specified in $ST:\n$line\n";
		}
		$ST_attr_size = (scalar @data) - $ST_type_size;
		if ($ST_attr_size > 0) {
		    foreach ($ST_type_size .. $#data) {
			$unknown .= "\t";
		    }
		}
	    } else {
		die "ERROR: header line for ST schema must begin with ST tab or ST(#) tab - not:\n$line\n";
	    }
	    print $new_fh "$line\n";
	}
	# Otherwise, store the ST number in schema as a value and the line (the ST's corresponding alleles) as a key.
	else {

	    if($max_st_num){
		$max_st_num = $ST if($ST > $max_st_num);
	    }else{
		$max_st_num = $ST;
	    }

	    my $size = scalar @data;
	    if (($size != ($ST_type_size + $ST_attr_size)) && ($size != $ST_type_size)) {
		die "ERROR: ST $ST does not have the same number of columns as the header line or the alleles only: $size versus $ST_type_size alleles and $ST_attr_size attributes.\n";
	    }

	    my $st_type = join("\t",@data[0 .. ($ST_type_size - 1)]);

	    unless(exists $schema->{$st_type}){
		$schema->{$st_type} = $ST;

		if ($ST_attr_size > 0) {

		    if ($size != $ST_type_size) {
			$stAttributes->{$ST} = join("\t", @data[$ST_type_size .. $#data]);
		    } else {
			$stAttributes->{$ST} = $unknown;
		    }
		}
	    } else {
		print $lfh "|WARNING: duplicated ST entry: ($st_type) in $orig_st with ST values $schema->{$st_type} and $ST: $schema->{$st_type} will be used and $ST ignored,\n";
		print STDERR "WARNING: duplicated ST entry: ($st_type) in $orig_st with ST values $schema->{$st_type} and $ST: $schema->{$st_type} will be used and $ST ignored,\n";
	    }
	}
    }


    #Parses new ST types and stores in hash
    open(my $new_st_fh, "<", $new_st_file);

    while(<$new_st_fh>){

	unless($_ =~ /Sample/){

	    chomp $_;

	    my @values = split(/\t/,$_);
	    (my $sample, my $st_num) = splice(@values,0,2);
	    my $st_type = join("\t",@values[0 .. ($ST_type_size - 1)]);

	    if ($st_num eq "UNKNOWN") {

		my $skip_bad = 0;

		foreach my $value (@values) {
		    if (($value eq "MISSING") || ($value eq "SHORT") || ($value eq "3'PRTL") || ($value eq "5'PRTL") || ($value eq "PSEUDO")) {
			$skip_bad = 1;
		    } elsif ($value eq "NEW") {
			$skip_bad = 1;

			print $lfh "|WARNING: no NEW alleles should be found for append_schema option but $sample had type $st_num for ($st_type).\n";
			print STDERR "WARNING: no NEW alleles should be found for append_schema option but $sample had type $st_num for ($st_type).\n";

		    }

		}

		next if ($skip_bad);

		unless(exists $schema->{$st_type}){
		    $max_st_num++;
		    $schema->{$st_type} = $max_st_num;
		    $stAttributes->{$max_st_num} = $unknown;
		}

	    } else {

		if ((!defined $schema->{$st_type}) || ($st_num != $schema->{$st_type})) {
		    die "ERROR: $st_num assigned for ($st_type) not found in $orig_st.\n";
		}

	    }

	}

    }

    #Prints the new ST stypes to the new appended scheme
    foreach my $st_type (sort {$schema->{$a} <=> $schema->{$b}} keys %$schema){

	print $new_fh "$schema->{$st_type}\t$st_type";

	if ($ST_attr_size > 0) {
	    print $new_fh "\t$stAttributes->{$schema->{$st_type}}";
	}

	print $new_fh "\n";
    }

    close $new_fh;

    #Write new ST_out w/ new schema numbers
    open ($new_st_fh, "<", $new_st_file);

    while(<$new_st_fh>){

	my $line = $_;
	$line =~ s/\s+$//;

	if($line =~ /UNKNOWN/){

	    my @values  = split(/\t/,$line);

	    (my $sample, my $st_num) = splice(@values,0,2);
	    my $st_type = join("\t",@values[0 .. ($ST_type_size - 1)]);
	    my $st_output = join("\t",@values);

	    if (defined $schema->{$st_type}) {
		$st_num = $schema->{$st_type};
	    }

	    print $final_st_fh "$sample\t$st_num\t$st_output\n";

	}else{

	    print $final_st_fh "$line\n";

	}

    }

    close $new_st_fh;
}

sub run_seq_type{
    print $lfh "|Step: run sequence typer\n";

    my ($input_file,$new_alleles) = @_;
    my %genomes_seen = ();

    # initialize ST and allele maps for ST typing later
    my %stMap;                              # hash of Sequence Type number to corresponding allele numbers.
    my %stAttributes;                       # hash of Sequence Type number to corresponding columns of attributes to be output.
    my %alleleMap;                          # hash of allele sequences to corresponding allele name
    my @allelesOrdered;                     # array of allele names in order expected by ST-schema.
    my $append = (defined $new_alleles) ? 1 : 0;

    &init_st_finder($opts{scheme}, (defined $new_alleles) ? $new_alleles : $opts{alleles}, \%stMap, \%alleleMap, \@allelesOrdered, \%stAttributes, $append) unless($opts{novel_schema});

    #Open input file and loop through genomes
    my @lines = read_file($input_file);
    my (@st_files,@fa_files,@top_seqs_files,@multi_copy);
    my $genomeHeader = "";
    my $header = 1;

    my $tfh = path("$OUTPUT/logs/tophits.log")->filehandle(">");;

    foreach my $line (@lines){

	chomp $line;

	if ($header) {#check if first line is a header

	    if ($line =~ /^GENOME\tPATH/) {
		$genomeHeader = $line;
		$genomeHeader =~ s/^GENOME\tPATH//;
		$header = 0;
		next;
	    }

	    $header = 0;

	}

	#Make individual genome directory
	my @line_values = split(/\t/,$line,3);
	my ($genome,$path,$genomeAttributes);

	#Check that the correct values were in the input file and that formatting is correct
	if(scalar @line_values < 2){

	    print $lfh "ERROR: Formatting in $input_file is not correct. Parsing died at $line. Genome and path must be seperated by a tab.\n";
	    print STDERR "ERROR: Formatting in $input_file is not correct.";
	    print STDERR "Parsing died at $line. Genome and path must be seperated by a tab.\n";
	    exit(1);

	}else{

	    $genome = $line_values[0];
	    $path = $line_values[1];
	    $genomeAttributes = $line_values[2] if $line_values[2];

	}

	$genome =~ s/\s+$//; #removing trailing spaces
	$path =~ s/\s+$//; #removing trailing spaces

	if(exists $genomes_seen{$genome}){

	    print STDERR "WARNING: a genome is allowed to occur only once in the genome list file - ignoring the second line.\n";
	    print STDERR "$line\n$genomes_seen{$genome}\n";

	    next;

	}else{

	    $genomes_seen{$genome} = $line;

	}

	print $lfh "|INFO:Processing genome $genome\n";
	my $logfile = "$OUTPUT/$genome/makeblastdb.log";

	unless((-s $logfile && $opts{skip_blast}) || $new_alleles){
	    make_directory($genome) unless(-d $genome);

	    #Copy over genome sequences
	    #Only copy over if given input list, ncbi downloads have no need to copy over
	    copy_genome_sequences($path,$genome);

	}

	#Run blastall
	my $blast_file = run_blastall($genome,$new_alleles); #for $new_alleles this is just used to return the file name

	#Die if blast results were zero
	unless(-s $blast_file){

	    print $lfh "|WARNING: Blast results for $genome where empty. Please check genome and remove from analysis if necessary\n";
	    next;

	}

	#Run MLST Scripts
	my ($top_hits_file,$th_log)= run_top_hits($blast_file,$genome,$new_alleles,$tfh); #for $new_alleles this is just used to return the file name
	my $top_seqs_file = run_pullseqs($top_hits_file,$genome,$new_alleles); #for $new_alleles this is just used to return the file name

	push(@top_seqs_files, $top_seqs_file);
	push(@multi_copy, $th_log) if(-s $th_log);
	
	my ($st_out,$fa_out);

	if($opts{novel_schema}){

	    push(@fa_files,$top_seqs_file);

	}else{

	    ($st_out,$fa_out) = run_st_finder($top_seqs_file, $genome, $new_alleles, \%stMap, \%alleleMap, \@allelesOrdered, \%stAttributes, $genomeAttributes, $genomeHeader) ;
	    push(@fa_files,$fa_out);
	    push(@st_files,$st_out);

	}

    }


    #Add previous genome's list output files
    unless($opts{retype}){
	if($opts{original_input_file}){

	    foreach (@ORIG_GENOME_LIST){

		my @values = split(/\t/,$_);

		if($opts{novel_schema}){

		    push(@fa_files, "$ORIG_DIR/" . $values[0] . "/" . $values[0] . "_hits_top_seqs.fa");

		}else{

		    push(@fa_files,"$ORIG_DIR/" . $values[0] . "/" . $values[0] . "_new.fa");
		    push(@st_files,"$ORIG_DIR/" . $values[0] . "/" . $values[0] . "_ST.out");

		}

	    }

	}
    }

    return(\@st_files,\@fa_files,\@top_seqs_files,\@multi_copy);
}

sub clean_fasta_file{
    print $lfh "|Step: Cleaning fasta files and make them non redundant\n";
    my ($old_alleles,$new_alleles,$dir) = @_;

    my $allele_numbers = {};
    my $new_file = "$dir/appended_alleles_nr.fa";

    my $ofh = path($old_alleles)->filehandle("<");
    my $nfh = path($new_alleles)->filehandle("<");
    my $nrfh = path($new_file)->filehandle(">");

    while(<$ofh>){

	my $line = $_;

	if($line =~ /^>/){

	    my $allele = $line;
	    $allele =~ s/^>//;
	    $allele =~ s/\s.*//;
	    $allele =~ s/\s+//;

	    my ($gene,$identifier) = split(/\_/,$allele, 2);
	    
	    if($identifier =~ /^NOVEL\d+$/){

		$identifier =~ s/^NOVEL//;

		if(exists $allele_numbers->{$gene}){

		    if($identifier > $allele_numbers->{$gene}){

			$allele_numbers->{$gene} = $identifier;

		    }

		}else{

		    $allele_numbers->{$gene} = $identifier;

		}

	    }else{

		$allele_numbers->{$gene} = $identifier;
	    }

	}

    }

    my $current_sequence;
    my $unique_sequences; #hsh ref of unique sequences of new alleles
    my ($p_allele,$c_allele);

    while(<$nfh>){

	my $line = $_;
	$line =~ s/s+$//;

	if($line =~ /^>/){

	    my ($h,$temp_allele) = split(/>/,$line); #split def line
	    my @alleles = split(/\_/,$temp_allele, 2);
	    $c_allele = $alleles[0];
	    $c_allele =~ s/\s+$//;
	    $p_allele = $c_allele unless($p_allele);

	    if($current_sequence){

		#Removes newlines and other whiespace to make for one string
		$current_sequence =~ s/\s+//g;

		unless (exists $unique_sequences->{$current_sequence}->{$p_allele}){

		    if(exists $allele_numbers->{$p_allele}){

			$allele_numbers->{$p_allele}++;

		    }else{

			$allele_numbers->{$p_allele} = 1;

		    }

		    $unique_sequences->{$current_sequence}->{$p_allele} = "NOVEL" . $allele_numbers->{$p_allele}; #store unique seq

		}

		$current_sequence = ""; #reset seq variable
		$p_allele = $c_allele;

	    }

	}else{

	    $current_sequence .= $line;

	}

    }

    #Account for final sequence
    if($current_sequence){

	$current_sequence =~ s/\s+//g;

	unless(exists $unique_sequences->{$current_sequence}->{$p_allele}){

	    if(exists $allele_numbers->{$p_allele}){

		$allele_numbers->{$p_allele}++;

	    }else{

		$allele_numbers->{$p_allele} = 1;

	    }

	    $unique_sequences->{$current_sequence}->{$p_allele} = "NOVEL" . $allele_numbers->{$p_allele}; #store unique seq

	}

    }

    #Print Unique Sequences to File
    #Cat with original allele
    #Clean Fasta file
    foreach my $sequence (keys %$unique_sequences){

	foreach my $allele (keys %{$unique_sequences->{$sequence}}){

	    print $nrfh ">$allele" . "_" . $unique_sequences->{$sequence}->{$allele} . "\n";
	    print $nrfh "$sequence\n";

	}
    }

    close $nrfh;

    return($new_file);
}

sub init_st_finder{
    my($scheme_file, $alleles_file, $stMap, $alleleMap, $allelesOrdered, $stAttributes, $append) = @_;

    # Parse through ST-schema file.
    my $stSchemaFile_fh = path($scheme_file)->filehandle("<");
    my $ST_type_size;
    my $ST_attr_size;
    my $alleles_seen_ST = {};
    my $alleles_seen_alleles = {};
    my $header =1;

    while (my $line = <$stSchemaFile_fh>) {
	chomp $line;
        my @data = split("\t",$line);
	my $ST = shift(@data);       # Remove the first element (Sequence Type) from the line and store it.

	# If this line is the header, store the line's array, which is the ordered list of allele names, in allelesOrdered.
	if ($header) {
	    $header = 0;
	    if ($ST eq "ST") {
		@{$allelesOrdered} = @data;
		$ST_type_size = scalar @data;
		$ST_attr_size = 0;
	    } elsif ($ST =~ /^ST\((\d+)\)$/) {
		$ST_type_size = $1;
		if ($ST_type_size > (scalar @data)) {
		    die "ERROR: fewer header columns than specified in $ST:\n$line\n";
		}
		@{$allelesOrdered} = @data[0 .. ($ST_type_size - 1)];
		$ST_attr_size = (scalar @data) - $ST_type_size;
		if ($ST_attr_size > 0) {
		    $stAttributes->{'HEADER'} = join("\t", @data[$ST_type_size .. $#data]);
		    my $unknown = "";
		    foreach ($ST_type_size .. $#data) {
			$unknown .= "\t";
		    }
		    $stAttributes->{'UNKNOWN'} = $unknown;
		}
	    } else {
		die "ERROR: header line for ST schema must begin with ST tab or ST(#) tab - not:\n$line\n";
	    }
	}
	# Otherwise, store the ST number in stMap as a value and the line (the ST's corresponding alleles) as a key.
	else {
	    my $size = scalar @data;
	    if (($size != ($ST_type_size + $ST_attr_size)) && ($size != $ST_type_size)) {
		die "ERROR: ST$ST does not have the same number of columns as the header line or the alleles only: $size versus $ST_type_size alleles and $ST_attr_size attributes.\n";
	    }
	    my $st_type = join("\t",@data[0 .. ($ST_type_size - 1)]);
	    unless(exists $stMap->{$st_type}){
		$stMap->{$st_type} = $ST;
		if ($ST_attr_size > 0) {
		    if ($size != $ST_type_size) {
			$stAttributes->{$ST} = join("\t", @data[$ST_type_size .. $#data]);
		    } else {
			$stAttributes->{$ST} = $stAttributes->{'UNKNOWN'};
		    }
		}
	    } else {
		print STDERR "WARNING: duplicated ST entry: ($st_type) in $scheme_file with ST values $stMap->{$st_type} and $ST: $stMap->{$st_type} will be used and $ST ignored,\n";
	    }
	    foreach my $index (0 .. ($ST_type_size - 1)) {
		$alleles_seen_ST->{$allelesOrdered->[$index]}{$data[$index]} = 1;
	    }
	}
    }

    close($stSchemaFile_fh);

    # Open Bio::SeqIO object for mlstAllelesFile.
    my $inMlstAlleles = Bio::SeqIO->new(-file => "<$alleles_file", -format => "fasta",);

    # Parse through the MLST alleles and create the allele mapping.
    while (my $mlstAllele = $inMlstAlleles->next_seq) {

	my($allele,$id) = split(/\_/,$mlstAllele->primary_id,2);

	$alleleMap->{lc($mlstAllele->seq)}->{$allele} = $id;

	if ($mlstAllele->primary_id !~ /\_/) {
	    die "FATAL: allele $mlstAllele->primary_id is not in expected format of GeneName_Identifier where GeneName the name of a gene with no special characters such as _ and Identifier is typically a number.\n";
	}

	$alleles_seen_alleles->{$allele}{$id} = 1;
    }


    #perform sanity checking
    foreach my $allele (keys %{$alleles_seen_ST}) {
	foreach my $id (keys $alleles_seen_ST->{$allele}) {
	    if (!defined $alleles_seen_alleles->{$allele}{$id}) {
		print STDERR "WARNING: Allele $allele\_$id is in $scheme_file but not in $alleles_file.\n";
	    }
	}
    }

    if (!$append) {
	foreach my $allele (keys %{$alleles_seen_alleles}) {
	    foreach my $id (keys $alleles_seen_alleles->{$allele}) {
		if (!defined $alleles_seen_ST->{$allele}{$id}) {
		    print STDERR "WARNING: Allele $allele\_$id is in $alleles_file but not in $scheme_file.\n";
		}
	    }
	}
    }

    return;
}
sub run_st_finder{
    my($top_seqs_file, $identifier, $new_alleles, $stMap, $alleleMap, $allelesOrdered, $stAttributes, $genomeAttributes, $genomeHeader) = @_;
    my %allelesFound;                  # hash of allele name to matching MLST allele.
    my %query_sequences;               # hash ref of query sequences

    my $outputfile;

    #Open File for NEW Allele sequences
    my($new_alleles_file,$new_alleles_fh);

    if($new_alleles){
	$outputfile = "$OUTPUT/$identifier/$identifier" . "_NEW_alleles_ST.out";
    }else{
	$outputfile = "$OUTPUT/$identifier/$identifier" . "_ST.out";
	$new_alleles_file = "$OUTPUT/$identifier/$identifier" . "_new.fa";

	$new_alleles_fh = path($new_alleles_file)->filehandle(">");
    }

    # Open filehandle for output file.
    my $output_fh = path($outputfile)->filehandle(">");

    print $output_fh "Sample\tST\t".join("\t", @{$allelesOrdered});

    if (keys %{ $stAttributes }) {
	print $output_fh "\t$stAttributes->{'HEADER'}";
    }
    if ($genomeHeader) {
	print $output_fh "\t$genomeHeader";
    }
    print $output_fh "\n";

    # Open Bio::SeqIO object for top_seqs_file.
    my $inQueries = Bio::SeqIO->new(-file => "<$top_seqs_file", -format => "fasta",);
    my $query_count = 1;

    # Parse through the query.
    while (my $query = $inQueries->next_seq) {
	my $queryName = $query->primary_id;                      # Query name.
	my $querySeq = $query->seq;                              # Query sequence.

	#Split name to get queryAllele
	my($queryAllele,$scheme) = split(/\_/,$queryName,2);

	my $unique_id = $queryAllele . ":" . $query_count;
	
	$allelesFound{$queryAllele}{$unique_id} = "NEW";
	$query_sequences{$queryAllele}{$unique_id} = $querySeq;

       	# If the query's sequence begins with "SHORT", declare the queryAllele's hit as SHORT.
	if ($querySeq =~ /^SHORT/) {
	    $allelesFound{$queryAllele}{$unique_id} = "SHORT";
	} elsif ($querySeq =~ /^PSEUDO/) {
	    $allelesFound{$queryAllele}{$unique_id} = "PSEUDO";
	} elsif ($querySeq =~ /^5'PRTL/) {
	    $allelesFound{$queryAllele}{$unique_id} = "5'PRTL";
	} elsif ($querySeq =~ /^3'PRTL/) {
	    $allelesFound{$queryAllele}{$unique_id} = "3'PRTL";
	} elsif (defined $alleleMap->{lc($querySeq)}->{$queryAllele}) {

	    #Note: Different alleles COULD have the same sequence
	    foreach my $mlstAlleleName(keys $alleleMap->{lc($querySeq)}){

		#If the query allele is a member of the MLST allele
		#(i.e. the query allele name forms the beginning of the MLST allele name),
		#it is a proper match.
		if ($mlstAlleleName eq $queryAllele) {
		    $allelesFound{$queryAllele}{$unique_id} = $alleleMap->{lc($querySeq)}->{$mlstAlleleName};
		}
	    }
	}

	$query_count++;
    }
    
    # Modify allelesFound to add _MC where the is a multi copy
    foreach my $allele(keys %allelesFound){

	my $count_alleles = scalar(keys %{$allelesFound{$allele}});
	my $labels;

	if($count_alleles > 1){ #means MC

	    foreach my $unique(keys %{$allelesFound{$allele}}){

		my $label = $allelesFound{$allele}{$unique};

		if(exists $labels->{$label}){

		    $allelesFound{$allele}{$unique} = $label . "_MC";

		    my $orig_label = $labels->{$label};
		    $allelesFound{$allele}{$orig_label} = $label . "_MC";

		}else{

		    $labels->{$label} = $unique;
		}

	    }

	}
    }
    
    #Determine which alleles have multiple variants
    my @multi_variants;
    my @base_string;
    my %locations;

    #Stores alleles with multi variants
    #Creates base string, with the alleles that do not have variants
    #Variants will be 'plugged' in later
    my $location_count = 0;

    foreach my $key (@$allelesOrdered){

	$locations{$key} = $location_count;
	$location_count++;

	my @v;

	if(exists $allelesFound{$key}){
	    @v = keys $allelesFound{$key};
	}else{
	    $allelesFound{$key}{$key} = "MISSING";
	}

	my $count = scalar (@v);

	if($count > 1){
	    
	    push(@multi_variants,\@v);
	    push(@base_string,"");

	}else{

	    if(@v){
		my $schema_num = $allelesFound{$key}{$v[0]};
		push(@base_string,$schema_num);
	    }else{
		push(@base_string,"MISSING");
	    }
	}
    }

    my @variants;
    push @variants, clone(\@base_string);

    foreach my $gene(@multi_variants){

	my($base_allele,$unique_num) = split(/:/,$gene->[0]);
      
	my $values_to_add = scalar(@variants);

	my $size = scalar(@$gene);
	grow_variants( \@variants, $size );

	my $allele_index = $locations{$base_allele};

	my $values_added = 0;
	my $index = 0;  #increment this until it's time to move on.

	for my $allele ( @$gene ) {
     
	    my $allele_schema = $allelesFound{$base_allele}{$allele};

	    while ( $values_added < $values_to_add ) {
		# Fill the location for allele at this entry
		$variants[$index][$allele_index] = $allele_schema;

		# then increment both index and values_added
		$index++;
		$values_added++;
	    }

	    $values_added = 0; # reset values added between alleles.

	}

    }

    # stKey - tab-delimited string of allelesFound in order of allelesOrdered.
    my $stKey = "";

    #Goal: To print all combinations of ST types, with all
    #possible allele variations shown

    foreach my $allele (@{$allelesOrdered}) {

	#Print new allele sequences to fasta file
	if(exists $allelesFound{$allele}){

	    foreach my $a(keys $allelesFound{$allele}){
	
		if($allelesFound{$allele}{$a} eq 'NEW'){

		    print $new_alleles_fh ">$allele\n" unless $new_alleles;
		    print $new_alleles_fh "$query_sequences{$allele}{$a}\n" unless $new_alleles;

		}
		
	    }
	    
	}
	
    }

    # Clean up @variants to remove duplicate type string
    my $unique_variants;
    foreach (@variants){
	my $variant_value = $_;
	my $variant_string = join("\t",@$variant_value);

	$unique_variants->{$variant_string} = 1;

    }

   # If the stKey exists in stMap, print both ST found and stKey to output file.
    foreach my $stKey (keys %$unique_variants){

	my $ST;

	if (exists $stMap->{$stKey}) {

	    $ST = $stMap->{$stKey};
	}

	# Otherwise, print "UNKNOWN" and stKey to output file.
	# Note: this will happen if any of the alleles in stKey are SHORT or NEW.
	else {
	    $ST = "UNKNOWN";
	}


	print $output_fh "$identifier\t$ST\t$stKey";

	if (keys %{ $stAttributes }) {
	    print $output_fh "\t$stAttributes->{$ST}";
	}

	if ($genomeAttributes) {
	    print $output_fh "\t$genomeAttributes";
	}

	print $output_fh "\n";
    }

    return($outputfile,$new_alleles_file);
}
sub grow_variants {

    my ( $variants, $size ) = @_;

    # Make a local copy of 1 'unit', being the current @variants array.
    # This allows us to add arbitrary number of copies to itself.
    # Also, we need to make sure these aren't copies of the same reference, or we'll be
    # overwriting the same locations when making our fills, so we'll use Clone.pm
    # to make 'deep copies'. However, we also need to make sure to clone the entries in the
    # copy, or, again, we'll be accessing the same location when making different fills... yeesh.
    my $copy = clone($variants);

    while ( $size > 1 ) {
	# Needs to be 1, not 0, since we're starting with 1 copy of @variants in @variants.

	for my $entry ( @$copy ) {
	    push @$variants, clone($entry);
	}
	$size--; # Don't forget to decrement size :)

    }
}

sub run_pullseqs{
    my ($top_hits,$genome,$new_alleles) = @_;

    my $file = "$OUTPUT/$genome/$genome" . "_hits_top_seqs.fa";

    my $cmd = "perl $Bin/pullseqs.pl";
    $cmd .= " -fastafile $OUTPUT/$genome/$genome.fasta";
    $cmd .= " -blastfile $top_hits";
    $cmd .= " -seeds $opts{seed_file}";
    $cmd .= " -max_mismatch $opts{blast_length}" if $opts{blast_length};

    unless(-s $file && $opts{skip_blast} || $new_alleles){
	print $lfh "Running: $cmd\n";
	system($cmd) == 0 || die("ERROR: $cmd failed");
    }

    return($file);
}
sub run_top_hits{
    my ($blast,$genome,$new_alleles,$log_file) = @_;

    my $file = "$OUTPUT/$genome/$genome" . "_hits_top.txt";
    my $multi_file = "$OUTPUT/$genome/$genome" . "_multi_copy_hits.txt";
   
    my $cmd = "perl $Bin/tophits.pl";
    $cmd .= " $blast";
    $cmd .= " $multi_file";
    $cmd .= " multi" if $opts{multi_copy};
    
    #skips running the top hits if file already
    #exists and the option --skip_blast is set

    unless(-s $file && $opts{skip_blast} || $new_alleles){
	print $lfh "Running: $cmd\n";

	capture_merged{

	    system($cmd) == 0 || die("ERROR running tophits!", __LINE__);

	} stdout => $log_file;
    }

    return ($file,$multi_file);

}
sub run_blastall{
    my ($genome,$new_alleles) = @_;

    my $fasta_file = "$OUTPUT/$genome/$genome.fasta";
    my $output_blast = "$OUTPUT/$genome/$genome" . "_hits.txt";

    #Runs blast unless the blast file is already there AND
    #the user specified the skip_blast option
    #Need to customize blast tabular output to include the query length and subject length as the final columns - first twelve columns remain the same
    unless((-s $output_blast && $opts{skip_blast}) || $new_alleles){

	if(-s $fasta_file){

	    my $cmd = $BLAST_CMD . " -task blastn -db $fasta_file -query $opts{seed_file} -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -out $output_blast -qcov_hsp_perc 30 -max_target_seqs 5 -num_threads 4";
	    print $lfh "Running: $cmd\n";
	    system($cmd) == 0 || die("ERROR: $cmd failed");

	}else{

	    print $lfh "WARNING: $fasta_file is size zero or doesn't exist so no blast could
be run\n";

	}

    }

    return $output_blast;
}

sub make_directory{
    my $genome = shift;
    my $dir = "$OUTPUT/$genome";

    mkdir($dir) unless (-d $dir);
    die("ERROR: Error creating $dir\n") unless (-d $dir);
}

sub copy_genome_sequences{
    my ($location,$genome) = @_;

    $location = _trim($location);
    $genome = _trim($genome);

    my $genome_dir = "$OUTPUT/$genome/";
    my $fasta = "$genome.fasta";
    my $logfile = "$OUTPUT/$genome/makeblastdb.log";

    #Create sym link to fasta file if not downloaded
    if($opts{input_file} || $opts{input_path}){
	my $real_path = path($location)->realpath;
	symlink($real_path,"$OUTPUT/$genome/$fasta");
    }else{
	$location = "$OUTPUT/$genome/$fasta";
    }

    if(-s $location){

	#Format fasta file into blast database
	my $cmd = $FORMATDB_EXEC . " -input_type fasta -in $location -dbtype nucl -logfile $logfile -out $genome_dir/$fasta";
	system($cmd) == 0 || die("ERROR: $cmd failed");

    }

}

sub create_itol_file{
	my $typer_file = shift;
	my $cmd = "perl $Bin/create_itol.pl";
	$cmd .= " --input_file $typer_file";
	print $lfh "|Step: Creating ITOL Annotation file.";
	print $lfh "Running $cmd\n";
	system($cmd) == 0 || die("ERROR: $cmd failed");
}

sub _trim{
    my $value = shift;

    $value =~ s/\s+$//;
    $value =~ s/^\s+//;

    return $value;

}
sub cleanup_files{

    my $download_dir = "$OUTPUT/logs/download";

    #make log file and move all logs to there

    mkdir($download_dir) unless (-d $download_dir);

    #move download files
    my @accession_download = glob("$OUTPUT/accession.*");
    my @biosample_download = glob("$OUTPUT/biosample.*");
    my @org_download = glob("$OUTPUT/org_search.*");
    my @curl = glob("$OUTPUT/*curl*");
    my @logs = glob("$OUTPUT/*log");

    unlink glob "$OUTPUT/*assembly_stats.txt";
    unlink "$OUTPUT/gb.list";
    unlink "$OUTPUT/org.list";

    foreach(@curl){
	move($_,$download_dir);
    }
    foreach(@logs){
	move($_,$log_dir);
    }
    foreach(@accession_download){
	move($_,$download_dir);
    }
    foreach(@biosample_download){
	move($_,$download_dir);
    }
    foreach(@org_download){
	move($_,$download_dir);
    }
}
sub check_params{

    my $error = "";

    if($opts{config}){
	if ($opts{config} !~ /^\//) {
	    $opts{config} = "$START_CWD/$opts{config}";
	}
	$error .= "ERROR: $opts{config} does not exist or is size zero\n" unless (-s $opts{config});
    }elsif("$Bin/default_config"){
	$opts{config} = "$Bin/default_config.ini";
    }else{
        $error .= "ERROR: Could not find valid config file in current working directory. Please provide one using --config\n";
    }

    if($opts{input_file}){
	if ($opts{input_file} !~ /^\//) {
	    $opts{input_file} = "$START_CWD/$opts{input_file}";
	}
	$error .= "ERROR: $opts{input_file} does not exist or is size zero\n" unless (-s $opts{input_file});

	if($opts{hmm_model}){
	    $error .= "ERROR: Must provide --gb_list if running HMM models without NCBI downloading\n" unless($opts{gb_list});
	}

	if($opts{org_search} || $opts{biosample_list} || $opts{accession_list} || $opts{input_path}){
	    $error .= "ERROR: Please provide only one input option: --input_file, --org_search, --biosample_list, --accession_list, --input_path\n";
	}

    }elsif(!($opts{org_search} || $opts{biosample_list} || $opts{accession_list} || $opts{input_path})){
	$error .= "ERROR: Option --input_file/-i is required\n";
    }

    if($opts{seed_file}){

	if ($opts{seed_file} !~ /^\//) {
	    $opts{seed_file} = "$START_CWD/$opts{seed_file}";
	}
	$error .= "ERROR: $opts{seed_file} does not exist or is size zero\n" unless (-s $opts{seed_file});
    }elsif($opts{incrememnt}){
	$error .= "ERROR: Option --seed_file is required when using the --incrememnt option\n";
    }
		if ($opts{download_schema}) {
			if ($opts{scheme}) {
				$error .= "ERROR: Cannot use the download schema option when providing schema\n";
				die("ERROR: Cannot use the download schema option when providing schema.\n");
			}
			if ($opts{alleles}) {
				$error .= "ERROR: Cannot use the download schema option when providing alleles file\n";
				die("ERROR: Cannot use the download schema option when providing alleles file.\n");
			}
			if ($opts{seed_file}) {
				$error .= "ERROR: Cannot use the download schema option when providing a seed file\n";
				die("ERROR: Cannot use the download schema option when providing a seed file.\n");
			}
		}
    unless($opts{download_schema}){
	if($opts{scheme}){
	    if ($opts{scheme} !~ /^\//) {
		$opts{scheme} = "$START_CWD/$opts{scheme}";
	    }
	    $error .= "ERROR: $opts{scheme} does not exist or is size zero\n" unless (-s $opts{scheme});
	}else{
	    $error .= "ERROR: Option --scheme/-m is required, unless doing --novel_schema\n" unless($opts{novel_schema});
	}

	if($opts{alleles}){
	    if ($opts{alleles} !~ /^\//) {
		$opts{alleles} = "$START_CWD/$opts{alleles}";
	    }
	    $error .= "ERROR: $opts{alleles} does not exist or is size zero\n" unless (-s $opts{alleles});
	}else{
	    $error .= "ERROR: Option --alleles/-a is required, unless doing --novel_schema\n" unless($opts{novel_schema});
	}
    }

    if($opts{novel_schema} && $opts{append_schema}){
	$error .= "ERROR: Can only use either --novel_schema OR --append_schema. Can not use both options\n";
    }

    if($opts{tree}){

	$error .= "ERROR: the arugment value '$opts{tree}' for tree option is not valid. The valid entries are 'raxml' or 'fasttree'\n" unless(lc($opts{tree}) eq 'raxml' || lc($opts{tree}) eq 'fasttree');

    }

    if($opts{genome_type}){

	$error .= "ERROR: the argument value '$opts{genome_type}' for genome_type is not valid. The valid entries are 'cg' or 'wgs'\n" unless(lc($opts{genome_type}) eq 'cg' || lc($opts{genome_type}) eq 'wgs');
    }

    my $output;

    if($opts{output}){

	if(-d $opts{output}){
	    $output = $opts{output};
	}else{
	    #Make directory path if it does not exist
	    $output = $opts{output};
	    mkdir($opts{output});
	}
	if ($output !~ /^\//) {
	    $output = "$START_CWD/$output";
	}
    }else{

	$output = $START_CWD;
    }


    if($opts{previous_output}){
	$output = $START_CWD;
    }

    if($opts{previous_output}){
	if($opts{original_input_file}){
	    $error .= "$opts{original_input_file} does not exist or is size zero\n" unless(-s $opts{original_input_file});
	}else{
	    $error .= "Must provide --original_input_file if using --previous_output\n";
	}
    }

    my $orig_dir = $opts{previous_output} // $START_CWD;

    if($opts{hmm_model}){
	$error .= "ERROR:$opts{hmm_model} does not exist or is size zero\n" unless(-s $opts{hmm_model});
    }

    if($error){
	die($error);
    }else{
	return($output,$orig_dir);
    }

}
sub _cat {
    # Given a list of file names, concatonate the first through n minus one-th
    # onto the nth.

    my ( $output_fh, $input ) = ( @_ );
    my $print_header = 0;

    for ( @$input ) {

	my $file = $_;

	if(-s $file){

	    my $ifh = path($file)->filehandle("<");

	    while ( <$ifh> ) {

		my $line = $_;

		if($line =~ /^Sample/){

		    print $output_fh $line unless($print_header);
		    $print_header = 1;

		}else{

		    print $output_fh $line;

		}
	    }

	}

    }

}

sub _cat_fa_files {
    # Given a list of file names, concatonate the first through n minus one-th
    # onto the nth.

    my ( $output_fh, $input ) = ( @_ );
    my $print_header = 0;

    for ( @$input ) {
	if(-s $_){

	    my $ifh = path($_)->filehandle("<");

	    while ( <$ifh> ) {
		if($_ =~ /Sample/){
		    print $output_fh $_ unless($print_header);
		    $print_header = 1;
		}else{
		    print $output_fh $_;
		}
	    }
	}
    }
}
