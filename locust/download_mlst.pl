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
use warnings;

=head1 NAME

download_mlst.pl -  A downloader for MLST profiles and allles from pubmlst.com

=head1 SYNOPSIS

  USAGE: download_mlst.pl --organism <Species Search Term>
                         --loci_num <Integer>
                         --help

=head1 OPTIONS

B<--organism, o>   : Term to search for MLST profile.

B<--loci_num, l>    : Number of expected loci in MLST schema.

B<--help, h>         : Display this help message.

=head1  DESCRIPTION

This program downloads the allele files and schema for use in LOCUST.

=head1  CONTACT

    Erin Beck
    ebeck@jcvi.org

=cut

use LWP::UserAgent;
use XML::Simple;
use Getopt::Long;
use Time::localtime;
use Bio::SeqIO;
use Pod::Usage;

sub download_files{
  my (%new_url_mapping, $schema_filename, @allele_files, $removal_string);
  my $ua = user_agent_command();
  my ($found_hash) = @_;
  my $directory = "logs/";
  unless (-d $directory){mkdir $directory or die $!};
  check_file_existence("logs/download.log");
  open (my $log_file, '>', "logs/download.log");
  my $tm = localtime;
  my ($day,$month,$year)=($tm->mday,($tm->mon + 1),($tm->year + 1900));
  print $log_file "$month-$day-$year\n";
  #found_hash should have one key-- found species
  while( my ($species, $xml_hash) = each %$found_hash){
    my $urls = $xml_hash ->{"mlst"}->{"database"}->{"loci"}->{"locus"};
    my $schema_file = $xml_hash ->{"mlst"}->{"database"}->{"profiles"}->{"url"};
    chomp $species;
    print $log_file "MLST Profile: $species\n";
    print "Schema: $schema_file\n";
    print $log_file "Schema: $schema_file\n";
    my %url_mapping;
    for (my $i = 0; $i < @{$urls}; $i ++){
      my $url_hash = @{$urls}[$i];
      my $url = $url_hash->{"url"};
      my $file = $url_hash->{"content"};
      chomp $file;
      $url_mapping{$file} = $url;
    }
    $removal_string = check_allele_name(\%url_mapping, $log_file);
    if ($removal_string eq ''){
      %new_url_mapping = %url_mapping;
    } else {
      %new_url_mapping = change_allele_name(\%url_mapping, $removal_string);
    }
    my $schema_content = download_content($schema_file);
    $schema_filename = $species . "_schema.txt";
    check_file_existence($schema_filename);
    open(FH, ">", $schema_filename) or die $!;
    print FH "$schema_content\n";
    close(FH);
  }
    while (my ($key, $value) = each %new_url_mapping){
      my $allele_content = download_content($value);
      my $allele_filename = $key . ".fa";
      check_file_existence($allele_filename);
      push (@allele_files, $allele_filename);
      open(FH, ">", $allele_filename) or die $!;
      print FH "$allele_content\n";
      close(FH);
  }
  return ($removal_string, \@allele_files, $schema_filename);

}

sub download_content(){
  my $ua = user_agent_command();
  my ($url) = @_;
  my $new_req = HTTP::Request->new(
    GET => $url,
  );
  my $res = $ua->request($new_req);
  my $content = $res->content;
  chomp $content;
  return $content;
}

sub check_allele_name(){
  my ($url_mapping, $log_file) = @_;
  my $underscore_count = 0;
  my $final_split_string = '';
  while (my ($key, $value) = each %$url_mapping){
    my $indv_count = ($key =~ tr/_//);
    $underscore_count = $underscore_count + $indv_count;
  }
  #Allels have prefix or suffix that needs to be fixed.
  if ($underscore_count == keys %$url_mapping) {
    print $log_file "Renaming files into the proper format for Locust.\n";
    #Generate sets to see whether it is prefix (left of split == 1 key) or suffix (right of split == 1 key)
    my (%left_of_underscore, %right_of_underscore, @str_to_remove, $log_info);
    while (my ($filename, $url_content) = each %$url_mapping){
      my ($left_split, $right_split) = split(_, $filename);
      $left_of_underscore{$left_split} = 1; #Left of split is prefix-- want right of _
      $right_of_underscore{$right_split} = 0; #Right of split is suffix-- want left of _
    }
    if (scalar keys %left_of_underscore == 1){
      @str_to_remove = keys %left_of_underscore;
      $final_split_string = $str_to_remove[0] . "_";
      $log_info = "prefix";
    }
    if (scalar keys %right_of_underscore == 1){
      @str_to_remove = keys %right_of_underscore;
      $final_split_string = "_" . $str_to_remove[0];
      $log_info = "suffix";
    }
    #Return $str_to_remove and remove first occurrence from FASTA file names, schema header, and sequence entries within FASTA file
    print $log_file "Reformatting files by removing the $final_split_string $log_info.\n";
    return $final_split_string;
  }  elsif ($underscore_count > 0){
    #At least one allele needs to be modified so it works in locust. Not suffix/prefix.
      print $log_file "Reformatting files by removing the $final_split_string _.\n";
      $final_split_string = "_";
      return $final_split_string;
  } else {
      return $final_split_string;
  }
}

sub change_allele_name{
  my ($url_mapping, $str_for_removal) = @_;
  my %new_mapping_file;
  while (my ($allele_name, $url) = each %$url_mapping){
    $allele_name =~ s/$str_for_removal//;
    $new_mapping_file{$allele_name} = $url;
  }
  return %new_mapping_file;
}

sub change_fasta_allele_names{
  (my $allele_file, my $str_for_removal) = @_;
  my $in = Bio::SeqIO->new( -file => $allele_file, -format => 'Fasta' );
  check_file_existence($allele_file);
  open(FH, ">", $allele_file);
  while ( my $seq = $in->next_seq() ) {
    my $oldSeq = $seq->display_id;
    $oldSeq =~ s/$str_for_removal//;
    my $sequence = $seq->seq;
    print FH ">$oldSeq\n$sequence\n";
  }
  close(FH);
}

sub change_schema_header{
  (my $schema_file, my $str_for_removal, my $loci_count) = @_;
  my $tmp = 'schema_temp.txt';
  check_file_existence($tmp);
  open (OUT, ">", $tmp);
  open (IN, '<', $schema_file) or die "Couldn't open schema file. $!";

  while (my $line = <IN>){
    chomp $line;
    if ($. == 1){
      my @header_columns = split(/\t/, $line);
      for (my $i = 0; $i < scalar @header_columns; $i++){
        if (not $loci_count == 0){
        if ($i == 0){
          $header_columns[$i] = $header_columns[$i] . "($loci_count)";
        }
      }
        $header_columns[$i] =~ s/$str_for_removal//;
      }
      my $final_header = join("\t", @header_columns);
      print OUT "$final_header\n";
    } else {
    print OUT "$line\n";
    }
  }
  rename($tmp, $schema_file);
}


sub user_agent_command{
  my $ua = LWP::UserAgent->new(
      ssl_opts => { verify_hostname => 0 },
      protocols_allowed => ['https'],
  );
  return $ua;
}

sub check_file_existence{
  my ($new_filename) = @_;
  if (-e $new_filename){unlink($new_filename) or die $!};
}

sub get_xml_file{
  my $ua = user_agent_command();
  my $req = HTTP::Request->new(
      GET => 'https://pubmlst.org/data/dbases.xml',
  );

  my $res = $ua->request($req);
  my $xmlContent = $res->content;
  return $xmlContent;
}

sub parse_xml_file{
  my ($search_term) = shift;
  my ($xml_file) = shift;
  my $parser = new XML::Simple;
  my $dom = $parser->XMLin($xml_file);
  my @entries = @{ $dom->{'species'} };
  my @found_entries;
  my %found_hash;
  foreach my $entry (@entries){
    my $species = $entry->{'content'};
    if ( $species =~ m/\Q$search_term/i){
      push (@found_entries, $species);
      my $key = join("_", split(/ /, $species));
      $found_hash{$key} = $entry;
    }
  }
    if (scalar @found_entries > 1){
      foreach my $entry (@found_entries){
        print "\n$entry";
      }
      die "\nMore than 1 MLST Profile matched your species of interest. Please choose one of the above.\n"
    }
    elsif (scalar @found_entries eq 0) {
      die "\nYou didn't match any MLST profile name. Check your spelling. Otherwise, check here for available schema's: https://pubmlst.org/data/\n"
    } else {
      my $found_entry = $found_entries[0];
      chomp $found_entry;
      print "Using MLST Profile: $found_entry\n";
    }
    return %found_hash;
}

sub make_concatenate_allele_file{
  my ($allele_file, $species, $final_split_string) = @_;
  my $concat_allele_file = $species . "_alleles.fa";
  check_file_existence($concat_allele_file);
  open(CAT_ALLELE_FILE, ">>", $concat_allele_file);
  foreach my $file (@$allele_file){
    change_fasta_allele_names($file, $final_split_string);
    if (open my $input_fasta, "<", $file){
      while (my $line = <$input_fasta>){
        chomp $line;
        print CAT_ALLELE_FILE "$line\n";
      }
      close $file;
    } else {
      die "Couldn't open individual allele file $file for concatenation.\n";
    }
  }
  close CAT_ALLELE_FILE;
  }




my %opts;
GetOptions( \%opts,
  'help|h',
  'organism|o=s',
  'loci_num|l=i',
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

#Main Function
my $search_term = $opts{organism} || die "\nEnter an Organism Name for MLST Download with '-o'.\n\n";
my $loci_num = $opts{loci_num} // 0;
my ($xmlContent) = get_xml_file();
my (%found_hash) = parse_xml_file($search_term, $xmlContent);
my $species = (keys %found_hash)[0]; #Will only be one key in hash.
chomp $species;
my ($final_split_string, $allele_files, $schema_filename) = download_files(\%found_hash);
change_schema_header($schema_filename, $final_split_string, $loci_num);
make_concatenate_allele_file($allele_files, $species, $final_split_string);
