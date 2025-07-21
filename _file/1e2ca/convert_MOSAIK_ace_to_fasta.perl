#!/usr/bin/perl

###############
# written by shyamala iyer
# University of Washington, Seattle

## MODIFIED: November 2012

# Script to convert a Mosaik ace file into fasta alignment file
# All the the sequences are converted into fasta format
# A single file is created and sequences that fall below a min read length are removed

## Individual read coverage for al positions is also listed
## Read coverage is the number of reads that cover a certain position ( Inserted bases will not be listed for read coverage)
###############

use strict;

my $usage = "perl convert_MOSAIK_ace_to_fasta.pl inAceFile output_fasta_file output_stat_file\n";

my $inFile = shift or die $usage;
my $output_fasta_file = shift or die $usage;
my $output_stat_file = shift or die $usage;

my ($inFileName) = split /\.ace/, $inFile;

our @in = ();


our @read_names = ();

our %all_read_start = ();  #read name, key == id, value is padded start position of read

our %all_read_seq = (); #just sequence of all reads key== id, value is sequence as a string

our $num_reads = 0;

our $ref_seq = "";
our $ref_name = "";
our $ref_len = 0;

our $padded_ref_len = 0;
our $padded_ref_start = 0;
our $padded_ref_end = 0;

our %all_read_start_pos = (); # store the start positions of all reads
our %all_read_length = ();

our %ref_gap_pos = (); ## gap positions in the reference sequence
our %ref_pos = ();

our %ref_G_pos = ();

our %possible_hyper_mutated_seq = ();
our %discarded_seq = ();

our %read_coverage = ();
our %duplicate_reads = ();

our $MIN_READ_LEN = 100;

open (INPUT_FILE, $inFile) or die ("couldn't open $inFile: $!\n");
@in = <INPUT_FILE>;
close INPUT_FILE;

foreach my $line (@in) {
  chomp $line;
  next if $line =~ /^\s*$/; # move to next line if empty
  
  if($line =~ m/^AS\s+(\S+)\s+(\S+)/){
    $num_reads = $2;
    $num_reads = $num_reads -1; # remove ref
  }
  elsif($line =~ /^CO\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/){
    $ref_name = $1;
    $padded_ref_len = $2;
    print "ref name = $ref_name\n";
    print "Padded length = $padded_ref_len\n";
  }
  
  # info about all the reads
  elsif ($line =~ /^AF\s+(\S+)\s+(\S+)\s+(\S+)$/) {
    
    my $read_str = $1;
    my $read_UC = $2; #read orientation 
    my $read_start_pos = $3;

    my @cols = split(/\./, $read_str);
    my $read_name = shift @cols;
    if((length $read_name > 0) && (!($read_name =~ /MosaikReference/)) ){
	
      $all_read_start{$read_name} = $read_start_pos;
           
      if(exists $all_read_start_pos{$read_start_pos}){
	# dont do anything
      }
      else{
	$all_read_start_pos{$read_start_pos} = 1;
      }
    } # read not mosaik ref
  }# end of elsif AF
}

my @sorted_pos = sort {$a <=> $b} keys %all_read_start_pos;
#print "@sorted_pos\n";

# Now get to the sequence part
for(my $i = 0; $i < scalar @in; $i++){
  
  if($in[$i] =~ m/^RD\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/){
    my $read_str = $1;
    my $seq_len = $2;
    
    if($read_str =~ /MosaikReference/){
      # reference   
      for (my $j = $i+1; $j < scalar @in; $j++){
	
	if($in[$j] =~ m/^QA*/){
	  # all ref sequence lines have been read
	  $j = scalar @in;
	}
	elsif($in[$j] =~ m/^\s+/){
	  # do nothing, empty lines
	}
	else{
	  chomp $in[$j];
	  #in this string, substitute '*' with '-'
	  my $curr_line = $in[$j];
	  $curr_line =~ s/\*/\-/g;
	  $curr_line =~ s/\s//g; $curr_line =~ s/\t//g;
	  $ref_seq = "$ref_seq"."$curr_line";
	  my @ref = split(//, $ref_seq);
	  $ref_len = length $ref_seq;
	  

	  %ref_pos = get_nt_pos(\@ref); # 1 based index for positions
	}# end of else
      }
    }
    else{
      # read line that is not reference
      # for each of the reads, the following needs to be done
      # 1) substitute "*s" with "-"
      # 2) pad the starting and trailing regions of the read with -s
      # 3) store the read starting and ending position (wrt to ref with gaps)

      my @cols = split(/\./, $read_str);
      my $curr_read_id = shift @cols;
      my $curr_read_seq = "";
      # for this read ID, get the starting position of the read
      my $read_pad_start = $all_read_start{$curr_read_id};

      for (my $j = $i+1; $j < scalar @in; $j++){
	if($in[$j] =~ m/^QA*/){
	  # all read sequence lines have been read
	  $j = scalar @in;
	}
	elsif($in[$j] =~ m/^\s+/){
	  # do nothing
	}
	else{
	  chomp $in[$j];
	  my $curr_line = $in[$j];
	  $curr_line =~ s/\*/\-/g;
	  $curr_line =~ s/\s//g; $curr_line =~ s/\t//g;
	  $curr_read_seq = "$curr_read_seq"."$curr_line";
	}
      }# end of for loop

      my $len = length $curr_read_seq;

      $all_read_length{$curr_read_id} = $len;
      my $tmp = $curr_read_seq;
      $tmp =~ s/\-//g;
      my $len_wo_gaps = length $tmp;
      

      if($len_wo_gaps >= $MIN_READ_LEN){
	
	# At this point, need to pad the read sequence with correct number of leading and trailing dashes
	my $lead_dashes = "";
	my $trailing_dashes = "";
	
	for(my $k = 1; $k < $read_pad_start; $k++){
	  $lead_dashes = "$lead_dashes"."-";
	}
	# get the number of trailing dashes that need to be filled
	my $num_dashes = ($padded_ref_len - ( ($read_pad_start-1) + ($seq_len)) );
	for (my $l =0; $l < $num_dashes; $l++){
	  $trailing_dashes = "$trailing_dashes"."-";
	}
	$curr_read_seq = "$lead_dashes"."$curr_read_seq"."$trailing_dashes";
	if(exists $all_read_seq{$curr_read_id}){
	  print "This read $curr_read_id is present more than once in the .ace file. This duplicate will be placed in a seperate file.\n";
	  $duplicate_reads{$curr_read_id} = $curr_read_seq;
	}
	else{
	  $all_read_seq{$curr_read_id} = $curr_read_seq;
	  push(@read_names, $curr_read_id);
	}
      }
      else{
	$discarded_seq{$curr_read_id} = $tmp;
      }
    }
  }
}


foreach my $p(sort {$a <=> $b} keys %ref_pos){
  my $coverage = get_read_coverage($p);
  $read_coverage{$p} = $coverage;
}

print "Reference sequence: \n$ref_seq\n";
print "reference alignment length: $ref_len\n";
my $num = scalar keys %all_read_seq;
print "Number of read clusters in the ace file $num\n";
print "writing reads to output files ...\n";

########## print the reads correctly #####

open (OUTPUT_FILE, ">$output_stat_file");
print OUTPUT_FILE "reference alignment length: $ref_len\n";
print OUTPUT_FILE "Number of read clusters in the ace file $num\n";
print OUTPUT_FILE "read coverage for positions in alignment\n";
foreach my $p(sort {$a <=> $b} keys %read_coverage){
  print OUTPUT_FILE "$p\t$read_coverage{$p}\n";
}
close OUTPUT_FILE;

open(OUTPUT_FILE,">$output_fasta_file");
print OUTPUT_FILE ">$ref_name\n";
print OUTPUT_FILE "$ref_seq\n";

foreach my $record(@read_names){
	print OUTPUT_FILE ">$record\n";
	print OUTPUT_FILE "$all_read_seq{$record}\n";  
}
close OUTPUT_FILE;

print "Done\n";
exit;

############## subroutines #################

##############
# for a given ref sequence position, find the number of reads that cover that position,
##############
sub get_read_coverage{
  my($position) = @_;
  my $out = 0;
  
  # cycle through all reads
  foreach my $rd(keys %all_read_seq){
    
    #get start and end positions for this read
    my $s = $all_read_start{$rd};
    my $e = ($all_read_length{$rd} + $all_read_start{$rd} - 1) ;
    
    if( ($position >= $s) && ($position <= $e) ){
      ## Since the red names includes frequencies, include freq information
   
      $out = $out + 1;
    }
  }
  return $out;
}

####################
# given ch with gaps, get the ch w/o gaps
####################
sub remove_gaps{
  my ($ref) = @_;
  my %out = ();
  
  my %hsh = %$ref;
  my $ct = 0;
  foreach my $p(sort {$a <=> $b} keys %hsh){
    
    my $ch = $hsh{$p};
    if($ch =~ /-/){
      # do nothing
    }
    else{
      $ct = $ct+1;
      $out{$ct} = $ch;
    }
  }
  return %out;
}


###################
#
####################
sub get_rev_arr_values{
  my ($aref) = @_;
  my @out = ();
  
  my @tmp = @$aref;
  
  my $num = scalar @tmp;
  for(my $i = $num-1; $i >= 0; $i--){
    push(@out, $tmp[$i]);
  }
  return @out;
}

###################
#
####################
sub get_rev_complement_for_seq{
  my($str) = @_;
  #print "$str\n";
  my $len = length $str;
  my $rev_str = "";
  for(my $j = $len-1; $j >= 0; $j--){
    
    my $char = substr($str, $j, 1);
    
    #print "J = $j, char = $char\n";
    
    my $comp = "";
    
    if ($char =~ /A/i){
      $comp = 'T';
    }
    elsif($char =~ /G/i){
      $comp = 'C';
    }
    elsif($char =~ /T/i){
      $comp = 'A';
    }
    elsif($char =~ /C/i){
      $comp = 'G';
    }
    elsif($char =~ /\-/){
      $comp = '-';
    }
    $rev_str = "$rev_str"."$comp";
  }
  #print "FWD = $str\n";
  #print "REV = $rev_str\n";
  return $rev_str;
}

###############
# given an array, calculate the average of all values in the aray
###############
sub calc_array_avg{
  my($aref) = @_;
  my $count = 0;
  my $sum = 0.0;
  my $avg = 0.0;
  foreach my $val(@$aref){
    $sum = $sum+$val;
    $count++;
  }
  
  if($count > 0){
    $avg = $sum/$count;
  }
  else{
    $avg = 0.0;
  }
  return $avg;
}




#############
# Given an array, split it into hash, with key being 1 based index
# and value being the array element
#############
sub get_nt_pos{
    my($aref) = @_;
    my %out = ();
    my $count = 1;
    foreach my $rec(@$aref){
	$rec =~ s/\s//g;
	$out{$count} = $rec;
	$count++;
    }
    return %out;
}


