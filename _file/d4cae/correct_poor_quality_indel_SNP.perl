#!/usr/bin/perl

###############
# written by shyamala iyer
# University of Washington, Seattle
# modified: November 15, 2012

## This script does the following:
## Uses a Multiple sequence alignment --> gets a list of all variants, gets the quality values
## of all variants, lists the variant orientation, checks if the variant is present in a homopolymeric
## or non-homopolymeric region, lists the frequency of the variant

## The following variants are flagged for correction: Singleton insertions and deletions, and substitutions
## Variants present only in one direction are also corrected with the consensus character
## This script also tracks which reads have been gap corrected
## Poor quality substitutions are also corrected
## SPECIAL CASES:
## 1) When substitutions, Multi-base pair indels are found in regions of poor coverage, epscially amplicon ends, if the INDEL is part of a multi-base INDEL, keep the
## INDEL intact even if it found only in one orientation, Keep the substitution intact if it is found in area of poor coverage 

## Carry forward errors: This is a special case where: first the reference sequence is checked for all single bp indels near homopolyer regions (after a homopolymer stretch)
## First the positions that can have carry forward errors are listed, at each of these positions, the type of base that would result in carry forward errors is listed
## Carry forward error check is carried out after checking for poor quality bases 
###############

use strict;

my $usage = "perl correct_poor_quality_indel_SNP.pl input_ace_file input_alignment_file orig_read_fasta_file orig_qual_file sample_id fold_coverage multiple_indel_num_check\n";

our $inFile = shift or die $usage;
our $input_alignfile = shift or die $usage;
our $orig_read_file = shift or die $usage;
our $orig_qual_file = shift or die $usage;
our $sample_id = shift or die $usage;
our $fold_coverage = shift;
our $indel_num_check = shift;


our $cutoff = 1;

# this is a value that can be changed, but what it means is that insertions that have greater than 5% occurance
# will not be removed
our $freq_cutoff = $cutoff/100; # variations (undercalls and overcalls) that fall over this range are not stripped

if(length $indel_num_check <= 0){
  $indel_num_check = 3;
} 

if(length $fold_coverage <= 0){
    $fold_coverage = 10;
}
our @in = ();

my @alignment = ();

our @orig_read = ();
our @qual = ();
our @read_names = ();

our %all_read_start = ();  #read name, key == id, value is padded start position of read
our %all_read_seq = (); #just sequence of all reads key== id, value is sequence as a string
our %all_read_orient = (); # store read orientation, fwd reads U rev reads = C
our %read_complement = (); # store complement for reads that are in rev direction, keeping correct leading and trailing gaps
our %read_fwd = ();
our %all_read_seq = ();

our $num_reads = 0;

our $ref_seq = "";
our $ref_name = "";
our $ref_len = 0;

our %all_read_end = ();

our %qual_values = ();
our %orig_read_seq = ();

our %ref_gap_pos = (); ## gap positions in the reference sequence
our %ref_deletion_pos = (); # store all the positions where certain reads have gaps 

our %all_indel_pos = ();
our %hp_indels = (); # positions that fall in homopolymeric regions of the alignment
our %non_hp_indels = (); #positions that dont fall in homopolymeric regions

our %all_sub_pos = (); # all positions with substitutions
our %hp_sub = ();
our %non_hp_sub = ();
our %all_read_len = ();

our %carry_forward_pos = (); ## store positions that could have carry forward errors

our %consensus_char_for_pos = (); # for each of the positions, store the consensus char

our %ref_pos = ();
our %ref_wo_gap = ();

our %qual_scores_for_indels = (); # store all the quality score information
our %high_freq_indels = (); # store the high freq indel values

our %all_sub_scores = (); # qualiy scores for substitutions

our %flagged_indel_pos = (); # positions in the alignment that have only a single read, or only in one orientation
our %flagged_sub_pos = (); # positions with only a single read with substitutions, or substitutions found in one orientation

our %all_pos_freq = (); # Store freq for all non-gap substitutions

our %indel_quality_scores_hp = ();
our %indel_quality_scores_nhp = ();

our %multiple_insert = ();
our %multiple_deletion = ();

our %part_of_multiple_insert = ();
our %part_of_multiple_deletion = ();

our %all_pos_read_coverage = ();

our $max_coverage = 0;
our $min_coverage = 0;

our %poor_coverage_pos = (); ## Store positions that have poor coverage, treat differently

our $log_file = "$sample_id". "_run.log";

my $annotation_file = "$sample_id" . "_sample_annotations.out";


our $total_ins = 0;
our $total_del = 0;
our $total_sub = 0;

our $single_sub = 0;
our $single_indel = 0;
our $one_dir_indel = 0;
our $one_dir_sub = 0;

our $poor_quality_indel = 0;

open (INPUT_FILE, $inFile) or die ("couldn't open $inFile: $!\n");
@in = <INPUT_FILE>;
close INPUT_FILE;

foreach my $line (@in) {
  chomp $line;
  next if $line =~ /^\s*$/; # move to next line if empty
  
  if ($line =~ /^AF\s+(\S+)\s+(\S+)\s+(\S+)$/) {
    my $read_str = $1;
    my $read_UC = $2; #read orientation 
    
    my @cols = split(/\./, $read_str);
    my $read_name = shift @cols;
    if((length $read_name > 0) && (!($read_name =~ /MosaikReference/)) ){
	if($read_name =~ /reversed/){
	    $all_read_orient{$read_name} = 'C';
	}
	else{
	    $all_read_orient{$read_name} = 'U';
	}
	
    } # read not mosaik ref
  }# end of elsif AF
}

# Now get to the sequence part

open (INPUT_FILE, $input_alignfile) or die ("couldn't open $input_alignfile: $!\n");
@alignment = <INPUT_FILE>;
close INPUT_FILE;

open (LOG_FILE, ">>$log_file");
print LOG_FILE "****************************\n";


print LOG_FILE "Base quality check Running\n";

print LOG_FILE "input ace file: $inFile, inputalignment file: $input_alignfile\n read fasta file: $orig_read_file, read base quality file: $orig_qual_file\n sample name: $sample_id, num of consecutive INDELS to check for: $indel_num_check \n";

for( my $i = 0; $i < scalar @alignment; $i++){
  chomp $alignment[$i];
  $alignment[$i] =~ s/\s//g; $alignment[$i] =~ s/\r//g; $alignment[$i] =~ s/\n//g;

  if($alignment[$i] =~ />(\S+)/){
    my $header = $1;
    $header =~ s/\s//g; $header =~ s/\t//g;
    $header =~ s/\r//g; $header =~ s/\n//g;

    if($i== 0){
      $ref_name = $header;
      # deal with the reference sequence
      for(my $k=$i+1; $k < scalar @alignment; $k++){
	chomp $alignment[$k];
	$alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g;
	if($alignment[$k] =~ />(\S+)/){
	  #next header
	  $k = scalar @alignment;
	}
	else{
	  $alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g; $alignment[$k] =~ s/\t//g;
	  $ref_seq = "$ref_seq" . "$alignment[$k]";
	}
      }
      my @ref = split(//, $ref_seq);
      $ref_len = length $ref_seq;
      %ref_pos = get_nt_pos(\@ref); # 1 based index for positions
    }
    else{
      my $curr_read = "";
      my $read_name = $header;
      $read_name =~ s/\s//g; $read_name =~ s/\t//g;
      if($read_name =~ /extra/){
	print LOG_FILE "INCORRECT READ NAME: $read_name\n";
	print LOG_FILE "Program will exit due to errors\n";
	exit;
      }
      if(! exists $all_read_orient{$read_name}){
	print LOG_FILE "INCORRECT READ NAME: $read_name  Read does not exist in the .ace file, $inFile \n";
	print LOG_FILE "Program will exit due to errors\n";

	print "INCORRECT READ NAME: $read_name  Read does not exist in the .ace file, $inFile \n";
	print "Program will exit due to errors\n";
	
	exit;
      }
      else{
	for(my $k=$i+1; $k < scalar @alignment; $k++){
	  chomp $alignment[$k];
	  $alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g;
	  if($alignment[$k] =~ />(\S+)/){
	    #next header
	    $k = scalar @alignment;
	  }
	  else{
	    $alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g; $alignment[$k] =~ s/\t//g;
	    $curr_read = "$curr_read" . "$alignment[$k]";
	  }
	}
	my @read = split(//, $curr_read);
	my ($start, $end) = get_start_end(\@read);
	my $tmp = $curr_read;
	$tmp =~ s/\-//g;

	my $rdlen = length $tmp;

	$all_read_len{$read_name} = $rdlen;

	$all_read_start{$read_name} = $start;
	$all_read_end{$read_name} = $end;
	
	
	if($all_read_orient{$read_name} eq 'C'){
	    $read_complement{$read_name} = $curr_read;
	}
	else{
	    $read_fwd{$read_name} = $curr_read;
	}
	if(exists $all_read_seq{$read_name}){
	  print LOG_FILE "THIS READ $read_name SEEMS TO BE PRESENT IN THIS ALIGNMENT MORE THAN ONCE.\n";
	  print LOG_FILE "Please fix this error and then run the program again. Program will exit.\n";
	  exit;
	}
	else{
	  $all_read_seq{$read_name} = $curr_read;
	  push(@read_names, $read_name);
	}
      }
    }
  }
}

#print "Reference sequence: \n$ref_seq\n";
print LOG_FILE "reference alignment length: $ref_len\n";
print "reference alignment length: $ref_len\n";

my $num = scalar keys %all_read_seq;
print LOG_FILE "length of alignment before starting the base quality check program: $ref_len\n";
print LOG_FILE "Number of read clusters in the alignment file $num\n";

print LOG_FILE "Getting Insertion position positions....\n";

print "Number of read clusters in the alignment file $num\n";

print "Getting Insertion position positions....\n";

my @coverage_vals = (); # store all coverage values

### Get the read coverage for all positions and get positions that have gaps 
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  my $cov = get_read_coverage($p);
  push(@coverage_vals, $cov);
  
  $all_pos_read_coverage{$p} = $cov;
  
  if($ref_pos{$p} eq '-'){
    $ref_gap_pos{$p} = 1;
    #print "Gap found in position $p\n";
  }
}


my @sortarr = sort {$b <=> $a} @coverage_vals;

$max_coverage = shift @sortarr;
$min_coverage = pop @sortarr;

my $num = scalar keys %ref_gap_pos;
print LOG_FILE "Total positions with insertions: $num\n";
print "Total positions with insertions: $num\n";

print "Max read coverage: $max_coverage  Min read coverage: $min_coverage\n";
print LOG_FILE "Max read coverage: $max_coverage  Min read coverage: $min_coverage\n";

$total_ins = $num;

## Now go through all the reads without insertions and looks for reads that have deletions w.r.t reference and put these reads in the file with insertions
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  
  my $refch = $ref_pos{$p};
 
  ## if the ref character is not a gap, 
  if($refch ne '-'){
    # go through all reads and check if that read has a gap at that positions
    foreach my $rd(keys %all_read_seq){
      
      my $st = $all_read_start{$rd};
      my $en = $all_read_end{$rd};
      
      # check if this read falls within the position being checked
      if ( ($p >= $st) && ($p <= $en)){
	my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
	if($ch eq '-'){
	  $ref_deletion_pos{$p} = 1; # some reads at this position in the alignment have deletions
	}
      }
    }
  }
}

$num = scalar keys %ref_deletion_pos;
print "Total positions that have deletions: $num\n";
print LOG_FILE "Total positions that have deletions: $num\n";
$total_del = $num;

%ref_wo_gap = get_ref(\%ref_pos);   ## STORE THE REFERENCE BASED POSITION

##########################################################################
## FIND GAP-ONLY POSITIONS #########
## ignore gap-only pose
my %gap_only_pos = ();
foreach my $p(sort {$a <=> $b} keys %ref_gap_pos){
  my $gap_only = -1;
  foreach my $rd(keys %all_read_seq){
    my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
    if($ch ne '-'){
      $gap_only = 0;
      last; # found a read with a non gap char, move to next gap position
    }
    else{
      $gap_only = 1;
    }
  }
  if($gap_only > 0){
    $gap_only_pos{$p} = 1;
  }
}

############### DEAL WITH SUBSTITUTIONS ##################################
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  if( ! exists $ref_gap_pos{$p}){ ## this position should not be an insertion
    my %char_ct = ();
    foreach my $rd(keys %all_read_seq){
      if( ($p >= $all_read_start{$rd})  && ($p <= $all_read_end{$rd})){
	my $rdch = substr($all_read_seq{$rd}, ($p-1), 1);
	if($rdch ne '-'){
	  if(exists $all_sub_pos{$p}{$rdch}){
	    my $ct = $all_sub_pos{$p}{$rdch};
	    $all_sub_pos{$p}{$rdch} = $ct + 1;
	    $char_ct{$rdch} = $all_sub_pos{$p}{$rdch};
	  }
	  else{
	    $all_sub_pos{$p}{$rdch} = 1;
	    $char_ct{$rdch} = $all_sub_pos{$p}{$rdch};
	  }
	}
      }
    } # done going through all reads for this position
    
    ## Get the consensus character for this position
    my $conchar = get_key_with_highest_val(\%char_ct);
    $consensus_char_for_pos{$p} = $conchar;
  }
}

#### Based on all sub positions, now get the positions that have only a singleton sub ##
my $ct = 0;

foreach my $p(sort {$a <=> $b} keys %all_sub_pos){
  
  foreach my $ch(keys %{$all_sub_pos{$p}}){
    if($ch ne '-'){
    
      if($all_sub_pos{$p}{$ch} == 1){
	#print "position $p char $ch singleton sub\n";
	$flagged_sub_pos{$p}{$ch} = 1;
	$ct++;
      }
    }
  }
}
$single_sub = $ct;

print "Reads with SNPs found only in a single read: $single_sub\n";

foreach my $pos(sort {$a <=> $b} keys %ref_gap_pos){
  if(! exists $gap_only_pos{$pos}){
    $all_indel_pos{$pos} = 'I';
  }
}
foreach my $pos(sort {$a <=> $b} keys %ref_deletion_pos){
  if(exists $ref_gap_pos{$pos}){
    print "Error: the same position $pos has both insertions and deletions\n";
  }
  else{
    $all_indel_pos{$pos} = 'D';
  }
}

if(scalar keys %all_indel_pos == 0){
  print LOG_FILE "please check your alignment... something seems to be incorrect, this alignment seems to have has zero INDELS\n";
  print LOG_FILE "Check if the first sequence in your alignment, the reference sequence is aligned correctly to the rest of the alignment\n";
  exit;
}

################### CHECK POSITIONS THAT HAVE ONLY INDELS in a single read ############

foreach my $p(sort{$a <=> $b} keys %all_indel_pos){
  if($all_indel_pos{$p} eq 'I'){
    my $insert_ct = 0;
    foreach my $rd(keys %all_read_seq){
      my $ch = substr($all_read_seq{$rd},($p-1), 1);
      if($ch ne '-'){
	$insert_ct++;
      }
    }
    if($insert_ct == 1){
      $single_indel++;
      $flagged_indel_pos{$p} = $all_indel_pos{$p};
    }
  }
  else{
    my $del_ct = 0;
    foreach my $rd(keys %all_read_seq){
      my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
      if($ch eq '-'){
	$del_ct++;
      }
    }
    if($del_ct == 1){
      $single_indel++;
      $flagged_indel_pos{$p} = $all_indel_pos{$p};
    }
  }
}
print "Number of indels found only in a single read: $single_indel\n";

## store positions for indels that fall in hp vs nonhp regions

######## for all indels ### check if the indel position falls in a region of homopolymers or not ######
## for a region to be considered hp vs non hp, atleast 2 consensus bases before and after the indel base must be the same
#################################################################################################
#print "Geting weights for the indel positions\n";

foreach my $p(sort {$a <=> $b} keys %all_indel_pos){
  if(! exists $flagged_indel_pos{$p}){
    my $prev1 = get_prev(\%ref_wo_gap, $p);
    my $prev1ch = $ref_pos{$prev1};
    
    my $prev2 = get_prev(\%ref_wo_gap, $prev1);
    my $prev2ch = $ref_pos{$prev2};
    
    my $wt = get_hp_weight($p);
    
    #print "$p position is an indel with weight $wt\n";
    
    
    if($prev1ch ne $prev2ch){
      
      my $after1 = get_next(\%ref_wo_gap, $p);
      my $after1ch = $ref_pos{$after1};
      
      my $after2 = get_next(\%ref_wo_gap, $after1);
      my $after2ch = $ref_pos{$after2};
      #print "position after $p $after1 $after2\n";
      
      if($after1ch eq $after2ch){
	$hp_indels{$p} = $wt;
      }
      else{
	$non_hp_indels{$p} = 1; ## for a non homopolymeric region, the weight is 1
      }
    }
    else{
      $hp_indels{$p} = $wt;
    }
  }
}

#### Check all positions that have substitutions and see if that substitution falls in a homopolymeric or non homopolymeric region
foreach my $p(sort {$a <=> $b} keys %all_sub_pos){
  my $prev1 = get_prev(\%ref_wo_gap, $p);
  my $prev1ch = $ref_pos{$prev1};
  
  my $prev2 = get_prev(\%ref_wo_gap, $prev1);
  my $prev2ch = $ref_pos{$prev2};
  
  if($prev1ch ne $prev2ch){

    my $after1 = get_next(\%ref_wo_gap, $p);
    my $after1ch = $ref_pos{$after1};

    my $after2 = get_next(\%ref_wo_gap, $after1);
    my $after2ch = $ref_pos{$after2};
  
    #print "position after $p $after1 $after2\n";

    if($after1ch eq $after2ch){
      $hp_sub{$p} = 1;
    }
    else{
      $non_hp_sub{$p} = 1; ## for a non homopolymeric region, the weight is 1
    }
  }
  else{
    $hp_sub{$p} = 1;
  }
}
####################################################################################

##################################################################################  

## Calculating substitution frequencies
foreach my $p(sort {$a <=> $b} keys %all_sub_pos){
  
  my $coverage = get_read_coverage($p);

  foreach my $ch(keys %{$all_sub_pos{$p}}){
    
    ## at this point the flagged_sub pos only has substitutions that are present in a single read
    if(! exists $flagged_sub_pos{$p}{$ch}){
      
      if($ch ne '-'){

	my $nt_freq = $all_sub_pos{$p}{$ch}/$coverage;
	$all_pos_freq{$p}{$ch} = $nt_freq;
      }
    }
  }
}

print "Getting BASE QUALITY NUMBERS for each of the INSERTION, DELETION positions.";
print " Depending on the number of positions with gaps(insertions), this might take a while...\n";

my $added = 0;
## Since these files can be big
open(INPUT_FILE, $orig_qual_file);
my $str = "";
my $header = "";

while(my $line = <INPUT_FILE>){
  if($line =~ />(\S+)/){

    ## Starting a new sequence
    $str = "";
    $header = $1;
    #print $header;
    $header =~ s/\s//g; $header =~ s/\t//g;
    
    if(exists $all_read_seq{$header}){
      #print "$header\n";
      
      $added = 1;
    }
    else{
      $added = 0;
    }
  }
  else{
    if($added == 1){
      $str = "$str" . "$line";
      $qual_values{$header} = $str;
    }
  }
}

close INPUT_FILE;
print "Done reading quality values from $orig_qual_file\n";

$added = 0;
open(INPUT_FILE, $orig_read_file);

my $str = "";
my $header = "";

while(my $line = <INPUT_FILE>){
  if($line =~ />(\S+)/){

    ## Starting a new sequence
    $str = "";

    $header = $1;
    #print $header;
    $header =~ s/\s//g; $header =~ s/\t//g;
    
    if(exists $all_read_seq{$header}){
      #print "$header\n";
      $added = 1;
    }
    else{
      $added = 0;
    }
  }
  else{
    if($added == 1){
      $str = "$str" . "$line";
      $orig_read_seq{$header} = $str;
    }
  }
}

close INPUT_FILE;

print "Done reading sequences from $orig_read_file\n";

########################################
# Get INDEL positions that are part of multiple gaps

for(my $i = 1; $i <= scalar keys %ref_wo_gap; $i++){
  
  if( (exists $all_indel_pos{$i}) && (! exists $flagged_indel_pos{$i}) ){
    
    if($all_indel_pos{$i} eq 'I'){
      
      my %ref_start_end = ();
      my $reforig_pos = $ref_wo_gap{$i};
      my @cols = split(/\-/, $reforig_pos);
      my $st = $cols[0];
      my @vals = split(/\./, $cols[1]);
      my $en = $vals[0];

      $ref_start_end{$st}{$en} = 1;
      #print "I\t$i\t$st\t$en\n";

      my $multi_start = $i;
      my $multi_end = 0;
     
      for(my $j = $i+1; $j < scalar keys %ref_wo_gap; $j++){
	
	if( (exists $all_indel_pos{$j}) && ($all_indel_pos{$j} eq 'I') ){
	  my $curr_ref_pos = $ref_wo_gap{$j};
	  my @cols = split(/\-/, $curr_ref_pos);
	  
	  my $currst = $cols[0];
	  
	  if($currst == $st){
	    
	    my @vals = split(/\./, $cols[1]);
	    
	    my $curren = $vals[0];
	    if($curren == $en){
	      ## If this is a consecutive insertion, then the start and ed would be the same
	      if(exists $ref_start_end{$currst}{$curren}){
		my $count = $ref_start_end{$currst}{$curren};
		$count = $count + 1;
		$ref_start_end{$currst}{$curren} = $count;
	      }
	    }
	  }
	}
	
	else{
	  $i = $j;
	  last;
	}
      }
      #print "I\t$i\t$st\t$en\t current count = $ref_start_end{$st}{$en}\n";

      ## now getting the count of how many consecutive insertions are there, 
      if($ref_start_end{$st}{$en} >= $indel_num_check){
	$multi_end = $multi_start + ($ref_start_end{$st}{$en} - 1);
	#print "$multi_start\t$multi_end\n";
	$multiple_insert{$multi_start} = $multi_end;
      }
    }
  }
}

for(my $i = 1; $i <= scalar keys %ref_wo_gap; $i++){
  
  if( (exists $all_indel_pos{$i}) && (! exists $flagged_indel_pos{$i}) ){
    if($all_indel_pos{$i} eq 'D'){
      
      my %ref_start_end = ();
      
      my $reforig_pos = $ref_wo_gap{$i};

      my $multi_start = $i;
      my $multi_end = 0;
      my $ct = 1;
      #print "D\t$i\t$reforig_pos\n";
     
      
      my $prev_pos = $reforig_pos;
      for(my $j = $i+1; $j < scalar keys %ref_wo_gap; $j++){

	if( (exists $all_indel_pos{$j}) && ($all_indel_pos{$j} eq 'D')){
	  my $curr_ref_pos = $ref_wo_gap{$j};
	  
	  if($curr_ref_pos - $prev_pos == 1){
	    $ct = $ct + 1;
	    $prev_pos = $curr_ref_pos;
	    $multi_end = $j;
	  }
	  else{
	    $i = $j;
	    last;
	  }
	}
      }
      #print "D\t$i\t$reforig_pos\t$prev_pos\t current count = $ct\n";
      if($ct >= $indel_num_check){
	#print "$multi_start\t$multi_end\n";
	$multiple_deletion{$multi_start} = $multi_end;
      }
    }
  }
}

foreach my $p(sort {$a <=> $b} keys %all_indel_pos){
 
  if(! exists $flagged_indel_pos{$p}){
  
    if(exists $multiple_insert{$p}){
      my $end = $multiple_insert{$p};
      
      for(my $i = $p; $i <= $end; $i++){
	if( (exists $all_indel_pos{$i}) && ($all_indel_pos{$i} eq 'I')){
	  $part_of_multiple_insert{$i} = 1;
	}
      }
    }
  }
}

print LOG_FILE "Getting regions that the program considers part of multi-base INDELs\n";
print LOG_FILE "The positions are alignment positions in the $input_alignfile\n";
print LOG_FILE "Please check these regions for accuracy\n";

foreach my $p(sort {$a <=> $b} keys %part_of_multiple_insert){
  print "part of multiple insert\t$p\n";
  print LOG_FILE "part of multiple insert\t$p\n";
}

foreach my $p(sort {$a <=> $b} keys %all_indel_pos){

  if(! exists $flagged_indel_pos{$p}){
    
    if(exists $multiple_deletion{$p}){
      my $end = $multiple_deletion{$p};
      for(my $i = $p; $i <= $end; $i++){
	$part_of_multiple_deletion{$i} = 1;
      }
    }
  }
}

foreach my $p(sort {$a <=> $b} keys %part_of_multiple_deletion){
  print "part of multiple deletion\t$p\n";
  print LOG_FILE "part of multiple deletion\t$p\n";
}

############################
%qual_scores_for_indels = get_qualscores_for_indels();

## Get qallity values for substitutions
get_qualscores_for_subs();
#############################

### With these quality score values for all Indels, get the average score values for deletions, Insertions in hp and nhp regions
my $hp_del_avg = 0; my @hp_del = ();
my $nhp_del_avg = 0; my @nhp_del = ();

my $hp_ins_avg = 0; my @hp_ins = ();
my $nhp_ins_avg = 0; my @nhp_ins = ();

foreach my $p(sort {$a <=> $b} keys %indel_quality_scores_hp){
  if($all_indel_pos{$p} eq 'I'){
    push(@hp_ins, $indel_quality_scores_hp{$p});
  }
  else{
    push(@hp_del, $indel_quality_scores_hp{$p});
  }
}

foreach my $p(sort {$a <=> $b} keys %indel_quality_scores_nhp){
  if($all_indel_pos{$p} eq 'I'){
    push(@nhp_ins, $indel_quality_scores_nhp{$p});
  }
  else{
    push(@nhp_del, $indel_quality_scores_nhp{$p});
  }
}

$hp_del_avg = calc_array_avg(\@hp_del);
$nhp_del_avg = calc_array_avg(\@nhp_del);

$hp_ins_avg = calc_array_avg(\@hp_ins);
$nhp_ins_avg = calc_array_avg(\@nhp_ins);

print "average hp insertion q score drop: $hp_ins_avg\n";
print "average nhp insertion q score drop: $nhp_ins_avg\n";

print "average hp del q score drop $hp_del_avg\n";
print "average nhp del q score drop $nhp_del_avg\n";

##### based on these averages, flag the positions that have "below average quality values" ########
###### DONOT FLAG POSITIONS THAT HAVE MULTIPLE BASE INSERTION AND DELTIONS EVEN IF THEY FALL BELOW AVERAGE Q score values #####

foreach my $p(sort {$a <=> $b} keys %all_indel_pos){

  if($all_indel_pos{$p} eq 'I'){
    
    if ( (! exists $part_of_multiple_insert{$p}) && (! exists $flagged_indel_pos{$p}) ){
      
      if($indel_quality_scores_hp{$p} > $hp_ins_avg){
	# flag position
	$poor_quality_indel++;
	$flagged_indel_pos{$p} = $all_indel_pos{$p};
      }
    }
  }
  else{
    if( (! exists $part_of_multiple_deletion{$p}) && (! exists $flagged_indel_pos{$p}) ){
      
      if($indel_quality_scores_hp{$p} > $hp_del_avg){
	# flag position
	$poor_quality_indel++;
	$flagged_indel_pos{$p} = $all_indel_pos{$p};
      }
    }
  }
}


foreach my $p(sort {$a <=> $b} keys %all_indel_pos){
  if($all_indel_pos{$p} eq 'I'){
    
    if( (! exists $part_of_multiple_insert{$p}) && (! exists $flagged_indel_pos{$p}) ){
      if($indel_quality_scores_nhp{$p} > $nhp_ins_avg){

	$poor_quality_indel++;
	$flagged_indel_pos{$p} = $all_indel_pos{$p};
      }
    }
  }
  else{
    if( (! exists $part_of_multiple_deletion{$p}) && (! exists $flagged_indel_pos{$p}) ){
      if($indel_quality_scores_nhp{$p} > $nhp_del_avg){
	
	$poor_quality_indel++;
	$flagged_indel_pos{$p} = $all_indel_pos{$p};
      }
    }
  }
}


#### using the quality scores for all SNP variants, get the distribution average, get the variants that have higher than average quality score reduction
#### flag positions with higher than average score reduction

my $snp_average = 0;
my @all_snp_scores = ();
my $poor_quality_snp = 0;

foreach my $p(sort {$a <=> $b} keys %all_sub_scores){
	foreach my $snpch(keys %{$all_sub_scores{$p}}){
		push(@all_snp_scores, $all_sub_scores{$p}{$snpch});
	}	
}
$snp_average = calc_array_avg(\@all_snp_scores);

## Based on this average distribution, get the SNPs falling above distribution 

foreach my $p(sort {$a <=> $b} keys %all_sub_scores){
	foreach my $snpch(keys %{$all_sub_scores{$p}}){
		if($all_sub_scores{$p}{$snpch} > $snp_average){
			$flagged_sub_pos{$p}{$snpch} = 1;
			$poor_quality_snp++;
		}	
	}	
}


print "Total insertions: $total_ins\n";
print "Total deletions: $total_del\n";
print "total substitutions $total_sub\n";

print "Substitutions found only in a single read: $single_sub\n";
print "Indels found only in a single read: $single_indel\n";

print "Indels  found only in one dir: $one_dir_indel\n";

print "Substitutions found only in one dir: $one_dir_sub\n";

print "Poor quality indels: $poor_quality_indel\n";

print "Poor quality snps: $poor_quality_snp\n";

print LOG_FILE "Total insertions: $total_ins\n";
print LOG_FILE "Total deletions: $total_del\n";
print LOG_FILE "total substitutions $total_sub\n";

print LOG_FILE "Substitutions found only in a single read: $single_sub\n";
print LOG_FILE "Indels found only in a single read: $single_indel\n";

print LOG_FILE "Indels  found only in one dir: $one_dir_indel\n";


print LOG_FILE "Substitutions found only in one dir: $one_dir_sub\n";

print LOG_FILE "Poor quality indels: $poor_quality_indel\n";
print LOG_FILE "Poor quality snps: $poor_quality_snp\n";

print LOG_FILE "Done quality check program\n";

my $flagged_sub_file = "$sample_id" . "_" . "flagged_sub_positions.out";
my $output_stat_file = "$sample_id" . "_" . "variant_base_quality_table.out";
my $flagged_indel_file = "$sample_id" . "_" . "flagged_indel_positions.out";

print LOG_FILE "printing base quaity values in $output_stat_file\nprinting flagged indel positions in $flagged_indel_file\nprinting flagged sub positions in $flagged_sub_file\n";

close LOG_FILE;

##### print the quality score information into the output file

print "Printing position and quality score information to $output_stat_file ...\n";

open(OUTPUT_FILE, ">$output_stat_file");
# print OUTPUT_FILE "Baseline homopolymeric insertion rate = $hp_insert_freq , Baseline homopolymeric deletion rate = $hp_del_freq\n";
# print OUTPUT_FILE "Baseline non homopolymeric insertion rate = $nhp_insert_freq , Baseline non homopolymeric deletion rate = $nhp_del_freq\n";
# print OUTPUT_FILE "Baseline substitution rate for homopolymeric regions = $sub_hp_freq\n";
# print OUTPUT_FILE "Baseline substitution rate for non homopolymeric regions = $sub_nhp_freq\n";

print OUTPUT_FILE "\n-----------------------------------------\n";

print OUTPUT_FILE "Total insertions: $total_ins\n";
print OUTPUT_FILE "Total deletions: $total_del\n";
#print OUTPUT_FILE "total substitutions $total_sub\n";

print OUTPUT_FILE "Singleton indels: $single_indel\n";
print OUTPUT_FILE "Singleton sub: $single_sub\n";

print OUTPUT_FILE "Indels  found only in one dir: $one_dir_indel\n";
print OUTPUT_FILE "Substitutions found only in one dir: $one_dir_sub\n";

print OUTPUT_FILE "Poor quality indels: $poor_quality_indel\n";
#print OUTPUT_FILE "Poor quality sub: $poor_quality_sub\n";

print OUTPUT_FILE "\n\n---------------- INSERTIONS DELETION VALUES ---------------------\n";

print OUTPUT_FILE "Position in alignment\tInsertion OR Deletion\t 1% FREQ\t hp, nhp region\tIndel Freq\tRead_coverage_for_pos\tprev_col_valid reads\tprev_avg Quality value avg\tIndel col valid reads\t #fwd reads\t #rev reads\tIndel col Q value avg\tNext col valid reads\tNext col Q value avg\tQ score diff\t weighted Q score\n";
foreach my $k(keys %qual_scores_for_indels){
  print OUTPUT_FILE "$qual_scores_for_indels{$k}";
}

print OUTPUT_FILE "\n\n-------------- SUBSTITUTION SCORES ---------------\n";
print OUTPUT_FILE "Position\tSub char\t 1 % FREQ\t hp,nhp region \tconsensus char\t Sub Freq\t Read coverage\t #reads with substitution char\t #fwd reads with sub \t #rev reads with sub \tSub q value avg\tConsensus q value avg\t Q score diff\t weighted Q score\n";
foreach my $p(sort {$a <=> $b} keys %all_sub_scores){
  foreach my $ch(keys %{$all_sub_scores{$p}}){
    print OUTPUT_FILE "$all_sub_scores{$p}{$ch}";
  }
}

print OUTPUT_FILE "average hp insertion q score drop: $hp_ins_avg\n";
print OUTPUT_FILE "average nhp insertion q score drop: $nhp_ins_avg\n";

print OUTPUT_FILE "average hp del q score drop $hp_del_avg\n";
print OUTPUT_FILE "average nhp del q score drop $nhp_del_avg\n";

close OUTPUT_FILE;


open(OUTPUT_FILE, ">$flagged_sub_file");
foreach my $p(sort {$a <=> $b} keys %flagged_sub_pos){
  foreach my $ch(keys %{$flagged_sub_pos{$p}}){
    print OUTPUT_FILE "$p\t$ch\n";
  }
}
close OUTPUT_FILE;


open(OUTPUT_FILE, ">$flagged_indel_file");
foreach my $p(sort {$a <=> $b} keys %flagged_indel_pos){
  print OUTPUT_FILE "$p\t$flagged_indel_pos{$p}\n";
}
close OUTPUT_FILE;


#### Using the input alignment file and flagged positions, update the alignment 
## Updates: 
## 1) Insertion columns are removed
## 2) Deletes bases are updated with the consensus char
## 3) Substituted bases are updated with the consensus char

my %annotation = ();
my %high_anno = ();


open(LOG_FILE, ">>$log_file");
print LOG_FILE "updating the alignment to correct flagged INDEL and flagged substitution positions\n";




######### UPDATE DELETIONS ###########

foreach my $p(sort {$a <=> $b} keys %flagged_indel_pos){
  if(exists $ref_pos{$p}){
    
    if( ($ref_pos{$p} ne '-') && ($flagged_indel_pos{$p} eq 'D')){
      ## Reference has a valid character, cycle through all the reads
      ## reads that have a gap in this position should be updated to reflect the char
      ## at the reference positon
      my $ch = $consensus_char_for_pos{$p};
      
      foreach my $rd(keys %all_read_seq){

	my $st = $all_read_start{$rd};
	my $en = $all_read_end{$rd};
	if( ($p >= $st) && ($p <= $en)){

	  my $readch = substr($all_read_seq{$rd}, ($p-1), 1);
	  
	  if($readch eq '-'){
	    
	    my @readseq = split(//, $all_read_seq{$rd});
	    my %read_pos = get_nt_pos(\@readseq); ## 1 based index
	    
	    if($read_pos{$p} eq '-'){
	      $read_pos{$p} = $ch;
	      
	      if(exists $annotation{$rd}){
		my $ant = $annotation{$rd};
		$ant = "$ant" . "," . "$p" . "D";
		$annotation{$rd} = $ant;
	      }
	      else{
		my  $ant = "$p" . "D";
		$annotation{$rd} = $ant;
	      }
	    }
	    my $str = "";
	    foreach my $ct(sort {$a <=> $b} keys %read_pos){
	      $str = "$str" . "$read_pos{$ct}";
	    }
	    $all_read_seq{$rd} = $str;

	    
	  }
	}
      }
    }
  }
  else{
    print "$p is not a valid position\n";
    print "Program will exit.\n";
    #exit;
  }
}
print "All flagged deletions have been updated with the consensus character at that position.\n";
print LOG_FILE "All flagged deletions have been updated with the consensus character at that position.\n";

######## BEFORE UPDATING GAPS, UPDATE SUBSTITUTIONS #####
foreach my $p(sort {$a <=> $b} keys %flagged_sub_pos){

  if(exists $ref_pos{$p}){

    foreach my $sub_char(keys %{$flagged_sub_pos{$p}}){
    
      if($ref_pos{$p} ne '-'){

	## Reference has a valid character, cycle through all the reads
	## reads that have a gap in this position should be updated to reflect the char
	## at the reference positon
	
	my $conch = $consensus_char_for_pos{$p};
	
	foreach my $rd(keys %all_read_seq){
	  my $st = $all_read_start{$rd};
	  my $en = $all_read_end{$rd};
	  
	  if( ($p >= $st) && ($p <= $en)){
	    my $readch = substr($all_read_seq{$rd}, ($p-1), 1);
	    if($readch eq $sub_char){

	 
	      my @readseq = split(//, $all_read_seq{$rd});
	      my %read_pos = get_nt_pos(\@readseq); ## 1 based index
	      if($read_pos{$p} eq $sub_char){
		$read_pos{$p} = $conch; # update with the reference or consensus char

		## Updating character
		if(exists $annotation{$rd}){
		  my $ant = $annotation{$rd};
		  $ant = "$ant" . "," . "$p" . "S";
		  $annotation{$rd} = $ant;
		  
		}
		else{
		  my $ant = "$p" ."S";
		  $annotation{$rd} = $ant;
		}
	      }
	      my $str = "";
	      foreach my $ct(sort {$a <=> $b} keys %read_pos){
		$str = "$str" . "$read_pos{$ct}";
	      }
	      $all_read_seq{$rd} = $str;
	    }
	  }
	}
      }
    }
  }
  else{
    print "$p is not a valid position\n";
  }
}

print "All flagged substitutions have been updated with the consensus character at that position.\n";
print LOG_FILE "All flagged substitutions have been updated with the consensus character at that position.\n";

#### BEFORE UPDATING ALL INSERTIONS CHECK EACH READ AND SEE IF THAT READ HAS A VALID BASE IN THE POSITION THAT NEEDS UPDATING
foreach my $rd(keys %all_read_seq){
  foreach my $p(sort {$a <=> $b} keys %flagged_indel_pos){

    if($flagged_indel_pos{$p} eq 'I'){
      
      my $rdch = substr($all_read_seq{$rd}, ($p-1), 1);
      if($rdch ne '-'){
	if(exists $annotation{$rd}){
	  my $ant = $annotation{$rd};
	  
	  $ant = "$ant" . "," . "$p" ."I";
	  $annotation{$rd} = $ant;
	  
	}
	else{
	  my $ant = "$p" ."I";
	  $annotation{$rd} = $ant;
	}
      }
    }
  }
}


## Remove all Insertion columns that have been flagged for singleton, orientation and base quality
foreach my $rd(keys %all_read_seq){
  my $ct = 0;

  foreach my $p(sort {$a <=> $b} keys %ref_pos){
    if( ( (exists $flagged_indel_pos{$p}) && ($flagged_indel_pos{$p} eq 'I')) || (exists $gap_only_pos{$p})){
      
      #print "Insertion Postion $p being updated\n";
      my $curr_pos = $p-$ct;
      
      my $rdlen = length $all_read_seq{$rd};
      my $prefix = substr($all_read_seq{$rd}, 0, ($curr_pos-1));
      my $suffix = substr($all_read_seq{$rd}, $curr_pos, ($rdlen - $curr_pos));
      my $str = "$prefix" . "$suffix";
      $all_read_seq{$rd} = $str;
      $ct = $ct+1;
    }
  }
}

print "poor quality flagged insertions have been removed from  sequences.\n";
print LOG_FILE "Poor quality flagged insertions have been removed from  sequences.\n";
## do the same for the ref position

my $ct = 0;
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  if( ( (exists $flagged_indel_pos{$p}) && ($flagged_indel_pos{$p} eq 'I')) || (exists $gap_only_pos{$p})){
    my $curr_pos = $p-$ct;
    #print "Reference sequence for position $p being updated\n";
    
    my $curr_ref_len = length $ref_seq; ## get the updated ref seq length

    my $prefix = substr($ref_seq, 0, ($curr_pos-1));
    my $suffix = substr($ref_seq, $curr_pos, ($curr_ref_len - $curr_pos));
    $ref_seq = "$prefix" . "$suffix";
    $ct = $ct+1;
  }
}

my $new_ref_length = length $ref_seq;

print "Check the updated sequences for carry forward errors.\n";
print LOG_FILE "Check the updated sequences for carry forward errors.\n";

############# IDENTIFY POSITIONS WITH POTENTIAL CARRY FORWARD ERRORS #####

## Now with the updated sequences, we need to get an updated list of indels, 
my @ref = split(//, $ref_seq);
%ref_pos = get_nt_pos(\@ref); # 1 based index for positions

my %updated_indel_pos = ();
my %updated_flagged_pos = ();

foreach my $p(sort {$a <=> $b} keys %ref_pos){
  if($ref_pos{$p} eq '-'){
  
    $updated_indel_pos{$p} = 1;
    #print "Gap found in position $p\n";
  }
}


my $num = scalar keys %updated_indel_pos;
print "Total positions with insertions after base quality checks: $num\n";
print LOG_FILE "Total positions with insertions after base quality checks: $num\n";
%ref_wo_gap = get_ref(\%ref_pos);   ## STORE THE REFERENCE BASED POSITION

my %single_base_insert = ();


# 1) get a list of single base insertions
# 2) of these single base insert positions, go back two positions (consensus positions), and check if that position is part of a homopolymer region, if so then get the homopolymer base, and store
# the single base insertion position as a potential carryforward position

foreach my $pos(sort  {$a <=>$b} keys %updated_indel_pos){
  my $prev = $pos -1;
  my $nex = $pos +1;

  if( (! exists $updated_indel_pos{$prev} ) && (! exists $updated_indel_pos{$nex}) ){
    $single_base_insert{$pos} = 1;
    #print "single base insert $pos\n";
  }
}


foreach my $p(sort {$a <=> $b} keys %single_base_insert){
  my $prev1 = get_prev(\%ref_wo_gap, $p);
  my $prev1ch = $ref_pos{$prev1};
  
  my $prev2 = get_prev(\%ref_wo_gap, $prev1);
  my $prev2ch = $ref_pos{$prev2};

  if($prev1ch ne $prev2ch){
    my $prev3 = get_prev(\%ref_wo_gap, $prev2);
    my $prev3ch = $ref_pos{$prev3};
    if($prev3ch eq $prev2ch){
      ## Found a potential carry forward position,
      $carry_forward_pos{$p} = $prev2ch;
      print "Potential carry forward position $p  with char $prev2ch\n";
    }
  }
  ## Similarly do the same for rev reads by looking at the next position
    my $next1 = get_next(\%ref_wo_gap, $p);
    my $next1ch = $ref_pos{$next1};
    my $next2 = get_next(\%ref_wo_gap, $next1);
    my $next2ch = $ref_pos{$next2};

    if($next1ch ne $next2ch){
        my $next3 = get_next(\%ref_wo_gap, $next2);
        my $next3ch = $ref_pos{$next3};
        if($next3ch eq $next2ch){
            ## Found a potential carry forward position,                                                                                          
            $carry_forward_pos{$p} = $next2ch;
            print "Potential carry forward position $p  with char $next2ch\n";
        }
    }


}
###############################################################

## Go through all reads, "correct the reads that have a carry forward error"
my $ct = 0;
my $corr_ct = 0;

foreach my $p(sort {$a <=> $b} keys %carry_forward_pos){
  $ct++;
  foreach my $rd(keys %all_read_seq){

    my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
    if($ch eq $carry_forward_pos{$p}){
      
      $corr_ct = $corr_ct + 1;
      print "read $rd has a carry forward error that was corrected at position $p\n";
      print LOG_FILE "read $rd has a carry forward error that was corrected at position $p\n";
      ## make a new read seq that now has a gap in the carry foward position
      my @readseq = split(//, $all_read_seq{$rd});
      my %read_pos = get_nt_pos(\@readseq); ## 1 based index

      ## chnge this position to a gap
      $read_pos{$p} = '-';

      my $str = "";
      foreach my $ct(sort {$a <=> $b} keys %read_pos){
	$str = "$str" . "$read_pos{$ct}";
      }
      $all_read_seq{$rd} = $str;
    }
  }
}

print "Total positions with potential carry forward errors: $ct\n";
print "After checking all reads, $corr_ct reads have positions that have these carry forward errors\n";

print LOG_FILE "Total positions with potential carry forward errors: $ct\n";
print LOG_FILE "After checking all reads, $corr_ct reads have positions that have these carry forward errors\n";

################### CHECK POSITIONS THAT HAVE ONLY INDELS in a single read ############
$num = 0;

foreach my $p(sort{$a <=> $b} keys %updated_indel_pos){

  my $insert_ct = 0;

  foreach my $rd(keys %all_read_seq){
    my $ch = substr($all_read_seq{$rd},($p-1), 1);
    if($ch ne '-'){
      $insert_ct = $insert_ct + 1;
    }
  }
  if($insert_ct == 1){
    $num++;
    $updated_flagged_pos{$p} = 1;
  }
}
print "Number of indels found only in a single read after removing carry forward errors: $num\n";

### get gap only positions after correcting carry forward errors ####
foreach my $p(sort {$a <=> $b} keys %updated_indel_pos){
 
  my $gap_only = -1;
  
  foreach my $rd(keys %all_read_seq){
    my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
    if($ch ne '-'){
      $gap_only = 0;
      last; # found a read with a non gap char, move to next gap position
    }
    else{
      $gap_only = 1;
    }
  }
  if($gap_only > 0){
    $updated_flagged_pos{$p} = 1;
  }
}



## Correct all other flagged insertion positions ###

foreach my $rd(keys %all_read_seq){
  
  my $ct = 0;

  foreach my $p(sort {$a <=> $b} keys %updated_flagged_pos){

    #print "Insertion Postion $p being updated\n";
    my $curr_pos = $p-$ct;
    
    my $rdlen = length $all_read_seq{$rd};
    my $prefix = substr($all_read_seq{$rd}, 0, ($curr_pos-1));
    my $suffix = substr($all_read_seq{$rd}, $curr_pos, ($rdlen - $curr_pos));
    my $str = "$prefix" . "$suffix";
    $all_read_seq{$rd} = $str;
    $ct = $ct+1;
  }
}

## do the same for the ref position

my $ct = 0;
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  
  if(exists $updated_flagged_pos{$p}){
    my $curr_pos = $p-$ct;
    #print "Reference sequence for position $p being updated\n";
    
    my $curr_ref_len = length $ref_seq; ## get the updated ref seq length
    
    my $prefix = substr($ref_seq, 0, ($curr_pos-1));
    my $suffix = substr($ref_seq, $curr_pos, ($curr_ref_len - $curr_pos));
    $ref_seq = "$prefix" . "$suffix";
    $ct = $ct+1;
  }
}

$new_ref_length = length $ref_seq;



##################### CHECK WHICH READS HAVE MORE THAN 5% UPDATING
foreach my $rd(keys %annotation){

  chomp $annotation{$rd};

  my @cols = split(/,/, $annotation{$rd});
  my $num = 0;
  foreach my $val(@cols){
    if($val =~ /S/){

    }
    else{
      $num = $num+ 1;
    }
  }
  if($num >= ( ($all_read_len{$rd}*5)/100) ){
    print LOG_FILE "$rd has corrected INDELS more than 5% of positions\n";
    print "$rd has corrected INDELS in  more than 5% of positions\n";
    $high_anno{$rd} = $num;
  }
}

my $updated_file = "$sample_id" . "_" . "base_qual_updated_alignment.fasta";
print "Printing the updated alignment to $updated_file\n";
print LOG_FILE "Printing the updated alignment to $updated_file\n";
print LOG_FILE "The alignment length after correcting for singleton, one orientation Indels, substitutions, and poor quality indels in homopolymer regions: $new_ref_length\n";

print "The alignment length after correcting for singleton, one orientation Indels, substitutions, and poor quality indels in homopolymer regions: $new_ref_length\n";
print LOG_FILE "*************************\n";
close LOG_FILE;

##### WRTE ANNOTATIONS TO ANNOTATIONS FILE ####
open(ANNT, ">$annotation_file");
foreach my $rd(keys %annotation){
  print ANNT "$rd\t$annotation{$rd}\n";
}
close ANNT;



########### write updated alignment to output file 
open(OUTPUT_FILE, ">$updated_file");

print OUTPUT_FILE ">$ref_name\n";
print OUTPUT_FILE "$ref_seq\n";
foreach my $rd(@read_names){
  
  if( exists $all_read_seq{$rd} ){
    if(! exists $high_anno{$rd}){
      print OUTPUT_FILE ">$rd\n";
      print OUTPUT_FILE "$all_read_seq{$rd}\n";
    }
  }
}
close OUTPUT_FILE;

print "Done\n";



exit;

############## subroutines #################

###################
# Get read consensus based on given alignment
###################
sub get_key_with_highest_val{
  my ($href1) = @_;
  my %vals = %$href1;

  my %curr_ct = (); # store all freq for this position and their

  foreach my $ch(keys %vals){
  
    my $ct = $vals{$ch};
    $curr_ct{$ct} = $ch;
  }

  my @sorted_cts = sort ({$b <=> $a} keys %curr_ct);
  my $high_ct = shift @sorted_cts;
  my $high_char = $curr_ct{$high_ct};
  return $high_char;
}


#####################
# For a given homopolymer region , get the weights
## For insertion: an average weight for a region is selected
## For deletion, the number of chars that match the deleted base is selected as the weight
#####################
sub get_hp_weight{
  my ($curr_pos) = @_;

  my $out = 0; ## the weight of the homopolymer region

  my $type = $all_indel_pos{$curr_pos};
  
  if($type eq 'I'){

    ## for an insertion sometimes alignment is not perfect hence similar charcters both before and after the 
    
    my $prevpos = get_prev(\%ref_wo_gap, $curr_pos);
    my $prev_refch = $ref_pos{$prevpos};

    my $prevch_ct = 0;
    
    for(my $i = $prevpos; $i >= 1; $i--){
      my $prev = get_prev(\%ref_wo_gap, $i);
      my $prevch = $ref_pos{$i};
      if($prevch eq $prev_refch){
	$prevch_ct++;
      }
      else{
	last; ## exit loop
      }
    }
    
    my $after_ch_ct = 0;

    my $after_pos = get_next(\%ref_wo_gap, $curr_pos);
    my $after_ref_ch = $ref_pos{$after_pos};
    
    for(my $i = $after_pos; $i <= $ref_len; $i++){
      my $after = get_next(\%ref_wo_gap, $i);
      my $afterch = $ref_pos{$i};
      if($afterch eq $after_ref_ch){
	$after_ch_ct++;
      }
      else{
	last; ## exit loop
      }
    }

    ## for an insertion, return the greater weight
    if($prevch_ct > $after_ch_ct){
      $out = $prevch_ct;
    }
    elsif($prevch_ct <= $after_ch_ct){
      $out = $after_ch_ct;
    }
  }
  else{

    ## for deletion need to look only at the previous characters, since by delafult the alignment should have been done so that the undercall
    ## is the last base in the hp region
    my $refch = $ref_pos{$curr_pos};
    
    my $char_ct = 0;
    for(my $i = $curr_pos; $i >= 1; $i--){
      my $prev = get_prev(\%ref_wo_gap, $i);
      my $prevch = $ref_pos{$i};
      if($prevch eq $refch){
	$char_ct++;
      }
      else{
	$out = $char_ct;
	last;
      }
    }
  }
  return $out;
}


###################
# 
###################
sub get_prev{
  my($ref, $pos) = @_;
  my %hsh = %$ref;

  my $out = 0;

  for(my $i = $pos-1; $i > 0; $i--){
    my $refconpos = $hsh{$i};
    my @cols = split(/\-/, $refconpos);
    if(scalar @cols == 1){
      $out = $i;
      last;
    }
  }
  return $out;
}


##################
#
##################
sub get_next{
  my($ref, $pos) = @_;
  my %hsh = %$ref;
  my $out = 0;
  for(my $i = $pos+1; $i < scalar keys %hsh; $i++){
    my $refconpos = $hsh{$i};
    my @cols = split(/\-/, $refconpos);
    if(scalar @cols ==1){
      $out = $i;
      last;
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
# Given an alignment, get the positions with substitutions and calculate scores for each of these
## substitutions
##################
sub get_qualscores_for_subs{
  
  foreach my $p(sort {$a <=> $b} keys %all_sub_pos){
    
    my @consensus_scores = ();
    my %other_scores = ();
    my %qscore_diff = (); # store difference between all non-con characters and avg scores for con characters
    my $coverage = get_read_coverage($p);
    
    #print "Read coverage $p = $coverage\n";
    # Get consensus char for this
    # Go through all the scoring ONLY IF THERE IS MORE THAN 1 CHARACTER AT A POSITION
    
    
    ## for this positionn , get the number of reads in fwd and rev directions for each of the substitutions
    my %char_fwd = ();
    my %char_rev = ();
    
    foreach my $ch(keys %{$all_sub_pos{$p}}){
      # for this char, at this position, track the number of fwd and rev reads that have the substitution 

      my $num_fwd = 0;
      my $num_rev = 0;

      if( ($all_pos_freq{$p}{$ch} > 0) && ($all_pos_freq{$p}{$ch} < 1) ){
	
	#print "Calculating substitution score for $p char $ch\n";
	if($ch eq $consensus_char_for_pos{$p}){
	  #print "$p consensus char $ch\n";
	  foreach my $rd(keys %all_read_seq){
	    my $st = $all_read_start{$rd};
	    my $en = $all_read_end{$rd};

	    my $freq = 1;

	    if ( ($p >= $st) && ($p <= $en)){
	      my $rdch = substr($all_read_seq{$rd}, ($p-1), 1);
	      if($rdch eq $ch){
		
		# get the orig read sequence 
		my $orig_rd_seq = $orig_read_seq{$rd}; # will match the complemented read if rev reads
		my $qlstr = ""; ## make sure the quality scores match correctly
		$qlstr = $qual_values{$rd};
		
		my $score = match_char_score($rd, $orig_rd_seq, $qlstr, $p);
		#print "consensus $rd $ch $p score = $score\n";
		#print "$rd $ch $p score = $score\n";
		if($score >= 0){
		  for(my $i = 0; $i < $freq; $i++){
		    push(@consensus_scores, $score);
		  }
		}
		elsif($score < 0){
		  print "Negative score for $ch at position $p, This read sequence has been altered\n";
		}
	      }
	    }
	  }
	}
	else{
	  # Char not a consensus score
	  #print "char $ch not consensus\n";
	  my @tmp = (); # array to store the scores of non-consensus characters
	 
	  ## go through all reads, get scorea from all reads that have this character at this position
	  foreach my $rd(keys %all_read_seq){
	    my $st = $all_read_start{$rd};
	    my $en = $all_read_end{$rd};

	    my $freq = 1;

	    if ( ($p >= $st) && ($p <= $en)){

	      my $rdch = substr($all_read_seq{$rd}, ($p-1), 1);
	      #print "$rdch  $ch\n";

	      if($rdch eq $ch){ ## Pick the reads that have this character at that position since we want scores from these particular reads
		
		#print "\n";
		# get the orig read sequence 
		my $orig_rd_seq = $orig_read_seq{$rd}; # will match the complemented read if rev reads
		my $qlstr = ""; ## make sure the quality scores match correctly
		$qlstr = $qual_values{$rd};
		
		my $score = match_char_score($rd, $orig_rd_seq, $qlstr, $p);
		print "$rd at position $p has a substitution $ch with quality score = $score\n";
		if($score >= 0){
		  for(my $i = 0; $i < $freq; $i++){
		    push(@tmp, $score);
		  }
		}
		elsif($score < 0){
		  print "Negative score for $ch at position $p, This read sequence has been altered\n";
		}

		if(exists $read_fwd{$rd}){
		  $num_fwd = $num_fwd + 1;
		}
		elsif(exists $read_complement{$rd}){
		  $num_rev = $num_rev + 1;
		}
	      }
	    }
	  }
	  $other_scores{$ch} = \@tmp;

	  $char_fwd{$ch} = $num_fwd;
	  $char_rev{$ch} = $num_rev;
	  
	  if(($char_fwd{$ch} == 0) || ($char_rev{$ch} == 0)){
	    
	    $one_dir_sub++;
	    ## Check the frequency of this substitution ONLY IF THE COVERAGE OF THIS REGION FALLS TO LESS THAN 1% OF THE MAX COVERAGE
	    my $curr_coverage = get_read_coverage($p);
	    print "$p\t$curr_coverage\n";
	    ## OVER 10 FOLD DIFFERENCE BETWEEN MAX COVERAGE AND CURRENT READ COVERAGE
	    if($max_coverage/$curr_coverage >= $fold_coverage ){
 	      ## the coverage in this region is too low to make a call, keep the substitution, dont flag this position
 	      print "Substitution character $ch at $p is present in a region of very poor read coverage\n";
 	      print "This substitution is present only in ONE orientation. This SUBSTITUTION IS BEING RETAINED\n";

	      print LOG_FILE "Substitution character $ch at $p is present in a region of very poor read coverage\n";
 	      print LOG_FILE "This substitution is present only in ONE orientation. This SUBSTITUTION IS BEING RETAINED\n";
 	    }
 	    else{
 	      ## The coverage at this region is high enough to include both fwd and rev reads,
 	      ## Flag this position for correction
 	      $flagged_sub_pos{$p}{$ch} = 1;
 	    }
	  }
	}
      }
    }# done going through all the characters at this position
    
    # Find score averages, score differences, substitution scores
    my $con_avg = calc_array_avg(\@consensus_scores);
    
    if(scalar keys %other_scores > 0){
      
      foreach my $char(keys %other_scores){
	my $aref = $other_scores{$char};
	my $sub_num = scalar @$aref;
	if($sub_num > 1){ 
	  my $avg = calc_array_avg($aref);
	  $qscore_diff{$char} = abs($con_avg - $avg);

	  if(exists $hp_sub{$p}){
	    if(exists $all_sub_scores{$p}{$char}){
	      print "Score already exists for this $p and $char\n";
	    }
	    else{
	      if($all_pos_freq{$p}{$char} > $freq_cutoff){
		$all_sub_scores{$p}{$char} = "$p\t$char\t >1 % FREQ\t homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
		print "$p\t$char\t >1 % FREQ\t homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
	      }
	      else{
		$all_sub_scores{$p}{$char} = "$p\t$char\t <1 % FREQ\t homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
		print "$p\t$char\t <1 % FREQ\t homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
	      }
	    }
	  }
	  elsif(exists $non_hp_sub{$p}){
	    if(exists $all_sub_scores{$p}{$char}){
	      print "Score already exists for this $p and $char\n";
	    }
	    else{
	      if($all_pos_freq{$p}{$char} > $freq_cutoff){
		$all_sub_scores{$p}{$char} = "$p\t$char\t >1 % FREQ\t non homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
		print "$p\t$char\t >1 % FREQ\t non homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
	      }
	      else{
		$all_sub_scores{$p}{$char} = "$p\t$char\t <1 % FREQ\t non homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
		print "$p\t$char\t <1 % FREQ\t non homopolymeric region \t$consensus_char_for_pos{$p}\t$all_pos_freq{$p}{$char}\t$coverage\t$sub_num\t$char_fwd{$char}\t$char_rev{$char}\t$avg\t$con_avg\t$qscore_diff{$char}\n";
	      }
	    }
	  }
	}
      } ## done going through all characters for the substitutions in this position
    }
  }# Done going through all positions
  return;
}

######################
# Given a position, original read, read id, original quality scores, get the correct score for the character
######################
sub match_char_score{
  my ($rd_id, $rd, $ql, $pos) = @_;
  
  my %corrected_read_start = (); # this is the corrected read positions
  my %pos_char_score = ();;
  
  my $out = 0;
  my $corr_end = 0;
  my $corr_start = 0;
  
  my $corr_pos = 0;
  
  ############### the given alignment, start and end positions correspond to the full alignment and not the start and end positions of that individual read
  my $rdstart = $all_read_start{$rd_id};
  my $rdend = $all_read_end{$rd_id};
  
  # we need to check if the read is fwd or rev orientation and for this we need the actual read sequence without 
  ## leading or trailing dashes
  my $curr_read = "";
  if(exists $read_complement{$rd_id}){
    $curr_read = $read_complement{$rd_id};
    $curr_read =~ s/\-//g;
    #print "rev\n";
  }
  else{
    $curr_read = $read_fwd{$rd_id};
    $curr_read =~ s/\-//g;
    #print "fwd\n";
  }
  
  ## IF THIS READ EXISTS IN THE ORIGINAL SEQUENCE FILE, GET THE CORRECT POSITIONS THAT MATCH TO THE ORIGINAL SEQUENCE
  ## return the first instance of this substr from the original string
  ## if the read is in rev orientation, then we need to compare th rev complement to the orig read seq to get correct "start" position for this read
  ## this "start" position is 0 index based and the position is w.r.t to the original read sequence in the sample fasta file and will determine the quality scores
  
  $rd =~ s/\r//g;
  $rd =~ s/\n//g;
  $curr_read =~ s/\r//g;
  $curr_read =~ s/\n//g;

  my $res = index($rd, $curr_read);
  print "Found this read $rd_id amongst the original read and qulity files $res\n";

  my $str_len = length $curr_read;
  # res has to be a positive number or 0
  if($res >= 0){
    
    ## NOW GET THE CORRECT QUALITY SCORES CORRESPONDING TO THE ORIG SEQUENCE POSITION INFORMATION
    ## These are the condensed original sequence and original qual score values that start at the first read position
    #print "Found $res\n";
    
    my @qual_vals = split(/\s/, $ql);  ### split qual values based on gap, for convenience to add to string
    my @cond_qual_vals = ();
    for(my $i = $res; $i < ($res+$str_len); $i++){
      push(@cond_qual_vals, $qual_vals[$i]);
    }

    ## The quality has already been reversed, for us, no need to reverse it
    #if(exists $read_complement{$rd_id}){
      ## READ IN REV ORIENTATION
      ## Quality scores need to be dealt with seperately in this case
     # my @revquals = get_rev_arr_values(\@cond_qual_vals);
     # @cond_qual_vals = @revquals;
    #}

    my %qualvals = get_nt_pos(\@cond_qual_vals); ## 1 based index, starting at actual read start
    my $curr_rd_str = $all_read_seq{$rd_id}; ## ALIGNED READ THAT HAS A DEF START AND END POSITION
    my $ct = 1;
    for(my $i = $rdstart; $i <= $rdend; $i++){
      my $ch = substr($curr_rd_str, ($i-1), 1);
      $corrected_read_start{$ct} = $ch;
      $ct++;
    }
    
    # Now we have to figure out how far along the alignment the start position was and where do the column, start and end position figure out along the alignment
    $ct = 1;
    foreach my $p(sort {$a <=> $b} keys %corrected_read_start){
      my $ch = $corrected_read_start{$p};
      if($ch ne '-'){
	#print "pos $p\t$ct\t$ch\t$qualvals{$ct}\n";
	$pos_char_score{$p}{$ch} = $qualvals{$ct};
	$ct++;
      }
      else{
	$pos_char_score{$p}{$ch} = -1;
      }
    }
    
    $ct = 1;
    my $ch = "";
    # the actual read sequence and corresponding quality scores start at this position
    $corr_pos = ($pos - $rdstart + 1);
    $ch = $corrected_read_start{$corr_pos};
    $out = $pos_char_score{$corr_pos}{$ch};
  }
  else{
    print "Read $rd_id with the sequence \n$curr_read\n was not found in the original fasta file $orig_read_file\n";
  }
  return $out;
}

####################
# Given an alignment sequence, original read sequence, quality scores for the orig read
# since the alignment could trim beginning and ends of reads, the correct position of the 
# aligned sequence relative to original should be established
####################
sub get_qualscores_for_indels{

  my %out = (); # store quality score information for indels
 
  ## Go through all the gap positions in the alignment
  ## go through all the reads, calculate q score stats for any reads that have a valid char in this position
  #print " processing indels in alignment\n";

  foreach my $p(sort {$a <=> $b} keys %all_indel_pos){
    
    if(! exists $flagged_indel_pos{$p}){
      my $num_fwd = 0;
      my $num_rev = 0;
      
      my $type = $all_indel_pos{$p};
      my $region = "";
      if(exists $hp_indels{$p}){
	$region = 'homopolymer region';
      }
      else{
	$region = 'non homopolymer region';
      }
      
      my %read_match = ();
      my $indel_freq = 0.0;
      
      ## For each position that has an indel, collect qscore values for column, and previous, after columns
      ## for this position get all the reads that have a valid char at this position
      ## if the previous position is also a gap, then get the average quality score of all valid chars in that column
      my @scores = ();
      my @prev = ();
      my @after = ();
      
      my ($prev, $after) = get_flanking_consensus_pos($p);

      print "Scores for these positions are needed: $p\t, previous:  $prev, after: $after\n";

      ## Go through all the reads, get sequence, associated quality scores for the reads
      ## All reads that are cover this position will be used
      foreach my $rd(keys %all_read_seq){
	
	my $bases = 0;
	my $st = $all_read_start{$rd};
	my $en = $all_read_end{$rd};
	
	if ( ($p >= $st) && ($p <= $en)){
	
	
	  # get the orig read sequence 
	  my $orig_rd_seq = $orig_read_seq{$rd}; # will match the complemented read if rev reads
	  my $qlstr = ""; ## make sure the quality scores match correctly
	  $qlstr = $qual_values{$rd};
	  
	  my %match_aligned_seq_to_scores = match_seq_score($rd, $orig_rd_seq, $qlstr, $prev, $p, $after);
	  
	  if (scalar keys %match_aligned_seq_to_scores > 0){
	    $read_match{$rd} = \%match_aligned_seq_to_scores;
	  }
	  
	  if($type eq 'I'){
	    ## Insertion
	    # does  this read have a valid non gap base at this position, if check whether the read is fwd or rev orientation
	    my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
	    if($ch ne '-'){
	      
	      if(exists $read_fwd{$rd}){
		$num_fwd = $num_fwd + 1;
	      }
	      if(exists $read_complement{$rd}){
		$num_rev = $num_rev + 1;
	      }
	      
	    }
	  }
	  else{
	    my $ch = substr($all_read_seq{$rd}, ($p-1), 1);
	    if($ch eq '-'){
	      if(exists $read_fwd{$rd}){
		$num_fwd = $num_fwd + 1;
	      }
	      if(exists $read_complement{$rd}){
		$num_rev = $num_rev + 1;
	      }
	    }
	  }
	}
      } ## done going through all reads
     
      ## At the end of going through all reads, if the Indel is present only in one orientation
      ## add to flagged position
      if (($num_fwd == 0) || ($num_rev == 0)){
	$one_dir_indel++;
      }
      ###########################################
      ## get read coverage for this gap position, this means all reads that span this region will be included in read coverage, though only
      ## reads that have a valid char in that gap position are used for q score calculation
      my $read_coverage_for_pos = get_read_coverage($p);
      
      my $valid_base_ct = get_valid_char_ct($p); ## for a given position, get the valid bases
      
      ## for each of the reads that have a valid char in the gap position, get qscore stats
      
      if(scalar keys %read_match > 0){
	foreach my $rd(keys %read_match){
	  my $ref = $read_match{$rd};
	  
	  my %hsh =%$ref;
	  my $freq = 1;
	  
	  ## each of the reads has three positions: , before variant column, var column and after variant column
	  ## Only non-gap characters will have quality scores associated and can be included in counts and calculating 
	  ## average quality scores
	  
	  my $ct = 1;
	  foreach my $ch(keys %{$hsh{$ct}}){
	    if($hsh{$ct}{$ch} >= 0){
	      
	      for(my $i = 0; $i < $freq; $i++){
		push(@prev, $hsh{$ct}{$ch});
	      }
	    }
	  }
	  $ct = 2;
	  foreach my $ch(keys %{$hsh{$ct}}){
	    if($hsh{$ct}{$ch} >= 0){
	      for(my $i = 0; $i < $freq; $i++){
		push(@scores, $hsh{$ct}{$ch});
	      }
	    }
	  }
	  $ct = 3;
	  foreach my $ch(keys %{$hsh{$ct}}){
	    if($hsh{$ct}{$ch} >= 0){
	      for(my $i = 0; $i < $freq; $i++){
		push(@after, $hsh{$ct}{$ch});
	      }
	    }
	  }
	} ## gone through all reads
	
	# calculate average of the 3 columns, check if the var column has reduced quality compared to previous and subsequenct
	# columns
	my $out_str = "";
	my $num_scores_in_col = scalar @scores;
	#print "$p\n @scores\n";
	
	if($type eq 'I'){
	  if($read_coverage_for_pos > 0){
	    $indel_freq = $num_scores_in_col/$read_coverage_for_pos;
	  }
	  else{
	    $indel_freq = 0.0;
	  }
	}
	elsif($type eq 'D'){
	  if($read_coverage_for_pos > 0){
	    $indel_freq = ($read_coverage_for_pos - $num_scores_in_col)/$read_coverage_for_pos;
	  }
	  else{
	    $indel_freq = 0.0;
	  }
	}
	my @sorted = sort{$a <=> $b} @prev;
	my $prev_avg = calc_array_avg(\@prev);
	my $prev_num = scalar @prev;
	@sorted = sort{$a <=> $b} @scores;
	my $col_avg = calc_array_avg(\@sorted);
	my $col_num = scalar @scores;
	@sorted = sort{$a <=> $b} @after;
	my $after_avg = calc_array_avg(\@after);
	my $after_num = scalar @after;
	my $valid_bases = 0;
	
	## This means that valid bases actually refer to bases where quality scores could have been found
	## 
	if($type eq 'I'){
	  $valid_bases = $col_num;
	}
	elsif($type eq 'D'){
	  $valid_bases = ($read_coverage_for_pos - $col_num);
	}
	my $diff1 = abs($prev_avg - $col_avg);
	my $diff2 = abs($col_avg - $after_avg);
	
	my @tmp = ();
	push(@tmp, $diff1);
	push(@tmp, $diff2);
	
	my $diff = calc_array_avg(\@tmp);
	
	#my $align_freq = 0.0;
	## A SCORING SYSTEM FOR INDELS
	if ($indel_freq < $freq_cutoff){
	  $out_str = "$p\t$type\t <1 % FREQ\t$region\t$indel_freq\t$read_coverage_for_pos\t$prev_num\t$prev_avg\t$valid_bases\t$num_fwd\t$num_rev\t$col_avg\t$after_num\t$after_avg\t$diff\n";
	  print "$out_str\n";
	}
	else{
	    
	  $out_str = "$p\t$type\t >1 % FREQ\t$region\t$indel_freq\t$read_coverage_for_pos\t$prev_num\t$prev_avg\t$valid_bases\t$num_fwd\t$num_rev\t$col_avg\t$after_num\t$after_avg\t$diff\n";
	  print "$out_str\n";
	}
	$out{$p} = $out_str;
	if(exists $hp_indels{$p}){
	  $indel_quality_scores_hp{$p} = $diff;
	}
	else{
	  $indel_quality_scores_nhp{$p} = $diff;
	}
      }
    }
  } ## move to next position
  return %out;
}

#####################
# Given a position, cycle through all reads that cover that position and get the count of valid bases
#####################
sub get_valid_char_ct{
  my($currpos) = @_;
  
   my $out = 0;
  
  # cycle through all reads
  foreach my $rd(keys %all_read_seq){
    
    #get start and end positions for this read
    my $s = $all_read_start{$rd};
    my $e = $all_read_end{$rd};
    
    if( ($currpos >= $s) && ($currpos <= $e) ){
      ## Since the red names includes frequencies, include freq information
      # my @vals = split(/\_/, $rd);
#       my $rdfreq = pop @vals;
      my $char = substr($all_read_seq{$rd}, ($currpos-1), 1);
      if($char ne '-'){
	$out = $out + 1;
      }
    }
  }
  return $out;
}

###################
# given an alignment position, get the previous consensus position and the next consensus position
###################
sub get_flanking_consensus_pos{
  my($alignpos) = @_;

  my $s = 0;
  my $e = 0;

  my $rp = $ref_wo_gap{$alignpos};
 

  my @cols = split(/\-/, $rp);
 
  if(scalar @cols == 1){
    # deletion, deal with it
    # aligment pos is the same,
    for(my $i = $alignpos-1; $i > 0; $i--){
      $rp = $ref_wo_gap{$i};
      my @arr = split(/\-/, $rp);
      if(scalar @arr == 1){
	$s = $i;
	last;
      }
    }
    
    for(my $i = $alignpos+1; $i < scalar keys %ref_wo_gap; $i++){
      $rp = $ref_wo_gap{$i};
      my @arr = split(/\-/, $rp);
      if(scalar @arr == 1){
	$e = $i;
	last;
      }
    }
  }
  else{
    ## insertion column in alignment
    my @tmp = split(/\./, $cols[1]);

    my $r1 = $cols[0];
    my $r2 = $tmp[0];


    foreach my $rec(sort {$a <=> $b} keys %ref_wo_gap){
      
      my $val = $ref_wo_gap{$rec};
      if($val eq $r1){
	$s = $rec;
	last;
      }
    }
    foreach my $rec(sort {$a <=> $b}  keys %ref_wo_gap){
      my $val =  $ref_wo_gap{$rec};
      if($val eq $r2){
	$e =$rec;
      }
    }
  }
  return($s, $e);
}

####################
# script to check for reference region: does the indel fall in homopolymeric or nonhomopolymeric region
####################
sub check_hp_region{
  my ($str, $position) = @_;
  my $hp_flag = 0;
  
  my $prevch = "";
  my $afterch = "";
  my $ch = substr($str, ($position-1), 1);

  
  for(my $i = $position-1; $i > 0; $i--){
    $prevch = substr($str, ($i-1), 1);
    if ($prevch ne '-'){
      last;
    }
  }
  
  for(my $i = $position+1; $i < length $str; $i++){
    $afterch = substr($str, ($i-1), 1);
    if($afterch ne '-'){
      last;
    }
  }
  #print "prev char $prevchar next char $afterch char $ch\n";
  if( ($ch eq $prevch) || ($ch eq $afterch)){
    $hp_flag = 1;
  }
  else{
    $hp_flag = 0;
  }
  return $hp_flag;
}

###################
# Get the nt characters in ref sequence 
####################
sub get_ref{
  my($href) = @_;
  my %align_pos = %$href;
  my %out = ();
  
  my @sorted_pos = sort {$a <=> $b} keys %align_pos;
  my $align_st = shift @sorted_pos;
  my $align_end = pop @sorted_pos;
  
  print "align start $align_st  align end $align_end\n"; 
  my $ct = 1; # count for orig pos

  #####################################
  my $insertion_count = 0;
  for(my $k = $align_st; $k <= $align_end; $k++){

    my $char = $align_pos{$k};
    $char =~ s/\s//g;
    if($char ne '-'){
      # not a gap, process
      $insertion_count = 0; # reset the insertion counter;
      if($ct <= 1){
	$out{$k} = "1";
	$ct = $ct+1;
      }
      else{
	$out{$k} = "$ct";
	$ct = $ct+1;
      }
    }
    else{
      # place in alignment with possible insertions
      my $st_pos = "";
      my $end_pos = "";

      # increment the insertion counter for each position where there is a gap
      $insertion_count = $insertion_count + 1;

      # find start of insertion
      for(my $i = $k; $i >= $align_st; $i--){
	if($i == 1){
	  $st_pos = "1";
	}
	else{
	  if($align_pos{$i} ne '-'){
	    if(exists $out{$i}){
	      $st_pos = $out{$i};
	      last;
	    }
	  }
	}
      }
      for(my $i = $k; $i <= $align_end; $i++){
	if($align_pos{$i} ne '-'){
	  # found the end point
	  $end_pos = "$ct";
	  last;	
	}
      }
      $out{$k} = "$st_pos" . "-" . "$end_pos" . "." . "$insertion_count";
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

####################
# given a read with positions and quality scores, get the correct position that this read belongs to
# 
####################
sub match_seq_score{
  my ($rd_id, $rd, $ql, $beg ,$pos, $end) = @_;

  my %corrected_read_start = (); # this is the corrected read positions
  my %pos_char_score = ();

  # if a read has a valid non gap character at the given position pos, then move to the previous non gap character in that read to get the prev score
  # do the same for character after gap column

  # if a read has a gap character in the given position, then q score for that position = -1, but the prev and subsequenct

  my %out = ();

  ## This is the corrected poitions w.r.t to the read in the alignment 
  my $corr_start = 0;
  my $corr_pos = 0;
  my $corr_end = 0;

  ############### the given alignment, start and end positions correspond to the full alignment and not the start and end positions of that individual read
  my $rdstart = $all_read_start{$rd_id};
  my $rdend = $all_read_end{$rd_id};
  
  #print "$rd_id ";
  # we need to check if the read is fwd or rev orientation and for this we need the actual read sequence without 
  ## leading or trailing dashes
  
  my $curr_read = "";
  if(exists $read_complement{$rd_id}){
    $curr_read = $read_complement{$rd_id};
    $curr_read =~ s/\-//g;
   #  print "rev\t";
  }
  else{
    $curr_read = $read_fwd{$rd_id};
    $curr_read =~ s/\-//g;
    # print "fwd\t";
  }

  ## IF THIS READ EXISTS IN THE ORIGINAL SEQUENCE FILE, GET THE CORRECT POSITIONS THAT MATCH TO THE ORIGINAL SEQUENCE
  ## return the first instance of this substr from the original string
  ## if the read is in rev orientation, then we need to compare th rev complement to the orig read seq to get correct "start" position for this read
  ## this "start" position is 0 index based and the position is w.r.t to the original read sequence in the sample fasta file and will determine the quality scores
  $rd =~ s/\r//g;
  $rd =~ s/\n//g;
  $rd =~ s/\s//g;

  $curr_read =~ s/\r//g;
  $curr_read =~ s/\n//g;

  my $res = index($rd, $curr_read);
  
  #print "current read: $curr_read\n";
  #print "Original read: $rd\n";
  print "Found this read $rd_id amongst the original read and qulity files $res\n";
  my $str_len = length $curr_read;

  # res has to be a positive number or 0
  if($res >= 0){
   
    ## NOW GET THE CORRECT QUALITY SCORES CORRESPONDING TO THE ORIG SEQUENCE POSITION INFORMATION
    ## These are the condensed original sequence and original qual score values that start at the first read position
    #print "Found $res\n";

    my @qual_vals = split(/\s/, $ql);  ### split qual values based on gap, for convenience to add to string

    my $qual_len = scalar @qual_vals;
    my $rd_len = length $rd;

    my @cond_qual_vals = ();

    
    ## when the quality length is not equal to the read length
    ## we need to take care of this
    ## for forward reads:
    # res is the starting point of the string match
    ## the end of teh string match will be the end of 
    ## for rev reads
    
     if($res+$str_len > scalar @qual_vals){
      
       print "$rd_id original qual length not equal to read length\n"; 
       my $difflen = ($res+$str_len - scalar @qual_vals);
      
       for(my $i = $res; $i < scalar @qual_vals; $i++){
 	push(@cond_qual_vals, $qual_vals[$i]);
       }
       for(my $j = 0; $j < $difflen; $j++){
 	push(@cond_qual_vals, 0);
	
       }
     }
    else{
     #  print "$rd_id\n";
      for(my $i = $res; $i < ($res+$str_len); $i++){
	push(@cond_qual_vals, $qual_vals[$i]);
      }
      
      ## debugging:
   #    foreach my $v(@cond_qual_vals){
# 	print "$v ";
#       }
    }
    
    ## Since the amplicon trimmed files have already been reversed for both seq and quality,
    ## no need to reverse quality
    
    #if(exists $read_complement{$rd_id}){
      ## READ IN REV ORIENTATION
      ## Quality scores need to be dealt with seperately in this case
     # my @revquals = get_rev_arr_values(\@cond_qual_vals);
     # @cond_qual_vals = @revquals;

       ## debugging:
    #   print "Rev completement\n";
#       foreach my $v(@cond_qual_vals){
# 	print "$v ";
#       }
#       print "\n";

    #}
    my %qualvals = get_nt_pos(\@cond_qual_vals); ## 1 based index, starting at actual read start

    
    my $curr_rd_str = $all_read_seq{$rd_id}; ## ALIGNED READ THAT HAS A DEF START AND END POSITION

    my $ct = 1;
    for(my $i = $rdstart; $i <= $rdend; $i++){
      my $ch = substr($curr_rd_str, ($i-1), 1);
      $corrected_read_start{$ct} = $ch;
      $ct++;
    }
    
    # Now we have to figure out how far along the alignment the start position was and where do the column, start and end position figure out along the alignment
    $ct = 1;

    foreach my $p(sort {$a <=> $b} keys %corrected_read_start){
      my $ch = $corrected_read_start{$p};

      if($ch ne '-'){
	#print "pos $p\t$ct\t$ch\t$qualvals{$ct}\n";
	$pos_char_score{$p}{$ch} = $qualvals{$ct};
	$ct++;
      }
      else{
	$pos_char_score{$p}{$ch} = -1;
      }
    }

    ## debug ##
 #    foreach my $k1(sort {$a <=> $b} keys %pos_char_score){
#       foreach my $k2(keys %{$pos_char_score{$k1}}){
# 	print "$k1\t$k2\t$pos_char_score{$k1}{$k2}\n";
#       }
#     }
    
    $ct = 1;
    my $ch = "";
    #print "id: $rd_id\t Found position: $res\n";
    if($pos > $rdstart){
      
      # the actual read sequence and corresponding quality scores start at this position
      $corr_start = ($beg - $rdstart + 1);
      $ct = 1;
      $ch = $corrected_read_start{$corr_start};
      $out{$ct}{$ch} = $pos_char_score{$corr_start}{$ch};
      #print "posi $pos Correct start: $corr_start char $ch $pos_char_score{$corr_start}{$ch}\n";
      
      $ct = 2;
      $corr_pos = ($pos - $rdstart + 1);
      $ch = $corrected_read_start{$corr_pos};
      $out{$ct}{$ch} = $pos_char_score{$corr_pos}{$ch};
      #print "posi $pos Correct pos: $corr_pos char $ch $pos_char_score{$corr_pos}{$ch}\n";
      
      $ct = 3;
      $corr_end = ($end -$rdstart + 1);
      $ch = $corrected_read_start{$corr_end};
      $out{$ct}{$ch} = $pos_char_score{$corr_end}{$ch};
      #print "posi $pos Correct end: $corr_end char $ch $pos_char_score{$corr_end}{$ch}\n";
    }
    elsif($pos == $rdstart){
      ## the start position is the position with the 
      $corr_pos = ($pos - $rdstart + 1);
      $corr_start = -1; # no previous position since the indel column is the start of the read
      
      $ct = 1;
      $ch = '-';
      $out{$ct}{$ch} = -1;
      
      $ct = 2;
      $ch = $corrected_read_start{$corr_pos};
      $out{$ct}{$ch} = $pos_char_score{$corr_pos}{$ch};
      
      $ct = 3;
      $corr_end = ($end - $rdstart +1);
      $ch = $corrected_read_start{$corr_end};
      $out{$ct}{$ch} = $pos_char_score{$corr_end}{$ch};
     #  print "posi $pos Correct start: $corr_start char $ch\n";
#       print "posi $pos Correct pos: $corr_pos char $ch $pos_char_score{$corr_pos}{$ch}\n";
#       print "posi $pos Correct end: $corr_end char $ch $pos_char_score{$corr_end}{$ch}\n";
    }
  }
  else{
    print "Read $rd_id with the sequence $curr_read was not found in the original fasta file $orig_read_file\n";
  }
  #print "$pos\t$qualvals{$pos}\t$prevpos\t$qualvals{$prevpos}\t$after_pos\t$qualvals{$after_pos}\n";
  return %out;
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
    my $e = $all_read_end{$rd};
    
    if( ($position >= $s) && ($position <= $e) ){
      ## Since the red names includes frequencies, include freq information
     #  my @vals = split(/\_/, $rd);
#       my $rdfreq = pop @vals;
      #print "$rd read freq: $rdfreq\n";
      $out = $out + 1;
    }
  }
  return $out;
}


#############
# Given an array, get the start and end position
#############
sub get_start_end{
    my($aref) = @_;

    my @arr = @$aref;

    # Both start and end positions are based based on 1-based indeces
    my $st = 0;
    my $end = 0;
    # start position
    for(my $i =0; $i < scalar @arr; $i++){
      $arr[$i] =~ s/\s//g; $arr[$i] =~ s/\t//g; $arr[$i] =~ s/\r//g; $arr[$i] =~ s/\n//g;
      if($arr[$i] ne '-'){
	$st = $i+1;
	$i = scalar @arr;
      }
    }
    # end position
    for(my $j = (scalar @arr -1); $j > 0; $j--){
      $arr[$j] =~ s/\s//g; $arr[$j] =~ s/\t//g; $arr[$j] =~ s/\r//g; $arr[$j] =~ s/\n//g;
      if($arr[$j] ne '-'){
	$end = $j+1;
	$j = 0;
      }
    }
    return($st, $end);
}
