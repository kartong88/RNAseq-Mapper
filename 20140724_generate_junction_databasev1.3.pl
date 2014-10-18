#!/usr/bin/perl

use strict;
use warnings;


# This script will generate the junction data base based
# based on the known transcript annotations


#my $ucsc_anno		= "/home/Kar\ Tong/annotations/Transcript_annotation/hg19_UCSCGenes_20140718.txt";
#my $refseq_anno		= "/home/Kar\ Tong/annotations/Transcript_annotation/hg19_refseq_20140718.txt";
#my $genecode_anno	= "/home/Kar\ Tong/annotations/Transcript_annotation/hg19_GENCODEgenesV19_20140718.txt";
#my $ensembl_anno	= "/home/Kar\ Tong/annotations/Transcript_annotation/hg19_EnsemblGenes_20140718.txt";


my $ucsc_anno		= "/home/kartong/Annotations/Transcript_annotation/hg19_UCSCGenes_20140718.txt";
my $refseq_anno		= "/home/kartong/Annotations/Transcript_annotation/hg19_refseq_20140718.txt";
my $genecode_anno	= "/home/kartong/Annotations/Transcript_annotation/hg19_GENCODEgenesV19_20140718.txt";
my $ensembl_anno	= "/home/kartong/Annotations/Transcript_annotation/hg19_EnsemblGenes_20140718.txt";


my $read_length			= $ARGV[0];
my $overlap_requirement	= 1;
my $required_length		= $read_length - $overlap_requirement;
my %anno_hash;
my $anno_hash_ref 		= \%anno_hash;
my $hg19_file			= "/home/kartong/genome/hg19/icgc_genome_broadvariant/Homo_sapiens_assembly19.fasta";

#my $hg19_file			= "/home/kartong/Downloads/temp.genome.fa";



open(my $error_junction, ">", "junctions_err.txt") || die $!;


parse_annotation($anno_hash_ref, $ucsc_anno, 8, 9, 1);
parse_annotation($anno_hash_ref, $refseq_anno, 9, 10, 2);
parse_annotation($anno_hash_ref, $genecode_anno, 9, 10, 2);
parse_annotation($anno_hash_ref, $ensembl_anno, 9, 10, 2);

#while( my ($k, $v) = each(%anno_hash)){
#	print $k . " " . $v . "\n";
#}


my $hg19_hash = parse_genome($hg19_file);


foreach my $junction_coords ( keys( %{ $anno_hash_ref }) ){
	my $seq = generate_junction_seq($junction_coords, $hg19_hash);
	
	# Check if there are
	if($seq eq "NULL"){
		print $error_junction $junction_coords . "\n";
	}
	else{
		print '>' .$junction_coords . "\n";
		print $seq . "\n";
	}
}


close($error_junction);


sub parse_annotation{
	my $result_hash		= $_[0];
	my $anno			= $_[1];
	my $exonstart_col	= $_[2];
	my $exonend_col		= $_[3];
	my $chr_col			= $_[4];
	

	open(my $anno_file, $anno) || die $!;

	my $header = <$anno_file>;

	while(my $line = <$anno_file>){
		$line =~ s/\r\n//;
		chomp($line);
		my @line_array 			= split(/\t/, $line);
		my @exonstart_array  	= split(/\,/, $line_array[$exonstart_col]); # Exon start coords is zero-based
		my @exonend_array  		= split(/\,/, $line_array[$exonend_col]); # Exon end coords is 1-based
		
		# Remove the last element because the string ends with
		# comma and then nothing
		#pop(@exonstart_array);
		#pop(@exonend_array);
		
		my $chr = $line_array[$chr_col];
		
		
		# Calculate the size of each exon
		my @exonsize_array;
		for(my $i=0; $i<scalar(@exonstart_array); $i++){
			my $exon_size = $exonend_array[$i] - $exonstart_array[$i];
			push(@exonsize_array, $exon_size);
		}
		
		# Loop through each junction and generate the coordiantes needed
		# Here we loop to scalar(@exonend_array) - 1 because there is one
		# less junction as comapred to the num of exon ends
		for(my $i=0; $i<scalar(@exonend_array) - 1; $i++){
			
			my @coords;
			my $current_posn = $i;
			
			# Identify the coords needed for LHS
			my $left_size; # total size of exons left of junction
			my $left_flag = 1;
			while($current_posn>=0 && $left_flag){
				my $end_posn	= $exonend_array[$current_posn];
				my $start_posn	= $exonstart_array[$current_posn] + 1;
				my $exon_size	= $exonsize_array[$current_posn];
				
				$left_size += $exon_size;
				
				unshift(@coords, $end_posn);
				unshift(@coords, $start_posn);
				
				# Check if using this exon is sufficient to cover 
				# the size that we need
				if($left_size > $required_length){
					# If including this exon is more than sufficient to cover 
					# the size that we need, then we now need to find the 
					# exact position that we need to stop
					
					my $excess_distance = $left_size - $required_length;
					my $start_posn_new = $start_posn + $excess_distance;
					
					# Since the exon size is already too big, we need to adjust 
					# the coords
					$coords[0] 	= $start_posn_new;
					$left_flag 	= 0;
				}
				$current_posn--;
			}
			
			
			$current_posn = $i + 1;
			
			# Identify the coords needed for RHS
			my $right_size; # total size of exons left of junction
			my $right_flag = 1;
			
			# print $current_posn . " curr\n";
			# print scalar(@exonstart_array) . "\n";
			# print join("\t", @exonstart_array) . "\n";
			while($current_posn< scalar(@exonstart_array) && $right_flag){
				my $start_posn	= $exonstart_array[$current_posn] + 1;
				my $end_posn	= $exonend_array[$current_posn];
				my $exon_size	= $exonsize_array[$current_posn];
				
				$right_size += $exon_size;
				
				push(@coords, $start_posn);
				push(@coords, $end_posn);
				
				# Check if using this exon is sufficient to cover 
				# the size that we need
				
				# print $right_size . "\n";
				# print $required_length . "\n";
				
				if($right_size > $required_length){
					# If including this exon is more than sufficient to cover 
					# the size that we need, then we now need to find the 
					# exact position that we need to stop
					
					my $excess_distance = $right_size - $required_length;
					my $end_posn_new = $end_posn - $excess_distance;
					
					# Since the exon size is already too big, we need to adjust 
					# the coords
					$coords[-1] 	= $end_posn_new;
					$right_flag 	= 0;
					# print "apple\n"
				}
				$current_posn++;
			}
			
			my $junction = $chr . "\:" . join("-", @coords);
			
			$result_hash ->{$junction} += 1;
		}
	}
	close($anno_file);
}




sub parse_genome{
# Parse the genome of interest and return it 
# as a hash.
	my $file = $_[0];
	
	open(my $genomefile, $file) || die $!;
	
	my %genome_hash;
	my $current_chr;	
	
	while(my $line = <$genomefile>){
		chomp($line);
		$line =~ s/\r\n//;
		
		# Check if it is header
		if ($line =~ m/^\>/ ){	
			my ($chr) = ($line =~ m/^\>(\S+)\s/);		
			$current_chr = $chr;
		}
		else{
			$genome_hash{$current_chr} .= $line
		}
	}	
	close($genomefile);

	return \%genome_hash;
}


sub generate_junction_seq{
	my $junction_coords = $_[0];
	my $genome_hash		= $_[1];
	

	my ($chrom, $coords) = split(/\:/, $junction_coords);
	
	
	$chrom =~ s/^chr//;
	
	my @coords_array = split(/\-/, $coords);
	my $seq;
	
	if( exists($genome_hash->{$chrom}) ){	
		for(my $i=0; $i<scalar(@coords_array); $i+=2){
			my $start 	= $coords_array[$i];
			my $end 	= $coords_array[$i+1];
			
			my $seqsize = $end - $start + 1;
			my $start_offsetted = $start -1;
			
			# Extract the seq needed
			my $subseq = substr($genome_hash->{$chrom}, $start_offsetted ,$seqsize);
			$seq .= $subseq;
		}
		

	}
	else{
		$seq = "NULL";
	}
	
	return $seq;
}




