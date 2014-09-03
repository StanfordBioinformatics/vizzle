#!/usr/bin/perl
use strict;

my $in = shift;
my $genes = shift;
my $out = shift;
my $chr = shift;
#my $in = "/srv/gsfs0/projects/gbsc/benchmark/giab/NA12878.clia/hugeseq-1.3/hc2/manual-QC/coverage/acmg_coverage_hist.txt";
#my $genes = "/srv/gsfs0/projects/gbsc/Clinical_Service/dbases/ACMG/ACMG_genes.bed";
#my $out = "/srv/gsfs0/projects/gbsc/benchmark/giab/NA12878.clia/hugeseq-1.3/hc2/manual-QC/coverage/acmg_coverage_hist_genes.txt";

open (IN, "<", $in) or die "cannot open input sample file!";
open (OUT, ">", $out) or die "cannot open output file!";
#open (OUT, ">", "test.txt") or die "cannot open output file!";
open (GENE, "<", $genes) or die "cannot open input gene list!";

my %genes;

my $l = <GENE>;
while (my $l = <GENE>){
	chomp $l;
	my @rec = split("\t", $l);
	#print "@rec\n";
	push(@{$genes{$rec[0]}}, "$rec[1]\t$rec[2]\t$rec[3]");  
}
#print "@{$genes{'chr1'}}\n";

my $line=<IN>;
chomp $line;
#print OUT "$line\tgene\n";


#my $chr="chr1";
my @intervals;
foreach my $rec (@{$genes{$chr}}){
	my @rec = split("\t", $rec);
	my $beg = $rec[0];
	my $end = $rec[1];
	my $gene = $rec[2];
	push(@intervals, $beg);	
	push(@intervals, $end);	
	push(@intervals, $gene);	
}
#print "@intervals\n";
#die;

my @lines;
while (my $line=<IN>){
	if ($line=~m/$chr\:\d+/){
		chomp $line;
		while ($line=~m/$chr\:\d+/){
			push(@lines, $line);
			$line=<IN>;
			chomp $line;
		}
		goto END;
	}	
}
END:

my @olines;
my $search = 0;
my $gene_found = 0;
my $gene_beg = -1;
my $gene_end = -1;
my $gene_name;
#my @line = split("\t", $lines[0]);
#my $loc = $1 if ($line[0] =~ m/chr[\dXYM]{1,2}:(\d+)/); 
#$loc = $loc-1;
foreach my $line (@lines){
	#$loc = $loc + 1;
	my @l = split("\t", $line);
	my $loc = $1 if ($l[0] =~ m/chr[\dXYM]{1,2}:(\d+)/); 
	#print ">$loc\t$gene_end\n";
	#die;
	if ($gene_found>0){
		if ($loc <= $gene_end){
			#print OUT ">>>$loc\t$search\t$gene_name\n";
			push(@olines,  "$chr\t$loc\t$l[1]\t$l[2]\t$l[3]\t$gene_name\n");
		} else {
			foreach my $l (@olines){
				print OUT "$l";
			}
			#print OUT "$loc\n";
			$gene_found=0;
			$search=0;
			$gene_beg = -1;
			$gene_end = -1;
			$gene_name = "";
			@olines=();
			#die;
		}
	} else {
		foreach my $rec (@intervals){
			#print "$loc\t$rec\t$search\n";
			#if ($loc >= $rec && $search==1 && $rec !~ m/[A-Z]/i){
			#	$gene_beg = $rec;
			#	print "1:$gene_name\t$gene_beg\t$gene_end\n";
			#	next;
			#} elsif ($loc <= $rec && $search==1 && $rec !~ m/[A-Z]/i){
			if ($loc <= $rec && $search==1 && $rec !~ m/[A-Z]/i){
				$gene_end = $rec;
				#print "2:$gene_name\t$gene_beg\t$gene_end\n";
				next;	
			} elsif ($loc > $rec && $search==1 && $rec !~ m/[A-Z]/i){
				#print "3:$gene_name\t$gene_beg\t$gene_end\n";
				$search=0;
				next;	
			} elsif ($search == 1 && $rec =~ m/[A-Z]/i){ 
				push(@olines,  "$chr\t$loc\t$l[1]\t$l[2]\t$l[3]\t$rec\n");
				$gene_found = 1;
				$gene_name = $rec;
				print "4:$gene_name\t$gene_beg\t$gene_end\n";
				#die;
				goto GENEFOUND;
			} elsif ($loc >= $rec && $search==0 && $rec !~ m/[A-Z]/i){
				$gene_beg = $rec;
				$search = 1;
				#print "5:$rec\t$gene_name\t$gene_beg\t$gene_end\n";
			}
		}
	}
	GENEFOUND:
	#print "@olines\n";
}
#print OUT "$chr\t$loc\t$line[1]\t$line[2]\t$line[3]\t$gene\n";
close IN;
close OUT;
