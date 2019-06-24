#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


########################################################
### Concatenation of treated peaks and control peaks ###
########################################################

my $inputA = shift;
my $inputB =shift;
my $output = shift;


my %temoin;

open (F1,"<",$inputA) or die "ouverture $inputA impossible, cause : $!";
while (my $li = <F1>)
{
	chomp $li;
	my @T = split("\t",$li);
	my @final;
	for (my $i = 0 ; $i < 3 ; $i++)
	{
		push(@final, $T[$i]);
	}
	my $line = join ("\t",@final);
	$temoin{$T[3]} = $line;
}
close F1;


my %treatment;

open (F1,"<",$inputB) or die "ouverture $inputB impossible, cause : $!";
while (my $li = <F1>)
{
	chomp $li;
	my @A = split("\t",$li);
	my @final;
	for (my $i = 0 ; $i < 3 ; $i++)
	{
		push(@final, $A[$i]);
	}
	my $line = join ("\t",@final);
	$treatment{$A[3]} = $line;
}
close F1;

# creation of a new file which contain all peaks

open (FE1,">",$output) or die "création impossible, cause $!";
print FE1 "#./$inputA\n";
my @nA = split("\.",$inputA);
my $nameA = "$nA[0]";
foreach my $k (keys (%temoin))
{
	print FE1 "$temoin{$k}\t$k\t$nameA\n";
}

print FE1 "#./$inputB\n";
my @nB = split("\.",$inputB);
my $nameB = "$nB[0]";
foreach my $k (keys (%treatment))
{ 
	print FE1 "$treatment{$k}\t$k\t$nameB\n";
}
close FE1;