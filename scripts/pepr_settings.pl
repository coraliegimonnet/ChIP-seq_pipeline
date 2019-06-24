###########################################################
### Script to create a settings file for PePr execution ###
###########################################################

use strict;
use warnings;
use Data::Dumper;

my $input = shift;
my $cpus = shift;
my $windowSize = shift;


my %h;

# read file create by echo bash command
open (F1,"<",$input) or die "cannot open $input, cause : $!\n";
while (my $li = <F1>)
{
	chomp $li;
	my @T = split (",",$li);
	# create a hash - key = experiment - value = line
	if (!defined $h{$T[3]})
	{
		$h{$T[3]} = [];
	}
	push (@{$h{$T[3]}}, $li);
}
close F1;

# write settings files for pepr execution
foreach my $k (keys (%h))
{
	my $output = "$k\.setting_PePr.txt";
	open (FE1, ">", $output) or die "cannot create $output, cause : $!\n";
	print FE1 "#filetype\tfilename\n";
	for (my $i = 0 ; $i < scalar(@{$h{$k}}) ; $i++)
	{
		my @T = split (',',${$h{$k}}[$i]);
		my @C = split (/\./,$T[0]);
		my @I = split (/\./, $T[1]);
		print FE1 "chip1\t$C[0].dedup.bam\n";
		print FE1 "input1\t$I[0].dedup.bam\n";
	}
	print FE1 "file-format\tbam\n";
	print FE1 "peaktype\tsharp\n";
	print FE1 "difftest\tFALSE\n";
	print FE1 "name\t$k\n";
	if ($windowSize ne "false")
	{
		print FE1 "windowsize\t$windowSize\n";
	}
	if (defined $cpus)
	{
		print FE1 "num-processors\t$cpus\n";
	}
	close FE1;
}
