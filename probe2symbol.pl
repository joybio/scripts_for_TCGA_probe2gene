use strict;
use warnings;
print STDERR "gene symbol column number: ";
my $geneSymbolCol=<STDIN>;
chomp($geneSymbolCol);
$geneSymbolCol--;
my $expFile="probeMatrix.txt";
my $gplFile="ann.txt";
my $expFileWF="geneMatrix.txt";
my %hash=();
my @sampleName=();

open(EXP,"$expFile") or die $!;
while(my $exp=<EXP>)
{
	next if ($exp=~/^(\n|\!)/);
	chomp($exp);
	if($.==1)
	{
		my @expArr=split(/\t/,$exp);
		for(my $i=0;$i<=$#expArr;$i++)
		{
			my $singleName=$expArr[$i];
			$singleName=~s/\"//g;
			if($i==0)
			{
				push(@sampleName,"ID_REF");
			}
			else
			{
				my @singleArr=split(/\_|\./,$singleName);
				push(@sampleName,$singleArr[0]);
			}
		}
	}
	else
	{
		my @expArr=split(/\t/,$exp);
		for(my $i=0;$i<=$#sampleName;$i++)
		{
			$expArr[$i]=~s/\"//g;
			push(@{$hash{$sampleName[$i]}},$expArr[$i]);
		}
	}
}
close(EXP);

my %probeGeneHash=();

open(GPL,"$gplFile") or die $!;
while(my $gpl=<GPL>)
{
	next if($gpl=~/^(\#|ID|\!|\n)/);
	chomp($gpl);
	my @gplArr=split(/\t/,$gpl);
	if((exists $gplArr[$geneSymbolCol]) && ($gplArr[$geneSymbolCol] ne '') && ($gplArr[$geneSymbolCol] !~ /.+\s+.+/))
	{
		$gplArr[$geneSymbolCol]=~s/(.+?)\/\/\/(.+)/$1/g;
		$gplArr[$geneSymbolCol]=~s/\"//g;
		$probeGeneHash{$gplArr[0]}=$gplArr[$geneSymbolCol];
	}
}
close(GPL);

my @probeName=@{$hash{"ID_REF"}};
delete($hash{"ID_REF"});

my %geneListHash=();
my %sampleGeneExpHash=();
foreach my $key (keys %hash)
{
	my %geneAveHash=();
	my %geneCountHash=();
	my %geneSumHash=();
	my @valueArr=@{$hash{$key}};
	for(my $i=0;$i<=$#probeName;$i++)
	{
		if(exists $probeGeneHash{$probeName[$i]})
		{
			my $geneName=$probeGeneHash{$probeName[$i]};
			$geneListHash{$geneName}++;
			$geneCountHash{$geneName}++;
			$geneSumHash{$geneName}+=$valueArr[$i];
		}
	}
	foreach my $countKey (keys %geneCountHash)
	{
		$geneAveHash{$countKey}=$geneSumHash{$countKey}/$geneCountHash{$countKey};
	}
	$sampleGeneExpHash{$key}=\%geneAveHash;
}

open(WF,">$expFileWF") or die $!;
$sampleName[0]="geneNames";
print WF join("\t",@sampleName) . "\n";
foreach my $probeGeneValue (sort(keys %geneListHash))
{
	print WF $probeGeneValue . "\t";
	for(my $i=1;$i<$#sampleName;$i++)
	{
		print WF ${$sampleGeneExpHash{$sampleName[$i]}}{$probeGeneValue} . "\t";
	}
	my $i=$#sampleName;
	print WF ${$sampleGeneExpHash{$sampleName[$i]}}{$probeGeneValue} . "\n";
}
close(WF);
