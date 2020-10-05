#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my %opts = ();
GetOptions(
	"threshold=s"		=> \$opts{"threshold"},
	"rename_tbl=s"		=> \$opts{"rename_table"},
	"dir=s"				=> \$opts{"dir"},
	"fasta=s"			=> \@{$opts{"fasta"}}, 
	"stats=s"			=> \$opts{"stats"},
	"help" 				=> \$opts{help},
);

$opts{fasta_count}=scalar @{$opts{"fasta"}};
validate_options(\%opts);


my %stats=(
	"number of fasta files"=>0,
	"pass completeness"=>0,
	"fail completeness"=>0,
	"completeness threshold"=>$opts{threshold}
);

if($opts{stats}){
	(open $opts{STATS},">",$opts{stats}) || usage("unable to open table for completeness stats at $opts{stats}");
	print {$opts{STATS}} "completeness\theader\nfile\n";
}

my %rename_tbl;
if($opts{rename_table}){
	(open my $TBL,"<",$opts{rename_table}) || usage("unable to open rename table $opts{rename_table}");
	while(my $rec=<$TBL>){
		chomp $rec;
		my ($oldname,$newname)=split /\t/,$rec;
		$rename_tbl{$oldname}=$newname;
	}
}


for my $fa(@{$opts{fasta}}){
	
	(open my $FA,"<",$fa) || usage("unable to open fasta file $fa");
	
	$stats{"number of fasta files"}++;
	
	my $headerline=<$FA>;chomp $headerline;
	### cleanup the headerline
	$headerline=~s/>Consensus_/>/;
	$headerline=~s/\.primertrimmed.*//;
	
	### rename the consensus, if there is a rename table
	if($opts{rename_table}){
		my ($name)=$headerline=~/>(.*)/;
		if(my $newname=$rename_tbl{$name}){
			$headerline=">$newname";
		}
	}
	
	
	my @seqlines=<$FA>;chomp @seqlines;
	my $seq=join("",@seqlines);
	
	my $seqlength=length($seq);
	my $ncount = () = $seq=~/(n)/ig;
	my $completeness=($seqlength-$ncount)/$seqlength;
	
	if($opts{STATS}){
		print {$opts{STATS}} "$completeness\t$headerline\t$fa\n";
	}
	
	
	my $rec = "$headerline\n$seq";
	if($completeness>=$opts{threshold}){
		print STDOUT "$rec\n";
		$stats{"pass completeness"}++;
		
	}else{
		print STDERR "failed completeness, $completeness $fa\n";
		$stats{"fail completeness"}++;
	}
}


print STDERR "Completeness Counts\n";
print STDERR join("\n",map{ "$_ : $stats{$_} "} ("completeness threshold","number of fasta files","pass completeness","fail completeness")) ."\n";

close $opts{STATS} if($opts{STATS});


sub validate_options{
	my ($opts)=@_;
	

	usage("Help requested.") if($opts{help});
	usage("completness threshold is required") unless($opts{threshold});
	usage("completness threshold must be between 0 and 1") unless( $opts{threshold}>=0 && $opts{threshold}<=1  );
	
	

	
	if($opts{fasta_count} && $opts{dir}){
		usage("only one of --fasta or --dir can be used");
	}
	if( !$opts{fasta_count} && !$opts{dir} ){
		usage("either --fasta or --dir must be provided");
	}
	
	### get the list of fasta files to process, store in opts
	if($opts{dir}){
		opendir(my $DIR,$opts{dir}) || usage("directory $opts{dir} not found");
		my @fa=grep{ /\.fa$/ } readdir $DIR;
		@{$opts{fasta}}=map{"$opts{dir}/$_"} @fa;
		
		closedir $DIR;
		$opts{fasta_count}=scalar @{$opts{"fasta"}};
		if(! $opts{fasta_count}){
			usage("no fasta files found in $opts{dir}");
		}
		if($opts{fasta_count}<2){
			usage("less than 2 fasta files found in $opts{dir}");
		}
	}
}

sub usage{

	print "\npreprocess_fasta.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--threshold.  Completeness threshold, required\n";
	print "\t--rename. Table to rename consensus sequence names\n";
	print "\tinput data, only one of --dir or --fasta can be used.  Fasta files must have .fa suffix\n";
	print "\t  --dir. Directory of consensus fasta files\n";
	print "\t  --fasta. Consensus fasta files, can take multiple values\n";
	print "\t--stats.  Path to a file to completeness statistics\n";
	print "\t--help displays this usage message.\n";
	die "\n@_\n\n";
}
