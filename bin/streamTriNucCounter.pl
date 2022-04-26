# Count trinucleotide occurances on the forward strand
# Soft-masked and N bases count result in XXX counts.
# Martin Taylor, based on TriNucClassifierCount.pl but modified June 2016 to
# accept a streamed output of tabluated sequence from bedtools getfasta

use strict;


my ($tripos) = classes();

my @outOrder = sort {$a cmp $b} keys %{$tripos};
my $name = "NULL";
($name) = @ARGV;

my $process=0;

my @header = ("SeqName",@outOrder);
print join "\t",@header;
print "\n";

my ($triposSeqLocal) = classes();
while (<STDIN>){
    chomp;
    my @sp = split /\t/;
    my $ln = length($sp[1]);
    next if ($ln < 3);
    
    my ($classed) = triNuc($triposSeqLocal,uc($sp[1]));
    #my @outvec = ($name);
    #foreach my $n (@outOrder){
#	$tripos->{$n}+=$triposSeqLocal->{$n};
 #   }
}
my @outvec = ($name);
foreach my $n (@outOrder){
    push @outvec, $triposSeqLocal->{$n};
}
print join "\t", @outvec;
print "\n";
    
sub classes {
    my %cls;
    my @nucs = qw{A C G T};
    my $id=1;
    foreach my $n (@nucs){
	foreach my $m (@nucs){
	    foreach my $o (@nucs){
		my $tri = $n.$m.$o;
		if (!defined $cls{$tri}){
		    $cls{$tri} = 0;
		}
	    }
	}
    }
    $cls{'XXX'} = 0;
    return(\%cls);
}

sub triNuc {
    my ($tref,$sqstream) = @_;
    my $ln = length($sqstream);
    for my $i (0 .. ($ln-3)){
	my $trip=substr($sqstream,$i,3);
	if (!defined $tref->{$trip}){
	    $tref->{'XXX'}++;
	}
	else {
	    $tref->{$trip}++;
	}
    }
    return(1);
}
