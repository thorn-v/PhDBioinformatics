#!/usr/bin/perl
##Comparison.pl
## Modified from Yue Wang 2023
##IMPORTANT!!! This script now prints to stdout so you need to redirect
my %hash;
my $r='';

# my $file = "/Users/veronica/Documents/School/McMaster/OtherSoftware/mt_all_noRoot.min4.fasta";
my $file = $ARGV[0]; # takes bash argument
print STDERR "Processing file = $file";
open VCF, $file or die "Can't open file '$file' for read. $!"; #would be nice if not hardcoded
for $rec (<VCF>)
{
    chomp($rec);
    if ($rec =~ />(.*)$/) #gets the record ex. ">foobar"
    {
        $id = $1;
        # print $id,"\n";
    }
    else
    {
        # print $id . "\n";
        $hash{$id} = {
            record => $rec, 
            letters => split(undef, $rec) #splits the ssequence into a letter array for that id
        };
        # print $rec . "\n";
    }
}
close(VCF);

my @arr = keys(%hash); #gets all the ids

$len = length($hash{$arr[0]}->{record}); #gets length of the sequence, I think

# $file = "/Users/veronica/Documents/School/McMaster/OtherSoftware/distance_pairwise.txt";
# open O,">>", $file or die "Can't open file '$file' for write. $!";

print "\t"; #shifts over so the col names line up correctly
my $num_ids = scalar(@arr);

foreach my $id (@arr) #prints the col names
{
    print $id . "\t";
}

print "\n";
my @brr = (); #will dynamically build

for my $i (0..$#arr)
{
    print STDERR "."; #loading bar
    my $a_id = $arr[$i];
    my @a_letters = $hash{$a_id}->{letters};
    for my $j ($i..$#arr)
    {
        my $b_id = $arr[$j];
        my @b_letters = $hash{$b_id}->{letters};

        my $count = 0;
        for my $p (0..$len-1)
        {
            if ($a_letters[$p] ne $b_letters[$p])
            {
                $count++;
            }
        }
        
        push @brr, $count; #appends count to the end of brr 
    }

    my $str = join "\t", @brr;
    my @brrp=('x') x $i;
    my $s = join "\t", @brrp;
    @brr = ();
    if ($i == 0)
    {
        print $arr[$i],"\t",$str,"\n";
    }
    else
    {
        print $arr[$i],"\t",$s,"\t",$str,"\n";
    }
}
