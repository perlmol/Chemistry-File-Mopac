#!/home/ivan/bin/perl

use warnings;
use strict;
use Chemistry::File qw(MDLMol XYZ);
use blib;
use Chemistry::File::Mopac;

#my $mol = Chemistry::Mol->read($ARGV[0] || "test.mop", format => "mop");
#printf "%s\n", $mol->formula;
#print $_->coords->stringify("%10.4f%10.4f%10.3f\n") for $mol->atoms;
my $mol = Chemistry::Mol->read(shift);
#print $mol->print;
print $mol->print(format => "mop", coords => "internal", rebuild => 1);
