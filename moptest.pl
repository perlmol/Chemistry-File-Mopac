#!/home/ivan/bin/perl

use warnings;
use strict;
use Chemistry::Mol;
use blib;
use Chemistry::File::Mopac;

my $mol = Chemistry::Mol->read($ARGV[0] || "test.mop", format => "mop");

printf "%s\n", $mol->formula;
print $_->coords for $mol->atoms;
