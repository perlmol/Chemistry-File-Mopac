package Chemistry::File::Mopac;

$VERSION = '0.01';

use Chemistry::Mol 0.07;
use Carp;
use 5.006001;
use strict;
use warnings;
use base "Chemistry::File";
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

=head1 NAME

Chemistry::File::Mopac

=head1 SYNOPSIS

    use Chemistry::File::Mopac;

=cut

=head1 DESCRIPTION

This module reads Mopac input files. 

=cut

Chemistry::Mol->register_type(
    "mop" => __PACKAGE__,  
);

=head1 FUNCTIONS

=over 4

=item $my_mol = mop_read($fname, option => value...)

=cut

my ($C_sym, $C_1, $C_o1, $C_2, $C_o2, $C_3, $C_o3, $C_len, 
    $C_ang, $C_dih) = 0 .. 9;

sub parsefile {
    my $class = shift;
    my $fname = shift;
    my %options = @_; 
    local $_;

    open F, $fname or croak "Could not open file $fname";

    my $mol = $options{mol} || Chemistry::Mol->new();

    # read header
    my @header;
    for (1..3) {
        my $line = <F>;
        chomp $line;
        push @header, $line;
    }
    my @keys;

    $keys[0] = $header[0];
    if ($keys[0] =~ /&/) {
        $keys[1] = $header[1];
        if ($keys[1] =~ /&/) {
            $keys[2] = $header[2];
        }
    } elsif ($keys[0] =~ /\+/) {
        $keys[1] = $header[1];
        push @header, scalar <F>;
        if ($keys[1] =~ /\+/) {
            $keys[2] = $header[2];
            push @header, scalar <F>;
        }
    }

    $mol->attr(keys_line => join "\n", @keys);
    $mol->attr(text_line => join "\n", @header[@keys..$#header]);
    
    my @coords;
    my $mode;
    # read coords
    while (<F>) {
      # O    1.232010  1  128.812332  1  274.372818  1    3   2   1
        last if /^\s*$/; #blank line
        push @coords, [split];
    }
    
    if (@coords <= 3 or $coords[3][$C_len]) { # Internal coords
        $mode = 'internal';
    } else { # Cartesian coords
        $mode = 'cartesian';
    }

    my @int_coords;
    for my $coord (@coords) {
        my $atom = $mol->new_atom(symbol => $coord->[$C_sym]);
        if ($mode eq 'internal') {
            $atom->attr("int/len_val" => $coord->[$C_1]);
            $atom->attr("int/ang_val" => $coord->[$C_2]);
            $atom->attr("int/dih_val" => $coord->[$C_3]);
            $atom->attr("mop/len_opt" => $coord->[$C_o1]);
            $atom->attr("mop/ang_opt" => $coord->[$C_o2]);
            $atom->attr("mop/dih_opt" => $coord->[$C_o3]);
            $atom->attr("int/len_ref" => $coord->[$C_len]);
            $atom->attr("int/ang_ref" => $coord->[$C_ang]);
            $atom->attr("int/dih_ref" => $coord->[$C_dih]);
            push @int_coords, [@{$coord}[$C_1, $C_len, $C_2, $C_ang, 
                                        $C_3, $C_dih]]
        }
    }
    close F;

    my $i = 1;
    for my $v (int_to_cart(@int_coords)) {
        $mol->atoms($i++)->coords($v);
    }
    return $mol;
}

=item is_mop($fname)

Returns true if the specified file is a Mopac file. 

=cut

sub isfile {
    my $class = shift;
    my $fname = shift;
    
    return 1 if $fname =~ /\.(?:mop|zt)$/i;

    open F, $fname or croak "Could not open file $fname";
    
    my $line = <F>;
    close F;
    return 1 if $line =~ /am1|pm3|mndo|mdg|pdg/i;
    return 0;
}


# Expects a list of array references, where each array lists
# [ distance, atom#, angle(deg), atom#, dihedral(deg), atom# ]
# returns a list of vectors with the cartesian coordinates.
sub int_to_cart {
    #use Data::Dumper;
    #local $Data::Dumper::Indent = 0;
    use Math::VectorReal ":all";
    #local $Math::VectorReal::FORMAT = "[ %.5f %.5f %.5f ]\n";
    my (@ints) = @_; #internal coordinates
    my @carts = (X, Y, Z); # base cartesians
    my $base = 2; # point 3 is the first 'real' one
    my $pi180 = 3.1415926535/180;
    
    return () unless @ints; # make sure we have something
    shift @ints; # throw away first point
    push @carts, vector(0,0,0); # first point at origin
    @{$ints[0]}[2..5] = (90, -1, -90, 0) if @ints; # refer 1st atom to YZ
    @{$ints[1]}[4..5] = (-90, 0) if @ints > 1; # make second atom refer to Z

    for my $int (@ints) {
        my $v1 = $carts[$base+$int->[5]]; # 'oldest' point
        my $v2 = $carts[$base+$int->[3]];
        my $v3 = $carts[$base+$int->[1]]; # 'newest' point
        my $d1 = $v1 - $v2;
        my $d2 = $v3 - $v2;

        # $xp = normal to atoms 1 2 3 
        my $xp = $d1 x $d2;
        my $yp = $d2 x $xp;

        my $ang1 = $int->[4] * $pi180;   # dihedral
        # $r = normal to atoms 2 3 4 (where 4 is the new atom)
        #      obtained by rotating $xp through $d2
        my $r = $xp->norm * cos($ang1) + $yp->norm * sin($ang1);

        my $ypp = $d2 x $r; # complete new frame of reference
        my $ang2 = $int->[2] * $pi180;   # angle
        my $d3 = -$d2->norm * cos($ang2) + $ypp->norm * sin($ang2);

        $d3 = $d3 * $int->[0]; # mult by distance to $v3
        my $v4 = $v3 + $d3; # define new point
        push @carts, $v4;
        
        #print "INT: ", Dumper ($int), "\n";
        #print "v1:$v1 v2:$v2 v3:$v3 d1:$d1 d2:$d2 xp:$xp yp:$yp r:$r ypp:$ypp d3:$d3 v4:$v4\n";
        #print "v4: $v4";
    }
    @{$ints[0]}[2..5] = (0,0,0,0) if @ints; # clean up
    @{$ints[1]}[4..5] = (0,0) if @ints > 1;         
    splice @carts, 0, 3; # throw away XYZ

    return @carts;
}

1;


=back

=head1 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>

The PDB format description at 
L<http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=cut

