package Chemistry::File::Mopac;

$VERSION = '0.10';

use 5.006001;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol 0.10;
use Carp;

=head1 NAME

Chemistry::File::Mopac

=head1 SYNOPSIS

    use Chemistry::File::Mopac;

    my $mol = Chemistry::Mol->read('file.mop');

=cut

=head1 DESCRIPTION

This module reads Mopac 6 input files. It can handle both internal coordinates
and cartesian coordinates.

This module registers the C<mop> format with Chemistry::Mol. For detection
purposes, it assumes that filenames ending in .mop or .zt have the Mopac 
format.

=cut

Chemistry::Mol->register_format("mop");

my ($C_sym, $C_1, $C_o1, $C_2, $C_o2, $C_3, $C_o3, $C_len, 
    $C_ang, $C_dih) = 0 .. 9;
my %Pos = (
    sym     => 0,
    x       => 1,
    x_opt   => 2,
    y       => 3,
    y_opt   => 4,
    z       => 5,
    z_opt   => 6,
    l       => 1,
    l_opt   => 2,
    a       => 3,
    a_opt   => 4,
    d       => 5,
    d_opt   => 6,
    l_ref   => 7,
    a_ref   => 8,
    d_ref   => 9,
);

sub parse_string {
    my $class = shift;
    my $string = shift;
    my %opts = @_; 
    my $mol_class = $opts{mol_class} || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || "Chemistry::Atom";
    my $bond_class = $opts{bond_class} || "Chemistry::Bond";

    my $mol = $mol_class->new();
    my @lines = split "\n", $string;

    # read header
    my @header = splice @lines, 0, 3;
    my @keys;

    $keys[0] = $header[0];
    if ($keys[0] =~ /&/) {
        $keys[1] = $header[1];
        if ($keys[1] =~ /&/) {
            $keys[2] = $header[2];
        }
    } elsif ($keys[0] =~ /\+/) {
        $keys[1] = $header[1];
        push @header, shift @lines;
        if ($keys[1] =~ /\+/) {
            $keys[2] = $header[2];
            push @header, shift @lines;
        }
    }

    $mol->attr(keys_line => join "\n", @keys);
    $mol->attr(text_line => join "\n", @header[@keys..$#header]);
    
    my @coords;
    my $mode;
    # read coords
    for (@lines) {
      # Sample line below
      # O    1.232010  1  128.812332  1  274.372818  1    3   2   1
        last if /^\s*$/; #blank line
        push @coords, [split];
    }
    
    if (@coords <= 3 or $coords[3][$Pos{l_ref}]) { # Internal coords
        $mode = 'internal';
    } else { # Cartesian coords
        $mode = 'cartesian';
    }

    my @int_coords;
    for my $coord (@coords) {
        my $atom = $mol->new_atom(symbol => $coord->[$Pos{sym}]);
        if ($mode eq 'internal') {
            $atom->attr("int/len_val" => $coord->[$Pos{l}]);
            $atom->attr("int/ang_val" => $coord->[$Pos{a}]);
            $atom->attr("int/dih_val" => $coord->[$Pos{d}]);
            $atom->attr("mop/len_opt" => $coord->[$Pos{l_opt}]);
            $atom->attr("mop/ang_opt" => $coord->[$Pos{a_opt}]);
            $atom->attr("mop/dih_opt" => $coord->[$Pos{d_opt}]);
            $atom->attr("int/len_ref" => $coord->[$Pos{l_ref}]);
            $atom->attr("int/ang_ref" => $coord->[$Pos{a_ref}]);
            $atom->attr("int/dih_ref" => $coord->[$Pos{d_ref}]);
            push @int_coords, [@$coord[$C_1, $C_len, $C_2, $C_ang, 
                                        $C_3, $C_dih]]
        } else { # cartesian
            croak "Cartesian not yet implemented";
        }
    }

    my $i = 1;
    for my $v (int_to_cart(@int_coords)) {
        $mol->atoms($i++)->coords($v);
    }
    return $mol;
}

sub file_is {
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
    @{$ints[0]}[2..5] = (90, -1, -90, 0) if @ints; # refer 1st atom to Y,Z
    @{$ints[1]}[4..5] = (-90, 0) if @ints > 1; # make second atom refer to Z

    for my $int (@ints) {
        my $v1 = $carts[$base+$int->[5]]; # 'oldest' point
        my $v2 = $carts[$base+$int->[3]];
        my $v3 = $carts[$base+$int->[1]]; # 'newest' point
        my $d1 = $v1 - $v2;
        my $d2 = $v3 - $v2;

        # $xp = normal to atoms 1 2 3 
        my $xp = $d1 x $d2;
        # $yp = normal to xp and atoms 2 3
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


=head1 SEE ALSO

L<Chemistry::Mol>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=cut

