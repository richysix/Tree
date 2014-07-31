## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package Tree::GenomicIntervalTree;
## use critic

# ABSTRACT: GenomicIntervalTree for efficiently retrieving genomic intervals

use namespace::autoclean;
use Moose;
use autodie;
use Carp qw(cluck confess);
use Readonly;

use Set::IntervalTree;

=method genomic_tree

  Usage       : $tree->genomic_tree;
  Purpose     : Getter for genomic_tree attribute
  Returns     : HashRef
  Parameters  : None
  Throws      : If input is given
  Comments    : The keys of the returned hashref are seq_region_names
                (chromosome/scaffolds)
                The values are Set::IntervalTree objects

=cut

has 'genomic_tree' => (
    is      => 'ro',
    isa     => 'HashRef',
    default => sub { {} },
);

=method add_intervals_to_genomic_tree_from_hash

  Usage       : $tree->add_intervals_to_genomic_tree_from_hash( $hashref );
  Purpose     : add intervals to a genomic interval tree
  Returns     : 1 on Success
  Parameters  : Hash or HashRef
  Throws      : If Hash(Ref) does not match expected data structure
  Comments    : TO DO: Accepting a hash as well as a hashref

=cut

sub add_intervals_to_genomic_tree_from_hash {
    my ( $self, $hash ) = @_;

    # check $hash
    # $hash is expected to have structure:
    Readonly my $HASH_STRUCTURE => <<'END_HASH_STRUCTURE';
$hash = {
    CHR => {
        'START1-END1' => $object,
        'START2-END2' => $object,
    }
}
END_HASH_STRUCTURE
    if ( ref $hash ne 'HASH' ) {
        confess "Supplied object does not match required structure!\n\n",
          $HASH_STRUCTURE;
    }
    foreach my $chr ( keys %{$hash} ) {
        if ( ref $hash->{$chr} ne 'HASH' ) {
            confess "Supplied object does not match required structure!\n\n",
              $HASH_STRUCTURE;
        }
        foreach my $interval ( keys %{ $hash->{$chr} } ) {
            if ( $interval !~ m/\d+\-\d+/xms ) {
                confess
                  "Supplied object does not match required structure!\n\n",
                  $HASH_STRUCTURE;
            }
        }
    }

    foreach my $chr ( keys %{$hash} ) {
        foreach my $interval ( keys %{ $hash->{$chr} } ) {
            my ( $start, $end ) = split /-/xms, $interval;
            my $object = $hash->{$chr}->{$interval};
            $self->insert_interval_into_tree( $chr, $start, $end, $object, );
        }
    }
    return 1;
}

=method fetch_overlapping_intervals

  Usage       : $tree->fetch_overlapping_intervals( $chr, $start, $end );
  Purpose     : fetch intervals from the tree that match the supplied interval
  Returns     : ArrayRef
  Parameters  : CHR:    String
                START:  Integer
                END:    Integer
  Throws      : 
  Comments    : TO DO: Checking input is good

=cut

sub fetch_overlapping_intervals {
    my ( $self, $chr, $q_start, $q_end ) = @_;

    if ( exists $self->genomic_tree->{$chr} ) {
        my $results =
          $self->genomic_tree->{$chr}->fetch( $q_start, $q_end + 1 );
        return $results;
    }
    else {
        return [];
    }
}

=method insert_interval_into_tree

  Usage       : $tree->insert_interval_into_tree( $chr, $start, $end, $object );
  Purpose     : inserts an interval to a genomic interval tree
  Returns     : 1 on Success
  Parameters  : CHR:    String
                START:  Integer
                END:    Integer
                OBJECT: Any
  Throws      : 
  Comments    : 

=cut

sub insert_interval_into_tree {
    my ( $self, $chr, $start, $end, $object ) = @_;
    if ( !exists $self->genomic_tree->{$chr} ) {
        $self->_make_new_subtree_for_chr($chr);
    }
    $self->genomic_tree->{$chr}->insert( $object, $start, $end + 1 );

    return 1;
}

=method _make_new_subtree_for_chr

  Usage       : $tree->_make_new_subtree_for_chr( $chr );
  Purpose     : creates an new empty Set::InetrvalTree for the supplied
                seq_region
  Returns     : Set::InetrvalTree object
  Parameters  : CHR:    String
  Throws      : 
  Comments    : 

=cut

sub _make_new_subtree_for_chr {
    my ( $self, $chr ) = @_;

    # check tree doesn't already exist
    if ( exists $self->genomic_tree->{$chr} ) {
        warn "Tree already exists!\n";
        return $self->genomic_tree->{$chr};
    }
    else {
        my $sub_tree = Set::IntervalTree->new;
        $self->genomic_tree->{$chr} = $sub_tree;
        return $sub_tree;
    }
}

__PACKAGE__->meta->make_immutable;
1;

__END__

=pod

=head1 SYNOPSIS

    use Tree::GenomicIntervalTree;
    my $genomic_tree = Tree::GenomicIntervalTree->new()

    my $intervals = {
        1 => {
            '1-10' => 'exon1.1',
            '11-29' => 'intron1.1',
            '30-40' => 'exon1.2',
        },
        2 => {
            '10-20' => 'exon2.1',
            '21-25' => 'intron2.1',
            '26-45' => 'exon2.2',
            '21-35' => 'intron2.2',
            '36-45' => 'exon2.3',
        },
        'Zv9_scaffold1345' => {
            '10-20' => 'exon3.1',
            '21-29' => 'intron3.1',
            '30-50' => 'exon3.2',
        },
    };

    # add intervals to tree
    $genomic_tree->add_intervals_to_genomic_tree_from_hash( $test_intervals );

    # add a single interval
    $genomic_tree->insert_interval_into_tree( $chr, $start, $end, $object );

    # fetch overlapping intervals
    my $results_arrayref = $genomic_tree->fetch_overlapping_intervals( '1', 8,
        35 );


=head1 DESCRIPTION

    Tree::GenomicIntervalTree is an implementation of an interval tree for
    genomic intervals. Basically, it contains a separate interval tree for each
    chromosome/scaffold. Intervals used in this implementation are closed. i.e.
    The interval 10-20 contains both the endpoints (10 and 20). [10,20]

=cut
