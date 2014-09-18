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
  Returns     : HashRef of HashRefs
  Parameters  : None
  Throws      : If input is given
  Comments    : The keys of the returned hashref are seq_region_names
                (chromosome/scaffolds)
                The values are hashrefs keyed by strand (-1, 0 or 1) of
                Set::IntervalTree objects

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
        'START1-END1'         => $object,
        'START2-END2:STRAND2' => $object,
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
            if ( $interval !~ m/ \A \d+ [-] \d+ :? [\d+-]* \z /xms ) {
                confess
                  "Supplied object does not match required structure!\n\n",
                  $HASH_STRUCTURE;
            }
        }
    }

    foreach my $chr ( keys %{$hash} ) {
        foreach my $interval ( keys %{ $hash->{$chr} } ) {
            my ( $start, $end, $strand ) =
              $interval =~ m/ \A (\d+) [-] (\d+) :? ([\d+-]*) \z /xms;
            my $object = $hash->{$chr}->{$interval};

            $strand = $self->_normalise_strand($strand);

            $self->insert_interval_into_tree( $chr, $start, $end, $strand,
                $object );
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
                STRAND: Integer
  Throws      : 
  Comments    : TO DO: Checking input is good

=cut

sub fetch_overlapping_intervals {
    my ( $self, $chr, $q_start, $q_end, $q_strand ) = @_;

    $q_strand = $self->_normalise_strand($q_strand);

    ## no critic (ProhibitMagicNumbers)
    my @q_strands = ( -1, 0, 1 );
    if ( defined $q_strand && ( $q_strand == 1 || $q_strand == -1 ) ) {
        ## use critic
        @q_strands = ( $q_strand, 0 );
    }

    my $results = [];

    foreach my $strand (@q_strands) {
        if ( exists $self->genomic_tree->{$chr}->{$strand} ) {
            push @{$results},
              @{ $self->genomic_tree->{$chr}->{$strand}
                  ->fetch( $q_start, $q_end + 1 ) };
        }
    }

    return $results;
}

=method insert_interval_into_tree

  Usage       : $tree->insert_interval_into_tree( $chr, $start, $end, $object );
  Purpose     : inserts an interval to a genomic interval tree
  Returns     : 1 on Success
  Parameters  : CHR:    String
                START:  Integer
                END:    Integer
                STRAND: Integer
                OBJECT: Any
  Throws      : 
  Comments    : 

=cut

sub insert_interval_into_tree {    ## no critic (ProhibitManyArgs)
    my ( $self, $chr, $start, $end, $strand, $object ) = @_;

    # Support old interface (i.e. no strand parameter)
    if ( !defined $object ) {
        $object = $strand;
        $strand = 0;
    }

    $strand = $self->_normalise_strand($strand);

    if ( !exists $self->genomic_tree->{$chr}->{$strand} ) {
        $self->_make_new_subtree( $chr, $strand );
    }

    $self->genomic_tree->{$chr}->{$strand}->insert( $object, $start, $end + 1 );

    return 1;
}

=method _make_new_subtree

  Usage       : $tree->_make_new_subtree( $chr, $strand );
  Purpose     : creates an new empty Set::IntervalTree for the supplied
                seq_region and strand
  Returns     : Set::InetrvalTree object
  Parameters  : CHR:    String
                STRAND: Integer
  Throws      : 
  Comments    : 

=cut

sub _make_new_subtree {
    my ( $self, $chr, $strand ) = @_;

    # 0 = unstranded interval
    if ( !defined $strand ) {
        $strand = 0;
    }

    # check tree doesn't already exist
    if ( exists $self->genomic_tree->{$chr}->{$strand} ) {
        warn "Tree already exists!\n";
        return $self->genomic_tree->{$chr}->{$strand};
    }
    else {
        my $sub_tree = Set::IntervalTree->new;
        $self->genomic_tree->{$chr}->{$strand} = $sub_tree;
        return $sub_tree;
    }
}

=method _normalise_strand

  Usage       : $strand = $tree->_normalise_strand( $strand );
  Purpose     : change strands like + or - to 1 or -1
  Returns     : Integer
  Parameters  : STRAND: Integer
  Throws      : 
  Comments    : 

=cut

sub _normalise_strand {
    my ( $self, $strand ) = @_;

    if ( $strand && ( $strand eq q{+} || $strand eq q{1} || $strand eq '+1' ) )
    {
        $strand = 1;
    }
    elsif ( $strand && ( $strand eq q{-} || $strand eq '-1' ) ) {
        $strand = -1;    ## no critic (ProhibitMagicNumbers)
    }
    else {
        $strand = 0;
    }

    return $strand;
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
            '1-10:1' => 'exon1.1',
            '11-29:1' => 'intron1.1',
            '30-40:1' => 'exon1.2',
        },
        2 => {
            '10-20:-1' => 'exon2.1',
            '21-25:-1' => 'intron2.1',
            '26-45:-1' => 'exon2.2',
            '21-35:-1' => 'intron2.2',
            '36-45:-1' => 'exon2.3',
        },
        'Zv9_scaffold1345' => {
            '10-20' => 'gc',
            '30-50' => 'gc',
        },
    };

    # add intervals to tree
    $genomic_tree->add_intervals_to_genomic_tree_from_hash( $test_intervals );

    # add a single interval
    $genomic_tree->insert_interval_into_tree( $chr, $start, $end, $strand,
        $object );

    # fetch overlapping intervals
    my $results_arrayref = $genomic_tree->fetch_overlapping_intervals( '1', 8,
        35, 1 );


=head1 DESCRIPTION

    Tree::GenomicIntervalTree is an implementation of an interval tree for
    genomic intervals. Basically, it contains a separate interval tree for each
    chromosome/scaffold. Intervals used in this implementation are closed. i.e.
    The interval 10-20 contains both the endpoints (10 and 20). [10,20]

=cut
