package Tree::GenomicIntervalTree;
use namespace::autoclean;
use Moose;
use Carp qw(cluck confess);
use Readonly;

use Set::IntervalTree;

has 'genomic_tree' => (
	is => 'ro',
    isa => 'HashRef',
    default => sub { {} },
);


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
    if( ref $hash ne 'HASH' ){
        confess "Supplied object does match required structure!\n\n", $HASH_STRUCTURE;
    }
    foreach my $chr ( keys %{$hash} ){
        if( ref $hash->{$chr} ne 'HASH' ){
            confess "Supplied object does match required structure!\n\n", $HASH_STRUCTURE;
        }
        foreach my $interval ( keys %{$hash->{$chr}} ){
            if( $interval !~ m/[0-9]+\-[0-9]+/ ){
                confess "Supplied object does match required structure!\n\n", $HASH_STRUCTURE;
            }
        }
    }
    
    foreach my $chr ( keys %{$hash} ){
        foreach my $interval ( keys %{$hash->{$chr}} ){
            my ( $start, $end ) = split /-/, $interval;
            my $object = $hash->{$chr}->{$interval};
            $self->insert_interval_into_tree( $chr, $start, $end, $object, );
        }
    }
}

sub fetch_overlapping_intervals {
    my ( $self, $chr, $q_start, $q_end ) = @_;
    
    if( exists $self->genomic_tree->{$chr} ){
        my $results = $self->genomic_tree->{$chr}->fetch( $q_start, $q_end );
        return $results;
    }
    else{
        return [];
    }
}

sub insert_interval_into_tree {
    my ( $self, $chr, $start, $end, $object ) = @_;
    if( !exists $self->genomic_tree->{$chr} ){
        $self->_make_new_subtree_for_chr( $chr );
    }
    $self->genomic_tree->{$chr}->insert( $object, $start - 1, $end + 1);
}

sub _make_new_subtree_for_chr {
    my ( $self, $chr ) = @_;
    
    # check tree doesn't already exist
    if( exists $self->genomic_tree->{$chr} ){
        warn "Tree already exists!\n";
        return $self->genomic_tree->{$chr};
    }
    else{
        my $sub_tree = Set::IntervalTree->new;
        $self->genomic_tree->{$chr} = $sub_tree;
        return $sub_tree;
    }
}


__PACKAGE__->meta->make_immutable;
1;
