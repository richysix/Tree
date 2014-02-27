package Tree::AnnotationTree;
use namespace::autoclean;
use Moose;

extends 'Tree::GenomicIntervalTree';

sub add_annotations_from_annotation_file {
    my ( $self, $annotation_file ) = @_;
    
    # open annotation file - should be tab-separated and of the form: chr    start    end   annotation
    open my $anno_fh, '<', $annotation_file or die "Couldn't open annotation file, $annotation_file:$!\n";
    while(<$anno_fh>){
        chomp;
        my( $chr, $start, $end, $annotation, ) = split /\t/, $_;
        if( !exists $self->genomic_tree->{$chr} ){
            $self->_make_new_subtree_for_chr( $chr );
        }
        $self->insert_annotation_into_genomic_tree( $chr, $start, $end, $annotation );
    }
}

sub add_annotations_from_gff {
    my ( $self, $annotation_file ) = @_;
    
    # open annotation file - should be gff format:  chr    source annotation  start    end   score    strand  frame   attribute
    open my $anno_fh, '<', $annotation_file or die "Couldn't open annotation file, $annotation_file:$!\n";
    while(<$anno_fh>){
        chomp;
        my( $chr, undef, $annotation, $start, $end, undef, ) = split /\t/, $_;
        if( !exists $self->genomic_tree->{$chr} ){
            $self->_make_new_subtree_for_chr( $chr );
        }
        $self->insert_annotation_into_genomic_tree( $chr, $start, $end, $annotation );
    }
}

sub fetch_overlapping_annotations {
    my ( $self, $chr, $q_start, $q_end ) = @_;
    return $self->fetch_overlapping_intervals( $chr, $q_start, $q_end );
}

sub insert_annotation_into_genomic_tree {
    my ( $self, $chr, $start, $end, $annotation ) = @_;
    $self->insert_interval_into_tree( $chr, $start, $end, $annotation, );
}

__PACKAGE__->meta->make_immutable;
1;
