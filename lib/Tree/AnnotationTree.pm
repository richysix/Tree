## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package Tree::AnnotationTree;
## use critic

# ABSTRACT: AnnotationTree object - used for holding and efficiently retrieving genomic annotation.

use namespace::autoclean;
use Moose;

extends 'Tree::GenomicIntervalTree';

=method new

  Usage       : my $annotation_tree = Tree::AnnotationTree->new();
  Purpose     : Constructor for creating AnnotationTree objects
  Returns     : Tree::AnnotationTree object
  Parameters  : 
  Throws      : If parameters are not the correct type
  Comments    : TO DO: Add support for strandedness

=cut

=method add_annotations_from_annotation_file

  Usage       : $annotation_tree->add_annotations_from_annotation_file( $annotation_file );
  Purpose     : adding annotations to the tree from a file
  Returns     : 1 on Success
  Parameters  : annotation_file String
  Throws      : If file cannot be opened
  Comments    : File must be tab-separated and of the form CHR  START   END ANNOTATION

=cut

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
    
    return 1;
}

=method add_annotations_from_gff

  Usage       : $annotation_tree->add_annotations_from_gff( $gff_file );
  Purpose     : adding annotations to the tree from a gff file
  Returns     : 1 on Success
  Parameters  : gff_file String
  Throws      : If file cannot be opened
  Comments    : File must be in gff format

=cut

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
    
    return 1;
}

=method fetch_overlapping_annotations

  Usage       : $annotation_tree->fetch_overlapping_annotations( $chr, $start, $end );
  Purpose     : retrieve overlapping annotations
  Returns     : ArrayRef
  Parameters  : CHR:    String
                START:  Integer
                END:    Integer
  Throws      : 
  Comments    : 

=cut

sub fetch_overlapping_annotations {
    my ( $self, $chr, $q_start, $q_end ) = @_;
    return $self->fetch_overlapping_intervals( $chr, $q_start, $q_end );
}

=method insert_annotation_into_genomic_tree

  Usage       : $annotation_tree->insert_annotation_into_genomic_tree( $chr, $start, $end, $annotation );
  Purpose     : add a single annotations to the tree
  Returns     : 1 on Success
  Parameters  : CHR:    String
                START:  Integer
                END:    Integer
                OBJECT: Any
  Throws      : 
  Comments    : 

=cut

sub insert_annotation_into_genomic_tree {
    my ( $self, $chr, $start, $end, $annotation ) = @_;
    $self->insert_interval_into_tree( $chr, $start, $end, $annotation, );
}

__PACKAGE__->meta->make_immutable;
1;

__END__

=pod

=head1 SYNOPSIS
 
    use Tree::AnnotationTree;
    my $annotation_tree = Tree::AnnotationTree->new()
    
    # add annotation to tree
    $annotation_tree->add_annotations_from_annotation_file( $annotation_file );
    $annotation_tree->add_annotations_from_gff( $gff_file );
    
    # add a single interval
    $annotation_tree->insert_annotation_into_genomic_tree( $chr, $start, $end, $annotation );
    
    # fetch overlapping intervals
    my $results_arrayref = $annotation_tree->fetch_overlapping_annotations( '1', 8, 35 );
  
  
=head1 DESCRIPTION

    Tree::AnnotationTree is an implementation of an interval tree for genomic annotations.
    It is a subclass of Tree::GenomicIntervalTree.
    Basically, it contains a separate interval tree for each chromosome/scaffold. The objects stored in the interval tree are the annotation.
    Intervals used in this implementation are closed. i.e. The interval 10-20 contains both the endpoints (10 and 20). [10,20]

=cut
