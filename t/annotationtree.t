#!/usr/bin/env perl
# annotationtree.t

use strict;
use warnings;

use Test::More;

use List::MoreUtils qw{ any };
use Data::Dumper;

plan tests => 1 + 9 + 3 + 6 + 10 + 2;

use Tree::AnnotationTree;

my $annotation_tree = Tree::AnnotationTree->new();

# 1 test
isa_ok( $annotation_tree, 'Tree::AnnotationTree' );

# check methods - 9 tests
my @methods = qw(
    genomic_tree
    add_intervals_to_genomic_tree_from_hash
    fetch_overlapping_intervals
    insert_interval_into_tree
    _make_new_subtree_for_chr
    add_annotations_from_annotation_file
    add_annotations_from_gff
    fetch_overlapping_annotations
    insert_annotation_into_genomic_tree
);

foreach my $method (@methods) {
    can_ok( $annotation_tree, $method );
}

my $test_annotation = << "END_ANNOTATION";
1	1	10	exon1.1
1	11	29	intron1.1
1	30	40	exon1.2
2	10	20	exon2.1
2	21	25	intron2.1
2	26	45	exon2.2
2	21	35	intron2.2
2	36	45	exon2.3
3	10	20	exon3.1
3	21	29	intron3.1
3	30	50	exon3.2
END_ANNOTATION

my $test_annotation_file = 'test_annotation.txt';
open my $test_fh, '>', $test_annotation_file
  or die "Couldn't open test annotation file!\n";
print $test_fh $test_annotation;
close $test_fh;

# build interval tree
$annotation_tree->add_annotations_from_annotation_file($test_annotation_file);

# check some overlaps - 3 tests
my $results = $annotation_tree->fetch_overlapping_annotations( '1', 8, 35 );
my @expected_results = qw( exon1.1 intron1.1 exon1.2 );

foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 1' );
}

# check ends of interval - 6 tests
$results = $annotation_tree->fetch_overlapping_annotations( '1', 10, 10 );
foreach my $anno ( @{$results} ) {
    ok( $anno eq 'exon1.1', 'testing results 2' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '1', 11, 11 );
foreach ( @{$results} ) {
    ok( $_ eq 'intron1.1', 'testing results 3' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '1', 29, 29 );
foreach ( @{$results} ) {
    ok( $_ eq 'intron1.1', 'testing results 4' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '1', 30, 30 );
foreach ( @{$results} ) {
    ok( $_ eq 'exon1.2', 'testing results 5' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '1', 10, 11 );
@expected_results = qw( exon1.1 intron1.1 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 6' );
}

# chr2 - 10 tests
$results = $annotation_tree->fetch_overlapping_annotations( '2', 1, 1 );
ok( !@{$results}, 'testing results 7' );

$results = $annotation_tree->fetch_overlapping_annotations( '2', 10, 10 );
foreach ( @{$results} ) {
    ok( $_ eq 'exon2.1', 'testing results 8' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '2', 21, 21 );
@expected_results = qw( intron2.1 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 9' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '2', 25, 25 );
@expected_results = qw( intron2.1 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 10' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '2', 26, 26 );
@expected_results = qw( exon2.2 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 11' );
}

$results = $annotation_tree->fetch_overlapping_annotations( '2', 45, 45 );
@expected_results = qw( exon2.2 exon2.3 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 12' );
}

my @test_gff = qw{ 1 lncRNAs lincRNA-1 5 15 0.9 + . . };

my $test_gff_file = 'test_annotation.gff';
open my $gff_fh, '>', $test_gff_file or die "Couldn't open test gff file!\n";
print $gff_fh join( "\t", @test_gff, ), "\n";
close $gff_fh;

# build interval tree
$annotation_tree->add_annotations_from_gff($test_gff_file);

# check some overlaps - 2 tests
$results = $annotation_tree->fetch_overlapping_annotations( '1', 10, 10 );
@expected_results = qw( exon1.1 lincRNA-1 );

foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 13' );
}

unlink( $test_annotation_file, $test_gff_file, );
