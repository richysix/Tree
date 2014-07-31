#!/usr/bin/env perl
# genomicintervaltree.t

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warn;

use List::MoreUtils qw{ any };
use Data::Dumper;

plan tests => 1 + 5 + 3 + 6 + 10 + 2 + 1 + 5;

use Tree::GenomicIntervalTree;

my $genomic_tree = Tree::GenomicIntervalTree->new();

# 1 test
isa_ok( $genomic_tree, 'Tree::GenomicIntervalTree' );

# check methods - 5 tests
my @methods = qw(
  genomic_tree
  add_intervals_to_genomic_tree_from_hash
  fetch_overlapping_intervals
  insert_interval_into_tree
  _make_new_subtree_for_chr
);

foreach my $method (@methods) {
    can_ok( $genomic_tree, $method );
}

my $test_intervals = {
    1 => {
        '1-10'  => 'exon1.1',
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
    3 => {
        '10-20' => 'exon3.1',
        '21-29' => 'intron3.1',
        '30-50' => 'exon3.2',
    },
};

# build interval tree
$genomic_tree->add_intervals_to_genomic_tree_from_hash($test_intervals);

# check some overlaps - 3 tests
my $results = $genomic_tree->fetch_overlapping_intervals( '1', 8, 35 );
my @expected_results = qw( exon1.1 intron1.1 exon1.2 );

foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 1' );
}

# check ends of interval - 6 tests
$results = $genomic_tree->fetch_overlapping_intervals( '1', 10, 10 );
foreach my $anno ( @{$results} ) {
    ok( $anno eq 'exon1.1', 'testing results 2' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '1', 11, 11 );
foreach ( @{$results} ) {
    ok( $_ eq 'intron1.1', 'testing results 3' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '1', 29, 29 );
foreach ( @{$results} ) {
    ok( $_ eq 'intron1.1', 'testing results 4' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '1', 30, 30 );
foreach ( @{$results} ) {
    ok( $_ eq 'exon1.2', 'testing results 5' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '1', 10, 11 );
@expected_results = qw( exon1.1 intron1.1 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 6' );
}

# chr2 - 10 tests
$results = $genomic_tree->fetch_overlapping_intervals( '2', 1, 1 );
ok( !@{$results}, 'testing results 7' );

$results = $genomic_tree->fetch_overlapping_intervals( '2', 10, 10 );
foreach ( @{$results} ) {
    ok( $_ eq 'exon2.1', 'testing results 8' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '2', 21, 21 );
@expected_results = qw( intron2.1 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 9' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '2', 25, 25 );
@expected_results = qw( intron2.1 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 10' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '2', 26, 26 );
@expected_results = qw( exon2.2 intron2.2 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 11' );
}

$results = $genomic_tree->fetch_overlapping_intervals( '2', 45, 45 );
@expected_results = qw( exon2.2 exon2.3 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 12' );
}

# insert new interval for a current and a new chr
$genomic_tree->insert_interval_into_tree( '1', 5, 15, 'lincRNA-1' );
$genomic_tree->insert_interval_into_tree( 'Zv9_scaffold4', 1, 100,
    'promoter-1' );

#check new intervals - 2 tests
$results = $genomic_tree->fetch_overlapping_intervals( '1', 10, 10 );
@expected_results = qw( exon1.1 lincRNA-1 );
foreach my $anno ( @{$results} ) {
    ok( ( any { $anno eq $_ } @expected_results ), 'testing results 13' );
}

# check _make_new_subtree_for_chr throws correctly - 1 test
warning_like { $genomic_tree->_make_new_subtree_for_chr('1') }
qr/Tree\salready\sexists!/xms, '_make_new_subtree warning correct';

# check add_intervals_to_genomic_tree_from_hash throws correctly - 5 tests
throws_ok { $genomic_tree->add_intervals_to_genomic_tree_from_hash('scalar') }
qr/Supplied\sobject\sdoes\snot\smatch\srequired\sstructure/xms,
  'Scalar input to add_intervals_to_genomic_tree_from_hash';
my @array;
throws_ok { $genomic_tree->add_intervals_to_genomic_tree_from_hash(@array) }
qr/Supplied\sobject\sdoes\snot\smatch\srequired\sstructure/xms,
  'Array input to add_intervals_to_genomic_tree_from_hash';
throws_ok { $genomic_tree->add_intervals_to_genomic_tree_from_hash( \@array ) }
qr/Supplied\sobject\sdoes\snot\smatch\srequired\sstructure/xms,
  'ArrayRef input to add_intervals_to_genomic_tree_from_hash';
throws_ok {
    $genomic_tree->add_intervals_to_genomic_tree_from_hash( { 1 => 'scalar' } );
}
qr/Supplied\sobject\sdoes\snot\smatch\srequired\sstructure/xms,
  'Chr entry not hash input to add_intervals_to_genomic_tree_from_hash';
throws_ok {
    $genomic_tree->add_intervals_to_genomic_tree_from_hash(
        { 1 => { 'scalar' => 'scalar' }, } );
}
qr/Supplied\sobject\sdoes\snot\smatch\srequired\sstructure/xms,
  'Chr entry not correct format to add_intervals_to_genomic_tree_from_hash';
