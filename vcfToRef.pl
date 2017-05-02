#!/usr/bin/perl
#version 1.01
#Changed fasta heading format

=head1 vcfToRef

vcfToRef - A simple program to build reference list based on input VCF file. Internet connection required. Reference list is returned to STDOUT in Fasta format.

=head1 SYNOPSIS

  vcfToRef [options] --vcf <inputFile>
  
  Options:
   --help       Print this message and exit.
   --span       Specify spanning length and respctively reference length. Default 250.
   --build      Specify genome build version of VCF file. Default hg19.

=cut

use strict;
use warnings;
use String::Util qw(trim);
use List::Util qw(max);
use XML::Simple;
use LWP::UserAgent;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

$| = 1;
my $lengthConst	= 250;
my $build	= 'hg19';
my $inputVcf;
my $help	= 0;

GetOptions(
	'vcf=s'	=> \$inputVcf,
	'build=s' => \$build,
	'help!' => \$help,
	'span=i' => \$lengthConst) or pod2usage(1);
pod2usage(1) if $help;
pod2usage(1) unless defined $inputVcf;
main();

sub vcfSortRule {
	return -1 if (($a->[1] cmp $b->[1]) eq -1);
	return 1 if (($a->[1] cmp $b->[1]) eq 1);
	return -1 if $a->[2] < $b->[2];
	return 1 if $a->[2] > $b->[2];
	return 0;
	}

sub getSeq {
	my $build = shift;
	my $chr = shift;
	my $pos1 = shift;
	my $pos2 = shift;
	my $try = 0;
	my $ua = LWP::UserAgent->new;
	FETCH:
		my $response = $ua->get("http://genome.ucsc.edu/cgi-bin/das/$build/dna?segment=$chr:$pos1,$pos2");
		unless ($response && $response->is_success) {
			++$try;
			die "Cant reach http://genome.ucsc.edu\n" if $try > 100;
			goto FETCH;
			}
	my $xmlString = $response->content or return '-1';
	my @options = ();
	my $ref = XMLin($xmlString, @options) or die "$chr:$pos1-$pos2 falls out from $build.";
	my $return = trim((((($ref)->{"SEQUENCE"})->{"DNA"})->{"content"}));
	$return =~ s/\n//g;
	return uc($return);
	}

sub mergeBed {
	my $vcf = shift;
	my $bed = shift;
	return ($bed, $vcf) if scalar @{$bed} eq 0;
	my $current = 0;
	my $bedResult;
	push (@{$bedResult}, $bed->[$current]);
	for (my $i = 1; $i < scalar @{$bed}; $i++) {
		if ($bed->[$i]->[1] eq $bed->[$current]->[1]) {
			if ($bed->[$i]->[2] <= $bed->[$current]->[3]) {
				$bed->[$current]->[3] = $bed->[$i]->[3];
				$vcf->[$i]->[0] = $current;
				next;
				}
			}
		$current = $i;
		push(@{$bedResult}, $bed->[$current]);
		}
	return ($vcf, $bedResult);
	}

sub mutate {
	my $seq = shift;
	my $seqPos = shift;
	my $mutPos = shift;
	my $mutRef = shift;
	my $mutAlt = shift;
	my $mutName = shift;
	my $result = substr $seq, $mutPos - $seqPos, length($mutRef), $mutAlt;
	print STDERR "Warning: reference allele for mutation '$mutName' does not correspond to the genome sequence. Check build version\n" unless $result eq $mutRef;
	return $seq;
	}


sub main {
	open (READ, "<$inputVcf") or die "$inputVcf: no such file\n";
	if (($build ne 'hg19')and($build ne 'hg38')) {die "Only GRCh37 and GRCh38 genome builds are supported. Use either 'hg19' or 'hg38' to specify build"}
	my $vcf;
	my $bed;
	while (<READ>) {
		next if m!^#!;
		my @mas = split/\t/;
		die "Wrong VCF file" unless defined $mas[4];
		my @var = split/,/, $mas[4];
		foreach my $arg (@var) {
			if ((length($mas[3]) > $lengthConst)or(length($arg) > $lengthConst)) {
				$lengthConst = max(length($mas[3]), length($arg)) + 10;
				print STDERR "Warning: spanning length was adjusted according to the substitution length: $lengthConst\n";
				}
			push (@{$vcf}, [0, $mas[0], $mas[1], $mas[3], $arg, $mas[2]]);
			}
		}
	@{$vcf} = sort vcfSortRule @{$vcf};
	for (my $i = 0; $i < scalar @{$vcf}; $i++) {
		$vcf->[$i]->[0] = $i;
		push (@{$bed},
			[$i, 
			$vcf->[$i]->[1], 
			$vcf->[$i]->[2] - $lengthConst, 
			$vcf->[$i]->[2] + $lengthConst]);
		}
	($vcf, $bed) = mergeBed($vcf, $bed);
	my $ref;
	for (my $i = 0; $i < scalar(@{$bed}); $i++) {
		my $seq = getSeq($build, $bed->[$i]->[1], $bed->[$i]->[2], $bed->[$i]->[3]);
		if ($seq eq '-1') {
			print STDERR  "Warning: ",$bed->[$i]->[1],":",$bed->[$i]->[2],"-",$bed->[$i]->[3]," falls out from $build.\n";
			}
		$ref->{$bed->[$i]->[0]} = [@{($bed->[$i])}, $seq];
		}
	my $mut;
	for (my $i = 0; $i < scalar(@{$vcf}); $i++) {
		next if $ref->{$vcf->[$i]->[0]}->[4] eq '-1';
		$mut->{$vcf->[$i]->[0] . "|mut|" . 
			$vcf->[$i]->[5] . "|" . 
			$vcf->[$i]->[1] . "|" . 
			$vcf->[$i]->[2] . "|" . 
			$vcf->[$i]->[3] . "|" . 
			$vcf->[$i]->[4]} = mutate($ref->{$vcf->[$i]->[0]}->[4],
							$ref->{$vcf->[$i]->[0]}->[2],
							$vcf->[$i]->[2],
							$vcf->[$i]->[3],
							$vcf->[$i]->[4],
							$vcf->[$i]->[5]);
		}
	foreach my $key (keys %{$ref}) {
		next if $ref->{$key}->[4] eq '-1';
		print ">AMPL",$key,"|wt|",$ref->{$key}->[1],"|",$ref->{$key}->[2],"|",$ref->{$key}->[3],"\n",$ref->{$key}->[4],"\n";
		}
	foreach my $key (keys %{$mut}) {
		print ">AMPL",$key,"\n",$mut->{$key},"\n";
		}
	}


















