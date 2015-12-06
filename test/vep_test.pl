use strict;
use warnings;
 
use HTTP::Tiny;
use JSON;
use Data::Dumper;
 
my $http = HTTP::Tiny->new();

=begin comment

  Example output from VEP:
  [
    {
      'colocated_variants' => [
	{
	  'phenotype_or_disease' => 1,
	  'seq_region_name' => '9',
	  'eas_allele' => 'C',
	  'amr_maf' => '0.4553',
	  'strand' => 1,
	  'sas_allele' => 'C',
	  'id' => 'rs1333049',
	  'allele_string' => 'G/C',
	  'sas_maf' => '0.4908',
	  'amr_allele' => 'C',
	  'minor_allele_freq' => '0.4181',
	  'afr_allele' => 'C',
	  'eas_maf' => '0.5367',
	  'afr_maf' => '0.2133',
	  'end' => 22125504,
	  'eur_maf' => '0.4722',
	  'eur_allele' => 'C',
	  'minor_allele' => 'C',
	  'pubmed' => [
	    24262325,
	  ],
	  'start' => 22125504
	}
      ],
      'assembly_name' => 'GRCh38',
      'end' => 22125504,
      'seq_region_name' => '9',
      'transcript_consequences' => [
	{
	  'gene_id' => 'ENSG00000240498',
	  'variant_allele' => 'C',
	  'distance' => 4932,
	  'biotype' => 'antisense',
	  'gene_symbol_source' => 'HGNC',
	  'consequence_terms' => [
	    'downstream_gene_variant'
	  ],
	  'strand' => 1,
	  'hgnc_id' => 'HGNC:34341',
	  'gene_symbol' => 'CDKN2B-AS1',
	  'transcript_id' => 'ENST00000584020',
	  'impact' => 'MODIFIER'
	},
      ],
      'strand' => 1,
      'id' => 'rs1333049',
      'most_severe_consequence' => 'downstream_gene_variant',
      'allele_string' => 'G/C',
      'start' => 22125504
    }
  ]

=cut
 
my $server = "http://rest.ensembl.org";
my $ext = "/vep/human/id/rs1333049?content-type=application/json;Conservation=1";
my $response = $http->get($server.$ext);
die "Failed!\n" unless $response->{success};

if(length $response->{content}) {
  my $hash = decode_json($response->{content});
  local $Data::Dumper::Terse = 1;
  local $Data::Dumper::Indent = 1;
  print Dumper $hash;
  print "\n";
}
