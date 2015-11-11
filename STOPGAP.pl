#!/usr/bin/env perl

use strict;
use JSON;
use Data::Dumper;

my $DEBUG = 0;
my $FANTOM5_LOCATION = undef;

# TODO: OMIM http://omim.org/api
# TODO: MeSH: http://ws.ncibi.org/g2m.html
# TODO: Regulome
# TODO: PCHIC
# TODO: DHS
# TODO: GRASP: https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip
# TODO: Phewas Catalog: http://phewas.mc.vanderbilt.edu/phewas-catalog.csv
# TODO: GWAS DB: http://jjwanglab.org:8080/gwasdb/GWASdb_snp_v4.zip
# TODO: CATO score (Maurano 2015) and gkm-SVM (Ghandi et al 2014)
# TODO: Enhance PIC with epigenomic priors

main();

sub main {
  # TODO
}

=head2 diseases_to_genes

  Associates genes from a list of diseases
  Args: Listref of Bio::EnsEMBL::Variation::Phenotype
  Returntype: Listref of hashrefs {
                gene => Bio::EnsEMBL::Gene
                disease => Bio::EnsEMBL::Variation::Phenotype
		score => scalar
                evidence => listref of hashrefs returned by markers_to_genes
              }

=cut

sub diseases_to_genes {
  my ($diseases, $populations, $tissues, $parameters) = @_;
  markers_to_genes(diseases_to_markers($diseases), $populations, $tissues);
}

=head2 diseases_to_markers

  Associates markers from a list of diseases
  Args: Listref of Bio::EnsEMBL::Variation::Phenotype
  Returntype: Listref of hashrefs {
                pvalue => scalar
                marker => Bio::EnsEMBL::Variation::Variation
                disease => Bio::EnsEMBL::Variation::Phenotype
                source => string
              }

=cut

sub diseases_to_markers {
  my ($diseases) = @_;

  # Extract crude data
  my $markers = scan_disease_databases($diseases);

  # Filter by p-value
  @filtered_markers = grep { $_->{pvalue} < $PVALUE_CUTOFF} @$markers;

  return \@filtered_markers;
}

=head2 scan_disease_databases 

  Associates markers from a list of diseases
  Args: Listref of Bio::EnsEMBL::Variation::Phenotype
  Returntype: Listref of hashrefs {
                pvalue => scalar
                marker => Bio::EnsEMBL::Variation::Variation
                disease => Bio::EnsEMBL::Variation::Phenotype
                source => string
              }

=cut

sub scan_disease_databases {
  my ($diseases) = @_;

  my @res = ();
  foreach my $disease (@$diseases) {
    push @res, scan_disease_databases_2($disease);
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." markers associated to all diseases\n";
  } 

  return \@res;
}

=head2 scan_disease_databases_2

  Associates markers from a disease
  Args: Bio::EnsEMBL::Variation::Phenotype
  Returntype: Listref of hashrefs {
                pvalue => scalar
                marker => Bio::EnsEMBL::Variation::Variation
                disease => Bio::EnsEMBL::Variation::Phenotype
                source => string
              }

=cut

my @database_functions = (\&GWASCatalog);

sub scan_disease_databases_2 {
  my ($disease) = @_;

  my @res = ();
  foreach my $scan_disease_database (@database_functions) {
    push @res, scan_disease_database($disease);
  }
    
  return \@res;
}

=head2 GWASCatalog

  Returns all markers associated to a disease in GWAS Catalog
  Arg: Trait name or EFO ID
  Returntype: listref of Bio::EnsEMBL::Variation::Variation objects

=cut

sub GWASCatalag {
  my ($disease) = @_
  my $server = "http://www.ebi.ac.uk/gwas/";
  my $ext = "search?";
  my $response = $http->get($server.$ext, {headers => { 'Content-type' => 'application/json' }});
  die "Failed!\n" unless $response->{success};
  if(length $response->{content}) {
    my $hash = decode_json($response->{content});
    # TODO
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." markers associated to disease $disease in GWAS Catalog\n";
  } 

  return \@res;
}

=head2 markers_to_genes

  Associates Genes to markers of interest
  Args: 
  * Listref produced by diseases_to_markers
  * Listref of Bio::EnsEMBL::Variation::Population objects
  * Listref of strings (tissue names)
  Returntype: Listref of hashrefs {
                gene => Bio::EnsEMBL::Gene
                marker_evidence => hashref returned by disease_to_markers
                score => scalar
                association_evidence => listref of hashrefs returned by snp_to_genes
              }

=cut

sub markers_to_genes {
  my ($markers, $populations, $tissues) = @_;

  # Must set the population settings before separating out the markers
  if (! defined $populations) { 
    $populations = markers_to_populations($markers);
  }

  # Must set the tissue settings before separating out the markers
  if (! defined $tissues) {
    $tissues = markers_to_tissues($markers);
  }

  # Compute for each marker separately, concatenate results
  my @res = map { @$_ } map { marker_to_genes($_, $populations, $tissues) } @$markers;

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to all markers\n";
  } 

  return \@res;
}

=head2 markers_to_populations

  Associates list of populations to list of markers
  Args: 
  * Hashref produced by diseases_to_markers
  Returntype: Listref of Bio::EnsEMBL::Variation::Population objects

=cut

sub markers_to_populations {
  my ($markers) = @_;
  return undef; # TODO
}

=head2 markers_to_tissues

  Associates list of tissues to list of markers
  Args: 
  * Hashref produced by diseases_to_markers
  Returntype: Listref of strings 

=cut

sub markers_to_tissues {
  my ($markers) = @_;
  return undef; # TODO See FORGE??
}

=head2 marker_to_genes

  Associates Genes to marker of interest
  Args: 
  * Hashref produced by diseases_to_markers
  * Listref of Bio::EnsEMBL::Variation::Population objects
  * Listref of strings (tissue names)
  Returntype: Listref of hashrefs {
                gene => Bio::EnsEMBL::Gene
                marker_evidence => hashref returned by disease_to_markers
                score => scalar
                association_evidence => listref of hashrefs returned by snp_to_genes
              }

=cut

sub marker_to_genes {
  my ($marker, $populations, $tissues) = @_;
  my %res = ();
  foreach my $snp_gene_pair (snps_to_genes(marker_to_snps($marker, $populations), $tissues)) {
    my $key = $snp_gene_pair->{gene}->stable_id;
    $res{ $key }->gene = $snp_gene_pair->{gene};
    $res{ $key }->marker_evidence = $marker;
    $res{ $key }->score = # TODO ???
    my @evidence_list = @{$res{ $key }->association_evidence};
    push @evidence_list, $snp_gene_pair;
  }
  my @values = values %res;

  if ($DEBUG) {
    printf "Found ".scalar @values." genes associated to marker ".$marker->{marker}->name."\n";
  } 

  return \@values;
}

=head2 marker_to_snps

  Associates LD linked SNPs to markers of interest
  Args: 
  * Hashref produced by diseases_to_markers
  * Listref of Bio::EnsEMBL::Variation::Population objects
  Returntype: Listref of hashrefs {
                snp => Bio::EnsEMBL::Variation::Variation
                ld => scalar
                marker_evidence => hashref returned by diseases_to_markers 
              }

=cut

sub marker_to_snps {
  my ($marker, $populations) = @_;
  my $rsquares = $marker->get_all_LD_values->get_all_r_square_values();
  my @filtered_rsquares = grep {exists $populations->{$_->{population_id}}} $rsquares;
  
  my @res = ();
  foreach my $filtered_rsquare (@filtered_rsquares) {
    if ($filtered_rsquare->{variation1} == $marker) {
      push @res, {
        snp => $filtered_rsquare->{variation2}, 
        ld => $filtered_rsquare->{r_square},
        marker_evidence => $marker,
      }
    } else {
      push @res, {
        snp => $filtered_rsquare->{variation1}, 
        ld => $filtered_rsquare->{r_square},
        marker_evidence => $marker,
      }
    }
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." LD SNPs associated to marker ".$marker->{marker}->name."\n";
  } 

  return \@res;
}

sub cluster_snps {
  my ($snps, $populations) = @_;
  # TODO
  return \@ld_snps;
}

=head2 snps_to_genes

  Associates genes to listref of LD linked SNPs
  Args: 
  * Listref of hashrefs returned by marker_to_snps
  * Listref of strings (tissues)
  Returntype: Listref of hashrefs {
    snp => Hashref returned by marker_to_snps,
    gene => Bio::EnsEMBL::Gene,
    tissue => tissue,
    score => scalar,
    source => string
  }

=cut

sub snps_to_genes {
  my ($snps, $tissues) = @_;
  my @res = map { $_ } map { snp_to_genes->($_, $tissues) } @$snps;

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to all SNPs\n";
  } 

  return \@res;
}

=head2 snp_to_genes

  Associates genes to LD linked SNP
  Args: 
  * Listref of hashrefs returned by marker_to_snps
  * Listref of strings (tissues)
  Returntype: Listref of hashrefs {
    snp => Hashref returned by marker_to_snps,
    gene => Bio::EnsEMBL::Gene,
    tissue => tissue,
    score => scalar,
    evidence => listref of hashrefs returned by database queries
  }

=cut

my @snp_to_gene_functions = (\&GTEx, \&fantom5, \&VEP);

sub snp_to_genes {
  my ($snp, $tissues) = @_;
  # @snp_to_gene_functions is an array of function pointers
  # We compute the result of each of these functions, which creates an array of arrayrefs
  # Each arrayref is flattened into an array, in effect concatenating all arrays into one.
  my @evidence = map { @$_ } map { $_->($snp, $gene, $tissues) } @snp_to_gene_functions;

  my %res = ();
  foreach my $record (@evidence) {
    my $gene = $record->{gene};
    $res{$gene}{snp} = $snp;
    $res{$gene}{gene} = $gene;
    $res{$gene}{v2g_score} += $record->{score};
    $res{$gene}{tissue} = $tissues; # TODO separate by tissues eventually
    push @{$res{$gene}{evidence}}, $record;
  }

  foreach my $gene (keys %res) {
    $res{$gene}{gene_score} = $snp->{ld} * $res{$gene}{v2g_score};
    $res{$gene}{pics_gene_score} = total_score(XXX, $res{$gene}{gene_score}); # TODO PICS?
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP ".$snp->{snp}->name."\n";
  } 

  return \@res;
}

=head2 total_score

  Computes a weird mean function from snp PICs score and Gene/SNP association score
  Args: 
  * PICS, scalar
  * gene_score, scalar

=cut

sub total_score {
  my ($pics, $gene_score) = @_;
  
  my $A = $pics * ($pics ** (1/3));
  my $B = $gene_score * ($gene_score ** (1/3));
  return (($A + $B) / 2) ** 3;
}

=head2 GTEx

  Returns all genes associated to a SNP in GTEx 
  Arg: Hashref returned by marker_to_snps
  Returntype: listref of hashrefs {
    snp => Hashref returned by marker_to_snps,
    gene => Bio::EnsEMBL::Gene,
    tissue => tissue,
    score => scalar,
    source => string
  }

=cut

sub GTEx {
  my ($snp) = @_
  my $server = "http://rest.ensembl.org";
  my $ext = "eqtl?"; # TODO
  my $response = $http->get($server.$ext, {headers => { 'Content-type' => 'application/json' }});
  die "Failed!\n" unless $response->{success};
  if(length $response->{content}) {
    my $hash = decode_json($response->{content});
    # TODO
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$snp->name.' in GTEx\n";
  } 

  return \@res;
}

=head2 VEP

  Returns all genes associated to a SNP in VEP 
  Arg: Hashref returned by marker_to_snps
  Returntype: listref of hashrefs {
    snp => Hashref returned by marker_to_snps,
    gene => Bio::EnsEMBL::Gene,
    tissue => tissue,
    score => scalar,
    source => string
  }

=cut

sub VEP {
  my ($snp) = @_
  my $server = "http://rest.ensembl.org";
  my $ext = "eqtl?"; # TODO
  my $response = $http->get($server.$ext, {headers => { 'Content-type' => 'application/json' }});
  die "Failed!\n" unless $response->{success};
  if(length $response->{content}) {
    my $hash = decode_json($response->{content});
    # TODO
  }

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$snp->name.' in VEP\n";
  } 

  return \@res;
}

=head2 Fantom5

  Returns all genes associated to a SNP in Fantom5 
  Arg: Hashref returned by marker_to_snps
  Returntype: listref of hashrefs {
    snp => Hashref returned by marker_to_snps,
    gene => Bio::EnsEMBL::Gene,
    tissue => tissue,
    score => scalar,
    source => string
  }

=cut

sub Fantom5 {
  my ($snp) = @_

  if ($DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$snp->name.' in Fantom5\n";
  } 

  return \@res;
}
