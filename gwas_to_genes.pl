#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use Bio::EnsEMBL::Registry;
use JSON;
use XML::Simple qw(:strict);

our $registry = 'Bio::EnsEMBL::Registry';
our $SPECIES = 'Human';
our $DEBUG = 1;
our $DATABASES_DIR = "databases";
our %phenotype_cache = ();

main();

=begin comment

  Datasets to integrate:
    Cis-regulatory annotations:
      Regulome (STOPGAP Scoring: 1-2: +2, 3: +1, 4: +0)
      PCHIC (STOPGAP Scoring: Single cell line: +1, multiple cell lines: +2)
      DHS associations (FPR 0-0.6: +2, 0.6-0.85: +1,0.85-1: +0)
      PhyloP (STOPGAP Scoring: FPR 0-0.6: +2, 0.6-0.85: +1,0.85-1: +0)

    Epigenetic activity:
      DHS
      Fantom5

  Missing details:
    DHScor FDR % distance

  Improvements:
    Ontological search for phenotypes? EFO?
    Replace PICS
    Epigenetic priors

=cut 

=head2 main

  Reads commandline parameters, prints corresponding associated genes with evidence info

=cut

sub main {
  my $options = get_options();
  print encode_json(diseases_to_genes($options->{diseases}, $options->{population_weights}, $options->{tissues}));
}

=head2 get_options

  Reads commandline parameters
  Returntype: {
                diseases => [ string ],
                population_weights => [ Bio::EnsEMBL::Variation::Population ],
                tissues => [ string ],
              }

=cut

sub get_options {
  my %options = ();

  GetOptions(\%options, "help=s", "disease|d=s@", "population|p=s@", "tissue|t=s@", "databases|D=s", "host=s", "user=s", "port=s", "species=s", "debug|g");

  $options{diseases} = $options{disease};
  $options{population_weights} = map { $population_adaptor->fetch_by_name($_) } @$options{population};
  $options{tissues} = $options{tissue};
  
  # Connecting to Ensembl registry
  $options->{host} ||= 'ensembldb.ensembl.org';
  $options->{user} ||= 'anonymous';
  our $registry;
  $registry->load_registry_from_db(
    -host => $options->{host},
    -user => $options->{user},
    -port => $options->{port},
  ); 

  if (defined $options{databases}) {
    our $DATABASES_DIR = $options{databases};
  }

  if (defined $options->{species}) {
    our $SPECIES = $options->{species};
  }

  if (defined $options->{$debug}) {
    our $DEBUG = 1;
  }

  return \%options;
}

=head2 diseases_to_genes

  Associates genes from a list of diseases
  Args: [ string ] 
  Returntype: [  
                {
                  gene => Bio::EnsEMBL::Gene,
                  score => scalar, # aka pics.gene.score
                  MeSH => [ string ],
                  gene_phenotype_associations => [ Bio::EnsEMBL::Variation::PhenotypeFeature ],
                  gwas_snp_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      pvalue => scalar,
                      evidence => [ 
                        {
                          pvalue => scalar,
                          disease => string,
                          source => string,
                          study => string,
                        }
                      ],
                    }
                  ],
                  association_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      score => scalar, # aka v2g.score
                      evidence => [
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ],
                    }
                  ],
                }
              ]

=cut

sub diseases_to_genes {
  my ($diseases, $population_weights, $tissues) = @_;
  my $res = gwas_snps_to_genes(diseases_to_gwas_snps($diseases), $population_weights, $tissues);

  # Add some gene/disease external annotations
  map { $_->{MeSH} = disease_to_MeSH($diseases) } @$res;
  map { $_->{gene_phenotype_association} = gene_to_phenotypes($_->{gene})  } @$res;

  return $res;
}

=head2 diseases_to_gwas_snps

  Associates gwas_snps from a list of diseases
  Args: [ string ]
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::Variation,
                  pvalue => scalar, (min of pvalues - could be improved)
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

sub diseases_to_gwas_snps {
  my ($diseases) = @_;

  # Extract crude data
  my $gwas_snps = scan_disease_databases($diseases);

  # Filter by p-value
  @filtered_gwas_snps = grep { $_->{pvalue} < $PVALUE_CUTOFF } @$gwas_snps;

  return \@filtered_gwas_snps;
}

=head2 scan_disease_databases 

  Associates gwas_snps from a list of diseases
  Args: [ Bio::EnsEMBL::Variation::Phenotype ]
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::Variation,
                  pvalue => scalar, (min of pvalues - could be improved)
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

sub scan_disease_databases {
  my ($diseases) = @_;

  my @hits = map { @$_ } map( scan_disease_databases_2, @$diseases)

  my %hash = ();
  foreach my $hit (@hits) {
    if (exists $hash{$hits->{snp}->name}) {
      my $record = $hash{$hits->{snp}->name};
      push @{$record->{evidence}}, $hit;
      if ($record->{pvalue} > $hit->{pvalue}) {
        $record->{pvalue} = $hit->{pvalue};
      }
    } else {
       $hash{$hits->{snp}->name}->{evidence} = {
         snp => $hit->{snp},
         pvalue => $hit->{pvalue},
         evidence => [ $hit ],
       }
    }
  }
  
  @res = values %hash;

  if (our $DEBUG) {
    printf "Found ".scalar @res." gwas_snps associated to all diseases\n";
  } 

  return \@res;
}

=head2 scan_disease_databases_2

  Associates gwas_snps from a disease
  Args: string 
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::Variation,
                  pvalue => scalar, (min of pvalues - could be improved)
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

my @database_functions = (\&GWASCatalog, \&GRASP, \&GWAS_DB, \&Phewas_Catalog);

sub scan_disease_databases_2 {
  my ($disease) = @_;

  our $registry;
  our $SPECIES;
  my $variation_adaptor = $registry->get_adaptor($SPECIES, 'Variation', 'Variation');

  my @hits = map { @$_ } map { $_->($disease, $variation_adaptor) } @database_functions;

  my %hash = ();
  foreach my $hit (@hits) {
    if (exists $hash{$hits->{snp}->name}) {
      my $record = $hash{$hits->{snp}->name};
      push @{$record->{evidence}}, $hit;
      if ($record->{pvalue} > $hit->{pvalue}) {
        $record->{pvalue} = $hit->{pvalue};
      }
    } else {
       $hash{$hits->{snp}->name}->{evidence} = {
         snp => $hit->{snp},
         pvalue => $hit->{pvalue},
         evidence => [ $hit ],
       }
    }
  }
  
  my @res = values %hash;
    
  return \@res;
}

=head2 GWASCatalog

  Returns all gwas_snps associated to a disease in GWAS Catalog
  Arg: Trait name or EFO ID (string)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  source => "GWAS Catalog",
                  study => string,
                }
              ]

=cut

=begin comment

  GWAS Catalog flat file format (waiting for REST API...)

  1.  DATE ADDED TO CATALOG: Date added to catalog
  2.  PUBMEDID: PubMed identification number
  3.  FIRST AUTHOR: Last name of first author
  4.  DATE: Publication date (online (epub) date if available)
  5.  JOURNAL: Abbreviated journal name
  6.  LINK: PubMed URL
  7.  STUDY: Title of paper (linked to PubMed abstract)
  8.  DISEASE/TRAIT: Disease or trait examined in study
  9.  INITIAL SAMPLE DESCRIPTION: Sample size for Stage 1 of GWAS
  10.  REPLICATION SAMPLE DESCRIPTION: Sample size for subsequent replication(s)
  11.  REGION: Cytogenetic region associated with rs number (NCBI)
  12.  CHR_ID: Chromosome number associated with rs number (NCBI)
  13.  CHR_POS: Chromosomal position associated with rs number (dbSNP Build 144, Genome Assembly GRCh38.p2, NCBI)
  14.  REPORTED GENE (S): Gene(s) reported by author
  15.  MAPPED GENE(S): Gene(s) mapped to the strongest SNP (NCBI). If the SNP is located within a gene, that gene is listed. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen.
  16.  UPSTREAM_GENE_ID: Entrez Gene ID for nearest upstream gene to rs number, if not within gene (NCBI)
  17.  DOWNSTREAM_GENE_ID: Entrez Gene ID for nearest downstream gene to rs number, if not within gene (NCBI)
  18.  SNP_GENE_IDS: Entrez Gene ID, if rs number within gene; multiple genes denotes overlapping transcripts (NCBI)
  19.  UPSTREAM_GENE_DISTANCE: distance in kb for nearest upstream gene to rs number, if not within gene (NCBI)
  20.  DOWNSTREAM_GENE_DISTANCE: distance in kb for nearest downstream gene to rs number, if not within gene (NCBI)
  21.  STRONGEST SNP-RISK ALLELE: SNP(s) most strongly associated with trait + risk allele (? for unknown risk allele). May also refer to a haplotype.
  22.  SNPS: Strongest SNP; if a haplotype is reported above, may include more than one rs number (multiple SNPs comprising the haplotype)
  23.  MERGED: denotes whether the SNP has been merged into a subsequent rs record (0 = no; 1 = yes; NCBI) SNP_ID_CURRENT: current rs number (will differ from strongest SNP when merged = 1)
  24.  CONTEXT: SNP functional class (NCBI)
  25.  INTERGENIC: denotes whether SNP is in intergenic region (0 = no; 1 = yes; NCBI)
  26.  RISK ALLELE FREQUENCY: Reported risk allele frequency associated with strongest SNP
  27.  P-VALUE: Reported p-value for strongest SNP risk allele (linked to dbGaP Association Browser)
  28.  PVALUE_MLOG: -log(p-value)
  29.  P-VALUE (TEXT): Information describing context of p-value (e.g. females, smokers). Note that p-values are rounded to 1 significant digit (for example, a published pvalue of 4.8 x 10-7 is rounded to 5 x 10-7).
  30.  OR or BETA: Reported odds ratio or beta-coefficient associated with strongest SNP risk allele
  31.  95% CI (TEXT): Reported 95% confidence interval associated with strongest SNP risk allele
  32.  PLATFORM (SNPS PASSING QC): Genotyping platform manufacturer used in Stage 1; also includes notation of pooled DNA study design or imputation of SNPs, where applicable
  33.  MAPPED_TRAIT: Mapped Experimental Factor Ontology trait for this study
  34.  MAPPED_TRAIT_URI: URI of the EFO trait

=cut

sub GWASCatalog {
  my ($disease, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GWASCatalog.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if ($items[7] == $disease) {
      foreach my $snp (split /,/, $items[21]) {
        push @res, (
          pvalue => $items[26],
          snp => $variation_adaptor->fetch_by_name($snp),
          disease => $disease,
          source => 'GWAS Catalog',
        )
      }
    }
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gwas_snps associated to disease $disease in GWAS Catalog\n";
  } 

  return \@res;
}

=head2 GRASP

  Returns all gwas_snps associated to a disease in GRASP
  Arg: Trait name or EFO ID (string)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  source => "GRASP",
                  study => string,
                }
              ]

=cut

=begin comment

  GRASP file format:
  1. NHLBIkey
  2. HUPfield
  3. LastCurationDate
  4. CreationDate
  5. SNPid(dbSNP134)
  6. chr(hg19)
  7. pos(hg19)
  8. PMID
  9. SNPid(in paper)
  10. LocationWithinPaper
  11. Pvalue
  12. Phenotype
  13. PaperPhenotypeDescription
  14. PaperPhenotypeCategories
  15. DatePub
  16. InNHGRIcat(as of 3/31/12)
  17. Journal
  18. Title
  19. IncludesMale/Female Only Analyses
  20. Exclusively Male/Female
  21. Initial Sample Description
  22. Replication Sample Description
  23. Platform [SNPs passing QC]
  24. GWASancestryDescription
  25. TotalSamples(discovery+replication)
  26. TotalDiscoverySamples
  27. European Discovery
  28. African Discovery
  29. East Asian Discovery
  30. Indian/South Asian Discovery
  31. Hispanic Discovery
  32. Native Discovery
  33. Micronesian Discovery
  34. Arab/ME Discovery
  35. Mixed Discovery
  36. Unspecified Discovery
  37. Filipino Discovery
  38. Indonesian Discovery
  39. Total replication samples
  40. European Replication
  41. African Replication
  42. East Asian Replication
  43. Indian/South Asian Replication
  44. Hispanic Replication
  45. Native Replication
  46. Micronesian Replication
  47. Arab/ME Replication
  48. Mixed Replication
  49. Unspecified Replication
  50. Filipino Replication
  51. Indonesian Replication
  52. InGene
  53. NearestGene
  54. InLincRNA
  55. InMiRNA
  56. InMiRNABS
  57. dbSNPfxn
  58. dbSNPMAF
  59. dbSNPalleles/het/se
  60. dbSNPvalidation
  61. dbSNPClinStatus
  62. ORegAnno
  63. ConservPredTFBS
  64. HumanEnhancer
  65. RNAedit
  66. PolyPhen2
  67. SIFT
  68. LS-SNP
  69. UniProt
  70. EqtlMethMetabStudy

=cut

sub GRASP {
  my ($disease, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GRASP.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if ($items[11] eq $disease) {
      push @res, {
        pvalue => $items[10],
        snp => $variation_adaptor->fetch_by_name($items[4]),
        disease => $disease,
        source => "GRASP",
        study => $items[7],
      }
    }
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gwas_snps associated to disease $disease in GRASP\n";
  } 

  return \@res;
}

=head2 PhewasCatalog

  Returns all gwas_snps associated to a disease in PhewasCatalog
  Arg: Trait name or EFO ID (string)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  source => "Phewas Catalog",
                  study => string,
                }
              ]

=cut

=begin comment

  Phewas Catalog format:
  1. chromosome
  2. snp
  3. phewas phenotype
  4. cases
  5. p-value
  6. odds-ratio
  7. gene_name
  8. phewas code
  9. gwas-associations
  
=cut

sub PhewasCatalog {
  my ($disease, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/PhewasCatalog.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if ($items[2] eq $disease) {
      push @res, {
          pvalue => $items[4],
          snp => $variation_adaptor->fetch_by_name($items[1]),
          disease => $disease,
          source => "Phewas Catalog",
          study => undef,
      }
    }
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gwas_snps associated to disease $disease in PhewasCatalog\n";
  } 

  return \@res;
}

=head2 GWAS_DB

  Returns all gwas_snps associated to a disease in GWAS_DB
  Arg: Trait name or EFO ID (string)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  source => "GWAS DB",
                  study => string,
                }
              ]

=cut

=begin comment

  GWAS DB data
  1. CHR
  2. POS
  3. SNPID
  4. REF
  5. ALT
  6. ORI_SNPID
  7. PMID
  8. P_VALUE
  9. P_VALUE_TEXT
  10. OR/BETA
  11. CI95_TEXT
  12. GWAS_INITIAL_SAMPLE_SIZE
  13. SUB_POPULATION
  14. SUPER_POPULATION
  15. GWAS_TRAIT
  16. HPO_ID
  17. HPO_TERM
  18. DO_ID
  19. DO_TERM
  20. MESH_ID
  21. MESH_TERM
  22. EFO_ID
  23. EFO_TERM
  24. DOLITE_TERM
  25. RISK_ALLELE
  26. PUBLICATION_TYPE
  27. AA
  28. GENE_SYMBOL
  29. TYPE
  30. REFGENE

=cut 

sub GWAS_DB {
  my ($disease, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GWAS_DB.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if ($items[14] eq $disease) {
      push @res, {
        pvalue => $items[7],
        snp => $variation_adaptor->fetch_by_name($items[2]),
        disease => $disease,
        source => "GWAS DB",
        study => $items[6],
      }
    }
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gwas_snps associated to disease $disease in GWAS_DB\n";
  } 

  return \@res;
}



=head2 gwas_snps_to_genes

  Associates Genes to gwas_snps of interest
  Args: 
  * [
      {
        snp => Bio::EnsEMBL::Variation::Variation,
        pvalue => scalar,
        evidence => [
          {
            pvalue => scalar,
            disease => string,
            source => string,
            study => string,
          }
        ]
      }
    ]
  * { $population_name => scalar (weight) }
  * { $tissue_name => $scalar (weight) }
  Returntype: [
                {
                  gene => Bio::EnsEMBL::Gene,
                  score => scalar, # aka pics.gene.score
                  gwas_snp_evidence => [
                    {
                      pvalue => scalar,
                      snp => Bio::EnsEMBL::Variation::Variation,
                      disease => string,
                      source => string,
                      study => string,
                    }
                  ]
                  association_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      score => scalar, # aka v2g.score
                      evidence => [
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ]
                    }
                  ]
                }
              ]

=cut

sub gwas_snps_to_genes {
  my ($gwas_snps, $populations, $tissue_weights) = @_;

  # Must set the tissue settings before separating out the gwas_snps
  if (! defined $tissue_weights) {
    $tissue_weights = gwas_snps_to_tissue_weights($gwas_snps);
  }

  my $clusters = cluster_gwas_snps(\@gwas_snp_locations, $populations);
  my @res = map { @$_ } map { cluster_to_genes($_, $tissue_weights, $population_weights) } @$clusters;

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to all gwas_snps\n";
  } 

  my @sorted_res = sort { $b->{score} <=> $a->{score} } @res;

  return [ $sorted_res[0] ];
}

=head2 gwas_snps_to_tissue_weights

  Associates list of tissues to list of gwas_snps
  Args: 
  * [
      {
        snp => Bio::EnsEMBL::Variation::Variation,
        pvalue => scalar,
        evidence => [
          {
            pvalue => scalar,
            disease => string,
            source => string,
            study => string,
          }
        ]
      }
    ]
  Returntype: [ string ]

=cut

sub gwas_snps_to_tissue_weights {
  my ($gwas_snps) = @_;
  return; # See FORGE??
}


=head2 cluster_gwas_snps

  Bundle together gwas_snps within LD threshold
  * [
      {
        snp => Bio::EnsEMBL::Variation::Variation,
        pvalue => scalar,
        evidence => [
          {
            pvalue => scalar,
            disease => string,
            source => string,
            study => string,
          }
        ],
      }
    ]
  * [ Bio::Ensembl::Variation::Population ]
  Returntype: [
                {
                  gwas_snps => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      pvalue => scalar,
                      evidence => [
                        {
                          pvalue => scalar,
                          disease => string,
                          source => string,
                          study => string,
                        }
                      ],
                    }
                  ],
                  ld_snps => [ Bio::EnsEMBL::Variation::VariationFeature ],
                }
              ]

=cut

sub cluster_gwas_snps {
  my ($gwas_snps, $populations) = @_;

  # Obtain locations of SNPs:
  my @gwas_snp_locations = map { @$_ } map(gwas_snp_to_locations, @$gwas_snps); 

  my @preclusters = map { gwas_snp_to_precluster($_, $populations) } @$gwas_snp_locations;

  return merge_clusters(@preclusters);
}

=head2 gwas_snp_to_locations

  Convert rsID association to location association:
  Args
  * {
      snp => Bio::EnsEMBL::Variation::Variation,
      pvalue => scalar,
      evidence => [
        {
          pvalue => scalar,
          disease => string,
          source => string,
          study => string,
        }
      ],
    }
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  pvalue => scalar,
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      source => string,
                      study => string,
                    }
                  ],
                }
              ]

=cut

sub gwas_snp_to_locations {
  my ($gwas_snp) = @_;

  my $vfs = $gwas_snp->{snp}->get_all_VariationFeatures;
  my $gwas_snp_locations = map { 
    {
      snp => $_,
      pvalue => $gwas_snp->{pvalue},
      evidence => $gwas_snp->{evidence},
    } 
  } @$vfs;

  return \@res;
}

=head2 gwas_snp_to_precluster

  Extract neighbourhood of GWAS snp
  Args:
  * {
      snp => Bio::EnsEMBL::Variation::VariationFeature,
      pvalue => scalar,
      evidence => [
        {
          pvalue => scalar,
          disease => string,
          source => string,
          study => string,
        }
      ],
    }
  * [ Bio::Ensembl::Variation::Population ]
  Returntype: {
                gwas_snps => [
                  {
                    snp => Bio::EnsEMBL::Variation::VariationFeature,
                    pvalue => scalar,
                    evidence => [
                      {
                        pvalue => scalar,
                        disease => string,
                        source => string,
                        study => string,
                      }
                    ],
                  }
                ],
                ld_snps => [ Bio::EnsEMBL::Variation::VariationFeature ],
              }
  
=cut

sub gwas_snp_to_precluster {
  my ($gwas_snp, $populations) = @_;

  # Get all LD values around SNP
  my $ld_feature_container = $hits->{snp}->get_all_LD_values;

  # TODO Configurable filter on population 
  my @filtered_ls_r_squared_values = grep { $_->{population} == $ld_feature_container->{'_default_population'} } @{$ld_feature_container->get_all_r_square_values};

  # Filter on R2 value
  my @filtered_ls_r_squared_values2 = grep { $_->{r2} > 0.5 } @$filtered_ld_r_squared_values;

  # Filter on distance
  my @filtered_ls_r_squared_values3 =  grep { 
    $_->{variation1}->seq_region_name eq $_->{variation2}->seq_region_name && abs($_->{variation1}->seq_region_start - $_->{variation2}->seq_region_start) < 500000 
  } @$ld_r_squared_values;

  # Reduce to list of interacting SNPs
  my @ld_snps = map { $_->{variation2} } grep { $_->{variation1}->name eq $gwas_snp->{snp}->name } @filtered_ld_r_squared_values3;
  push @ld_snps, map { $_->{variation1} } grep { $_->{variation2}->name eq $gwas_snp->{snp}->name } @filtered_ld_r_squared_values3;
  push @ld_snps, $gwas_snp->{snp};

  return {
    gwas_snps => [ $gwas_snp ],
    ld_snps => \@ld_snps,
  }
}

=head2 merge_clusters

  Bundle together clusters that share one LD snp
  * [
      {
        gwas_snps => [
          {
            snp => Bio::EnsEMBL::Variation::VariationFeature,
            pvalue => scalar,
            evidence => [
              {
                pvalue => scalar,
                disease => string,
                source => string,
                study => string,
              }
            ],
          }
        ],
        ld_snps => [ Bio::EnsEMBL::Variation::VariationFeature ],
      }
    ]
  Returntype: (same in input)

=cut

sub merge_clusters {
  my ($clusters) = @_;
 
  my %hash = ();
  foreach my $cluster (@$clusters) {
    foreach my $ld_snp (@{$cluster->{ld_snps}}) {
      if (exists $hash{$ld_snp->name}) {
        my $other_cluster = $hash{$ld_snp->name};

        # Merge data from current cluster into previous cluster
        push @{$other_cluster->{gwas_snps}}, @{$cluster->{gwas_snps}};
        $other_cluster->{ld_snps} = unique_variation_features( (@{$other_cluster->{ld_snps}}, @{$cluster->{ld_snps}}) );

        # Mark for deletion
        $cluster->{delete} = 1;

        # Exit from that cluster
        last;

      } else {
        $hash{$ld_snp->name} = $cluster;
      }
    }
  }
 
  my @new_clusters = grep { ! exists $_->{delete} } @$clusters;

  return \@new_clusters;
}

=head2 unique_variation_features

  Remove duplicates from list of VariationFeature objects
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  Return type: same as input 

=cut

sub unique_variation_features {
  my ($features) = @_;
  my %hash = map { $feature->name => $feature } @$features;
  my @unique_features = values %hash;
  return \@unique_features;
}

=head2 cluster_to_genes

  Associated Genes to a cluster of gwas_snps
  Args: 
  * [
      {
        gwas_snps => [
          {
            pvalue => scalar,
            snp => Bio::EnsEMBL::Variation::VariationFeature,
            disease => string,
            source => string,
            study => string,
          }
        ],
        ld_snps => [ Bio::EnsEMBL::Variation::VariationFeature ],
      }
    ]
  * { $tissue_name => scalar (weights) }
  * { $population_name => scalar (weight) }
  Returntype: [
                {
                  gene => Bio::EnsEMBL::Gene,
                  score => scalar, # aka pics.gene.score
                  gwas_snp_evidence => [
                    {
                      pvalue => scalar
                      snp => Bio::EnsEMBL::Variation::VariationFeature
                      disease => string 
                      source => string
                      study => string,
                    }
                  ],
                  association_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      score => scalar, # aka v2g.score
                      evidence => [ 
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ]
                    }
                  ]
                }
              ]

=cut

sub cluster_to_genes {
  my ($cluster, $tissues, $population_weights) = @_;

  # Obtain interaction data from LD snps
  my $hits = ld_snps_to_genes($cluster->{ld_snps}, $tissues);

  # Compute gene.score
  my $top_gwas_snp = get_top_gwas_snp($cluster);
  my $ld = get_lds_from_top_gwas($cluster->{ld_snps}, $top_gwas_snp->{snp});
  my %gene_scores = map { $_->{snp}->name => $_->{score} * $ld->{$_->{snp}->name} } @$hits;
  my $max_score = (sort { $b <=> $a } values %gene_scores)[0];
 
  # OMIM exception
  foreach my $hit (@$hits) {
    if (scalar @{gene_to_phenotypes($hit->{gene})}) {
      my $gene_score = $max_score;
    }
  }

  my $pics = PICS($ld, $top_gwas_snp->{pvalue});

  foreach my $hit (@$hits) {
    push @res, { 
      gene => $hit->{gene},
      score => total_score($pics->{$key}, $ld->{$hit->{snp}->name} * $hit->{score}), 
      gwas_snp_evidence => $cluster->{gwas_snp_evidence},
      association_evidence => $hit,
    }
    delete $hit->{gene};
  }

  if (our $DEBUG) {
    printf "Found ".scalar @values." genes associated around gwas_snp ".$top_gwas_snp->name."\n";
  } 


  # Pick the association with the highest score
  my $sorted_res = sort { $b->score <=> $a->score } @res;
  return [ $sorted_res[0] ];
}

=head2 get_top_gwas_snp

  Get GWAS SNP with highest P-value in cluster
  Args:
  * [
      {
        gwas_snps => [
          {
            pvalue => scalar,
            snp => Bio::EnsEMBL::Variation::VariationFeature,
            disease => string,
            source => string,
            study => string,
          }
        ],
        ld_snps => [ Bio::EnsEMBL::Variation::VariationFeature ],
      }
    ]
  Returntype: {
                pvalue => scalar,
                snp => Bio::EnsEMBL::Variation::VariationFeature,
                disease => string,
                source => string,
                study => string,
              }

=cut 

sub get_top_gwas_snp {
  my ($cluster) = @_;

  my @sorted_gwas_snps = sort { $a->{pvalue} <=> $b->{pvalue} } @{$cluster->{gwas_snps}};
  return $sorted_gwas_snps[0];
}

=head2 get_lds_from_top_gwas

  Compute LD between top GWAS hit and all LD snps in list.
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * Bio::EnsEMBL::Variation::VariationFeature
  Returntype: {
                $rsId => scalar (ld),
              }

=cut 

sub get_lds_from_top_gwas {
  my ($ld_snps, $gwas_snp, $population_weights) = @_;
  my $ld_container = $gwas_snp->get_all_r_square_values;
  my %hash = map { $_->name => $ld_container->get_r_square($_, $gwas_snp, $population_weights) } @$ld_snps;
  return \%hash; 
}

=head2 PICS

  PICS score presented in http://pubs.broadinstitute.org/pubs/finemapping/
  Args: 
  * { $rs_id => scalar } (LD)
  * scalar (pvalue)
  Returntype: { $rsID => scalar } (PICS score)

=cut

sub PICS {
  my ($ld, $pvalue) = @_;
 
  my $minus_log_pvalue = - log($pvalue) / log(10); 
  my %SD = ();
  my %Mean = ();
  my %prob = ();
  my $sum = 0;

  foreach $snp (keys %$ld) {
    if (defined $ld->{$snp}) {
      # Calculate the standard deviation of the association signal at the SNP 
      $SD->{$snp} = sqrt(1 - sqrt($ld->{$snp}) ** 6.4) * sqrt($minus_log_pvalue) / 2; 

      # calculate the expected mean of the association signal at the SNP 
      $Mean->{$snp} = $ld->{$snp} * $minus_log_pvalue; 
    } else {
      # Defaults for remote SNPs
      $SD->{$snp} = 0;
      $Mean->{$snp} = 1 + $minus_log_pvalue;
    }

    # Calculate the probability of each SNP
    if ($SD->{$snp} != 0) {
      $prob->{$snp} = 1 - pnorm($minus_log_pvalue, $Mean, $SD);
    } else {
      $prob->{$snp} = 1;
    }

    # Normalisation sum
    $sum += $prob->{$snp};

  }

  # Normalize the probabilies so that their sum is 1.
  my %res = map { $_ => $prob->{$_} / $sum } keys %prob;

  return \%res;
}

=head pnorm

  Normal distribution PDF
  Args:
  * scalar: variable
  * scalar: mean
  * scalar: standard deviation
  Return type: scalar (probability density)

=cut

sub pnorm {
  my ($x, $mu, $sd) = @_;
  return exp ( - (($x - $mu) / $sd) ** 2 / 2 ) / ($sd * 2.5);
}

=head2 total_score

  Computes a weird mean function from ld_snp PICs score and Gene/SNP association score
  Args: 
  * PICS: scalar
  * gene_score: scalar
  Returntype: scalar

=cut

sub total_score {
  my ($pics, $gene_score) = @_;
  
  my $A = $pics * ($pics ** (1/3));
  my $B = $gene_score * ($gene_score ** (1/3));
  return (($A + $B) / 2) ** 3;
}

=head2 ld_snps_to_genes

  Associates genes to LD linked SNP
  Args: 
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  Returntype: [ 
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  gene => Bio::EnsEMBL::Gene,
                  score => scalar, # aka v2g.score
                  evidence => [
                    {
                      tissue => tissue,
                      score => scalar,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

my @ld_snp_to_gene_functions = (\&GTEx, \&Fantom5, \&VEP);

sub ld_snps_to_genes {
  my ($ld_snps, $tissues) = @_;
  our $registry;
  our $SPECIES;
  my $variation_adaptor = $registry->get_adaptor($SPECIES, 'Core', 'Gene');
  my @evidence = map { @$_ } map { $_->($ld_snps, $tissues, $gene_adaptor) } @ld_snp_to_gene_functions;

  my %hash = ();
  foreach my $record (@evidence) {
    if ( !($record->{gene}->biotype eq "protein_coding")) { 
      next;
    }
    my $gene_name = $record->{gene}->name;
    $hash{$gene_name}{snp} = $ld_snp;
    $hash{$gene_name}{gene} = $gene;
    $hash{$gene_name}{score} += $record->{score};
    delete $record{gene};
    delete $record{snp};
    push @{$hash{$gene}{evidence}}, $record;
  }

  my @res = values %hash;
  return \@res;
}

=head2 GTEx

  Returns all genes associated to a set of SNPs in GTEx 
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  gene => Bio::EnsEMBL::Gene,
                  tissue => tissue,
                  score => scalar,
                  source => "GTEx"
                  study => string,
                }
              ]

=cut

sub GTEx {
  my ($ld_snp, $tissues) = @_
  
  # TODO

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$ld_snp->name.' in GTEx\n";
  } 

  return \@res;
}

=head2 VEP

  Returns all genes associated to a set of SNPs in VEP 
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  gene => Bio::EnsEMBL::Gene,
                  tissue => tissue,
                  score => scalar,
                  source => "VEP"
                  study => string,
                }
              ]

=cut

sub VEP {
  my ($ld_snp, $tissues) = @_

  # TODO

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$ld_snp->name.' in VEP\n";
  } 

  return \@res;
}

=head2 Fantom5

  Returns all genes associated to a set of SNPs in Fantom5 
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  gene => Bio::EnsEMBL::Gene,
                  tissue => tissue,
                  score => scalar,
                  source => "Fantom5"
                  study => string,
                }
              ]

=cut

=begin comment

  Format of Fantom5 correlation files:
  1.  chrom
  2.  chromStart
  3.  chromEnd
  4.  name ";" separated: 
    1. chrom:start-end
    2. Refseq
    3. HGNC
    4. R:$r2
    5. FDR:$fdr
  5.  score
  6.  strand
  7.  thickStart
  8.  thickEnd
  9.  itemRgb
  10.  blockCount
  11.  blockSizes
  12.  chromStarts


=cut

sub Fantom5 {
  my ($ld_snps, $tissues, $gene_adaptor) = @_;

  # Store rsID -> VariationFeature hash
  my %ld_hash = map { $_->name => $_ } @$ld_snps;

  # Dump LD SNP coords in temporary BED file
  my ($fh, $filename) = tempfile;
  foreach my $ld_snp (@$ld_snps) {
    print $fh join("\t", ())."\n";
  }    

  # Search for overlaps using bedtools
  my ($fh2, $filename2) = tempfile;
  our $DATABASES_DIR;
  system "bedtools intersect -wa -wb -a $DATABASES_DIR/Fantom5.txt $filename > $filename2 ";

  # Parse output: first 12 columns are from Fantom5 file, the next 4 are LD SNP coords
  while (<$fh2>) {
    chomp;
    my @items = split;
    my @association_data = split /;/, $items[3];
    my $new = {
      snp => $ld_hash{$items[15]},
      gene => $gene_adaptor->fetch_all_by_name($association_data[2])[0];
      tissue => undef,
      source => "Fantom5",
    }

    my @FDR = split(":", $association_data[4]); # Do we really need to re-do FDR calcs???
    if ($FDR[1] < .6) {
      $new->{score} = 2;
    } elsif ($FDR[1] < .85) {
      $new->{score} = 1;
    } else {
      $new->{score} = 0;
    }

    push @res, $new;
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gene associations in Fantom5\n";
  } 

  return \@res;
}

=head2 disease_to_MeSH

  Look up MeSH annotations for diseases
  Args:
  * [ string ] (disease names)
  Return type: [ string ] (annotations)

=cut

sub disease_to_MeSH {
  my ($diseases) = @_;
 
  my @res = ();
  my $server = "http://gene2mesh.ncibi.org/";

  foreach my $disease (@$diseases) {
    my $ext = "fetch?mesh="$disease;
    my $response = $http->get($server.$ext);
    die "Failed!\n" unless $response->{success};
    if(length $response->{content}) {
      my $hash = XMLin([$response->{content}]); # TODO Probably incorrect
      my $results = $hash->{NCIBI}{Gene2MeSH}{Response}{ResultSet};
      foreach my $result (@$results) {
        push @res, $result->{MeSH}{Descriptor}{Identifier};
      }
    }
  }

  return \@res; 
}

=head2 gene_to_phenotypes

  Look up phenotype annotations for gene
  Args:
  * Bio::EnsEMBL::Gene
  Return type: [ Bio::EnsEMBL::Variation::PhenotypeFeature ]

=cut

sub gene_to_phenotypes {
  my ($gene) = @_;

  if (!exists $phenotype_cache{$gene}) {
    $phenotype_cache{$gene->stable_id} = $phenotype_feature_adaptor->fetch_all_by_Gene($gene);
  }

  return $phenotype_cache{$gene->stable_id}; 
}

1;
