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

our @database_functions = (\&GWASCatalog, \&GRASP, \&GWAS_DB, \&Phewas_Catalog);
our @ld_snp_to_gene_functions = (\&GTEx, \&Fantom5, \&VEP);
our @snp_regulatory_functions = (\&GERP, \&Regulome);

our %phenotype_cache = ();
our $VEP_impact_to_score = {
  HIGH => 4,
  MEDIUM => 3,
  LOW => 2,
  MODIFIER => 1,
}

main();

=begin comment

  Development TODO list:

  A. Datasets to integrate:
    -Cis-regulatory annotations:
      PCHIC (STOPGAP Scoring: Single cell line: +1, multiple cell lines: +2)

    -Epigenetic activity:
      PhyloP (STOPGAP Scoring: FPR 0-0.6: +2, 0.6-0.85: +1,0.85-1: +0)
      DHS
      Fantom5

  B. Code improvements:
    Pathways analysis (Downstream)
    Take into account population composition in LD calcs

  C. Model improvements:
    Replace PICS with Bayesian model
    Fine mapping of summary data
    Tissue selection

=cut 

=head2 main

  Reads commandline parameters, prints corresponding associated genes with evidence info

=cut

sub main {
  my $options = get_options();
  print encode_json(diseases_to_genes($options->{diseases}, $options->{efos}, $options->{populations}, $options->{tissues}));
}

=head2 get_options

  Reads commandline parameters
  Returntype: {
                diseases => [ string ],
                populations => [ Bio::EnsEMBL::Variation::Population ],
                tissues => [ string ],
              }

=cut

sub get_options {
  my %options = ();

  GetOptions(\%options, "help=s", "efo|o=s@", "disease|d=s@", "population|p=s@", "tissue|t=s@", "databases|D=s", "host=s", "user=s", "port=s", "species=s", "debug|g");

  $options{diseases} = $options{disease};
  $options{populations} = map { $population_adaptor->fetch_by_name($_) } @$options{population};
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

  our $DEBUG ||= defined $options->{$debug});

  if (!defined $options->{efo}) {
    my @efos = grep { defined $_ } map { efo_suggest($_) } @{$options->{diseases}};
    $options->{efo} = \@efos;
  }

  # Expand list of EFOs to children, concatenate, remove duplicates
  my $options->{efos} = keys map { $_ => $_ } map { @$_ } map(efo_children, @{$options->{efo}});

  return \%options;
}

=head2 efo_suggest
  
  Find most appropriate EFO term for arbitrary string
  Arg:
  * string
  Returntype: string (EFO ID)

=cut

=begin comment

  Example search result:
  {
    status: "/api/status/ok",
    result: [
      {
	'mid' => '3804279AF8462F3A01EAEE2589C94781F248F9D7',
	'notable' => {
	  'name' => 'disease; EFO_0000311',
	  'id' => 'summary_from_disease_to_EFO_0000311'
	},
	'name' => 'cancer',
	'score' => '86.11212'
      }
    ]  
  }

=cut

=begin comment

  See info in preprocessing./EFO_suggest.pl

=cut

sub efo_suggest {
  my ($term) = shift;
  my $server = 'http://www.ebi.ac.uk/spot/zooma/v2/api';
  my $url_term = $term;
  $url_term =~ s/ /%20/g;
  my $ext = "/summaries/search?query=$url_term";
  my $response = $http->get($server.$ext);
   
  die "Failed!\n" unless $response->{success};

  if(length $response->{content}) {
    my $hash = decode_json($response->{content});
    my $result = $hash->{result};
    my @hits = sort { $b->{score} <=> $a->{score} } grep { $_->{notable}{name} =~ /EFO_[0-9]+/ } @$result;
    my $selection = $hits[0]{notable}{name};
    if (!defined $selection) {
      return;
    }
    $selection =~ s/.*(EFO_[0-9]*)$/$1/;
    return $selection;
  } else {
    return;
  }
}

=head2 efo_children

  Return list of children EFO IDs
  Arg:
  * string (EFO ID)
  Returntype: [ string ] (EFI IDs)

=cut

=begin comment

  OWL Output format:
  {
    "_links" : {
      "first" : {
	"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0001071/descendants?page=0&size=10"
      },
      "prev" : {
	"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0001071/descendants?page=0&size=10"
      },
      "self" : {
	"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0001071/descendants?page=1&size=10"
      },
      "last" : {
	"href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_0001071/descendants?page=1&size=10"
      }
    },
    "_embedded" : {
      "terms" : [ {
	"iri" : "http://www.ebi.ac.uk/efo/EFO_1000333",
	"label" : "Lung Inflammatory Myofibroblastic Tumor",
	"description" : [ "An intermediate fibroblastic neoplasm arising from the lung. It is characterized by the presence of spindle-shaped fibroblasts and myofibroblasts, and a chronic inflammatory infiltrate composed of eosinophils, lymphocytes and plasma cells." ],
	"annotation" : {
	  "NCI_Thesaurus_definition_citation" : [ "NCIt:C39740" ]
	},
	"synonyms" : null,
	"ontology_name" : "efo",
	"ontology_prefix" : "EFO",
	"ontology_iri" : "http://www.ebi.ac.uk/efo",
	"is_obsolete" : false,
	"is_defining_ontology" : true,
	"has_children" : false,
	"is_root" : false,
	"short_form" : "EFO_1000333",
	"obo_id" : "EFO:1000333",
	"_links" : {
	  "self" : {
	    "href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_1000333"
	  },
	  "parents" : {
	    "href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_1000333/parents"
	  },
	  "ancestors" : {
	    "href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_1000333/ancestors"
	  },
	  "jstree" : {
	    "href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_1000333/jstree"
	  },
	  "graph" : {
	    "href" : "http://www.ebi.ac.uk/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252FEFO_1000333/graph"
	  }
	}
      } ]
    },
    "page" : {
      "size" : 10,
      "totalElements" : 20,
      "totalPages" : 2,
      "number" : 1
    }
  }

=cut

sub efo_children {
  my ($efo) = @_;

  my ($term) = shift;
  my $server = 'http://www.ebi.ac.uk';
  my $page = 1;
  my @res = ($efo);

  while (1) {
    my $ext = "/ols/beta/api/ontologies/efo/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252F$efo/descendants?page$page&size=10";
    my $response = $http->get($server.$ext);
     
    die "Failed!\n" unless $response->{success};

    if(length $response->{content}) {
      my $hash = decode_json($response->{content});
      my $results = $hash->{_embedded}{terms};
      push @res, map { $_->{short_form} } @$results;
      if ($page > $hash->{page}{totalPages}) {
        last;
      }
    } 

    $page++;
  }
  
  return \@res;
}

=head2 diseases_to_genes

  Associates genes from a list of diseases
  Args: 
  * [ string ] (trait descriptions - free strings)
  * [ string ] (trait EFO identifiers)
  * [ string ] (population names)
  * [ string ] (tissue names)
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
                          efo => string,
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
                      cis_regulatory_evidence => [
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ],

                      regulatory_evidence => [
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
  my ($diseases, $efos, $populations, $tissues) = @_;
  my $res = gwas_snps_to_genes(diseases_to_gwas_snps($diseases, $efos), $populations, $tissues);

  foreach my $candidate (@$res) {
    $candidate->{MeSH} = gene_to_MeSH($candidate->{gene});
    $candidate->{gene_phenotype_association} = gene_to_phenotypes($candidate->{gene});
  }

  return $res;
}

=head2 diseases_to_gwas_snps

  Associates gwas_snps from a list of diseases
  Args: 
  * [ string ] (trait descriptions - free strings )
  * [ string ] (trait EFO identifiers)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::Variation,
                  pvalue => scalar, (min of pvalues - could be improved)
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      efo => string,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

sub diseases_to_gwas_snps {
  my ($diseases, $efos) = @_;

  # Extract crude data
  my $gwas_snps = scan_disease_databases($diseases, $efos);

  # Filter by p-value
  @filtered_gwas_snps = grep { $_->{pvalue} < $PVALUE_CUTOFF } @$gwas_snps;

  return \@filtered_gwas_snps;
}

=head2 scan_disease_databases 

  Associates gwas_snps from a list of diseases
  Args: 
  * [ string ] (trait descriptions)
  * [ string ] (trait EFO identifiers)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::Variation,
                  pvalue => scalar, (min of pvalues - could be improved)
                  evidence => [
                    {
                      pvalue => scalar,
                      disease => string,
                      efo => string,
                      source => string,
                      study => string,
                    }
                  ]
                }
              ]

=cut

sub scan_disease_databases {
  my ($diseases, $efos) = @_;

  our $registry;
  our $SPECIES;
  my $variation_adaptor = $registry->get_adaptor($SPECIES, 'Variation', 'Variation');

  my @hits = map { @$_ } map { $_->($diseases, $efos, $variation_adaptor) } @database_functions;

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

=head2 GWASCatalog

  Returns all gwas_snps associated to a disease in GWAS Catalog
  Args:
  * [ string ] (trait descriptions)
  * [ string ] (trait EFO identifiers)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  efo => string,
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
  my ($diseases, $efos, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GWASCatalog.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if (grep(/^$items[7]$/, @$diseases) || grep( /http:\/\/www.ebi.ac.uk\/efo\/$items[33]$/, @$efos)) {
      foreach my $snp (split /,/, $items[21]) {
        push @res, (
          pvalue => $items[26],
          snp => $variation_adaptor->fetch_by_name($snp),
          disease => $items[7],
	  efo => $items[34],
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
  Args:
  * [ string ] (trait descriptions)
  * [ string ] (trait EFO identifiers)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  efo => string,
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
  71. EFO string 

=cut

sub GRASP {
  my ($diseases, $efos, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GRASP.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if ($items[11] eq $disease) {
    if (grep(/^$items[11]$/, @$diseases) || grep(/^$items[70]$/, @$efos)) {
      push @res, {
        pvalue => $items[10],
        snp => $variation_adaptor->fetch_by_name($items[4]),
        disease => $items[11],
	efo => $items[70],
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
  Args:
  * [ string ] (trait descriptions)
  * [ string ] (trait EFO identifiers)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  efo => string,
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
  10. [Inserte] EFO identifier (or N/A)
  
=cut

sub PhewasCatalog {
  my ($diseases, $efos, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/PhewasCatalog.txt";
  my @res = ();

  while (<$fh>) {
    chomp;
    my @items = split;
    if (grep(/^$items[2]$/, @$diseases) || grep(/^$items[9]$/, @$efos)) {
      push @res, {
          pvalue => $items[4],
          snp => $variation_adaptor->fetch_by_name($items[1]),
          disease => $items[3],
	  snp => $items[9],
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
  Args:
  * [ string ] (trait descriptions)
  * [ string ] (trait EFO identifiers)
  Returntype: [ 
                {
                  pvalue => scalar,
                  snp => Bio::EnsEMBL::Variation::Variation,
                  disease => string,
                  efo => string,
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
  my ($diseases, $efos, $variation_adaptor) = @_
  our $DATABASES_DIR;
  open my $fh, "<", $DATABASES_DIR."/GWAS_DB.txt";

  # Note the format of the EFO strings is modified in this file format, so we need to change the queries
  my @efos2 = map { (my $v = $_) =~ s/_/ID:/; $v} @$efos;

  my @res = ();
  while (<$fh>) {
    chomp;
    my @items = split;
    if (grep(/^$items[14]$/, @$diseases) || grep(/^$items[21]$/, @efos2)) {
      push @res, {
        pvalue => $items[7],
        snp => $variation_adaptor->fetch_by_name($items[2]),
        disease => $items[14],
	efo => $items[21],
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
            efo => string,
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
                      efo => string,
                      source => string,
                      study => string,
                    }
                  ]
                  association_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      score => scalar, # aka v2g.score
                      cis_regulatory_evidence => [
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ]
                      regulatory_evidence => [
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
  my @res = map { @$_ } map { cluster_to_genes($_, $tissue_weights, $populations) } @$clusters;

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
            efo => string,
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
            efo => string,
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
                          efo => string,
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
          efo => string,
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
                      efo => string,
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
          efo => string,
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
                        efo => string,
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

  # Configurable filter on population 
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
                efo => string,
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
            efo => string,
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
                      efo => string,
                      source => string
                      study => string,
                    }
                  ],
                  association_evidence => [
                    {
                      snp => Bio::EnsEMBL::Variation::VariationFeature,
                      score => scalar, # aka v2g.score
                      cis_regulatory_evidence => [ 
                        {
                          tissue => tissue,
                          score => scalar,
                          source => string,
                          study => string,
                        }
                      ]
                      regulatory_evidence => [ 
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
  my ($cluster, $tissues, $populations) = @_;

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
            efo => string,
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
                efo => string,
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
  my ($ld_snps, $gwas_snp, $populations) = @_;
  my $ld_container = $gwas_snp->get_all_r_square_values;
  my %hash = map { $_->name => $ld_container->get_r_square($_, $gwas_snp, $populations) } @$ld_snps;
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
                  cis_regulatory_evidence => [
                    {
                      tissue => tissue,
                      score => scalar,
                      source => string,
                      study => string,
                    }
                  ]
                  regulatory_evidence => [
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

sub ld_snps_to_genes {
  my ($ld_snps, $tissues) = @_;
  our $registry;
  our $SPECIES;
  my $gene_adaptor = $registry->get_adaptor($SPECIES, 'Core', 'Gene');
  my @evidence = map { @$_ } map { $_->($ld_snps, $tissues, $gene_adaptor) } @ld_snp_to_gene_functions;

  # Group by (gene,snp) pair:
  foreach my $record (@evidence) {
    if ( !($record->{gene}->biotype eq "protein_coding")) { 
      next;
    }

    my $gene_name = $record->{gene}->name;
    my $rsID = $hit->{snp}->name;

    # Basic info
    $hash{$gene_name}{$rsID}{snp} = $hit->{snp};
    $hash{$gene_name}{$rsID}{gene} = $hit->{gene};

    # Integrate cis-interaction data
    $hash{$gene_name}{$rsID}{score} += $record->{score};
    push @{$hash{$gene}{$rsID}{cis_regulatory_evidence}}, $record;
    delete $record{gene};
    delete $record{snp};
  }

  # Collapse these hashrefs of hashrefs into a single list
  my @res = map { values %$_ } values %hash;

  # Integrate SNP specific info
  my %selected_snp_hash = map { $_->{snp}->name => $_->{snp} } @res;
  my @selected_snps = values %selected_snp_hash; 
  my %regulatory_evidence_hash = regulatory_annotation_at_snps(\@selected_snps, $tissues);
  foreach my $record (@res) {
    $record->{regulatory_evidence} = $regulatory_evidence_hash->{$snp->name}; 
    map { $record->{score} += $_->{score} } @{$record->{regulatory_evidence}};
  }

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
  my ($ld_snps, $tissues, $gene_adaptor) = @_
  
  # Find all genes with 1Mb
  my @starts = sort { $a <=> $b } map { $_->seq_region_start } @$ld_snps;
  my @ends = sort { $b <=> $a } map { $_->seq_region_end } @$ld_snps;
  my $slice = $ld_snps->slice;
  $slice->seq_region_start($starts->[0] - 1000000);
  $slice->seq_region_end($ends->[0] + 1000000);
  my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

  my %snp_hash = map { $_->name => $_ } @$ld_snps;
  my @res = map { @$_ } map { GTEx_gene($gene, $tissues, \%snp_hash) };

  return \@res;
}

=head2 GTEx_gene

  Returns all SNPs associated to a gene in GTEx 
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  * { $rsID => Bio::EnsEMBL::Variation::VariationFeature }
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

sub GTEx_gene {
  my ($gene, $tissues, $snp_hash) = @_
  
  my @res = map { @$_ } map { GTEx_gene_tissue($gene, $tissue) } @$tissues;

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP '.$gene->external_name.' in GTEx\n";
  } 

  return \@res;
}

=head2 GTEx_gene_tissue

  Returns all SNPs associated to a gene in GTEx in a given tissue
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  * { $rsID => Bio::EnsEMBL::Variation::VariationFeature }
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

=begin comment

  Example return object:
  [
    {
      'value' => '0.804108648395327',
      'snp' => 'rs142557973'
    },
  ]

=cut 

sub GTEx_gene_tissue {
  my ($gene, $tissue, $snp_hash) = @_
  
  my $server = "http://193.62.54.30:5555";
  my $ext = "/eqtl/id/homo_sapiens/$gene->stable_id?content-type=application/json;statistic=p-value;tissue=$tissue"; 
  my $response = $http->get($server.$ext);
  die "Failed!\n" unless $response->{success};

  if(length $response->{content}) {
    my $list = decode_json($response->{content});
    my @res = ();
    foreach my $hit (@$list) {
      if (exists $snp_hash->{$hit->{snp}}) {
        push @res, {
	  snp => $snp_hash->{$hit->{snp}}, 
	  gene => $gene,
	  tissue => $tissue,
	  score => $hit->{value},
	  source => "GTEx",
	  study => undef,
        };
      }
    }
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP $gene->external_name in tissue $tissue in GTEx\n";
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
  my ($ld_snps, $tissues) = @_

  # TODO Batch requests to the REST server
  my @res = map { @$_ } map { VEP_snp($snp, $tissues, $gene_adaptor) } @$ld_snps;

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated in cluster in VEP\n";
  } 

  return \@res;
}

=head2 VEP_snp

  Returns all genes associated to a SNP in VEP 
  Args:
  * Bio::EnsEMBL::Variation::VariationFeature
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

sub VEP_snp {
  my ($ld_snp, $tissues, $gene_adaptor) = @_
  our $VEP_impact_to_score;

  my $server = "http://rest.ensembl.org";
  my $ext = "/vep/human/id/$ld_snp->name?content-type=application/json";
  my $response = $http->get($server.$ext);
  die "Failed!\n" unless $response->{success};
  
  my $list = decode_json($response->{content});
  my @res = map {
    {
      snp => $ld_snp,
      gene => $gene_adaptor->fetch_by_stable_id($_->{gene_id}),
      tissue => undef,
      score => $VEP_impact_to_score->{$_->{impact}},
      source => "VEP"
      study => undef,
    }
  } @{$list->{transcript_consequences}};

  if (our $DEBUG) {
    printf "Found ".scalar @res." genes associated to SNP $ld_snp->name in cluster in VEP\n";
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

  my $fdr_model = retrieve("$DATABASES_DIR/Fantom5.fdrs");

  # Parse output: first 12 columns are from Fantom5 file, the next 4 are LD SNP coords
  while (<$fh2>) {
    chomp;
    my @items = split;
    my @association_data = split /;/, $items[3];

    my $gene = $gene_adaptor->fetch_all_by_name($association_data[2])[0];
    my $snp = $ld_hash{$items[15]};
    my $score = STOPGAP_FDR($snp, $gene, $fdr_model);

    if ($score == 0) {
      next;
    }

    my $new = {
      snp => $snp,
      gene => $gene,
      tissue => undef,
      source => "Fantom5",
      score => $score,
    }

    push @res, $new;
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gene associations in Fantom5\n";
  } 

  return \@res;
}

=head2 DHS

  Returns all genes associated to a set of SNPs in DHS 
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ] (tissues)
  Returntype: [
                {
                  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  gene => Bio::EnsEMBL::Gene,
                  tissue => tissue,
                  score => scalar,
                  source => "DHS"
                  study => string,
                }
              ]

=cut

=begin comment

  Format of DHS correlation files:
  1. chrom
  2. chromStart
  3. chromEnd
  4. HGNC
  5. Correlation

=cut

sub DHS {
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
  system "bedtools intersect -wa -wb -a $DATABASES_DIR/DHS.txt $filename > $filename2 ";

  my $fdr_model = retrieve("$DATABASES_DIR/DHS.fdrs");

  # Parse output: first 5 columns are from DHS file, the next 4 are LD SNP coords
  while (<$fh2>) {
    chomp;
    my @items = split;

    my $gene = $gene_adaptor->fetch_all_by_name($items[3])[0];
    my $snp = $ld_hash{$items[8]};
    my $score = STOPGAP_FDR($snp, $gene, $fdr_model);

    if ($score == 0) {
      next;
    }

    my $new = {
      snp => $snp,
      gene => $gene,
      tissue => undef,
      source => "DHS",
      score => $score,
    }

    push @res, $new;
  }

  if (our $DEBUG) {
    printf "Found ".scalar @res." gene associations in DHS\n";
  } 

  return \@res;
}

=head2 STOPGAP_FDR

  Special function for cis-regulatory interactions 
  Args:
  * Bio::EnsEMBL::Variation::VariationFeature
  * Bio::EnsEMBL::Gene
  * 

  Returntype: scalar

=cut

sub STOPGAP_FDR {
  my ($snp, $gene, $fdr_model) = @_;

  if (! ($gene->seq_region_name eq $snp->seq_region_name)) {
    return 0;
  }

  my $pos = ($start + $end) / 2;

  my $tss;
  if ($gene->strand) {
    $tss = $gene->seq_region_start;
  } else {
    $tss = $gene->seq_region_end;
  }

  my $distance = abs($pos - $tss);

  if ($distance > $fdr_model->MAX_DISTANCE) {
    return 0;
  }

  my $FDR = $fdr_model->FDR[int($distance / $FDR->BIN_WIDTH)];

  if (!defined $FDR) {
    return 0;
  }

  my $score;
  if ($FDR < .6) {
    return 2;
  } elsif ($FDR < .85) {
    return 1;
  } else {
    return 0;
  }
}

=head2 regulatory_annotation_at_snps

  Extract regulatory evidence linked to SNPs and stores them in a hash
  * [ Bio::EnsEMBL::Variation::VariationFeature ]
  * [ string ]
  Returntype: {
                $dbID => [
		  {
		    tissue => tissue,
		    score => scalar,
		    source => string,
		    study => string,
		  }
		]
              } 

=cut 

sub regulatory_annotation_at_snps {
  my ($snp, $tissues) = @_;

  my @res = map { @$_ } map { $_->($snps, $tissues) } @snp_regulatory_functions;

  # Group by SNP
  my %hash = ();
  foreach my $hit (@res) {
    push @{$hash{$hit->{snp}->name}}, $hit;
    delete $hit->{snp};
  }

  return \%hash;
}

=head2 GERP

  Extract GERP score at position
  Args:
  * [ Bio::EnsEMBL::Variation::VariationFeature ] 
  Returntype: [
                {
		  snp => Bio::EnsEMBL::Variation::VariationFeature,
                  tissue => undef,
                  score => scalar,
                  source => 'GERP',
                  study => undef,
                }
              ]

=cut 

sub GERP {
  my ($snps, $tissues) = @_;

  my %hash = map { 
    {
      snp => $snp,
      tissue => undef,
      score => GERP_at_snp($_),
      source => 'GERP',
      study => undef,
    }
  } @$snps;

  return \%hash;
}

=head2 GERP_at_snp

  Extract GERP score at position
  Args:
  * Bio::EnsEMBL::Variation::VariationFeature
  Returntype: scalar 
=cut 

sub GERP {
  my ($snp) = @_;

  my $server = "http://rest.ensembl.org";
  my $ext = "/vep/human/id/$ld_snp->name?content-type=application/json;Conservation=1";
  my $response = $http->get($server.$ext);
  die "Failed!\n" unless $response->{success};
  
  my $hash = decode_json($response->{content});

  return \%hash;
}

=head2 Regulome

  Extract Regulome score at sns of interest 
  Args:
  * Bio::EnsEMBL::Variation::VariationFeature
  Returntype: [
                {
                  tissue => undef,
                  score => scalar,
                  source => 'Regulome',
                  study => undef,
                }
              ]

=cut 

=begin comment

  Regulome file format:
  1. chrom
  2. start
  3. end
  4. category

=cut

sub Regulome {
  my ($snps, $tissues) = @_;

  # Store rsID -> VariationFeature hash
  my %hash = map { $_->name => $_ } @$snps;

  # Dump LD SNP coords in temporary BED file
  my ($fh, $filename) = tempfile;
  foreach my $snp (@$snps) {
    print $fh join("\t", ($snp->seq_region_name, $snp->seq_region_start, $snp->seq_region_end, $snp->name))."\n";
  }    

  # Search for overlaps using bedtools
  my ($fh2, $filename2) = tempfile;
  our $DATABASES_DIR;
  system "bedtools intersect -wa -wb -a $DATABASES_DIR/Regulome.txt $filename > $filename2 ";

  my $fdr_model = retrieve("$DATABASES_DIR/Regulome.fdrs");

  # Parse output: first 4 columns are from Regulome file, the next 4 are LD SNP coords
  my @res = ();
  while (<$fh2>) {
    chomp;
    my @items = split;

    my $snp = $hash{$items[7]};
    my $category = $hash{$items[3]};

    my $score;
    if ($category =~ /^1/ || $category =~ /^2/) {
      $score = 2;
    } else {
      $score = 1;
    }

    push @res, {
      snp => $snp,
      tissue => undef,
      source => "DHS",
      score => $score,
    }
  }

  return \@res;
}

=head2 gene_to_MeSH

  Look up MeSH annotations for gene
  Args:
  * [ string ] (gene names)
  Return type: [ string ] (annotations)

=cut

=begin comment

 Example MeSH output:
  {
    'Gene2MeSH' => {
      'Request' => {
	'ParameterSet' => {
	  'Tool' => 'none',
	  'GeneSymbol' => 'csf1r',
	  'TaxonomyID' => '9606',
	  'Limit' => '1000',
	  'Email' => 'anonymous'
	},
	'type' => 'fetch'
      },
      'Response' => {
	'ResultSet' => {
	  'Result' => [
	    {
	      'FisherExact' => {
		'content' => '1.8531319238671E-230',
		'type' => 'p-value'
	      },
	      'ChiSquare' => '112213.6506462',
	      'Fover' => '1498.1813411401',
	      'MeSH' => {
		'Qualifier' => {
		  'Name' => 'metabolism'
		},
		'Descriptor' => {
		  'TreeNumber' => [
		    'D08.811.913.696.620.682.725.400.500',
		    'D12.776.543.750.060.492',
		    'D12.776.543.750.705.852.150.150',
		    'D12.776.543.750.750.400.200.200',
		    'D12.776.624.664.700.800',
		    'D23.050.301.264.035.597',
		    'D23.101.100.110.597'
		  ],
		  'Identifier' => 'D016186',
		  'Name' => 'Receptor, Macrophage Colony-Stimulating Factor',
		  'UMLSID' => {}
		}
	      },
	      'DocumentSet' => {
		'type' => 'pubmed',
		'PMID' => [
		]
	      },
	      'Gene' => {
		'Taxonomy' => {
		  'Identifier' => '9606'
		},
		'Identifier' => '1436',
		'type' => 'Entrez',
		'Description' => 'colony stimulating factor 1 receptor',
		'Symbol' => 'CSF1R'
	      }
	    },
	  ],
	  'sort' => 'FisherExact',
	  'count' => '94',
	  'order' => 'ascending'
	},
	'Copyright' => {
	  'Details' => 'http://nlp.ncibi.org/Copyright.txt',
	  'Year' => '2009',
	  'Statement' => 'Copyright 2009 by the Regents of the University of Michigan'
	},
	'Support' => {
	  'Details' => 'http://www.ncibi.org',
	  'GrantNumber' => 'U54 DA021519',
	  'Statement' => 'Supported by the National Institutes of Health as part of the NIH\\\'s National Center for Integrative Biomedical Informatics (NCIBI)'
	}
      }
    }
  }

=cut 

sub gene_to_MeSH {
  my ($gene) = @_;

  my %NCBI_Taxon_ID = (
    'Human' => 9606;
  );
 
  my $server = "http://gene2mesh.ncibi.org";
  my $ext = "/fetch?genesymbol=$gene->external_name&taxid=$NCBI_Taxon_ID{$SPECIES}";
  my $response = $http->get($server.$ext);
  die "Failed!\n" unless $response->{success};

  if(length $response->{content}) {
    my $hash = XMLin($response->{content});
    my $hits = $hash->{Response}{ResultSet}{Result};
    my @res = map { $_->{MeSH}{Descriptor}{Name} } @$hits;
    return \@res;
  } else {
    return;
  }
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
