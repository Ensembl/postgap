use strict;
use warnings;
 
use HTTP::Tiny;
use XML::Simple;
use Data::Dumper;
local $Data::Dumper::Terse = 1;
local $Data::Dumper::Indent = 1;
 
my $http = HTTP::Tiny->new();

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

my $server = "http://gene2mesh.ncibi.org";
my $ext = "/fetch?genesymbol=csf1r&taxid=9606";
my $response = $http->get($server.$ext);
die "Failed!\n" unless $response->{success};

if(length $response->{content}) {
  my $hash = XMLin($response->{content});
  print Dumper $hash;
  print "\n";
}
