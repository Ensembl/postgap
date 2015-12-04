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
	  'TaxonomyID' => '9606',
	  'MeSHDescriptor' => 'diabetes',
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
		'content' => '0',
		'type' => 'p-value'
	      },
	      'ChiSquare' => '10267.441327824',
	      'Fover' => '41.772344933059',
	      'MeSH' => {
		'Qualifier' => {
		  'Name' => 'genetics'
		},
		'Descriptor' => {
		  'TreeNumber' => [
		    'C18.452.394.750.124',
		    'C19.246.267',
		    'C20.111.327'
		  ],
		  'Identifier' => 'D003922',
		  'Name' => 'Diabetes Mellitus, Type 1',
		  'UMLSID' => {}
		}
	      },
	      'DocumentSet' => {
		'type' => 'pubmed',
		'PMID' => [
		  '2187469',
		]
	      },
	      'Gene' => {
		'Taxonomy' => {
		  'Identifier' => '9606'
		},
		'Identifier' => '3119',
		'type' => 'Entrez',
		'Description' => 'major histocompatibility complex, class II, DQ beta 1',
		'Symbol' => 'HLA-DQB1'
	      }
	    },
	  ],
	  'sort' => 'FisherExact',
	  'count' => '800',
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

my $server = "http://gene2mesh.ncibi.org/";
my $ext = "fetch?mesh=diabetes&taxid=9606";
my $response = $http->get($server.$ext);
die "Failed!\n" unless $response->{success};

if(length $response->{content}) {
  my $hash = XMLin($response->{content});
  print Dumper $hash;
  print "\n";
}
