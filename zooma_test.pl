use strict;
use warnings;
 
use HTTP::Tiny;
use JSON;
use Data::Dumper;
 
my $http = HTTP::Tiny->new();

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

  Example annotation result
  {
    'annotationSummaryTypeID' => 'summary_from_disease_to_EFO_0000311',
    'semanticTags' => [
      'http://www.ebi.ac.uk/efo/EFO_0000311'
    ],
    'uri' => 'http://rdf.ebi.ac.uk/resource/zooma/annotation_summary/3804279AF8462F3A01EAEE2589C94781F248F9D7',
    'annotatedPropertyValue' => 'cancer',
    'annotatedPropertyType' => 'disease',
    'annotationSummaryTypeName' => 'disease; EFO_0000311',
    'annotatedPropertyUri' => 'http://rdf.ebi.ac.uk/resource/zooma/8963601691028A7A9AA7BC747B04E586',
    'annotationSourceURIs' => [
      'http://www.ebi.ac.uk/gxa',
      'http://www.ebi.ac.uk/efo/efo.owl'
    ],
    'annotationURIs' => [
      'http://rdf.ebi.ac.uk/resource/zooma/gxa/EDC7A5E2985C660FBFDF64A0E0D1397A',
    ],
    'id' => '3804279AF8462F3A01EAEE2589C94781F248F9D7',
    'quality' => '86.11212'
  }

=cut
 
my $server = 'http://www.ebi.ac.uk/spot/zooma/v2/api';
my $ext = '/summaries/search?query=cancer';
my $response = $http->get($server.$ext);
 
die "Failed!\n" unless $response->{success};

my $hit;
 
if(length $response->{content}) {
  my $hash = decode_json($response->{content});
  my $result = $hash->{result};
  my @hits = sort { $b->{score} <=> $a->{score} } @$result;
  local $Data::Dumper::Terse = 1;
  local $Data::Dumper::Indent = 1;
  $hit = $hits[0];
  print Dumper $hit;
  print "\n";
}

$ext = "/summaries/" . $hit->{mid};
$response = $http->get($server.$ext);

die "Failed!\n" unless $response->{success};
 
if(length $response->{content}) {
  my $hash = decode_json($response->{content});
  local $Data::Dumper::Terse = 1;
  local $Data::Dumper::Indent = 1;
  print Dumper $hash;
  print "\n";
}
