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
use warnings;
 
use HTTP::Tiny;
use JSON;
 
my $http = HTTP::Tiny->new();

main();

sub main() {
  my $phenotype_column = $ARGV[0];
  while (my $line = <STDIN>) {
    chomp $line;
    my @columns = split(/\t/, $line);
    my $efo = efo_suggest($columns[$phenotype_column]);
    if (defined $efo) {
      print STDOUT $line."\t".$efo."\n";
    } else {
      print STDERR "Could not find EFO for $columns[$phenotype_column]\n";
      print STDOUT $line."\tN/A\n";
    }
  }
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

sub efo_suggest {
  my ($term) = shift;
  my $server = 'http://www.ebi.ac.uk/spot/zooma/v2/api';
  my $url_term = $term;
  $url_term =~ s/ /%20/g;
  my $ext = "/summaries/search?query=$url_term";
  my $response = $http->get($server.$ext);
   
  die "Failed!\n" unless $response->{success};

  my $hit;
   
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

1;
