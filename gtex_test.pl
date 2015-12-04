use strict;
use warnings;
 
use HTTP::Tiny;
use JSON;
use Data::Dumper;
 
=begin comment

  Example return object:
  [
    {
      'value' => '0.804108648395327',
      'snp' => 'rs142557973'
    },
  ]

=cut 
my $http = HTTP::Tiny->new();
 
my $server = "http://193.62.54.30:5555";
my $ext = "/eqtl/id/homo_sapiens/ENSG00000227232?content-type=application/json;statistic=p-value;tissue=Whole_Blood"; 
my $response = $http->get($server.$ext);
die "Failed!\n" unless $response->{success};

if(length $response->{content}) {
  my $hash = decode_json($response->{content});
  local $Data::Dumper::Terse = 1;
  local $Data::Dumper::Indent = 1;
  print Dumper $hash;
  print "\n";
}
