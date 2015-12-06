use strict;
use Storable;
use Data::Dumper;
local $Data::Dumper::Terse = 1;
local $Data::Dumper::Indent = 1;

my $BIN_SIZE = 10000;
my $MAX_DISTANCE = 5e5;

use Bio::EnsEMBL::Registry;
our $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
  -host => 'ens-livemirror',
  -user => 'ensro',
); 
my $gene_adaptor = $registry->get_adaptor('Homo sapiens', 'Core', 'Gene');

my %distribution = ();

while (<STDIN>) {
  process_line($_, \%distribution);
}

my %res = {
  FDR => map { $_ => $distribution{99} * .99 / $distribution{$_} } keys %distribution;
  MAX_BIN => $MAX_DISTANCE,
  BIN_SIZE => $BIN_SIZE,
};

print Dumper $res;

store_fd $res, *STDOUT;

sub process_line {
  my ($line, $distribution) = shift;
  chomp $line;
  my @columns = split /\t/, $line;

  my $chr = $columns[0];
  if ($chr =~ /^#/) {
    return;
  }
  $chr =~ s/chr//;
  my $start = $columns[1];
  my $end = $columns[2];
  my $genes = $gene_adaptor->fetch_all_by_external_name($columns[3]);
  if (scalar @$genes == 0) {
    return;
  } 

  my $pos = ($start + $end) / 2;
  my $gene = $genes->[0];

  if (! ($gene->seq_region_name eq $chr)) {
    return; 
  }

  my $tss;
  if ($gene->strand) {
    $tss = $gene->seq_region_start;
  } else {
    $tss = $gene->seq_region_end;
  }

  my $distance = abs($pos - $tss);

  if ($distance > $MAX_DISTANCE) {
    return;
  }

  $distribution->{int($distance / $BIN_SIZE)}++;
}
