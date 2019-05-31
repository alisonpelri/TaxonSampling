use strict;
use warnings;

my %speciesName2sequenceIDs;

my %taxonID2speciesName;

open(IN, $ARGV[0]); #file containing species IDs in first column and sequence IDs, in second, one per line

while (my $line = <IN>) {
  chomp $line;
  my @aux = split(/\t/, $line);
  my $tmp = $aux[1];
  $tmp =~ s/\s\+/_/g;
  $speciesName2sequenceIDs{$aux[0]} = $tmp;
}

close IN;

open(IN, $ARGV[1]);

while (my $line = <IN>) {
  next if ($line !~ /^Query:/);
  $line =~ s/^Query://;
  chomp $line;
  my @aux = split(/\t/, $line);
  next if ($#aux < 1);
  $taxonID2speciesName{$aux[1]} = $aux[0];
}

close IN;

foreach my $key (keys %taxonID2speciesName) {
  my $species = $taxonID2speciesName{$key};
  if (defined $speciesName2sequenceIDs{$species}) {
    print "$key\t$speciesName2sequenceIDs{$species}\n";
    #print "Here!\n";
  }
}
