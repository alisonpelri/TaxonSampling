use strict;
use warnings;
use Bio::SeqIO;

my $file = $ARGV[0];

my $in = new Bio::SeqIO(-format => 'genbank',
                           -file => $file);
my $seq;

my $name = $file;

my $i = 1;

$file =~ s/.gb$//;

open (OUT1, ">$file.fasta");
open (OUT2, ">$file\_species_list.txt");
open (OUT3, ">$file\_species_name_2_sequence_id.txt");

while($seq = $in->next_seq()){
  my $sequence = $seq->seq();
  my $id = $seq->display_id();
  print OUT1 (">$id\n$sequence\n");
  my $start;
  my $end;
  my $organism;
    for my $feat_object ($seq->get_all_SeqFeatures) {
    print "primary tag: ", $feat_object->primary_tag, "\n";
    for my $tag ($feat_object->get_all_tags) {
      for my $value ($feat_object->get_tag_values($tag)) {
      print ("$tag\t$value\n");
      }
    }
    if ($feat_object->primary_tag eq "source") {
#       print "primary tag: ", $feat_object->primary_tag, "\n";
      for my $tag ($feat_object->get_all_tags) {
#          print ("$tag\n");
        if ($tag eq "organism") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $organism = $value;
            print OUT2 ("$organism\n");
            print OUT3 ("$organism\t$id\n");
          }
        }
      }
    }
  }
}

close OUT1;
close OUT2;
close OUT3;
