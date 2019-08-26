#!/bin/perl -w
#Used the following websites for references
#https://www.perlmonks.org/?node_id=124970
#https://www.perlmonks.org/?node_id=1037314
($bamfile)=@ARGV;
#open(SAM_HEAD, "samtools view -H ../190BD_bowtie.bam |");
open(SAM_HEAD, "samtools view -H $bamfile |");
while (<SAM_HEAD>){
    chomp;
    my @tmp = split /\t/, $_;
    #print "@tmp\n";
    $seq = $tmp[0];
    $contig = $tmp[1];
    $length = $tmp[2];
    if ($contig =~ /SN:(\S+)/){
        $contig =$1;
       #  print "$contig\n";
    }

    if ($length =~ /LN:(\d+)/){
        $length =$1;
#        print "$length\n";

    $l_c{$contig} = $length;
#    print "$length\t$l_c{$contig}\t$contig\n";
        print "$tmp[0]\t$tmp[1]\t$tmp[2]\n";
    }
}
#

##For very strict selection (-f3 included:
#open(SAM_FILE, "samtools view -F 1548 -f 3 -q 1 $bamfile | ");
open(SAM_FILE, "samtools view -F 1548 -q 1 $bamfile | ");

while (<SAM_FILE>){
    chomp;
    my @tmp2 = split /\t/, $_;
    $QNAME=$tmp2[0];
    $FLAG=$tmp2[1];
    $trinity= $tmp2[2];
    $position=$tmp2[3];
    $MAPQ=$tmp2[4];
    $cigar=$tmp2[5];
    $RNEXT= $tmp2[6];
    $PNEXT= $tmp2[7];
    $TLEN= $tmp2[8];
    $SEQ=$tmp2[9];
    $QUAL=$tmp2[10];
 #   print "\~\~\-\-\-\-\-\-";
 #   print "$trinity\t$cigar\tPosition\: $position\n";
 #   print "\~\~\-\-\-\-\-\-\n";
##-------------------

    @parts = split /(\d+\D)/, $cigar;
    @segs = grep { $_ ne '' } @parts;

    $sumseqs = digitify(@segs);
 #   print "Total of whole cigar\: $sumseqs\n";

##For the matching ones
    @mcig = grep { /\d+M/ } @segs;
    $sum_mcig = digitify(@mcig);
#   print "Total only matching bases in cigar: $sum_mcig\n";

    @m_d_cig = grep { /[DMN]/ } @segs;
    $sum_m_d_cig= digitify(@m_d_cig);
#    print "MI_cigars \: $sum_m_d_cig\n";

#---------
#####
    $subject_margin = ($l_c{$trinity} +1 - ($position + $sum_m_d_cig));
    $left = ($segs[0] =~ /(\d+)S/) ? $1 : 0;
    $right = ($segs[-1] =~ /(\d+)S/) ? $1 : 0;

    ## modify for position and length of the template
    $left = ($left < $position) ? $left : $position;
    $right = ($right < $subject_margin) ? $right : $subject_margin;

    ## sanity check.. left and right should never be negative...
    if($left < 0 || $right < 0){
        die "impossible margins: $left / $right\n";
    }

    $score = $sum_mcig / ($sum_mcig + $left + $right);
    if ($score > 0.95){
#       print "$trinity\t$position\t$cigar\n";
        print "$QNAME\t$FLAG\t$trinity\t$position\t$MAPQ\t$cigar\t$RNEXT\t$PNEXT\t$TLEN\t$SEQ\t$QUAL\n";
#       print join "\t", @tmp2, "\n";
#    print "The score is : $score\n";
    }
}
close(SAM_HEAD);
close(SAM_FILE);

sub add {
    my $sum = 0;
    foreach ( @_ ) {
        $sum += $_;
    }
    return $sum;
}

sub digitify {
    @dn = @_;
    @seq_dn = grep { s/(\d+)\D/$1/g } @dn;
    $sum_dn= add(@seq_dn);
    return $sum_dn;
}


