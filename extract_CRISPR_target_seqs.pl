#! usr/bin/perl -w
#use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2020.04.23	
# Modified:			 
############################################################

my $barcode_table = shift;
my $merged_fq = shift;
my $out = shift;

die  "
Usage:
perl   xxx.pl     barcode_table     merged.fq.gz     out_name

Note:
barcode_table.xls
Y1_FDDP202327169-1a_merged.fq.gz
Y1_FDDP202327169-1a_merged_targets_out

" if !defined $out;

my %reverse_comp = ( 
"A" => "T",
"T" => "A",
"C" => "G",
"G" => "C",
"N" => "N",
"a" => "t",
"t" => "a",
"c" => "g",
"g" => "c",
"n" => "n"
);

`date` =~ /\S+\s+(\S+)\s+(\d+)\s+(\d+):(\d+):\d+\s+\S+\s+(\d+)/;
my $date = $5.$1.$2.$3.$4;

my (%list, %head1, %head2, %check_duplicates, $rc_head1, $rc_head2);
my ($count, @t, $tail1, $tail2, $i, $j, $seq_head, $seq_tail);
my (%total_count, $rc_seq, $order, $percent);
############################################################


## Record barcode
open (IN, "< $barcode_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	if($t[0] !~ /^#/){
		$list{$t[0]} ++;
		$head1{$t[0]} = $t[1];
		$head2{$t[0]} = $t[2];
		
		$check_duplicates{$t[1]} ++;
		$check_duplicates{$t[2]} ++;
		$rc_head1 = base_reverse ($t[1]);
		$check_duplicates{$rc_head1}++;
		$rc_head2 = base_reverse ($t[2]);
		$check_duplicates{$rc_head2}++;
	}
}
close IN;


foreach $i (keys %check_duplicates){
	if($check_duplicates{$i} > 1){
		print "WARNING! Repeat Barcode!\t".$i."\n";
	}
}


## Count top3 abundant target seqs
if($merged_fq =~ /\.fq\.gz$/ || $merged_fq =~ /\.fastq\.gz$/){		## ".fq.gz"	or	".fastq.gz"
	open (IN, "gunzip -c $merged_fq |") or die $!;
}elsif($merged_fq =~ /\.fq$/ || $merged_fq =~ /\.fastq$/){			## ".fq"  	or	".fastq"
	open (IN, "< $merged_fq") or die $!;
}else{
	die "Check the fastq file name!\t".$merged_fq."\n"; 
}
while (<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if(/^@\w+:\d+:\w+/){
		$i = 0;
	}else{
		$i ++;
	}
	$t[$i] = $_;
	
	if($i == 1){
		$count ++;
        #print "Analysing merged reads: $count\r" if $count % 1000 == 0;
		
		foreach $j (sort keys %list){
			$tail1 = base_reverse ($head2{$j});
			$tail2 = base_reverse ($head1{$j});

			$seq_head = substr($t[$i], 0, 8);
			#print $seq_head."\t";
			$seq_tail = substr($t[$i], -8);
			#print $seq_tail."\n";
			
			if($seq_head eq $head1{$j} && $seq_tail eq $tail1){
				$total_count{$j}++;
				${"count_".$j}{$t[$i]} ++;
			}elsif($seq_head eq $head2{$j} && $seq_tail eq $tail2){
				$total_count{$j}++;
				$rc_seq = base_reverse ($t[$i]);
				${"count_".$j}{$rc_seq} ++;
			}
		}
	}
}
close IN;
print "\n";
undef $count;


## Output
open (OUT, "> $out"."_".$date.".xls") or die $!;
open (OUT2, "> $out"."_".$date."_top.fa") or die $!;
print OUT "#\tHead1\tHead2\tTotal_Number\tPercent\tTop1_seq\tTop1_Num\tPercent\tTop2_seq\tTop2_Num\tPercent\tTop3_seq\tTop3_Num\tPercent\n";
foreach $i (sort keys %list){
	
	if ($total_count{$i} && $total_count{$i} >= 100){
		print OUT $i."\t".$head1{$i}."\t".$head2{$i}."\t".$total_count{$i}."\t100%";
		
		undef $order;
		foreach $j (sort {${"count_".$i}{$b} <=> ${"count_".$i}{$a}} keys %{"count_".$i}){
			$order ++;
			
			undef $percent;
			$percent = int((${"count_".$i}{$j} / $total_count{$i})*10000) /100;
			
			if($order == 1){
				print OUT "\t".$j."\t".${"count_".$i}{$j}."\t".$percent."%";
				print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";
				$top1_percent = $percent;
			}elsif($order > 1){
				if($percent >= 5){
					print OUT "\t".$j."\t".${"count_".$i}{$j}."\t".$percent."%";
				}
				if( ($percent/$top1_percent) >= 0.2){
					print OUT2 ">".$i."_top".$order."_".$percent."%\n".$j."\n";	
				}
			}
		}
		print OUT "\n";
	}
}
close OUT;
close OUT2;
print "\n";


sub base_reverse{
	my @bre = reverse(split ('', $_[0]));
	my $baseout;
	foreach $i(@bre){
	     $baseout .= $reverse_comp{$i};
	}
	$baseout;
}
