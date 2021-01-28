#step1:extract read pairs for which the distance between the mapped locations of the two reads was between 1kb and 20kb(from sorted SAM file)

open(FI,"${file}.mem.bwa.sort.sam");
open(SE,">${file}_large_tlen.sam");
while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
  
   if(($all[6] eq "=") && ($all[2] ne "mitochondrion")&& ($all[4] >= 20)&& (abs($all[8]) >= 1000)&& (abs($all[8]) <= 20000)){     
	    print SE $line,"\n";
	    
    }
}

#step2: Remove the redundant readout from SAM file

open(FI,"${file}_large_tlen.sam");
 while($line=<FI>){
           chomp $line;
           @all=split/\t/,$line;
             $hash{$all[0]}++;
      }
 close(FI);
 
open(FI,"${file}_large_tlen.sam");
open(SE,">${file}_large_tlen_unique.sam");
 while($line=<FI>){
   chomp $line;
   @all=split/\t/,$line;
    if(($hash{$all[0]}==1) &&($line =~ m/SA:Z/)) {
	    print SE $line,"\n";
	}
    elsif(($hash{$all[0]} >= 2) &&($all[8] > 0)) {
	    print SE $line,"\n";
	}
}
 
#step3: cluster read pairs supporting the same event.  
 
open(FI,"${file}_large_tlen_unique.sam");
open(SE,">${file}_paired_cluster_1");
while($line=<FI>){
    chomp $line;
    @all=split/\t/,$line;
	if($all[8]>0){
       $left = int($all[3]/1000)*1000;
	   $right = int($all[7]/1000)*1000 + 1000;
	}
	if($all[8]<0){
	   $left = int($all[7]/1000)*1000;
	   $right = int($all[3]/1000)*1000 + 1000;
	}
	print SE $all[2],"\t",$left,"\t",$right,"\t",$line,"\n";
}
close(FI);
close(SE);

open(FI,"${file}_paired_cluster_1");
open(SE,">${file}_paired_cluster_2");
%hash=();
while($line=<FI>){
    chomp $line;
    @all=split/\t/,$line;
	$name = $all[0]."_".$all[1]."_".$all[2];
	$hash{$name}++;
}
for $a(keys %hash){
   @all=split/_/,$a;
   print SE "sample_name","\t",$all[0],"\t",$all[1],"\t",$all[2],"\t",$hash{$a},"\n";
}
close(FI);
close(SE);

#step4: Define excision event
 open(FI,"${file}_paired_cluster_2");
 open(SE,">${file}_paired_cluster_TE");
 while($line=<FI>){
          chomp $line;
          @all=split/\t/,$line;
	      $i=0;
            open(IN,"sc_reference_sequence.fa.out.gff");
            while($line2=<IN>){
              chomp $line2;
              @all2=split/\t/,$line2;
			  if (($all[1] eq $all2[0])&& ($all2[3] - $all[2] >= 0)&&($all2[4] - $all[3] <= 0)){
			    $i++;
			  }
			}		
		    close(IN);
			if($i>3){
				 print SE $line,"\n";
            }
        }
 


