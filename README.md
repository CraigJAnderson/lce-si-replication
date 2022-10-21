# lce-si-replication

set up the working environment
<pre>
conda env create -f lce-si.yml
conda activate lce
</pre>

run the following to combine strands
<pre>
snakemake --snakefile WG_snakefile -pr -j 2
</pre>

run the following to run for watson and crick strands, separately
<pre>
snakemake --snakefile WVC_snakefile -pr -j 2
</pre>

If you're working with nodMat format from the other LCE work, use the following to reformat the list of variants
<pre>
while read line ; do awk -F"," '{print $1"\t"$2"\t"$2+1"\t"$10"\tNA"}' /exports/igmm/eddie/taylor-lab/lceCORE/data/nodules/${line}.nodMat | sed 1d | awk '$4 ~ /[A-Z]/' > /exports/igmm/eddie/taylor-lab/craig/bin/lce-si-replication/variants/${line}.filtered.nodMat ; done < src/sample_list.txt
</pre>

Below is some rough code for stripping out the results in unix and plotting in R
<pre>
##now sumarise results, needs to have input for CI calc, as in you tell it how many samples it needs to head and tail off as 5%
cd out

##comparing watson and crick in genic regions
for z in C W ;
 do for y in $(seq 1 1 ${NUMBER_OF_RFD_BINS}) 
  do WIN=$(sed -n ${y}p src/rfd_list.txt | cut -f3)
  for x in {1..520}
   do VAR=$(grep $z ./genic/bs${x}/bs${x}.${WIN}.mu | awk '{ total += $5 } END { print total/NR }') 
   echo $WIN $z $VAR >> wc_all_mean_mu.txt
  done
  AVE=$(grep $z genic/bs0/bs0.${WIN}.mu | awk '{ total += $5 } END { print total/NR }')
  MAX=$(grep $z wc_all_mean_mu.txt | grep -w "^$WIN" | sort -gk2,2 | tail -n 13 | head -n 1)
  MIN=$(grep $z wc_all_mean_mu.txt | grep -w "^$WIN" | sort -gk2,2 | head -n 13 | tail -n 1)
  echo $z ${WIN} ${AVE} ${MIN} ${MAX} >> wc_all_mean_mu_95CI.txt
 done
 rm wc_all_mean_mu.txt
done

##comparing watson and crick in nongenic regions
for z in C W ;
 do for y in $(seq 1 1 ${NUMBER_OF_RFD_BINS}) 
  do WIN=$(sed -n ${y}p src/rfd_list.txt | cut -f3)
  for x in {1..520}
   do VAR=$(grep $z ./nongenic/bs${x}/bs${x}.${WIN}.mu | awk '{ total += $5 } END { print total/NR }') 
   echo $WIN $z $VAR >> wc_all_mean_mu.txt
  done
  AVE=$(grep $z nongenic/bs0/bs0.${WIN}.mu | awk '{ total += $5 } END { print total/NR }')
  MAX=$(grep $z wc_all_mean_mu.txt | grep -w "^$WIN" | sort -gk2,2 | tail -n 13 | head -n 1)
  MIN=$(grep $z wc_all_mean_mu.txt | grep -w "^$WIN" | sort -gk2,2 | head -n 13 | tail -n 1)
  echo $z ${WIN} ${AVE} ${MIN} ${MAX} >> wc_all_mean_mu_95CI.txt
 done
 rm wc_all_mean_mu.txt
done

###combined analysis
##combined genic
for y in $(seq 1 1 ${NUMBER_OF_RFD_BINS}) 
 do WIN=$(sed -n ${y}p src/rfd_list.txt | cut -f3)
 for x in {1..520}
  do VAR=$(awk '{ total += $3 } END { print total/NR }' genic/bs${x}/bs${x}.${WIN}.mu) 
  echo $WIN $VAR >> combined_genic_all_mean_mu.txt
 done
 AVE=$(awk '{ total += $3 } END { print total/NR }' genic/bs0/bs0.${WIN}.mu)
 MAX=$(grep -w "^$WIN" combined_genic_all_mean_mu.txt | sort -gk2,2 | tail -n 13 | head -n 1 | cut -f2)
 MIN=$(grep -w "^$WIN" combined_genic_all_mean_mu.txt | sort -gk2,2 | head -n 13 | tail -n 1 | cut -f2)
 echo genic ${WIN} ${AVE} ${MIN} ${MAX} >> combined_genic_all_mean_mu_95CI.txt
done
rm combined_genic_all_mean_mu.txt

##combined nongenic
for y in $(seq 1 1 ${NUMBER_OF_RFD_BINS})
 do WIN=$(sed -n ${y}p src/rfd_list.txt | cut -f3)
 for x in {1..520}
  do VAR=$(awk '{ total += $3 } END { print total/NR }' nongenic/bs${x}/bs${x}.${WIN}.mu) 
  echo $WIN $VAR >> combined_nongenic_all_mean_mu.txt
 done
 AVE=$(awk '{ total += $3 } END { print total/NR }' nongenic/bs0/bs0.${WIN}.mu)
 MAX=$(grep -w "^$WIN" combined_nongenic_all_mean_mu.txt | sort -gk2,2 | tail -n 13 | head -n 1 | cut -f2)
 MIN=$(grep -w "^$WIN" combined_nongenic_all_mean_mu.txt | sort -gk2,2 | head -n 13 | tail -n 1 | cut -f2)
 echo nongenic ${WIN} ${AVE} ${MIN} ${MAX} >> combined_nongenic_all_mean_mu_95CI.txt
done
rm combined_nongenic_all_mean_mu.txt

###plotting
##now plot in R- this example is for watson v crick in the genic and nongenic regions
c <- read.table("wc_all_mean_mu_95CI.txt",header=T,sep="\t")
c$mu <- c$mu*1000000
c$u <- c$u*1000000
c$l <- c$l*1000000
par(oma=c(2,1,0,4),mar=c(3,3,3,7),xpd=TRUE)
plot(mu ~ rfd,data=c,pch=19,ylim=c(6.8,16),cex=0,xlim=c(-0.7,0.7), ylab="mu/Mb", xlab="RFD")
genic <- subset(c, type == "Genic")
nongenic <- subset(c, type == "Nongenic")
watson <- subset(c, type == "Watson")
crick <- subset(c, type == "Crick")

polygon(c(genic$rfd,rev(genic$rfd)),c((genic$u),rev((genic$l))),col=rgb(0, 0, 0,0.15),border=NA)
polygon(c(nongenic$rfd,rev(nongenic$rfd)),c((nongenic$u),rev((nongenic$l))),col=rgb(0, 0, 0,0.15),border=NA)
polygon(c(watson$rfd,rev(watson$rfd)),c((watson$u),rev((watson$l))),col=rgb(0, 0, 0,0.15),border=NA)
polygon(c(crick$rfd,rev(crick$rfd)),c((crick$u),rev((crick$l))),col=rgb(0, 0, 0, 0.15),border=NA)

lines(genic$rfd,genic$mu,col="black",lwd=3)
lines(nongenic$rfd,nongenic$mu,col="black",lwd=3)
lines(watson$rfd,watson$mu,col="black",lwd=3)
lines(crick$rfd,crick$mu,col="black",lwd=3)

#axis(side=1,at=c(-0.3,-0.15,0,0.15,0.3),labels=c("-0.3","-0.15","0","0.15","0.3"),lwd=1,col="black",) ##fewer ticks
#axis(side=2,line=-1,at=c(7,9,11,13,15),labels=c("7","9","11","13","15"),lwd=1,col="black")
legend(0.35,12, legend=c("Nongenic","Genic","Watson","Crick"),col=c("black","black","black","black"), lty=c(3,1),box.lty=0)
title(ylab="mu/Mb", line=1.5, cex.lab=1)
title(xlab="RFD", line=2.2, cex.lab=1)
</pre>
