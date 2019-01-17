#!/bin/bash

#wd=$(pwd)

# #Code for how parameter matrix is generated:
# R
# fineness=.005
# efseq <- seq(1.05,1.3,by=fineness)
# eureffectsize <- c("matching","null","weaker","stronger","euronly")
# afmaf <- c(0.1,0.2, 0.3,0.4)
# eurmaf <- c("smaller","larger","same")
# admixturegroup <- c(0.5,0.8)
# N <- c(12000)
# prevalence <- c(0.1,0.2)
# parameter_matrix  <- expand.grid(efseq,eureffectsize,afmaf,eurmaf,admixturegroup,N,prevalence, KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)
# names(parameter_matrix) <- c("afeffectsize","eureffectsize","afmaf","eurmaf","admixturegroup","N","prevalence")
# write.csv(parameter_matrix,file="parameter_matrix.csv",row.names=F)


#Do multiple runs. Takes about 5 runs of to make really stable results

# for runno in {1..5}
# do
# qsub simulate_power.sh -lwalltime=3:00:00 -d $wd -e errandout/ -o errandout/ -F "-n 7 -r $runno"
# done

while getopts n:r: option
do
  case "${option}"
    in
      n) nodeuse=${OPTARG};;
      r) run=${OPTARG};;
    esac
done

#want to split N commands into K jobs, run L jobs at a time

ncommands=$(wc -l parameter_matrix.csv   | awk '{print $1}')
nodeuse=7

nodesize=$nodeuse
totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))


for i in $(seq 1 $nodesize $totjobs)
 do
  #Run jobs K..K+node num
  jstart=$i
  jstop=$(($i + $nodesize - 1))
  min=$([ $ncommands -le $jstop ] && echo "$ncommands" || echo "$jstop")
  jstop=$min
  
  for j in $(seq  $jstart 1 $jstop) 
  do
  linestart=$((($j-1)*$nodeuse +1))
  linestop=$(($j*$nodeuse))

  Rscript /mnt/sdb/genetics/elizabeth_power_simulation/lanc_simulation_v5_loopsims2.txt parameter_matrix.csv $linestart $linestop $run &
  
  echo $linestart $linestop
  
  done
  wait
  echo "batch $i"
done
