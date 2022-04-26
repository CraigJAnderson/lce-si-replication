#!/bin/bash
##pt3_mut_count.sh
##script to annotate all the windows phasable through lesion segregation from each nodule
##use as: ./pt3_mut_count.sh 0.15 0.10 0.125 0 1 2 watson_v_crick ${SAMPLE_LIST} input.anno.bed

POS1=${1}
POS2=${2}
POS3=${3}
MODE=${4} ##if 1, get random selection of mutations from cohort for bootstrap. If 0, run on observed mutations
COMBINE=${5} ##if 1 then combine watson and crick, else if 0, then report for watson and crick, separately. 
BS=${6} ##bootstrap number
FOLDER=${7}
NOD=${8}
BIN=${9}
PTH=${10}
DIR=${11}
GENOME=${12}

##check the folders for output exist.
if [ ! -d ${FOLDER} ]
 then mkdir ${FOLDER}
fi

if [ ! -d ${FOLDER}/bs${BS} ]
 then mkdir ${FOLDER}/bs${BS}
fi

if [[ $COMBINE -eq 1 ]] 
then
 if [[ $MODE -eq 1 ]]
 then
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.mut
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.W.mut
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.C.tnt_count
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.W.tnt_count
   if [ ! -f ${FOLDER}/bs${BS}/${NOD}.${BS}.bed ]
    then NUM=$(wc -l ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed | cut -d " " -f 1)
    python ${BIN}/random_select2.py ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed ${FOLDER}/bs${BS}/${NOD}.${BS}.bed ${NUM}
   fi
  for x in C W
   do X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${FOLDER}/bs${BS}/${NOD}.${BS}.bed | bedtools intersect -a ${PTH}/${DIR}/variants/${NOD}.filtered.nodMat -b stdin | Z=$x awk '{ if ($5 == ENVIRON["Z"]) print $1"\t"$2"\t"$3"\t"substr($4,1,3)}' | sort -gk1,1 -k2,2g | grep -v X | grep -v N | cut -f 4 | cat - ${BIN}/tnt_list.txt | sort | uniq -c | awk '{print $2"\t"$1-1}' | datamash transpose | tail -n 1 | paste <(echo ${POS3}) - >> ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
   VAR=$(X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${FOLDER}/bs${BS}/${NOD}.${BS}.bed | bedtools getfasta -s -tab -name+ -fi ${GENOME} -bed stdin | perl ${BIN}/streamTriNucCounter.pl | tail -n 1 | cut -f2-65) 
   echo $POS3 $VAR | sed 's/ /\t/g' >> ${FOLDER}/bs${BS}/${NOD}.rfd.${x}.tnt_count
  done
  python ${BIN}/combine_tnt_counts.py ${FOLDER}/bs${BS}/${NOD}.rfd.C.tnt_count ${FOLDER}/bs${BS}/${NOD}.rfd.W.tnt_count > ${FOLDER}/bs${BS}/${NOD}.rfd.tnt_count
  python ${BIN}/combine_mut_counts.py ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.mut ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.W.mut > ${FOLDER}/bs${BS}/${NOD}.rfd.mut
  python ${BIN}/matrix_mu_rate.weighted_norm.py ${FOLDER}/bs${BS}/${NOD}.rfd.tnt_count ${FOLDER}/bs${BS}/${NOD}.rfd.mut ${BIN}/rfd_tntcount.txt ${POS3} >> ${FOLDER}/bs${BS}/bs${BS}.${NOD}.mu
  rm ${FOLDER}/bs${BS}/${NOD}.rfd.mut
  rm ${FOLDER}/bs${BS}/${NOD}.rfd.tnt_count
 else
  for x in C W
   do cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
   X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed | bedtools intersect -a ${PTH}/${DIR}/variants/${NOD}.filtered.nodMat -b stdin | Z=$x awk '{ if ($5 == ENVIRON["Z"]) print $1"\t"$2"\t"$3"\t"substr($4,1,3)}' | sort -gk1,1 -k2,2g | grep -v X | grep -v N | cut -f 4 | cat - ${BIN}/tnt_list.txt | sort | uniq -c | awk '{print $2"\t"$1-1}' | datamash transpose | tail -n 1 | paste <(echo ${POS3}) - >> ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
  done
 python ${BIN}/combine_mut_counts.py ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.mut ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.W.mut > ${FOLDER}/bs${BS}/${NOD}.rfd.mut
 python ${BIN}/matrix_mu_rate.weighted_norm.py ${FOLDER}/${NOD}.rfd.COM.tnt_count ${FOLDER}/bs${BS}/${NOD}.rfd.mut ${BIN}/rfd_tntcount.txt ${POS3} >> ${FOLDER}/bs${BS}/bs${BS}.${NOD}.mu
 fi
else
 if [[ $MODE -eq 1 ]]
 then 
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.mut
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.W.mut
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.C.tnt_count
   cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.W.tnt_count
   if [ ! -f ${FOLDER}/bs${BS}/${NOD}.${BS}.bed ]
    then NUM=$(wc -l ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed | cut -d " " -f 1)
    python ${BIN}/random_select2.py ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed ${FOLDER}/bs${BS}/${NOD}.${BS}.bed ${NUM}
   fi
  for x in C W
   do X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${FOLDER}/bs${BS}/${NOD}.${BS}.bed | bedtools intersect -a ${PTH}/${DIR}/variants/${NOD}.filtered.nodMat -b stdin | Z=$x awk '{ if ($5 == ENVIRON["Z"]) print $1"\t"$2"\t"$3"\t"substr($4,1,3)}' | sort -gk1,1 -k2,2g | grep -v X | grep -v N | cut -f 4 | cat - ${BIN}/tnt_list.txt | sort | uniq -c | awk '{print $2"\t"$1-1}' | datamash transpose | tail -n 1 | paste <(echo ${POS3}) - >> ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
   VAR=$(X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${FOLDER}/bs${BS}/${NOD}.${BS}.bed | bedtools getfasta -s -tab -name+ -fi ${GENOME} -bed stdin | perl ${BIN}/streamTriNucCounter.pl | tail -n 1 | cut -f2-65)
   echo $POS3 $VAR | sed 's/ /\t/g' >> ${FOLDER}/bs${BS}/${NOD}.rfd.${x}.tnt_count
  done
 else
  for x in C W
   do cat ${BIN}/tnt_header.txt > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
   X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' ${PTH}/${DIR}/sample_specific_rfd_windows/${NOD}.input.anno.${FOLDER}.bed | bedtools intersect -a ${PTH}/${DIR}/variants/${NOD}.filtered.nodMat -b stdin | Z=$x awk '{ if ($5 == ENVIRON["Z"]) print $1"\t"$2"\t"$3"\t"substr($4,1,3)}' | sort -gk1,1 -k2,2g | grep -v X | grep -v N | cut -f 4 | cat - ${BIN}/tnt_list.txt | sort | uniq -c | awk '{print $2"\t"$1-1}' | datamash transpose | tail -n 1 | paste <(echo ${POS3}) - >> ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.${x}.mut
  done
 fi
 python ${BIN}/revcom_mut_counts.py ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.mut > ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.revcom.mut
 MU1=$(python ${BIN}/matrix_mu_rate.weighted_norm.py ${FOLDER}/${NOD}.rfd.C.tnt_count ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.revcom.mut ${BIN}/rfd_tntcount.txt ${POS3})
 MU2=$(python ${BIN}/matrix_mu_rate.weighted_norm.py ${FOLDER}/${NOD}.rfd.W.tnt_count ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.W.mut ${BIN}/rfd_tntcount.txt ${POS3})
 echo "C" $MU1 $'\nW' $MU2 >>  ${FOLDER}/bs${BS}/bs${BS}.${NOD}.mu
 rm ${FOLDER}/bs${BS}/${NOD}.rfd.${BS}.C.revcom.mut
fi
