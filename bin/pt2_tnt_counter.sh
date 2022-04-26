#!/bin/bash
##pt2_tnt_counter.sh
##script to count trinucleotide occurrence throughout all the windows phasable through lesion segregation from each nodule

POS1=$1
POS2=$2
POS3=$3
COMBINE=$4 ##if 0, then report for watson and crick, separately. 
REGIONS=$5 ##analysis subset
NOD=$6
BIN=$7
GENOME=$8

##check the folders for output exist.
if [ ! -d ${REGIONS} ]
 then mkdir ${REGIONS}
fi

for x in C W
 do cat ${BIN}/tnt_header2.txt > ${REGIONS}/tmp_${NOD}.rfd.${x}.tnt_count
 VAR=$(X=${POS1} Y=${POS2} Z=$x awk '{if (($4 < ENVIRON["X"]) && ($4 > ENVIRON["Y"]) && ($5 == ENVIRON["Z"])) print $0}' sample_specific_rfd_windows/${NOD}.input.anno.${REGIONS}.bed | bedtools getfasta -s -tab -name+ -fi ${GENOME} -bed stdin | perl ${BIN}/streamTriNucCounter.pl | sed '/^$/d' | tail -n 1 | cut -f2-65)
 echo $POS3 $VAR | sed 's/ /\t/g' >> ${REGIONS}/tmp_${NOD}.rfd.${x}.tnt_count
done

if [[ $COMBINE == 1 ]]
then
 if [ ! -e ${REGIONS}/${NOD}.rfd.COM.tnt_count ]
  then python ${BIN}/combine_tnt_counts.py ${REGIONS}/tmp_${NOD}.rfd.C.tnt_count ${REGIONS}/tmp_${NOD}.rfd.W.tnt_count > ${REGIONS}/${NOD}.rfd.COM.tnt_count
 else python ${BIN}/combine_tnt_counts.py ${REGIONS}/tmp_${NOD}.rfd.C.tnt_count ${REGIONS}/tmp_${NOD}.rfd.W.tnt_count | tail -n 1 >> ${REGIONS}/${NOD}.rfd.COM.tnt_count
 fi
else
 if [ ! -e ${REGIONS}/${NOD}.rfd.C.tnt_count ]
  then python ${BIN}/revcom_tnt_counts.py ${REGIONS}/tmp_${NOD}.rfd.C.tnt_count > ${REGIONS}/${NOD}.rfd.C.tnt_count
 else
  python ${BIN}/revcom_tnt_counts.py ${REGIONS}/tmp_${NOD}.rfd.C.tnt_count | sed '/^$/d' | tail -n 1 >> ${REGIONS}/${NOD}.rfd.C.tnt_count
 fi
 if [ ! -e ${REGIONS}/${NOD}.rfd.W.tnt_count ]
  then python ${BIN}/tnt_autosome_corrector.py ${REGIONS}/tmp_${NOD}.rfd.W.tnt_count > ${REGIONS}/${NOD}.rfd.W.tnt_count
 else
  python ${BIN}/tnt_autosome_corrector.py ${REGIONS}/tmp_${NOD}.rfd.W.tnt_count | sed '/^$/d' | tail -n 1 >> ${REGIONS}/${NOD}.rfd.W.tnt_count
 fi
fi

rm ${REGIONS}/tmp_${NOD}.rfd.C.tnt_count
rm ${REGIONS}/tmp_${NOD}.rfd.W.tnt_count
