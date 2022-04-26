#!/bin/bash

PTH=$1
DIR=$2
RFD_VALUES=$3
NOD=$4
VARIANT_PATH=$5

#this will generate mutation sets for each strand, for each nodule, as well as orientating RFD windows with resepct to watson and crick.
cat ${PTH}/${DIR}/wac/${NOD}.wac.bed | awk  '{if ($4 == "W") print $1"\t"$2"\t"$3"\t"$4}' | bedtools intersect -a ${VARIANT_PATH}/${NOD}.filtered.nodMat -b stdin | sort -k1,1V -k2,2g | awk '{print $1"\t"$2"\t"$3"\t"$4"\tW"}' > ${PTH}/${DIR}/variants/${NOD}.filtered.W.nodMat
cat ${PTH}/${DIR}/wac/${NOD}.wac.bed | awk  '{if ($4 == "C") print $1"\t"$2"\t"$3"\t"$4}' | bedtools intersect -a ${VARIANT_PATH}/${NOD}.filtered.nodMat -b stdin | sort -k1,1V -k2,2g | awk '{print $1"\t"$2"\t"$3"\t"$4"\tC"}' > ${PTH}/${DIR}/variants/${NOD}.filtered.C.nodMat
cat ${PTH}/${DIR}/variants/${NOD}.filtered.C.nodMat ${PTH}/${DIR}/variants/${NOD}.filtered.W.nodMat | sort -k1,1V -k2,2g > ${PTH}/${DIR}/variants/${NOD}.filtered.nodMat
rm ${PTH}/${DIR}/variants/${NOD}.filtered.C.nodMat
rm ${PTH}/${DIR}/variants/${NOD}.filtered.W.nodMat
#set up okseq search space, keeping only windows that liftedover into those greater than 500 bp. The following also has an awk script to invert rfd values if they're on the watson.
awk '{if (($3-$2) > 500) print $0 }' ${RFD_VALUES} | bedtools intersect -a stdin -b ${PTH}/${DIR}/wac/${NOD}.wac.bed  -wb | cut -f1-4,8 | awk '{ if (($5 == "W") && ($4 > 0)) {print $1"\t"$2"\t"$3"\t"0-$4"\t"$5} else if (($5 == "W") && ($4 < 0)) {print $1"\t"$2"\t"$3"\t"$4-(2*$4)"\t"$5} else if ($5 == "C") {print $0}}' > sample_specific_rfd_windows/${NOD}.input.anno.bed

