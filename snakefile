#snakefile: /hps/nobackup2/flicek/user/cander21/strandAsym/release/snakefile
import pandas as pd
import os

configfile : "config_lcesi.yaml"
#include : "rules/test.smk"

if not os.path.exists(config["PTH"]+"/"+config["DIR"]):
 os.makedirs(config["PTH"]+"/"+config["DIR"])

os.chdir(config["PTH"]+"/"+config["DIR"])

if not os.path.exists("variants"):
 os.mkdir("variants")

if not os.path.exists("sample_specific_rfd_windows"):
 os.mkdir("sample_specific_rfd_windows")

if not os.path.exists("wac"):
 os.mkdir("wac")

sample_names = pd.read_table(config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['nod']
samples= list(sample_names.nod)

boots = list(range(1,int(config["BOOTSTRAP"])+1))

rule all:
 input:
  pt1 = expand("{pth}/{dir}/wac/{nod}.wac.bed", pth=config["PTH"], dir=config["DIR"], nod=samples), 
  pt2 = expand("{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.bed", pth=config["PTH"], dir=config["DIR"], nod=samples), 
  pt3 = expand("{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.genic.bed", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  pt4 = expand("{pth}/{dir}/genic/{nod}.rfd.COM.tnt_count", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  pt5 = expand("{pth}/{dir}/genic/bs0/bs0.{nod}.mu", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  #pt6 = expand("{pth}/{dir}/genic/bs{bs}/bs{bs}.{nod}.mu", nod=samples, bs=boots, pth=config["PTH"], dir=config["DIR"]),
  pt7 = expand("{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.nongenic.bed", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  pt8 = expand("{pth}/{dir}/nongenic/{nod}.rfd.COM.tnt_count", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  pt9 = expand("{pth}/{dir}/nongenic/bs0/bs0.{nod}.mu", nod=samples, pth=config["PTH"], dir=config["DIR"]),
  #pt10 = expand("{pth}/{dir}/nongenic/bs{bs}/bs{bs}.{nod}.mu", nod=samples, bs=boots, pth=config["PTH"], dir=config["DIR"]),

  
ruleorder: make_wac > pt1_RFD_window_sorter > subset_beds > pt2_tnt_counter > pt3_mut_count_obs #> pt3_mut_count_BS_G > pt3_mut_count_BS_NG

rule make_wac:
 output:
  "{pth}/{dir}/wac/{nod}.wac.bed"
 #input:
 # "{fnod}/{nod}.fnodDrcrMa"
 params:
  PTH = config["PTH"],
  DIR = config["DIR"],
  FNOD = config["FNOD_DRCR_MA_LOC"]
 shell:
  """awk -F"," '{{print $3"\\t"$6"\\t"$7"\\t"$15}}' {params.FNOD}/{wildcards.nod}.fnodDrcrMa | sed '/^$/d' | awk '{{if ($4 >= 0.33) {{print $1"\\t"$2"\\t"$3"\\tW"}} else if ($4 <= -0.33) print $1"\\t"$2"\\t"$3"\\tC"}}' | sort -k1,1V -k2,2g | grep -v chr > {output}"""

rule pt1_RFD_window_sorter:
 output:
  "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.bed"
 input:
  "{pth}/{dir}/wac/{nod}.wac.bed"
 params:
  BIN = config["BIN"],
  PTH = config["PTH"],
  DIR = config["DIR"],
  RFD_VALUES = config["RFD_VALUES"],
  VARIANT_PATH = config["VARIANT_PATH"]
 shell :
  """{params.BIN}/pt1_RFD_window_sorter.sh {params.PTH} {params.DIR} {params.RFD_VALUES} {wildcards.nod} {params.VARIANT_PATH}"""

rule subset_beds:
 input:
  IN_BED1 = "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.bed"
 output:
  OUT_BED1 = "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.genic.bed",
  OUT_BED2 = "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.nongenic.bed"
 params:
  GENE_ANNO = config["GENE_ANNO"]
 shell :
  """ bedtools intersect -a {input.IN_BED1} -b {params.GENE_ANNO} > {output.OUT_BED1} ; bedtools subtract -a {input.IN_BED1} -b {params.GENE_ANNO} > {output.OUT_BED2} """

rule pt2_tnt_counter:
 output:
  "{pth}/{dir}/nongenic/{nod}.rfd.COM.tnt_count",
  "{pth}/{dir}/genic/{nod}.rfd.COM.tnt_count"
 input:
  "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.genic.bed",
  "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.nongenic.bed"
 params:
  BIN = config["BIN"],  
  DIR = config["DIR"],
  RFD = config["RFD_LIST"],
  GENOME = config["GENOME"]
 shell :
  """ while read line ; do bits=($line) ; {params.BIN}/pt2_tnt_counter.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 1 genic {wildcards.nod} {params.BIN} {params.GENOME} ; done < {params.RFD} ; while read line ; do bits=($line) ; {params.BIN}/pt2_tnt_counter.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 1 nongenic {wildcards.nod} {params.BIN} {params.GENOME} ; done < {params.RFD} """

rule pt3_mut_count_obs:
 output:
  "{pth}/{dir}/genic/bs0/bs0.{nod}.mu",
  #temp("{pth}/{dir}/genic/bs0/{nod}.rfd.0.C.mut"),
  #temp("{pth}/{dir}/genic/bs0/{nod}.rfd.0.W.mut"),
  "{pth}/{dir}/nongenic/bs0/bs0.{nod}.mu",
  #temp("{pth}/{dir}/nongenic/bs0/{nod}.rfd.0.C.mut"),
  #temp("{pth}/{dir}/nongenic/bs0/{nod}.rfd.0.W.mut")
 input:
  "{pth}/{dir}/genic/{nod}.rfd.COM.tnt_count",
  "{pth}/{dir}/nongenic/{nod}.rfd.COM.tnt_count"
  #expand("{pth}/{dir}/{reg}/{nod}.rfd{com}.tnt_count", reg=regions, nod=samples, com=combine, pth=config["PTH"], dir=config["DIR"]),
 params:
  BIN = config["BIN"],
  DIR = config["DIR"],
  RFD = config["RFD_LIST"],
  PTH = config["PTH"],
  GENOME = config["GENOME"]
 shell:
  """while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 0 1 0 genic {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} ; while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 0 1 0 nongenic {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} """

#rule pt3_mut_count_BS_G:
# output:
#  "{pth}/{dir}/genic/bs{bs}/bs{bs}.{nod}.mu",
#  #temp("{pth}/{dir}/genic/bs{bs}/{nod}.{bs}.bed"),
#  #temp("{pth}/{dir}/genic/bs{bs}/{nod}.rfd.W.tnt_count"),
#  #temp("{pth}/{dir}/genic/bs{bs}/{nod}.rfd.C.tnt_count"),
#  #temp("{pth}/{dir}/genic/bs{bs}/{nod}.rfd.{bs}.C.mut"),
#  #temp("{pth}/{dir}/genic/bs{bs}/{nod}.rfd.{bs}.W.mut"),
# input:
#  "{pth}/{dir}/genic/bs0/bs0.{nod}.mu",
# params:
#  BIN = config["BIN"],
#  DIR = config["DIR"],
#  RFD = config["RFD_LIST"],
#  PTH = config["PTH"],
#  GENOME = config["GENOME"]
# shell:
#  """ for x in {{1..3}} ; do while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 1 1 $x genic {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} ; done """
#
#rule pt3_mut_count_BS_NG:
# output:
#  "{pth}/{dir}/nongenic/bs{bs}/bs{bs}.{nod}.mu",
#  #temp("{pth}/{dir}/nongenic/bs{bs}/{nod}.{bs}.bed"),
#  #temp("{pth}/{dir}/nongenic/bs{bs}/{nod}.rfd.W.tnt_count"),
#  #temp("{pth}/{dir}/nongenic/bs{bs}/{nod}.rfd.C.tnt_count"),
#  #temp("{pth}/{dir}/nongenic/bs{bs}/{nod}.rfd.{bs}.C.mut"),
#  #temp("{pth}/{dir}/nongenic/bs{bs}/{nod}.rfd.{bs}.W.mut")
# input:
#  "{pth}/{dir}/nongenic/bs0/bs0.{nod}.mu"
# params:
#  BIN = config["BIN"],
#  DIR = config["DIR"],
#  RFD = config["RFD_LIST"],
#  PTH = config["PTH"],
#  GENOME = config["GENOME"]
# shell:
#  """ for x in {{1..3}} ; do while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 1 1 $x nongenic {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} ; done """