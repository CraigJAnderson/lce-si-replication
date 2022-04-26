#snakefile: /hps/nobackup2/flicek/user/cander21/strandAsym/release/snakefile
import pandas as pd

configfile : "config_lcesi.yaml"
#include : "rules/test.smk"

os.mkdir(config["PTH"]+"/"+config["DIR"])
os.chdir(config["PTH"]+"/"+config["DIR"])
os.mkdir("variants")
os.mkdir("sample_specific_rfd_windows")
os.mkdir("wac")

sample_names = pd.read_table(config["SAMPLE_LIST"],sep=" ",header=None)
sample_names.columns = ['nod']
samples= list(sample_names.nod)

boots = list(range(1,int(config["BOOTSTRAP"])+1))

if (config["COMBINE"] == 0):
 combine = ['.W']
elif (config["COMBINE"] == 1):
 combine = ['.COM']
 print("combining!")

regions = ["genic","nongenic"]

rule all:
 input:
  pt1 = expand("{pth}/{dir}/wac/{nod}.wac.bed", pth=config["PTH"], dir=config["DIR"], nod=samples), 
  pt2 = expand("{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.bed", pth=config["PTH"], dir=config["DIR"], nod=samples), 
  pt3 = expand("{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.{reg}.bed", reg=regions, nod=samples, pth=config["PTH"], dir=config["DIR"]),
  pt4 = expand("{pth}/{dir}/{reg}/{nod}.rfd{com}.tnt_count", nod=samples, com=combine, reg=regions, pth=config["PTH"], dir=config["DIR"]),
  pt5 = expand("{pth}/{dir}/{reg}/bs0/bs0.{nod}.mu", nod=samples, pth=config["PTH"], reg=regions, dir=config["DIR"]),
  pt6 = expand("{pth}/{dir}/{reg}/bs{bs}/bs{bs}.{nod}.mu", nod=samples, bs=boots, reg=regions, pth=config["PTH"], dir=config["DIR"]),
  
ruleorder: make_wac > pt1_RFD_window_sorter > subset_beds > pt2_tnt_counter > pt3_mut_count_obs > pt3_mut_count_BS

rule make_wac:
 output:
  "{pth}/{dir}/wac/{nod}.wac.bed"
 input:
  INPUT = expand("{fnod}/{nod}.fnodDrcrMa", nod=samples, fnod=config["FNOD_DRCR_MA_LOC"])
 params:
  PTH = config["PTH"],
  DIR = config["DIR"],
 shell:
  """awk -F"," '{{print $3"\\t"$6"\\t"$7"\\t"$15}}' {input} | sed '/^$/d' | awk '{{if ($4 >= 0.33) {{print $1"\\t"$2"\\t"$3"\\tW"}} else if ($4 <= -0.33) print $1"\\t"$2"\\t"$3"\\tC"}}' | sort -k1,1V -k2,2g | grep -v chr > {output}"""

rule pt1_RFD_window_sorter:
 output:
  "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.bed"
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
  OUT_BED1 = "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.{reg}.bed"
 params:
  GENE_ANNO = config["GENE_ANNO"]
 shell :
  """ bedtools subtract -a {input.IN_BED1} -b {params.GENE_ANNO} > {output.OUT_BED1} """

rule pt2_tnt_counter:
 output:
  "{pth}/{dir}/{reg}/{nod}.rfd{com}.tnt_count"
 input:
  "{pth}/{dir}/sample_specific_rfd_windows/{nod}.input.anno.{reg}.bed"
 params:
  BIN = config["BIN"], 
  COMBINE = config["COMBINE"], 
  DIR = config["DIR"],
  RFD = config["RFD_LIST"],
  GENOME = config["GENOME"]
 shell :
  """ while read line ; do bits=($line) ; {params.BIN}/pt2_tnt_counter.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} {params.COMBINE} {wildcards.reg} {wildcards.nod} {params.BIN} {params.GENOME} ; done < {params.RFD} """

rule pt3_mut_count_obs:
 output:
  "{pth}/{dir}/{reg}/bs0/bs0.{nod}.mu",
  temp("{pth}/{dir}/{reg}/bs0/{nod}.rfd.0.C.mut"),
  temp("{pth}/{dir}/{reg}/bs0/{nod}.rfd.0.W.mut")
 input:
  expand("{pth}/{dir}/{reg}/{nod}.rfd{com}.tnt_count", reg=regions, nod=samples, com=combine, pth=config["PTH"], dir=config["DIR"]),
 params:
  BIN = config["BIN"],
  COMBINE = config["COMBINE"],
  DIR = config["DIR"],
  RFD = config["RFD_LIST"],
  PTH = config["PTH"],
  GENOME = config["GENOME"]
 shell:
  """while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 0 {params.COMBINE} 0 {wildcards.reg} {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} """

rule pt3_mut_count_BS:
 output:
  "{pth}/{dir}/{reg}/bs{bs}/bs{bs}.{nod}.mu",
  temp("{pth}/{dir}/{reg}/bs{bs}/{nod}.{bs}.bed"),
  temp("{pth}/{dir}/{reg}/bs{bs}/{nod}.rfd.W.tnt_count"),
  temp("{pth}/{dir}/{reg}/bs{bs}/{nod}.rfd.C.tnt_count"),
  temp("{pth}/{dir}/{reg}/bs{bs}/{nod}.rfd.{bs}.C.mut"),
  temp("{pth}/{dir}/{reg}/bs{bs}/{nod}.rfd.{bs}.W.mut")
 input:
  expand("{pth}/{dir}/{reg}/{nod}.rfd{com}.tnt_count", reg=regions, nod=samples, com=combine, pth=config["PTH"], dir=config["DIR"]),
 params:
  BIN = config["BIN"],
  COMBINE = config["COMBINE"],
  DIR = config["DIR"],
  RFD = config["RFD_LIST"],
  PTH = config["PTH"],
  GENOME = config["GENOME"]
 shell:
  """while read line ; do bits=($line) ; {params.BIN}/pt3_mut_count.sh ${{bits[0]}} ${{bits[1]}} ${{bits[2]}} 1 {params.COMBINE} {wildcards.bs} {wildcards.reg} {wildcards.nod} {params.BIN} {params.PTH} {params.DIR} {params.GENOME} ; done < {params.RFD} """
