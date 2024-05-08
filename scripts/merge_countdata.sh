#!/bin/bash
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }

while getopts "c:h:o:" opt; do
  case ${opt} in
    c ) # cntrl file in 
      in_cntrl=${OPTARG} ;;
    h ) # chx file in 
      in_chx=${OPTARG} ;;
    o ) # out gene or exon 
      out=${OPTARG} ;;
    \? | h | * )
      usage ;;
  esac
done
shift $((OPTIND -1))

#if diff <(awk '{print $1}' ${in_cntrl}) <(awk '{print $1}' ${in_chx}) >& /dev/null;
if diff <(cut -f1 ${in_cntrl}) <(cut -f1 ${in_chx}) >& /dev/null;
then
  cut -f2- ${in_chx} > tmp.tsv
  paste <(cat ${in_cntrl}) <(cat tmp.tsv) > umcu_rnaseq_${out}s_counts.tsv
  rm tmp.tsv
  echo "done"
else
  echo "first column not equal"
fi