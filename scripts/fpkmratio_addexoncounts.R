suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser(description = "Process some integers")
parser$add_argument("ratioresult", metavar = "ratio fpkm outrider result file", nargs = "+",
                    help = "The ratio fpkm outrider result file.")
parser$add_argument("exonresult", metavar = "exon outrider result file", nargs = "+",
                    help = "The exon outrider result file.")
parser$add_argument("-o","--output_path", metavar = "output_path", nargs = "+",
                    help = "Path where output file with replaced normalized and mean correct values.", default="./")
args <- parser$parse_args()

#url <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_ratiofpkm/umcu_rnaseq_fib_untreated_res_outrider_exons_counts.tsv"
umc_res_fpkm <- read_delim(args$ratioresult,show_col_types = FALSE)

#url <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_cntrl/umcu_rnaseq_res_outrider_exon_counts_L00022_pheno_added.tsv"
umc_res_e <- read_delim(args$exonresult,show_col_types = FALSE)

umc_res_e <- umc_res_e[,c(1,2,7:9)]
umc_res_merge <- right_join(umc_res_e, umc_res_fpkm, by=c("geneID","sampleID"))
res_merge_out <- umc_res_merge[,c(1:9,13:25)]
colnames(res_merge_out)[3:5] <- c("rawcounts","normcounts","meanCorrected")
#write_csv(res_merge_out, "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_ratiofpkm/umcu_rnaseq_fib_untreated_res_outrider_exons_counts_countsadded.tsv", append=FALSE, col_names = TRUE)
write_csv(res_merge_out, paste0(args$output_path,"umcu_rnaseq_fib_untreated_res_outrider_exons_counts_countsadded.tsv"), append=FALSE, col_names = TRUE)
