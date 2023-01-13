#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

preds = read.table(args[1], header=TRUE, sep="\t")
truths = read.table(args[2], header=TRUE, sep="\t")

#truths$Pair_name_a = paste0("ped1_",truths$V1,"-","ped1_", truths$V2)
#truths$Pair_name_b = paste0("ped1_",truths$V2,"-","ped1_", truths$V1)
merged = merge(preds, truths, by="Pair_name")

#output   = rbind(merged_a[,c("Most_Likely_rel", "V4")], merged_b[,c("Most_Likely_rel", "V4")])

#head(merged)
#colnames(merged)
cat(paste0(merged$Most_Likely_rel.x, " ", merged$Most_Likely_rel.y,"\n"), sep = "")
