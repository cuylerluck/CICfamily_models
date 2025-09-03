#Script to perform DE analysis on EGAD00001003121 CIC-DUX4 and CIC-NUTM1 samples
#With questions, reach out to Cuyler Luck cuyler.luck@ucsf.edu or Ross Okimoto ross.okimoto@ucsf.edu
#This is not intended to be plug-and-play, but rather a guide for other users.
#It will likely not work out of the box to directly reproduce the data. At the least, changing file paths where necessary will be required.

#The DE performed uses a classic edgeR analysis with the exact test.


#Package versions included at the time of writing documentation
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(sva) #3.46.0
library(edgeR) #3.40.2
library(biomaRt) #2.54.1
library(ggrepel) #0.9.3
library(tidyr) #1.3.0

#set the working directory as needed
setwd("/Volumes/cuyler/ucsf_okimoto_lab/nutm1_dux4")

#define a list of sample names, this will be used for reading in data
samples = c("EWN1069T", "EWN1072T", "INI42", "RNA003_16_008", "RNA012_16_073", "RNA012_16_074", "SARC013", "SARC038", "SARC051", "SARC084", "SARC101")


#gather gene-level read counts into one matrix, using column 4 (counts for the 2nd read strand aligned with RNA, AKA htseq-count option -s reverse)
#this is the correct option for this data, which was generated with Illumina TruSeq Stranded mRNA protocol.
#see https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
#see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

#also gather uniquely mapped read percentages into one vector, for QC purposes
percents = c()
names = c()
for (sample in samples){
  temp = fread(paste(sample,"ReadsPerGene.out.tab",sep=""))
  genes_only = temp[grep("ENSG",temp$V1),]
  gene_stranded = dplyr::select(genes_only, V1, V4)
  colnames(gene_stranded) = c("gene",sample)
  
  if(!exists("master")){
    master=gene_stranded
  }
  else{
    master=dplyr::inner_join(master, gene_stranded, by="gene")
  }
  
  log_final_out = readLines(paste(sample,"Log.final.out",sep=""))
  percents = c(percents, log_final_out[grep("Uniquely mapped reads %",log_final_out)])
  names = c(names, sample)
  
}

#by looking at the vectors "names" and "percents" you can see what the % uniquely mapped reads was for each sample
#note that the RNAxxx... samples have a lower % uniquely mapped reads than the rest (~76% vs. low 90s% for the rest)

#retrieve metadata, which I have already cleaned
metadata = fread("NUTM1_vs_DUX4_RNAseq_metadata.csv", header = T)
metadata = metadata[,-1]

#the above metadata duplicates info for each sample since its based on the fastq file names, so let's simplify by taking info for just one line per sample
#I also only really care about partner gene and sequencer model here
simplified_meta = metadata[seq(1,21,2),]
simplified_meta = dplyr::select(simplified_meta,"sampleID","partner","sequencer_model")


#the three samples with a lower % uniquely mapped reads were sequenced on a different machine than the rest. 

### Batch correcting with ComBat_seq 
#This starts from the completely raw data again, with all genes.
#Batch used is sequencer model, groups are defined as the DUX4 and NUTM1 samples. 
batch_input_data = master
#Make into a matrix
batch_input_data = tibble::column_to_rownames(batch_input_data, "gene")
batch_input_data = as.matrix(batch_input_data)

#order metadata per matrix and define batch / group variables
simplified_meta_batch_corrected = tibble::column_to_rownames(simplified_meta,"sampleID")
simplified_meta_batch_corrected = simplified_meta_batch_corrected[colnames(batch_input_data),]

batches = ifelse(simplified_meta_batch_corrected$sequencer_model == "Illumina HiSeq 2500", 1, 2)
groups_toPreserve = ifelse(simplified_meta_batch_corrected$partner == "DUX4", 1, 2)

#providing sequencer model as the batch, and telling ComBat_seq what the biological group of interest is (fusion partner)
corrected = ComBat_seq(batch_input_data, batch=batches, group = groups_toPreserve)



### Standard edgeR analysis using batch corrected data from above
group = factor(simplified_meta_batch_corrected$partner)
colnames(corrected)    #use to double check that samples and metadata are in the same order
rownames(simplified_meta_batch_corrected)      #use to double check that samples and metadata are in the same order
y = DGEList(counts=corrected,group=group)
keep = filterByExpr(y)
y = y[keep,,keep.lib.sizes=F]
y = calcNormFactors(y)
y = estimateDisp(y)
et = exactTest(y)
#move results into an easier-to-use dataframe
results = et$table
results = tibble::rownames_to_column(results,"gene")
#manually adjust p-values using the FDR method
results = dplyr::mutate(results, p_fdr_adj = p.adjust(PValue, method = "fdr"))
#and also calculate -log10(p_fdr_adj) for use in volcano plot
results = dplyr::mutate(results, negLog10_p_fdr_adj = -log10(p_fdr_adj))


###Translate gene names from ENSG IDs to gene symbols using biomaRt
#Using Ensembl genes for GRCh38.p13.
#This will almost definitely lose some genes for which multiple symbols match or for which no symbols exist

ensembl = useEnsembl(biomart = "ensembl")
listDatasets(ensembl) #this shows that at the time of analysis, the 'hsapiens_gene_ensembl' dataset is for GRCh38.p13 which is the same version I used to align to.
ensembl = useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')

listFilters(ensembl) #I will filter by ensembl_gene_id
listAttributes(ensembl) #I will retrieve hgnc_symbol and ensembl_gene_id
#actually get the key from biomaRt
symbol_key = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id",values=results$gene, mart=ensembl)
#trim out entries that have no gene symbol
symbol_key = symbol_key[symbol_key$hgnc_symbol!="",]

#make sure there are unique matches
length(symbol_key$hgnc_symbol)
length(unique(symbol_key$hgnc_symbol))
length(unique(symbol_key$ensembl_gene_id))
#there are duplicate symbols & IDs, let's get rid of them
symbol_key = symbol_key[!duplicated(symbol_key$hgnc_symbol),]
symbol_key = symbol_key[!duplicated(symbol_key$ensembl_gene_id),]

#make a new DE results table that just includes genes that have symbols instead of ENSG IDs.
colnames(symbol_key) = c("symbol","gene")
symbol_results = left_join(symbol_key, results, by="gene")
#this loses about 3500 genes for which there were ENSG IDs but which did not map nicely to gene symbols.

#write this table to a file
write.table(symbol_results, file = "DE_DUX4_NUTM1.tsv",sep = "\t", quote = F, row.names = F)


#volcano plot for FOXB1
pdf(file = "for_Sep2025_revisions/foxb1_volcano.pdf", height = 5, width= 5)
ggplot(symbol_results, aes(x=logFC, y=negLog10_p_fdr_adj)) + 
  geom_point(color = "dark gray", size = 4, alpha = 0.25) +
  geom_point(data = symbol_results[symbol_results$symbol %in% c('FOXB1','ETV4','ETV5','ETV1','NUTM1'),], color = "blue", size = 4) +
  geom_hline(yintercept=1, linetype = "dashed", alpha = 0.5) + 
  geom_vline(xintercept=1, linetype = "dashed", alpha = 0.5) + 
  geom_vline(xintercept=-1, linetype = "dashed", alpha = 0.5) +
  geom_label_repel(data =  symbol_results[symbol_results$symbol %in% c('FOXB1','ETV4','ETV5','ETV1','NUTM1'),], aes(label = symbol), nudge_y = 7) + #labels
  theme_classic(base_size = 20) +
  xlab("log2 Fold Change") +
  ylab("-log10(q)")
dev.off()

#Pull out cpm for FOXB1 from batch corrected data, and also from the original raw read counts, and plot.
cpm = as.data.frame(cpm(y))
cpm_foxb1 = cpm['ENSG00000171956',]

cpm_foxb1_pivot = pivot_longer(cpm_foxb1, cols = c(1:11), names_to = "sample", values_to = "batch_corrected_cpm")
cpm_foxb1_pivot = dplyr::mutate(cpm_foxb1_pivot, partner = ifelse(sample %in% c("SARC084","SARC101","RNA003_16_008","RNA012_16_073","RNA012_16_074"), "NUTM1", "DUX4"))

#calculate the mean cpm per partner too, for plotting
mean_foxb1_dux4 = mean(cpm_foxb1_pivot$batch_corrected_cpm[cpm_foxb1_pivot$partner == "DUX4"])
mean_foxb1_nutm1 = mean(cpm_foxb1_pivot$batch_corrected_cpm[cpm_foxb1_pivot$partner == "NUTM1"])
means_foxb1 = data.frame(partner = c("DUX4","NUTM1"), means = c(mean_foxb1_dux4,mean_foxb1_nutm1))


raw_counts_foxb1 = master[master$gene == 'ENSG00000171956',]
raw_counts_foxb1 = tibble::column_to_rownames(raw_counts_foxb1, "gene")

raw_counts_foxb1_pivot = pivot_longer(raw_counts_foxb1, cols = c(1:11), names_to = "sample", values_to = "raw_counts")
raw_counts_foxb1_pivot = dplyr::mutate(raw_counts_foxb1_pivot, partner = ifelse(sample %in% c("SARC084","SARC101","RNA003_16_008","RNA012_16_073","RNA012_16_074"), "NUTM1", "DUX4"))

#make simple plots with these data to show FOXB1 expression in individual samples
set.seed(924152) #setting a seed so that jitter is always the same

pdf(file = "for_Sep2025_revisions/batchCorr_foxb1_cpm.pdf", width = 5, height = 4)
ggplot(data = cpm_foxb1_pivot, aes(x = partner, y = batch_corrected_cpm, color = partner)) +
  geom_crossbar(data = means_foxb1, aes(y=means,ymin=means,ymax=means), color = "black", width = 0.5) +
  geom_jitter(height = 0, width = 0.2, size = 4) + 
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks=seq(0,300,50)) +
  ylab("Batch Corrected, TMM-Normalized \n Counts Per Million Reads") +
  xlab("3' Fusion Partner") +
  scale_color_manual(values = c('#fc558e','#31a851'))
dev.off()
  
set.seed(924152) #setting a seed so that jitter is always the same

pdf(file = "for_Sep2025_revisions/raw_counts_foxb1.pdf", width = 5, height = 4)
ggplot(data = raw_counts_foxb1_pivot, aes(x = partner, y = raw_counts, color = partner)) +
  geom_jitter(height = 0, width = 0.2, size = 4) + 
  theme_classic(base_size = 16) +
  ylab("Raw Gene Counts") +
  xlab("3' Fusion Partner") +
  scale_color_manual(values = c('#fc558e','#31a851'))
dev.off() 


#let's also do one for ETV5 using the batch corrected, TMM-normalized CPM values
cpm_etv5 = cpm['ENSG00000244405',]

cpm_etv5_pivot = pivot_longer(cpm_etv5, cols = c(1:11), names_to = "sample", values_to = "batch_corrected_cpm")
cpm_etv5_pivot = dplyr::mutate(cpm_etv5_pivot, partner = ifelse(sample %in% c("SARC084","SARC101","RNA003_16_008","RNA012_16_073","RNA012_16_074"), "NUTM1", "DUX4"))

#calculate the mean cpm per partner too, for plotting
mean_etv5_dux4 = mean(cpm_etv5_pivot$batch_corrected_cpm[cpm_etv5_pivot$partner == "DUX4"])
mean_etv5_nutm1 = mean(cpm_etv5_pivot$batch_corrected_cpm[cpm_etv5_pivot$partner == "NUTM1"])
means_etv5 = data.frame(partner = c("DUX4","NUTM1"), means = c(mean_etv5_dux4,mean_etv5_nutm1))

set.seed(816152) #setting a seed so that jitter is always the same

pdf(file = "for_Sep2025_revisions/batchCorr_etv5_cpm.pdf", width = 5, height = 4)
ggplot(data = cpm_etv5_pivot, aes(x = partner, y = batch_corrected_cpm, color = partner)) +
  geom_crossbar(data = means_etv5, aes(y=means,ymin=means,ymax=means), color = "black", width = 0.5) +
  geom_jitter(height = 0, width = 0.2, size = 4) + 
  theme_classic(base_size = 16) +
  ylab("Batch Corrected, TMM-Normalized \n Counts Per Million Reads") +
  xlab("3' Fusion Partner") +
  scale_color_manual(values = c('#fc558e','#31a851'))
dev.off()

  