#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu

#This analysis is on stranded RNA-seq data from 293T cells ~48h after transient transfection with one of three conditions:
# A = empty vector 
# B = HA-ATXN1::DUX4
# C = HA-ATXN1::DUX4 with AXH domain deleted

#I have triplicates for all conditions, hence the numbers in sample names

#First load packages and set working directory:
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(edgeR) #3.40.2
library(ggplot2) #3.5.2
library(pheatmap) #1.0.12
library(biomaRt) #2.54.1
library(tidyr) #1.3.0
library(ggrepel) #0.9.3
library(patchwork) #1.3.2
library(gprofiler2) #0.2.3


#set wd
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_ATXN1_DUX4/hg38_data")

#Next read in all data
#Skip the first four lines, they just give info on odd-mapping statistics
A1 = fread(file="ReadsPerGene/A1ReadsPerGene.out.tab", skip = 4)
A2 = fread(file="ReadsPerGene/A2ReadsPerGene.out.tab", skip = 4)
A3 = fread(file="ReadsPerGene/A3ReadsPerGene.out.tab", skip = 4)
B1 = fread(file="ReadsPerGene/B1ReadsPerGene.out.tab", skip = 4)
B2 = fread(file="ReadsPerGene/B2ReadsPerGene.out.tab", skip = 4)
B3 = fread(file="ReadsPerGene/B3ReadsPerGene.out.tab", skip = 4)
C1 = fread(file="ReadsPerGene/C1ReadsPerGene.out.tab", skip = 4)
C2 = fread(file="ReadsPerGene/C2ReadsPerGene.out.tab", skip = 4)
C3 = fread(file="ReadsPerGene/C3ReadsPerGene.out.tab", skip = 4)

#This is stranded data generated using the NEBNext Ultra II Directional Library Prep kit for Illumina
#so, we can take column 4 from STAR as the output (equivalent to HTseq -stranded "reverse")
#see https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
#see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
A1_strand = dplyr::select(A1, c(1,4))
A2_strand = dplyr::select(A2, c(1,4))
A3_strand = dplyr::select(A3, c(1,4))
B1_strand = dplyr::select(B1, c(1,4))
B2_strand = dplyr::select(B2, c(1,4))
B3_strand = dplyr::select(B3, c(1,4))
C1_strand = dplyr::select(C1, c(1,4))
C2_strand = dplyr::select(C2, c(1,4))
C3_strand = dplyr::select(C3, c(1,4))

#Let's put names on the gene count columns that correspond to the sample so that we can merge these into one big data frame
colnames(A1_strand) = c("ENSG", "A1")
colnames(A2_strand) = c("ENSG", "A2")
colnames(A3_strand) = c("ENSG", "A3")
colnames(B1_strand) = c("ENSG", "B1")
colnames(B2_strand) = c("ENSG", "B2")
colnames(B3_strand) = c("ENSG", "B3")
colnames(C1_strand) = c("ENSG", "C1")
colnames(C2_strand) = c("ENSG", "C2")
colnames(C3_strand) = c("ENSG", "C3")


#join into one master data.frame
master = inner_join(A1_strand, A2_strand, by = "ENSG")
master = inner_join(master, A3_strand, by = "ENSG")
master = inner_join(master, B1_strand, by = "ENSG")
master = inner_join(master, B2_strand, by = "ENSG")
master = inner_join(master, B3_strand, by = "ENSG")
master = inner_join(master, C1_strand, by = "ENSG")
master = inner_join(master, C2_strand, by = "ENSG")
master = inner_join(master, C3_strand, by = "ENSG")


#move ENSG to rownames for edgeR syntax
master = tibble::column_to_rownames(master, "ENSG")

#Let's remove ATXN1 and DUX4 just to make sure the overexpression of fusions doesn't change any clustering.
master = master[!(rownames(master) %in% c("ENSG00000124788", "ENSG00000260596")),]

#We can now pass this master data frame into an edgeR pipeline
#I am choosing to use a GLM for this analysis. So I will use the non-classic pipeline.


groups = c("EV", "EV", "EV", 
           "ATXN1DUX4", "ATXN1DUX4", "ATXN1DUX4", 
           "delAXH", "delAXH", "delAXH")
dg = DGEList(counts=master, group = groups)
keep = filterByExpr(dg)
dg = dg[keep, , keep.lib.sizes = FALSE]
dg = calcNormFactors(dg)
#reordering the levels of the group factor so "EV" becomes the reference
dg$samples$group = relevel(dg$samples$group, ref="EV")
#checking the order with
#dg$samples$group
#reveals that the present levels are:
#EV, ATXN1DUX4, delAXH
design = model.matrix(~dg$samples$group, data = dg$samples)
#checking the design matrix shows that now the intercept (baseline condition) is EV.
#coefficient 2 = ATXN1DUX4
#coefficient 3 = delAXH
#this information will be useful in doing pairwise comparisons soon.
dg = estimateDisp(dg, design)
#plotBCV(dg) #to see what the BCV looks like if desired


#and we can ask edgeR to plot an MDS plot to see how well samples cluster

#set color palette for all plots in same order as the samples are in master
all_palette = rep(c("#a9a9a9", "#a020f0", "#8455fa"), each = 3)

pdf(file = "plots/plotMDS_all_samples.pdf", width = 5, height = 4)
plotMDS(dg, col = all_palette, 
        labels = c("EV1","EV2","EV3","AD4.1","AD4.2","AD4.3","delAXH.1","delAXH.2","delAXH.3"))
dev.off()


pdf(file = "plots/plotMDS_all_symbols.pdf", width = 6, height = 4)
par(mar=c(6.1,5.1,5.1,11.1))
plotMDS(dg, col = all_palette,
        pch = rep(c(15, 16, 17), each = 3), 
        cex = 1.5)
legend("topright", inset = c(-.80,0), 
       legend = c("EV", "HA-ATXN1::DUX4", "delAXH"), 
       pch = c(15, 16, 17),
       col = all_palette[c(1,3,5)],
       xpd = T,
       bty = "n")

dev.off()


#This suggests that something is off with the HA-ATXN1::DUX4 sample from replicate 3.
#To be safe, let's exclude the samples from replicate 3 and move forward with just duplicates for each sample (reps 1 and 2).

master_all9 = master #saving as a backup
master = master[,c(1:2,4:5,7:8)]

#now do a new GLM analysis with just the duplicates
#first, remove the original objects associated with processing all samples to avoid accidental reuse

rm(groups)
rm(keep)
rm(dg)
rm(design)


groups2 = c("EV", "EV", 
           "ATXN1DUX4", "ATXN1DUX4", 
           "delAXH", "delAXH")
dg2 = DGEList(counts=master, group = groups2)
keep2 = filterByExpr(dg2)
dg2 = dg2[keep2, , keep.lib.sizes = FALSE]
dg2 = calcNormFactors(dg2)
#reordering the levels of the group factor so "EV" becomes the reference
dg2$samples$group = relevel(dg2$samples$group, ref="EV")
#checking the order with
#dg2$samples$group
#reveals that the present levels are:
#EV, ATXN1DUX4, delAXH
design2 = model.matrix(~dg2$samples$group, data = dg2$samples)
#checking the design matrix shows that now the intercept (baseline condition) is EV.
#coefficient 2 = ATXN1DUX4
#coefficient 3 = delAXH
#this information will be useful in doing pairwise comparisons soon.
dg2 = estimateDisp(dg2, design2)
#plotBCV(dg2) #to see what the BCV looks like if desired


#now we can do pairwise comparisons between groups of interest
#with this design matrix, providing glmQLFTest a number for coef means "compare this coefficient to the baseline [EV]"
#first we get the fit
fit2 = glmQLFit(dg2, design2)

#compare ATXN1DUX4 to EV
qlf.ATXN1DUX4.EV_dups = glmQLFTest(fit2, coef= 2)

#compare delAXH to EV
qlf.delAXH.EV_dups = glmQLFTest(fit2, coef= 3)

#compare delAXH to ATXN1DUX4
qlf.delAXH.ATXN1DUX4_dups = glmQLFTest(fit2, contrast=c(0, -1, 1)) #i.e. 3 coefficient with 2 coefficient


#now we can pull out tables of results from these tests for use in volcano plots etc.
ATXN1DUX4.EV_results_dups = qlf.ATXN1DUX4.EV_dups$table

delAXH.EV_results_dups = qlf.delAXH.EV_dups$table

delAXH.ATXN1DUX4_results_dups = qlf.delAXH.ATXN1DUX4_dups$table

#Let's add FDR-adjusted p-values to each of these
ATXN1DUX4.EV_results_dups = dplyr::mutate(ATXN1DUX4.EV_results_dups, fdr_adj = p.adjust(PValue, method = "fdr"))

delAXH.EV_results_dups = dplyr::mutate(delAXH.EV_results_dups, fdr_adj = p.adjust(PValue, method = "fdr"))

delAXH.ATXN1DUX4_results_dups = dplyr::mutate(delAXH.ATXN1DUX4_results_dups, fdr_adj = p.adjust(PValue, method = "fdr"))


#let's also pull out the TMM-normalized log(cpm) values from the edgeR pipeline
#this will be useful later for looking at expression of specific genes
logcpm2 = as.data.frame(cpm(dg2, log = T))


#and we can ask edgeR to plot an MDS plot to see how well samples cluster

#set color palette for all plots in same order as the samples are in master
all_palette = rep(c("#a9a9a9", "#a020f0", "#8455fa"), each = 2)

pdf(file = "plots/plotMDS_dups_samples.pdf", width = 5, height = 4)
plotMDS(dg2, col = all_palette, 
        labels = c("EV1","EV2","AD4.1","AD4.2","delAXH.1","delAXH.2"))
dev.off()


pdf(file = "plots/plotMDS_dups_symbols.pdf", width = 6, height = 4)
par(mar=c(6.1,5.1,5.1,11.1))
plotMDS(dg2, col = all_palette,
        pch = rep(c(15, 16, 17), each = 2), 
        cex = 1.5)
legend("topright", inset = c(-.80,0), 
       legend = c("EV", "HA-ATXN1::DUX4", "delAXH"), 
       pch = c(15, 16, 17),
       col = all_palette[c(1,3,5)],
       xpd = T,
       bty = "n")

dev.off()



#before we do anything more, I want to add on the gene symbols (not just ENSG IDs) for simplicity

#Translate gene names from ENSG IDs to gene symbols using biomaRt
#Using Ensembl genes for GRCh38.p14, which is one version newer than what I used to align samples, but should be OK (same overall genome)
#This will almost definitely lose some genes in the key for which multiple symbols match or for which no symbols exist.

#ensembl = useEnsembl(biomart = "ensembl")
#listDatasets(ensembl) #this shows that at the time of analysis, the 'hsapiens_gene_ensembl' dataset is for GRCh38.p14.
ensembl = useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')

#listFilters(ensembl) #I will filter by ensembl_gene_id
#listAttributes(ensembl) #I will retrieve hgnc_symbol and ensembl_gene_id
#actually get the key from biomaRt
#I'm getting symbols for all genes that came out of STAR originally, using A1_strand
symbol_key = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id",values=A1_strand$ENSG, mart=ensembl)
#trim out entries that have no gene symbol
symbol_key = symbol_key[symbol_key$hgnc_symbol!="",]

#make sure there are unique matches
length(symbol_key$hgnc_symbol)
length(unique(symbol_key$hgnc_symbol))
length(unique(symbol_key$ensembl_gene_id))
#there are duplicate symbols and keys, let's get rid of them
symbol_key = symbol_key[!duplicated(symbol_key$hgnc_symbol),]
symbol_key = symbol_key[!duplicated(symbol_key$ensembl_gene_id),]
#check again
length(symbol_key$hgnc_symbol)
length(unique(symbol_key$hgnc_symbol))
length(unique(symbol_key$ensembl_gene_id))
#now all lengths match, good.


#now we can add gene symbol to the results DFs
#for these I will left join -- this will save rows even if they have no symbol
ATXN1DUX4.EV_results_dups = tibble::rownames_to_column(ATXN1DUX4.EV_results_dups, "ensembl_gene_id")
ATXN1DUX4.EV_results_dups = left_join(ATXN1DUX4.EV_results_dups, symbol_key, by = "ensembl_gene_id")

delAXH.EV_results_dups = tibble::rownames_to_column(delAXH.EV_results_dups, "ensembl_gene_id")
delAXH.EV_results_dups = left_join(delAXH.EV_results_dups, symbol_key, by = "ensembl_gene_id")

delAXH.ATXN1DUX4_results_dups = tibble::rownames_to_column(delAXH.ATXN1DUX4_results_dups, "ensembl_gene_id")
delAXH.ATXN1DUX4_results_dups = left_join(delAXH.ATXN1DUX4_results_dups, symbol_key, by = "ensembl_gene_id")


#and lets add the names of genes in the logcpm dataframe. for this I only want genes that have symbols because I'm probably going to need symbols later
#so I will use inner_join

logcpm2 = tibble::rownames_to_column(logcpm2, "ensembl_gene_id")
logcpm2 = inner_join(logcpm2, symbol_key, by = "ensembl_gene_id")
logcpm2 = tibble::column_to_rownames(logcpm2, "hgnc_symbol")
#then get rid of ENSMBL IDs (currently in column 1)
logcpm2 = logcpm2[,-1]



#first, I will look at some volcano plots, and use these to define significantly upregulated genes.

#this plots ATXN1::DUX4 vs. EV results.
#FDR cutoff plotted is at -log10(q) = 3, or q < 0.001
#log2FC cutoffs are at +/- 1

ATXN1DUX4vsEVplot = ggplot(data = ATXN1DUX4.EV_results_dups, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1],size=4) +
  geom_point(data = ATXN1DUX4.EV_results_dups[ATXN1DUX4.EV_results_dups$logFC > 1 & -log10(ATXN1DUX4.EV_results_dups$fdr_adj) > 3,], 
             color = all_palette[3],size=4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = ATXN1DUX4.EV_results_dups[ATXN1DUX4.EV_results_dups$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4","ZSCAN4"),],
                  aes(label = hgnc_symbol), box.padding = 1, size = 8 ) +
  ggtitle("ATXN1DUX4 vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-5, 12), ylim = c(-2, 16))

pdf(file = "plots/AD4_vs_EV.pdf",width = 8, height = 6)
ATXN1DUX4vsEVplot
dev.off()

#this plots delAXH vs. EV results.
#FDR cutoff plotted is at -log10(q) = 3, or q < 0.001
#log2FC cutoffs are at +/- 1

delAXHvsEVplot = ggplot(data = delAXH.EV_results_dups, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1],size=4) +
  geom_point(data = delAXH.EV_results_dups[delAXH.EV_results_dups$logFC > 1 & -log10(delAXH.EV_results_dups$fdr_adj) > 3,], 
             color = all_palette[5], size=4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = delAXH.EV_results_dups[delAXH.EV_results_dups$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4","ZSCAN4"),],
                  aes(label = hgnc_symbol), box.padding = 1, size = 8 ) +
  ggtitle("delAXH vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-5, 12), ylim = c(-2, 16))

pdf(file = "plots/delAXH_vs_EV.pdf",width = 8, height = 6)
delAXHvsEVplot
dev.off()

pdf(file = "plots/patchwork_volcano.pdf", width = 8, height = 12)
ATXN1DUX4vsEVplot / delAXHvsEVplot
dev.off()


#now for heatmaps

#first, heatmap of all genes on in either the full length vs EV or delAXH vs EV comparisons
#cutoff for a gene being turned on is same as in volcano plots, ie log2FC > 1, -log10q > 3

on_full = ATXN1DUX4.EV_results_dups[ATXN1DUX4.EV_results_dups$logFC > 1 & -log10(ATXN1DUX4.EV_results_dups$fdr_adj) > 3 & !is.na(ATXN1DUX4.EV_results_dups$hgnc_symbol),]$hgnc_symbol
on_delAXH = delAXH.EV_results_dups[delAXH.EV_results_dups$logFC > 1 & -log10(delAXH.EV_results_dups$fdr_adj) > 3 & !is.na(delAXH.EV_results_dups$hgnc_symbol),]$hgnc_symbol
union_on = union(on_full, on_delAXH)

#separate logcpm with just the on genes
logcpm_union_on = logcpm2[union_on,]

#annotation dfs
sample_annotations = data.frame(condition = rep(c("EV", "HA-ATXN1::DUX4", "delAXH"), each = 2))
rownames(sample_annotations) = colnames(logcpm2)

heatmap_colors = list(
  condition = all_palette
)

names(heatmap_colors$condition) = sample_annotations$condition

pdf(file = "plots/heatmap_union_on.pdf", width = 6, height = 8)
pheatmap(logcpm_union_on, 
         cluster_cols = F, 
         scale = "row", 
         show_rownames = F,
         angle_col = 45,
         annotation_col = sample_annotations,
         annotation_colors = heatmap_colors)
dev.off()



#how about log2cpm plots

logcpm_pivot = tibble::rownames_to_column(logcpm2, "gene")
logcpm_pivot = pivot_longer(logcpm_pivot, cols = c(2:7), names_to = "sample", values_to = "log2cpm")

#add a "condition" variable to group each sample by
logcpm_pivot = mutate(logcpm_pivot, condition = ifelse(sample %in% c("A1","A2"), "EV", 
                                                       ifelse(sample %in% c("C1", "C2"), "delAXH", "ATXN1::DUX4")))

#make a simplified condition-color key to use in a second
condition_colors = heatmap_colors$condition[c(1,3,5)]
names(condition_colors) = c("EV", "ATXN1::DUX4", "delAXH")

#write a little function to make a log(cpm) plot given a gene and the samples to use (in the order you want them, as a vector)
make_cpm_plot = function(gene, groups){
  logcpm_pivot_subset = logcpm_pivot[logcpm_pivot$gene == gene & logcpm_pivot$condition %in% groups,]
  return(ggplot(logcpm_pivot_subset, 
                aes(x = factor(condition, levels = groups), y = log2cpm, fill = condition))+
           stat_summary(fun = mean, geom = "bar", show.legend = F) +
           geom_point(show.legend = F) +
           ggtitle(paste(gene)) +
           theme_classic() +
           scale_fill_manual(values = condition_colors[groups]) +
           xlab("Condition") +
           theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 12)))
}

#plot some interesting genes

selected = c("ETV5","ETV4","VGF","SPRY4","ZSCAN4")

plots_selected = lapply(selected, make_cpm_plot, groups = c("EV","ATXN1::DUX4","delAXH"))

#plot using patchwork
patchwork_plots_selected=wrap_plots(plots_selected[1:length(selected)]) + plot_layout(ncol=5)

pdf(file = "plots/patchwork_logcpm_selected.pdf",width = 10, height = 3)
patchwork_plots_selected
dev.off()



#one more volcano plot to see what genes are sig. different in delAXH vs full-length ATXN1::DUX4
#this plots delAXH vs. ATXN1::DUX4 results.
#FDR cutoff plotted is at -log10(q) = 3, or q < 0.001
#log2FC cutoffs are at +/- 1

delAXHvsATXN1DUX4plot = ggplot(data = delAXH.ATXN1DUX4_results_dups, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1],size=5) +
  geom_point(data = delAXH.ATXN1DUX4_results_dups[delAXH.ATXN1DUX4_results_dups$logFC < -1 & -log10(delAXH.ATXN1DUX4_results_dups$fdr_adj) > 3,], 
             color = all_palette[3],size=5) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = delAXH.ATXN1DUX4_results_dups[delAXH.ATXN1DUX4_results_dups$hgnc_symbol %in% c("ETV4","ETV5","MAFF","DUSP4","SPRY4"),],
                  aes(label = hgnc_symbol), box.padding = 1, size = 8 ) +
  ggtitle("delAXH vs. ATXN1::DUX4") +
  xlab("log2FC") +
  ylab("-log10(q)") 

pdf(file = "plots/delAXH_vs_AD4.pdf",width = 8, height = 6)
delAXHvsATXN1DUX4plot
dev.off()



#write out datasets to share

write.table(tibble::rownames_to_column(logcpm2, "Gene"), file = "outputs/logcpm_Feb2025.txt", sep = "\t", quote = F, row.names = F)
write.table(tibble::rownames_to_column(master, "Gene"), file = "outputs/counts_Feb2025.txt", sep = "\t", quote = F, row.names = F)
write.table(ATXN1DUX4.EV_results_dups, file = "outputs/ATXN1DUX4_vs_EV_Feb2025.txt", sep = "\t", quote = F, row.names = F)
write.table(delAXH.EV_results_dups, file = "outputs/delAXH_vs_EV_Feb2025.txt", sep = "\t", quote = F, row.names = F)
write.table(delAXH.ATXN1DUX4_results_dups, file = "outputs/delAXH_vs_ATXN1DUX4_Feb2025.txt", sep = "\t", quote = F, row.names = F)


#Reviewers asked for gene ontology analysis. Doing this with gProfiler, using the list of genes sig up in ATXN1::DUX4 vs EV.

#ATXN1DUX4 vs EV
ATXN1_gprofiler = gost(on_full, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
ATXN1_GO = ATXN1_gprofiler$result
trimmed_ATXN1_GO = dplyr::select(ATXN1_GO, p_value, term_id, term_name)
write.table(trimmed_ATXN1_GO, file = "outputs/ATXN1DUX4_vs_EV_gProfiler.txt", sep = "\t", quote = F, row.names = F)

#also writing out the list of genes sig. up in full length ATXN1::DUX4 vs EV
write.table(on_full, file = "outputs/up_ATXN1DUX4_vs_EV.txt", sep = "\t", quote = F, row.names = F, col.names = F)


