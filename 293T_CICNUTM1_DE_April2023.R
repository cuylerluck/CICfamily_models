#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu


#This analysis is on stranded RNA-seq data from 293T cells ~48h after transient transfection with one of three conditions:
# A = empty vector (same backbone as CIC-NUTM1 plasmids)
# B = CIC(ex20)-NUTM1(ex6)
# D = CIC-DUX4 (different backbone from the other 3 plasmids, but same promoter)

#note there was originally a condition C, but it was discarded because it was likely subject to issues with transfection efficiency in this experiment

#I have triplicates for all three conditions, which are included

#hex colors:
#EV - #a9a9a9
#CIC::DUX4 - #fc558e
#CICex20::NUTM1ex6 - #31a851

#First load packages and set working directory:
library(data.table) #1.14.8
library(dplyr) #1.1.1
library(edgeR) #3.40.2
library(ggplot2) #3.4.1
library(pheatmap) #1.0.12
library(biomaRt) #2.54.1
library(tidyr) #1.3.0
library(ggrepel) #0.9.3
library(patchwork) #1.1.2
library(gprofiler2) #0.2.3


setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_CICNUTM1") #this will change by user & location of data

#Next read in all data
#Skip the first four lines, they just give info on odd-mapping statistics
A1 = fread(file="Counts/A1ReadsPerGene.out.tab", skip = 4)
A2 = fread(file="Counts/A2ReadsPerGene.out.tab", skip = 4)
A3 = fread(file="Counts/A3ReadsPerGene.out.tab", skip = 4)
B1 = fread(file="Counts/B1ReadsPerGene.out.tab", skip = 4)
B2 = fread(file="Counts/B2ReadsPerGene.out.tab", skip = 4)
B3 = fread(file="Counts/B3ReadsPerGene.out.tab", skip = 4)
D1 = fread(file="Counts/D1ReadsPerGene.out.tab", skip = 4)
D2 = fread(file="Counts/D2ReadsPerGene.out.tab", skip = 4)
D3 = fread(file="Counts/D3ReadsPerGene.out.tab", skip = 4)


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
D1_strand = dplyr::select(D1, c(1,4))
D2_strand = dplyr::select(D2, c(1,4))
D3_strand = dplyr::select(D3, c(1,4))

#Let's put names on the gene count columns that correspond to the sample so that we can merge these into one big data frame
colnames(A1_strand) = c("ENSG", "A1")
colnames(A2_strand) = c("ENSG", "A2")
colnames(A3_strand) = c("ENSG", "A3")
colnames(B1_strand) = c("ENSG", "B1")
colnames(B2_strand) = c("ENSG", "B2")
colnames(B3_strand) = c("ENSG", "B3")
colnames(D1_strand) = c("ENSG", "D1")
colnames(D2_strand) = c("ENSG", "D2")
colnames(D3_strand) = c("ENSG", "D3")

master = inner_join(A1_strand, A2_strand, by = "ENSG")
master = inner_join(master, A3_strand, by = "ENSG")
master = inner_join(master, B1_strand, by = "ENSG")
master = inner_join(master, B2_strand, by = "ENSG")
master = inner_join(master, B3_strand, by = "ENSG")
master = inner_join(master, D1_strand, by = "ENSG")
master = inner_join(master, D2_strand, by = "ENSG")
master = inner_join(master, D3_strand, by = "ENSG")

#move ENSG to rownames for edgeR syntax
master = tibble::column_to_rownames(master, "ENSG")

#Let's remove CIC, NUTM1, and DUX4 just to make sure the overexpression of fusions doesn't change anything.
master = master[!(rownames(master) %in% c("ENSG00000079432", "ENSG00000260596", "ENSG00000184507")),]

#We can now pass this master data frame into an edgeR pipeline
#I am choosing to use a GLM for this analysis. So I will use the non-classic pipeline.

groups = c("EV", "EV", "EV", "CICex20_NUTM1ex6", "CICex20_NUTM1ex6", "CICex20_NUTM1ex6", "CICDUX4", "CICDUX4", "CICDUX4")
dg = DGEList(counts=master, group = groups)
keep = filterByExpr(dg)
dg = dg[keep, , keep.lib.sizes = FALSE]
dg = calcNormFactors(dg)
#reordering the levels of the group factor so "EV" becomes the reference
dg$samples$group = relevel(dg$samples$group, ref="EV")
#checking the order with
#dg$samples$group
#reveals that the present levels are EV CICDUX4 CICex20_NUTM1ex6
design = model.matrix(~dg$samples$group, data = dg$samples)
#checking the design matrix shows that now the intercept (baseline condition) is EV.
#coefficient 2 = CICDUX4
#coefficient 3 = CICex20_NUTM1ex6
#this information will be useful in doing pairwise comparisons soon.
dg = estimateDisp(dg, design)
#plotBCV(dg) #to see what the BCV looks like if desired

#now we can do pairwise comparisons between groups of interest
#with this design matrix, providing glmQLFTest a number for coef means "compare this coefficient to the baseline [EV]"
#first we get the fit
fit = glmQLFit(dg, design)

#compare CICDUX4 to EV
qlf.CICDUX4.EV = glmQLFTest(fit, coef=2)

#compare CICex20_NUTM1ex6 to EV
qlf.CICex20_NUTM1ex6.EV = glmQLFTest(fit, coef=3)

#compare CICex20_NUTM1ex6 to CICDUX4
qlf.CICex20_NUTM1ex6.CICDUX4 = glmQLFTest(fit, contrast=c(0,-1,1)) #i.e. third coefficient with second coefficient

#now we can pull out tables of results from these tests for use in volcano plots etc.
CICDUX4.EV_results = qlf.CICDUX4.EV$table

CICex20_NUTM1ex6.EV_results = qlf.CICex20_NUTM1ex6.EV$table

CICex20_NUTM1ex6.CICDUX4 = qlf.CICex20_NUTM1ex6.CICDUX4$table


#Let's add FDR-adjusted p-values to each of these
CICDUX4.EV_results = dplyr::mutate(CICDUX4.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

CICex20_NUTM1ex6.EV_results = dplyr::mutate(CICex20_NUTM1ex6.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

CICex20_NUTM1ex6.CICDUX4 = dplyr::mutate(CICex20_NUTM1ex6.CICDUX4, fdr_adj = p.adjust(PValue, method = "fdr"))


#let's also pull out the TMM-normalized log(cpm) values from the edgeR pipeline
#this will be useful later for looking at expression of specific genes
logcpm = as.data.frame(cpm(dg, log = T))


#and we can ask edgeR to plot an MDS plot to see how well samples cluster
pdf(file = "plots/plotMDS_EV_DUX_NUTM1_samples.pdf", width = 5, height = 4)
plotMDS(dg, labels = c("EV1", "EV2", "EV3", "C20N6_1", "C20N6_2", "C20N6_3", "CD1", "CD2", "CD3"))
dev.off()

#here is a version with symbols, not sample names
pdf(file = "plots/plotMDS_EV_DUX_NUTM1_symbols.pdf", width = 6, height = 4)
par(mar=c(5.1,5.1,5.1,11.1))
plotMDS(dg, col = c(rep("#a9a9a9",3), rep("#31a851",3), rep("#fc558e",3)),
        pch = rep(c(15, 16, 17), each = 3), 
        cex = 2)
legend("topright", inset = c(-0.2,-0.53), 
       legend = c("EV", "HA-CICex20::NUTM1ex6", "HA-CIC::DUX4"), 
       pch = c(15, 16, 17),
       col = c("#a9a9a9", "#31a851", "#fc558e"),
       xpd = T,
       bty = "n",
       cex = 1.5)
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
CICDUX4.EV_results = tibble::rownames_to_column(CICDUX4.EV_results, "ensembl_gene_id")
CICDUX4.EV_results = left_join(CICDUX4.EV_results, symbol_key, by = "ensembl_gene_id")

CICex20_NUTM1ex6.EV_results = tibble::rownames_to_column(CICex20_NUTM1ex6.EV_results, "ensembl_gene_id")
CICex20_NUTM1ex6.EV_results = left_join(CICex20_NUTM1ex6.EV_results, symbol_key, by = "ensembl_gene_id")

CICex20_NUTM1ex6.CICDUX4 = tibble::rownames_to_column(CICex20_NUTM1ex6.CICDUX4, "ensembl_gene_id")
CICex20_NUTM1ex6.CICDUX4 = left_join(CICex20_NUTM1ex6.CICDUX4, symbol_key, by = "ensembl_gene_id")


#and lets add the names of genes in the logcpm dataframe. for this I only want genes that have symbols because I'm probably going to need symbols later
#so I will use inner_join

logcpm = tibble::rownames_to_column(logcpm, "ensembl_gene_id")
logcpm = inner_join(logcpm, symbol_key, by = "ensembl_gene_id")
logcpm = tibble::column_to_rownames(logcpm, "hgnc_symbol")
#then get rid of ENSMBL IDs (currently in column 1)
logcpm = logcpm[,-1]



#now we have all the base data structures necessary to make some plots with.



#first, I will look at some volcano plots, and use these to define significantly upregulated genes.


#this plots CIC-DUX4 vs. EV results.
#FDR cutoff plotted is at -log10(q) = 3, or q < 0.001.
#log2FC cutoffs are at +/- 1.5

CICDUX4vsEVplot = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "#a9a9a9") +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$logFC > 1.5 & -log10(CICDUX4.EV_results$fdr_adj) > 3,], 
             color = "#fc558e") +
  theme_classic() +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1.5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1.5, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5", "VGF", "SPRY4", "DUSP6"),],
    aes(label = hgnc_symbol), nudge_y = 1.5, nudge_x = -1, box.padding = 0.5) +
  ggtitle("CICDUX4 vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  xlim(c(-5,17)) +
  ylim(c(0,18))

#this plots CIC(ex20)-NUTM1(ex6) vs. EV results.
#FDR cutoff plotted is at -log10(q) = 3, or q < 0.001.
#log2FC cutoffs are at +/- 1.5

CICNUTM1vsEVplot = ggplot(data = CICex20_NUTM1ex6.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = "#a9a9a9") +
  geom_point(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$logFC > 1.5 & -log10(CICex20_NUTM1ex6.EV_results$fdr_adj) > 3,],
             color = "#31a851") +
  theme_classic() +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1.5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1.5, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5", "VGF", "SPRY4", "DUSP6"),],
                  aes(label = hgnc_symbol), nudge_y = 1, nudge_x = -1, box.padding = 0.5) +
  ggtitle("CIC(ex20)-NUTM1(ex6) vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  xlim(c(-5,17)) +
  ylim(c(0,18))


patchwork_volcano = CICDUX4vsEVplot / CICNUTM1vsEVplot 
pdf(file = "plots/volcanoes.pdf", width = 6, height = 8)
patchwork_volcano
dev.off()


#Now I'd like to do some heatmaps to get another way of looking at how gene expression varies across all conditions simultaneously.

#To look at significant genes changed relative to EV with either fusion,
# let's get IDs of all genes sig. diff over-expressed (q<0.001, logFC > 1.5) for either CIC-NUTM1 or CIC-DUX4, & hgnc_symbol is not NA

# Note that this restricts the data for the heatmaps to just genes being turned ON relative to EV.
# This is because we are looking at the fusions as activators of transcription.

de_CICex20_NUTM1ex6_EV = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$fdr_adj<0.001 & CICex20_NUTM1ex6.EV_results$logFC > 1.5 & !is.na(CICex20_NUTM1ex6.EV_results$hgnc_symbol),]
de_CICDUX4_EV = CICDUX4.EV_results[CICDUX4.EV_results$fdr_adj<0.001 & CICDUX4.EV_results$logFC > 1.5 & !is.na(CICDUX4.EV_results$hgnc_symbol),]
de_in_CN1_or_CD4 = unique(c(de_CICex20_NUTM1ex6_EV$hgnc_symbol, de_CICDUX4_EV$hgnc_symbol))

#Just making a duplicate of logcpm for use in heatmapping
heatmap_cpm = logcpm

#and defining vectors for annotating columns by experimental condition
annotation_col = data.frame(
  Condition = c(rep("EV",3), rep("CICNUTM1",3), rep("CICDUX4",3))
)
rownames(annotation_col) = colnames(heatmap_cpm)

annotation_colors = list(
  Condition = c(EV = "#a9a9a9", CICNUTM1 = "#31a851", CICDUX4 = "#fc558e")
)

#We should use rowmax/min to look at these... just coloring by log(cpm) will make it hard to see within-row differences. so, scale = "row" will do this

#Here is a heatmap for all genes DE in at least one of the two comparisons, per the above definition
pdf(file = "plots/de_in_at_least_one_comparison.pdf", width = 6, height = 5)
heatmap_on_at_least_one = pheatmap(heatmap_cpm[rownames(heatmap_cpm) %in% de_in_CN1_or_CD4,], 
         scale = "row", 
         show_rownames = F,
         cluster_cols = F,
         labels_col = c("EV 1","EV 2","EV 3","CICNUTM1 1","CICNUTM1 2","CICNUTM1 3","CICDUX4 1","CICDUX4 2","CICDUX4 3"),
         angle_col = 45,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors)
dev.off()



#it's interesting that there are a few blocks of genes that are specifically turned on in one of the two fusions only. 
#what are they?

#we already have dataframes containing the significantly turned ON (positive log2FC) genes for both fusions vs. EV.
#now which do not intersect (are exclusive to one condition)?

on_CICDUX4_only = de_CICDUX4_EV[!(de_CICDUX4_EV$hgnc_symbol %in% de_CICex20_NUTM1ex6_EV$hgnc_symbol),] 
#60 genes sig. on in CICDUX4 exclusively
on_CICex20_NUTM1ex6_only = de_CICex20_NUTM1ex6_EV[!(de_CICex20_NUTM1ex6_EV$hgnc_symbol %in% de_CICDUX4_EV$hgnc_symbol),] 
#178 genes sig. on in CIC-NUTM1 exclusively
on_both_fusions_temp = c(de_CICDUX4_EV$hgnc_symbol, de_CICex20_NUTM1ex6_EV$hgnc_symbol)
on_both_fusions = on_both_fusions_temp[duplicated(on_both_fusions_temp)] 
#226 genes significantly on in both

#let's write some text files for the fusion-specific and shared targets to do some GO with or share.

write.table(on_both_fusions, file = "outputs/genes_on_both_fusions.txt", row.names = F, quote = F, col.names = F)
write.table(on_CICDUX4_only$hgnc_symbol, file = "outputs/genes_on_CICDUX4_only.txt", row.names = F, quote = F, col.names = F)
write.table(on_CICex20_NUTM1ex6_only$hgnc_symbol, file = "outputs/gene_on_CICNUTM1_only.txt", row.names = F, quote = F, col.names = F)


#One concern about just using heatmaps to draw conclusions about genes being turned on is that
# the heatmaps are scaled to row max/min. This means that even objectively small log(cpm) expression values can look large.
#Let's test this directly by plotting logcpm itself for some targets of interest.

#hex colors:
#EV - #a9a9a9
#CIC::DUX4 - #fc558e
#CICex20::NUTM1ex6 - #31a851


#I need to pivot_longer to be able to plot the data properly
logcpm_pivot = tibble::rownames_to_column(logcpm, "ID")
logcpm_pivot = pivot_longer(logcpm_pivot, cols = c("A1","A2","A3","B1","B2","B3","D1","D2","D3"), names_to = "condition")
#need to also make a new variable just simplifying to A/B/D for the replicates, will make it easier to plot with ggplot
logcpm_pivot = mutate(logcpm_pivot, group = ifelse(condition %in% c("A1","A2","A3"), "EV", ifelse(condition %in% c("B1","B2","B3"),"CIC(ex20)-NUTM1(ex6)", ifelse(condition %in% c("D1","D2","D3"), "CIC-DUX4","other"))))

#change the levels of the group factor to set a better order for the x axis
logcpm_pivot$group = factor(logcpm_pivot$group, levels = c("EV", "CIC-DUX4", "CIC(ex20)-NUTM1(ex6)"))
#also adding a color variable, just to make ggplotting easier
logcpm_pivot = mutate(logcpm_pivot, color = ifelse(group == "EV", "#a9a9a9", ifelse(group == "CIC-DUX4", "#fc558e", ifelse(group == "CIC(ex20)-NUTM1(ex6)", "#31a851", "black"))))

#FOXD3
logcpm_pivot_FOXD3 = logcpm_pivot[logcpm_pivot$ID == "FOXD3",]
FOXD3_pivot = ggplot(data = logcpm_pivot_FOXD3, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("FOXD3") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

#FOXC2
logcpm_pivot_FOXC2 = logcpm_pivot[logcpm_pivot$ID == "FOXC2",]
FOXC2_pivot = ggplot(data = logcpm_pivot_FOXC2, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("FOXC2") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))


#FOXB1
logcpm_pivot_FOXB1 = logcpm_pivot[logcpm_pivot$ID == "FOXB1",]
FOXB1_pivot = ggplot(data = logcpm_pivot_FOXB1, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("FOXB1") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

#FOXG1
logcpm_pivot_FOXG1 = logcpm_pivot[logcpm_pivot$ID == "FOXG1",]
FOXG1_pivot = ggplot(data = logcpm_pivot_FOXG1, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("FOXG1") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

#SOX2
logcpm_pivot_SOX2 = logcpm_pivot[logcpm_pivot$ID == "SOX2",]
SOX2_pivot = ggplot(data = logcpm_pivot_SOX2, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("SOX2") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

patchwork_forkhead = FOXC2_pivot | FOXD3_pivot | FOXB1_pivot | FOXG1_pivot | SOX2_pivot

pdf(file = "plots/logcpm_forkheads_SOX2.pdf", width = 9, height = 3)
patchwork_forkhead
dev.off()


#and never hurts to look at ETV1/4/5

#ETV1
logcpm_pivot_ETV1 = logcpm_pivot[logcpm_pivot$ID == "ETV1",]
ETV1_pivot = ggplot(data = logcpm_pivot_ETV1, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("ETV1") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

#ETV4
logcpm_pivot_ETV4 = logcpm_pivot[logcpm_pivot$ID == "ETV4",]
ETV4_pivot = ggplot(data = logcpm_pivot_ETV4, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("ETV4") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

#ETV5
logcpm_pivot_ETV5 = logcpm_pivot[logcpm_pivot$ID == "ETV5",]
ETV5_pivot = ggplot(data = logcpm_pivot_ETV5, aes(x = group, y = value, fill = color)) + 
  geom_point() +
  stat_summary(fun = mean, geom = "bar", alpha = 0.3) +
  theme_classic() +
  ggtitle("ETV5") +
  ylab("log2(cpm)") +
  xlab("Group") +
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, hjust = 1, angle = 30),
        legend.position = "none") +
  scale_fill_manual(values = c("#31a851","#a9a9a9","#fc558e"))

patchwork_ETV = ETV1_pivot | ETV4_pivot | ETV5_pivot

pdf(file = "plots/logcpm_ETVs.pdf", width = 6, height = 3)
patchwork_ETV
dev.off()



#writing some data files to share
#logcpm, edgeR output

write.table(tibble::rownames_to_column(logcpm,"Gene"), file = "outputs/logcpm_April2023.txt", sep = "\t", quote = F, row.names = F)
write.table(CICDUX4.EV_results, file = "outputs/CICDUX4vsEV_April2023.txt", sep = "\t", quote = F, row.names = F)
write.table(CICex20_NUTM1ex6.EV_results, file = "outputs/CICex20NUTM1ex6vsEV_April2023.txt", sep = "\t", quote = F, row.names = F)
write.table(tibble::rownames_to_column(master, "Gene"), file = "outputs/counts_April2023.txt", sep = "\t", quote = F, row.names = F)




