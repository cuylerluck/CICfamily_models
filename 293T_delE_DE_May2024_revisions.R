#Written by Cuyler Luck
#Contact: cuyler.luck@ucsf.edu / cuylerluck@gmail.com or ross.okimoto@ucsf.edu

#This analysis is on stranded RNA-seq data from 293T cells ~48h after transient transfection with one of seven conditions:
# A = empty vector (same backbone as CIC-NUTM1 plasmids)
# B = HA-CIC::DUX4 (different backbone than others, but same promoter)
# C = HA-CIC(ex20)::NUTM1(ex6)
# D = HA-CIC(ex20)::NUTM1(ex6) dTAD
# E = HA-CIC(ex20)::NUTM1(ex6) delE
# F = HA-CIC(ex18)::NUTM1(ex3)
# G = HA-CIC::LEUTX

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
library(UpSetR) #1.4.0

setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_293T_delE_CICNUTM1") #this will change by user & location of data

#Next read in all data
#Skip the first four lines, they just give info on odd-mapping statistics
A1 = fread(file="Counts/A1ReadsPerGene.out.tab", skip = 4)
A2 = fread(file="Counts/A2ReadsPerGene.out.tab", skip = 4)
A3 = fread(file="Counts/A3ReadsPerGene.out.tab", skip = 4)
B1 = fread(file="Counts/B1ReadsPerGene.out.tab", skip = 4)
B2 = fread(file="Counts/B2ReadsPerGene.out.tab", skip = 4)
B3 = fread(file="Counts/B3ReadsPerGene.out.tab", skip = 4)
C1 = fread(file="Counts/C1ReadsPerGene.out.tab", skip = 4)
C2 = fread(file="Counts/C2ReadsPerGene.out.tab", skip = 4)
C3 = fread(file="Counts/C3ReadsPerGene.out.tab", skip = 4)
D1 = fread(file="Counts/D1ReadsPerGene.out.tab", skip = 4)
D2 = fread(file="Counts/D2ReadsPerGene.out.tab", skip = 4)
D3 = fread(file="Counts/D3ReadsPerGene.out.tab", skip = 4)
E1 = fread(file="Counts/E1ReadsPerGene.out.tab", skip = 4)
E2 = fread(file="Counts/E2ReadsPerGene.out.tab", skip = 4)
E3 = fread(file="Counts/E3ReadsPerGene.out.tab", skip = 4)
F1 = fread(file="Counts/F1ReadsPerGene.out.tab", skip = 4)
F2 = fread(file="Counts/F2ReadsPerGene.out.tab", skip = 4)
F3 = fread(file="Counts/F3ReadsPerGene.out.tab", skip = 4)
G1 = fread(file="Counts/G1ReadsPerGene.out.tab", skip = 4)
G2 = fread(file="Counts/G2ReadsPerGene.out.tab", skip = 4)
G3 = fread(file="Counts/G3ReadsPerGene.out.tab", skip = 4)

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
D1_strand = dplyr::select(D1, c(1,4))
D2_strand = dplyr::select(D2, c(1,4))
D3_strand = dplyr::select(D3, c(1,4))
E1_strand = dplyr::select(E1, c(1,4))
E2_strand = dplyr::select(E2, c(1,4))
E3_strand = dplyr::select(E3, c(1,4))
F1_strand = dplyr::select(F1, c(1,4))
F2_strand = dplyr::select(F2, c(1,4))
F3_strand = dplyr::select(F3, c(1,4))
G1_strand = dplyr::select(G1, c(1,4))
G2_strand = dplyr::select(G2, c(1,4))
G3_strand = dplyr::select(G3, c(1,4))

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
colnames(D1_strand) = c("ENSG", "D1")
colnames(D2_strand) = c("ENSG", "D2")
colnames(D3_strand) = c("ENSG", "D3")
colnames(E1_strand) = c("ENSG", "E1")
colnames(E2_strand) = c("ENSG", "E2")
colnames(E3_strand) = c("ENSG", "E3")
colnames(F1_strand) = c("ENSG", "F1")
colnames(F2_strand) = c("ENSG", "F2")
colnames(F3_strand) = c("ENSG", "F3")
colnames(G1_strand) = c("ENSG", "G1")
colnames(G2_strand) = c("ENSG", "G2")
colnames(G3_strand) = c("ENSG", "G3")


master = inner_join(A1_strand, A2_strand, by = "ENSG")
master = inner_join(master, A3_strand, by = "ENSG")
master = inner_join(master, B1_strand, by = "ENSG")
master = inner_join(master, B2_strand, by = "ENSG")
master = inner_join(master, B3_strand, by = "ENSG")
master = inner_join(master, C1_strand, by = "ENSG")
master = inner_join(master, C2_strand, by = "ENSG")
master = inner_join(master, C3_strand, by = "ENSG")
master = inner_join(master, D1_strand, by = "ENSG")
master = inner_join(master, D2_strand, by = "ENSG")
master = inner_join(master, D3_strand, by = "ENSG")
master = inner_join(master, E1_strand, by = "ENSG")
master = inner_join(master, E2_strand, by = "ENSG")
master = inner_join(master, E3_strand, by = "ENSG")
master = inner_join(master, F1_strand, by = "ENSG")
master = inner_join(master, F2_strand, by = "ENSG")
master = inner_join(master, F3_strand, by = "ENSG")
master = inner_join(master, G1_strand, by = "ENSG")
master = inner_join(master, G2_strand, by = "ENSG")
master = inner_join(master, G3_strand, by = "ENSG")

#move ENSG to rownames for edgeR syntax
master = tibble::column_to_rownames(master, "ENSG")

#Let's remove CIC, NUTM1, DUX4, and LEUTX just to make sure the overexpression of fusions doesn't change any clustering.
master = master[!(rownames(master) %in% c("ENSG00000079432", "ENSG00000260596", "ENSG00000184507", "ENSG00000213921")),]

#We can now pass this master data frame into an edgeR pipeline
#I am choosing to use a GLM for this analysis. So I will use the non-classic pipeline.

groups = c("EV", "EV", "EV", 
           "CICDUX4", "CICDUX4", "CICDUX4", 
           "CICex20_NUTM1ex6", "CICex20_NUTM1ex6", "CICex20_NUTM1ex6",
           "dTAD", "dTAD", "dTAD",
           "delE", "delE", "delE",
           "CICex18_NUTM1ex3", "CICex18_NUTM1ex3", "CICex18_NUTM1ex3",
           "CICLEUTX", "CICLEUTX", "CICLEUTX")
dg = DGEList(counts=master, group = groups)
keep = filterByExpr(dg)
dg = dg[keep, , keep.lib.sizes = FALSE]
dg = calcNormFactors(dg)
#reordering the levels of the group factor so "EV" becomes the reference
dg$samples$group = relevel(dg$samples$group, ref="EV")
#checking the order with
#dg$samples$group
#reveals that the present levels are:
#EV, CICDUX4, Cex18Nex3, Cex20Nex6, CICLEUTX, delE, dTAD
design = model.matrix(~dg$samples$group, data = dg$samples)
#checking the design matrix shows that now the intercept (baseline condition) is EV.
#coefficient 2 = CICDUX4
#coefficient 3 = CICex18-NUTM1ex3
#coefficient 4 = CICex20-NUTM1ex6
#coefficient 5 = CIC-LEUTX
#coefficient 6 = delE
#coefficient 7 = dTAD
#this information will be useful in doing pairwise comparisons soon.
dg = estimateDisp(dg, design)
#plotBCV(dg) #to see what the BCV looks like if desired

#now we can do pairwise comparisons between groups of interest
#with this design matrix, providing glmQLFTest a number for coef means "compare this coefficient to the baseline [EV]"
#first we get the fit
fit = glmQLFit(dg, design)

#compare CICDUX4 to EV
qlf.CICDUX4.EV = glmQLFTest(fit, coef= 2)

#compare CICex20_NUTM1ex6 to EV
qlf.CICex20_NUTM1ex6.EV = glmQLFTest(fit, coef= 4)

#compare dTAD to EV
qlf.dTAD.EV = glmQLFTest(fit, coef= 7)

#compare delE to EV
qlf.delE.EV = glmQLFTest(fit, coef= 6)

#compare CICex18_NUTM1ex3 to EV
qlf.CICex18_NUTM1ex3.EV = glmQLFTest(fit, coef= 3)

#compare CICLEUTX to EV
qlf.CICLEUTX.EV = glmQLFTest(fit, coef= 5)


#compare delE to CICex20_NUTM1ex6
qlf.delE.CICex20_NUTM1ex6 = glmQLFTest(fit, contrast=c(0, 0, 0, -1, 0, 1, 0)) #i.e. 6 coefficient with 4 coefficient

#compare dTAD to CICex20_NUTM1ex6
qlf.dTAD.CICex20_NUTM1ex6 = glmQLFTest(fit, contrast=c(0, 0, 0, -1, 0, 0, 1)) #i.e. 7 coefficient with 4 coefficient

#compare delE to dTAD
qlf.delE.dTAD = glmQLFTest(fit, contrast=c(0, 0, 0, 0, 0, 1, -1)) #i.e. 6 coefficient with 7 coefficient



#now we can pull out tables of results from these tests for use in volcano plots etc.
CICDUX4.EV_results = qlf.CICDUX4.EV$table

CICex20_NUTM1ex6.EV_results = qlf.CICex20_NUTM1ex6.EV$table

dTAD.EV_results = qlf.dTAD.EV$table

delE.EV_results = qlf.delE.EV$table

CICex18_NUTM1ex3.EV_results = qlf.CICex18_NUTM1ex3.EV$table

CICLEUTX.EV_results = qlf.CICLEUTX.EV$table

delE.CICex20_NUTM1ex6_results = qlf.delE.CICex20_NUTM1ex6$table

dTAD.CICex20_NUTM1ex6_results = qlf.dTAD.CICex20_NUTM1ex6$table

delE.dTAD_results = qlf.delE.dTAD$table


#Let's add FDR-adjusted p-values to each of these
CICDUX4.EV_results = dplyr::mutate(CICDUX4.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

CICex20_NUTM1ex6.EV_results = dplyr::mutate(CICex20_NUTM1ex6.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

dTAD.EV_results = dplyr::mutate(dTAD.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

delE.EV_results = dplyr::mutate(delE.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

CICex18_NUTM1ex3.EV_results = dplyr::mutate(CICex18_NUTM1ex3.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

CICLEUTX.EV_results = dplyr::mutate(CICLEUTX.EV_results, fdr_adj = p.adjust(PValue, method = "fdr"))

delE.CICex20_NUTM1ex6_results = dplyr::mutate(delE.CICex20_NUTM1ex6_results, fdr_adj = p.adjust(PValue, method = "fdr"))

dTAD.CICex20_NUTM1ex6_results = dplyr::mutate(dTAD.CICex20_NUTM1ex6_results, fdr_adj = p.adjust(PValue, method = "fdr"))

delE.dTAD_results = dplyr::mutate(delE.dTAD_results, fdr_adj = p.adjust(PValue, method = "fdr"))


#let's also pull out the TMM-normalized log(cpm) values from the edgeR pipeline
#this will be useful later for looking at expression of specific genes
logcpm = as.data.frame(cpm(dg, log = T))


#and we can ask edgeR to plot an MDS plot to see how well samples cluster

#set color palette for all plots in same order as the samples are in master
all_palette = rep(c("#a9a9a9", "#fc558e", "#31a851", "#62f5e6", "#1e4291","#03611c", "#e2bc29"), each = 3)


pdf(file = "plots/plotMDS_all_symbols.pdf", width = 6, height = 4)
par(mar=c(6.1,5.1,5.1,11.1))
plotMDS(dg, col = all_palette,
        pch = rep(c(15, 16, 17, 18, 8, 3, 4), each = 3), 
        cex = 1.5)
legend("topright", inset = c(-.80,0), 
       legend = c("EV", "HA-CIC::DUX4", "HA-CICex20::NUTM1ex6", "dTAD", "delE", "HA-CICex18::NUTM1ex3", "HA-CIC::LEUTX"), 
       pch = c(15, 16, 17, 18, 8, 3, 4),
       col = all_palette[c(1,4,7,10,13,16,19)],
       xpd = T,
       bty = "n")

dev.off()


pdf(file = "plots/plotMDS_all_samples.pdf", width = 5, height = 4)
plotMDS(dg, col = all_palette, 
        labels = rep(c("EV", "HCD", "HC20N6", "dTAD", "delE", "HC18N3", "HCL"), each = 3))
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

dTAD.EV_results = tibble::rownames_to_column(dTAD.EV_results, "ensembl_gene_id")
dTAD.EV_results = left_join(dTAD.EV_results, symbol_key, by = "ensembl_gene_id")

delE.EV_results = tibble::rownames_to_column(delE.EV_results, "ensembl_gene_id")
delE.EV_results = left_join(delE.EV_results, symbol_key, by = "ensembl_gene_id")

CICex18_NUTM1ex3.EV_results = tibble::rownames_to_column(CICex18_NUTM1ex3.EV_results, "ensembl_gene_id")
CICex18_NUTM1ex3.EV_results = left_join(CICex18_NUTM1ex3.EV_results, symbol_key, by = "ensembl_gene_id")

CICLEUTX.EV_results = tibble::rownames_to_column(CICLEUTX.EV_results, "ensembl_gene_id")
CICLEUTX.EV_results = left_join(CICLEUTX.EV_results, symbol_key, by = "ensembl_gene_id")

delE.CICex20_NUTM1ex6_results = tibble::rownames_to_column(delE.CICex20_NUTM1ex6_results, "ensembl_gene_id")
delE.CICex20_NUTM1ex6_results = left_join(delE.CICex20_NUTM1ex6_results, symbol_key, by = "ensembl_gene_id")

dTAD.CICex20_NUTM1ex6_results = tibble::rownames_to_column(dTAD.CICex20_NUTM1ex6_results, "ensembl_gene_id")
dTAD.CICex20_NUTM1ex6_results = left_join(dTAD.CICex20_NUTM1ex6_results, symbol_key, by = "ensembl_gene_id")

delE.dTAD_results = tibble::rownames_to_column(delE.dTAD_results, "ensembl_gene_id")
delE.dTAD_results = left_join(delE.dTAD_results, symbol_key, by = "ensembl_gene_id")

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
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICDUX4vsEVplot = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$logFC > 2 & -log10(CICDUX4.EV_results$fdr_adj) > 5,], 
             color = all_palette[4], size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4","FOXB1"),],
    aes(label = hgnc_symbol), nudge_y = 3, nudge_x = -1, box.padding = 0.5, size = 5 ) +
  ggtitle("CICDUX4 vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CD4_vs_EV.pdf",width = 8, height = 6)
CICDUX4vsEVplot
dev.off()

#this plots CN73 vs. EV results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICex20_NUTM1ex6vsEVplot = ggplot(data = CICex20_NUTM1ex6.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$logFC > 2 & -log10(CICex20_NUTM1ex6.EV_results$fdr_adj) > 5,], 
             color = all_palette[7], size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4", "FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 3, nudge_x = -1, box.padding = 0.5, size = 5) +
  ggtitle("CICex20_NUTM1ex6 vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CN73_vs_EV.pdf",width = 8, height = 6)
CICex20_NUTM1ex6vsEVplot
dev.off()



#this plots CN84 vs. EV results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2


CICex18_NUTM1ex3vsEVplot = ggplot(data = CICex18_NUTM1ex3.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex18_NUTM1ex3.EV_results[CICex18_NUTM1ex3.EV_results$logFC > 2 & -log10(CICex18_NUTM1ex3.EV_results$fdr_adj) > 5,], 
             color = all_palette[16], size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICex18_NUTM1ex3.EV_results[CICex18_NUTM1ex3.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4", "FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 3, nudge_x = -1, box.padding = 0.5, size= 5) +
  ggtitle("CICex18_NUTM1ex3 vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CN84_vs_EV.pdf",width = 8, height = 6)
CICex18_NUTM1ex3vsEVplot
dev.off()


#this plots CIC-LEUTX vs. EV results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICLEUTXvsEVplot = ggplot(data = CICLEUTX.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICLEUTX.EV_results[CICLEUTX.EV_results$logFC > 2 & -log10(CICLEUTX.EV_results$fdr_adj) > 5,], 
             color = all_palette[19], size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICLEUTX.EV_results[CICLEUTX.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4", "FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 0, nudge_x = 0, box.padding = 0.5, size = 5) +
  ggtitle("CICLEUTX vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CL_vs_EV.pdf",width = 8, height = 6)
CICLEUTXvsEVplot
dev.off()



#this plots delE vs. dTAD results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

delEvsdTADplot = ggplot(data = delE.dTAD_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = delE.dTAD_results[delE.dTAD_results$logFC > 2 & -log10(delE.dTAD_results$fdr_adj) > 5,], 
             color = all_palette[13], size = 4) +
  geom_point(data = delE.dTAD_results[delE.dTAD_results$logFC < -2 & -log10(delE.dTAD_results$fdr_adj) > 5,], 
             color = all_palette[10], size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = delE.dTAD_results[delE.dTAD_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4", "FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 4, nudge_x = -5, box.padding = 0.5, size = 6) +
  ggtitle("delE vs. dTAD") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/delE_vs_dTAD.pdf",width = 8, height = 6)
delEvsdTADplot
dev.off()


#this plots delE vs. EV results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2


delEvsEVplot = ggplot(data = delE.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = delE.EV_results[delE.EV_results$logFC > 2 & -log10(delE.EV_results$fdr_adj) > 5,], 
             color = all_palette[13], size = 4) +
  geom_point(data = delE.EV_results[abs(delE.EV_results$logFC) <= 1 | -log10(delE.EV_results$fdr_adj) <= 2,], 
             color = "dark grey") +
  theme_classic() +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = delE.EV_results[delE.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4","FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 1.5, nudge_x = -1, box.padding = 0.5) +
  ggtitle("delE vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/delE_vs_EV.pdf",width = 8, height = 6)
delEvsEVplot
dev.off()

#this plots dTAD vs. EV results.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2


dTADvsEVplot = ggplot(data = dTAD.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = dTAD.EV_results[dTAD.EV_results$logFC > 2 & -log10(dTAD.EV_results$fdr_adj) > 5,], 
             color = all_palette[10], size = 4) +
  theme_classic() +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = dTAD.EV_results[dTAD.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5","VGF","DUSP6","SPRY4","FOXB1"),],
                  aes(label = hgnc_symbol), nudge_y = 1.5, nudge_x = -1, box.padding = 0.5) +
  ggtitle("dTAD vs. EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/dTAD_vs_EV.pdf",width = 8, height = 6)
dTADvsEVplot
dev.off()


### Reviewers asked for a more thorough analysis of where previously nominated CIC::DUX4 gene targets
# align in this transient 293T system. An easy way to evaluate this is by highlighting them on a volcano
# plot, for the CIC::DUX4 vs EV comparison.

# There are two gene lists I will overlay.

# First, from Watson et al bioRxiv 2019: https://www.biorxiv.org/content/10.1101/517722v2
# this list aggregates genes that were identified iin separate CDS experiments evaluating 
# fish tumors, human tumors, and CIC::DUX4 knockdown in a CDS cell line (IB120)

watson_list = c("VGF", "ETV4", "LBH", "SCG3", "SPRY4", "DUSP4", "ETV5", "SHC3", "ETV1", "HMGA2",
                "CRH", "SPRED1", "ZNF423", "PCDH18", "MAP2", "ETS1", "HMGB2", "KIF11", "RAD51",
                "CDH4", "TGFB3", "PSAT1", "NRARP", "PRKAR2B", "ANGPT2")

#how many of these are also in the results from the CIC::DUX4 vs EV comparison?

watson_list_retained = watson_list[watson_list %in% CICDUX4.EV_results$hgnc_symbol]

#all but one - CRH is not there. 24 are.

#volcano plot with the watson genes overlaid

#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICDUX4vsEVplot_watson = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% watson_list_retained,], color = "orange", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_label_repel(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% watson_list_retained,],
                  aes(label = hgnc_symbol), box.padding = 0.5) +
  ggtitle("CICDUX4 vs. EV, 24 Watson genes in orange") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CD4_vs_EV_watson.pdf",width = 8, height = 6)
CICDUX4vsEVplot_watson
dev.off()


# Similar but using genes from Ringnalda et al, Nat Comm, 2025: nature.com/articles/s41467-025-62673-2
# In particular, I am using the genes from Fig 4 B and C which are in the block of genes specific to
# CDS tumoroids over other tumoroids.
# excluding CIC and DUX4 because they are excluded from my RNA-seq data.

ringnalda = c("WT1", "DUSP4", "ETS1", "ETV1", "ETV4","ETV5", "MCL1", "MUC5AC","DUSP6", "HMGA2", "SPRED1",
              "SPRY4", "VGF", "POLE", "SHC3", "VGF", "NME3", "SPON2", "PTPN9", 
              "LINC01299", "SPACA5", "FOXL1", "ACSF2", "IRS2", "PER2", "KLF9", "SHC4", "LINC00865",
              "LRRC4C", "FOXN3", "SCNN1G", "COLEC11", "TMEM132E", "IL31RA")

#which are in my data too?

ringnalda_retained = ringnalda[ringnalda %in% CICDUX4.EV_results$hgnc_symbol]

#most - 30. Doesn't include SPACA5, LINC00865, LRRC4C, or IL31RA.

#Overlay on volcano plot again

#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICDUX4vsEVplot_ringnalda = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% ringnalda_retained,], color = "turquoise", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_label_repel(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% ringnalda_retained,],
                   aes(label = hgnc_symbol), box.padding = 0.2) +
  ggtitle("CICDUX4 vs. EV, 30 Ringnalda genes in turquoise") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CD4_vs_EV_ringnalda.pdf",width = 8, height = 6)
CICDUX4vsEVplot_ringnalda
dev.off()




### Let's make some heatmaps too.

#Let's use all genes significantly upregulated (log2FC > 2 and q < 0.00001) in any of the four full-length fusions
#i.e. CIC-DUX4, CICex20-NUTM1ex6, CICex18-NUTM1ex3, and CIC-LEUTX
#note that I am choosing to only use genes with HGNC symbols

up_CICDUX4_vs_EV = CICDUX4.EV_results[CICDUX4.EV_results$logFC > 2 & CICDUX4.EV_results$fdr_adj < 0.00001 & !is.na(CICDUX4.EV_results$hgnc_symbol),]$hgnc_symbol
up_CICex20_NUTM1ex6_vs_EV = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$logFC > 2 & CICex20_NUTM1ex6.EV_results$fdr_adj < 0.00001 & !is.na(CICex20_NUTM1ex6.EV_results$hgnc_symbol),]$hgnc_symbol
up_CICex18_NUTM1ex3_vs_EV = CICex18_NUTM1ex3.EV_results[CICex18_NUTM1ex3.EV_results$logFC > 2 & CICex18_NUTM1ex3.EV_results$fdr_adj < 0.00001 & !is.na(CICex18_NUTM1ex3.EV_results$hgnc_symbol),]$hgnc_symbol
up_CICLEUTX_vs_EV = CICLEUTX.EV_results[CICLEUTX.EV_results$logFC > 2 & CICLEUTX.EV_results$fdr_adj < 0.00001 & !is.na(CICLEUTX.EV_results$hgnc_symbol), ]$hgnc_symbol

#take the union of all of these gene lists
up_any_full_fusion = union(c(c(up_CICDUX4_vs_EV, up_CICex20_NUTM1ex6_vs_EV), up_CICex18_NUTM1ex3_vs_EV), up_CICLEUTX_vs_EV)

#now we can shorten logcpm to just these genes to make the heatmap more manageable.
logcpm_up_any = logcpm[up_any_full_fusion,]

#Plot a heatmap for all samples and all genes significantly up in at least one full length fusion

#Set column annotations and colors

sample_annotations = data.frame(condition = rep(c("EV", "HA-CIC::DUX4", "HA-CICex20::NUTM1ex6", "dTAD", "delE", "HA-CICex18::NUTM1ex3", "HA-CIC::LEUTX"), each = 3))
rownames(sample_annotations) = colnames(logcpm)

heatmap_colors = list(
  condition = all_palette
)

names(heatmap_colors$condition) = sample_annotations$condition

pdf(file = "plots/heatmap_all.pdf", width = 8, height = 10)
pheatmap(logcpm_up_any, 
         cluster_cols = F, 
         scale = "row", 
         show_rownames = F,
         angle_col = 45,
         annotation_col = sample_annotations,
         annotation_colors = heatmap_colors)
dev.off()



#It's a little confusing to do all of them together. Let's do some in smaller combinations

#Just CIC-DUX4 and CICex20_NUTM1ex6 -- this should show if the initial findings are reproducible
up_DUX4_CN73 = union(up_CICDUX4_vs_EV, up_CICex20_NUTM1ex6_vs_EV)

logcpm_up_DUX4_CN73 = logcpm[up_DUX4_CN73,]
logcpm_up_DUX4_CN73_only_full_fusions = logcpm_up_DUX4_CN73[,-c(10:21)] #removes dTAD, delE, CICex18_NUTM1ex3, and CIC-LEUTX

#run pheatmap with just EV, CIC-DUX4, and CICex20-NUTM1ex6
pdf(file = "plots/de_EV_CD4_CN73.pdf", width = 8, height = 10)
pheatmap(logcpm_up_DUX4_CN73_only_full_fusions, 
                                                  cluster_cols = F, 
                                                  scale = "row", 
                                                  show_rownames = F,
                                                  angle_col = 45,
                                                  annotation_col = sample_annotations,
                                                  annotation_colors = heatmap_colors)
dev.off()





#Just CICex20_NUTM1ex6 and delE
logcpm_up_CN73_delE = logcpm[up_CICex20_NUTM1ex6_vs_EV,]
logcpm_up_CN73_delE = logcpm_up_CN73_delE[,c(1:3,7:9, 13:15)]

#run pheatmap with just EV, CICex20-NUTM1ex6, and delE
pdf(file = "plots/de_EV_CN73_delE.pdf", width = 8, height = 8)
pheatmap(logcpm_up_CN73_delE, 
         cluster_cols = F, 
         scale = "row", 
         show_rownames = F,
         angle_col = 45,
         annotation_col = sample_annotations,
         annotation_colors = heatmap_colors)
dev.off()

#To better see which genes are potentially turned on without the delE region,
#let's see what genes are sig. up in BOTH HCN73 vs EV AND delE vs EV.
#sig. on = log2FC > 2, q < 0.00001, and has an HGNC symbol. 

#first, genes sig. up in delE vs EV
up_delE_vs_EV = delE.EV_results[delE.EV_results$logFC > 2 & delE.EV_results$fdr_adj < 0.00001 & !is.na(delE.EV_results$hgnc_symbol),]$hgnc_symbol

#then, intersect with genes sig. up in HCN73 vs EV
delE_ignorers = intersect(up_delE_vs_EV, up_CICex20_NUTM1ex6_vs_EV)

#Let's make some annotations for the mutant heatmap to display which genes fall into this bucket
mutant_anno_row = data.frame(gene = row.names(logcpm_up_CN73_delE))
mutant_anno_row = dplyr::mutate(mutant_anno_row, delE_ignore = gene %in% delE_ignorers)
mutant_anno_row = tibble::column_to_rownames(mutant_anno_row, "gene")
#change FALSE --> NO and TRUE --> YES for annotation purposes
mutant_anno_row$delE_ignore[mutant_anno_row$delE_ignore] = "YES"
mutant_anno_row$delE_ignore[mutant_anno_row$delE_ignore == "FALSE"] = "NO"


#similarly, we can narrow down genes that ARE delE responsive by getting those with the following characteristics:
#sig. on in HCN73 vs EV
#not sig. different, using stricter criteria, in delE vs EV
#sig. on = log2FC > 2, q < 0.00001, and has an HGNC symbol. 
#not sig different in either direction = |log2FC| <= 1, q >= 0.01, has an HGNC symbol

#first, not sig different between delE and EV
not_sig_delE_vs_EV = delE.EV_results[(abs(delE.EV_results$logFC) <= 1 | delE.EV_results$fdr_adj >= 0.01) & !is.na(delE.EV_results$hgnc_symbol),]$hgnc_symbol

#intersect with genes up in HCN73 vs EV
delE_responders = intersect(not_sig_delE_vs_EV, up_CICex20_NUTM1ex6_vs_EV)

#add to the annotation DF and plot heatmap again from above
mutant_anno_row = dplyr::mutate(mutant_anno_row, delE_respond = rownames(mutant_anno_row) %in% delE_responders)
mutant_anno_row$delE_respond[mutant_anno_row$delE_respond] = "YES"
mutant_anno_row$delE_respond[mutant_anno_row$delE_respond == "FALSE"] = "NO"


#make factors
mutant_anno_row$delE_ignore = as.factor(mutant_anno_row$delE_ignore)
mutant_anno_row$delE_respond = as.factor(mutant_anno_row$delE_respond)


#choose colors for each annotation type

heatmap_colors$delE_ignore = c(NO="light gray",YES="dark red")
heatmap_colors$delE_respond = c(NO="light gray", YES="navy blue")


pdf(file = "plots/de_EV_CN73_delE_ignorers_responders.pdf", width = 8, height = 8)
pheatmap(logcpm_up_CN73_delE, 
                  cluster_cols = F, 
                  scale = "row", 
                  show_rownames = F,
                  annotation_row = mutant_anno_row, 
                  annotation_colors = heatmap_colors,
                  annotation_col = sample_annotations,
                  angle_col = 45)
dev.off()


#let's write the delE_ignorers and delE_responders to text files for GO or other analysis
write.table(delE_ignorers, file = "outputs/delE_ignorers.txt", quote = F, row.names = F, col.names = F)
write.table(delE_responders, file = "outputs/delE_responders.txt", quote = F, row.names = F, col.names = F)


#what does it look like if we take just the delE ignorers + responders and do the same heatmap of EV, CICDUX4, and CICex20-NUTM1ex6 
#as above, but with them labeled?

#then again make an annotation df for this
EV_CD4_CN73_anno_row = data.frame(gene = row.names(logcpm_up_DUX4_CN73_only_full_fusions))
EV_CD4_CN73_anno_row = dplyr::mutate(EV_CD4_CN73_anno_row, delE_ignore = gene %in% delE_ignorers)
EV_CD4_CN73_anno_row = tibble::column_to_rownames(EV_CD4_CN73_anno_row, "gene")
#change FALSE --> NO and TRUE --> YES for annotation purposes
EV_CD4_CN73_anno_row$delE_ignore[EV_CD4_CN73_anno_row$delE_ignore] = "YES"
EV_CD4_CN73_anno_row$delE_ignore[EV_CD4_CN73_anno_row$delE_ignore == "FALSE"] = "NO"
#same for delE_responders
EV_CD4_CN73_anno_row = dplyr::mutate(EV_CD4_CN73_anno_row, delE_respond = rownames(EV_CD4_CN73_anno_row) %in% delE_responders)
EV_CD4_CN73_anno_row$delE_respond[EV_CD4_CN73_anno_row$delE_respond] = "YES"
EV_CD4_CN73_anno_row$delE_respond[EV_CD4_CN73_anno_row$delE_respond == "FALSE"] = "NO"


#make factors
EV_CD4_CN73_anno_row$delE_ignore = as.factor(EV_CD4_CN73_anno_row$delE_ignore)
EV_CD4_CN73_anno_row$delE_respond = as.factor(EV_CD4_CN73_anno_row$delE_respond)

#using same color annotation df is fine.


pdf(file = "plots/delE_EV_HCD_HCN73.pdf", width = 8, height = 8)
delE_heatmap = pheatmap(logcpm_up_DUX4_CN73_only_full_fusions,
         cluster_cols = F,
         cluster_rows = T,
         scale = "row",
         show_rownames = F,
         angle_col = 45,
         annotation_row = EV_CD4_CN73_anno_row,
         annotation_colors = heatmap_colors,
         annotation_col = sample_annotations)
dev.off()





#Some GO analysis for the delE_responders and delE_ignorers using gProfiler2:

delE_responders_gprofiler = gost(delE_responders, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
delE_responders_GO = delE_responders_gprofiler$result
#can do interactive plot if desired
#gostplot(delE_responders_gprofiler, interactive = T)
#dev.off()

delE_ignorers_gprofiler = gost(delE_ignorers, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
delE_ignorers_GO = delE_ignorers_gprofiler$result
#can do interactive plot if desired
#gostplot(delE_ignorers_gprofiler, interactive = T)
#dev.off()


multi_delE_gprofiler = gost(list("delE_responders" = delE_responders, "delE_ignorers" = delE_ignorers),
                            organism = "hsapiens",
                            correction_method = "gSCS",
                            multi_query = T,
                            sources = "GO:BP")
multi_delE_GO = multi_delE_gprofiler$result
p1 = gostplot(multi_delE_gprofiler, interactive = F)
pdf(file = "plots/delE_responders_ignorers_gprofiler2.pdf", width = 10, height = 8)
publish_gostplot(p1, highlight_terms = c("GO:0032502",
                                         "GO:0048856",
                                         "GO:0070372",
                                         "GO:0007399", 
                                         "GO:0030182",
                                         "GO:0048880",
                                         "GO:0009790"))
dev.off()
  

#Now that we have lists of delE_ignorers and delE_responders, I'd like to look at where they fall on
#some of the volcano plots

#this plots CIC-DUX4 vs. EV results, navy blue for delE responders
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICDUX4vsEVplot_delE_responders = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% delE_responders,], 
             color = "navy blue", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICDUX4 vs. EV, delE_responders") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

#this plots CIC-DUX4 vs. EV results, dark red for delE ignorers
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICDUX4vsEVplot_delE_ignorers = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% delE_ignorers,], 
             color = "dark red", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICDUX4 vs. EV, delE_ignorers") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CD4_vs_EV_delE.pdf",width = 8, height = 10)
CICDUX4vsEVplot_delE_responders / CICDUX4vsEVplot_delE_ignorers
dev.off()


#this plots CN73 vs. EV results, navy blue for delE_responders
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICex20_NUTM1ex6vsEVplot_delE_responders = ggplot(data = CICex20_NUTM1ex6.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$hgnc_symbol %in% delE_responders,], 
             color = "navy blue", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICex20_NUTM1ex6 vs. EV, delE_responders") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

#this plots CN73 vs. EV results, dark red for delE_ignorers
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICex20_NUTM1ex6vsEVplot_delE_ignorers = ggplot(data = CICex20_NUTM1ex6.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex20_NUTM1ex6.EV_results[CICex20_NUTM1ex6.EV_results$hgnc_symbol %in% delE_ignorers,], 
             color = "dark red", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICex20_NUTM1ex6 vs. EV, delE_ignorers") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CN73_vs_EV_delE.pdf",width = 8, height = 10)
CICex20_NUTM1ex6vsEVplot_delE_responders / CICex20_NUTM1ex6vsEVplot_delE_ignorers
dev.off()




#this plots CN84 vs. EV results, navy blue for delE_responders.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2


CICex18_NUTM1ex3vsEVplot_delE_responders = ggplot(data = CICex18_NUTM1ex3.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex18_NUTM1ex3.EV_results[CICex18_NUTM1ex3.EV_results$hgnc_symbol %in% delE_responders,], 
             color = "navy blue", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICex18_NUTM1ex3 vs. EV, delE_responders") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

#this plots CN84 vs. EV results, dark red for delE_ignorers.
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2


CICex18_NUTM1ex3vsEVplot_delE_ignorers = ggplot(data = CICex18_NUTM1ex3.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICex18_NUTM1ex3.EV_results[CICex18_NUTM1ex3.EV_results$hgnc_symbol %in% delE_ignorers,], 
             color = "dark red", size = 4) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  ggtitle("CICex18_NUTM1ex3 vs. EV, delE_ignorers") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))


pdf(file = "plots/CN84_vs_EV_delE.pdf",width = 8, height = 10)
CICex18_NUTM1ex3vsEVplot_delE_responders / CICex18_NUTM1ex3vsEVplot_delE_ignorers
dev.off()



#what proportion of delE_responders and delE_ignorers are still sig. up for these fusions?
CN73_prop_delE_responders_up = length(up_CICex20_NUTM1ex6_vs_EV[up_CICex20_NUTM1ex6_vs_EV %in% delE_responders]) / length(delE_responders)
CN73_prop_delE_ignorers_up = length(up_CICex20_NUTM1ex6_vs_EV[up_CICex20_NUTM1ex6_vs_EV %in% delE_ignorers]) / length(delE_ignorers)
#for HA-CICex20::NUTM1ex6 these are 1.00 by definition (those genes were chosen with being up in HCN73 vs EV as a criteria)

CICDUX4_prop_delE_responders_up = length(up_CICDUX4_vs_EV[up_CICDUX4_vs_EV %in% delE_responders]) / length(delE_responders)
CICDUX4_prop_delE_ignorers_up = length(up_CICDUX4_vs_EV[up_CICDUX4_vs_EV %in% delE_ignorers]) / length(delE_ignorers)
#So about 90% of delE_ignorers are up in HCD vs EV, but only 14% of delE_responders!

CN84_prop_delE_responders_up = length(up_CICex18_NUTM1ex3_vs_EV[up_CICex18_NUTM1ex3_vs_EV %in% delE_responders]) / length(delE_responders)
CN84_prop_delE_ignorers_up = length(up_CICex18_NUTM1ex3_vs_EV[up_CICex18_NUTM1ex3_vs_EV %in% delE_ignorers]) / length(delE_ignorers)
#For HA-CICex18::NUTM1ex3, it's 35% of ignorers, 24% of responders. Consistent with both groups being affected similarly
#i.e. that HCN84 is overall not as good at turning these genes on with this specific context and these specific thresholds


#The next informative thing will be plotting log2(cpm) values and comparing across conditions.

logcpm_pivot = tibble::rownames_to_column(logcpm, "gene")
logcpm_pivot = pivot_longer(logcpm_pivot, cols = c(2:22), names_to = "sample", values_to = "log2cpm")

#add a "condition" variable to group each sample by
logcpm_pivot = mutate(logcpm_pivot, condition = ifelse(sample %in% c("A1","A2","A3"), "EV", 
                                                      ifelse(sample %in% c("C1", "C2", "C3"), "HCN73",
                                                             ifelse(sample %in% c("D1","D2","D3"), "dTAD",
                                                                    ifelse(sample %in% c("E1","E2","E3"), "delE", 
                                                                           ifelse(sample %in% c("F1","F2","F3"), "HCN84", 
                                                                                  ifelse(sample %in% c("B1","B2","B3"), "HCD", "HCL")))))))

#make a simplified condition-color key to use in a second
condition_colors = heatmap_colors$condition[c(1,4,7,10,13,16,19)]
names(condition_colors) = c("EV", "HCD", "HCN73", "dTAD", "delE", "HCN84", "HCL")

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

  
#let's make a bunch of these for delE_ignorers and delE_responders, and compare how they look with each other
#including dTAD as a general impaired control, and HCD as an indicator of if the gene is HCN73-skewed or not
selected_delE_ignorers = c("DUSP4","DUSP6","ETV1","ETV4","ETV5","SHC3","SPRED1","SPRED2","SPRED3","SPRY3","SPRY4","VGF")
selected_delE_responders = c("DLL1","FLRT3","FOXB1","FOXD3","HMX3","JAG1","NEUROD4","NPAS3","NR4A2","NR4A3","PAX5","RUNX2","SHOX","SOX2")

plots_selected_delE_ignorers = lapply(selected_delE_ignorers, make_cpm_plot, groups = c("EV","HCD","HCN73","dTAD","delE"))
plots_selected_delE_responders = lapply(selected_delE_responders, make_cpm_plot, groups = c("EV","HCD","HCN73","dTAD","delE"))

#plot using patchwork
patchwork_plots_selected_delE_ignorers=wrap_plots(plots_selected_delE_ignorers[1:length(selected_delE_ignorers)]) + plot_layout(ncol=5)

patchwork_plots_selected_delE_responders=wrap_plots(plots_selected_delE_responders[1:length(selected_delE_responders)]) + plot_layout(ncol=5)

pdf(file = "plots/patchwork_logcpm_selected_delE_ignorers.pdf",width = 13, height = 8)
patchwork_plots_selected_delE_ignorers
dev.off()

pdf(file = "plots/patchwork_logcpm_selected_delE_responders.pdf",width = 13, height = 8)
patchwork_plots_selected_delE_responders
dev.off()



#Let's do a heatmap of EV, CICDUX4, and CICLEUTX
#with any genes up in CICDUX4 vs EV or CICLEUTX vs EV
up_CICDUX4_or_CICLEUTX = union(up_CICDUX4_vs_EV, up_CICLEUTX_vs_EV)

logcpm_EV_CICDUX4_CICLEUTX = logcpm[up_CICDUX4_or_CICLEUTX, c(1:6, 19:21)]


pdf(file = "plots/EV_HCD_HCL.pdf", width = 8, height = 6)
EV_HCD_HCL_heatmap = pheatmap(logcpm_EV_CICDUX4_CICLEUTX,
                        cluster_cols = F,
                        cluster_rows = T,
                        scale = "row",
                        show_rownames = F,
                        angle_col = 45,
                        annotation_colors = heatmap_colors,
                        annotation_col = sample_annotations)
dev.off()


#How about a volcano plot of CICLEUTX vs EV where genes ON in CICDUX4 vs EV are maroon?
#FDR cutoff plotted is at -log10(q) = 5, or q < 0.00001.
#log2FC cutoffs are at +/- 2

CICLEUTXvsEV_CICDUX4_on_genes_plot = ggplot(data = CICLEUTX.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICLEUTX.EV_results[CICLEUTX.EV_results$hgnc_symbol %in% up_CICDUX4_vs_EV,], 
             color = "maroon", size = 4) +
  theme_classic() +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICLEUTX.EV_results[CICLEUTX.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5", "VGF", "DUSP6","FOXB1","DUSP4","DUSP6","SPRED1","SPRED2","SPRED3","SPRY4"),],
                  aes(label = hgnc_symbol), nudge_y = 0, nudge_x = 0, box.padding = 0.5) +
  ggtitle("CICLEUTX vs. EV, maroon = genes on in CICDUX4 vs EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CL_vs_EV_maroon_CICDUX4_on.pdf",width = 8, height = 6)
CICLEUTXvsEV_CICDUX4_on_genes_plot
dev.off()

#and the opposite, CICDUX4 vs EV where genes ON in CICLEUTX vs EV are maroon
#FDR cutoff is plotted at -log10(q) = 5, or q < 0.00001
#log2FC cutoffs are +/-2

CICDUX4vsEV_CICLEUTX_on_genes_plot = ggplot(data = CICDUX4.EV_results, aes(x = logFC, y = -log10(fdr_adj))) +
  geom_point(color = all_palette[1], size = 4) +
  geom_point(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% up_CICLEUTX_vs_EV,], 
             color = "maroon", size = 4) +
  theme_classic() +
  geom_hline(yintercept = 5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -2, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = CICDUX4.EV_results[CICDUX4.EV_results$hgnc_symbol %in% c("ETV1","ETV4","ETV5", "VGF", "DUSP6"),],
                  aes(label = hgnc_symbol), nudge_y = 3, nudge_x = 0, box.padding = 0.5) +
  ggtitle("CICDUX4 vs. EV, maroon = genes on in CICLEUTX vs EV") +
  xlab("log2FC") +
  ylab("-log10(q)") +
  coord_cartesian(xlim = c(-10, 20), ylim = c(-2, 35))

pdf(file = "plots/CD4_vs_EV_maroon_CICLEUTX_on.pdf",width = 8, height = 6)
CICDUX4vsEV_CICLEUTX_on_genes_plot
dev.off()


#make some logcpm plots for CICDUX4 and CICLEUTX
CD_CL = c("ETV1","ETV4","ETV5","SHC3","VGF","DUSP6","MUC5AC")

plots_selected_HCL = lapply(CD_CL, make_cpm_plot, groups = c("EV","HCD","HCL"))


#plot using patchwork
patchwork_plots_selected_HCL=wrap_plots(plots_selected_HCL[1:length(CD_CL)]) + plot_layout(ncol=5)


pdf(file = "plots/patchwork_logcpm_HCD_HCL.pdf",width = 13, height = 6)
patchwork_plots_selected_HCL
dev.off()




#write out datasets to share

write.table(tibble::rownames_to_column(logcpm, "Gene"), file = "outputs/logcpm_May2024.txt", sep = "\t", quote = F, row.names = F)
write.table(CICDUX4.EV_results, file = "outputs/CICDUX4_vs_EV_May2024.txt", sep = "\t", quote = F, row.names = F)
write.table(CICex20_NUTM1ex6.EV_results, file = "outputs/CICex20-NUTM1ex6_vs_EV_May2024.txt", sep = "\t", quote = F, row.names = F)
write.table(CICex18_NUTM1ex3.EV_results, file = "outputs/CICex18-NUTM1ex3_vs_EV_May2024.txt", sep = "\t", quote = F, row.names = F)
write.table(CICLEUTX.EV_results, file = "outputs/CICLEUTX_vs_EV_May2024.txt", sep = "\t", quote = F, row.names = F)

write.table(delE_ignorers, file = "outputs/E_ignorers_May2024.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(delE_responders, file = "outputs/E_responders_May2024.txt", sep = "\t", quote = F, row.names = F, col.names = F)

write.table(tibble::rownames_to_column(master, "Gene"), file = "outputs/counts_May2024.txt", sep = "\t", quote = F, row.names = F)


#Reviewers also asked for GO analysis of genes upregulated by each fusion.
#I will generate these analyses using gProfiler on the lists of genes upregulated in each fusion vs EV comparison.

#CICDUX4 vs EV
DUX4_gprofiler = gost(up_CICDUX4_vs_EV, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
DUX4_GO = DUX4_gprofiler$result
trimmed_DUX4_GO = dplyr::select(DUX4_GO, p_value, term_id, term_name)
write.table(trimmed_DUX4_GO, file = "outputs/CICDUX4_vs_EV_gProfiler.txt", sep = "\t", quote = F, row.names = F)

#CN73 [CICex20_NUTM1ex6] vs EV
CN73_gprofiler = gost(up_CICex20_NUTM1ex6_vs_EV, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
CN73_GO = CN73_gprofiler$result
trimmed_CN73_GO = dplyr::select(CN73_GO, p_value, term_id, term_name)
write.table(trimmed_CN73_GO, file = "outputs/CICex20_NUTM1ex6_vs_EV_gProfiler.txt", sep = "\t", quote = F, row.names = F)

#CN84 [CICex18_NUTM1ex3] vs EV
CN84_gprofiler = gost(up_CICex18_NUTM1ex3_vs_EV, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
CN84_GO = CN84_gprofiler$result
trimmed_CN84_GO = dplyr::select(CN84_GO, p_value, term_id, term_name)
write.table(trimmed_CN84_GO, file = "outputs/CICex18_NUTM1ex3_vs_EV_gProfiler.txt", sep = "\t", quote = F, row.names = F)

#CICLEUTX vs EV
LEUTX_gprofiler = gost(up_CICLEUTX_vs_EV, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
LEUTX_GO = LEUTX_gprofiler$result
trimmed_LEUTX_GO = dplyr::select(LEUTX_GO, p_value, term_id, term_name)
write.table(trimmed_LEUTX_GO, file = "outputs/CICLEUTX_vs_EV_gProfiler.txt", sep = "\t", quote = F, row.names = F)


#Read in ATXN1::DUX4 vs EV gProfiler results, to do an upset plot with the others
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_ATXN1_DUX4/hg38_data")
trimmed_ATXN1_GO = fread("outputs/ATXN1DUX4_vs_EV_gProfiler.txt")
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_293T_delE_CICNUTM1")

upSetList = list(CIC_DUX4 = trimmed_DUX4_GO$term_name, 
                 C20N6 = trimmed_CN73_GO$term_name,
                 C18N3 = trimmed_CN84_GO$term_name,
                 CIC_LEUTX = trimmed_LEUTX_GO$term_name,
                 ATXN1_DUX4 = trimmed_ATXN1_GO$term_name)

pdf(file = "plots/GO_UpSeT.pdf", width = 6, height = 5, onefile = F)
upset(fromList(upSetList), order.by = c("freq"), text.scale = c(2,2,2,2,2,1.5))
dev.off()

#I will also write out the list of significantly up genes for each condition, for the sake of completeness.
write.table(up_CICDUX4_vs_EV, file = "outputs/up_CICDUX4_vs_EV.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(up_CICex20_NUTM1ex6_vs_EV, file = "outputs/up_CICex20_NUTM1ex6_vs_EV.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(up_CICex18_NUTM1ex3_vs_EV, file = "outputs/up_CICex18_NUTM1ex3_vs_EV.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(up_CICLEUTX_vs_EV, file = "outputs/up_CICLEUTX_vs_EV.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#A reviewer also asked for a list of genes specifically regulated by ATXN1::DUX4, so let's load that in
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_ATXN1_DUX4/hg38_data")
up_ATXN1DUX4_vs_EV = fread("outputs/up_ATXN1DUX4_vs_EV.txt", header = F)
up_ATXN1DUX4_vs_EV = up_ATXN1DUX4_vs_EV$V1
setwd("/Volumes/cuyler/ucsf_okimoto_lab/cuyler_293T_delE_CICNUTM1")

#union of all genes sig. up by one of the CIC-fusions (this is the same as up_any_full_fusion, but recalculating for sanity)
up_all_CIC = union(c(c(up_CICDUX4_vs_EV, up_CICex20_NUTM1ex6_vs_EV), up_CICex18_NUTM1ex3_vs_EV), up_CICLEUTX_vs_EV)

up_ATXN1DUX4_only = up_ATXN1DUX4_vs_EV[!up_ATXN1DUX4_vs_EV %in% up_all_CIC]
write.table(up_ATXN1DUX4_only, file = "outputs/up_ATXN1DUX4_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#any GO for these genes?
AD_only_gprofiler = gost(up_ATXN1DUX4_only, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP")
AD_only_GO = AD_only_gprofiler$result
trimmed_AD_only_GO = dplyr::select(AD_only_GO, p_value, term_id, term_name)
write.table(trimmed_AD_only_GO, file = "outputs/ATXN1DUX4_only_gProfiler.txt", sep = "\t", quote = F, row.names = F)


#A reviewer also asked for genes upregulated by CIC::DUX4 but not CIC::NUTM1 
#(I will combine ex20-ex6 and ex18-ex3 genes as "CIC::NUTM1 responsive")

up_either_CICNUTM1 = union(up_CICex20_NUTM1ex6_vs_EV, up_CICex18_NUTM1ex3_vs_EV)

up_CICDUX4_not_CICNUTM1 = up_CICDUX4_vs_EV[!up_CICDUX4_vs_EV %in% up_either_CICNUTM1] 
write.table(up_CICDUX4_not_CICNUTM1, file = "outputs/up_CICDUX4_not_CICNUTM1.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#any GO for these genes?
CD4_not_CN_gprofiler = gost(up_CICDUX4_not_CICNUTM1, organism = "hsapiens", correction_method = "gSCS", source = "GO:BP", significant = F)
CD4_not_CN_GO = CD4_not_CN_gprofiler$result

#no - none of them have a corrected p value less than 0.1. The closest is animal organ morphogenesis at 0.15
  
