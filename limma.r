#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
totres <- data.frame()
col.stat <- c(Baseline="#fef0d9",
              Early="#fdbb84",
              Middle="#e34a33",
              Late="#b30000",
              Recovered="grey60")
col.regs <- c(Down="#377eb8", None="grey60", Up="#e41a1c")
scalefact <- 1.5

# helper functions
normalize.counts <- function(matx, minsamp=3, thld=10) {
  # remove genes (rows) with < minsamp samples with threshold counts
  matx[apply(matx, 1, function(i) {sum(i > thld)}) >= minsamp, ]%>%
    # format as DGElist and calculate normalization factors
    DGEList() %>%
    calcNormFactors()
}
diffexpr <- function(exprmat, modmat, contr, plim=0.05, llim=1) {
  # run differential expression and format output
  x <- exprmat %>%
       voom(modmat) %>%
       lmFit(design=modmat) %>%
       contrasts.fit(contr) %>%
       eBayes() %>%
       topTable(number=Inf) %>%
       rownames_to_column("gene") %>%
       dplyr::rename(padj=adj.P.Val,
                     pvalue=P.Value,
                     lfc=logFC,
                     avexpr=AveExpr) %>%
       select(gene, avexpr, lfc, padj, pvalue)
  # extract name of comparison
  x$Comparison <- names(which(contr[, 1]==1))
  # add significance and annotate regulation
  x$Significant <- (x$padj < plim & abs(x$lfc) > llim)
  x$Regulation <- "None"
  x[x$Significant & x$lfc > llim, "Regulation"] <- "Up"
  x[x$Significant & x$lfc < -llim, "Regulation"] <- "Down"
  # return full matrix
  return(x)
}
format.ipa <- function(res) {
  res %>%
    select(gene, Comparison, lfc, padj) %>%
    reshape2::melt(id.vars=c("gene", "Comparison"),
                   measure.vars=c("lfc", "padj")) %>%
    mutate(label=paste0(Comparison, ".", variable)) %>%
    reshape2::dcast(gene ~ label, value.var="value")
}

## inputs ----------------------------------------------------------------------
# metadata
meta <- read.csv("data/metadata.csv", 
                 na.strings="") %>%
        mutate(Sample.ID=str_replace(Sample.ID, "-", ".")) %>%
        # QC filter: AGM 5000 & two cyno samples with low [RNA]
        filter(Filter.QC)
# update status & DPI to factor
meta <- meta %>%
        mutate(Status=factor(Status, levels=c("Baseline", "Early", "Middle",
                                              "Late", "Recovered")),
               Species=factor(Species, levels=c("MCy", "AGM")))
# set row names
rownames(meta) <- meta$Sample.ID

# counts
cmat <- read.csv("data/counts-thresholded.csv",
                 row.names=1)

# align rows and columns
x <- intersect(colnames(cmat), rownames(meta))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## normalization ---------------------------------------------------------------
# remove survivors 
meta <- meta %>%
        filter(Outcome != "Survived") %>%
        select(-Outcome)
cmat <- cmat[ , rownames(meta)]

# remove positive and negative controls
cmat <- cmat[!str_detect(rownames(cmat), "POS_|NEG_"), ]

# format full matrix to get log2 CPM
cpm2 <- cmat %>%
        normalize.counts() %>%
        cpm()
# format CPM as log2 for remainder of analysis
cpm2 <- log2(cpm2)

## timeline --------------------------------------------------------------------
pA <- meta %>%
      select(Study.day, Status) %>%
      distinct() %>%
      ggplot(aes(Study.day, 1, group=1)) +
      geom_line(linewidth=2) +
      geom_point(aes(fill=Status), pch=21, size=5) +
      geom_text(aes(label=Study.day), nudge_y=0.15) +
      ylim(0.9, 1.2) +
      scale_fill_manual(values=col.stat) +
      labs(fill=element_blank()) +
      theme_void() +
      theme(legend.position="top")
pA
ggsave("analysis/timeline.png", scale=scalefact,
       units="in", width=3, height=0.5)

## individual PCA --------------------------------------------------------------
# cyno only
c <- filter(meta, Species=="MCy")
c <- cpm2[ , rownames(c)]
# run PCA
pca <- c %>%
       t() %>%
       prcomp()
# calculate PCs
pcs <- pca$sdev^2
pcs <- pcs/sum(pcs)
pcs <- round(100*pcs)
pcs <- c(paste0("PC1 (", pcs[1], "%)"),
         paste0("PC2 (", pcs[2], "%)"))
# extract PCA matrix
pca <- pca$x %>%
       as.data.frame() %>%
       rownames_to_column("Sample.ID") %>%
       select(Sample.ID, PC1, PC2) %>%
       left_join(meta, by="Sample.ID")
# plot!
pB <- pca %>%
      ggplot(aes(PC1, PC2, fill=Status)) +
      geom_point(pch=21, size=2, alpha=0.8) +
      geom_vline(xintercept=0, linetype=2) +
      geom_hline(yintercept=0, linetype=2) +
      scale_fill_manual(values=col.stat) +
      labs(x=pcs[1],
           y=pcs[2],
           fill=element_blank(),
           title="MCy") +
      theme(legend.position=c(0.15, 0.8))
pB
ggsave("analysis/pca-cyno.png", scale=scalefact, 
       units="in", width=3, height=2)
rm(c, pca, pcs)

# AGM only
a <- filter(meta, Species=="AGM")
a <- cpm2[ , rownames(a)]

# run PCA
pca <- a %>%
      t() %>%
      prcomp()
# calculate PCs
pcs <- pca$sdev^2
pcs <- pcs/sum(pcs)
pcs <- round(100*pcs)
pcs <- c(paste0("PC1 (", pcs[1], "%)"),
         paste0("PC2 (", pcs[2], "%)"))
# extract PCA matrix
pca <- pca$x %>%
        as.data.frame() %>%
        rownames_to_column("Sample.ID") %>%
        select(Sample.ID, PC1, PC2) %>%
        left_join(meta, by="Sample.ID")
# plot!
pC <- pca %>%
      ggplot(aes(PC1, PC2, fill=Status)) +
      geom_point(pch=21, size=2, alpha=0.8) +
      geom_vline(xintercept=0, linetype=2) +
      geom_hline(yintercept=0, linetype=2) +
      scale_fill_manual(values=col.stat) +
      labs(x=pcs[1],
           y=pcs[2],
           fill=element_blank(),
           title="AGM") +
      theme(legend.position=c(0.15, 0.8))
pC
ggsave("analysis/pca-agms.png", scale=scalefact, 
       units="in", width=3, height=2)
rm(a, pca, pcs)

## cyno DE ---------------------------------------------------------------------
# subset cyno
m <- filter(meta, Species=="MCy")
c <- cmat[ , rownames(m)] %>%
     normalize.counts()
x <- droplevels(m$Status)
# format model matrix
m <- model.matrix(~0+x)
colnames(m) <- levels(x)
rm(x)
# run contrasts
x <- diffexpr(c, m, makeContrasts(Early - Baseline, 
                                  levels=colnames(m)))
y <- diffexpr(c, m, 
              makeContrasts(Middle - Baseline, 
                            levels=colnames(m))) 
z <- diffexpr(c, m, 
              makeContrasts(Late - Baseline, 
                            levels=colnames(m))) 
dexp <- rbind(x, y) %>%
  rbind(z) %>%
  mutate(Comparison=factor(Comparison, levels=c("Early", "Middle", "Late")),
         Regulation=factor(Regulation, levels=c("Down", "None", "Up")))
rm(x, y, z)

# save for IPA
dexp %>%
  format.ipa() %>% 
  write.csv("data/ipa-input-cyno.csv",
            row.names=FALSE)

# volcano plot
# get top 10 genes by fold change
glist <- dexp %>%
         filter(Significant) %>%
         group_by(Comparison) %>%
         top_n(n=10, wt=abs(lfc)) %>%
         ungroup()
pD <- dexp %>%
      ggplot(aes(lfc, -log10(padj), size=Regulation, col=Regulation)) +
      geom_point(alpha=0.8) +
      ggrepel::geom_text_repel(data=glist, aes(label=gene), 
                               col="black", size=2, max.overlaps=100) +
      scale_size_manual(values=c(Down=1, None=0.5, Up=1)) +
      scale_color_manual(values=col.regs) +
      facet_wrap(~Comparison, nrow=1) +
      scale_x_continuous(limits=c(-3, 6), breaks=c(-2, 0, 2, 4, 6)) +
      scale_y_continuous(limits=c(0, 25), breaks=c(0, 5, 10, 15, 20, 25),
                         labels=c("1", "1e-5", "1e-10", "1e-15", 
                                  "1e-20", "1e-25")) +
      labs(x="Fold change (log2)",
           y="FDR-adjusted p-value",
           title="MCy")
pD
ggsave("analysis/volcano-cyno.png", scale=scalefact,
       units="in", width=4, height=2)

# add results to full list and save
totres <- dexp %>%
          mutate(Species="MCy") %>% 
          rbind(totres)

# clean up
rm(dexp, c, m, glist)

## AGM DE ----------------------------------------------------------------------
# subset AGM
m <- filter(meta, Species=="AGM")
c <- cmat[ , rownames(m)] %>%
     normalize.counts()
x <- droplevels(m$Status)
# format model matrix
m <- model.matrix(~0+x)
colnames(m) <- levels(x)
rm(x)
# run contrasts
x <- diffexpr(c, m, makeContrasts(Early - Baseline, 
                                  levels=colnames(m)))
y <- diffexpr(c, m, 
              makeContrasts(Middle - Baseline, 
                            levels=colnames(m))) 
z <- diffexpr(c, m, 
              makeContrasts(Late - Baseline, 
                            levels=colnames(m))) 
dexp <- rbind(x, y) %>%
  rbind(z) %>%
  mutate(Comparison=factor(Comparison, levels=c("Early", "Middle", "Late")),
         Regulation=factor(Regulation, levels=c("Down", "None", "Up")))
rm(x, y, z)

# save results for IPA
dexp %>%
  format.ipa() %>% 
  write.csv("data/ipa-input-agm.csv",
            row.names=FALSE)

# volcano plot
# get top 10 genes by fold change
glist <- dexp %>%
        filter(Significant) %>%
        group_by(Comparison) %>%
        top_n(n=10, wt=abs(lfc)) %>%
        ungroup()
pE <- dexp %>%
      ggplot(aes(lfc, -log10(padj), size=Regulation, col=Regulation)) +
      geom_point(alpha=0.8) +
      ggrepel::geom_text_repel(data=glist, aes(label=gene), 
                               col="black", size=2, max.overlaps=100) +
      scale_size_manual(values=c(Down=1, None=0.5, Up=1)) +
      scale_color_manual(values=col.regs) +
      facet_wrap(~Comparison, nrow=1) +
      scale_x_continuous(limits=c(-3, 6), breaks=c(-2, 0, 2, 4, 6)) +
      scale_y_continuous(limits=c(0, 25), breaks=c(0, 5, 10, 15, 20, 25),
                         labels=c("1", "1e-5", "1e-10", "1e-15", 
                                  "1e-20", "1e-25")) +
      labs(x="Fold change (log2)",
           y="FDR-adjusted p-value",
           title="AGM") 
pE
ggsave("analysis/volcano-agm.png", scale=scalefact,
       units="in", width=4, height=2)

# add results to full list
totres <- dexp %>%
          mutate(Species="AGM") %>%
          rbind(totres)
write.csv(totres, "data/de-results.csv", row.names=FALSE)

# plot volcano together
x <- cowplot::plot_grid(pA,
                        pB + theme(legend.position="none"), 
                        pC + theme(legend.position="none"), 
                        labels="AUTO", ncol=1, rel_heights=c(1, 3, 3))
y <- cowplot::plot_grid(pD, pE, labels=c("D", "E"), ncol=1)
cowplot::plot_grid(x, y, nrow=1)
ggsave("analysis/figure.png", scale=scalefact,
       units="in", width=7.5, height=4)

# clean up
rm(dexp, c, m, glist, pA, pB, pC, pD, pE, x, y)

## correlation -----------------------------------------------------------------
# get the intersecting DE genes for each comparison
x <- totres %>%
    reshape2::dcast(gene + Comparison ~ Species, 
                    value.var="lfc") %>%
    filter(abs(MCy) > 1 & abs(AGM) > 1)

# define the most divergent genes for each state
glist <- x %>%
         mutate(Dif=abs(AGM-MCy)) %>%
         filter(Dif > 1) %>%
         group_by(Comparison) %>%
         top_n(n=10, wt=Dif) %>%
         ungroup()

# plot it with correlation lines
x <- x %>%
     ggplot(aes(MCy, AGM)) +
     ggpmisc::stat_poly_line(col="darkgrey") +
     ggpmisc::stat_poly_eq(label.y=0.8) +
     geom_point(fill="grey40", col="grey30", pch=21, alpha=0.8) + 
     ggrepel::geom_text_repel(data=glist, aes(label=gene), size=2, force=50) +
     facet_wrap(~Comparison) +
     xlim(-3, 6) +
     ylim(-3, 6) +
     labs(x="MCy fold change (log2)",
          y="AGM fold change (log2)")
cowplot::plot_grid(x, labels="I")
ggsave("analysis/correlation.png", scale=scalefact,
       units="in", width=4.45, height=2)

## DE gene overlaps ------------------------------------------------------------
# early: both up and down-regulated
# define list
list(AGM=filter(totres, 
                Significant, 
                Species=="AGM",
                Comparison=="Early") %>%
          select(gene) %>%
          unlist(),
     Cyno=filter(totres, 
                 Significant, 
                 Species=="MCy",
                 Comparison=="Early") %>%
         select(gene) %>%
         unlist()) %>%
  VennDiagram::venn.diagram(filename="analysis/venn-early.png", 
                            imagetype="png",
                            disable.logging=TRUE, 
                            units="in", height=2, width=2, 
                            fontfamily="sans", 
                            cat.fontfamily="sans",
                            scaled=TRUE, 
                            cat.pos=c(-45, 45), cat.dist=0.1)

# middle
list(AGM=filter(totres, 
                Significant, 
                Species=="AGM",
                Comparison=="Middle") %>%
       select(gene) %>%
       unlist(),
     Cyno=filter(totres, 
                 Significant, 
                 Species=="MCy",
                 Comparison=="Middle") %>%
       select(gene) %>%
       unlist()) %>%
  VennDiagram::venn.diagram(filename="analysis/venn-middle.png", 
                            imagetype="png",
                            disable.logging=TRUE, 
                            units="in", height=2, width=2, 
                            fontfamily="sans", 
                            cat.fontfamily="sans",
                            scaled=TRUE, 
                            cat.pos=c(-45, 45), cat.dist=0.1)

# late
list(AGM=filter(totres, 
                Significant, 
                Species=="AGM",
                Comparison=="Late") %>%
       select(gene) %>%
       unlist(),
     Cyno=filter(totres, 
                 Significant, 
                 Species=="MCy",
                 Comparison=="Late") %>%
       select(gene) %>%
       unlist()) %>%
  VennDiagram::venn.diagram(filename="analysis/venn-late.png", 
                            imagetype="png",
                            disable.logging=TRUE, 
                            units="in", height=2, width=2, 
                            fontfamily="sans", cat.fontfamily="sans",
                            scaled=TRUE, 
                            cat.pos=c(-45, 45), cat.dist=0.1)

## IPA heatmap -----------------------------------------------------------------
col.hmap <- circlize::colorRamp2(c(0, 2.5, 5), c("white", "#fb6a4a", "#a50f15"))
png("analysis/ipa-heatmap.png", units="in", width=5, height=5.5, res=300)
read.csv("data/ipa-output.csv") %>%
  select(Pathway, 
         CynoEarly, AGMEarly, 
         CynoMiddle, AGMMiddle, 
         CynoLate, AGMLate) %>%
  column_to_rownames("Pathway") %>%
  as.matrix() %>% 
  ComplexHeatmap::Heatmap(col=col.hmap, 
                          name="z-score", 
                          border=TRUE, 
                          na_col="white", 
                          cluster_columns=FALSE,
                          column_labels=c("MCy early", "AGM early", 
                                          "MCy middle", "AGM middle",
                                          "MCy late", "AGM late"), 
                          column_names_rot=45)
dev.off()

# fin --------------------------------------------------------------------------
sessionInfo()
