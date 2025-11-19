#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
options(stringsAsFactors=FALSE)
theme_set(ggpubr::theme_pubr())

# load correlation data matrix
df <- read.csv("data/correlations.csv") %>%
      # log-transform with pseudo-count
      mutate(GEq=log10(GEq+0.1),
             PFU=log10(PFU+0.1),
             IgG=log10(IgG+0.1),
             Dose=log10(Dose+0.1)) %>%
      # set NHPs as rownames
      column_to_rownames("NHP")

# calculate correlation using spearman, which does not assume normality
# unfortunately we cannot include outcome since it's categorical
corr.mat <- df %>%
            select(-Outcome) %>%
            as.matrix() %>%
            Hmisc::rcorr(type="spearman")
pval <- corr.mat$P %>%
        as.data.frame() %>%
        rownames_to_column("A") %>%
        reshape2::melt(id.vars="A",
                       variable.name="B", 
                       value.name="p")
corr.mat <- corr.mat$r %>%
            as.data.frame() %>%
            rownames_to_column("A") %>%
            reshape2::melt(id.vars="A",
                           variable.name="B", 
                           value.name="r") %>%
            left_join(pval, by=c("A", "B")) %>%
            # format factors and p-value. Transform r to r2
            mutate(A=factor(A, levels=colnames(corr.mat$r)),
                   B=factor(B, levels=colnames(corr.mat$r)),
                   p.format=gtools::stars.pval(p),
                   r=r*r)
corr.mat %>%
  ggplot(aes(A, B)) +
  geom_tile(aes(fill=r), col="black") +
  scale_fill_gradientn(breaks=c(0, 1), 
                       colors=c("white", "black"),
                       limits=c(0, 1)) +
  geom_text(aes(label=p.format), col="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x=element_blank(),
       y=element_blank(),
       fill=bquote(r^2)) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1),
        axis.ticks=element_blank(), panel.border=element_rect(fill=NA, color="black", linewidth=1.5))
ggsave("analysis/correlations.png",
       units="in", width=3, height=3)

# clean up 
rm(pval, corr.mat)

## outcome predictions ---------------------------------------------------------
df$Outcome.tf <- as.numeric(df$Outcome=="Fatal")

# loop over predictors and make GLMs
predictors <- c("IgG", "GEq", "PFU", "Dose") 
glms <- predictors %>%
        paste0("Outcome.tf ~ ", .) %>%
        lapply(as.formula) %>%
        lapply(glm, data=df, family=binomial)
names(glms) <- predictors

# get ROCs
rocs <- glms %>%
        lapply(function(i) {
          pROC::roc(response=i$y, predictor=i$fitted.values, quiet=TRUE)
        })

# get threshold stats, AUC, and significance
stat.mat <- rocs %>%
            lapply(pROC::coords, 
                   x="best", 
                   ret=c("threshold", "specificity", "sensitivity"), 
                   best.method="closest.topleft") %>%
            do.call(rbind, .) %>%
            rownames_to_column("Predictor") %>%
            mutate(AUC=unlist(lapply(rocs, pROC::auc)),
                   Significance=unlist(lapply(glms, function(i) { 
                     summary(i)$coefficients[2, 4] 
            })))

# calculate FP, FN, and accuracy
stat.mat <- df %>%
            select(-Outcome.tf) %>%
            reshape2::melt(id.vars="Outcome",
                           variable.name="Predictor") %>%
            left_join(select(stat.mat, Predictor, threshold),
                      by="Predictor") %>%
            mutate(Predicted=(value >= threshold),
                   Predicted=factor(Predicted, levels=c(FALSE, TRUE),
                                    labels=c("Survived", "Fatal")),
                   # where do predicted & ground truth not match?
                   Agreement=factor(Predicted==Outcome,
                                    levels=c(FALSE, TRUE),
                                    labels=c("Disagree", "Agree"))) %>%
            group_by(Predictor, Agreement) %>%
            summarise(Cases=n(),
                      .groups="drop") %>%
            reshape2::dcast(Predictor ~ Agreement, 
                            value.var="Cases", fill=0) %>%
            mutate(Accuracy=Agree/(Agree + Disagree)) %>%
            select(-Disagree, -Agree) %>%
            left_join(stat.mat, by="Predictor") %>%
            mutate(threshold=(10^threshold)-0.1,
                   Significance=format.pval(as.numeric(Significance), digits=2)) %>%
            mutate_if(is.numeric, round, digits=2)
colnames(stat.mat) <- c("Predictor", "Accuracy", "Threshold", "Specificity", 
                        "Sensitivity", "AUC", "Significance")
write.csv(stat.mat, "analysis/outcome-predictors.csv", row.names=FALSE)
