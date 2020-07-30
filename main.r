#
packages.to.load = c("readxl", "tidyr", "VennDiagram", "dplyr", "edgeR", "preprocessCore", "ggplot2", "gridExtra", "lme4", "lmerTest", "ggrepel","tidyverse", "DT")
for(p in packages.to.load){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
rm(list = c("p", "packages.to.load"))


###########
############
############

dat <- read.table( "proteinGroups.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t", fill =  TRUE)
colnames(dat)

dat <- dat[ dat$Reverse != "+" , ]
dat <- dat[ dat$Potential.contaminant != "+" , ]


# optional
# dat <- dat[ dat$Razor...unique.peptides >= 2, ]
#

# data = data.frame(fc.raw.file = rep(c("file A", "file B", "file C"), each=81),
#                   RT = c(20:100),
#                   intensity = c(rnorm(81, mean=20), rnorm(81, mean=10), rnorm(81, mean=30)))


# optional QC
colnames(dat)
evidence <- read.table( "evidence.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t", fill =  TRUE)
evidence <- evidence[ evidence$Reverse != "+", ]
evidence <- evidence[ evidence$Potential.contaminant != "+", ]
evidence <- evidence[ c( grep("Cont_bio2_A", evidence$Raw.file ), grep("TNF_Bio1_A", evidence$Raw.file)) , ]

colnames(evidence)
data <- data.frame(
  fc.raw.file = evidence$Raw.file,
  RT = evidence$Retention.time, 
  intensity = log2(evidence$Intensity)
  
)
summary(data)
PTXQC::plot_TIC(data, c(10, 100), c(0, 35) )

#rm(evidence); rm(data)
gc()
# 
##
####

# 
##
####


colnames(dat)
dat$ID <- paste(dat$Majority.protein.IDs, dat$Gene.names, sep = ".")
dat <- dat[ , c(131, 82:93)]
dat[ dat == 0] <- NA
ids <- dat$ID
row.names(dat) <- dat$ID
dat$ID <- NULL
# VIM::matrixplot(dat, sortby = 2)
Amelia::missmap(dat[ , ])
###

## Remove rows with more than 50% NA
dat <- dat[which(rowMeans(!is.na(dat[ , ])) > 0.5), ]
Amelia::missmap(dat[ , ])
##

# 
# nPCs <- missMDA::estim_ncpPCA(dat[ , -1])
# nPCs$ncp
# dat.complete <- imputePCA(dat, ncp = 2, scale = TRUE)

dat <- log2(dat)

dat_pca_methods <- pcaMethods::pca(dat, nPcs=2, method="ppca", center = TRUE)
dat_imp_pcamethods <- pcaMethods::completeObs(dat_pca_methods)
##


##
# 
# cor.dat <- cor(dat, use="complete.obs", method="spearman")
# 
# cor.dat <- cor(dat, use="complete.obs", method="spearman") %>%
#   as.data.frame() %>%
#   mutate(var1 = rownames(.)) %>%
#   gather(var2, value, -var1) %>%
#   arrange(desc(value))
# 
# 
# for(i in 1:ncol(dat)){
#   name <- colnames(dat)[i]
#   
#   
#   cor.to.use <- cor.dat[ cor.dat$var2 == name, ][with(cor.dat[ cor.dat$var2 == name, ] , order(value,decreasing = T)),]
#   cor.to.use[ 2,3 ]
#   
#   
#   dat[ , i ]
#   
# }
##


# qc plots

pca_com <- prcomp((t(dat_imp_pcamethods)), scale=TRUE, center = TRUE)
pca_group  <- data.frame(Condition = sub("_.*", "", colnames(dat_imp_pcamethods)) ,
                         Name = colnames(dat_imp_pcamethods))

ggbiplot::ggbiplot(pca_com,
                   Conditions = pca_group$Condition, var.axes = FALSE)+
  geom_point(aes(color = pca_group$Condition), size = 3)+
  geom_text(aes(label=pca_group$Name), size=3,
             position=position_jitter(width=0.2, height=0.2)) +
  labs(title = "PCA - log2(Intensity)", color = "Condition")+
  # scale_shape_manual(values = c(15,16,7:9,17,11)) +
  theme_bw()
##



##
dat.long <- reshape2::melt(dat_imp_pcamethods)

## count number of missing values

# plots for data quality check
ggplot(data = dat.long, aes(x = Var2, y = (value))) +
  geom_boxplot() +
  labs(title = "Boxplot for raw data") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 

ggplot(dat.long, aes(x=value, color=Var2))+ 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for raw data") + 
  theme_classic()



#########################
#########################
# Normalization method: https://pwilmart.github.io/IRS_normalization/understanding_IRS.html

#####
## Sample loading
#####
# This assumes the same amount of protein labelled in each sample and the total signal in each channel summing to the same value within the batch.
dat.sl <- dat.long %>%
  group_by(Var2) %>%
  dplyr::summarise(SumEachSample = sum(value, na.rm = TRUE)) %>%
  mutate(MeanofSumWithinExp = mean(SumEachSample, na.rm = TRUE),
         correct.factor.sl = MeanofSumWithinExp/SumEachSample) %>%
  right_join(dat.long, by = c("Var2")) %>%
  mutate(Abundances.sl = value * correct.factor.sl) 



# plots for data quality check
ggplot(data = dat.sl, aes(x = Var2, y = (Abundances.sl))) +
  geom_boxplot() +
  labs(title = "Boxplot loading normalized data") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 

ggplot(dat.sl, aes(x=Abundances.sl, color=Var2))+ 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for loadng normalized data") + 
  theme_classic()
# 

# convert to data frame in wide format
dat.norm <- dat.sl %>%
  dplyr::select(Var1, Var2, Abundances.sl) %>%
  pivot_wider(id_cols = c("Var1"), names_from = "Var2", values_from = "Abundances.sl")

dat.norm <- as.data.frame(dat.norm)
row.names(dat.norm) <- dat.norm$Var1
dat.norm$Var1 <- NULL


pca_com <- prcomp((t(dat.norm)), scale=TRUE, center = TRUE)
pca_group  <- data.frame(Condition = sub("_.*", "", colnames(dat.norm)) ,
                         Name = colnames(dat.norm))

ggbiplot::ggbiplot(pca_com,
                   Conditions = pca_group$Condition, var.axes = FALSE)+
  geom_point(aes(color = pca_group$Condition), size = 3)+
  geom_text(aes(label=pca_group$Name), size=3,
            position=position_jitter(width=0.2, height=0.2)) +
  labs(title = "PCA - log2(Intensity)", color = "Condition")+
  # scale_shape_manual(values = c(15,16,7:9,17,11)) +
  theme_bw()
##




##################
##################
# plot limma results
plot.limma <- function(fit2, coeff){
  tab <- topTable(fit2, coef = coeff, adjust="BH", number=Inf, sort.by = "logFC", resort.by = "P")
  tab <- data.frame(tab) %>%  rownames_to_column(var = "Accession")
  p1 <- ggplot(data=tab, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(alpha=0.5)+
    geom_point(data= filter(tab, adj.P.Val < 0.05 ), #& Accession %in% peptides.filter
               aes(x=logFC, y=-log10(adj.P.Val)), col=AZcolor$Mulberry)+
    geom_hline(aes(yintercept = -log10(0.05)), linetype = 6, col = "black") + 
    labs(x= expression( ~ log[2] ~ FC), ylab = expression( ~ -log[10] ~ (Adjusted.p.value)),        title = paste(colnames(fit2$coefficients)[coeff])) + geom_text(aes(label=ifelse(adj.P.Val<0.05,as.character(Accession),'')),hjust=0,vjust=0, cex = 1.9)  + 
    theme_bw()
  
  p1
  ggsave(plot = p1, filename = paste(colnames(fit2$coefficients)[coeff], ".png",sep = "") )
  
}
#############
#########################
#########################
##
## limma
##
#########################
#########################

# construct model matrix
metadata <- data.frame(Condition = sub("_.*", "", colnames(dat.norm)) ,
           Name = colnames(dat.norm))
metadata$Condition <- sub("Intensity.", "", metadata$Condition)

mod.matrix1 <- model.matrix( ~ 0 + Condition,  data = metadata) 
##
colnames(mod.matrix1) <- sub("Condition", "", colnames(mod.matrix1))
colnames(mod.matrix1)

#########################
#########################
##
## run model
##
#########################
#########################
# Fit the model for the treatment effects to obtain the pooled variance
fit <- lmFit((dat.norm), mod.matrix1) # Coefficients not estimable: Time8h
# fit1 <- eBayes(fit)
names(fit)
head(fit$coef)
# 
colnames(mod.matrix1)

contr.mat <- makeContrasts(
  Ins - Cont,
  
  TNF - Cont, 
  
  TNF.Ins - Cont, 
  
  levels = mod.matrix1 )

#
fit2 <-contrasts.fit(fit, contr.mat)		
fit2 <- eBayes(fit2)
############
result <- decideTests(fit2, method="separate",adjust.method="BH",p.value=0.05,lfc=0)
summary(decideTests(fit2, method="separate",adjust.method="BH",p.value=0.05,lfc=0))
vennDiagram(result[ , 1:3], circle.col=c(1,2,3 ),counts.col =  "black", cex = 1.2 ) # include=c("up"



diff.ids <- which(result[,1]!=0 | result[,2]!=0 | result[,3]!=0)
# de.common <- which(result[,1]!=0 | result[,2]!=0)
#
pheatmap::pheatmap(fit2$coefficients[ (diff.ids), 1:3], scale = "none", cluster_rows = T, cluster_cols = F )

# Biomarkers
row.names(fit2$coefficients[ (diff.ids), ])

temp = melt(dat.norm[ "P12265.Gusb",] )
temp$variable  <- sub("_.*", "", (temp$variable))
lattice::bwplot( value ~ variable , data=temp)

temp = melt(dat.norm[ "Q9EQU5.Set",] )
temp$variable  <- sub("_.*", "", (temp$variable))
lattice::bwplot( value ~ variable , data=temp)

# Beta-glucuronidase - Plays an important role in the degradation of dermatan and keratan sulfates.
# Protein SET - Multitasking protein, involved in apoptosis, transcription, nucleosome assembly and histone chaperoning.
####
#####


coeff = 1 
tab <- topTable(fit2, coef = coeff, adjust="BH", number=Inf, sort.by = "logFC", resort.by = "P")
  tab <- data.frame(tab) %>%  rownames_to_column(var = "Accession")
  p1 <- ggplot(data=tab, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(alpha=0.5)+
    geom_point(data= filter(tab, adj.P.Val < 0.05 ), #& Accession %in% peptides.filter
               aes(x=logFC, y=-log10(adj.P.Val)), col=AZcolor$Mulberry)+
    geom_hline(aes(yintercept = -log10(0.05)), linetype = 6, col = "black") + 
    labs(x= expression( ~ log[2] ~ FC), ylab = expression( ~ -log[10] ~ (Adjusted.p.value)),        title = paste(colnames(fit2$coefficients)[coeff])) + geom_text(aes(label=ifelse(adj.P.Val<0.05,as.character(Accession),'')),hjust=0,vjust=0, cex = 1.9)  + 
    theme_bw()
  
  p1

##

result <- decideTests(fit2, method="separate",adjust.method="none",p.value=0.05,lfc=0)
summary(decideTests(fit2, method="separate",adjust.method="none",p.value=0.05,lfc=0))
vennDiagram(result[ , 1:3], circle.col=c(1,2,3 ),counts.col =  "black", cex = 1.2 ) # include=c("up"


#### 
# gsea
###

library("biomaRt")
# Select ensembl database and mouse dataset:
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl" ,host="www.ensembl.org")
# Map IDs to GO

mapGO <- getBM(attributes=c('uniprotswissprot', "name_1006"), mart = ensembl)
# Remove blanks ("")
mapGO <- mapGO[mapGO[,2]!="",]
# Check the 10 first rows to see what we got:
mapGO[1:10,]

##
myGsc <- loadGSC(mapGO)

P.Value <- tab$P.Value
names(P.Value) <- sub( "\\..*" , "", tab$Accession)
logFC <- tab$logFC
names(logFC) <- sub("\\..*", "", tab$Accession)

gsaRes <- runGSA(P.Value,logFC,gsc=myGsc,gsSizeLim=c(5,300))
exploreGSAres(gsaRes)

