###Load packages
library(ggpubr)
library(readxl)
##############------------------Dominance score for supplementary Figure S4A----------------------------------
snv <- read_xlsx("Supplementary Table 2.Somatic mutations in our SCLC cohort", sheet="Sheet1",skip=1)
snv <- as.data.frame(snv, stringsAsFactors = F)

sclc_driver <- read.delim("~/TCGA/Mutated_Genes_SCLC_239P_cBioportal.txt",header=T,fill=T,stringsAsFactors=F)
colnames(sclc_driver) <- c("gene","q","mut","mut_patient","Freq","cosmic")
sclc_driver$Freq <- as.numeric(gsub("%","",sclc_driver$Freq))
sclc_driver <- sclc_driver[which( sclc_driver$Freq >=5),]

cosmic <- read.delim("~/Census_allThu_Mar_19_09_31_24_2020.tsv",stringsAsFactors=F)

sclc_driver_snv <- snv[which((snv$Hugo_Symbol %in% sclc_driver$gene) | (snv$Hugo_Symbol %in% cosmic$Gene.Symbol)),]

sclc_driver_snv2 <- tidyr::separate(sclc_driver_snv, col="Tumor_Sample_Barcode", into=c("pid","sid"),sep="-")

pid_gene <- unique(sclc_driver_snv2[,c("pid","Hugo_Symbol")])

#Di is the dominance score for gene i, 
#Ni is the number of patients carrying non-silent mutations in gene i
#while dj is the number of LUAD driver mutations in patient j. Di=(Sum(1/dj))/Ni
#Ocurrence calculatoin	frequency of patients carrrying that driver mutation
##Gene dominant.score occurrence

##dj calculation
dj <- as.data.frame(table(pid_gene$pid))
dj$ratio <- 1/dj$Freq

##Ni calculation
gene.freq <- as.data.frame(table(pid_gene$Hugo_Symbol))
##
gene.freq$Occurrence <- gene.freq$Freq/40


##
Di <- NULL
for (i in gene.freq$Var1){
Ni <- gene.freq$Freq[gene.freq$Var1==i]
pid_tmp <- pid_gene$pid[pid_gene$Hugo_Symbol==i]

	dj_tmp <- dj$ratio[dj$Var1%in%pid_tmp]
	sum_dj <- sum(dj_tmp)
	sum_dj <- as.numeric(sum_dj)

	Di_tmp <- c(as.character(gene.freq$Var1[gene.freq$Var1==i]),sum_dj/Ni)
	Di <- rbind(Di,Di_tmp)
}

Di <- as.data.frame(Di,stringsAsFactors=F)
rownames(Di) <- Di[,1]
colnames(Di) <- c("Gene","Dominant.score")
dij <- merge(gene.freq, Di, by.x="Var1", by.y="Gene",all=T, no.dups=T)
dij$Dominant.score <- as.numeric(dij$Dominant.score)
colnames(dij)[1] <- "Gene"

dij.p <- dij[dij$Gene%in%cosmic$Gene.Symbol,]
dij.p$name <- ""
dij.p$name[(dij.p$Occurrence >=0.25) | (dij.p$Dominant.score>=0.03)] <- as.character(dij.p$Gene[(dij.p$Occurrence >=0.25) | (dij.p$Dominant.score>=0.03)])

##plot

ggscatter(dij.p, "Occurrence", "Dominant.score", color="Occurrence",label = "name", repel = TRUE, ylab="Driver dominant score",main="SCLC")

## To observe the drivers in LUAD
luad_driver <- read_xlsx("~/41586_2014_BFnature13385_MOESM21_ESM.xlsx",sheet="S_Table 4-MutSig2CV List",skip=2)
luad_sigs <- luad_driver$Gene[luad_driver$`Bonferroni p-value` < 0.025]

luad_snv <- read_xlsx("~/41586_2014_BFnature13385_MOESM21_ESM.xlsx",sheet="S_Table 5-Verification",skip=2)

luad_snv2 <- luad_snv[luad_snv$"Validation Judgement"==1,c("Hugo Symbol","Tumor Sample Barcode")]
luad_snv2 <- as.data.frame(luad_snv2)
pid_gene <- unique(luad_snv2) 
colnames(pid_gene) <- c("Hugo_Symbol","pid")
##dj
dj <- as.data.frame(table(pid_gene$pid))
dj$ratio <- 1/dj$Freq

##Ni calculation
gene.freq <- as.data.frame(table(pid_gene$Hugo_Symbol))
##
gene.freq$Occurrence <- gene.freq$Freq/length(unique(pid_gene$pid))


##
Di <- NULL
for (i in gene.freq$Var1){
Ni <- gene.freq$Freq[gene.freq$Var1==i]
pid_tmp <- pid_gene$pid[pid_gene$Hugo_Symbol==i]

	dj_tmp <- dj$ratio[dj$Var1%in%pid_tmp]
	sum_dj <- sum(dj_tmp)
	sum_dj <- as.numeric(sum_dj)

	Di_tmp <- c(as.character(gene.freq$Var1[gene.freq$Var1==i]),sum_dj/Ni)
	Di <- rbind(Di,Di_tmp)
}

Di <- as.data.frame(Di,stringsAsFactors=F)
rownames(Di) <- Di[,1]
colnames(Di) <- c("Gene","Dominant.score")
dij <- merge(gene.freq, Di, by.x="Var1", by.y="Gene",all=T, no.dups=T)
dij$Dominant.score <- as.numeric(dij$Dominant.score)
colnames(dij)[1] <- "Gene"

ggscatter(dij, "Occurrence", "Dominant.score", color="Occurrence", ylab="Driver dominant score",main="LUAD")

dij$name <- ""
dij$name[(dij$Occurrence >=0.3) | (dij$Dominant.score>=0.3)] <- as.character(dij$Gene[(dij$Occurrence >=0.3) | (dij$Dominant.score>=0.3)])
##plot
ggscatter(dij, "Occurrence", "Dominant.score", color="Occurrence",label = "name", repel = TRUE, ylab="Driver dominant score",main="LUAD")

### To observe the drivers in LUSC 
lusc_snv <- read_xls("~/nature11404-s2/data.file.S7.5.clinical.and.genomic.data.table.xls",sheet="LUSC_CpG_Filtered.patients.coun",range=cell_cols("N:BF"))
lusc_snv <- lusc_snv[-1,]
colnames(lusc_snv) <- lusc_snv[1,]
lusc_snv <- lusc_snv[-1,]
head(lusc_snv)
lusc_snv <- as.data.frame(lusc_snv)

ff <- function(x){sum(!is.na(x))}
gene.freq <- apply(lusc_snv,2,ff)
gene.freq <- as.data.frame(gene.freq)
gene.freq$Occurrence <- gene.freq$gene.freq/178
gene.freq <- gene.freq[-1,]
gene.freq$Var1 <- rownames(gene.freq)

rownames(lusc_snv) <- lusc_snv$`Tumor ID`
lusc_snv <- lusc_snv[,-1]
dj <- apply(lusc_snv,1,ff)
dj <- as.data.frame(dj)
dj$ratio <- 1/dj$dj
dj$pid <- rownames(dj)
###

Di <- NULL
for (i in gene.freq$Var1){
Ni <- gene.freq$gene.freq[gene.freq$Var1==i]

pid_tmp <- rownames(lusc_snv)[colnames(lusc_snv)==i]

	dj_tmp <- dj$ratio[dj$pid%in%pid_tmp]
	sum_dj <- sum(dj_tmp)
	sum_dj <- as.numeric(sum_dj)

	Di_tmp <- c(as.character(gene.freq$Var1[gene.freq$Var1==i]),sum_dj/Ni)
	Di <- rbind(Di,Di_tmp)
}

Di <- as.data.frame(Di,stringsAsFactors=F)
rownames(Di) <- Di[,1]
colnames(Di) <- c("Gene","Dominant.score")
dij <- merge(gene.freq, Di, by.x="Var1", by.y="Gene",all=T, no.dups=T)
dij$Dominant.score <- as.numeric(dij$Dominant.score)
colnames(dij)[1] <- "Gene"

ggscatter(dij, "Occurrence", "Dominant.score", color="Occurrence", ylab="Driver dominant score",main="LUSC")

dij$name <- ""
dij$name[(dij$Occurrence >=0.2) | (dij$Dominant.score>=1.0)] <- as.character(dij$Gene[(dij$Occurrence >=0.2) | (dij$Dominant.score>=1.0)])

ggscatter(dij, "Occurrence", "Dominant.score", color="Occurrence",label = "name", repel = TRUE, ylab="Driver dominant score",main="LUSC")

