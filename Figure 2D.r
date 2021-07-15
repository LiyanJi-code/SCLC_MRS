### custom R code for the publication: Huaqiang Zhou, Yi Hu, Rongzhen Luo et al. "Multi-region exome sequencing reveals the intratumoral heterogeneity of surgically resected small cell lung cancer"
###################
### Load R packages
  library(gridExtra)
  library(plyr)
  library(ggplot2)
  library(data.table)
###################

###read in file, include at least three columns: Tumor_Sample_Barcode	Gene_symbol	CNV_status
cnv_mat <- read.table("~/sourcedata/SCLC_cnv.csv",header=T, stringsAsFactors=F, sep=",")

cnv_mat$flag <- paste(cnv_mat$Gene_symbol,cnv_mat$CNV_status,sep=" ")
cnv_mat <- cnv_mat[,c("Tumor_Sample_Barcode","flag")]
cnv_mat <- unique(cnv_mat)

patients <- unique(gsub("-.*$","",cnv_mat$Tumor_Sample_Barcode))

### subclonal CNV 
ITH_CNV <- NULL
for (i in patients){
    tmp <- cnv_mat[grepl(i,cnv_mat$Tumor_Sample_Barcode),]
    if (nrow(tmp)==0){
      ITH_CNV <- rbind(ITH_CNV,c(i,0,0,0))
      next
    }
    all <- unique(tmp$flag)
    cnv <- as.data.frame(table(tmp$flag))
    overlap <- cnv[cnv$Freq==3,]
    ITH <- (length(all)-dim(overlap)[1])/length(all)
    ITH_CNV <- rbind(ITH_CNV,c(i,nrow(overlap),length(all)-nrow(overlap),ITH))
}
colnames(ITH_CNV) <- c("patientID","trunk","branch","ITH")
ITH_CNV <- data.frame(ITH_CNV,stringsAsFactors = F)
ITH_CNV$ITH <- as.numeric(ITH_CNV$ITH)

### matrix for plot Figure 2D  
ITH_CNV <- ITH_CNV[order(ITH_CNV$ITH),]
ITH_CNV$patientID <- factor(ITH_CNV$patientID,levels = ITH_CNV$patientID, ordered = T)
ITH_CNV_plot <- melt(ITH_CNV,id=c("patientID","ITH"),measure.vars=c("trunk","branch"))
ITH_CNV_plot$value <- as.numeric(ITH_CNV_plot$value)

df_cumsum <- ddply(ITH_CNV_plot, .(patientID),
                     transform, sum=sum(value))
 
df_cumsum$variable <- factor(df_cumsum$variable, levels=rev(levels(df_cumsum$variable)))
  	
CNV_ITH.p1 <- ggplot(df_cumsum, aes(x=patientID,y=value,fill=variable))+geom_bar(stat="identity",position = 'stack')+ggpubr::theme_pubr()+
    theme(axis.text.x = element_text(angle=60,hjust=1))+labs(title = "",x="",y="Number of CNV per patient",fill="Group")+ 
	ggsci::scale_fill_npg()+ theme(legend.position="right")
	
df_cumsum$percentage <- ifelse(df_cumsum$variable == "trunk", 1- df_cumsum$ITH, df_cumsum$ITH)

CNV_ITH.p2 <- ggplot() + geom_bar(aes(x = patientID, y=percentage, fill = variable), stat="identity",data = df_cumsum)+ggpubr::theme_pubr()+
    theme(axis.text.x = element_text(angle=60,hjust=1))+labs(title = "",x="", y="Percentage of subclonal CNV",fill="Group")+ 
	ggsci::scale_fill_npg() + theme(legend.position="right")
    