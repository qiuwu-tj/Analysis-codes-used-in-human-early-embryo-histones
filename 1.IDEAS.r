library(ggsci)
library(ggplot2)
library(reshape2)
ncol = c(pal_npg("nrc")(8),"grey")

############################################################## 0. Quality control ##############################################################
peak_cor <- read.delim('0.qualityMetric/peak_cor_N.txt',header=F);colnames(peak_cor)<-c("stage","histone","cor")
promoter_cor <- read.delim('0.qualityMetric/promoter_cor_N.txt',header=F);colnames(promoter_cor)<-c("stage","histone","cor")
peak_cor$stage <- factor(peak_cor$stage,unique(as.character(peak_cor$stage)))
promoter_cor$stage <- factor(promoter_cor$stage,unique(as.character(promoter_cor$stage)))

my_colors <- RColorBrewer::brewer.pal(8, "Blues")[3:7]

p <- ggplot(peak_cor, aes(x=stage, y=cor, fill=stage)) + 
     geom_bar(position="dodge", stat="identity") +
     scale_fill_manual(values=my_colors) + 
     scale_y_continuous(limits = c(0,1)) + 
     xlab("Stage")+ylab("Mean pearson correlation of replicates")+labs(title="Replicates correlation on peaks")+
     facet_wrap(~histone) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank(),
           legend.position = "none")
ggsave("peak_cor_N.pdf",p,width=4.5,height=3.5)

p <- ggplot(promoter_cor, aes(x=stage, y=cor, fill=stage)) + 
     geom_bar(position="dodge", stat="identity") + 
     scale_fill_manual(values=my_colors) + 
     scale_y_continuous(limits = c(0,1))+
     xlab("Stage")+ylab("Mean pearson correlation of replicates")+labs(title="Replicates correlation on promoters")+
     facet_wrap(~histone) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank(),
           legend.position = "none")
ggsave("promoter_cor_N.pdf",p,width=4.5,height=3.5)

############################################################## 1. Segment number and features ##############################################################
###### read segment file
IDEAS_segment <- read.delim('1.IDEAS/HPCOS_IDEAS_output/HPCOS.state',header=T,sep=" ",check.names=F)
IDEAS_segment <- cbind(IDEAS_segment, IDEAS_segment[,4]-IDEAS_segment[,3]);colnames(IDEAS_segment)[16] <- "Length"
IDEAS_segment <- IDEAS_segment[,c(1:4,11:12,5:10,13:16)]

IDEAS_segment_rename <- IDEAS_segment[,5:14]
IDEAS_segment_rename[IDEAS_segment_rename==0] <- 'Nonmarked';
IDEAS_segment_rename[IDEAS_segment_rename==1] <- 'WeakH3K9me3';
IDEAS_segment_rename[IDEAS_segment_rename==2] <- 'WeakH3K27me3'
IDEAS_segment_rename[IDEAS_segment_rename==3] <- 'WeakH3K4me3';
IDEAS_segment_rename[IDEAS_segment_rename==4] <- 'StrongH3K9me3';
IDEAS_segment_rename[IDEAS_segment_rename==5] <- 'Bivalent'
IDEAS_segment_rename[IDEAS_segment_rename==6] <- 'StrongH3K27me3'
IDEAS_segment_rename[IDEAS_segment_rename==7] <- 'StrongH3K4me3'
# IDEAS_segment_rename[IDEAS_segment_rename==8] <- 'Trivalent'
IDEAS_segment <- cbind(IDEAS_segment[,1:4],IDEAS_segment_rename,IDEAS_segment[,15:16])

###### read value file
IDEAS_value <- read.delim('../1.IDEAS/HPCOS_IDEAS_output/HPCOS_IDEAS_input_NBP.bedgraph',header=F)

IDEAS_value <- read.delim('1.IDEAS/HPCOS_IDEAS_output/HPCOS_IDEAS_input_NBP.bedgraph',header=F)
colnames(IDEAS_value) <- c("CHR","POSst","POSed",
                         "MIIoocyte_N_H3K4me3","MIIoocyte_N_H3K27me3","MIIoocyte_N_H3K9me3","MIIoocyte_P_H3K4me3","MIIoocyte_P_H3K27me3","MIIoocyte_P_H3K9me3",
                         "4cell_N_H3K4me3","4cell_N_H3K27me3","4cell_N_H3K9me3","4cell_P_H3K4me3","4cell_P_H3K27me3","4cell_P_H3K9me3",
                         "8cell_N_H3K4me3","8cell_N_H3K27me3","8cell_N_H3K9me3","8cell_P_H3K4me3","8cell_P_H3K27me3","8cell_P_H3K9me3",
                           "ICM_N_H3K4me3","ICM_N_H3K27me3","ICM_N_H3K9me3","ICM_P_H3K4me3","ICM_P_H3K27me3","ICM_P_H3K9me3",
                           "TE_N_H3K4me3","TE_N_H3K27me3","TE_N_H3K9me3","TE_P_H3K4me3","TE_P_H3K27me3","TE_P_H3K9me3")
IDEAS_combined <- cbind(IDEAS_segment, IDEAS_value[,4:ncol(IDEAS_value)]) ####IDEAS_combined used for histone_DE.r

ONS <- 5; ONV <- 17:19; OPS <- 6; OPV <- 20:22
FNS <- 7; FNV <- 23:25; FPS <- 8; FPV <- 26:28
ENS <- 9; ENV <- 29:31; EPS <- 10; EPV <- 32:34
INS <- 11; INV <- 35:37; IPS <- 12; IPV <- 38:40
TNS <- 13; TNV <- 41:43 ; TPS <- 14; TPV <- 44:46

IDEAS_MII_N <- IDEAS_combined[,c(ONS,ONV)];colnames(IDEAS_MII_N) <- c("Segment","H3K4me3","H3K27me3","H3K9me3");IDEAS_MII_N <- melt(IDEAS_MII_N)
IDEAS_4cell_N <- IDEAS_combined[,c(FNS,FNV)];colnames(IDEAS_4cell_N) <- c("Segment","H3K4me3","H3K27me3","H3K9me3");IDEAS_4cell_N <- melt(IDEAS_4cell_N)
IDEAS_8cell_N <- IDEAS_combined[,c(ENS,ENV)];colnames(IDEAS_8cell_N) <- c("Segment","H3K4me3","H3K27me3","H3K9me3");IDEAS_8cell_N <- melt(IDEAS_8cell_N)
IDEAS_ICM_N <- IDEAS_combined[,c(INS,INV)];colnames(IDEAS_ICM_N) <- c("Segment","H3K4me3","H3K27me3","H3K9me3");IDEAS_ICM_N <- melt(IDEAS_ICM_N)
IDEAS_TE_N <- IDEAS_combined[,c(TNS,TNV)];colnames(IDEAS_TE_N) <- c("Segment","H3K4me3","H3K27me3","H3K9me3");IDEAS_TE_N <- melt(IDEAS_TE_N)

IDEAS_MII_N$Segment <- factor(IDEAS_MII_N$Segment, c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3'))
IDEAS_4cell_N$Segment <- factor(IDEAS_4cell_N$Segment, c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3'))
IDEAS_8cell_N$Segment <- factor(IDEAS_8cell_N$Segment, c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3'))
IDEAS_ICM_N$Segment <- factor(IDEAS_ICM_N$Segment, c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3'))
IDEAS_TE_N$Segment <- factor(IDEAS_TE_N$Segment, c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3'))

###### check the histone modification level of each segment, use 8cell_N, ICM_N as an example
p <- ggplot(IDEAS_MII_N, aes(x=variable, y=value, fill=variable)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,6))+
     xlab("Histone Mark")+ylab("Normalized ChIP-seq Signal")+labs(title="MII_N")+
     facet_wrap(~Segment) + scale_fill_manual(values=ncol[1:3]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_value_MII_N.pdf",p,width=6,height=6)
p <- ggplot(IDEAS_8cell_N, aes(x=variable, y=value, fill=variable)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,6))+
     xlab("Histone Mark")+ylab("Normalized ChIP-seq Signal")+labs(title="8cell_N")+
     facet_wrap(~Segment) + scale_fill_manual(values=ncol[1:3]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_value_8cell_N.pdf",p,width=6,height=6)
p <- ggplot(IDEAS_TE_N, aes(x=variable, y=value, fill=variable)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,6))+
     xlab("Histone Mark")+ylab("Normalized ChIP-seq Signal")+labs(title="TE_N")+
     facet_wrap(~Segment) + scale_fill_manual(values=ncol[1:3]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_value_TE_N.pdf",p,width=6,height=6)
p <- ggplot(IDEAS_4cell_N, aes(x=variable, y=value, fill=variable)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,6))+
     xlab("Histone Mark")+ylab("Normalized ChIP-seq Signal")+labs(title="4cell_N")+
     facet_wrap(~Segment) + scale_fill_manual(values=ncol[1:3]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_value_4cell_N.pdf",p,width=6,height=6)
p <- ggplot(IDEAS_ICM_N, aes(x=variable, y=value, fill=variable)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,6))+
     xlab("Histone Mark")+ylab("Normalized ChIP-seq Signal")+labs(title="ICM_N")+
     facet_wrap(~Segment) + scale_fill_manual(values=ncol[1:3]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_value_ICM_N.pdf",p,width=6,height=6)

###### segment numbers
H3K4me3_base <- as.numeric(quantile(IDEAS_combined[which(IDEAS_combined[,ONS]=="WeakH3K4me3"),ONV[1]],0.05)) # 1% percentile MII_N Weak H3K4me3 as baseline 0.7345003
H3K27me3_base <- as.numeric(quantile(IDEAS_combined[which(IDEAS_combined[,ONS]=="WeakH3K27me3"),ONV[2]],0.05)) # 1% percentile MII_N Weak H3K27me3 as baseline 0.5852117
H3K9me3_base <- as.numeric(quantile(IDEAS_combined[which(IDEAS_combined[,ONS]=="WeakH3K9me3"),ONV[3]],0.05)) # 1% percentile MII_N H3K9me3 as baseline

H3K4me3_base <-  max(IDEAS_combined[which(IDEAS_combined[,ONS]=="Nonmarked"),ONV[1]]) # 0.7345003
H3K27me3_base <-  max(IDEAS_combined[which(IDEAS_combined[,ONS]=="Nonmarked"),ONV[2]]) # 0.5852117
H3K9me3_base <- max(IDEAS_combined[which(IDEAS_combined[,ONS]=="Nonmarked"),ONV[3]]) # 0.7748691

H3K4me3_fraction <- cbind(table(IDEAS_combined[which(IDEAS_combined[,ONV[1]]>H3K4me3_base),ONS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,OPV[1]]>H3K4me3_base),OPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,FNV[1]]>H3K4me3_base),FNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                             table(IDEAS_combined[which(IDEAS_combined[,FPV[1]]>H3K4me3_base),FPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,ENV[1]]>H3K4me3_base),ENS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,EPV[1]]>H3K4me3_base),EPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,INV[1]]>H3K4me3_base),INS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,IPV[1]]>H3K4me3_base),IPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,TNV[1]]>H3K4me3_base),TNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,TPV[1]]>H3K4me3_base),TPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')])/nrow(IDEAS_combined);colnames(H3K4me3_fraction) <- c("MII_N","MII_P","4cell_N","4cell_P","8cell_N","8cell_P","ICM_N","ICM_P","TE_N","TE_P")
H3K27me3_fraction <- cbind(table(IDEAS_combined[which(IDEAS_combined[,ONV[2]]>H3K27me3_base),ONS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,OPV[2]]>H3K27me3_base),OPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,FNV[2]]>H3K27me3_base),FNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,FPV[2]]>H3K27me3_base),FPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,ENV[2]]>H3K27me3_base),ENS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,EPV[2]]>H3K27me3_base),EPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,INV[2]]>H3K27me3_base),INS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,IPV[2]]>H3K27me3_base),IPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,TNV[2]]>H3K27me3_base),TNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                           table(IDEAS_combined[which(IDEAS_combined[,TPV[2]]>H3K27me3_base),TPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')])/nrow(IDEAS_combined);colnames(H3K27me3_fraction) <- c("MII_N","MII_P","4cell_N","4cell_P","8cell_N","8cell_P","ICM_N","ICM_P","TE_N","TE_P")
H3K9me3_fraction <- cbind(table(IDEAS_combined[which(IDEAS_combined[,ONV[3]]>H3K9me3_base),ONS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,OPV[3]]>H3K9me3_base),OPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,FNV[3]]>H3K9me3_base),FNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,FPV[3]]>H3K9me3_base),FPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,ENV[3]]>H3K9me3_base),ENS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,EPV[3]]>H3K9me3_base),EPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,INV[3]]>H3K9me3_base),INS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,IPV[3]]>H3K9me3_base),IPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,TNV[3]]>H3K9me3_base),TNS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')],
                          table(IDEAS_combined[which(IDEAS_combined[,TPV[3]]>H3K9me3_base),TPS])[c('Bivalent','Nonmarked','StrongH3K4me3','StrongH3K27me3','StrongH3K9me3','WeakH3K4me3','WeakH3K27me3','WeakH3K9me3')])/nrow(IDEAS_combined);colnames(H3K9me3_fraction) <- c("MII_N","MII_P","4cell_N","4cell_P","8cell_N","8cell_P","ICM_N","ICM_P","TE_N","TE_P")
for(i in c(ONV[1],OPV[1],FNV[1],FPV[1],ENV[1],EPV[1],INV[1],IPV[1],TNV[1],TPV[1]))
{
  write.table(IDEAS_combined[which(IDEAS_combined[,i]>H3K4me3_base),c(2:4,i,i+1,i+2)],paste0("H3K4me3_",i),sep="\t",quote=F,col.names=F,row.names=F)
  write.table(IDEAS_combined[which(IDEAS_combined[,i+1]>H3K27me3_base),c(2:4,i,i+1,i+2)],paste0("H3K27me3_",i),sep="\t",quote=F,col.names=F,row.names=F)
  write.table(IDEAS_combined[which(IDEAS_combined[,i+2]>H3K9me3_base),c(2:4,i,i+1,i+2)],paste0("H3K9me3_",i),sep="\t",quote=F,col.names=F,row.names=F)
}

IDEAS_fraction <- data.frame(rbind(melt(H3K4me3_fraction), melt(H3K27me3_fraction), melt(H3K9me3_fraction)), 
                             histone = c(rep("H3K4me3",nrow(melt(H3K4me3_fraction))),rep("H3K27me3",nrow(melt(H3K27me3_fraction))),rep("H3K9me3",nrow(melt(H3K9me3_fraction)))))
colnames(IDEAS_fraction) <- c("segment","stage","value","histone")
IDEAS_fraction$histone <- factor(IDEAS_fraction$histone, c('H3K4me3','H3K27me3','H3K9me3'))
IDEAS_fraction$segment <- factor(IDEAS_fraction$segment, c('Bivalent','StrongH3K4me3','WeakH3K4me3','StrongH3K27me3','WeakH3K27me3','StrongH3K9me3','WeakH3K9me3','Nonmarked'))
IDEAS_fraction <- na.omit(IDEAS_fraction)
IDEAS_fraction_N <- IDEAS_fraction[which(IDEAS_fraction[,2]=='MII_N'|IDEAS_fraction[,2]=='4cell_N'|IDEAS_fraction[,2]=='8cell_N'|IDEAS_fraction[,2]=='ICM_N'|IDEAS_fraction[,2]=='TE_N'),]

p <- ggplot(IDEAS_fraction_N, aes(x=stage, y=value, fill=segment)) + 
     geom_bar(position="stack", stat="identity") +
     scale_y_continuous(limits = c(0,0.05)) +
     xlab("Stage")+ylab("Fraction of the Genome-wide Bins")+labs(title="Genome-wide fraction of histone segments")+
     facet_wrap(~histone) + scale_fill_manual(values=ncol[c(2,8,1,4,6,3,7,9)]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_fraction_N.pdf",p,width=6,height=3.5)
p <- ggplot(IDEAS_fraction, aes(x=stage, y=value, fill=segment)) + 
     geom_bar(position="stack", stat="identity") +
     scale_y_continuous(limits = c(0,0.05)) + 
     xlab("Stage")+ylab("Fraction of the Genome-wide Bins")+labs(title="Genome-wide fraction of histone segments")+
     facet_wrap(~histone) + scale_fill_manual(values=ncol[c(2,8,1,4,6,3,7,9)]) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_genome_fraction_NP.pdf",p,width=10,height=3.5)

###### segment genomic distribution
segment_overlap <- read.table('1.IDEAS/IDEAS_state_enrichment_N.txt')
segment_overlap[,1:8] <- segment_overlap[,1:8]/segment_overlap[,9]
bg <- c(733269268,289459557,6744452833,4481639,409839127,663376365,277960976,38948936)/3088286376 # promoter, exon, intron, SVA, SINE, LINE, LTR, simple_repeats
segment_odds <- t(rbind(segment_overlap[,1:8], bg))
segment_odds[,1:40] <- segment_odds[,1:40]/segment_odds[,41];rownames(segment_odds)<-c("Promoter","Exon","Intron","SVA","SINE","LINE","LTR","SimpleRepeats")
segment_odds <- melt(log2(t(segment_odds[,1:40])))

p <- ggplot(segment_odds, aes(x=Var2, y=value, fill=value)) + 
     geom_bar(position="dodge", stat="identity") +
     scale_y_continuous(limits = c(-6,5))+
     xlab("Genomic elements")+ylab("log2(Observed/Expected)")+labs(title="Enrichment of histone segments")+
     facet_wrap(~Var1, ncol=8) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank(),
           legend.position = "none",
           strip.text.x = element_text(size = 6.5, color = "blue", face = "bold.italic"))
ggsave("IDEAS_genome_enrichment.pdf",p,width=10.5,height=9)

#####################pie chart of segment genomic distribution###################
segment_overlap_length <- read.table('1.IDEAS/IDEAS_state_enrichment_N_RepeatsSub.txt')
segment_overlap_K9 <- segment_overlap_length[grepl('H3K9me3$',rownames(segment_overlap_length)),]
segment_overlap_K9 <- segment_overlap_K9[,c(9:12,14:15,18:19,21)]
colnames(segment_overlap_K9) <- c('ERV1','ERVK','ERVL','ERVL_MaLR','L1','L2','Alu','MIR','total')

segment_overlap_K9$others <- segment_overlap_K9[,9]-apply(segment_overlap_K9[,1:8],1,sum)

segment_overlap_K9_all <- rbind(segment_overlap_K9[1,]+segment_overlap_K9[2,],segment_overlap_K9[3,]+segment_overlap_K9[4,],segment_overlap_K9[5,]+segment_overlap_K9[6,],segment_overlap_K9[7,]+segment_overlap_K9[8,],segment_overlap_K9[9,]+segment_overlap_K9[10,])
segment_overlap_K9_strong <- segment_overlap_K9[grepl('StrongH3K9me3$',rownames(segment_overlap_K9)),]

segment_overlap_K9_p <- segment_overlap_K9_all[,c(1:8,10)]/segment_overlap_K9_all[,9]
segment_overlap_StrongK9_p <- segment_overlap_K9_strong[,c(1:8,10)]/segment_overlap_K9_strong[,9]

my_col <- c('#E64B35FF','#DC0000B2','#E64B35B2','#F39B7FFF','#00A087FF','#00A087B2','#8491B4FF','#8491B4B2','lightgrey')

data <- data.frame(group=colnames(segment_overlap_StrongK9_p),value=t(segment_overlap_StrongK9_p[1,]))
data$group=factor(data$group, levels=rev(data$group), order=T)
p <- ggplot(data, aes(x="", y= MII_N_StrongH3K9me3, fill= group))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="MII")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=my_col,breaks = data$group, labels = as.character(data$group))+
  geom_text(aes(x = 1.7, y = cumsum(data$MII_N_StrongH3K9me3)-data$MII_N_StrongH3K9me3/2 , 
                label =paste(as.character(round(data$MII_N_StrongH3K9me3*100,2)),"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("MII_StrongH3K9me3_enrichment_repeats_RepeatsSub.pdf",p,width=5,height=5,useDingbats=FALSE)

data <- data.frame(group=colnames(segment_overlap_StrongK9_p),value=t(segment_overlap_StrongK9_p[2,]))
data$group=factor(data$group, levels=rev(data$group), order=T)
p <- ggplot(data, aes(x="", y= X4cell_N_StrongH3K9me3, fill= group))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="4-cell")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=my_col,breaks = data$group, labels = as.character(data$group))+
  geom_text(aes(x = 1.7, y = cumsum(data$X4cell_N_StrongH3K9me3)-data$X4cell_N_StrongH3K9me3/2 , label =paste(as.character(round(data$X4cell_N_StrongH3K9me3*100,2)),"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("4cell_StrongH3K9me3_enrichment_repeats_RepeatsSub.pdf",p,width=5,height=5,useDingbats=FALSE)

data <- data.frame(group=colnames(segment_overlap_StrongK9_p),value=t(segment_overlap_StrongK9_p[3,]))
data$group=factor(data$group, levels=rev(data$group), order=T)
p <- ggplot(data, aes(x="", y= X8cell_N_StrongH3K9me3, fill= group))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="8-cell")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=my_col,breaks = data$group, labels = as.character(data$group))+
  geom_text(aes(x = 1.7, y = cumsum(data$X8cell_N_StrongH3K9me3)-data$X8cell_N_StrongH3K9me3/2 , label =paste(as.character(round(data$X8cell_N_StrongH3K9me3*100,2)),"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("8cell_StrongH3K9me3_enrichment_repeats_RepeatsSub.pdf",p,width=5,height=5,useDingbats=FALSE)

data <- data.frame(group=colnames(segment_overlap_StrongK9_p),value=t(segment_overlap_StrongK9_p[4,]))
data$group=factor(data$group, levels=rev(data$group), order=T)
p <- ggplot(data, aes(x="", y= ICM_N_StrongH3K9me3, fill= group))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="ICM")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=my_col,breaks = data$group, labels = as.character(data$group))+
  geom_text(aes(x = 1.7, y = cumsum(data$ICM_N_StrongH3K9me3)-data$ICM_N_StrongH3K9me3/2 , label =paste(as.character(round(data$ICM_N_StrongH3K9me3*100,2)),"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("ICM_StrongH3K9me3_enrichment_repeats_RepeatsSub.pdf",p,width=5,height=5,useDingbats=FALSE)

data <- data.frame(group=colnames(segment_overlap_StrongK9_p),value=t(segment_overlap_StrongK9_p[5,]))
data$group=factor(data$group, levels=rev(data$group), order=T)
p <- ggplot(data, aes(x="", y= TE_N_StrongH3K9me3, fill= group))+
  geom_bar(stat="identity",width=1)+
  coord_polar(theta="y")+
  labs(x="", y="",title="TE")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  scale_fill_manual(values=my_col,breaks = data$group, labels = as.character(data$group))+
  geom_text(aes(x = 1.7, y = cumsum(data$TE_N_StrongH3K9me3)-data$TE_N_StrongH3K9me3/2 , label =paste(as.character(round(data$TE_N_StrongH3K9me3*100,2)),"%",sep="")), show.legend = FALSE, color="black")+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line= element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("TE_StrongH3K9me3_enrichment_repeats_RepeatsSub.pdf",p,width=5,height=5,useDingbats=FALSE)

##### segment CpG density
segment_CpG <- read.delim('1.IDEAS/HPCOS_IDEAS_output/HPCOS.CpG',header=F)
IDEAS_combined <- cbind(IDEAS_combined, segment_CpG[,4]);colnames(IDEAS_combined)[ncol(IDEAS_combined)] <- "CpG"
segment_CpG <- na.omit(rbind(data.frame(segment=IDEAS_combined[,ONS],stage="MIIOocyte_N",value=IDEAS_combined[,ncol(IDEAS_combined)]),
#                             data.frame(segment=IDEAS_combined[,OPS],stage="MIIOocyte_P",value=IDEAS_combined[,ncol(IDEAS_combined)]),
                             data.frame(segment=IDEAS_combined[,FNS],stage="4cell_N",value=IDEAS_combined[,ncol(IDEAS_combined)]), 
#                             data.frame(segment=IDEAS_combined[,FPS],stage="4cell_P",value=IDEAS_combined[,ncol(IDEAS_combined)]), 
                             data.frame(segment=IDEAS_combined[,ENS],stage="8cell_N",value=IDEAS_combined[,ncol(IDEAS_combined)]),
#                             data.frame(segment=IDEAS_combined[,EPS],stage="8cell_P",value=IDEAS_combined[,ncol(IDEAS_combined)]), 
                             data.frame(segment=IDEAS_combined[,INS],stage="ICM_N",value=IDEAS_combined[,ncol(IDEAS_combined)]), 
#                             data.frame(segment=IDEAS_combined[,IPS],stage="ICM_P",value=IDEAS_combined[,ncol(IDEAS_combined)]),
                             data.frame(segment=IDEAS_combined[,TNS],stage="TE_N",value=IDEAS_combined[,ncol(IDEAS_combined)])))
#                             data.frame(segment=IDEAS_combined[,TPS],stage="TE_P",value=IDEAS_combined[,ncol(IDEAS_combined)]))
segment_CpG[,1] <- factor(segment_CpG[,1], c('Bivalent','StrongH3K4me3','WeakH3K4me3','StrongH3K27me3','WeakH3K27me3','StrongH3K9me3','WeakH3K9me3','Nonmarked'))
segment_CpG[,2] <- factor(segment_CpG[,2], c("MIIOocyte_N","4cell_N","8cell_N","ICM_N","TE_N"))

p <- ggplot(segment_CpG, aes(x=segment, y=value, fill=value)) + 
     geom_boxplot(outlier.shape = NA) +
     scale_y_continuous(limits = c(0,10))+
     xlab("Segments")+ylab("CpG density in 100bp tiles")+labs(title="CpG density of histone segments")+
     facet_wrap(~stage,ncol=5) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
           panel.grid.major = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.ticks = element_blank())
ggsave("IDEAS_CpG_density.pdf",p,width=10,height=3.5)


















