####################################################
# Include codes to re-generate plots in Figure 1 and Supplementary Figure 1
# Necessary data can be found in ".../data/" directory
####################################################

datapath=".../data/"
# Load results for differential splicing test
library(readxl)
res_se <- read_excel("data/Test results of splicing events for 102 GB patients.xlsx", sheet = "SE")
res_mxe <- read_excel("data/Test results of splicing events for 102 GB patients.xlsx", sheet = "MXE")
res_a53ss <- read_excel("data/Test results of splicing events for 102 GB patients.xlsx", sheet = "A53SS") 
res_sum<-rbind(res_se[,c("geneSymbol","Diff_RvsP","P_value.adj")],res_mxe[,c("geneSymbol","Diff_RvsP","P_value.adj")])

# Plot Figure 1b
library(EnhancedVolcano)
res_sum$geneSymbol2<-NA
res_sum$geneSymbol2[intersect(which(abs(res_sum$Diff_RvsP)>0.1),which(res_sum$P_value.adj<0.03))]<-res_sum$geneSymbol[intersect(which(abs(res_sum$Diff_RvsP)>0.1),which(res_sum$P_value.adj<0.03))]
EnhancedVolcano(res_sum,lab = res_sum$geneSymbol2,
                labSize = 3.5,
                x="Diff_RvsP",y="P_value.adj",
                pCutoff = 0.05,
                FCcutoff = 0.05,
                pointSize = 1,
                xlim = c(-0.21,0.21),
                ylim = c(0,2.2),
                colAlpha = 1,
                cutoffLineCol = '#C6C6C8',
                col = c("#C6C6C8","#5E856B","#5D6AB7","#C795A7"),
                gridlines.major = F,
                gridlines.minor = F,
                border = 'full',
                borderWidth = 0.6,
                borderColour = 'black')

# Plot Figure 1c
res_sum<-rbind(res_se[,c("geneSymbol","Diff_RvsP","P_value.adj")],res_mxe[,c("geneSymbol","Diff_RvsP","P_value.adj")],res_a53ss[,c("geneSymbol","Diff_RvsP","P_value.adj")])
gene<-unique(res_sum$geneSymbol[intersect(which(res_sum$P_value.adj<0.05),which(abs(res_sum$Diff_RvsP)>0.05))])
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
go_enrich_up<-enrichGO(gene = gene,
                       universe =unique(res_sum$geneSymbol),
                       OrgDb = org.Hs.eg.db, 
                       keyType = 'SYMBOL',
                       readable = F,
                       ont = "BP",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.1)
categorys <- go_enrich_up@result$Description[c(1,3,5,8,18,22,23,33,46,47)]
library(ggpubr)
barplot(go_enrich_up, showCategory=categorys, label_format = 60)+theme_pubclean()

# Plot Figure 1d
GB_psi_info <- readRDS("data/Fig1_dse_psi_info.rds")
data<-GB_psi_info[,c("Patient","condition","MSE_30744")]
colnames(data)<-c("Patient","condition","PSI")
library(reshape2)
data_short<-dcast(data,Patient~condition,value.var = "PSI")
data_short$Diff<-(data_short$Recurrent-data_short$Initial)
data_short<-data_short[-which(is.na(data_short$Diff)),]
data_short<-data_short[order(data_short$Diff,decreasing = T),]
data$group<-ifelse(data$Patient%in%data_short$Patient[1:8],"Top","Other")
cols<-c("Top"="#C71585","Other"="light grey","Initial"="#008B8B","Recurrent"="#7F17F5")
ggplot(data, aes(condition, PSI,color=condition,fill=condition)) +
  geom_boxplot(width=0.5) + scale_fill_manual(values = c("#F3FCFF","#E6E6F8")) +
  geom_point(aes(color = group),size=2) + scale_color_manual(values = cols) +
  geom_point(data=data[data$group == "Top",],color="#C71585")+
  geom_line(aes(group=Patient),colour="light grey",alpha=0.8) + 
  geom_line(data=data[data$group == "Top",],aes(group=Patient),color="#C71585")+ stat_compare_means(paired = T) +
  theme_classic()+labs(y = "Inclusion level of exon 9&10 on CD47") + theme(axis.text.x = element_text(size = 10,face = "bold"))  

# Plot Supplementary Figure 1a
data<-GB_psi_info[,c("Patient","condition","SE_219424")]
colnames(data)<-c("Patient","condition","PSI")
data_short<-dcast(data,Patient~condition,value.var = "PSI")
data_short$Diff<-(data_short$Recurrent-data_short$Initial)
data_short<-data_short[-which(is.na(data_short$Diff)),]
data_short<-data_short[order(data_short$Diff,decreasing = T),]
data$group<-ifelse(data$Patient%in%data_short$Patient[1:8],"Top","Other")
cols<-c("Top"="#C71585","Other"="light grey","Initial"="#008B8B","Recurrent"="#7F17F5")
ggplot(data, aes(condition, PSI,color=condition,fill=condition)) +
  geom_boxplot(width=0.5) + scale_fill_manual(values = c("#F3FCFF","#E6E6F8")) +
  geom_point(aes(color = group),size=2) + scale_color_manual(values = cols) +
  geom_point(data=data[data$group == "Top",],color="#C71585")+
  geom_line(aes(group=Patient),colour="light grey",alpha=0.8) + 
  geom_line(data=data[data$group == "Top",],aes(group=Patient),color="#C71585")+ stat_compare_means(paired = T) +
  theme_classic()+labs(y = "Inclusion level of exon 9 on MAP2") + theme(axis.text.x = element_text(size = 10,face = "bold"))  

# Plot Figure 1f
library(pheatmap)
psi_mtx <- t(GB_psi_info[,6:159])
for(i in colnames(psi_mtx)){psi_mtx[,i][is.na(psi_mtx[,i])] <- 0}
col_annotation<-GB_psi_info[,2:5]
row.names(col_annotation)<-GB_psi_info$SampleID
colors<-list(splice_grp=c(C1="#EE9235",C2="#4A68DA"),
             condition=c(Initial="#40E0D0",Recurrent="#8A00FF"),
             subtype=c(Neural="#FF6347",Proneural="#9ACD32",Mesenchymal="#4682B4", Classical="#C71585"))
set.seed(123)
pheatmap(psi_mtx,annotation_col = col_annotation[colnames(psi_mtx),2:4],show_colnames = F,show_rownames = F,
         cutree_cols = 2,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_colors = colors,clustering_method = "ward.D2")

# Plot Supplementary Figure 1c
data<-GB_psi_info[,c("Patient","condition","SE_86199")]
colnames(data)<-c("Patient","condition","PSI")
data_short<-dcast(data,Patient~condition,value.var = "PSI")
data_short$Diff<-(data_short$Recurrent-data_short$Initial)
data_short<-data_short[-which(is.na(data_short$Diff)),]
data_short<-data_short[order(data_short$Diff,decreasing = T),]
data<-data[which(data$Patient%in%data_short$Patient),]
data$group<-ifelse(data$Patient%in%GB_psi_info$Patient[intersect(which(GB_psi_info$splice_grp=="C1"),which(GB_psi_info$condition=="Recurrent"))],"C1","C2")
cols<-c("C1"="#C71585","C2"="light grey","Initial"="#008B8B","Recurrent"="#7F17F5")
ggplot(data, aes(condition, PSI,color=condition,fill=condition)) +
  geom_boxplot(width=0.5) + scale_fill_manual(values = c("#F3FCFF","#E6E6F8")) +
  geom_point(aes(color = group),size=2) + scale_color_manual(values = cols) +
  geom_point(data=data[data$group == "C1",],color="#C71585")+
  geom_line(aes(group=Patient),colour="light grey",alpha=0.8) + 
  geom_line(data=data[data$group == "C1",],aes(group=Patient),color="#C71585")+ stat_compare_means(paired = T) +
  theme_classic()+labs(y = "Inclusion level of exon 17 on MAP4K4") + theme(axis.text.x = element_text(size = 10,face = "bold"))  

# Plot Supplementary Figure 1e
data<-GB_psi_info[,c("Patient","condition","SE_154929")]
colnames(data)<-c("Patient","condition","PSI")
data_short<-dcast(data,Patient~condition,value.var = "PSI")
data_short$Diff<-(data_short$Recurrent-data_short$Initial)
data_short<-data_short[order(data_short$Diff,decreasing = T),]
data$group<-ifelse(data$Patient%in%GB_psi_info$Patient[intersect(which(GB_psi_info$splice_grp=="C1"),which(GB_psi_info$condition=="Recurrent"))],"C1","C2")
cols<-c("C1"="#C71585","C2"="light grey","Initial"="#008B8B","Recurrent"="#7F17F5")
ggplot(data, aes(condition, PSI,color=condition,fill=condition)) +
  geom_boxplot(width=0.5) + scale_fill_manual(values = c("#F3FCFF","#E6E6F8")) +
  geom_point(aes(color = group),size=2) + scale_color_manual(values = cols) +
  geom_point(data=data[data$group == "C1",],color="#C71585")+
  geom_line(aes(group=Patient),colour="light grey",alpha=0.8) + 
  geom_line(data=data[data$group == "C1",],aes(group=Patient),color="#C71585")+ stat_compare_means(paired = T) +
  theme_classic()+labs(y = "Inclusion level of exon 25 on ITGA6") + theme(axis.text.x = element_text(size = 10,face = "bold"))  

# Plot Figure 1g
library(ggalluvial)
data<-GB_psi_info[,c("Patient","condition","subtype")]
grp4_ssgsea_short<-dcast(data,Patient~condition,value.var = "subtype")
grp4_ssgsea_short$group<-"C2.NEU-"
grp4_ssgsea_short$group[which(grp4_ssgsea_short$Recurrent=="Neural")]<-"C2.NEU+"
grp4_ssgsea_short$group[which(grp4_ssgsea_short$Patient%in%GB_psi_info$Patient[intersect(which(GB_psi_info$splice_grp=="C1"),which(GB_psi_info$condition=="Recurrent"))])]<-"C1.NEU+"
plot_table<-as.data.frame(table(grp4_ssgsea_short$Initial,grp4_ssgsea_short$Recurrent,grp4_ssgsea_short$group))
colnames(plot_table)<-c("Subtype.I","Subtype.R","group","Count")
ggplot(plot_table,
       aes(axis1 = Subtype.I,
           axis2 = Subtype.R,
           y = Count)) +
  geom_alluvium(aes(fill = group)) +
  scale_fill_manual(values = c(`C1.NEU+` = "#EE9238", `C2.NEU+` = "#3069F6", `C2.NEU-` = "#6C7FAE"))+
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("Primary", "Recurrent"),
                   expand = c(.1, .1)) +
  theme_minimal()
