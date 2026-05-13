####################################################
# Include codes to re-generate plots in Figure 5 and Supplementary Figure 8
# Necessary data can be found in ".../data/" directory
####################################################

# Plot Supplementary Figure 8a
library(enrichplot)
em<-readRDS("data/FigS8_GSEA_res.rds")
data<-em@result
colors<-c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
gseaplot2(em, geneSetID = data$ID[which(data$p.adjust<1e-6)],color = colors)

# Plot Supplementary Figure 8b
library(ggpubr)
SF_score<-readRDS("data/Fig5_regulatory_score.rds")
data<-SF_score
data$abs_FC<-abs(data$logFC_bulk)*abs(data$logFC_sc)
data<-data[order(data$abs_FC,decreasing = T),]
de_sf<-data$SF_gene[1:10]
ggscatter(data, x="logFC_sc", y="logFC_bulk", size = "abs_FC",color = "logFC_bulk", label = "SF_gene", label.select = de_sf)+
  scale_color_gradient2(low = "#191970", high = "#8B0000")

# Plot Figure 5b
dse_sf<-SF_score$SF_gene[which(SF_score$rank<=5)]
SF_score$log_score<-log2(SF_score$regulatory_score+1)
ggscatter(SF_score,x="rank",y="regulatory_score",color = "black",fill = "log_score",shape = 21,size = "regulatory_score",
          legend="right", label = "SF_gene", label.select = dse_sf, repel = T)+
  scale_fill_gradient2(high = "#FAEEB7",mid = "#DF6863",low="#140E36",midpoint=1)

# Plot Figure 5c
library(Seurat)
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data_new.rds")
sn_sum$celltype_subtype<-factor(sn_sum$celltype_subtype,levels=c("Oligodendrocyte","Neuron","Myeloid","Lymphocyte","Endothelial","Astrocyte","Tumor.Cycle","Tumor.OPC","Tumor.AC","Tumor.NPC","Tumor.MES","Tumor.Stem","Tumor.Fiber"))
cnv_res<-readRDS("data/infercnv_cn_50k.rds")
library(scCustomize)
sn_sum.sub<-sn_sum[,row.names(cnv_res)]
sn_sum.sub$chr7_cnv<-as.numeric(cnv_res$chr7)
sn_sum.sub$chr10_cnv<-as.numeric(cnv_res$chr10)
Stacked_VlnPlot(seurat_object = sn_sum.sub, features = c("chr7_cnv","chr10_cnv"), x_lab_rotate = TRUE, group.by = "celltype_subtype",
                colors_use = c("#91D1C2FF","#8491B4FF","#F39B7FFF","#3C5488FF","#00A087FF","#4DBBD5FF","#AC9C88","#C82B1E","#CD7692","#4A885C","#E7AC53","#7A634B","#8E39C5"))

library(viridis)
DotPlot(sn_sum,features = c("RBFOX3","CELF4","SRRM4","CELF5","RBFOX1"),group.by = "celltype_subtype") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  coord_flip()+rotate_x_text(angle = 45)

# Plot Figure 5e-5f
psi_res<-read.csv("data/RT_PCR_splic_res.csv")
ggviolin(psi_res,x="Group",y="PSI",add = "boxplot",fill = "Group",palette = c("#BED2DF","#D98986"), facet.by = "Event")+stat_compare_means(label = "p.format",method = "t.test")

# Plot Figure 5g-5h, Supplementary Figure 8d
pcr_res<-read.csv("data/RT_qPCR_res.csv")
data<-pcr_res[which(pcr_res$Gene=="RBFOX3"),]
ggbarplot(data,x="Group",y="Exp",add = "mean_se",fill = "Group",palette = c("#BED2DF","#D98986"))+stat_compare_means(label = "p.format",method = "t.test")
data<-pcr_res[which(pcr_res$Gene=="SOX11"),]
ggbarplot(data,x="Group",y="Exp",add = "mean_se",fill = "Group",palette = c("#BED2DF","#D98986"))+stat_compare_means(label = "p.format",method = "t.test")
data<-pcr_res[which(pcr_res$Gene=="SNAP25"),]
ggbarplot(data,x="Group",y="Exp",add = "mean_se",fill = "Group",palette = c("#BED2DF","#D98986"))+stat_compare_means(label = "p.format",method = "t.test")

# Plot Supplementary Figure 8f
tmz_dose<-read.csv("/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/experiment/Drug_response_curves.csv", header = T)

library(drc)
p <- ggplot(tmz_dose, aes(x = TMZ, y = Prop, colour = Group)) +
  geom_point() + 
  geom_smooth(
    method = "drm", 
    method.args = list(fct = LL.3()), # 改用 3 参数模型，默认下限为 0
    se = FALSE
  ) +
  coord_trans(x = "log10")+scale_x_continuous(breaks = c(0.001, 0.01, 0.1, 1, 10)) 

p+labs(
  x = "TMZ Concentration (\u00B5M)", 
  y = "Relative Cell Viability",     
  colour = "Group"                   
) +scale_color_manual(
  values = c(
    "NC" = "#0072B5", # 深蓝色 (你可以改成 "black" 或 "blue" 等)
    "OE" = "#BC3C29"  # 深红色 (你可以改成 "red" 等)
  )
) +theme_classic()
