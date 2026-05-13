####################################################
# Include codes to re-generate plots in Figure 3 and Supplementary Figure 3-6
# Necessary data can be found in ".../data/" directory
####################################################

# Plot Figure 3b
sn_meta<-read.csv("data/snRNA-seq.GBM.metadata.csv")
library(ggpubr)
data<-as.data.frame(table(sn_meta$splice_type,sn_meta$Stage,sn_meta$grp4_ssgsea))
colnames(data)<-c("group","Stage","subtype","count")
ggbarplot(data,x="Stage",y="count",facet.by = "group",fill = "subtype", color = "subtype",orientation = "horizontal",
          palette = c("#C71585","#4682B4","#FF6347","#9ACD32"), ncol = 1, label = TRUE, lab.pos = "in", lab.col = "white")+
  theme_classic2()

# Plot Figure 3c
library(reshape2)
grp4_ssgsea_short<-dcast(sn_meta,Pair.+grp_3type~Stage,value.var = "grp4_ssgsea")

library(ggalluvial)
plot_table<-as.data.frame(table(grp4_ssgsea_short$Primary,grp4_ssgsea_short$Recurrent,grp4_ssgsea_short$grp_3type))
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

# Plot Supplementary Figure 3b
load("data/FigS3_sn_splicing.RData")
library(pheatmap)
colors<-list(group=c(C1="#EE9235",C2="#4A68DA"),condition=c(I="#40E0D0",R="#8A00FF"),
             subtype=c(Neural="#FF6347",Proneural="#9ACD32",Mesenchymal="#4682B4", Classical="#C71585"))
data0<-nofil_psi_dse
for(i in colnames(data0)){data0[,i][is.na(data0[,i])] <- 0}
pheatmap(data0,show_rownames = F,annotation_col = anno_sn,scale = "row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_colors = colors,cutree_cols = 2,border_color = NA)

# Plot Supplementary Figure 3c
sn_meta<-sn_meta[order(sn_meta$SampleID),]
sn_meta<-sn_meta[order(sn_meta$Stage),]
sn_meta$grp_3type<-factor(sn_meta$grp_3type, levels = c("C1.NEU+", "C2.NEU+", "C2.NEU-"))
ggpaired(sn_meta,x="Stage",y="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",facet.by = "grp_3type",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",paired = T)+ theme_classic()+ ylab("Enrichment score of MES markers")
ggpaired(sn_meta,x="Stage",y="VERHAAK_GLIOBLASTOMA_NEURAL",facet.by = "grp_3type",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",paired = T)+ theme_classic()+ ylab("Enrichment score of NEU markers")

# Plot Supplementary Figure 3d
tpm_sn<-readRDS("data/snRNA_pseudo_bulk_tpm.rds")
data<-sn_meta
data$exp<-tpm_sn["EGR4", data$SampleID]
ggpaired(data,x="Stage",y="exp",facet.by = "grp_3type",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",paired = T)+ theme_classic()+ ylab("Normalized TPM of EGR4")
data$exp<-tpm_sn["ANXA3", data$SampleID]
ggpaired(data,x="Stage",y="exp",facet.by = "grp_3type",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",paired = T)+ theme_classic()+ ylab("Normalized TPM of ANXA3")
data$exp<-tpm_sn["MLH1", data$SampleID]
ggpaired(data,x="Stage",y="exp",facet.by = "grp_3type",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",paired = T)+ theme_classic()+ ylab("Normalized TPM of MLH1")

# Plot Figure 3d
library(Seurat)
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data_new.rds")
DimPlot(object = sn_sum, reduction = 'umap',group.by = "condition",cols = c("#73DDD0","#7F17F5"))
DimPlot(sn_sum,reduction = 'umap',group.by = "celltype_subtype",cols = 
          c("#91D1C2FF","#8491B4FF","#F39B7FFF","#3C5488FF","#00A087FF","#4DBBD5FF","#AC9C88","#C82B1E","#CD7692","#4A885C","#E7AC53","#7A634B","#8E39C5"))

# Plot Figure 3e
data<-sn_sum[,-which(is.na(sn_sum$chr7_cnv))]
data$chr7_cnv[which(data$chr7_cnv>1.2)]<-1.2
FeaturePlot(data, features = c("chr7_cnv"),order = T) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(0.8,1.2))
data<-sn_sum[,-which(is.na(sn_sum$chr10_cnv))]
data$chr10_cnv[which(data$chr10_cnv<0.8)]<-0.8
data$chr10_cnv[which(data$chr10_cnv>1.2)]<-1.2
FeaturePlot(data, features = c("chr10_cnv")) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(0.8,1.2))

# Plot Figure 3f
data<-as.data.frame(proportions(table(sn_sum$celltype,sn_sum$orig.ident2),margin = 2))
data$condition<-substr(data$Var2,1,1)
data$pair<-gsub("I","",data$Var2)
data$pair<-gsub("R","",data$pair)
data$grp_3type<-"C2.NEU-"
data$grp_3type[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C2.NEU+")])]<-"C2.NEU+"
data$grp_3type[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C1.NEU+")])]<-"C1.NEU+"
data$grp_3type<-factor(data$grp_3type,levels = c("C1.NEU+","C2.NEU+","C2.NEU-"))
data<-data[which(data$Var1%in%c("Myeloid","Oligodendrocyte","Neuron","Lymphocyte","Endothelial","Astrocyte")),]
data<-data[order(data$pair),]
data<-data[order(data$condition),]
data0<-data[which(data$Var1=="Neuron"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "grp_3type")+
  stat_compare_means(paired = T, label = "p.format",method = "t.test")+
  ylab("Proportion of Neurons")+theme_classic()
data0<-data[which(data$Var1=="Oligodendrocyte"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "grp_3type")+
  stat_compare_means(paired = T, label = "p.format",method = "t.test")+
  ylab("Proportion of Oligodendrocytes")+theme_classic()
data0<-data[which(data$Var1=="Myeloid"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "grp_3type")+
  stat_compare_means(paired = T, label = "p.format",method = "t.test")+
  ylab("Proportion of Myeloid cells")+theme_classic()

# Plot Figure 3g
cell_meta<-sn_sum@meta.data
data0<-cell_meta[which(cell_meta$condition=="Primary"),]
data0<-data0[which(data0$celltype=="Tumor"),]
data<-as.data.frame(proportions(table(data0$grp_3type,data0$celltype_subtype),1))
data$condition<-"Primary"
data0<-cell_meta[which(cell_meta$condition=="Recurrent"),]
data0<-data0[which(data0$celltype=="Tumor"),]
data0<-as.data.frame(proportions(table(data0$grp_3type,data0$celltype_subtype),1))
data0$condition<-"Recurrent"
data<-rbind(data,data0)
data$condition2<-ifelse(data$condition=="Primary",0,1)
tumor_sub = c(Tumor.Cycle="#C82B1E",Tumor.OPC="#7A634B",Tumor.AC="#AC9C88",Tumor.NPC="#E7AC53",
              Tumor.MES="#4A885C",Tumor.Stem="#8E39C5",Tumor.Fiber="#CD7692")
data$Var2<-factor(data$Var2,levels = c("Tumor.NPC","Tumor.Cycle","Tumor.OPC","Tumor.AC","Tumor.MES","Tumor.Stem","Tumor.Fiber"))
data$Var1<-factor(data$Var1,levels = c("C1.NEU+","C2.NEU+","C2.NEU-"))
ggplot(data, aes(x=condition2, y=Freq, fill=Var2)) + 
  geom_area(alpha=1 , size=0.5, colour="black")+scale_fill_manual(values = tumor_sub)+
  facet_wrap(~Var1)+theme_classic2()

# Plot Figure 3h
sn_tumor<-subset(sn_sum,subset=celltype=="Tumor")
de_tumor_npc<-FindMarkers(sn_tumor,group.by = "celltype_subtype",ident.1 = "Tumor.NPC")
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tibble)
de_tumor_npc$gene<-row.names(de_tumor_npc)
gsea_genes<-de_tumor_npc %>%
  arrange(desc(p_val_adj), desc(avg_log2FC)) %>%
  dplyr::select(gene,avg_log2FC)
ranks <- deframe(gsea_genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nPermSimple = 10000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>% arrange(desc(NES))
fgseaResTidy$pathway<-gsub("GOBP_","",fgseaResTidy$pathway)
# fgseaResTidy<-readRDS("data/Fig3_GSEA_res.rds")
data<-fgseaResTidy[which(fgseaResTidy$pathway%in%c("NEUROTRANSMITTER_SECRETION","CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC","REGULATION_OF_MEMBRANE_POTENTIAL",
                                                   "GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY","DENDRITE_DEVELOPMENT","EXOCYTOSIS","ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME",
                                                   "MRNA_SPLICE_SITE_SELECTION","CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION","EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",              
                                                   "VASCULATURE_DEVELOPMENT","BLOOD_VESSEL_MORPHOGENESIS","INNATE_IMMUNE_RESPONSE","INFLAMMATORY_RESPONSE",                                      
                                                   "POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION","POSITIVE_REGULATION_OF_CELL_CELL_ADHESION",                  
                                                   "CELL_SUBSTRATE_ADHESION","REGULATION_OF_MITOTIC_CELL_CYCLE")),]
ggplot(data, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) + scale_fill_manual(values = c("#B22222","#008080")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GOBP pathways") +
  theme_minimal()

# Plot Figure 3i
load("data/Fig3_sn_celltype_splicing_res.RData")
ggscatter(pca.res,x="PC1",y="PC2",color = "group",shape = "type",label = "label2",palette = c("#FD8D0E","#106BFF","black"))

# Plot Figure 3j
data<-sn_sum@meta.data
data<-data[which(data$celltype%in%c("Neuron","Tumor")),]
data$celltype2<-as.character(data$celltype)
data$celltype2[which(data$celltype_subtype=="Tumor.NPC")]<-"NPC"
data$celltype2<-factor(data$celltype2,levels = c("Neuron","NPC","Tumor"))
data<-data[-which(data$chr7_cnv>1.5),]
data<-data[-which(data$chr10_cnv>1.5),]
data<-data[-which(is.na(data$chr10_cnv)),]
data$ID<-row.names(data)
library(reshape2)
data_long<-melt(data, id.vars = c("ID","celltype2"), 
                measure.vars = c("chr7_cnv","chr10_cnv"),
                variable.name = c('type'),
                value.name = 'chr_cnv')
ggboxplot(data_long,x="celltype2",y="chr_cnv",fill = "type",palette = "jco",outlier.shape = NA)+coord_flip()

# Plot supplementary Figure 6a
data<-as.data.frame(proportions(table(sn_sum$celltype,sn_sum$orig.ident2),margin = 2))
data$condition<-substr(data$Var2,1,1)
data$pair<-gsub("I","",data$Var2)
data$pair<-gsub("R","",data$pair)
data$group<-"C2.NEU-"
data$group[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C1.NEU+")])]<-"C1.NEU+"
data$group[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C2.NEU+")])]<-"C2.NEU+"
data$group<-factor(data$group,levels = c("C1.NEU+","C2.NEU+","C2.NEU-"))
data<-data[which(data$Var1%in%c("Myeloid","Oligodendrocyte","Neuron","Lymphocyte","Endothelial","Astrocyte")),]
data<-data[order(data$pair),]
data<-data[order(data$condition),]
data0<-data[which(data$Var1=="Neuron"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(paired = T, label = "p.format",method = "t.test")+
  ylab("Proportion of Neurons")+theme_classic()

# Plot supplementary Figure 6b
data0<-sn_sum@meta.data
data0<-data0[which(data0$celltype=="Tumor"),]
data<-as.data.frame(proportions(table(data0$celltype_subtype,data0$orig.ident2),margin = 2))
data$condition<-substr(data$Var2,1,1)
data$pair<-gsub("I","",data$Var2)
data$pair<-gsub("R","",data$pair)
data$group<-"C2.NEU-"
data$group[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C1.NEU+")])]<-"C1.NEU+"
data$group[which(data$pair%in%sn_meta$Pair.[which(sn_meta$grp_3type=="C2.NEU+")])]<-"C2.NEU+"
data$group<-factor(data$group,levels = c("C1.NEU+","C2.NEU+","C2.NEU-"))
data<-data[order(data$pair),]
data<-data[order(data$condition),]
data0<-data[which(data$Var1=="Tumor.NPC"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(paired = T, label = "p.format")+
  ylab("Tumor abundance of Tumor.NPC")+theme_classic()