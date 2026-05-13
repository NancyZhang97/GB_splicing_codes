####################################################
# Include codes to re-generate plots in Figure 2 and Supplementary Figure 2
# Necessary data can be found in ".../data/" directory
####################################################

# Plot Supplementary Figure 2a
IDHwt_tmz_clinical<-readRDS("data/Fig2_patient_clinical_info.rds")
IDHwt_tmz_clinical$grp_3type<-factor(IDHwt_tmz_clinical$grp_3type,levels = c("C1.NEU+","C2.NEU+","C2.NEU-"))
data<-t(IDHwt_tmz_clinical[order(IDHwt_tmz_clinical$grp_3type,IDHwt_tmz_clinical$Patient,decreasing = F),])
data<-data[2:7,]
library(ComplexHeatmap)
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  CGGA = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#B22222", col = NA)),
  NG2016 = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#0066FFFF", col = NA)),
  NM2019 = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FF7F50", col = NA)),
  SMCnew = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6A5ACD", col = NA)),
  TCGA = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FFD700", col = NA)),
  Asian = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#006400", col = NA)),
  `non-Asian` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#DB7093", col = NA)),
  M = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6495ED", col = NA)),
  `F` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CD5C5C", col = NA)),
  Yes = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#BA55D3", col = NA)),
  No = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#00CED1", col = NA)),
  II = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E6E6FA", col = NA)),
  III = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#8B008B", col = NA)),
  IV = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#9370DB", col = NA))
)
oncoPrint(data, alter_fun = alter_fun,row_order = 1:nrow(data),column_order = 1:ncol(data),
          col = c(CGGA = "#B22222", NG2016 = "#0066FFFF",NM2019="#FF7F50",
                  SMCnew ="#6A5ACD", TCGA ="#FFD700", Asian="#006400", `non-Asian`="#DB7093",
                  M="#6495ED", `F`="#CD5C5C", Yes="#BA55D3", 
                  No="#00CED1", II="#E6E6FA", III="#8B008B", IV="#9370DB"))

# Plot Supplementary Figure 2b
library(ggpubr)
comparisons<-list(c("C1.NEU+","C2.NEU+"),c("C2.NEU+","C2.NEU-"),c("C1.NEU+","C2.NEU-"))
ggboxplot(IDHwt_tmz_clinical,x="grp_3type",y="age",fill = "grp_3type",palette = c("#EE9238","#3069F6","#6C7FAE"))+
  stat_compare_means(label = "p.format",comparisons = comparisons, method = "t.test")

# Plot Supplementary Figure 2c
mut_matrix<-readRDS("data/mutation_info_matrix.rds")
library(ComplexHeatmap)
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  # red rectangles
  Initial_only = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#B22222", col = NA)),
  # blue rectangles
  Recurrent_only = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#0066FFFF", col = NA)),
  # dots
  Both_pos = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FF7F50", col = NA)),
  Both_neg = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E6E6FA", col = NA)),
  `C1.NEU+` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#EE9238", col = NA)),
  `C2.NEU+` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#3069F6", col = NA)),
  `C2.NEU-` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6C7FAE", col = NA))
)
oncoPrint(mut_matrix, alter_fun = alter_fun,row_order = 1:nrow(mut_matrix),column_order = 1:ncol(mut_matrix),
          col = c(Initial_only = "#B22222", Recurrent_only = "#0066FFFF",Both_pos="#FF7F50",
                  Both_neg ="#E6E6FA", `C1.NEU+` ="#EE9238", `C2.NEU+`="#3069F6",`C2.NEU-`="#6C7FAE"),show_column_names = TRUE)

# Plot Figure 2a
mutation_info<-read.table("data/mutation_info_IDHwt_tmz.txt", header = T, sep = "\t")
mutation_info$group<-factor(mutation_info$group, levels = c("C1.NEU+", "C2.NEU+", "C2.NEU-"))
data<-mutation_info[which(mutation_info$condition=="Initial"),]
data0<-as.data.frame(proportions(table(data$PTENdel,data$group),margin = 2))
colnames(data0)<-c("Event","group","Freq")
library(ggplot2)
ggplot(data=data0, aes(x=" ", y=Freq, group=group, fill=Event)) + 
  geom_bar(width = 1, stat = "identity") + 
  #scale_color_manual(values = c(`TRUE` = "#AD5BCD", `FALSE` = "#3DCBCF")) +
  scale_fill_manual(values = c(`0` = "#D3D3D3", `1` = "#00CED1")) +
  coord_polar("y", start=0) +  
  facet_grid(.~ group) +theme_void() 

data0<-as.data.frame(proportions(table(data$PIK3CA,data$group),margin = 2))
colnames(data0)<-c("Event","group","Freq")
ggplot(data=data0, aes(x=" ", y=Freq, group=group, fill=Event)) + 
  geom_bar(width = 1, stat = "identity") + 
  #scale_color_manual(values = c(`TRUE` = "#AD5BCD", `FALSE` = "#3DCBCF")) +
  scale_fill_manual(values = c(`0` = "#D3D3D3", `1` = "#00CED1")) +
  coord_polar("y", start=0) +  
  facet_grid(.~ group) +theme_void() 

data<-mutation_info[which(mutation_info$condition=="Recurrent"),]
data0<-as.data.frame(proportions(table(data$CDKN2Adel,data$group),margin = 2))
colnames(data0)<-c("Event","group","Freq")
ggplot(data=data0, aes(x=" ", y=Freq, group=group, fill=Event)) + 
  geom_bar(width = 1, stat = "identity") + 
  #scale_color_manual(values = c(`TRUE` = "#AD5BCD", `FALSE` = "#3DCBCF")) +
  scale_fill_manual(values = c(`0` = "#D3D3D3", `1` = "#BA55D3")) +
  coord_polar("y", start=0) +  
  facet_grid(.~ group) +theme_void() 

# Plot Figure 2b
library(forcats)
gsea_res<-readRDS("data/GSEA_RvsP_grp3type.rds")
gsea_res %>%
  arrange(pathway) %>%
  ggplot(aes(NES, pathway)) +
  geom_line(aes(group = pathway)) +
  geom_point(aes(color = group, size = log_pval))+
  scale_color_manual(values = c(`C1.NEU+` = "#EE9238", `C2.NEU+` = "#3069F6", `C2.NEU-` = "#6C7FAE"))+theme_pubclean() + theme(legend.position = "right")

# Plot Figure 2c
subtype_score<-readRDS("data/Fig1_dse_psi_info.rds")
subtype_score<-subtype_score[,c(1:5,163:166)]
data<-IDHwt_tmz_clinical[,c("Patient","grp_3type")]
subtype_score<-merge(subtype_score, data, by = "Patient")
subtype_score<-subtype_score[order(subtype_score$Patient),]
subtype_score<-subtype_score[order(subtype_score$condition),]
ggpaired(subtype_score,x="condition",y="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of MES markers")
ggpaired(subtype_score,x="condition",y="VERHAAK_GLIOBLASTOMA_NEURAL",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of NEU markers")

# Plot Figure 2d, Supplementary Figure 2d
tpm_sum<-readRDS("data/TPM_102_GB.rds")
data<-subtype_score
data$exp<-tpm_sum["MGMT",data$SampleID]
ggpaired(data,x="condition",y="exp",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Normalized TPM of MGMT")
data$exp<-tpm_sum["MLH1",data$SampleID]
ggpaired(data,x="condition",y="exp",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Normalized TPM of MLH1")
data$exp<-tpm_sum["ANXA3",data$SampleID]
ggpaired(data,x="condition",y="exp",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T, label.y = 4.1)+ theme_classic()+ ylab("Normalized TPM of ANXA3")
data$exp<-tpm_sum["EGR4",data$SampleID]
ggpaired(data,x="condition",y="exp",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T, label.y = 6.1)+ theme_classic()+ ylab("Normalized TPM of EGR4")
data$exp<-tpm_sum["MSH6",data$SampleID]
ggpaired(data,x="condition",y="exp",facet.by = "grp_3type",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T, label.y = 5.9)+ theme_classic()+ ylab("Normalized TPM of MSH6")

# Plot Supplementary Figure 2e
data<-subtype_score
data<-cbind(data, t(tpm_sum[c("MGMT","ANXA3","EGR4","MLH1","MSH6"),data$SampleID]))
library(reshape2)
data_long<-melt(data, id.vars = c("SampleID", "grp_3type", "condition"), 
                measure.vars = c("MGMT","ANXA3","EGR4","MLH1","MSH6"),
                variable.name = c('TMZ_gene'),
                value.name = 'TPM')
data_long<-data_long[which(data_long$condition=="Recurrent"),]
comparisons<-list(c("C1.NEU+","C2.NEU+"),c("C2.NEU+","C2.NEU-"),c("C1.NEU+","C2.NEU-"))
ggboxplot(data_long,x="grp_3type",y="TPM",fill = "grp_3type",palette = c("#EE9238","#3069F6","#6C7FAE"), facet.by = "TMZ_gene")+
  stat_compare_means(comparisons = comparisons,label = "p.format",method = "t.test", label.y = c(7,7.2,8.1))+theme_classic2()+ylim(c(0,9))
