# for testing: 

library(imputeLCMD)

#dia_data <- read.table('for_heatmap_Dia_sum_med_6056_146_log2.txt')

#names:
#New_TMT_144_7083_for_imputing.txt
#New_TMT_144_7083_merged_with_second_strategy_for_imputing.txt
# run NMF for 12129 and 144:
#final_144_12129_for_imputing.txt
#final_144_10663_second_str_for_imputing.txt
#New_TMT_144_9412_cut_off_72_merged_second_strategy_for_imputing.txt
  
#for_imputing_8498_genes_949_cell_lines_pancancer.txt
#For_imputing_New_DIA_9564_247_test_train_tumors.txt


#For_imputing_Final_File_TMT_215_2763_Corrected_HarmonizR.txt

#for_imputing_Final_File_TMT_136_tumors_merged_second_str_12096_Corrected_HarmonizR.txt
dia_data <- read.table('for_imputing_Final_File_TMT_136_tumors_merged_second_str_12096_Corrected_HarmonizR.txt',header=T,sep="\t")

#second strategy:
#dia_data <- read.table('for_heatmap_Dia_sum_med_146_6056.txt')

dia_data<- matrix(unlist(dia_data), ncol = 136)


exprsData.imputed = impute.MinProb(dia_data,0.01,1)
exprsData.imputed
write.table(exprsData.imputed, 'Imputed_Final_File_TMT_136_tumors_merged_second_str_12096_Corrected_HarmonizR.txt.txt')

################## imputing for 122 TNBC tumor sample (2023):
#name of file: "For_imputing_log2_New_DIA_122_tumors_9564_genes.txt"
# "For_imputing_DIA_122_9564_normalized_median_log2.txt"
# for TMT data :   "for_imputing_TMT_12129_122.txt"
# for RNA data: "for_imputing_RNA_Log2_19657_121_tumors.txt"

#For_imputing_log2_DIA_225_TNBC_9564.txt
#For_imputing_DIA_IBCT_66_10888.txt
library(imputeLCMD)
dia_data <- read.table('For_imputing_TMT_data_Protegenomics_paper_11063_71.txt',header=T,sep="\t")
dia_data<- matrix(unlist(dia_data), ncol = 71)
exprsData.imputed = impute.MinProb(dia_data,0.01,1)
exprsData.imputed
write.table(exprsData.imputed, 'Imputed_TMT_data_Protegenomics_paper_11063_71.txt')

# QRILC approch:
dia_data <- read.table('For_imputing_log2_DIA_225_TNBC_9564.txt',header=T,sep="\t")
dia_data<- matrix(unlist(dia_data), ncol = 225)
exprsData.imputed = impute.QRILC(dia_data, tune.sigma = 1)
exprsData.imputed=exprsData.imputed[[1]]
write.table(exprsData.imputed, 'Imputed_log2_DIA_225_TNBC_9564_QRILC.txt')


#########################################################
  #    calculating quantile for each sample for TMT   #
#quantile(data_q, 0.01, na.rm = TRUE)

#data_q<-read.table("spearman_correlation_coloumns_RNA_Normalized_diaPASEF.txt",header=F,sep="\n")
#data_q<-as.vector(data_q)

data_q <- read.table('for_quantile_calculating.txt',header=T,sep="\t")
data_q<- matrix(unlist(data_q), ncol = 122)


print(quantile(as.vector(data_q[,1]), 0.01, na.rm = TRUE))


quantile(data_q, 0.01, na.rm = TRUE)




################ enrichment for DIA data: 
library(GSVA)
library(readxl)
library(RColorBrewer)


############## top annotation for DIA  #####

#DIA:
Data_file <- read_excel("Data_file_for_annotation_New_TMT.xlsx")
annot_df_DIA <- data.frame(NMF=Data_file$NMF5,
                          #Cons=Data_file$pearson5_New_TMT,
                          
                          PAM50=Data_file$PAM50,
                          TNBCtype=Data_file$TNBCtype,
                          TNBC =Data_file$TNBC,
                          #ESR1=as.double(Data_file$ESR1),
                          #AR=as.double(Data_file$AR),
                          #FOXA1 =as.double(Data_file$FOXA1),
                          ERBB2=as.double(Data_file$ERBB2_DIA_N),
                          EGFR =as.double(Data_file$EGFR_DIA_N),
                          #MET =as.double(Data_file$MET),
                          CLDN3=as.double(Data_file$CLDN3_DIA_N),
                          CDH1 =as.double(Data_file$CDH1_DIA_N),
                          VIM =as.double(Data_file$VIM_DIA_N),
                          #SOX10=as.double(Data_file$SOX10),
                          KRT5=as.double(Data_file$KRT5_DIA_N),
                          KRT14=as.double(Data_file$KRT14_DIA_N),
                          KRT17=as.double(Data_file$KRT17_DIA_N)
                          
                          
                          
                          
                          
)
col_fun = colorRamp2(c(-1,0,1), c("#4575b4", "white", "#d73027"))

# for Dia_TMT:
col_fun_new = colorRamp2(c(0.2,0.4,0.8), c("#4575b4", "white", "#d73027"))

# for Ki67:
col_fun_new_1 = colorRamp2(c(10,40,100), c("#4575b4", "white", "#d73027"))

# ASCAT_TUM_FRAC:
col_fun_new_2= colorRamp2(c(0.1,0.4,0.99), c("#4575b4", "white", "#d73027"))


col = list(PAM50 = c("Basal-like" = "#e31a1c", "Luminal A" = "#1f78b4", 
                     "Luminal B" = "#a6cee3", "HER2" = "#fb9a99",
                     "Normal-like" = "#33a02c", "NA"="lightgrey"),
           
           #NMF = c("1" = "blue", "2" = "red","3"="orange","4"="black","5" = "yellow","6"="purple","7"="lightblue"),
           NMF = c("1" = "#882255", "2" = "#DDCC77","3"="#999933","4"="#CC6677","5"="#44AA99"),
           #,"6"="black"),
           # "5"="yellow","6"="black")
           TNBC = c("Yes" = "black", "No" = "white"),
           Cons=c("1" = "green", "2" = "blue","3"="red","4"="yellow","5"="orange"),
           #,"6"="purple"),
           #,"5"="orange","6"="purple")
           ESR1 = col_fun,ERBB2 = col_fun, KRT5 = col_fun, AR = col_fun, CLDN3 = col_fun, VIM = col_fun, CDH1 = col_fun, EGFR = col_fun, MET = col_fun, FOXA1 = col_fun,KRT14 = col_fun,KRT17 = col_fun,SOX10=col_fun,
           Dia_TMT=col_fun_new,
           ER=c("0"="white", "<5"="#fc9272" ,"5-10"="#de2d26","NA"="lightgrey"),
           HER2=c("0+"="white","1+"="#fc9272","2+"="#de2d26","NA"="lightgrey"),
           Ki67=col_fun_new_ki,
           TNBCtype=c("M"="#E69F00","LAR"="#56B4E9","BL1"="#009E73","BL2"="#0072B2","IM"="#D55E00","MSL"="#CC79A7","NA"="lightgrey"),
           ASCAT_F=col_fun_new_2,
           HRDetect=c("high"="#d73027","intermed"="white","low"="#4575b4","NA"="lightgrey"),
           Grade=c( "2"="#f4a582","3"="#b2182b","NA"="lightgrey"),
           Age = c("<50" = "#fddbc7", "50-75" = "#ef8a62","75-95" = "#b2182b","NA"="lightgrey"),
           Chemo= c("Chemo" = "black", "Non-Chemo" = "white","NA"="lightgrey"),
           mpg = circlize::colorRamp2(c(17, 25), 
                                      c("lightblue", "purple")) )





ha_DIA <- HeatmapAnnotation(df = annot_df_DIA,
                        which="col", 
                        col = col,
                        simple_anno_size = unit(3, "mm"),
                        annotation_name_gp = gpar(fontsize = 8),
                        annotation_name_side = "left", 
                        na_col = "lightgrey", 
                        gap=unit(0.2, "mm"))
#draw(ha)
################################################




#dataa_DIA<-read.table("for_enrichment_DIA_data_all_subgroups_1 to 146.txt",header=TRUE,sep="\t")
dataa_DIA<-read.table("for_enrichment_DIA_data_Imputed_medianP_medianT.txt",header=TRUE,sep="\t")

title=tempdir()

dataa_DIA<- matrix(unlist(dataa_DIA), ncol = 146,
               dimnames=list(paste0("g", 1:6056), paste0("s", 1:146)))

#Data_file <- read_excel("For_enrichmenr_50_pathways.xlsx")
#Data_file_en <- read_excel("For_enrichmenr_50_pathways_New_TMT_7083_genes.xlsx")


Data_file_en_DIA <- read_excel("For_enrichmenr_50_pathways_DIA_data_6056_genes.xlsx")


gs_new_order_DIA <- list(Data_file_en_DIA$HALLMARK_APICAL_JUNCTION,
                         Data_file_en_DIA$HALLMARK_APICAL_SURFACE,
                         Data_file_en_DIA$HALLMARK_PEROXISOME,
                         Data_file_en_DIA$HALLMARK_ADIPOGENESIS,
                         Data_file_en_DIA$HALLMARK_ANGIOGENESIS,
                         Data_file_en_DIA$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
                         Data_file_en_DIA$HALLMARK_MYOGENESIS,
                         Data_file_en_DIA$HALLMARK_SPERMATOGENESIS,
                         Data_file_en_DIA$HALLMARK_PANCREAS_BETA_CELLS,
                         Data_file_en_DIA$HALLMARK_DNA_REPAIR,
                         Data_file_en_DIA$HALLMARK_UV_RESPONSE_DN,
                         Data_file_en_DIA$HALLMARK_UV_RESPONSE_UP,
                         Data_file_en_DIA$HALLMARK_ALLOGRAFT_REJECTION,
                         Data_file_en_DIA$HALLMARK_COAGULATION,
                         Data_file_en_DIA$HALLMARK_COMPLEMENT,
                         Data_file_en_DIA$HALLMARK_INTERFERON_ALPHA_RESPONSE,
                         Data_file_en_DIA$HALLMARK_INTERFERON_GAMMA_RESPONSE,
                         Data_file_en_DIA$HALLMARK_IL6_JAK_STAT3_SIGNALING,
                         Data_file_en_DIA$HALLMARK_INFLAMMATORY_RESPONSE,
                         Data_file_en_DIA$HALLMARK_BILE_ACID_METABOLISM,
                         Data_file_en_DIA$HALLMARK_CHOLESTEROL_HOMEOSTASIS,
                         Data_file_en_DIA$HALLMARK_FATTY_ACID_METABOLISM,
                         Data_file_en_DIA$HALLMARK_GLYCOLYSIS,
                         Data_file_en_DIA$HALLMARK_HEME_METABOLISM,
                         Data_file_en_DIA$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                         Data_file_en_DIA$HALLMARK_XENOBIOTIC_METABOLISM,
                         Data_file_en_DIA$HALLMARK_APOPTOSIS,
                         Data_file_en_DIA$HALLMARK_HYPOXIA,
                         Data_file_en_DIA$HALLMARK_PROTEIN_SECRETION,
                         Data_file_en_DIA$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                         Data_file_en_DIA$HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY,
                         Data_file_en_DIA$HALLMARK_E2F_TARGETS,
                         Data_file_en_DIA$HALLMARK_G2M_CHECKPOINT,
                         Data_file_en_DIA$HALLMARK_MYC_TARGETS_V1,
                         Data_file_en_DIA$HALLMARK_MYC_TARGETS_V2,
                         Data_file_en_DIA$HALLMARK_P53_PATHWAY,
                         Data_file_en_DIA$HALLMARK_MITOTIC_SPINDLE,
                         Data_file_en_DIA$HALLMARK_ANDROGEN_RESPONSE,
                         Data_file_en_DIA$HALLMARK_ESTROGEN_RESPONSE_EARLY,
                         Data_file_en_DIA$HALLMARK_ESTROGEN_RESPONSE_LATE,
                         Data_file_en_DIA$HALLMARK_IL2_STAT5_SIGNALING,
                         Data_file_en_DIA$HALLMARK_KRAS_SIGNALING_UP,
                         Data_file_en_DIA$HALLMARK_KRAS_SIGNALING_DN,
                         Data_file_en_DIA$HALLMARK_MTORC1_SIGNALING,
                         Data_file_en_DIA$HALLMARK_NOTCH_SIGNALING,
                         Data_file_en_DIA$HALLMARK_PI3K_AKT_MTOR_SIGNALING,
                         Data_file_en_DIA$HALLMARK_HEDGEHOG_SIGNALING,
                         Data_file_en_DIA$HALLMARK_TGF_BETA_SIGNALING,
                         Data_file_en_DIA$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
                         Data_file_en_DIA$HALLMARK_WNT_BETA_CATENIN_SIGNALING)

gsva_DIA.es <- gsva(dataa_DIA, gs_new_order_DIA, verbose=FALSE)
dim(gsva_DIA.es)

annot_df_2_DIA <- data.frame(gsva_DIA.es)

annot_df_2_2_DIA<- t(annot_df_2_DIA)
dim(annot_df_2_2_DIA)



###############

colnames(annot_df_2_2_DIA) <- c("APICAL_JUNCTION","APICAL_SURFACE","PEROXISOME","ADIPOGENESIS","ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION","MYOGENESIS","SPERMATOGENESIS","PANCREAS_BETA_CELLS","DNA_REPAIR","UV_RESPONSE_DN","UV_RESPONSE_UP","ALLOGRAFT_REJECTION","COAGULATION","COMPLEMENT","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE","BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","FATTY_ACID_METABOLISM","GLYCOLYSIS","HEME_METABOLISM","OXIDATIVE_PHOSPHORYLATION","XENOBIOTIC_METABOLISM","APOPTOSIS","HYPOXIA","PROTEIN_SECRETION","UNFOLDED_PROTEIN_RESPONSE","REACTIVE_OXIGEN_SPECIES_PATHWAY","E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","P53_PATHWAY","MITOTIC_SPINDLE","ANDROGEN_RESPONSE","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","IL2_STAT5_SIGNALING","KRAS_SIGNALING_UP","KRAS_SIGNALING_DN","MTORC1_SIGNALING","NOTCH_SIGNALING","PI3K_AKT_MTOR_SIGNALING","HEDGEHOG_SIGNALING","TGF_BETA_SIGNALING","TNFA_SIGNALING_VIA_NFKB","WNT_BETA_CATENIN_SIGNALING")



mycols_new <- colorRamp2(breaks = c(-0.2, 0,0.2), 
                         colors = c("#4575b4", "white", "#d73027")) # quant heatmap - blue, white, red
#col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))
col_2222 = list(APICAL_JUNCTION=mycols_new,APICAL_SURFACE=mycols_new,PEROXISOME=mycols_new,ADIPOGENESIS=mycols_new,ANGIOGENESIS=mycols_new,EPITHELIAL_MESENCHYMAL_TRANSITION=mycols_new,MYOGENESIS=mycols_new,SPERMATOGENESIS=mycols_new,PANCREAS_BETA_CELLS=mycols_new,DNA_REPAIR=mycols_new,UV_RESPONSE_DN=mycols_new,UV_RESPONSE_UP=mycols_new,ALLOGRAFT_REJECTION=mycols_new,COAGULATION=mycols_new,COMPLEMENT=mycols_new,INTERFERON_ALPHA_RESPONSE=mycols_new,INTERFERON_GAMMA_RESPONSE=mycols_new,IL6_JAK_STAT3_SIGNALING=mycols_new,INFLAMMATORY_RESPONSE=mycols_new,BILE_ACID_METABOLISM=mycols_new,CHOLESTEROL_HOMEOSTASIS=mycols_new,FATTY_ACID_METABOLISM=mycols_new,GLYCOLYSIS=mycols_new,HEME_METABOLISM=mycols_new,OXIDATIVE_PHOSPHORYLATION=mycols_new,XENOBIOTIC_METABOLISM=mycols_new,APOPTOSIS=mycols_new,HYPOXIA=mycols_new,PROTEIN_SECRETION=mycols_new,UNFOLDED_PROTEIN_RESPONSE=mycols_new,REACTIVE_OXIGEN_SPECIES_PATHWAY=mycols_new,E2F_TARGETS=mycols_new,G2M_CHECKPOINT=mycols_new,MYC_TARGETS_V1=mycols_new,MYC_TARGETS_V2=mycols_new,P53_PATHWAY=mycols_new,MITOTIC_SPINDLE=mycols_new,ANDROGEN_RESPONSE=mycols_new,ESTROGEN_RESPONSE_EARLY=mycols_new,ESTROGEN_RESPONSE_LATE=mycols_new,IL2_STAT5_SIGNALING=mycols_new,KRAS_SIGNALING_UP=mycols_new,KRAS_SIGNALING_DN=mycols_new,MTORC1_SIGNALING=mycols_new,NOTCH_SIGNALING=mycols_new,PI3K_AKT_MTOR_SIGNALING=mycols_new,HEDGEHOG_SIGNALING=mycols_new,TGF_BETA_SIGNALING=mycols_new,TNFA_SIGNALING_VIA_NFKB=mycols_new,WNT_BETA_CATENIN_SIGNALING=mycols_new,
                mpg = circlize::colorRamp2(c(17, 25), 
                                           c("lightblue", "purple")) )


split_2=c("Apical_junction","Apical_surface","Peroxisome","Adipogenesis","Angiogenesis","Epithelial_mesenchymal_transition","Myogenesis","Spermatogenesis","Pancreas_BETA_cells","DNA_repair","UV_response_DN","UV_response_UP","Allograft_rejection","Coagulation","Complement","Interferon_alpha_response","Interferon_gamma_response","IL6_JAK_STAT3_signaling","Inflammatory_response","Bile_Acid_metabolism","Cholesterol_homeostasis","Fatty_Acid_metabolism","Glycolysis","Heme_metabolism","Oxidative_phosphorylation","Xenobiotic_metabolism","Apoptosis","Hypoxia","Protein_secretion","Unfolded_protein_response","Reactive_oxigen_species_pathway","E2F_targets","G2M_Checkpoint","MYC_targets_V1","MYC_targets_V2","P53_pathway","Mitotic_spindle","Androgen_response","Estrogen_response_early","Estrogen_response_late","IL2_STAT5_signaling","KRAS_signaling_UP","KRAS_signaling_DN","MTORC1_signaling","NOTCH_signaling","PI3K_AKT_MTOR_signaling","Hedgehog_signaling","TGF_BETA_signaling","TNFA_signaling_via_NFKB","WNT_BETA_Catenin_signaling")



ha_2_DIA <- HeatmapAnnotation(df = annot_df_2_2_DIA, 
                          which="col", 
                          col = col_2222,
                          simple_anno_size = unit(1.5, "mm"),
                          #annotation_name_gp = gpar(fontsize = 5),
                          annotation_name_gp = gpar(col = c(rep("black", 3), rep("black", 6), rep("black", 3), rep("black", 7), rep("black", 7), rep("black", 5), rep("black", 6), rep("black", 13)),fontsize = 4.2),
                          annotation_name_side = "right", 
                          na_col = "lightgrey",
                          show_legend=FALSE,
                          annotation_label = split_2,
                          # gap = unit(1, "points"),
                          # space = unit(3, "mm"),
                          #column_split  = (rep(c("A", "B"), 25)),
                          #row_split  = split,
                          height = unit(9, "cm"),
                          #gap=unit(0.2, "mm")
                          gap=unit(c(0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), "mm"))

################# heat map for dia data



#ddd_new_1_215_dia<-read.table("ddd_new_1_146_Dia_data_medianT_medianP.txt",header=TRUE,sep="\t")
ddd_new_1_215_dia<-read.table("for_enrichment_DIA_data_Imputed_medianP_medianT.txt",header=TRUE,sep="\t")

title=tempdir()

ddd_new_1_215_dia<- matrix(unlist(ddd_new_1_215_dia), ncol = 146)
#colnames(ddd_new_1_215) = paste0("RNr", seq_len(146))
colnames(ddd_new_1_215_dia) = c("RNr04","RNr05","RNr07","RNr08","RNr09","RNr10","RNr11","RNr20","RNr21","RNr23","RNr24","RNr25","RNr27","RNr28","RNr29","RNr30","RNr31","RNr32","RNr33","RNr34","RNr35","RNr36","RNr37","RNr38","RNr39","RNr40","RNr41","RNr42","RNr44","RNr45","RNr46","RNr47","RNr48","RNr49","RNr50","RNr55","RNr56","RNr60","RNr62","RNr64","RNr65","RNr66","RNr67","RNr68","RNr69","RNr72","RNr73","RNr74","RNr77","RNr79","RNr80","RNr81","RNr82","RNr83","RNr84","RNr85","RNr86","RNr87","RNr88","RNr90","RNr91","RNr94","RNr95","RNr96","RNr97","RNr99","RNr101","RNr102","RNr103","RNr104","RNr106","RNr107","RNr108","RNr109","RNr112","RNr114","RNr116","RNr120","RNr121","RNr122","RNr123","RNr124","RNr128","RNr129","RNr131","RNr132","RNr136","RNr140","RNr144","RNr145","RNr146","RNr147","RNr153","RNr154","RNr157","RNr159","RNr162","RNr163","RNr164","RNr212","RNr213","RNr214","RNr215","RNr216","RNr217","RNr218","RNr219","RNr220","RNr221","RNr222","RNr223","RNr224","RNr225","RNr226","RNr228","RNr229","RNr230","RNr231","RNr232","RNr233","RNr234","RNr235","RNr236","RNr237","RNr238","RNr239","RNr240","RNr241","RNr242","RNr243","RNr244","RNr245","RNr246","RNr247","RNr248","RNr249","RNr250","RNr251","RNr252","RNr262","RNr263","RNr264","RNr265","RNr266","RNr269","RNr274"
)
hr <- hclust(as.dist(1-cor(t(ddd_new_1_215_dia), method="pearson", use = "pairwise.complete.obs")), method="ward.D2")# cluster rows for
mycols_dia <- colorRamp2(breaks = c(-1, 0, 1), 
                         colors = c("#4575b4", "white", "#d73027")) # quant heatmap - blue, white, red


Ht77<- Heatmap(ddd_new_1_215_dia, 
               name = "Protein quant", #title of legend
               column_title = "Breast tumors", 
               row_title = "Proteins",
               cluster_columns =as.dendrogram(consensushc(estim.r)),
               column_names_gp = gpar(fontsize = 3.1),
               show_row_names = FALSE,
               show_column_names = T, 
               cluster_rows = hr,
               col = mycols_dia,
               column_dend_reorder = F, # makes quant data come in better order
               row_dend_reorder = T, # makes quant data come in better order
               top_annotation = ha_DIA,
               
               bottom_annotation = ha_2_DIA,
               #right_annotation = ha_R,
               width = unit(10.7, "cm"),
               height = unit(4.7, "cm"),
               # use_raster = FALSE,
               # raster_resize_mat = TRUE,
               heatmap_legend_param = list(legend_direction = "horizontal")) 
Ht77

lgd_gsva = Legend(title = "gsva score", col_fun = mycols_new, 
)

draw(lgd_gsva, x = unit(0.62, "npc"), y = unit(0.36, "npc"))

library(circlize)

decorate_annotation("APICAL_SURFACE", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(5.3, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("Cell component   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))


######################################################

decorate_annotation("EPITHELIAL_MESENCHYMAL_TRANSITION",{
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(9.35, "mm") , gp = gpar(fill = TRUE),
            hjust =c(1,1) ,vjust=c(0.58,0))
  grid.text(paste("Development   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

decorate_annotation("UV_RESPONSE_DN", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(5.3, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("DNA-damage    ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

decorate_annotation("INTERFERON_ALPHA_RESPONSE", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(11, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("Immune   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

decorate_annotation("GLYCOLYSIS", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(11, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("Metabolic   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

decorate_annotation("PROTEIN_SECRETION", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(7.8, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("Pathway   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))


decorate_annotation("MYC_TARGETS_V1",{
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(9.35, "mm") , gp = gpar(fill = TRUE),
            hjust =c(1,1) ,vjust=c(0.58,0))
  grid.text(paste("proliferation   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

decorate_annotation("MTORC1_SIGNALING", {
  grid.rect(x = 0, width = unit(1.1, "mm"),height= unit(21, "mm") , gp = gpar(fill = TRUE),
            just = "right")
  grid.text(paste("Signaling   ", collapse = "\n"),x = unit(1, "mm") , just = "right",gp = gpar(fontsize = 9))
}, slice = 1,
envir = new.env(parent = parent.frame()))

