######## for testing: 

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
 

