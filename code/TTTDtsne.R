library(openxlsx)
library(umap)
library(Rtsne)
library(ggplot2)
library(ggsci)

qe = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_QE_all7038_20230623.pg_matrix.tsv", sep = '\t')
dat480 = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_480all2700_20230624.pg_matrix.tsv", sep = '\t')

TTTDinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240821_TTTD_QE_480_unique_patient_info_no_NA.xlsx")
sinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240107_TTTD_sample_info.xlsx")

rownames(qe) = paste0(qe$Protein.Group,'_', qe$Genes)
qe[1:4,1:5]
qe = qe[6:ncol(qe)]
colnames(qe) = sapply(colnames(qe) , function(x){
  x = gsub('X..172.16.13.136.','', x, fixed = T)
  x
})
na_row = colSums(is.na(qe))/nrow(qe)
qe = qe[na_row!=1]
# for(i in TTTDinfo$raw_data){
#   i = paste0(i,'.')
#   index = grep(i, colnames(qe), fixed = T)
#   if(length(index)>1){
#     print(i)
#   }
#   # colnames(qe)[index] = i
# }

rownames(dat480) = paste0(dat480$Protein.Group,'_', dat480$Genes)
dat480 = dat480[6:ncol(dat480)]
colnames(dat480)[1:10]
colnames(dat480) = sapply(colnames(dat480) , function(x){
  y = strsplit(x, '.', fixed = T)[[1]][10]
  y
})
na_row = colSums(is.na(dat480))/nrow(dat480)
dat480 = dat480[na_row!=1]
dat480[1:10,1:10]
qe[1:10,1:2]

# temp = data.frame(intensity = na.omit(unlist(qe)))
# temp$intensity = as.numeric(temp$intensity)
# # temp= data.frame(intensity = sample(1:10000, 500))
# ggplot(temp , aes(x = intensity))+
#   # geom_histogram(bins = 20)
#   geom_density(color = 'black' )

# merge
qe480 = data.frame(row.names = unique(c(rownames(qe), rownames(dat480))))
qe480[rownames(qe) ,colnames(qe)] = qe
qe480[rownames(dat480) ,colnames(dat480)] = dat480
qe480 = qe480[!grepl(";", rownames(qe480)),]
qe480 = qe480[!grepl(":",rownames(qe480)),]
qe480 = qe480[!grepl("iRT",rownames(qe480)),]
na_row = colSums(is.na(qe480))/nrow(qe480)
na_row[na_row==1]

na_col = colSums(!is.na(qe480))
length(na_col[na_col<1000])
# na_col = rowSums(is.na(qe480))/ncol(qe480)
# na_col[na_col==1]
machine_info = data.frame(row.names = colnames(qe480), qe480 = c(rep('QE', ncol(qe)), rep('480', ncol(dat480))) )
machine_info$machine = NA
for(i in rownames(machine_info)){
  # A,B，CAA, CAB, CAC, O
  m='A'
  if(grepl('.A2', i, fixed = T)){
    m='A'
  }else if(grepl('.B2', i, fixed = T)){
    m='B'
  }else if(grepl('O2', i, fixed = T)){
    m='O'
  }else if(grepl('CAA2', i, fixed = T)){
    m='CAA'
  }else if(grepl('CAB2', i, fixed = T)){
    m='CAB'
  }else if(grepl('CAC2', i, fixed = T)){
    m='CAC'
  }else if(grepl('C2', i, fixed = T)){
    m='C'
  }else{
    print(i)
    next
  }
  machine_info[i, 'machine'] = m
}
unique(machine_info[ is.na(machine_info$machine), 'qe480'])
machine_info[machine_info$machine=='A', ]
dim(qe480)
qe480 = qe480[na_row!=1]
machine_info = machine_info[colnames(qe480),]
na_col = colSums(!is.na(qe480))
machine_info$protein_num = na_col
machine_info$sample = rownames(machine_info)
del = c('tpd.TPD_add_2021_mzML.', 'tpd.TPD_FNA_395.', 'tpd.TPD_retrospective_raw.', 'tpd.TPD_SG.',
  'rspa.RSPA1_FFPEmzML.', 'rspa.RSPA2_FFPEtest.', 'rspa.RSPA3_fna_mzML.', 'tttd.raw_qe.raw.', '.mzML', '.raw')
for (i in del){
  machine_info$sample = gsub(i, '', machine_info$sample, fixed = T)
}
unique(TTTDinfo$DataSet)
unique(TTTDinfo$SampleType)
unique(TTTDinfo$Histopathology_type)

# for(j in unique(TTTDinfo$DataSet)){
#   machine_info[grepl(j, machine_info$sample), c("DataSet" )] = j
# }

for(i in unique(machine_info$sample)){
  tid = sinfo[sinfo$raw_name == i, 'Thyroid_ID' ]
  if(length(tid)!=0){
    machine_info[machine_info$sample == i, "Thyroid_ID"] = tid
  }
  temp = TTTDinfo[TTTDinfo$Thyroid_ID == tid, c("DataSet", "SampleType", "Histopathology_type", "Gender","Age","Tissue_type","Hospital","Classificaiton_type")]
  if(nrow(temp) == 0 ){
    temp = TTTDinfo[TTTDinfo$raw_data == i, c("DataSet", "SampleType", "Histopathology_type", "Gender","Age","Tissue_type","Hospital","Classificaiton_type")]
  }
  if(nrow(temp) == 0 ){ next }
  
  machine_info[machine_info$sample == i, c( "DataSet" )] = unique(temp[, 1])
  machine_info[machine_info$sample == i, c( "SampleType" )] = unique(temp[,2])
  machine_info[machine_info$sample == i, c( "Histopathology_type" )] = unique(temp[,3])
  machine_info[machine_info$sample == i, c( "Gender" )] = unique(temp[,4])
  machine_info[machine_info$sample == i, c( "Age" )] = unique(temp[,5])
  machine_info[machine_info$sample == i, c( "Tissue_type" )] = unique(temp[, 6])
  machine_info[machine_info$sample == i, c( "Hospital" )] = unique(temp[, 7])
  machine_info[machine_info$sample == i, "Classificaiton_type"] = unique(temp[, 8])
}
for(i in unique(machine_info$sample)){
  tid = sinfo[sinfo$raw_name == i, 'Thyroid_ID' ]
  temp = TTTDinfo[TTTDinfo$Thyroid_ID == tid, c("DataSet", "Classificaiton_type")]
  if(nrow(temp) == 0 ){
    temp = TTTDinfo[TTTDinfo$raw_data == i, c("DataSet", "Classificaiton_type")]
  }
  if(nrow(temp) == 0 ){ next }

  machine_info[machine_info$sample == i, "Classificaiton_type"] = unique(temp[,2])
}

machine_info$Histopathology_type = gsub("PTC ", "PTC",machine_info$Histopathology_type)
unique(machine_info$DataSet)
unique(machine_info$SampleType)
unique(machine_info$Histopathology_type)
# 0-20,20-40,40-60,60-80，80-100
machine_info$Age_grade = NA
temp = subset(machine_info, Age<20)
machine_info[rownames(temp), 'Age_grade'] = '<20'
temp = subset(machine_info, Age<40 & Age>=20)
machine_info[rownames(temp), 'Age_grade'] = '20-40'
temp = subset(machine_info, Age<60& Age>=40)
machine_info[rownames(temp), 'Age_grade'] = '40-60'
temp = subset(machine_info, Age<80& Age>=60)
machine_info[rownames(temp), 'Age_grade'] = '60-80'
temp = subset(machine_info, Age>=80)
machine_info[rownames(temp), 'Age_grade'] = '≥80'

unique(machine_info$Gender)
machine_info$Gender = gsub(' F','F', machine_info$Gender)
machine_info$Gender = gsub('f','F', machine_info$Gender) 
machine_info$Gender = gsub('m','M', machine_info$Gender) 
machine_info$Gender = gsub('F ','F', machine_info$Gender) 
machine_info$Gender = gsub('na',NA, machine_info$Gender)


library(stringr)
pattern <- "\\d{8}"
for(i in 1:nrow(machine_info)){
  text = machine_info[i,'sample']
  match <- str_extract(text, pattern)
  if (!is.na(match)) {
    # print(paste("Extracted date:", match))
    # date_parsed <- as.Date(match, format = "%Y%m%d")
    # formatted_date <- format(date_parsed, "%Y-%m-%d")
    machine_info[i, 'date'] = match
  }
}
machine_info$date = as.numeric(machine_info$date)
machine_info[grepl(21180611, machine_info$date) , 'date'] = 20180611

# write.xlsx(machine_info, 'D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)
machine_info=read.xlsx('D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)
#### barplot
unique(machine_info$DataSet)

ggplot(machine_info , aes(DataSet, protein_num, color = SampleType )) +
  geom_boxplot()+ 
  #geom_jitter( size=1)+
  scale_color_lancet()+
  theme_classic()+
  labs(x = 'DataSet', y = 'protein_nums')+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5  ),
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black")
  )
# delete NA

# filter1 = machine_info[machine_info$DataSet != 'RSPA' & machine_info$SampleType != 'FNA_in_vivo', ]
filter1 = machine_info[!is.na(machine_info$DataSet), ]# 10429 sample

qe480_min0.8 = qe480[rownames(filter1)]
qe480_min0.8[is.na(qe480_min0.8)] = min(qe480_min0.8, na.rm = T) * 0.8

qe_min0.8 = qe[rownames(machine_info[machine_info$qe480 == 'QE' & !is.na(machine_info$DataSet),])]
qe_min0.8[is.na(qe_min0.8)] = min(qe_min0.8, na.rm = T) * 0.8

dat480_min0.8 = dat480[rownames(machine_info[machine_info$qe480 == '480' & !is.na(machine_info$DataSet),])]
dat480_min0.8[is.na(dat480_min0.8)] = min(dat480_min0.8, na.rm = T) * 0.8

# tsneinput = qe480_min0.8
# tsneinput = qe_min0.8
tsneinput = dat480_min0.8
set.seed(2024)
dfQE.tsne <- Rtsne(t(tsneinput), dims = 2, perplexity = 10,
                   partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
dim(dfQE.tsne$Y)
dim(machine_info )
# "DataSet", "SampleType", "Histopathology_type" "Gender","Age_grade","Tissue_type"

df1 = data.frame(dfQE.tsne$Y, machine_info[colnames(tsneinput), 'date'] )
names(df1) = c("tSNE1", "tSNE2", "Batch")
df1$Batch = as.numeric(df1$Batch)
dim(df1)
head(df1) 
# write.xlsx(df1, "tsne_filterInfoNA_qe.xlsx", rowNames = T)

ggplot(df1, aes(tSNE1,tSNE2, colour= Batch))+
  geom_point(size=2)+ # aes(color= Batch), 
  scale_color_gradient(low='cyan', high = 'red')+
  # scale_color_d3(palette = c( "category20" ))+# scale_color_jama()
  theme_bw()+
  #scale_colour_manual(values =  colors12)+ #+scale_colour_manual(values = c("green", "red"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13,color = "black"),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.25),
        axis.line.y = element_line(color="black", size = 0.25),
        plot.title  = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        axis.title  = element_text(size=13 ),
        axis.title.x = element_text( hjust = 0.5 ),
        axis.title.y  = element_text( hjust = 0.5 ),
        panel.background = element_blank(),
        #legend.position = c(0.9,0.1), # legend的位置信息
        #legend.background = element_rect(fill = "grey90", size = 1, colour = "white")
  ) + theme(plot.title = element_text(hjust = 0.5))

#
qe480[1:10, 1:3]
pool = qe480[grepl('ool',colnames(qe480))]
colnames(machine_info)
unique(machine_info$SampleType)
unique(machine_info$Tissue_type)
unique(machine_info$Histopathology_type)

delinfo = subset(machine_info, Tissue_type=='thyroid' & Histopathology_type!='uncertain' & SampleType!='FNA_in_vivo')
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
# write.xlsx(delinfo, 'D:/chh/2023workProject/20240821TTTD/QC/thyroid_9440sampleinfos.xlsx', rowNames = T)
unique(delinfo$SampleType)
unique(delinfo$Tissue_type)
unique(delinfo$Histopathology_type)
# deldat = qe480[rownames(delinfo)]# 9440
deldat = qe[rownames(delinfo[delinfo$qe480== 'QE',])]# 6448
# deldat = dat480[rownames(delinfo[delinfo$qe480 == '480',])]# 2992
dim(deldat)
na_col = colSums(!is.na(deldat))
deldat = deldat[na_col>=2000 ]
na_row = rowSums(is.na(deldat))/ncol(deldat)
length(na_row)
na_row[na_row==1 ]
deldat = deldat[na_row<=0.9, ]
dim(deldat)# QE:6359pro*6149sample,480:7499pro*2976
# write.xlsx(deldat, 'D:/chh/2023workProject/20240821TTTD/QC/thyroid_480.6642pros.2992samples_matrix.xlsx', rowNames = T)

Q = quantile(deldat, na.rm = T )# 
Q
Q3 = Q[4]
Q1 = Q[2]
IQR = Q3-Q1
IQR
up = Q3+1.5*IQR
down = Q1-1.5*IQR
deldat[1:20,1:2]
print(c(up, down))
deldat[deldat>up] = up
deldat[deldat<down] = down

qe_del_IQR = deldat
dat480_del_IQR = deldat
qe480IQR = data.frame(row.names = unique(c(rownames(qe_del_IQR), rownames(dat480_del_IQR))))
qe480IQR[rownames(qe_del_IQR), colnames(qe_del_IQR)] = qe_del_IQR
qe480IQR[rownames(dat480_del_IQR), colnames(dat480_del_IQR)] = dat480_del_IQR
dim(qe480IQR)# 7837 9125
# write.xlsx(qe480IQR, 'D:/chh/2023workProject/20240821TTTD/QC/thyroid_qe480IQR.7837pros.9125samples_matrix.xlsx', rowNames = T)
qe480IQR[1:3,1:3]
qe480IQR = read.xlsx('D:/chh/2023workProject/20240821TTTD/QC/thyroid_qe480IQR.7837pros.9125samples_matrix.xlsx', rowNames = T)

qe480IQR_min0.8 = qe480IQR
min(qe480IQR_min0.8, na.rm = T)
qe480IQR_min0.8[is.na(qe480IQR_min0.8)] = min(qe480IQR_min0.8, na.rm = T) * 0.8


library(sva)
dim(qe480IQR_min0.8)
qe480IQR_min0.8[1:20, 1:2]
colnames(machine_info)
qe480IQR_min0.8_combat = ComBat( qe480IQR_min0.8, batch = machine_info[colnames(qe480IQR_min0.8), 'qe480']) # , ref.batch = 'QE'
qe480IQR_min0.8_combat[1:10, 1:3]
min(qe480IQR_min0.8_combat, na.rm = T)
qe480IQR_min0.8_combat[qe480IQR_min0.8_combat<0] = 782.2736
tsneinput = qe480IQR_min0.8_combat
set.seed(2024)
dfQE.tsne <- Rtsne(t(tsneinput), dims = 2, perplexity = 10,
                   partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
dim(dfQE.tsne$Y)
dim(machine_info )
# "DataSet", "SampleType", "Histopathology_type" "Gender","Age_grade","Tissue_type"
df1 = data.frame(dfQE.tsne$Y, machine_info[colnames(tsneinput), 'qe480'] ) 
names(df1) = c("tSNE1", "tSNE2", "Batch")
dim(df1)
head(df1) 
# write.xlsx(df1, "tsne_filterInfoNA_qe.xlsx", rowNames = T)
ggplot(df1, aes(tSNE1,tSNE2, colour= Batch))+
  geom_point(aes(color= Batch), size=2)+
  scale_color_d3(palette = c( "category20" ))+# scale_color_jama()
  theme_bw()+
  #scale_colour_manual(values =  colors12)+ #+scale_colour_manual(values = c("green", "red"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13,color = "black"),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.25),
        axis.line.y = element_line(color="black", size = 0.25),
        plot.title  = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        axis.title  = element_text(size=13 ),
        axis.title.x = element_text( hjust = 0.5 ),
        axis.title.y  = element_text( hjust = 0.5 ),
        panel.background = element_blank(),
        #legend.position = c(0.9,0.1), # legend的位置信息
        #legend.background = element_rect(fill = "grey90", size = 1, colour = "white")
  ) + theme(plot.title = element_text(hjust = 0.5))


BiocManager::install("proBatch")
library(dplyr)
library(tibble)
library(ggplot2)
library(proBatch)
data('example_proteome', 'example_sample_annotation', 'example_peptide_annotation', package = 'proBatch')
feature_id_col = 'Protein_group_label'
measure_col = 'Intensity'
sample_id_col = 'FullRunName'
essential_columns = c(feature_id_col, measure_col, sample_id_col)
batch_col = 'qe480'


color_list <- sample_annotation_to_colors(machine_info,
                                          factor_columns = c('qe480', 'machine',"DataSet", "SampleType", "Histopathology_type",
                                                             'Tissue_type', "Gender"  ),
                                          numeric_columns = c("date" , "Age" ))
generated_sample_annotation <- date_to_sample_order(example_sample_annotation,
                                                    time_column = c('RunDate','RunTime'),
                                                    new_time_column = 'generated_DateTime',
                                                    dateTimeFormat = c('%b_%d', '%H:%M:%S'),
                                                    new_order_col = 'generated_order',
                                                    instrument_col = NULL)
library(knitr)
kable(generated_sample_annotation[1:5,] %>%
        select(c('RunDate', 'RunTime', 'order', 'generated_DateTime', 'generated_order')))
example_matrix <- long_to_matrix(example_proteome,
                                 feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')
log_transformed_matrix <- log_transform_dm(qe480IQR_min0.8, log_base = 2, offset = 1) # pro* file——>log2
log_transformed_matrix[1:3,1:3]
quantile_normalized_matrix = normalize_data_dm(as.matrix(log_transformed_matrix), normalize_func = 'quantile')
quantile_normalized_matrix[1:3,1:3]
machine_info$FullRunName = rownames(machine_info)

machine_info1 = machine_info[colnames(qe480IQR_min0.8),]
rownames(machine_info1) = 1:nrow(machine_info1)
machine_info1$order = 1;nrow(machine_info1)
# plot_PCA(quantile_normalized_matrix, machine_info, color_by='qe480', sample_id_col = "FullRunName",
#          plot_title= 'qe480',color_scheme=color_list[['qe480']])
plot_PVCA(quantile_normalized_matrix, machine_info1,
          technical_factors = c('qe480', 'machine'), biological_factors = c("Gender", "DataSet"))

quantile_normalized_long<-matrix_to_long(quantile_normalized_matrix)
colnames(quantile_normalized_long)
# colnames(quantile_normalized_long)[1] = "proteins"
head(quantile_normalized_long)
max(quantile_normalized_long$Intensity)
batch_col = 'qe480'
loess_fit_30=adjust_batch_trend_df(quantile_normalized_long, machine_info1, batch_col = batch_col, 
                                   #feature_id_col = "proteins", #order_col = "order",
                                   span=0.3)
head(loess_fit_30)
unique(loess_fit_30$fit)
plot_with_fitting_curve(feature_name = 'A0A0B4J2D5_GATD3B', #feature_id_col = "proteins",
                        fit_df=loess_fit_30, fit_value_col='fit',
                        df_long=quantile_normalized_long,
                        sample_annotation=machine_info1, batch_col = batch_col, 
                        color_by_batch=TRUE, color_scheme=color_list[[batch_col]],
                        plot_title = 'Span=30%')
example_sample_annotation$order


batch_corrected_df <- correct_batch_effects_df(df_long = quantile_normalized_long,
                                               sample_annotation = machine_info1,
                                               batch_col = batch_col,# feature_id_col = "proteins",
                                               discrete_func = 'MedianCentering',
                                               continuous_func = 'loess_regression',
                                               abs_threshold = 5, pct_threshold = 0.20)
head(batch_corrected_df)
unique(batch_corrected_df$fit)

temp = long_to_matrix(batch_corrected_df,
               feature_id_col = 'peptide_group_label',
               measure_col = 'preBatchCorr_Intensity',
               sample_id_col = 'FullRunName')
temp[1:3,1:3]
dfQE.tsne <- Rtsne(t(temp), dims = 2, perplexity = 10,
                   partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
dim(dfQE.tsne$Y)
dim(machine_info )
# "DataSet", "SampleType", "Histopathology_type" "Gender","Age_grade","Tissue_type"
df1 = data.frame(dfQE.tsne$Y, machine_info[colnames(temp), 'qe480'] ) 
names(df1) = c("tSNE1", "tSNE2", "Batch")
dim(df1)
head(df1) 
# write.xlsx(df1, "tsne_filterInfoNA_qe.xlsx", rowNames = T)
ggplot(df1, aes(tSNE1,tSNE2, colour= Batch))+
  geom_point(aes(color= Batch), size=2)+
  scale_color_d3(palette = c( "category20" ))+# scale_color_jama()
  theme_bw()+
  #scale_colour_manual(values =  colors12)+ #+scale_colour_manual(values = c("green", "red"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13,color = "black"),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.25),
        axis.line.y = element_line(color="black", size = 0.25),
        plot.title  = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        axis.title  = element_text(size=13 ),
        axis.title.x = element_text( hjust = 0.5 ),
        axis.title.y  = element_text( hjust = 0.5 ),
        panel.background = element_blank(),
        #legend.position = c(0.9,0.1), # legend的位置信息
        #legend.background = element_rect(fill = "grey90", size = 1, colour = "white")
  ) + theme(plot.title = element_text(hjust = 0.5))
