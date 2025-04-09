library(openxlsx)
library(dplyr)
TTTDinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240821_TTTD_QE_480_unique_patient_info_no_NA.xlsx")
sinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240107_TTTD_sample_info.xlsx")

machine_info=read.xlsx('D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)
machine_info$Histopathology_type = gsub('PTMC', 'PTC', machine_info$Histopathology_type)
#colnames(machine_info)
for(i in 1:nrow(machine_info)){
  tid = machine_info[i, 'Thyroid_ID']
  #pid = machine_info$
  if(is.na(tid)){next}
  sid = machine_info[i, 'sample']
  temp = subset(sinfo, Thyroid_ID == tid & raw_data == sid) #
  if (nrow(temp)==0){
    temp = subset(sinfo, Thyroid_ID == tid & raw_name == sid) #
  }
  if(nrow(temp)==1){
    machine_info[i, 'patient_ID'] = temp$patient_ID
  }
  if(nrow(temp)>1){print(i)}
}
#'tpd.TPD_add_2021_mzML.', 'tpd.TPD_SG.', 'rspa.RSPA2_FFPEtest.'

delinfo = subset(machine_info, Tissue_type=='thyroid' & Histopathology_type!='uncertain')
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
delinfo = subset(delinfo, Histopathology_type!='MEC' )
# FNA_in_vivo = subset(delinfo, !is.na(TrainTest) & SampleType=='FNA_in_vivo')
# write.xlsx(FNA_in_vivo, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_FNAinvivo_718samples.xlsx', rowNames = T)

delinfo = subset(delinfo, !(TrainTest=="Discovery" & SampleType=='FNA_in_vivo') )
delinfo = subset(delinfo, !( TrainTest=="Test2" & SampleType=='FNA_in_vivo') )
delinfo = subset(delinfo, !( TrainTest=="Test1" & SampleType=='FNA_in_vivo') )
dim(delinfo)
write.xlsx(delinfo, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_9624samples.xlsx', rowNames = T)

QE_deldat = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss.6359pros.6149samples_matrix.xlsx', rowNames = T)
#QE_deldat[1:3,1:3]
qeinfo = delinfo[colnames(QE_deldat), ]
qeinfoGender = subset(qeinfo, !is.na(qeinfo$Gender))
table(qeinfoGender$Classificaiton_type)
table(qeinfoGender$Histopathology_type)


dat480_deldat_min08 = read.xlsx('D:/chh/2023workProject/20240821TTTD/thyroid_480_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
dat480_deldat_min08[1:5, 1:3]

setdiff(rownames(dat480_deldat_min08), rownames(delinfo))
info480 = delinfo[rownames(dat480_deldat_min08),]
info480gender = subset(info480, !is.na(info480$Gender))

QETest3 = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/QETest3_log2min0.8.xlsx", rowNames = T)
QETest3[1:5, 1:3]
dat480Test3 = read.xlsx("D:/chh/2023workProject/20240821TTTD/dat480/dat480Test3_log2min0.8.xlsx", rowNames = T)

setdiff(rownames(QETest3), rownames(delinfo))
QETest3info = delinfo[rownames(QETest3), ]
dat480Test3info = delinfo[rownames(dat480Test3), ]

infoGender = rbind(qeinfoGender, info480gender)
infoGender = rbind(infoGender, QETest3info)
infoGender = rbind(infoGender, dat480Test3info)
infoGender =subset(infoGender,!is.na(TrainTest))
table(infoGender$Classificaiton_type)
infoGender = subset(infoGender, Classificaiton_type !='na')
write.xlsx(infoGender, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_reps.xlsx', rowNames = T)

length(unique(infoGender$patient_ID))
length(unique(infoGender$Thyroid_ID))
infoGenderPatient = infoGender %>% distinct( patient_ID, .keep_all = TRUE)
infoGenderPatient= subset(infoGenderPatient, !is.na(patient_ID))
dim(infoGenderPatient)
write.xlsx(infoGenderPatient, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6018patients.xlsx', rowNames = T)



# infoGenderPatient1 = subset(infoGenderPatient,  !(TrainTest!="Discovery" &  infoGenderPatient$qe480 == 'QE')) # "Discovery" "Test1"
# table(infoGenderPatient1$Classificaiton_type)
# table(infoGenderPatient1$Histopathology_type)

infoGenderThyroid = infoGender %>% distinct( Thyroid_ID, .keep_all = TRUE)
infoGenderThyroid= subset(infoGenderThyroid, !is.na(Thyroid_ID))
dim(infoGenderThyroid)
write.xlsx(infoGenderThyroid, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6109Thyroids.xlsx', rowNames = T)
infoGenderThyroid = read.xlsx( 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6109Thyroids.xlsx', rowNames = T)

# FNA_in_vivo
qe = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/QE_8509pros.7035samples_matrix.xlsx', rowNames = T)
qe[1:4,1:5]
dat480 = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_480all2700_20230624.pg_matrix.tsv", sep = '\t')
rownames(dat480) = paste0(dat480$Protein.Group,'_', dat480$Genes)
dat480[1:4,1:5]
dat480 = dat480[6:ncol(dat480)]
colnames(dat480) = sapply(colnames(dat480) , function(x){
  x = gsub('X..172.16.13.136.','', x, fixed = T)
  x = gsub('tttd.raw_480.MTC.','', x, fixed = T)
  x = gsub('tttd.raw_480.others.','', x, fixed = T)
  x = gsub('tttd.raw_480.TPM_ZY.','', x, fixed = T)
  x = gsub('tttd.raw_480.TPR.','', x, fixed = T)
  x = gsub('tttd.raw_480.SYF.','', x, fixed = T)
  x = gsub('.raw','', x, fixed = T)
  # x = gsub('tttd.raw_480.SYF.','', x, fixed = T)
  x
})
na_row = colSums(is.na(dat480))/nrow(dat480)
dat480 = dat480[na_row!=1]
dat480[1:3,1:3]
fnainvivo = rownames(FNA_in_vivo[ FNA_in_vivo$TrainTest !="Test3", ])

dat480_fnainvivo = data.frame(t(dat480[intersect(fnainvivo, colnames(dat480))]))
qe_fnainvivo = data.frame(t(qe[intersect(fnainvivo, colnames(qe))]))

dat480_fnainvivo = log2(dat480_fnainvivo)
dat480_fnainvivo[is.na(dat480_fnainvivo)] = 0.8*min(dat480_fnainvivo, na.rm = T)
dat480_fnainvivo[1:40,1:5]
write.xlsx(dat480_fnainvivo, "D:/chh/2023workProject/20240821TTTD/dat480/dat480FNA_in_vivo_DisTest12_log2min0.8.xlsx", rowNames = T)

qe_fnainvivo = log2(qe_fnainvivo)
qe_fnainvivo[is.na(qe_fnainvivo)] = 0.8*min(qe_fnainvivo, na.rm = T)
qe_fnainvivo[1:4,1:5]

write.xlsx(qe_fnainvivo, "D:/chh/2023workProject/20240821TTTD/dat480/qeFNA_in_vivo_DisTest12_log2min0.8.xlsx", rowNames = T)

#匹配，应该是参考Histology_ID，或者 将Patient_ID列，
#把HE_summary 里面的列：HE_image_file，
#ID加到 我们的sampleinfo 表格右边
he = read.xlsx('D:/chh/2023workProject/20240821TTTD/datas/HE_summary_20240613sunyt.xlsx')
for(i in 1:nrow(infoGenderThyroid)){
  pid = infoGenderThyroid[i, 'patient_ID']
  tid = infoGenderThyroid[i, 'Thyroid_ID']
  # tid = 'TTTD_T6567';pid = 'TTTD_P6567'
  temp = subset(TTTDinfo, patient_ID == pid & Thyroid_ID == tid)
  if(nrow(temp)!=1){prtint(i); break}
  if( is.na(temp$PatientID)){next}
  infoGenderThyroid[i, 'PATIENT_ID']=temp$PatientID
  
  if (temp$SampleType == 'FFPE_slide'){
    for (j in c(paste0('-A', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), paste0('-B', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), 
                paste0('-C', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), paste0('-D', c('1-1','1-2','2-1','2-2', '3-1', '3-2')))){
      temp$PatientID= gsub(j,'',temp$PatientID)
    }
  }
  
  temphe = subset(he, Patient_ID==temp$PatientID | Patient_ID== tolower(temp$PatientID) | Patient_ID== toupper(temp$PatientID))# PatientID
  if(nrow(temphe)==0){temphe = subset(he, Histology_ID==temp$PatientID | Histology_ID== tolower(temp$PatientID) | Histology_ID== toupper(temp$PatientID) )}#Histology_ID #
  
  if(nrow(temphe)==0 & grepl('^000', temp$PatientID)){
    tpid = as.character(as.numeric(temp$PatientID))
    temphe = subset(he, Patient_ID== tpid)# PatientID
    if(nrow(temphe)==0){temphe = subset(he, Histology_ID == tpid )}#Histology_ID
  }
  
  if(nrow(temphe)==2){temphe = subset(temphe, Project == temp$DataSet)}
  # if(nrow(temphe)==0 & grepl('^000', temp$PatientID)){
  #   temphe = subset(he, Patient_ID== gsub('^000','',temp$PatientID))# PatientID
  #   if(nrow(temphe)==0){temphe = subset(he, Histology_ID == gsub('^000','',temp$PatientID) )}#Histology_ID
  # }
  if(nrow(temphe)!=1){next}
  # print(nrow(temphe))
  infoGenderThyroid[i, 'HE_ID'] = temphe$ID
  infoGenderThyroid[i, 'HE_image_file'] = temphe$HE_image_file
}
write.xlsx(infoGenderThyroid, 'D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6109Thyroids_20250102.1.xlsx', rowNames = T)

library(stringi)
TTTD_QE_480_20240803 = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240803_TTTD_QE_480_unique_patient_info.xlsx")
#colnames(TTTD_QE_480_20240803)
for(i in 1:nrow(TTTD_QE_480_20240803)){
  # tid = 'TTTD_T6567';pid = 'TTTD_P6567'
  temp = TTTD_QE_480_20240803[i,] #subset(TTTD_QE_480_20240803, PatientID==pid)
  pid = temp$PatientID
  if (is.na(pid) | pid=='delete')next
  if (!is.na(temp$SampleType)){
    if (temp$SampleType == 'FFPE_slide'){
    for (j in c(paste0('-A', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), paste0('-B', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), 
                paste0('-C', c('1-1','1-2','2-1','2-2', '3-1', '3-2')), paste0('-D', c('1-1','1-2','2-1','2-2', '3-1', '3-2')))){
      temp$PatientID= gsub(j,'',temp$PatientID)
    }
  }
  }
  temphe = subset(he, Patient_ID==temp$PatientID | Patient_ID== tolower(temp$PatientID) | Patient_ID== toupper(temp$PatientID))# PatientID
  if(nrow(temphe)==0){temphe = subset(he, Histology_ID==temp$PatientID | Histology_ID== tolower(temp$PatientID) | Histology_ID== toupper(temp$PatientID) )}#Histology_ID #
  
  if(nrow(temphe)==0 & grepl('^000', temp$PatientID)){
    tpid = as.character(as.numeric(temp$PatientID))
    temphe = subset(he, Patient_ID== tpid)# PatientID
    if(nrow(temphe)==0){temphe = subset(he, Histology_ID == tpid )}#Histology_ID
  }
  
  if(nrow(temphe)==0 & temp$DataSet == 'MFT'){
    tpid = substr(pid,1,stri_length(pid)-1)
    temphe = subset(he, Patient_ID== tpid)# PatientID
    if(nrow(temphe)==0){temphe = subset(he, Histology_ID == tpid )}#Histology_ID
  }
  # if(nrow(temphe)==0 & temp$DataSet == 'MFT' & grepl('^R', pid)){
  #   tpid = strsplit(pid ,'-')[[1]][1]
  #   temphe = subset(he, Patient_ID== tpid)# PatientID
  #   if(nrow(temphe)==0){temphe = subset(he, Histology_ID == tpid )}#Histology_ID
  # }
  if(nrow(temphe)>1){temphe = subset(temphe, Project == temp$DataSet)}
  if(nrow(temphe)!=1){next}
  TTTD_QE_480_20240803[i, 'HE_ID'] = temphe$ID
  TTTD_QE_480_20240803[i, 'HE_image_file'] = temphe$HE_image_file
}
TTTD_QE_480_20240803$HE_image_file = gsub('^na', '', TTTD_QE_480_20240803$HE_image_file)
write.xlsx(TTTD_QE_480_20240803, 'D:/chh/2023workProject/20240821TTTD/datas/20240803_TTTD_QE_480_unique_patient_info.1.xlsx', rowNames = F)


########## 20250115 #######
xgb_results = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QE_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_results$TrainTest = machine_info[rownames(xgb_results), 'TrainTest']
xgb_results$samples = rownames(xgb_results)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/480_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_result$samples = rownames(xgb_result)
xgb_results = rbind(xgb_results, xgb_result)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QEtest3_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_result$samples = rownames(xgb_result)
xgb_results = rbind(xgb_results, xgb_result)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/480test3_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_result$samples = rownames(xgb_result)
xgb_results = rbind(xgb_results, xgb_result)

xgb_results$class = machine_info[rownames(xgb_results), 'Classificaiton_type']
xgb_results$his = machine_info[rownames(xgb_results), 'Histopathology_type']
xgb_results$Thyroid_ID = machine_info[rownames(xgb_results), 'Thyroid_ID']
xgb_results$qe480 = machine_info[rownames(xgb_results), 'qe480']
xgb_results$patient_ID = machine_info[rownames(xgb_results), 'patient_ID']

unique(xgb_results$his)
unique(xgb_results$class)
xgb_results$his = gsub('PTMC', 'PTC', xgb_results$his)
xgb_results$Age_gap = xgb_results$Pred - xgb_results$Age
xgb_results = subset(xgb_results, !is.na( Pred) )
unique(xgb_results$Gender)
for(i in 1:nrow(xgb_results)){
  tid = xgb_results[i, 'Thyroid_ID']
  xgb_results[i, 'Bethesda_Score']  = TTTDinfo[TTTDinfo$Thyroid_ID == tid, 'Bethesda_Score']
  xgb_results[i, 'BRAF']  = TTTDinfo[TTTDinfo$Thyroid_ID == tid, 'BRAF']
}
xgb_results$sid = machine_info[rownames(xgb_results), 'sample']

# uploadfile = read.xlsx('D:/chh/2023workProject/20240821TTTD/upload_file.xlsx' )

for(i in 1:nrow(xgb_results)){
  # i=1
  sid = xgb_results[i, 'sid']
  y = gsub("\\D", "", sid)
  xgb_results[i, 'year']= as.numeric(stri_sub( y, 1,4 ))
}

xgb_results$unpublished = ifelse(xgb_results$year > 2019, 'unpublished', NA)

# for(dt in unique(xgb_results$TrainTest)){
#   for (h in unique(xgb_results$Hospital)){
#     temp = subset(xgb_results, TrainTest==dt & Hospital==h)
#     if(nrow(temp)==0)next
#     p = unique(temp$patient_ID)
#     t = unique(temp$Thyroid_ID)
#     temp = temp %>% distinct( Thyroid_ID, .keep_all = TRUE)
#     unp = subset(temp, !is.na(unpublished))
#     print(c(dt, h, length(p), length(t), nrow(unp)))
#   }
# }


# patients = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6018patients.xlsx", rowNames = T)
# 
# differ = subset(xgb_results, !(xgb_results$patient_ID %in% patients$patient_ID))
# 
# intersect(rownames(differ), rownames(QE_deldat))
# patients = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6018patients.xlsx", rowNames = T)
# Thyroid = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/sampleinfo_6109Thyroids.xlsx", rowNames = T)


dat = infoGender
dat$sid = machine_info[rownames(dat), 'sample']
for(i in 1:nrow(dat)){
  # i=1
  sid = dat[i, 'sid']
  y = gsub("\\D", "", sid)
  dat[i, 'year']= as.numeric(stri_sub( y, 1,4 ))
}
dat$unpublished = ifelse(dat$year > 2019, 'unpublished', NA)
for(dt in unique(dat$TrainTest)){
  for (h in unique(dat$Hospital)){
    temp = subset(dat, TrainTest==dt & Hospital==h)
    if(nrow(temp)==0)next
    p = unique(temp$patient_ID)
    t = unique(temp$Thyroid_ID)
    temp = temp %>% distinct( Thyroid_ID, .keep_all = TRUE)
    unp = subset(temp, !is.na(unpublished))
    print(c(dt, h, length(p), length(t), nrow(unp)))
  }
}
