library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(reshape2)
# https://lehallib.github.io/DEswan/articles/DEswan.html
stdscale = function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}
pval_coefbind = function(pvalues, coefficients1, type1, dirs1){
  pvalues$factor = gsub('covariates$','',pvalues$factor, fixed = T)
  colnames(coefficients1) = gsub('covariates.','', colnames(coefficients1), fixed = T)
  colnames(coefficients1) = paste0(colnames(coefficients1), '_coef')
  pvalues_p = acast(pvalues[c("variable", "factor", "pvalue")], variable~factor)
  pvalues_p = data.frame(pvalues_p)
  pvalues_padjust = apply(pvalues_p, 2, function(x){p.adjust(x, method = 'BH')})
  pvalues_padjust = data.frame(pvalues_padjust)
  colnames(pvalues_p) = paste0(colnames(pvalues_padjust), '_pvalues')
  colnames(pvalues_padjust) = paste0(colnames(pvalues_padjust), '_padjust')
  coefficients1[rownames(pvalues_p) ,colnames(pvalues_p)] = pvalues_p
  coefficients1[rownames(pvalues_padjust) ,colnames(pvalues_padjust)] = pvalues_padjust
  coefficients1 = coefficients1[sort(colnames(coefficients1))]
  coefficients1$Uniprot = sapply(rownames(coefficients1), function(x){
    strsplit(x,'_')[[1]][1]
  })
  coefficients1$Gene = sapply(rownames(coefficients1), function(x){
    strsplit(x,'_')[[1]][2]
  })
  coefficients1 = cbind(coefficients1[grepl('qt',colnames(coefficients1))], coefficients1[!grepl('qt',colnames(coefficients1))])
  write.xlsx(coefficients1, paste0(dirs1, type1, '_glm_ageContinue.xlsx'), rowNames=T)
}


qe = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_QE_all7038_20230623.pg_matrix.tsv", sep = '\t')
# dat480 = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_480all2700_20230624.pg_matrix.tsv", sep = '\t')
# rownames(dat480) = paste0(dat480$Protein.Group,'_', dat480$Genes)
# dat480[1:4,1:5]
# dat480 = dat480[6:ncol(dat480)]
# colnames(dat480) = sapply(colnames(dat480) , function(x){
#   x = gsub('X..172.16.13.136.','', x, fixed = T)
#   x = gsub('tttd.raw_480.MTC.','', x, fixed = T)
#   x = gsub('tttd.raw_480.others.','', x, fixed = T)
#   x = gsub('tttd.raw_480.TPM_ZY.','', x, fixed = T)
#   x = gsub('tttd.raw_480.TPR.','', x, fixed = T)
#   x = gsub('tttd.raw_480.SYF.','', x, fixed = T)
#   x = gsub('.raw','', x, fixed = T)
#   # x = gsub('tttd.raw_480.SYF.','', x, fixed = T)
#   x
# })
# na_row = colSums(is.na(dat480))/nrow(dat480)
# dat480 = dat480[na_row!=1]

dat480_deldat_min08 = read.xlsx('D:/chh/2023workProject/20240821TTTD/thyroid_480_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
dat480_deldat_min08[1:5, 1:3]
dat480_deldat_min08
setdiff(rownames(dat480_deldat_min08), rownames(delinfo))
info480 = delinfo[rownames(dat480_deldat_min08),]
info480gender = subset(info480, !is.na(info480$Gender))

rownames(qe) = paste0(qe$Protein.Group,'_', qe$Genes)
qe[1:4,1:5]
qe = qe[6:ncol(qe)]
colnames(qe) = sapply(colnames(qe) , function(x){
  x = gsub('X..172.16.13.136.','', x, fixed = T)
  x
})
na_row = colSums(is.na(qe))/nrow(qe)
qe = qe[na_row!=1]
qe[1:5,1:2]
dim(qe)
write.xlsx(qe, 'D:/chh/2023workProject/20240821TTTD/QE/QE_8509pros.7035samples_matrix.xlsx', rowNames = T)
qe = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/QE_8509pros.7035samples_matrix.xlsx', rowNames = T)

# write.xlsx(deldat, 'D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss.6359pros.6149samples_matrix.xlsx', rowNames = T)
TTTDinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240821_TTTD_QE_480_unique_patient_info_no_NA.xlsx")
sinfo = read.xlsx("D:/chh/2023workProject/20240821TTTD/datas/20240107_TTTD_sample_info.xlsx")

machine_info=read.xlsx('D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)
machine_info$Histopathology_type = gsub('PTMC', 'PTC', machine_info$Histopathology_type)

colnames(machine_info)
machine_info[1:2,1:2]
unique(machine_info$Histopathology_type)
unique(machine_info$Classificaiton_type)
unique(machine_info$Tissue_type)
unique(machine_info$SampleType)
unique(machine_info$Hospital)

machine_info[machine_info$Hospital %in% c('SG', 'Y', 'N', 'Z', 'DS', 'I', 'SJ', 'Q', 'M' ), 'TrainTest'] = 'Discovery'
machine_info[machine_info$Hospital %in% c('S', 'ZE', 'F', 'B', 'DLschool', 'A' ), 'TrainTest'] = 'Test1'
machine_info[machine_info$Hospital %in% c('P', 'DL', 'J') & machine_info$qe480 == '480', 'TrainTest'] = 'Test2'
machine_info[machine_info$Hospital %in% c('D', 'W', 'C', 'H' ) & machine_info$SampleType == 'FNA_in_vivo', 'TrainTest'] = 'Test3'
unique(machine_info$TrainTest)
unique(machine_info[machine_info$TrainTest == 'Discovery', 'qe480'])

for(i in 1:nrow(machine_info)){
  tid = machine_info[i, 'Thyroid_ID']
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

write.xlsx(machine_info, 'D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)

delinfo = subset(machine_info, Tissue_type=='thyroid' & Histopathology_type!='uncertain' & SampleType!='FNA_in_vivo')
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
delinfo = subset(delinfo, Histopathology_type!='MEC' )

delinfo = subset(delinfo, qe480!='480' )
unique(delinfo$TrainTest)
QE_deldat = qe[,rownames(delinfo)]
na_col = colSums(!is.na(QE_deldat))
dim(QE_deldat) # 8509 6448
QE_deldat = QE_deldat[na_col>=2000 ] # 299
dim(QE_deldat)# 8509 6149
na_row = rowSums(is.na(QE_deldat))/ncol(QE_deldat)
length(na_row)
na_row[na_row==1 ]
QE_deldat = QE_deldat[na_row<=0.9, ]
dim(QE_deldat)# QE:6359pro*6149sample, 480:7499pro*2976
# write.xlsx(QE_deldat, 'D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss.6359pros.6149samples_matrix.xlsx', rowNames = T)

QE_deldat = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss.6359pros.6149samples_matrix.xlsx', rowNames = T)


QE_deldat = log2(QE_deldat)
QE_deldat[1:10,1:2]
delinfo = delinfo[colnames(QE_deldat), ]# important
unique(delinfo$qe480)

Q = quantile(QE_deldat, na.rm = T)
Q3 = Q[4]
Q1 = Q[2]
IQR = Q3-Q1
IQR
up = Q3+2*IQR
down = Q1-2*IQR
up # 26.33689
down # 13.47285 
QE_deldat_IQR = QE_deldat
QE_deldat_IQR[QE_deldat_IQR>up] = up
QE_deldat_IQR[QE_deldat_IQR<down] = down
temp = NA
plot(density(na.omit(unlist(QE_deldat_IQR))), main="2IQR density default")
sum(is.na(QE_deldat_IQR))/(ncol(QE_deldat_IQR) * nrow(QE_deldat_IQR))

min(QE_deldat_IQR, na.rm = T) * 0.8 # 10.77828

QE_deldat_IQR_min0.8 = QE_deldat_IQR
QE_deldat_IQR_min0.8[is.na(QE_deldat_IQR_min0.8)] = min(QE_deldat_IQR_min0.8, na.rm = T) * 0.8
QE_deldat_IQR_min0.8[1:20,1:2]
QE_deldat_IQR_min0.8 = data.frame(t(QE_deldat_IQR_min0.8))
dim(QE_deldat_IQR_min0.8)
plot(density(na.omit(unlist(QE_deldat_IQR_min0.8 ))), main="density default")
write.xlsx(data.frame(miss60 = colnames(QE_deldat_IQR_min0.8)), 'D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.6miss_pro.xlsx')

QE_deldat_IQR_min0.8$Age = as.numeric(machine_info[rownames(QE_deldat_IQR_min0.8), 'Age'])
QE_deldat_IQR_min0.8$SampleType = machine_info[rownames(QE_deldat_IQR_min0.8), 'SampleType']
QE_deldat_IQR_min0.8$DataSet = machine_info[rownames(QE_deldat_IQR_min0.8), 'DataSet']
QE_deldat_IQR_min0.8$Hospital = machine_info[rownames(QE_deldat_IQR_min0.8), 'Hospital']
QE_deldat_IQR_min0.8$Gender = machine_info[rownames(QE_deldat_IQR_min0.8), 'Gender']
QE_deldat_IQR_min0.8 = cbind(QE_deldat_IQR_min0.8[,6360:6364], QE_deldat_IQR_min0.8[,1:6359])
QE_deldat_IQR_min0.8[1:3,1:6]
write.xlsx(QE_deldat_IQR_min0.8, 'D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)

plot(density(na.omit(unlist(QE_deldat_IQR_min0.8[ ,6:ncol(QE_deldat_IQR_min0.8)]))), main="density default")

QE_deldat_IQR_min0.8 = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
QE_deldat_IQR_min0.8 = QE_deldat_IQR_min0.8[intersect(rownames(QE_deldat_IQR_min0.8), rownames(delinfo)),]
dim(QE_deldat_IQR_min0.8)
QE_deldat_IQR_min0.8[1:3,1:6]
delinfo = delinfo[rownames(QE_deldat_IQR_min0.8), ]
unique(delinfo$TrainTest)
table(delinfo$TrainTest)

unique(QE_deldat_IQR_min0.8$Gender)
unique(delinfo$Classificaiton_type)# N, B, M ,Borderline
unique(delinfo$Histopathology_type)# "PTC", "FTC", "PDTC"  "ATC"   "MTC"

QE_deldat_IQR_min0.8_zscore = QE_deldat_IQR_min0.8
QE_deldat_IQR_min0.8_zscore[,6:ncol(QE_deldat_IQR_min0.8_zscore)] = apply(QE_deldat_IQR_min0.8_zscore[,6:ncol(QE_deldat_IQR_min0.8_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})
dim(QE_deldat_IQR_min0.8_zscore)
###### 不校正 ######
QE_deldat[1:10,1:2]
QE_deldat_min0.8 = QE_deldat
min(QE_deldat_min0.8, na.rm = T) * 0.8 # 7.946766
QE_deldat_min0.8[is.na(QE_deldat_min0.8)] = min(QE_deldat_min0.8, na.rm = T) * 0.8
QE_deldat_min0.8 = data.frame(t(QE_deldat_min0.8))
QE_deldat_min0.8[1:10,1:6]
dim(QE_deldat_min0.8)
QE_deldat_min0.8$Age = as.numeric(machine_info[rownames(QE_deldat_min0.8), 'Age'])
QE_deldat_min0.8$SampleType = machine_info[rownames(QE_deldat_min0.8), 'SampleType']
QE_deldat_min0.8$DataSet = machine_info[rownames(QE_deldat_min0.8), 'DataSet']
QE_deldat_min0.8$Hospital = machine_info[rownames(QE_deldat_min0.8), 'Hospital']
QE_deldat_min0.8$Gender = machine_info[rownames(QE_deldat_min0.8), 'Gender']
dim(QE_deldat_min0.8)
QE_deldat_min0.8 = cbind(QE_deldat_min0.8[,6360:6364], QE_deldat_min0.8[,1:6359])

# write.xlsx(QE_deldat_min0.8, 'D:/chh/2023workProject/20240821TTTD/QE/thyroid_QEunQuantile_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
QE_deldat_min0.8 = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QEunQuantile_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
QE_deldat_min0.8[1:5,1:6]

QE_deldat_min0.8_zscore = QE_deldat_min0.8
QE_deldat_min0.8_zscore[,6:ncol(QE_deldat_min0.8_zscore)] = apply(QE_deldat_min0.8_zscore[,6:ncol(QE_deldat_min0.8_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})
dim(QE_deldat_min0.8_zscore) # 6149samp*6359pro
QE_deldat_min0.8_zscore[1:5,1:6]

QE_deldat_min0.8[1:5,1:7]

############ DEswan 1 window  ############
type1 = 'PDTC'
# dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
glm_overlap =list()
for(type1 in c('N', 'B', 'M')){#  "PTC" "FTC" "ATC" "MTC" "PDTC", 'Borderline'
  dat = QE_deldat_min0.8[rownames(delinfo[delinfo$Classificaiton_type %in% c(type1) & delinfo$TrainTest == 'Discovery', ]),]
  # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  if(type1 == 'PTC'){
    dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type %in% c('PTC', "PTMC") & delinfo$TrainTest == 'Discovery', ]),]
    type1 = 'PTC.PTMC' 
  }
  dat = subset(dat, !is.na(Gender))
  #dat = subset(dat, Gender = 'F')
  # dat[6:ncol(dat)] = apply(dat[6:ncol(dat)],2,stdscale)
  dat[1:10,1:7]
  dat = cbind(dat[1:5], class = machine_info[rownames(dat), 'Classificaiton_type'] , dat[,-c(1:5)])
  cor_result = cor(dat[,1], dat[, -c(1:6)], method = 'spearman')
  cor_result[1:10]
  cor_cut = 0
  #temp1 = colnames(cor_result)[abs(cor_result)>cor_cut]
  # dat1 = dat[,5+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>cor_cut]))]
  dat1 = dat[, 6+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>cor_cut]))]
  dat1[1:3,1:3]
  covariates = dat[,c(2,4,5)]
  #covariates = dat[,c(2,4,5, 6)]
  if(type1 == 'PDTC'){
    covariates = dat[c( 5)]
  }
  qt = dat[,1] # Age
  pvalues  <-  NULL
  coefficients  <-  NULL
  qt <- factor(qt)
  for(i in 1:ncol(dat1)){
    deswan.formula  <-  paste(c("dat1[, i] ~ qt", paste("covariates$",colnames(covariates), collapse  =  " + ", sep = "")), collapse = " + ", sep = "")
    # deswan.formula = "dat1[, i] ~ qt + covariates$Gender " # + covariates$SampleType + covariates$Hospital
    test.glm  <- try(glm.fit <- glm(as.formula(deswan.formula), family  =  gaussian), silent=TRUE)
    if(class(test.glm)[1] !=  "try-error"){
      print(c(type1, i, colnames(dat1)[i]))
      glm.res <- car::Anova(glm.fit,  type = "2")
      pvalues <- rbind(pvalues, data.frame(variable = colnames(dat1)[i], factor = rownames(glm.res), pvalue=glm.res$`Pr(>Chisq)`, stringsAsFactors  =  F))
      coefficients <- rbind(coefficients, data.frame(variable  =  colnames(dat1)[i] , factor = names(coefficients(glm.fit)), 
               coefficient=coefficients(glm.fit), stringsAsFactors  =  F))
    }
  }
  pvalues = data.frame(pvalues)
  pvalues_p = data.frame(acast(pvalues[c("variable", "factor", "pvalue")], variable~factor))
  pvalues_padjust = data.frame(apply(pvalues_p, 2, function(x){p.adjust(x, method = 'BH')}))
  temp = subset(pvalues_padjust, qt < 0.05 )
  # if(nrow(temp)>0){
  #   glm_overlap[[type1]] = rownames(temp)
  # }
  head(pvalues, 10)
  coefficients = data.frame(coefficients)
  head(coefficients)
  coefficients1 = acast(coefficients, variable~factor)
  coefficients1 = data.frame(coefficients1)
  coefficients1$qt = rowMeans(coefficients1[grepl('qt', colnames(coefficients1))], na.rm = T)
  coefficients1 = coefficients1[c('qt', colnames(coefficients1)[!grepl('qt', colnames(coefficients1))]) ]
  # coefficients(glm.fit), summary(glm.fit)$coefficients, glm.fit$coefficients
  # pval_coefbind(pvalues, coefficients1, type1)
  pval_coefbind(pvalues, coefficients1, type1, 'D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_')
}

# length(temp1)
# write.xlsx(pvalues, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QE_Discovery_',type1,'_glm_ageContinue_pvalue.xlsx'))
# write.xlsx(coefficients1, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QE_Discovery_',type1,'_glm_ageContinue_coefficients.xlsx'), rowNames  =T)
# write.csv(do.call(cbind, lapply(glm_overlap, function(x) c(x, rep(NA, max(lengths(glm_overlap)) - length(x))))), 
#           "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/glm_Discovery_padjust0.05.csv", row.names = T)

names(glm_overlap)
glm_overlap_dat = do.call(cbind, lapply(glm_overlap, function(x) c(x, rep(NA, max(lengths(glm_overlap)) - length(x)))))
glm_overlap_dat = data.frame(glm_overlap_dat)
write.xlsx(glm_overlap_dat, 'D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QEzscore_Discovery_All_glm_ageContinue_q0.05_upset.xlsx')
############ DEswan multi window ############
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(ggsci)
source('D:/chh/2023workProject/20240821TTTD/DEswan-master/R/DEswan.R')
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/reshape.DEswan.R")
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/nsignif.DEswan.R")
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/q.DEswan.R")

window_Qnum = data.frame(row.names = paste0('Age',unique(QE_deldat_IQR_min0.8_zscore$Age)))
windowDiffer_Qnum = data.frame(row.names = paste0('Age',unique(QE_deldat_IQR_min0.8_zscore$Age)))

type1 = 'N'# N, B, M ,Borderline,,"PTC", "FTC", "PDTC"  "ATC"   "MTC"

DEsan_slide = function(dat, slide, buckets.size, covariates, type1, cor_result){
  res.DEswan=DEswan(data.df = dat[,6+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>.1]))],
   qt = dat[,1], # age
   window.center = seq(min(qt)+1, max(qt), slide), # seq(35,95,10)
   buckets.size = buckets.size,
   covariates = covariates)# SampleType DataSet Hospital Gender
  Sys.time()
  pvals = res.DEswan$p
  
  write.xlsx(res.DEswan$p, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_',slide,'/QEzscore_Discovery_',type1, '_glm_ageContinue_p.xlsx'), rowNames = T)
  # coeff = res.DEswan$coeff
  write.xlsx(res.DEswan$coeff, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_',slide,'/QEzscore_Discovery_',type1, '_glm_ageContinue_coef.xlsx'), rowNames = T)
  res.DEswan.wide.coeff=reshape.DEswan(res.DEswan, parameter = 2, factor = "qt")
  colnames(res.DEswan.wide.coeff) = gsub('X','Age', colnames(res.DEswan.wide.coeff))
  
  res.DEswan.wide.p=reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")# "SampleType" "Hospital"   "Gender" 
  colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
  res.DEswan.wide = res.DEswan.wide.p
  rownames(res.DEswan.wide) = res.DEswan.wide$variable
  colnames(res.DEswan.wide) = paste0(colnames(res.DEswan.wide), '_p')
  for(i in colnames(res.DEswan.wide)[2:ncol(res.DEswan.wide)]){
    res.DEswan.wide[gsub('_p','_padjust',i)] = p.adjust(res.DEswan.wide[,i], method = 'BH')
    res.DEswan.wide[gsub('_p','_coeff',i) ] = res.DEswan.wide.coeff[,gsub('_p','',i)]
    # colnames(res.DEswan.wide)[i] = paste0(colnames(res.DEswan.wide)[i], '_pvalue')
  }
  res.DEswan.wide = res.DEswan.wide[sort(colnames(res.DEswan.wide)[2:ncol(res.DEswan.wide)])]
  write.xlsx(res.DEswan.wide, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_',slide,'/QEzscore_Discovery_',type1 ,'_glm_ageContinue_p.q.coef.xlsx'), rowNames = T)
  res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
  # window_Qnum[colnames(res.DEswan.wide.q.signif), type1] = as.numeric(res.DEswan.wide.q.signif[1,])
  # windowDiffer_Qnum[colnames(res.DEswan.wide.q.signif), buckets.size] = as.numeric(res.DEswan.wide.q.signif[1,])
  
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q )
  toPlot=res.DEswan.wide.q.signif#[1:3,]
  toPlot = melt(toPlot)
  colnames(toPlot) = c('cutoff', 'Age', 'Num')
  toPlot$Age = as.numeric(gsub('Age','', toPlot$Age))
  toPlot$cutoff = paste0('padjust<', toPlot$cutoff)
  p = ggplot(toPlot,aes(Age, Num, color = cutoff))+
    geom_line(size = 1)+
    labs(title = type1, y =paste0( "Significant Num"), x ="qt(Age)")+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),
  legend.position="right",
  axis.text.x = element_text(size = 13, color = "black"),
  axis.text = element_text(size = 13, color = "black"),
  #legend.title = element_blank(),
  legend.text = element_text(size = 13),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank())
  p
  ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_',slide,'/QEzscore_Discovery_',type1 ,'_glm_ageContinue_differQ_Num.pdf'),p)
  
  temp = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]))
  p = ggplot(temp, aes(Age, Num))+
    geom_line(colour = 'red')+
    geom_point(colour = 'black')+
    geom_text_repel( aes( label = paste0('Age',Age,'_',Num)),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
  segment.color = "black", show.legend = FALSE )+
    theme_bw()+labs(title = type1)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right",
          axis.text.x = element_text(size = 13, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          #legend.title = element_blank(),
          legend.text = element_text(size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  p
  ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_',slide,'/QEzscore_Discovery_',type1, '_glm_ageContinue_q0.05_Num.pdf'),p)
  
}
DEsan_slide_gender = function(dat, slide, buckets.size, covariates, type1, Gender){
  starttime = Sys.time()
  res.DEswan=DEswan(data.df = dat[,5+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>.1]))],
 qt = dat[,1], # age
 window.center = seq(min(qt)+1, max(qt), slide), # seq(35,95,10)
 buckets.size = buckets.size,
 covariates = covariates)# SampleType DataSet Hospital Gender
  Sys.time()
  pvals = res.DEswan$p
  # if(class(test.glm)[1] !=  "try-error"){}
  write.xlsx(res.DEswan$p, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/QEzscore_Discovery_',type1,'_sex',Gender,'_glm_ageContinue_p.xlsx'), rowNames = T)
  # coeff = res.DEswan$coeff
  write.xlsx(res.DEswan$coeff, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/QEzscore_Discovery_',type1,'_sex',Gender,'_glm_ageContinue_coef.xlsx'), rowNames = T)
  res.DEswan.wide.coeff=reshape.DEswan(res.DEswan, parameter = 2, factor = "qt")
  colnames(res.DEswan.wide.coeff) = gsub('X','Age', colnames(res.DEswan.wide.coeff))
  
  res.DEswan.wide.p=reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")# "SampleType" "Hospital"   "Gender" 
  colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
  res.DEswan.wide = res.DEswan.wide.p
  rownames(res.DEswan.wide) = res.DEswan.wide$variable
  colnames(res.DEswan.wide) = paste0(colnames(res.DEswan.wide), '_p')
  for(i in colnames(res.DEswan.wide)[2:ncol(res.DEswan.wide)]){
    res.DEswan.wide[gsub('_p','_padjust',i)] = p.adjust(res.DEswan.wide[,i], method = 'BH')
    res.DEswan.wide[gsub('_p','_coeff',i) ] = res.DEswan.wide.coeff[,gsub('_p','',i)]
    # colnames(res.DEswan.wide)[i] = paste0(colnames(res.DEswan.wide)[i], '_pvalue')
  }
  res.DEswan.wide = res.DEswan.wide[sort(colnames(res.DEswan.wide)[2:ncol(res.DEswan.wide)])]
  write.xlsx(res.DEswan.wide, paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/QEzscore_Discovery_',type1,'_sex',Gender,'_glm_ageContinue_p.q.coef.xlsx'), rowNames = T)
  res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
  # window_Qnum[colnames(res.DEswan.wide.q.signif), type1] = as.numeric(res.DEswan.wide.q.signif[1,])
  # windowDiffer_Qnum[colnames(res.DEswan.wide.q.signif), buckets.size] = as.numeric(res.DEswan.wide.q.signif[1,])
  
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q )
  toPlot=res.DEswan.wide.q.signif#[1:3,]
  toPlot = melt(toPlot)
  colnames(toPlot) = c('cutoff', 'Age', 'Num')
  toPlot$Age = as.numeric(gsub('Age','', toPlot$Age))
  toPlot$cutoff = paste0('padjust<', toPlot$cutoff)
  p = ggplot(toPlot,aes(Age, Num, color = cutoff))+
    geom_line(size = 1)+
    labs(title = paste0(type1, ': ', Gender), y =paste0( "Significant Num"), x ="qt(Age)")+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),
  legend.position="right",
  axis.text.x = element_text(size = 13, color = "black"),
  axis.text = element_text(size = 13, color = "black"),
  #legend.title = element_blank(),
  legend.text = element_text(size = 13),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank())
  p
  ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/QEzscore_Discovery_',type1,'_sex',Gender,'_glm_ageContinue_differQ_Num.pdf'),p)
  
  
  # ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QEzscore_Discovery_',type1,'_sex',Gender,'_glm_ageContinue_q0.05_Num.pdf'),p)
  # 
}

type1 = 
for(type1 in c('N', 'M','B', 'Borderline')){ #
  print(type1)
  type1 = 'NM'
  dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type %in% c(type1) & delinfo$TrainTest == 'Discovery', ]),]
  # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  dat = subset(dat, !is.na(Gender))
  # dat = subset(dat, Gender==sex) #'F' "M"
  dat = cbind(dat[1:5], class = machine_info[rownames(dat), 'Classificaiton_type'] , dat[,-c(1:5)])
  cor_result=cor(dat[,1], dat[,-c(1:6)], method = 'spearman')
  
  # cor_cut = 0.1
  # temp1 = colnames(cor_result)[abs(cor_result)>cor_cut]
  length(colnames(cor_result)[abs(cor_result)>.1])
  head(dat[,1:7])
  start_time <- Sys.time()
  qt = dat[,1]
  buckets.size = 10
  #covariates = dat[,c(2,4,5)]
  covariates = dat[,c(2, 4,5,6)]
  if(length(unique(dat$Gender))==1){
    covariates = dat[, c(2,4 )]
  }
  if(type1 == 'PDTC'){
    covariates = dat[ c(5)]
  }
  slide = 1
  DEsan_slide(dat, slide, buckets.size, covariates, type1, cor_result)
}

type1 = 'M'
for(type1 in c('N', 'B', 'Borderline')){ #
  print(type1)
  for(sex in c('F', 'M')){
    dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    dat = subset(dat, !is.na(Gender))
    dat = subset(dat, Gender==sex) #'F' "M"
    dat[1:5,1:7]
    cor_result=cor(dat[,1], dat[,-c(1:5)], method = 'spearman')
    # cor_cut = 0.1
    # temp1 = colnames(cor_result)[abs(cor_result)>cor_cut]
    length(colnames(cor_result)[abs(cor_result)>.1])
    head(dat[,1:7])
    start_time <- Sys.time()
    qt = dat[,1]
    buckets.size = 10
    covariates = dat[,c(2,4,5)]
    if(length(unique(dat$Gender))==1){
      covariates = dat[, c(2,4 )]
    }
    if(type1 == 'PDTC'){
      covariates = dat[ c( 5)]
    }
    slide = 1
    # DEsan_slide(dat, slide, buckets.size, covariates, type1)
    DEsan_slide_gender(dat, slide, buckets.size, covariates, type1, sex)
  }
  
}


for(type1 in c("PTC", "FTC", "PDTC", "ATC", "MTC")){ # "PTMC">PTC
  print(type1)
  # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  for(sex in c('F', 'M')){
    dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    if(type1 == 'PTC'){
      dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type %in% c('PTC', "PTMC") & delinfo$TrainTest == 'Discovery', ]),]
      type1 = 'PTC.PTMC' 
    }
    dat = subset(dat, !is.na(Gender))
    dat[1:10,1:6]
    cor_result=cor(dat[,1], dat[,-c(1:5)], method = 'spearman')
    # cor_cut = 0.1
    # temp1 = colnames(cor_result)[abs(cor_result)>cor_cut]
    length(colnames(cor_result)[abs(cor_result)>.1])
    head(dat[,1:7])
    start_time <- Sys.time()
    qt = dat[,1]
    buckets.size = 10
    covariates = dat[,c(2,4,5 )]
    if(length(unique(dat$Gender))==1){
      covariates = dat[, c(2,4 )]
    }
    if(type1 == 'PDTC'){
      covariates = dat[ c( 5)]
    }
    slide = 1
    # DEsan_slide(dat, slide, buckets.size, covariates,type1)
    DEsan_slide_gender(dat, slide, buckets.size, covariates, type1, sex)
  }
  
}


names(glm_overlap)
######## plot ######

files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/')
files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_5/')
files = files[grepl('_p.xlsx',files)]
files

for(n in c( 2, 10,5) ){
  f = files[n]
  toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/', f), rowNames = T)
  #toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/', f), rowNames = T)
  type1= gsub('QEzscore_Discovery_','',f)
  type1= gsub('_glm_ageContinue_p.xlsx','',type1)
  toReshape=toReshape[which(toReshape$factor=='qt'),-3]
  if (length(unique(toReshape$pvalue))==1){next}
  res.DEswan.wide.p=data.frame(reshape::cast(toReshape,variable~window.center,value = colnames(toReshape)[3]),stringsAsFactors = F)
  colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
  res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
  res.DEswan.wide.q.signif = data.frame(res.DEswan.wide.q.signif)
  if(n == 2){
    temp = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))
    next
  }
  temp1 = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))
  temp = rbind(temp, temp1)
}
unique(temp$Type)
temp$Type = gsub('PTC.PTMC', 'PTC',temp$Type, fixed = T)
temp$Type = gsub('_sexM', '',temp$Type, fixed = T)
temp$Type = gsub('_sexF', '',temp$Type, fixed = T)
unique(temp$Type)
temp1 = temp
temp = temp[temp$Type!='PDTC',]
temp$Type = factor(temp$Type, levels = c("N", "B", "Borderline", "FTC", "PTC", "ATC", "MTC", 'PDTC' ))
temp$Type = factor(temp$Type, levels = c("N", "B","Borderline","FTC","PTC", "ATC", "MTC" ))
temp$Type = factor(temp$Type, levels = c("N", "B","M" ))

p = ggplot(temp, aes(Age, Num, color = Type))+ # "ATC" "B"   "FTC" "M"   "MTC" "N"   "PTC"
  geom_line(size=1.2  )+
  scale_color_d3()+
  geom_point()+
  # geom_text_repel( aes( label = paste0('Age',Age,'_',Num)),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
  #                  segment.color = "black", show.legend = FALSE )+
  theme_bw()+labs(title = '', y = 'q<0.05, Proteins Num')+
  xlim(min(temp$Age), max(temp$Age))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 13, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
# ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QEzscore_Discovery_Nslide1_glm_ageContinue_q0.05_Num.pdf'),p)
p
ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_BMNBorderline_glm_ageContinue_q0.05_Num.1.pdf'),p)

colnames(temp)
unique(temp$Type)
temp_heatmap =data.frame(row.names = unique(temp$Type))
for(i in unique(temp$Type)){
  temp1 = temp[temp$Type == i, ]
  temp_heatmap[i, paste0('Age', temp1$Age)] = temp1$Num
}
temp_heatmap = temp_heatmap[c("N", "B","Borderline","FTC","PTC.PTMC", "ATC", "MTC" ),]
rownames(temp_heatmap) = c("N", "B","Borderline","FTC","PTC", "ATC", "MTC" )
temp_heatmap[temp_heatmap==0] = NA
temp_heatmap = temp_heatmap[sort(colnames(temp_heatmap)) ]
colnames(temp_heatmap) = gsub('Age','',colnames(temp_heatmap))
pheatmap(temp_heatmap, #scale="row",#对行进行归,一化
         color = colorRampPalette(c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A" ))(1000 ), # color参数自定义颜色
         # breaks = c(seq(-2, 2, length=1000)),
         # annotation_col = sampleinfo, # 样本信息
         # annotation_colors = ann_colors ,
         # border_color = NA,
         angle_col = 0,
         fontsize_col = 4, 
         fontsize_row = 15,  
         cluster_rows = F, 
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         fontsize = 10, 
         cellwidth= 5, 
         cellheight= 20)

library(RColorBrewer)
colors <- brewer.pal(9, "RdPu")
colors

###### pheatmap #####

glm_overlap = list()
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['N']] = rownames(subset(temp, qt_padjust<0.05 ))

temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_B_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['B']] = rownames(subset(temp, qt_padjust<0.05 ))

temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_M_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['M']] = rownames(subset(temp, qt_padjust<0.05 ))

olp = intersect(glm_overlap[['N']], glm_overlap[['B']])
olp = intersect(olp, glm_overlap[['M']])
length(olp)


Ntemp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]


temp = QE_deldat_min0.8_zscore[c('Age', olp)]
temp$class = machine_info[rownames(temp), 'Classificaiton_type']
unique(machine_info$TrainTest)
temp$TrainTest = machine_info[rownames(temp), 'TrainTest'] 
temp =  subset(temp, class =='N' & TrainTest== "Discovery" )
temp$Age = as.numeric(temp$Age)
temp = temp[order(temp[,1], decreasing = F), ]
temp =  temp[c( olp)]
temp = data.frame(t(temp))
temp$coeffient = Ntemp[rownames(temp), 'qt_coef']
temp = temp[order(temp['coeffient'], decreasing = F), ]
dim(temp)
temp = temp[1:285]

pheatmap(temp, #scale="row",#对行进行归,一化
         # color = colorRampPalette(c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A" ))(1000 ), # color参数自定义颜色
         color = colorRampPalette(c("blue", "white","red" ))(1000),
         breaks = c(seq(-5, 5, length=1000)),
         #breaks = c(seq(-5, 10, length=100)),
         # annotation_col = sampleinfo, # 样本信息
         # annotation_colors = ann_colors ,
         # border_color = NA,
         angle_col = 0,
         fontsize_col = 4, 
         fontsize_row = 15,  
         cluster_rows = F, 
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         fontsize = 10, 
         cellwidth= 1, 
         cellheight= .5)

#### pick 
files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/')
files = files[grepl('_p.xlsx',files)]
files

# c( 2, 10,5) )
f = files[2]
toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/', f), rowNames = T)
#toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/', f), rowNames = T)
type1= gsub('QEzscore_Discovery_','',f)
type1= gsub('_glm_ageContinue_p.xlsx','',type1)
toReshape=toReshape[which(toReshape$factor=='qt'),-3]
if (length(unique(toReshape$pvalue))==1){next}
res.DEswan.wide.p=data.frame(reshape::cast(toReshape,variable~window.center,value = colnames(toReshape)[3]),stringsAsFactors = F)
colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
# N Age23，Age27，Age80
df = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_N_glm_ageContinue_p.q.coef.xlsx", rowNames = T)
df$pros = rownames(df)
temp1 = subset(df, Age23_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age23_coeff']), ]
temp = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age23_coeff')]
temp$Age = 'Age23'
colnames(temp)[2] = 'coeff'
temp1 = subset(df, Age27_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age27_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age27_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age27'
temp = rbind(temp, temp1)
temp1 = subset(df, Age80_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age80_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age80_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age80'
temp = rbind(temp, temp1)
temp$Type = 'N'

# B Age24，Age58，Age72
df = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_B_glm_ageContinue_p.q.coef.xlsx", rowNames = T)
df$pros = rownames(df)
temp1 = subset(df, Age24_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age24_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age24_coeff')]
temp1$Age = 'Age24'
temp1$Type = 'B'
colnames(temp1)[2] = 'coeff'
temp = rbind(temp, temp1)
temp1 = subset(df, Age58_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age58_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age58_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age58'
temp1$Type = 'B'
temp = rbind(temp, temp1)
temp1 = subset(df, Age72_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age72_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age72_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age72'
temp1$Type = 'B'
temp = rbind(temp, temp1)

# M Age52，Age60，Age70
df = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_M_glm_ageContinue_p.q.coef.xlsx", rowNames = T)
df$pros = rownames(df)
temp1 = subset(df, Age52_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age52_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age52_coeff')]
temp1$Age = 'Age52'
temp1$Type = 'M'
colnames(temp1)[2] = 'coeff'
temp = rbind(temp, temp1)
temp1 = subset(df, Age60_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age60_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age60_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age60'
temp1$Type = 'M'
temp = rbind(temp, temp1)
temp1 = subset(df, Age70_padjust<0.05 )
temp1 = temp1[order(temp1[,'Age70_coeff']), ]
temp1 = temp1[c(1:10, nrow(temp1):(nrow(temp1)- 9)), c('pros', 'Age70_coeff')]
colnames(temp1)[2] = 'coeff'
temp1$Age = 'Age70'
temp1$Type = 'M'
temp = rbind(temp, temp1)
write.csv(temp, 'D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/NBM_absTop10_heatmap_pros.csv', row.names = F)
library(pheatmap)
temp = read.csv('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/NBM_absTop10_heatmap_pros.csv')

type1 = 'M'
n = 3
input1 = temp[temp$Type == type1, ]
samples = intersect(rownames(subset(machine_info, Classificaiton_type == type1)), rownames(QE_deldat_IQR_min0.8_zscore))
mat = QE_deldat_IQR_min0.8_zscore[samples, c( input1[input1$Age == unique(input1$Age)[n], 'pros'])]
info =data.frame( QE_deldat_IQR_min0.8_zscore[samples, 'Age'] ,row.names = samples)
colnames(info) = 'Age'
info$samples = rownames(info)
info = info[order(info[, 1]),]
info = data.frame(info['Age'])
mat = data.frame(t(mat))
mat[1:6,1:2]
pros = rownames(mat)
# pheatmap(mat[rev(pros), rownames(info)], # scale="row",#对行进行归,一化
#          color = colorRampPalette(c("blue", "white","red" ))(1000), 
#          #color = colorRampPalette(c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A" ))(1000 ), # color参数自定义颜色
#          breaks = c( seq(-max(mat)-1 , max(mat)+1 , length=1000)), # max(mat)
#          annotation_col = info, # 样本信息
#          # annotation_colors = ann_colors ,
#          # border_color = NA,
#          angle_col = 0,
#          fontsize_col = 6, 
#          fontsize_row = 4,  
#          cluster_rows = F, 
#          cluster_cols = F,
#          show_rownames = T,
#          show_colnames = F,
#          fontsize = 10, 
#          cellwidth= 0.5,  # B :0.12 N:0.5 M:0.05
#          cellheight= 6, 
#          main = paste0(type1,' ',  unique(input1$Age)[n]))
point_plot = mat[rev(pros), rownames(info)]
point_plot = data.frame(t(point_plot))
point_plot$Age = info$Age
point_plot = aggregate(point_plot, by = list(point_plot$Age), mean, na.rm = T)
point_plot$Age = paste0('Age', point_plot$Age)
rownames(point_plot) = point_plot$Age
point_plot = point_plot[2:ncol(point_plot) ]

point_plot = melt(point_plot, id.vars = 'Age')
head(point_plot)
# point_plot$value
class(point_plot)
typeof(point_plot)
point_plot$Age = gsub('Age', '', point_plot$Age)
point_plot$Age = as.numeric(point_plot$Age)

p = ggplot(point_plot, aes(Age, value, color = variable))+
  # geom_point()+
  geom_smooth(aes(  fill = variable), se = T, method = 'loess', alpha = 0.1)+ #  alpha = 0.15
  labs(y = 'log2 std Expressions', title = paste0(type1,' ',  unique(input1$Age)[n]))+
  theme_classic()+
  scale_color_d3(palette = "category20")+
  scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))
ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/', type1,'_',  unique(input1$Age)[n],'_abstop20.pdf'), p)


######## plot Gender ####
files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/')
files = files[grepl('_p.xlsx',files, fixed = T)]
files

files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/')
filesM = files[grepl('sexM_glm_ageContinue_p.xlsx', files)]
filesM

filesF = files[grepl('sexF_glm_ageContinue_p.xlsx', files)]
filesF

for(n in c(2,4,6)){
  f = filesM[n]
  toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/', f), rowNames = T)
  type1= gsub('QEzscore_Discovery_','',f)
  type1= gsub('_glm_ageContinue_p.xlsx','',type1)
  toReshape=toReshape[which(toReshape$factor=='qt'),-3]
  if (length(unique(toReshape$pvalue))==1){next}
  res.DEswan.wide.p=data.frame(reshape::cast(toReshape,variable~window.center,value = colnames(toReshape)[3]),stringsAsFactors = F)
  colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
  res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
  res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
  res.DEswan.wide.q.signif = data.frame(res.DEswan.wide.q.signif)
  if(n == 2){
    temp = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))
    next
  }
  temp1 = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))
  temp = rbind(temp, temp1)
}

f  = files[9]
toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/', f), rowNames = T)
type1= gsub('QEzscore_Discovery_','',f)
type1= gsub('_glm_ageContinue_p.xlsx','',type1)
toReshape=toReshape[which(toReshape$factor=='qt'),-3]
if (length(unique(toReshape$pvalue))==1){next}
res.DEswan.wide.p=data.frame(reshape::cast(toReshape,variable~window.center,value = colnames(toReshape)[3]),stringsAsFactors = F)
colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
res.DEswan.wide.q.signif = data.frame(res.DEswan.wide.q.signif)
temp = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))

f  = files[8]
toReshape = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1_Gender/', f), rowNames = T)
type1= gsub('QEzscore_Discovery_','',f)
type1= gsub('_glm_ageContinue_p.xlsx','',type1)
toReshape=toReshape[which(toReshape$factor=='qt'),-3]
if (length(unique(toReshape$pvalue))==1){next}
res.DEswan.wide.p=data.frame(reshape::cast(toReshape,variable~window.center,value = colnames(toReshape)[3]),stringsAsFactors = F)
colnames(res.DEswan.wide.p) = gsub('X','Age', colnames(res.DEswan.wide.p))
res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")
res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q, thresholds = 0.05)
res.DEswan.wide.q.signif = data.frame(res.DEswan.wide.q.signif)
temp1 = data.frame(Age = as.numeric(gsub('Age', '',colnames(res.DEswan.wide.q.signif))), Num = as.numeric(res.DEswan.wide.q.signif[1,]), Type = rep(type1, ncol(res.DEswan.wide.q.signif)))
temp = rbind(temp, temp1)
unique(temp$Type )
temp$Type = gsub('.PTMC', '',temp$Type, fixed = T)
sort(unique(temp$Type ), decreasing = T) 
temp$Type = factor(temp$Type, levels = sort(unique(temp$Type ), decreasing = T) )


unique(temp$Type)
temp$Type = gsub('PTC.PTMC', 'PTC',temp$Type, fixed = T)
temp$Type = gsub('_sexM', '',temp$Type, fixed = T)
temp$Type = gsub('_sexF', '',temp$Type, fixed = T)
unique(temp$Type)
temp1 = temp
temp = temp[temp$Type!='PDTC',]
# temp$Type = factor(temp$Type, levels = c("N", "B", "Borderline", "FTC", "PTC", "ATC", "MTC", 'PDTC' ))
# temp$Type = factor(temp$Type, levels = c("N", "B","Borderline","FTC","PTC", "ATC", "MTC" ))
temp$Type = factor(temp$Type, levels = c("N", "B", "M" )) # "Borderline",


p = ggplot(temp, aes(Age, Num, color = Type))+
  geom_line(size = 1 )+
  scale_color_d3()+theme_bw()+
  labs(title = 'M')+ # strsplit(type1, '_')[[1]][1] 'PTC'
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p
p+geom_point( )+
  geom_text_repel(aes( label = paste0('Age',Age,'_',Num)),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
                   segment.color = "black", show.legend = FALSE )
######### Age class barplot #####
classdat = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_NM_glm_ageContinue.xlsx", rowNames = T )
classdat = subset(classdat, class_padjust < 0.05)
Agedat = read.xlsx( "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_M_glm_ageContinue_p.q.coef.xlsx" , rowNames = T )
Agedat = Agedat[grepl('adjust', colnames(Agedat))]
temp = c()
temppros = list()
for(i in colnames(Agedat)){
  temp1 = rownames(Agedat[ Agedat[i]<0.05,])
  n = setdiff(temp1, rownames(classdat))
  m = intersect(temp1, rownames(classdat))
  #print(c(i, nrow(temp1), m))
  temp = rbind(temp, c(i, length(n), length(m)))
  if(length(n)!=0){temppros[[gsub('padjust','only',  i)]] = n}
  if(length(m)!=0){ temppros[[gsub('padjust','overlap', i )]] = m}
}

write.xlsx(data.frame(do.call(cbind, lapply(temppros, function(x) c(x, rep(NA, max(lengths(temppros)) - length(x)))))), 
           "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/QEzscore_Discovery_NM_glm_ageContinue_overlapPros.xlsx", rowNames = F)


temp = data.frame(temp)
colnames(temp) = c('Age', 'Age Only', 'overlap')
temp$variable
temp = melt(temp, id.vars = 'Age')
temp$value = as.numeric(temp$value)
temp$variable = factor(temp$variable, levels = c('overlap', 'Age Only' ))
temp$Age = gsub('_padjust', '', temp$Age )
temp$Age = as.numeric(gsub('Age', '', temp$Age ))
ggplot(temp,aes(Age, value, fill = variable))+
  #geom_line(size = 1)+
  scale_fill_d3()+
  geom_bar(stat = "identity")+
  labs(title = 'NM', y =paste0( ""), x ="Age")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   legend.position="right",
                   axis.text.x = element_text(size = 13, color = "black", vjust =  0.5),
                   axis.text = element_text(size = 13, color = "black"),
                   #legend.title = element_blank(),
                   legend.text = element_text(size = 13),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
p
######### upset #####
library(ggplot2)
library(ggsci)
library(UpSetR)
# overlap1:  N，B，Borderline，M 
# overlap2:  N，B，Borderline，FTC，PTC，PDTC，ATC, MTC
input = glm_overlap[c('B','N','Borderline','M')]
input = glm_overlap[c('B','N','Borderline','FTC','PTC','PDTC','ATC', 'MTC')]
#pdf(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QEzscore_Discovery_All_glm_ageContinue_q0.05_upset.pdf'), width = 13)

input = {}
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_All_glm_ageContinue_q0.05_upset.xlsx')
colnames(temp)
tempM = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_M_glm_ageContinue.xlsx', rowNames = T)
input[['N']] = temp[!is.na(temp$N), 'N']
input[['B']] = temp[!is.na(temp$B), 'B']
#input[['M']] = rownames(tempM[tempM$qt_padjust<0.05,  ] )
input[['FTC']] = temp[!is.na(temp$FTC), 'FTC']
input[['PTC']] = temp[!is.na(temp$PTC), 'PTC']
input[['ATC']] = temp[!is.na(temp$ATC), 'ATC']
input[['MTC']] = temp[!is.na(temp$MTC), 'MTC']


p1 = intersect(input[[1]], input[[2]])
p1 = intersect(p1, input[[3]])
data.frame(p1)

upset(fromList(input), 
      order.by = "freq",  # 主坐标系排序
      nsets = 16,
      number.angles = 0,  # 柱标倾角
      point.size = 3,  # 点大小
      line.size = 1,  # 线粗细
      #sets.x.label = "Datasets Size",  # x 标题
      set_size.show = T,
      #main.bar.color = "gray",
      #sets.bar.color = "gray",
      mainbar.y.label = "Count of Intersection",  # y 标题
      sets.x.label = "Datasets Size",  # x 标题
      text.scale = c(2, 1.5, 1, 1.5, 1.5, 1.5)
      # y 标题大小，y 刻度标签大小，datasetSize标题大小，datasetSize刻度标签大小，datasetSize分类标签大小，柱数字大小
)
dev.off()
dev.off()
dev.off()


files = list.files("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/")
files = files[grepl('xlsx',files, fixed = T)]
files

coef_pos = list()
coef_neg = list()
for(type1 in c('N', "B","FTC", "PTC", "ATC", "MTC")){
  #type1 = 'N'
  coefficients1 = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_',type1,'_glm_ageContinue.xlsx'), rowNames  =T)
  #p1 = glm_overlap[[type1]]
  temp = coefficients1#[p1, ]
  pros = rownames(subset(temp, qt_coef>0 & qt_padjust< 0.05 ))
  pros = sapply(pros, function(x){strsplit(x, '_')[[1]][1]})
  coef_pos[[type1]] = as.character(pros)
  
  pros = rownames(subset(temp, qt_coef<0 & qt_padjust< 0.05 ))
  pros = sapply(pros, function(x){strsplit(x, '_')[[1]][1]})
  coef_neg[[type1]] = pros
}


write.csv(do.call(cbind, lapply(coef_pos, function(x) c(x, rep(NA, max(lengths(coef_pos)) - length(x))))), 
          "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/glm_Discovery_padjust0.05_Coef_pos.csv", row.names = T)
write.csv(do.call(cbind, lapply(coef_neg, function(x) c(x, rep(NA, max(lengths(coef_neg)) - length(x))))), 
          "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/glm_Discovery_padjust0.05_Coef_neg.csv", row.names = T)


write.xlsx(data.frame(do.call(cbind, lapply(coef_pos, function(x) c(x, rep(NA, max(lengths(coef_pos)) - length(x)))))), 
          "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/NBFPAM_padjust0.05_allWindowUpCoef.xlsx", rowNames = F)
write.xlsx(data.frame(do.call(cbind, lapply(coef_neg, function(x) c(x, rep(NA, max(lengths(coef_neg)) - length(x)))))), 
          "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/NBFPAM_padjust0.05_allWindowdownCoef.xlsx", rowNames = F)


pros = intersect(glm_overlap[['PTC']], intersect(glm_overlap[['PDTC']], intersect(glm_overlap[['FTC']], glm_overlap[['B']])))
pros = setdiff(pros, glm_overlap[['N']])
pros = setdiff(pros, glm_overlap[['MTC']])
pros = setdiff(pros, glm_overlap[['ATC']])
pros = setdiff(pros, glm_overlap[['Borderline']])
length(pros)
data.frame(pros)

type1 = 'N'
coefficients1 = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_',type1,'_glm_ageContinue.xlsx'), rowNames  =T)
temp1 = coefficients1[p1, c('qt_coef', 'qt_padjust')]
temp1$type = type1
temp1$pro = rownames(temp1)
temp1 = temp1[order(temp1[,1]),]
temp = rbind(temp, temp1)
colnames(temp)
temp$type = factor(temp$type, levels = c('N','B','M'))
temp$pro = factor(temp$pro, levels = temp1$pro)

temppathway = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_NBM_upset41_overlap_pathways.xlsx', rowNames = T)
pros = sapply(temp1$pro, function(x){
  x = strsplit(x, '_')[[1]][2]
  x = gsub('.','-',x, fixed = T)
})
# as.character(pros)
pros = rev(as.character(pros))
pros[pros%in% rownames(temppathway)]
temppathway = temppathway[pros,]

custom_palette <- colorRampPalette(c('#BCF2F6', '#08C2FF','#478CCF','#006BFF','black', "#AE017E" ,"#DD3497", "#F768A1", "#FA9FB5", "#FCC5C0"))(1000)

head(temp)
ggplot(temp, aes(x = type, y = pro, color= qt_coef , size = -log10(qt_padjust))) +
  geom_point(  alpha=0.7 )+ # stroke = 0: 设置边框宽度为 0，以确保没有边框显示
  scale_size(range = c(3, 8) )+
  # scale_color_gradient(low = "#546de5", high = "#ff4757") +
  scale_color_gradientn(colors = custom_palette) +
  # scale_color_gradient2(low = "#546de5", mid = 'black', high = "yellow")+#ff4757
  #scale_color_viridis_c()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

library(pheatmap)
#temppathway[temppathway==0] = NA
pheatmap(temppathway, # scale="row",#对行进行归一化
         color = c('white', "black"), # color参数自定义颜色#ff4757
         # border_color = NA,
         # breaks = c(seq(-2, 2, length=1000)),
         fontsize_col = 12, 
         legend = F,
         # fontsize_row = 1,  
         angle_col = 45,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         fontsize = 10, 
         cellwidth= 16, 
         cellheight=  10)
####### volcano #####
library(ggrepel)
unique(delinfo$Classificaiton_type)# N, B, M ,Borderline
unique(delinfo$Histopathology_type)# "PTC", "FTC", "PDTC"  "ATC"   "MTC"
type1 = 'N'
# dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = subset(dat, !is.na(Gender))
# 
files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/')
files = files[grepl('_glm_ageContinue.xlsx',files)]
files

glm_overlap = list()
for(type1 in c('N', 'B', 'M')){ # "PTC.PTMC", "FTC", "PDTC",  "ATC",   "MTC"
  #coefficients1 = read.xlsx( paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_',type1,'_glm_ageContinue.xlsx'), rowNames  =T)
  type1 = 'M'
  type1 = 'B'
  # coefficients1 = read.xlsx( paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_',type1,'_glm_ageContinue.xlsx'), rowNames  =T)
  coefficients1 = read.xlsx( paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_',type1,'_glm_ageContinue.xlsx'), rowNames  =T)
  
  colnames(coefficients1)
  # temp = data.frame(coefficient = coefficients1$GenderM_coef, FDR = coefficients1$Gender_padjust, 
  temp = data.frame(coefficient = coefficients1$classN_coef, FDR = coefficients1$class_padjust, 
  proteins = sapply(rownames(coefficients1), function(x){
     x = strsplit(x, '_')[[1]][2]
     x = gsub('.','-',x, fixed = T)
     x
   })
    )
  temp$coefficient = 0 - temp$coefficient
  temp = subset(temp, FDR<0.05)
  # glm_overlap[[paste0(type1,'Gender')]] = rownames(temp )
  temp = data.frame(coefficient = coefficients1$qt_coef, FDR = coefficients1$qt_padjust,
                    proteins = sapply(rownames(coefficients1), function(x){
                      x = strsplit(x, '_')[[1]][2]
                      x = gsub('.','-',x, fixed = T)
                      x
                    })
  )
  temp = subset(temp, FDR<0.05)
  glm_overlap[[ type1 ]] = rownames(temp )
  
  
  temp$Sig = ifelse(is.na(temp$coefficient), "None", ifelse(temp$coefficient> 0 & temp$FDR<0.05, "Up",
        ifelse( temp$coefficient< 0& temp$FDR<0.05, "Down", "None")))
  table(temp$Sig)
  p = ggplot(temp, aes(x = coefficient, y = -log10(FDR), colour=Sig)) +
    geom_point(alpha=0.7, size=2 ) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
    geom_text_repel( data = temp[temp$FDR<0.05&abs(temp$coefficient)> 0, ],
  aes(x = coefficient, y = -log10(FDR), label = proteins),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
  segment.color = "black", show.legend = FALSE )+
    labs(x="coefficient", y="-log10 (P-adjust)", title = type1  )+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right",
          axis.text.x = element_text(size = 15, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          #legend.title = element_blank(),
          legend.text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  p
  ggsave(paste0('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_glm_',type1,'_volcano.pdf'),p)
  
}
names(glm_overlap)
library(UpSetR)
upset(fromList(glm_overlap[c(3:6)]), # 
      order.by = "freq",  # 主坐标系排序
      nsets = 16,
      number.angles = 0,  # 柱标倾角
      point.size = 3,  # 点大小
      line.size = 1,  # 线粗细
      #sets.x.label = "Datasets Size",  # x 标题
      set_size.show = T,
      #main.bar.color = "gray",
      #sets.bar.color = "gray",
      mainbar.y.label = "Count of Intersection",  # y 标题
      sets.x.label = "Datasets Size",  # x 标题
      text.scale = c(2, 1.5, 1, 1.5, 1.5, 1.5)
      # y 标题大小，y 刻度标签大小，datasetSize标题大小，datasetSize刻度标签大小，datasetSize分类标签大小，柱数字大小
)
p1 = intersect(glm_overlap[[1]], glm_overlap[[2]])
p1 = intersect(p1, glm_overlap[[3]])
p1 = intersect(p1, glm_overlap[[4]])
p1 


library(ggVennDiagram)
x = glm_overlap[c(4,3,12, 11 )] # c(1,2,9,10 )
ggVennDiagram(x, #category.names = c("A","B","C"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "dashed", 
              edge_size = 1) + labs(title = 'NM')+
  # labs(title = names(x)[1])+
  scale_fill_gradient(low="white",high = "#b9292b",name = "count")#library(ggplot2)
######## pheatmap ######
library(pheatmap)
pros = read.csv("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/glm_Discovery_padjust0.05.csv", row.names = 1)
pros = unique(na.omit(unlist(pros)))
pros = pros[!is.na(pros)]


# N, B, M ,Borderline,,"PTC", "FTC", "PDTC"  "ATC"   "MTC"
cla = subset(delinfo, Classificaiton_type == c('N', 'B', 'M' ,'Borderline') & TrainTest == 'Discovery')
cla =  rownames(cla)

histp = subset(delinfo, Histopathology_type == c("PTC", "FTC", "PDTC",  "ATC",   "MTC") & TrainTest == 'Discovery')
histp =  rownames(histp)

QE_deldat_IQR_min0.8[1:3,1:6]

dat = QE_deldat_IQR_min0.8[ unique(cla, histp), 6:ncol(QE_deldat_IQR_min0.8)]
dat = QE_deldat_IQR_min0.8[ 6:ncol(QE_deldat_IQR_min0.8)]

sampleinfo =delinfo[rownames(dat), c('Classificaiton_type','Histopathology_type', 'SampleType','Age','Gender')] #QE_deldat_IQR_min0.8_zscore[ unique(cla, histp), c(1,2,3 ,5)]
sampleinfo$Histopathology_type = gsub('PTMC', 'PTC', sampleinfo$Histopathology_type)
sampleinfo = sampleinfo[sampleinfo$Histopathology_type %in% c("N","HT","FA","MNG","UMP","NIFTP","FTC","PTC", "PDTC", "ATC", "MTC"),]
# sampleinfo$Histopathology_type = factor(sampleinfo$Histopathology_type, levels = )
sampleinfo$Age = as.numeric(sampleinfo$Age)
sampleinfo = sampleinfo[order(sampleinfo[,4]),]
sampleinfo = subset(sampleinfo, !is.na(Gender))
temp = sampleinfo[sampleinfo$Histopathology_type == 'N',]
for(i in c( "HT","FA","MNG","UMP","NIFTP","FTC","PTC", "PDTC", "ATC", "MTC")){
  temp = rbind(temp, sampleinfo[sampleinfo$Histopathology_type == i,])
}
sampleinfo = temp
mycol = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#8f7f7f',
          '#bcbd22', '#17becf',"#708090",'#68A180',
          '#F3B1A0', '#D6E7A3',"#99CCCC","#66CC99",
          "yellow", "blue", "gray")

SampleType <- mycol[1:length(unique(sampleinfo$SampleType))]
names(SampleType) <- unique(sampleinfo$SampleType)
Classificaiton_type <- mycol[1:length(unique(sampleinfo$Classificaiton_type))]
names(Classificaiton_type) <- unique(sampleinfo$Classificaiton_type)
Histopathology_type <- mycol[1:length(unique(sampleinfo$Histopathology_type))]
names(Histopathology_type) <- unique(sampleinfo$Histopathology_type)
Gender <- mycol[1:length(unique(sampleinfo$Gender))]
names(Gender) <- unique(sampleinfo$Gender)

ann_colors <- list( Classificaiton_type = Classificaiton_type, SampleType =  SampleType,
 Histopathology_type = Histopathology_type, Gender = Gender )
ann_colors


pheatmap(t(dat[rownames(sampleinfo),]), # scale="row",#对行进行归一化
         color = colorRampPalette(c("blue", "white","red" ))(1000), # color参数自定义颜色
         # breaks = c(seq(-2, 2, length=1000)),
         annotation_col = sampleinfo, # 样本信息
         annotation_colors = ann_colors ,
         fontsize_col = 0.55, 
         fontsize_row = 1,  
         cluster_rows = T,
         labels_col = as.character(sampleinfo$Age),
         cluster_cols = F,
         show_rownames = F,
         show_colnames = T,
         fontsize = 7, 
         cellwidth= 0.1, 
         cellheight= 0.05)

######### mfuzz ######
type1 = 'N'
coefficients = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_N_glm_ageContinue.xlsx", rowNames = T)

pros = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_N_glm_ageContinue.xlsx", rowNames = T)
pros = rownames(pros[pros$qt_padjust<0.05,])

dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
dat = subset(dat, !is.na(Gender))

library(Mfuzz)
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycolor <- colorRampPalette(mycol)(100)

mfuzzInput = aggregate(dat[c(1, 6:ncol(dat))], by = list(dat$Age), mean, na.rm = T)
dim(mfuzzInput)
mfuzzInput[1:4,1:3]
rownames(mfuzzInput) = paste0('Age', mfuzzInput$Age)
mfuzzInput = data.frame( t(mfuzzInput[3:ncol(mfuzzInput)]))


mat <- as.matrix(mfuzzInput[pros,])
dt <- new("ExpressionSet", exprs = mat)
# dt <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt , mode="mean")
tmp <- filter.std(dt.f, min.std=0)
dt.s <- standardise(tmp)
m1 <- mestimate(dt.s)
m1
c =6
set.seed(2024)
temp = dt.s@assayData
cl = mfuzz(dt.s, c = c, m = m1 )
# pdf(paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_",type1,'_',c,"cluster_noslide_glm_ageContinue_Mfuzz.pdf"), height = 15,width = 10)
# mfuzz.plot(dt.s,cl,mfrow=c(5,4), new.window= FALSE,
#            time.labels = colnames(dt.s), colo = mycolor)

mfuzz.plot2(dt.s, cl, #mfrow=c(5,4), 
            mfrow=c(2,3), 
            time.labels = colnames(mfuzzInput),
            centre=TRUE, centre.col="black", # 中心线
            x11 = F, # 新窗口
            xlab = "Age")
dev.off()
dev.off()
dev.off()

protein_cluster <- cl$cluster
protein_cluster <- cbind( protein_cluster, mat[names(protein_cluster), ] )
protein_cluster = data.frame(protein_cluster)
# write.xlsx(protein_cluster, paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_", type1,'_',c, "_MfuzzCluster.xlsx"),rowNames = T)
write.csv(protein_cluster, paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_", type1,'_',c, "_MfuzzCluster.csv") )
member = data.frame(cl$membership)
write.xlsx(member, paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_", type1,'_',c, "_MfuzzMembership.xlsx"),rowNames = T)

dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
dat = subset(dat, !is.na(Gender))
nrow(member[member$X2>0.5,])
protein_cluster = data.frame(protein_cluster)
for(clust in 1:6){
  pros_clust = rownames(protein_cluster[protein_cluster$protein_cluster == clust, ])
  corinput = dat[c('Age', pros_clust)]
  cor_result = cor(corinput[,'Age'], corinput[2:ncol(corinput)], method = 'spearman')
  max(cor_result)
  # names(cor_result)[1:10]
  boxplot_input = data.frame(cor_result)
  boxplot_input = data.frame(t(boxplot_input))
  colnames(boxplot_input) = 'spearmanr'
  boxplot_input$coef = coefficients[rownames(boxplot_input), 'qt_coef']
  # median(boxplot_input$spearmanr)
  # median(boxplot_input$coef)
  meddata = data.frame(apply(boxplot_input, 2, function(x){median(x)}))
  meddata$variable = rownames(meddata)
  colnames(meddata) = c('value', 'variable' )
  # boxplot_input$pros = rownames(boxplot_input)
  boxplot_input_melt = melt(boxplot_input)
  p = ggplot(boxplot_input_melt, aes(variable, value))+
    geom_boxplot()+
    # geom_jitter()+
    geom_text(aes(label = round(value, 3)), data = meddata, 
              vjust= -2.5, hjust = 0.5, size = 4, color = 'red')+
    labs(title = paste0('Cluster',clust), x = '')+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5) )
  ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_N_noslide_Mfuzz6_cluster",clust,".pdf"), p)
  
}

clust = 1
pros_clust = rownames(protein_cluster[protein_cluster$protein_cluster == clust, ])

dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
dat = subset(dat, !is.na(Gender))
# dat
# scale(x, center = TRUE, scale = TRUE)

ind = which.max(coefficients[rownames(member) , 'qt_coef'])
coef = coefficients[rownames(member)[ind] , 'qt_coef']
coef
point_plot = dat[c('Age', rownames(member)[ ind])] # 
colnames(point_plot) = c('Age', "Protein")

ggplot(point_plot, aes(Age, Protein))+
  # geom_line()
  # geom_point()+
  geom_smooth()+ 
  labs(y = 'log2 std Expressions')+
  theme_classic()+ 
  theme(axis.text = element_text(size = 13),
        text = element_text(size = 13))

point_plot1 = melt(dat[c('Age', rownames(member) )],id.vars = 'Age' )
colnames(point_plot1) = c('Age', "Protein", 'Intensity')
ggplot(point_plot1, aes(Age, Intensity, fill = Protein))+
  geom_smooth( se = FALSE, color = '#424242', alpha= 1, linewidth= 0.3 )+#method = "glm", aes( fill = Protein)
  theme_classic()+labs(y = 'log2 std Expressions')+
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        text = element_text(size = 13))

########### 甲状腺相关蛋白重点 #########
pro_list = c('P01266','P07858','P07202', 'P12830','P08670','Q14195', #甲状腺
             'P10828', 'P20396', 'P01222', 'P16473')# TSH

pro_list = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/Aging_related_overlap125_paper.xlsx", colNames = F)
pro_list = pro_list$X1
pro_lists = c()
for(i in pro_list){
  pro_lists = c(pro_lists, colnames(QE_deldat_IQR_min0.8_zscore)[ grepl(i, colnames(QE_deldat_IQR_min0.8_zscore))])
}
pro_lists = unique(c(pro_lists, "P01266_TG","P07858_CTSB","P07202_TPO","P12830_CDH1","P08670_VIM","Q14195_DPYSL3", "P16473_TSHR"))

files = list.files("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/")
files = files[grepl('ageContinue.xlsx',files, fixed = T)]
files
temp2 = data.frame(row.names = pro_lists)
for(f in files){
  temp = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/", f), rowNames = T)
  type1  = gsub('QEzscore_Discovery_', '', gsub('_glm_ageContinue.xlsx', '', f, fixed = T))
  temp1 = temp[pro_lists, grepl('qt', colnames(temp))]
  temp1 = subset(temp1, !is.na( qt_coef))
  if(nrow(temp1)==0){next}
  colnames(temp1) = gsub('qt',type1, colnames(temp1))
  temp2[rownames(temp1), colnames(temp1)] = temp1
}
temp2 = temp2[rowSums(is.na(temp2))!=ncol(temp2), ]
write.xlsx(temp2, "D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QEzscore_Discovery_impPros_Age_p.q.coef.xlsx", rowNames = T)


pro_lists = c('O76074_PDE5A','P06396_GSN','P07942_LAMB1','P31939_ATIC','P52655_GTF2A1','P56937_HSD17B7')

for(pro in pro_lists){
  type1 = 'N'
  dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  dat = subset(dat, !is.na(Gender))
  point_plot = dat[c('Age', pro)] # 
  colnames(point_plot) = c('Age', "Protein")
  point_plot$Type = type1
  dim(point_plot)
  for(type1 in c('B', 'Borderline', 'M')){
    dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    dat = subset(dat, !is.na(Gender))
    temp = dat[c('Age', pro)] # 
    colnames(temp) = c('Age', "Protein")
    temp$Type = type1
    point_plot = rbind(point_plot, temp)
  }
  # for(type1 in c("PTC", "FTC", "PDTC", "ATC", "MTC")){
  #   dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  #   if (type1 == "PTC"){
  #     dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Histopathology_type %in% c("PTC", 'PTMC') & delinfo$TrainTest == 'Discovery', ]),]
  #   }
  #   dat = subset(dat, !is.na(Gender))
  #   temp = dat[c('Age', pro)] # 
  #   colnames(temp) = c('Age', "Protein")
  #   temp$Type = type1
  #   point_plot = rbind(point_plot, temp)
  # }
  # point_plot$Type = factor(point_plot$Type, levels = c('N', 'B', 'Borderline', "PTC", "FTC", "PDTC", "ATC", "MTC"))
  point_plot$Type = factor(point_plot$Type, levels = c('N', 'B', 'Borderline', "M"))
  p = ggplot(point_plot, aes(Age, Protein, color = Type))+
    geom_smooth(aes(fill = Type), method = 'loess',alpha = 0.15)+ #  se = F
    # scale_fill_d3()+
    labs(y = 'log2 std Expressions', title = pro)+
    theme_classic()+ 
    theme(axis.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15))
  p
  # ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QEzscore_Discovery_impPros_",pro,'.pdf'), p)
  
  ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QEzscore_Discovery_impPros_NBBM_",pro,'.pdf'), p)
}


#'TG','TPO', 'TSHB'
#colnames(dat480_deldat_min08)[ grepl(paste0( 'TSHR$'), colnames(dat480_deldat_min08) )]
pro_lists = c("P01266_TG","P07202_TPO","P16473_TSHR")

QE_deldat[1:10,1:2]

point_plots = c()
for(pro in pro_lists){
  type1 = 'N'
  # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  # dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  dat = QE_deldat_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
  dat = subset(dat, !is.na(Gender))
  point_plot = dat[c('Age', pro)] # 
  colnames(point_plot) = c('Age', "value")
  if(sum(point_plot$value == min(point_plot$value ))>5){
    point_plot[point_plot$value == min(point_plot$value ) , 'value'] = NA
  }
  # point_plot$value = as.numeric(QE_deldat_min0.8_zscore[rownames(point_plot), pro]) # as.numeric(QE_deldat_min0.8_zscore[pro, rownames(point_plot)])
  point_plot$Type = type1
  point_plot$Protein = pro
  # point_plot[point_plot$value == unique(point_plot[table(point_plot$value)>5, 'value']), 'value'] = NA
  dim(point_plot)
  point_plots = rbind(point_plots, point_plot)
  for(type1 in c('B',  'M')){
    # dat = QE_deldat_IQR_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    # dat = QE_deldat_IQR_min0.8[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    dat = QE_deldat_min0.8_zscore[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
    dat = subset(dat, !is.na(Gender))
    point_plot = dat[c('Age', pro)] # 
    colnames(point_plot) = c('Age', "value")
    # point_plot$value = as.numeric(QE_deldat_min0.8_zscore[rownames(point_plot), point_plots[point_plots$Type == type1 & point_plots$Protein == p_cur , ]]) # as.numeric(QE_deldat_min0.8_zscore[pro, rownames(point_plot)])
    point_plot$Type = type1
    point_plot$Protein = pro
    if(sum(point_plot$value == min(point_plot$value ))>5){
      point_plot[point_plot$value == min(point_plot$value ) , 'value'] = NA
    }
    # point_plot[point_plot$value == unique(point_plot[table(point_plot$value)>10, 'value']), 'value'] = NA
    point_plots = rbind(point_plots, point_plot)
  }
}



colors = c('#46AA6E' , "#FF7F0E" , '#9C7822')
names(colors ) = c('N','B', 'M')
colors
type1 ='N' # 'B' # 'M'# 
for (p_cur in pro_lists){
  p = ggplot(point_plots[point_plots$Type == type1 & point_plots$Protein == p_cur , ], aes(Age, value, color = Protein))+
    geom_smooth(aes(fill = Protein), method = 'loess',alpha = 0.15)+ #  se = F
    geom_point(alpha=0.2)+
    xlim(c(min(point_plots$Age), max(point_plots$Age)))+
    #ylim(c(min(point_plots$value), max(point_plots$value)))+
    # scale_fill_d3()+
    scale_fill_manual(values = c(colors[[type1]]))+ #'#46AA6E' , "#FF7F0E" , '#9C7822'
    scale_color_manual(values = c(colors[[type1]]))+
    labs(y = 'log2 std Expressions', title = paste0(type1,' ', p_cur ))+
    # labs(y = 'Log2 Expressions', title = paste0(type1,' ', p_cur ))+
    theme_classic()+ 
    #scale_y_continuous(n.breaks = 8)+
    scale_x_continuous(n.breaks = 8)+
    theme(axis.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15))
  # ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QE_Discovery_impPros_",type1,"_",p_cur,'.pdf'), p)
  ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QEunnormalZscore_Discovery_impPros_",type1,"_",p_cur,'.pdf'), p)
}

for (j in pro_lists){
  N = point_plots[point_plots$Type == 'N' & point_plots$Protein == j , 'value']
  M = point_plots[point_plots$Type == 'M' & point_plots$Protein == j, 'value']
  B = point_plots[point_plots$Type  == 'B' & point_plots$Protein == j, 'value']
  NM = wilcox.test(N, M)
  NB= wilcox.test(N, B)
  BM= wilcox.test(B, M)
  print(c(j, NM$p.value, NB$p.value, BM$p.value))
}

glm_overlap = list()
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['N']] = rownames(subset(temp, qt_padjust<0.05 ))

temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_B_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['B']] = rownames(subset(temp, qt_padjust<0.05 ))

temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_M_glm_ageContinue.xlsx', rowNames = T)
temp[pro_lists,1:3]
glm_overlap[['M']] = rownames(subset(temp, qt_padjust<0.05 ))

write.xlsx(data.frame(do.call(cbind, lapply(glm_overlap, function(x) c(x, rep(NA, max(lengths(glm_overlap)) - length(x)))))),
          "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_NBM_padjust0.05.xlsx", rowNames = F)

library(VennDiagram)
df_inter <- get.venn.partitions(glm_overlap)
for (i in 1:nrow(df_inter)) df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ', ')
# df_inter[-c(5, 6)]
write.xlsx(data.frame(df_inter), 
           "D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_NBM_padjust0.05.xlsx", rowNames = F)



olp = intersect(glm_overlap[['N']], glm_overlap[['B']])
olp = intersect(olp, glm_overlap[['M']])
length(olp)

ggVennDiagram(glm_overlap, #category.names = c("A","B","C"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "dashed", 
              edge_size = 1) + labs(title = 'NBM')+
  # labs(title = names(x)[1])+
  scale_fill_gradient(low="white",high = "#b9292b",name = "count")#library(ggplot2)

# ggsave(paste0("D:/chh/2023workProject/20240821TTTD/QE/ImpPros/QEzscore_Discovery_impPros_",pro,'.pdf'), p)

########### gsva ############
anova_func <- function(exp, type){
  Pvalue = c()
  for(i in 1:nrow(exp)){
    #if (  sum(!(is.na(exp[i,  1:a ]))) <2 | sum(!(is.na(exp[i,  (a+1):(a+b)]))) <2 | sum(!(is.na(exp[i,  (a+b+1):(a+b+c)]))) <2 | sum(!(is.na(exp[i,  (a+b+c+1):(a+b+c+d)]))) <2){
    #  Pvalue = c(Pvalue, NA)
    #  next
    #}
    y = try(aov(as.numeric(exp[i,  ]) ~ type), silent=FALSE)
    if('try-error' %in% class(y)) # 判断当前循环的try语句中的表达式是否运行正确
    {
      Pvalue = c(Pvalue, NA)  # 此处可以对运行错误的情况进行处理应对
    }else{
      y = aov(as.numeric(exp[i,  ]) ~ type)# 默认var.equal = FALSE
      if (dim(summary(y)[[1]])[2] != 5){
        Pvalue[i]<- NA
        next
      }
      Pvalue[i]<- summary(y)[[1]][,5][1]
    }
  }
  FDR=p.adjust(Pvalue, "BH")
  out<-cbind(exp, Pvalue, FDR )
  out
}


library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSVA)
CP <- msigdbr(species = "Homo sapiens", category = "H" ) #%>% dplyr::select(gs_subcat, gs_name, entrez_gene, gene_symbol, human_gene_symbol )
head(CP)
unique(CP$gs_subcat)
unique(CP$gs_description)
#CP = subset(CP, gs_subcat!="CGP")

entriz = CP[ c("gs_name","gene_symbol")] %>% as.data.frame() # CP$human_gene_symbol %in% ttestout$symbol,
entriz <- split(entriz$gene_symbol, entriz$gs_name)

glm_overlap = list()
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
glm_overlap[['N']] = rownames(subset(temp, qt_padjust<0.05 ))
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_B_glm_ageContinue.xlsx', rowNames = T)
glm_overlap[['B']] = rownames(subset(temp, qt_padjust<0.05 ))
temp = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormal_Discovery_M_glm_ageContinue.xlsx', rowNames = T)
glm_overlap[['M']] = rownames(subset(temp, qt_padjust<0.05 ))
library(VennDiagram)
df_inter <- get.venn.partitions(glm_overlap)
for (i in 1:nrow(df_inter)) df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ', ')

NBMolp = strsplit(df_inter[1, "values"], ', ')[[1]]
NBMolp = data.frame(NBMolp, row.names = NBMolp)
NBMolp$symbol = as.character(lapply(rownames(NBMolp), function(x){
  x = strsplit(x, '_')[[1]][2]
  x = gsub('.', '-',x,fixed = T)
  x
}))

dat = QE_deldat_min0.8_zscore
dim(dat)# 6149*6364
dat[1:4,1:7]
dat$class = machine_info[rownames(dat), 'Classificaiton_type']
dat$TrainTest = machine_info[rownames(dat), 'TrainTest']
dat = subset(dat, class %in% c('N', 'B', 'M') &!is.na(Gender) & TrainTest=="Discovery")
dat = dat[order(dat['class']), ]
class = dat$class
dat = dat[6:6364]
# dat[] = lapply(dat, as.numeric)
dim(dat) # 6034 6359
dat = data.frame(t(dat))
dat$symbol = as.character(lapply(rownames(dat), function(x){
  x = strsplit(x, '_')[[1]][2]
  x = gsub('.', '-',x,fixed = T)
  x
}))
dim(dat)# 6359 3837


dat$symbol[grepl('NA$',dat$symbol,fixed = T)]

dat = subset(dat,  symbol !='NA' & !is.na(symbol))
rownames(dat) = dat$symbol

dat1 = aggregate(dat, by = list(dat$symbol), FUN = function(x){  mean(x, na.rm = T)} )
dim(dat )
dat1[1:3, 1:2]
# dat1$symbol = dat1$Group.1
rownames(dat1) = dat1$Group.1
dim(dat1 )
dat1[1:3, c(1,2,3837,3838)]
dat1 = dat1[2:3837]
dat1[1:3, c(1,2 )]


gsvaPar1 = gsva(expr = as.matrix(dat1[intersect(rownames(dat1), NBMolp$symbol),]) , 
                gset.idx.list = entriz,
                kcdf = "Gaussian" , # ("Gaussian", "Poisson", "none")
                method = "gsva",  # ("gsva", "ssgsea", "zscore", "plage")
                min.sz = 2,
                ssgsea.norm=T,
                verbose=F, 
                parallel.sz = parallel::detectCores())
gsvaPar1 = data.frame(gsvaPar1)
dim(gsvaPar1)
gsvaPar1[1:3, 1:3]

# write.xlsx(gsvaPar1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/QEunnormal_Discovery_NBMpadjust_gsva_C2-CP.xlsx', rowNames = T)
# write.xlsx(data.frame(t(gsvaPar1)), 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/QEunnormalzscore_Discovery_NBM_gsva_hallmark.xlsx', rowNames = T)
write.xlsx(data.frame(t(gsvaPar1)), 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/QEunnormalzscore_Discovery_NBMpadjust_gsva_hallmark.xlsx', rowNames = T)

# gsvaPar1_anova = anova_func(gsvaPar1, class)
# gsvaPar1_anova[c('Pvalue', 'FDR' )]
# 
# write.xlsx(gsvaPar1_anova[c('Pvalue', 'FDR')], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/QEunnormal_Discovery_NBMpadjust_gsva_C2-CP_anova.xlsx', rowNames = T)

# dim(subset(gsvaPar1_anova, Pvalue<0.05) )#& abs( logFC)>log2(1.5)
# dim(subset(gsvaPar1_anova, FDR<0.05))
# 
# gsvaPar1_ttest1 = subset(gsvaPar1_anova, FDR<0.05 )
# gsvaPar1_ttest1$pre = rowMeans( gsvaPar1_ttest1[sensitive$samp_ID] , na.rm = TRUE)
# gsvaPar1_ttest1$post = rowMeans( gsvaPar1_ttest1[nonRes$samp_ID] , na.rm = TRUE)
# gsvaPar1_ttest1$sub = gsvaPar1_ttest1$pre-gsvaPar1_ttest1$post
# dim(gsvaPar1_ttest1)
# gsvaPar1_ttest1 = gsvaPar1_ttest1[order(gsvaPar1_ttest1[,'sub']),]
# # gsvaPar1_ttest1= gsvaPar1_ttest1[c(1:15, 210:224),]
# gsvaPar1_ttest1= gsvaPar1_ttest1[c(1:15, 131:145), ]
# dim(gsvaPar1_ttest1)
# colnames(ptv3_Hsamp_sampinfo)
# 
# colnames(machine_info)


heatmapInfo = machine_info[colnames(gsvaPar1), c("Age", 'Classificaiton_type', 'Gender',"SampleType") ]
# rownames(heatmapInfo) = heatmapInfo$samp_ID
# heatmapInfo = heatmapInfo[colnames(gsvaPar1_ttest1)[1:192],]
colnames(heatmapInfo)
heatmapInfo$Age = as.numeric(heatmapInfo$Age)
heatmapInfo = heatmapInfo[order(heatmapInfo[,1]),]
heatmapInfo = heatmapInfo[order(heatmapInfo[,2]),]
dim(heatmapInfo)
heatmapInfo$Classificaiton_type <- factor(heatmapInfo$Classificaiton_type, levels = c('N', 'B', 'M'))
heatmapInfo <- heatmapInfo[order(heatmapInfo$Classificaiton_type), ]
dim(heatmapInfo)

# heatmapInfo = heatmapInfo[2:3]
# heatmapInfo  = subset(heatmapInfo , !is.na(heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO))
# heatmapInfo1 = data.frame(PRISM = heatmapInfo$PRISM1st_label_total, CCK8 = heatmapInfo$CCK8_cell_viability_before,
#                           PRISM1st_logFC = heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO,
#                           pert_time = paste0(heatmapInfo$pert_time,'h'),
#                           row.names = rownames(heatmapInfo))
# 
# rownames(gsvaPar1_ttest1) = tolower(rownames(gsvaPar1_ttest1))

rownames(gsvaPar1) = tolower(rownames(gsvaPar1))
# pert_time <- c("#007694", "#7DDFFF")
# names(pert_time) <- c(  "24h", "6h")
# ann_colors <- list(pert_time = pert_time, PRISM1st_logFC =c( "#366092","#C4D6E9") )
# ann_colors


#heatmapInfo1 = subset(heatmapInfo1, !is.na(heatmapInfo1$PRISM1st_logFC))
library(pheatmap)
pheatmap(gsvaPar1[rownames(heatmapInfo)], scale = 'row', #rownames(heatmapInfo)
         # color = colorRampPalette(c("blue", "white","red" ))(1000),
         color = colorRampPalette(c("#00B0F0", "white","yellow" ))(1000), 
         # breaks = c(seq(-2, 2, length=100)),
         annotation_col = heatmapInfo , # 样本信息
         # annotation_colors = ann_colors , 
         # labels_col = heatmapInfo$drugname, 
         # angle_col = 45,
         cluster_cols = F, 
         cluster_rows = T, 
         #clustering_method = 'median',
         cellwidth = 0.1,
         cellheight = 5,#2, # 
         fontsize_col = 1.5,
         fontsize_row = 4,
         show_rownames = T,
         show_colnames = F,
         name = ' ')

######### bar plot #######
# tempinfo = QE_deldat_IQR_min0.8[1:5]
# tempinfo$TrainTest = machine_info[rownames(tempinfo),  'TrainTest']
library(dplyr)
tempinfo = machine_info[c("Age", "SampleType", "DataSet", "Hospital", "Gender", 'TrainTest', "Thyroid_ID")]
dim(tempinfo)
tempinfo = tempinfo %>% distinct( Thyroid_ID, .keep_all = TRUE)

unique(tempinfo$Age)
temp = data.frame()
for(i in unique(tempinfo$Age)){
  if(is.na(i)){next}
  gendernum = table(tempinfo[tempinfo$Age == i, "Gender"])
  datanum = table(tempinfo[tempinfo$Age == i, "TrainTest"])
  temp1 = data.frame(gendernum)
  if(nrow(temp1)==0){next}
  temp1$Age = i
  temp = rbind(temp, temp1)
  # 
}
colnames(temp)
temp$Age = as.numeric(temp$Age)
ggplot(temp, aes(Age, Freq, fill = Var1))+
  geom_bar(stat = 'identity', position ="stack", width=0.5)+ # dodge
  scale_fill_d3(alpha = 0.7)+
  # geom_text(aes(label = Freq), position = position_stack(), vjust=-0.5, hjust = 0.5, size =2 )+
  # geom_text( aes(label = Age ), data = temp[temp$Var1 == 'Discovery',], #position = position_stack( ), 
  #            vjust= -0.5, hjust = 0.5, size =2 )+
  theme_classic()+ 
  theme(axis.text = element_text(size = 13),
        text = element_text(size = 13))
# + scale_y_continuous(expand = c(0,0) )
################### Test3 ####
delinfo = subset(machine_info, TrainTest=='Test3' & qe480 == 'QE' )
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
QETest3 = qe[rownames(delinfo)]
QETest3 = log2(QETest3)
QETest3_min0.8 = QETest3
QETest3_min0.8[is.na(QETest3_min0.8)] = min(QETest3_min0.8, na.rm = T) * 0.8
QETest3_min0.8 = data.frame(t(QETest3_min0.8))
write.xlsx(QETest3_min0.8, "D:/chh/2023workProject/20240821TTTD/QE/QETest3_log2min0.8.xlsx", rowNames = T)
############## ML result ####
# temp = data.frame(real = c(47, 39, 47, 31, 57, 68, 50, 42, 31, 31, 42, 47, 50, 65, 47, 65, 68, 54, 39, 64, 57, 64),
#                   Feat55 = c(44.96246337890625, 44.71418380737305, 44.160282135009766, 51.35565948486328, 51.61313247680664, 51.971988677978516, 42.35352325439453, 43.856536865234375, 36.79928207397461, 39.08518600463867, 44.24528884887695, 40.35720443725586, 50.20903396606445, 48.959312438964844, 41.01551818847656, 45.96669387817383, 51.624576568603516, 44.830726623535156, 43.33466720581055, 48.12007141113281, 57.341793060302734, 45.14720916748047 ),
#                   Feat60 = c(47.69786834716797, 44.9770393371582, 47.49500274658203, 54.25905227661133, 50.65324401855469, 51.255794525146484, 45.24101638793945, 40.7911491394043, 33.287227630615234, 34.96190643310547, 37.862525939941406, 38.64509963989258, 50.95775604248047, 46.749610900878906, 41.089237213134766, 43.447898864746094, 52.3726692199707, 45.29849624633789, 42.441471099853516, 50.60638427734375, 56.0914192199707, 49.86932373046875)
#                   )
# cor(temp, method = 'spearman')
# temp1 = melt(temp, id.vars = 'real')
# colnames(temp1)
# ggplot(temp1[temp1$variable == 'Feat55',], aes(real, value, color = variable))+
#   geom_smooth(aes(fill = variable), method = 'lm',alpha = 0.15)+ #  se = F
#   geom_point()+
#   labs(y = 'predict',x = 'real')+
#   theme_classic()+ 
#   annotate('text',x = 40, y= 55, label = 'Spearman r = 0.666' )+
#   theme(axis.text = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 15))
# Age_gap = temp$Feat55- temp$real
# temp1 = data.frame(Age_gap, Type = 'N')
# ggplot(temp1, aes(Type, Age_gap ))+
#   geom_boxplot()+
#   theme_classic()+ 
#   theme(axis.text = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 15))




temp2 = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/TTTD_ML_discoveryResult.xlsx')
colnames(temp2)
temp2$Age_gap = temp2$predictVal- temp2$realVal
temp2$Type = gsub('PTMC', 'PTC', temp2$Type)
unique(temp2$Type)
temp2$Type1 = temp2$Type
temp2[temp2$Type1 %in% c("FTC","PTC",'PDTC', "ATC", "MTC" ), 'Type1'] = 'M'

temp2$Type = factor(temp2$Type, levels = c( "B","Borderline","FTC","PTC",'PDTC', "ATC", "MTC" ))

unique(temp2$Type1)
temp2$Type1 = factor(temp2$Type1, levels = c( "B","Borderline","M" ))
ggplot(temp2, aes(Type1, Age_gap, color = Type1))+
  geom_boxplot()+
  scale_color_d3()+
  theme_classic()+ labs(title = 'Discovery')+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))

# 'Test1','Discovery'
# xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/models/xgb_predict.xlsx', rowNames = T)
# xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore/QE_xgbmae_imp_50feat_predict20241021.xlsx', rowNames = T)
# xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore/480_xgbmae_imp_50feat_predict20241021.xlsx', rowNames = T)
# xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QEtest3_xgbmae_imp_50feat_predict20241021.xlsx', rowNames = T)
# xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/480test3_xgbmae_imp_50feat_predict20241021.xlsx', rowNames = T)

xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QE_xgbmae_shap_95feat_predict20241031.xlsx', rowNames = T)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/480_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QEtest3_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)
xgb_result = read.xlsx('D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/480test3_xgbmae_shap_95feat_predict20241021.xlsx', rowNames = T)


xgb_result$class = machine_info[rownames(xgb_result), 'Classificaiton_type']
xgb_result$his = machine_info[rownames(xgb_result), 'Histopathology_type']
xgb_result$TrainTest = machine_info[rownames(xgb_result), 'TrainTest']
xgb_result$SampleType = machine_info[rownames(xgb_result), 'SampleType']
unique(xgb_result$his)
unique(xgb_result$class)
unique(xgb_result$SampleType)
xgb_result$his = gsub('PTMC', 'PTC', xgb_result$his)

xgb_result[is.na(xgb_result$TrainTest), 'TrainTest'] = 'Test2'
temp = subset(xgb_result, class == 'B' & !is.na(Pred))
temp = temp[c("Age", "Pred")]
#p = t.test(temp$Age, temp$Pred)
temp$Type = 'B'
temp$TrainTest = xgb_result[rownames(temp), 'TrainTest']
temp$SampleType = xgb_result[rownames(temp), 'SampleType']
for(i in c( 'N', 'M')){
  temp1 = subset(xgb_result, class == i & !is.na(Pred))
  if(nrow(temp1)==0){next}
  temp1 = temp1[c("Age", "Pred")]
  # p = t.test(temp1$Age, temp1$Pred)
  # print(c(i, p$p.value))
  temp1$Type = i
  temp1$TrainTest = xgb_result[rownames(temp1), 'TrainTest']
  temp1$SampleType = xgb_result[rownames(temp1), 'SampleType']
  temp = rbind(temp, temp1)
}

for(i in c( "PTC", "FTC", "PDTC" , "ATC", "MTC" )){
  temp1 = subset(xgb_result, his == i & !is.na(Pred))
  if(nrow(temp1)==0){next}
  temp1 = temp1[c("Age", "Pred")]
  # p = t.test(temp1$Age, temp1$Pred)
  # print(c(i, p$p.value))
  temp1$Type = i
  temp1$TrainTest = xgb_result[rownames(temp1), 'TrainTest']
  temp1$SampleType = xgb_result[rownames(temp1), 'SampleType']
  temp = rbind(temp, temp1)
}
temp$Age_gap = temp$Pred-temp$Age

# temp = subset(xgb_results, class =='B' &  TrainTest == 'Discovery')
# tempcor = temp[c('Age', 'Pred', 'TrainTest' )]

unique(temp$TrainTest)
tempcor = temp[temp$Type %in% c('N') & temp$TrainTest == 'Discovery',]# Test1 Discovery
pearsonr = cor(tempcor[1:2], method = 'pearson')
spearmanr = cor(tempcor[1:2], method = 'spearman')

ggplot(tempcor, aes(Age, Pred, color = Type))+
  geom_smooth(aes(fill = Type), method = 'lm',alpha = 0.15)+ #  se = F
  geom_point()+
  labs(y = 'predict',x = 'real')+
  labs(title = unique(tempcor$TrainTest))+
  theme_classic()+ 
  annotate('text',x = 40, y= 55, label = paste0('Pearson r=', round(pearsonr[2], 3) ))+  # 
  annotate('text',x = 40, y= 58, label = paste0('Spearman r=', round(spearmanr[2], 3) ) )+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))



temp = subset(xgb_results, class =='N' &  TrainTest == 'Test1')
tempcor = temp[c('Age', 'Pred', 'TrainTest' )]
pearsonr = cor(tempcor[1:2], method = 'pearson')
spearmanr = cor(tempcor[1:2], method = 'spearman')
ggplot(tempcor, aes(Age, Pred, color = TrainTest))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, size = 1)+
  labs(y = 'predict',x = 'real',title = unique(tempcor$TrainTest))+ xlim(5,90)+ylim(5, 90)+
  theme_classic()+ 
  annotate('text',x = 40, y= 55, label = paste0('Pearson r=', round(pearsonr[2], 3) ))+  # 
  annotate('text',x = 40, y= 58, label = paste0('Spearman r=', round(spearmanr[2], 3) ) )+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))


head(temp)
table(temp$Type)
unique(temp$TrainTest)
unique(xgb_result$SampleType)
# 'N', 'B', 'M' ,'FTC','PTC', 'PDTC','ATC','MTC'
# 'Discovery' 'Test2', 'Test1' 'Test3'
# "FFPE_punch"  "FNA_ex_vivo" "FFPE_slide"  "FF"
# temp1 = subset(temp, Type %in% c( 'N', 'B', 'FTC','PTC', 'PDTC','ATC','MTC') &  TrainTest == 'Test3' ) 
temp1 = subset(temp, Type %in% c( 'N', 'B', 'M') &  TrainTest == 'Discovery' )
dim(temp1)
table(temp1$Type)
unique(temp1$Type)
vec = combn(unique(temp1$Type), 2)
vecs = list()
for(i in seq(1,length(vec), 2)){
  # print((i+ 1) /2)
  vecs[[(i+ 1) /2]] = vec[i:(i+1)]
}
head(temp1)
IQR = quantile(temp1$Age)
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
#temp1 = subset(temp1, Age <= 66)#IQR[4]
temp1$Type = factor(temp1$Type, levels = c('N', 'B', 'FTC','PTC', 'PDTC','ATC','MTC'))
ggplot(temp1, aes(Type, Age_gap, color = Type))+
  scale_color_d3()+
  # scale_color_manual(values = c('#46AA6E', "#FF7F0E" , '#9C7822'))+  
  geom_boxplot(lwd = 1, outlier.color = NA)+
  geom_jitter( size=1, alpha = 0.1)+
  theme_classic()+ 
  labs(title = unique(temp1$TrainTest))+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) + 
  geom_signif(comparisons= vecs, #vecs 
          step_increase = 0.3,
          test="wilcox.test",  # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
          map_signif_level=T   # 标签样式F为数字，T为*号
        )

#### qe 480 med IQR####
dt = c("Discovery", "Test1", "Test2" , "Test3")
for(q4 in c('QE', '480')){
  for(i in dt){
    # i = "Discovery"
    temp = subset(xgb_results, TrainTest == i & !is.na(Pred) & qe480== q4)
    # temp$type = temp$class
    # temp2 = subset(temp,  his %in% c("PTC", "FTC", "PDTC" , "ATC", "MTC"))
    # temp2$type = temp2$his
    # temp = rbind(temp, temp2)
    temp = subset(temp,  class %in% c('N', 'B', 'M',"PTC", "FTC", "PDTC" , "ATC", "MTC"))
    temp$Age_gap = (temp$Age_gap-mean(temp$Age_gap, na.rm = T))/sd(temp$Age_gap, na.rm = T)
    for(type1 in c('N', 'B', 'M')){
      # type1 = 'N'
      temp1 = subset(temp, class == type1 & !is.na(Pred))
      if( nrow(temp1)<2){
        #print(c(qe480, nrow(temp1), i, type1 ))
        next
      }
      # temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
      med = median(temp1$Age_gap)
      IQR = quantile(temp1$Age_gap)
      IQR = IQR[[4]]-IQR[[2]]
      print(c(q4, nrow(temp1), i, type1, med, IQR))
    }
    
    temp = subset(xgb_results, TrainTest == i & !is.na(Pred) & qe480== q4)
    temp = subset(temp,  his %in% c('N', 'B', 'M',"PTC", "FTC", "PDTC" , "ATC", "MTC"))
    temp$Age_gap = (temp$Age_gap-mean(temp$Age_gap, na.rm = T))/sd(temp$Age_gap, na.rm = T)
    for(type1 in c("PTC", "FTC", "PDTC" , "ATC", "MTC")){
      temp1 = subset(temp, his == type1 & !is.na(Pred))
      if( nrow(temp1)<2){
        #print(c(q4, nrow(temp1), i, type1 ))
        next
      }
      # temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
      med = median(temp1$Age_gap)
      IQR = quantile(temp1$Age_gap)
      IQR = IQR[[4]]-IQR[[2]]
      print(c(q4, nrow(temp1), i, type1, med, IQR))
    }
    
    temp = subset(xgb_results, TrainTest == i & !is.na(Pred) & qe480== q4)
    temp = subset(temp,  Bethesda_Score %in% c("II",  "I",   "IV",  "VI",  "III", "V"))
    temp$Age_gap = (temp$Age_gap-mean(temp$Age_gap, na.rm = T))/sd(temp$Age_gap, na.rm = T)
    for(type1 in c("II",  "I",   "IV",  "VI",  "III", "V")){
      temp1 = subset(temp, Bethesda_Score == type1 & !is.na(Pred))
      if( nrow(temp1)<2){
        #print(c(nrow(temp1), i, type1 ))
        next
      }
      # temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
      med = median(temp1$Age_gap)
      IQR = quantile(temp1$Age_gap)
      IQR = IQR[[4]]-IQR[[2]]
      print(c(qe480, nrow(temp1), i, type1, med, IQR))
    }
  }
}


dt = c("Discovery", "Test1", "Test2" , "Test3")
#"QE"  "480"
q4 = '480'
i = "Test3"
temp = subset(xgb_results, TrainTest == i & !is.na(Pred) & xgb_results$qe480== q4 )# 
temp1 = subset(temp,  class %in% c('N', 'B', 'M'))
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
unique(temp1$class)
p1 = wilcox.test(temp1[temp1$class == 'N', 'Age_gap'], temp1[temp1$class == 'B', 'Age_gap'] )
p2 = wilcox.test(temp1[temp1$class == 'N', 'Age_gap'], temp1[temp1$class == 'M', 'Age_gap'] )
p3 = wilcox.test(temp1[temp1$class == 'M', 'Age_gap'], temp1[temp1$class == 'B', 'Age_gap'] )
print(c(p1$p.value, p2$p.value, p3$p.value)) #

temp1 = subset(temp,  his %in% c("PTC", "FTC", "PDTC" , "ATC", "MTC"))
unique(temp1$his)
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
p1 = aov(as.numeric(temp1$Age_gap) ~ temp1$his)
summary(p1)[[1]][,5][1]

temp = subset(xgb_results, TrainTest == i & !is.na(Pred) & xgb_results$qe480== q4)# 
temp1 = subset(temp, Bethesda_Score %in% c("II",  "I",   "IV",  "VI",  "III", "V"))
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)
p1 = aov(as.numeric(temp1$Age_gap) ~ temp1$Bethesda_Score)
summary(p1)[[1]][,5][1]





temp = subset(xgb_results, class == 'B' & !is.na(Pred))
temp = temp[c("Age", "Pred", 'TrainTest', 'SampleType', 'samples', 'Bethesda_Score')]
temp$Type = 'B'

for(i in c( 'N', 'M')){
  temp1 = subset(xgb_results, class == i & !is.na(Pred))
  if(nrow(temp1)==0){next}
  temp1 = temp1[c("Age", "Pred", 'TrainTest', 'SampleType', 'samples', 'Bethesda_Score')]
  # p = t.test(temp1$Age, temp1$Pred)
  # print(c(i, p$p.value))
  temp1$Type = i
  temp = rbind(temp, temp1)
}

for(i in c( "PTC", "FTC", "PDTC" , "ATC", "MTC" )){
  temp1 = subset(xgb_results, his == i & !is.na(Pred))
  if(nrow(temp1)==0){next}
  temp1 = temp1[c("Age", "Pred", 'TrainTest', 'SampleType', 'samples', 'Bethesda_Score')]
  # p = t.test(temp1$Age, temp1$Pred)
  # print(c(i, p$p.value))
  temp1$Type = i
  temp = rbind(temp, temp1)
}
temp$qe480 = machine_info[temp$samples, 'qe480']
temp$Age_gap = temp$Pred-temp$Age
unique(temp$TrainTest)
# "Discovery" "Test1"     "Test2"
# 'N', 'B', 'FTC','PTC', 'PDTC','ATC','MTC'
# qe_FNA_NBM_Discovery_zscore
# qe_FNA_Discovery_zscore
temp1 = subset(temp , qe480 == '480' & TrainTest == "Test2" & Type %in% c('N', 'B', 'FTC','PTC', 'PDTC','ATC','MTC') & SampleType %in% c("FNA_ex_vivo", "FNA_in_vivo"))
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T)

vec = combn(unique(temp1$Type), 2)
vecs = list()
for(i in seq(1,length(vec), 2)){
  # print((i+ 1) /2)
  vecs[[(i+ 1) /2]] = vec[i:(i+1)]
}
vecs
temp1$Type = factor(temp1$Type, levels = c('N', 'B', 'FTC','PTC', 'PDTC','ATC','MTC'))

ggplot(temp1, aes(Type, Age_gap, color = Type))+
  scale_color_d3()+
  # scale_color_manual(values = c( "#FF7F0E" , '#9C7822'))+ # '#46AA6E',
  geom_boxplot(lwd = 1, outlier.color = NA)+
  geom_jitter( size=1, alpha = 0.1)+
  theme_classic()+ 
  labs(title = unique(temp1$TrainTest))+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15)) + 
  geom_signif(comparisons= vecs, #vecs 
              step_increase = 0.3,
              test="wilcox.test",  # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
              map_signif_level=T   # 标签样式F为数字，T为*号
  )


temp1 = temp1[temp1$Type == 'N',]
unique(temp1$Type)
quantile(temp1$Age_gap)

#lm
for(typenow in unique(temp$Type)){
  templm = subset(temp, !is.na(Age_gap) & Type==typenow & TrainTest == 'Discovery')
  templm$Type = as.numeric(factor(templm$Type))
  fit = glm(Type ~ Age_gap, data = templm, family  =  gaussian)#lm(Type ~ Age_gap, data = templm )
#fit
  sumfit = summary(fit)
  print(c(typenow, sumfit$coefficients[2,4]) )
}




ML95pros = c('Q8NI22_MCFD2', 'Q96A11_GAL3ST3', 'P60033_CD81', 'P26440_IVD', 'P02743_APCS', 'O14639_ABLIM1', 'Q9P1F3_ABRACL', 'O94811_TPPP', 'Q9UHL4_DPP7', 'Q8NBI5_SLC43A3', 'P02794_FTH1', 'P30038_ALDH4A1', 'Q86WU2_LDHD', 'P01011_SERPINA3', 'O75503_CLN5', 'P06727_APOA4', 'Q9UNL2_SSR3', 'O60888_CUTA', 'P0DJI8_SAA1', 'Q06136_KDSR', 'Q9Y624_F11R', 'P04004_VTN', 'Q92743_HTRA1', 'P06703_S100A6', 'P01023_A2M', 'Q8TAE6_PPP1R14C', 'P10253_GAA', 'Q06481_APLP2', 'P12830_CDH1', 'P08603_CFH', 'Q16134_ETFDH', 'O43684_BUB3', 'P37840_SNCA', 'P11182_DBT', 'P40121_CAPG', 'Q96DG6_CMBL', 'P01903_HLA.DRA', 'P59044_NLRP6', 'P61626_LYZ', 'Q9H910_JPT2', 'Q96SI9_STRBP', 'Q96RS6_NUDCD1', 'P51858_HDGF', 'Q15276_RABEP1', 'Q14746_COG2', 'P52655_GTF2A1', 'P00488_F13A1', 'Q9P0M6_MACROH2A2', 'O95810_CAVIN2', 'P04439_HLA.A', 'P10768_ESD', 'P21399_ACO1', 'O75955_FLOT1', 'P37802_TAGLN2', 'P62979_RPS27A', 'P52907_CAPZA1', 'Q08380_LGALS3BP', 'Q12841_FSTL1', 'Q96EE3_SEH1L', 'P00325_ADH1B', 'P02787_TF', 'P83876_TXNL4A', 'P20700_LMNB1', 'Q12906_ILF3', 'Q8IW45_NAXD', 'P80723_BASP1', 'P43251_BTD', 'Q92777_SYN2', 'Q9NUP9_LIN7C', 'P50226_SULT1A2', 'P22748_CA4', 'P07602_PSAP', 'O75173_ADAMTS4', 'P02671_FGA', 'Q8IYI6_EXOC8', 'P43243_MATR3', 'O43521_BCL2L11', 'P30711_GSTT1', 'P61956_SUMO2', 'Q9UKK3_PARP4', 'Q96FZ7_CHMP6', 'Q16822_PCK2', 'P02790_HPX', 'Q9BXN1_ASPN', 'P06576_ATP5F1B', 'Q92522_H1.10', 'P05155_SERPING1', 'P15088_CPA3', 'P05388_RPLP0', 'P07942_LAMB1', 'Q96I24_FUBP3', 'Q99523_SORT1', 'Q9Y3B2_EXOSC1', 'O15439_ABCC4', 'P53582_METAP1', 'P01876_IGHA1', 'P54289_CACNA2D1', 'Q5HYI8_RABL3', 'Q8NE62_CHDH', 'Q00796_SORD')
ML95pros = ML95pros[1:95]

library(pheatmap)
QE_deldat_IQR_min0.8[1:3,1:6]

dat = QE_deldat_IQR_min0.8_zscore[ML95pros]
dat = dat480_deldat_min08[1:3]
# dat = QE_deldat_IQR_min0.8[ML95pros]
temp = intersect(rownames(dat), rownames(xgb_results) )
sampleinfo = xgb_results[temp, c(1,12,2,4,5,9,10,13,14)]
dat = dat[temp, ]
colnames(sampleinfo)
dim(sampleinfo)

sampleinfo = sampleinfo[order(sampleinfo[,1]),]
sampleinfo = sampleinfo[order(sampleinfo[,2]),]


mycol = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#9467bd', '#8c564b', '#e377c2', '#8f7f7f',
          '#bcbd22', '#17becf',"#708090",'#68A180',
          '#F3B1A0', '#D6E7A3',"#99CCCC","#66CC99",
          "yellow", "blue", "gray")

SampleType <- mycol[1:length(unique(sampleinfo$SampleType))]
names(SampleType) <- unique(sampleinfo$SampleType)
class <- mycol[1:length(unique(sampleinfo$class))]
names(class) <- unique(sampleinfo$class)
Hospital <- mycol[1:length(unique(sampleinfo$Hospital))]
names(Hospital) <- unique(sampleinfo$Hospital)
his <- mycol[1:length(unique(sampleinfo$his ))]
names(his ) <- unique(sampleinfo$his )
Gender <- mycol[1:length(unique(sampleinfo$Gender))]
names(Gender) <- unique(sampleinfo$Gender)
Bethesda_Score <- mycol[1:length(unique(sampleinfo$Bethesda_Score))]
names(Bethesda_Score) <- unique(sampleinfo$Bethesda_Score)
BRAF <- mycol[1:length(unique(sampleinfo$BRAF))]
names(BRAF) <- unique(sampleinfo$BRAF)

ann_colors <- list( class = class, SampleType = SampleType, Hospital = Hospital,
                    his  = his , Gender = Gender, Bethesda_Score = Bethesda_Score, BRAF=BRAF )
ann_colors

pheatmap(t(dat[rownames(sampleinfo),]), # scale="row",#对行进行归一化
         color = colorRampPalette(c("blue", "white","red" ))(1000), # color参数自定义颜色
         breaks = c(seq(-6, 6, length=1000)),
         annotation_col = sampleinfo, # 样本信息
         annotation_colors = ann_colors ,
         fontsize_col = 0.55, 
         fontsize_row = 1.5,
         cluster_rows = T,
         cluster_cols = F,
         #labels_col = as.character(sampleinfo$Age),
         show_rownames = F,
         show_colnames = F,
         fontsize = 7, 
         cellwidth= 0.1, 
         cellheight= 2)

######## cibersort #######
# library('devtools')
# devtools::install_github("Moonerss/CIBERSORT")
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
library(dplyr) # function： %>%
library (tidyr) # function： gather
library(tibble) # function： rownames_to_column
library(ggplot2)
library(CIBERSORT)
sig_matrix = system.file("extdata", "LM22.txt", package = "CIBERSORT")

QE_deldat = qe[,rownames(delinfo)]
na_col = colSums(!is.na(QE_deldat))
dim(QE_deldat) # 8509 6448
QE_deldat = QE_deldat[na_col>=2000 ] # 299
dim(QE_deldat)# 8509 6149
na_row = rowSums(is.na(QE_deldat))/ncol(QE_deldat)
length(na_row)
na_row[na_row==1 ]
QE_deldat = QE_deldat[na_row<=0.9, ]
Q = quantile(QE_deldat, na.rm = T)
Q3 = Q[4]
Q1 = Q[2]
IQR = Q3-Q1
IQR
up = Q3+2*IQR
down = Q1-2*IQR
max(QE_deldat, na.rm = T)
plot(density(na.omit(unlist(QE_deldat))), main="density")

QE_deldat[QE_deldat>up] = up
QE_deldat[1:4,1:2]
dim(QE_deldat)
QE_deldat = QE_deldat[!grepl('_NA$', rownames(QE_deldat)),]
dim(QE_deldat)
Symbols = lapply(rownames(QE_deldat), function(x){strsplit(x, '_')[[1]][2]})
head(Symbols)
QE_deldat[is.na(QE_deldat)] = min(QE_deldat, na.rm = T) * 0.8
QE_deldat = cbind(as.character(Symbols), QE_deldat)
colnames(QE_deldat)[1] = 'Symbols'
QE_deldat = aggregate(QE_deldat, by = list(QE_deldat$Symbols), FUN = 'mean' )
QE_deldat[1:4,1:2]
write.table(QE_deldat, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE.txt", sep = "\t", row.names = F)


N_QE_deldat = QE_deldat[c('Group.1', intersect(colnames(QE_deldat), rownames(delinfo[delinfo$Classificaiton_type == 'N', ])) )]
B_QE_deldat = QE_deldat[c('Group.1', intersect(colnames(QE_deldat), rownames(delinfo[delinfo$Classificaiton_type == 'B', ])) )]
M_QE_deldat = QE_deldat[c('Group.1', intersect(colnames(QE_deldat), rownames(delinfo[delinfo$Classificaiton_type == 'M', ])) )]
write.table(N_QE_deldat, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_N.txt", sep = "\t", row.names = F)
write.table(B_QE_deldat, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_B.txt", sep = "\t", row.names = F)
write.table(M_QE_deldat, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_M.txt", sep = "\t", row.names = F)
N_QE_deldat[1:3,1:3]


N_QE_deldat = read.table( "D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_N.txt", sep = "\t", header = T, row.names = 1)
B_QE_deldat = read.table("D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_B.txt", sep = "\t", header = T, row.names = 1)
M_QE_deldat = read.table("D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_M.txt", sep = "\t", header = T, row.names = 1)
NBM_mean = data.frame(symbol = rownames(N_QE_deldat), Nmean = rowMeans(N_QE_deldat), Bmean = rowMeans(B_QE_deldat), Mmean = rowMeans(M_QE_deldat))
write.table(NBM_mean, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_NBMmean.txt", sep = "\t", row.names = F)


type1 = 'B'
results <- cibersort(sig_matrix, paste0("D:/chh/2023workProject/20240821TTTD/QE/cibersort_QE_",type1,".txt"), perm=1000, QN=F)
# Ncells = results
Bcells = results
# Mcells = results


results <- cibersort(sig_matrix, "D:/chh/2023workProject/20240821TTTD/QE/cibersort_NBMmean.txt", perm=1000, QN=F)
dim(results)
colnames(results)
results[1:10, 1:2]
results[1:10, 21:23]

results = Mcells

results_depso <- results[,-( (ncol(results)-2):ncol(results))]
info = data.frame(Sample = rownames(results_depso), Age = machine_info[rownames(results_depso), 'Age'])
info = info[order(info[, 2]), ]

dat <- results_depso %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)
head(dat)
#x <- factor(as.integer(rownames(dat)),labels= dat$Sample)

plotOrder <- c(colnames(results_depso))
dat$Cell_type <- factor(dat$Cell_type,levels = plotOrder)
#dat$Sample = factor(dat$Sample, levels = info$Sample)

ggplot(data = dat, aes(x, Proportion,fill = Cell_type))+ 
  geom_bar(stat ="identity",position = 'stack')+
  labs(fill = "Cell", x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme( title = element_text(size = 13),
         plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        
  ) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22)) +
  theme(axis.text.x = element_text(
    #angle = 90,
    # hjust =1,
    # vjust =0.5
    
  ))

dat$type = 'M'
dat1 = dat

dat1 = rbind(dat, dat1)
head(dat1)

x <- factor(as.integer(rownames(dat1)), labels= dat1$Sample)
dat1$type = factor(dat1$type, levels = c('N', 'B', 'M'))
ggplot(data = dat1, aes(Cell_type, Proportion, fill = type))+ 
  geom_boxplot()+
  scale_fill_d3()+
  # scale_fill_manual(values = mypalette(22)) +
  theme_classic() +
  labs(fill = "Cell", x = "",y = "Estiamted Proportion") +
  theme( title = element_text(size = 13),
         plot.title = element_text(hjust = 0.5),
         #axis.text.x = element_blank(),
         axis.text.x = element_text( angle = 60, hjust = 1
                                     ),
         #axis.ticks.x = element_blank(),
         legend.position = "bottom",
         text = element_text(size = 13),
         legend.text = element_text(size = 12),
         axis.text.y = element_text(size = 12),
         axis.title.x = element_text(size = 12),
         axis.title.y = element_text(size = 12),
         legend.title = element_text(size = 12),
         
  )# +scale_y_continuous(expand = c(0.01,0)) +theme(axis.text.x = element_text( ))
  
  
  
############ lasso ####
# https://zhuanlan.zhihu.com/p/159307747
# https://cloud.tencent.com/developer/article/1975780
# https://ayueme.github.io/R_clinical_model/feature-selection_lasso.html
library(glmnet)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggthemes)
coefficients1 = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/QE_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
padjust = rownames(coefficients1[coefficients1$qt_padjust<0.05,] )
delinfo = delinfo[rownames(QE_deldat_IQR_min0.8_zscore),]
unique(delinfo$TrainTest)
NB = subset(delinfo, Classificaiton_type %in% c('N', 'B') &  TrainTest == 'Discovery')
NBtest = subset(delinfo, Classificaiton_type %in% c('N', 'B') &  TrainTest == 'Test1')

NB = subset(delinfo, Classificaiton_type %in% c('N' ) &  TrainTest == 'Discovery')
NBtest = subset(delinfo, Classificaiton_type %in% c('N' ) &  TrainTest == 'Test1')

lassoX = as.matrix(QE_deldat_IQR_min0.8_zscore[rownames(NB),padjust])
y = QE_deldat_IQR_min0.8_zscore[rownames(NB), 'Age']

alpha_seq <- 10^seq(-4, 2, length.out = 10)
alpha_seq = c(alpha_seq,  5, 10, 50, 100 )
# 1. 每个alpha值进行一次交叉验证
# 返回结果：
# cvm：就是这10次交叉验证的错误度量平均值，常规线性模型默认使用Deviance，也就是MSE（平均标准误差）,logistics回归是使用Bionomical Deviance
# cvsd：10次交叉验证的错误度量标准差
# lambda： 尝试的lambda值
# index_min：最低错误率模型对应的lambda值
# index_1se：错误率可接受的最精简模型对应的lambda值
cv_model_ls <-
  alpha_seq %>%
  set_names(., .) %>% # 对向量添加names为自身，保证map返回的列表也是有names的
  map(function(alpha){
    cv.model <- cv.glmnet(lassoX, y, type.measure = "mae", family="gaussian", nfolds = 3, alpha = alpha)
    # alpha, min_lambda, 1se_lambda, error
    print(c(alpha, cv.model$lambda.min, cv.model$cvm[cv.model$index["1se",1]])) #  cv.model$lambda.1se, 
    list(
      cvm = cv.model$cvm,
      cvsd = cv.model$cvsd,
      lambda = cv.model$lambda,
      index_min = cv.model$index["min",1],
      index_1se = cv.model$index["1se",1]
    )
  }
  )

# 2. 手动绘制11个alpha值下的lambda分布曲线，见下图
cv_model_ls  %>%
  map2(names(.), function(model, alpha){
    labmda_min <- model$lambda[model$index_min] %>% log
    labmda_1se <- model$lambda[model$index_1se] %>% log
    
    data.frame(x = log(model$lambda),
               y = model$cvm,
               std = model$cvsd) %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(color = "Red") +
      geom_errorbar(aes(ymin = y - std, ymax = y + std), color = "gray") +
      geom_vline(xintercept = c(labmda_min, labmda_1se) , linetype = 3) + 
      labs(x = expression(paste("Log(", lambda, ")")), y = "MAE") +
      theme_bw() +
      ggtitle(paste0("alpha: ", alpha ))
  }) %>%
  wrap_plots(
    ncol = 4
  )

# 选择alpha, CV建模筛选特征,  CV训练
cv_lassofit <- cv.glmnet(lassoX, y, type.measure="mae", alpha= 1, family="gaussian", nfolds = 3 )
summary(cv_lassofit)
plot(cv_lassofit)
# min代表的是在所有的lambda值中，是mse最小的那一个值,1se是指在min一个方差范围内得到最简单模型的那一个lambda值
# 1se给出的是一个具备优良性能且自变量个数最少的模型
print(cv_lassofit)
cv_lassofit$lambda
cv_lassofit$lambda.min
cv_lassofit$lambda.1se
lassocoef = data.frame(coef(cv_lassofit, s=cv_lassofit$lambda.min)) # lambda.min
lassocoef$pros = rownames(lassocoef)
lassocoef = lassocoef[lassocoef$s1!=0, ]
lassocoef



trainDat = as.matrix(QE_deldat_IQR_min0.8_zscore[rownames(NB), rownames(lassocoef)[2:nrow(lassocoef)]])
dim(trainDat)
testdat = QE_deldat_IQR_min0.8_zscore[rownames(NBtest), rownames(lassocoef)[2:nrow(lassocoef)]] # 
testdat[1:3,1:3]
# lassofit = glmnet(trainDat, y, type.measure="mae", alpha=0.8, family="gaussian" )
# plot(lassofit, xvar="lambda", label=TRUE)
cv_lassofit1 <- cv.glmnet(trainDat, y, type.measure="mae", alpha=1, family="gaussian", nfolds = 3 )
res = predict(cv_lassofit1, newx = as.matrix(testdat), s = "lambda.min", type = "response") # s = "lambda.1se"
res = data.frame(res)
res$real = QE_deldat_IQR_min0.8_zscore[rownames(NBtest), 'Age']
res$Type = delinfo[rownames(NBtest) ,'Classificaiton_type']
res$Age_gap= res$lambda.min-res$real
cor(res$lambda.min, res$real)

library(ggsci)
ggplot(res, aes(Type , Age_gap, color = Type ))+
  geom_boxplot()+
  scale_color_d3()+
  theme_classic()+ labs(title = 'Discovery N&B')+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))
