library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggsignif)
library(reshape2)
library(dplyr)
library(pheatmap)
### info
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
# date match
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
####### scale method ########
stdscale = function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}
#### dat480 
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


delinfo = subset(machine_info, Tissue_type=='thyroid' & Histopathology_type!='uncertain' & SampleType!='FNA_in_vivo')
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
delinfo = subset(delinfo, Histopathology_type!='MEC' )
delinfo = subset(delinfo, qe480=='480' )
unique(delinfo$TrainTest)
dat480_deldat = dat480[,rownames(delinfo)]
na_col = colSums(!is.na(dat480_deldat))
dim(dat480_deldat) # 8509 6448
dat480_deldat = dat480_deldat[na_col>=2000 ] # 299
dim(dat480_deldat)# 8509 6149
na_row = rowSums(is.na(dat480_deldat))/ncol(dat480_deldat)
length(na_row)
na_row[na_row==1 ]
dat480_deldat = dat480_deldat[na_row<=0.9, ]


dat480_deldat = log2(dat480_deldat)
dat480_deldat[1:2,1:2]
delinfo = delinfo[colnames(dat480_deldat), ]# important
unique(delinfo$qe480)
# 480 IQR
Q = quantile(dat480_deldat, na.rm = T)
Q3 = Q[4]
Q1 = Q[2]
IQR = Q3-Q1
IQR
up = Q3+2*IQR
down = Q1-2*IQR
dat480_deldat_IQR = dat480_deldat
dat480_deldat_IQR[dat480_deldat_IQR>up] = up
dat480_deldat_IQR[dat480_deldat_IQR<down] = down
temp = NA
plot(density(na.omit(unlist(dat480_deldat_IQR))), main="2IQR density default")
sum(is.na(dat480_deldat_IQR))/(ncol(dat480_deldat_IQR) * nrow(dat480_deldat_IQR))

dat480_deldat_IQR_min0.8 = dat480_deldat_IQR
dat480_deldat_IQR_min0.8[is.na(dat480_deldat_IQR_min0.8)] = min(dat480_deldat_IQR_min0.8, na.rm = T) * 0.8
dat480_deldat_IQR_min0.8[1:20,1:2]
dat480_deldat_IQR_min0.8 = data.frame(t(dat480_deldat_IQR_min0.8))
dim(dat480_deldat_IQR_min0.8)
plot(density(na.omit(unlist(dat480_deldat_IQR_min0.8 ))), main="density default")

dim(dat480_deldat_IQR_min0.8)# 2975sample 7498pro
dat480_deldat_IQR_min0.8$Age = as.numeric(machine_info[rownames(dat480_deldat_IQR_min0.8), 'Age'])
dat480_deldat_IQR_min0.8$SampleType = machine_info[rownames(dat480_deldat_IQR_min0.8), 'SampleType']
dat480_deldat_IQR_min0.8$DataSet = machine_info[rownames(dat480_deldat_IQR_min0.8), 'DataSet']
dat480_deldat_IQR_min0.8$Hospital = machine_info[rownames(dat480_deldat_IQR_min0.8), 'Hospital']
dat480_deldat_IQR_min0.8$Gender = machine_info[rownames(dat480_deldat_IQR_min0.8), 'Gender']
dim(dat480_deldat_IQR_min0.8)
dat480_deldat_IQR_min0.8 = cbind(dat480_deldat_IQR_min0.8[,7499:7503], dat480_deldat_IQR_min0.8[,1:7498])
dat480_deldat_IQR_min0.8[1:3,1:6]
write.xlsx(dat480_deldat_IQR_min0.8, 'D:/chh/2023workProject/20240821TTTD/dat480/thyroid_480_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
# 480 zscore
dat480_deldat_IQR_min0.8_zscore = dat480_deldat_IQR_min0.8
dat480_deldat_IQR_min0.8_zscore[,6:ncol(dat480_deldat_IQR_min0.8_zscore)] = apply(dat480_deldat_IQR_min0.8_zscore[,6:ncol(dat480_deldat_IQR_min0.8_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})

##### dat480 Test3
delinfo = subset(machine_info, TrainTest=='Test3' & qe480 == '480' )
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
dat480[1:10,1:3]
dat480Test3 = dat480[rownames(delinfo)]
dat480Test3 = log2(dat480Test3)
dat480Test3_min0.8 = dat480Test3
dat480Test3_min0.8[is.na(dat480Test3_min0.8)] = min(dat480Test3_min0.8, na.rm = T) * 0.8
dat480Test3_min0.8 = data.frame(t(dat480Test3_min0.8))
write.xlsx(dat480Test3_min0.8, "D:/chh/2023workProject/20240821TTTD/dat480/dat480Test3_log2min0.8.xlsx", rowNames = T)
##### qe
qe = read.csv("D:/chh/2023workProject/20240821TTTD/datas/TTTD_QE_all7038_20230623.pg_matrix.tsv", sep = '\t')
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

QE_deldat = log2(QE_deldat)
QE_deldat[1:10,1:2]
delinfo = delinfo[colnames(QE_deldat), ]# important
unique(delinfo$qe480)
# qe IQR
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
# qe zscore
QE_deldat_IQR_min0.8_zscore = QE_deldat_IQR_min0.8
QE_deldat_IQR_min0.8_zscore[,6:ncol(QE_deldat_IQR_min0.8_zscore)] = apply(QE_deldat_IQR_min0.8_zscore[,6:ncol(QE_deldat_IQR_min0.8_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})
# fig6 fig1 age and Gender barplot
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
}
colnames(temp)
temp$Age = as.numeric(temp$Age)
ggplot(temp, aes(Age, Freq, fill = Var1))+
  geom_bar(stat = 'identity', position ="stack", width=0.5)+ # dodge
  scale_fill_d3(alpha = 0.7)+
  theme_classic()+ 
  theme(axis.text = element_text(size = 13),
        text = element_text(size = 13))

# qe tsne
tsneinput = QE_deldat_IQR_min0.8[rownames(QE_deldat_IQR_min0.8) %in% rownames(machine_info[machine_info$TrainTest %in% c("Discovery", "Test1") & machine_info$Classificaiton_type == 'N',]), 
                                 6:ncol(QE_deldat_IQR_min0.8)]
dim(tsneinput)
set.seed(2024)
dfQE.tsne <- Rtsne(tsneinput, dims = 2, perplexity = 10,
                   partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
dim(dfQE.tsne$Y)
dim(machine_info )
# 'machine', "DataSet", "SampleType", "Histopathology_type" "Gender", "Age_grade", "Tissue_type", "Hospital"
label = "TrainTest" 
df1 = data.frame(dfQE.tsne$Y, machine_info[rownames(tsneinput), label] ) 
names(df1) = c("tSNE1", "tSNE2", "Batch")
dim(df1)
head(df1) 
# write.xlsx(df1, "tsne_filterInfoNA_qe.xlsx", rowNames = T)
ggplot(df1, aes(tSNE1,tSNE2, colour= Batch))+
  geom_point(aes(color= Batch), size=2)+
  scale_color_d3(palette = c( "category20" ))+# scale_color_jama()
  # scale_color_gradient(low='cyan', high = 'red')+ # date
  theme_bw() + labs(title = label)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13,color = "black"),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.25),
        axis.line.y = element_line(color="black", size = 0.25),
        plot.title  = element_text(size=13, hjust = 0.5),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        axis.title  = element_text(size=13 ),
        axis.title.x = element_text( hjust = 0.5 ),
        axis.title.y  = element_text( hjust = 0.5 ),
        panel.background = element_blank(),
        #legend.position = c(0.9,0.1), # legend的位置信息
        #legend.background = element_rect(fill = "grey90", size = 1, colour = "white")
  )
####### qe N DESWAN########
source('D:/chh/2023workProject/20240821TTTD/DEswan-master/R/DEswan.R')
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/reshape.DEswan.R")
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/nsignif.DEswan.R")
source("D:/chh/2023workProject/20240821TTTD/DEswan-master/R/q.DEswan.R")

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
### DEswan 1 window 
# linear model
type1 = 'N'
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

# non-linear model
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

for(type1 in c('N', 'M','B', 'Borderline')){ #
  print(type1)
  # type1 = 'NM'
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

files = list.files('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/slide_1/')
files = files[grepl('_p.xlsx',files)]
files
for(n in files){
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
####### volcano #####
# volcano fig2 B
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
### fig4 bubble  
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
####### ML boxplot #####
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
temp1$Age_gap = (temp1$Age_gap-mean(temp1$Age_gap, na.rm = T))/sd(temp1$Age_gap, na.rm = T) # zscore

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

########### 甲状腺相关蛋白重点 #########
# figs3 log2 std Expressions
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
         
  )
