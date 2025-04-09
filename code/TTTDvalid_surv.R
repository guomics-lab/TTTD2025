library(openxlsx)
library(survminer)
library(survival)
# library(glmnet)
# library(tidyverse)
# library(timeROC)
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(2024)
setwd("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid")
top_features_ = c('Q8NI22_MCFD2', 'Q96A11_GAL3ST3', 'P60033_CD81', 'P26440_IVD', 'P02743_APCS', 'O14639_ABLIM1', 'Q9P1F3_ABRACL', 'O94811_TPPP', 'Q9UHL4_DPP7', 'Q8NBI5_SLC43A3', 'P02794_FTH1', 'P30038_ALDH4A1', 'Q86WU2_LDHD', 'P01011_SERPINA3', 'O75503_CLN5', 'P06727_APOA4', 'Q9UNL2_SSR3', 'O60888_CUTA', 'P0DJI8_SAA1', 'Q06136_KDSR', 'Q9Y624_F11R', 'P04004_VTN', 'Q92743_HTRA1', 'P06703_S100A6', 'P01023_A2M', 'Q8TAE6_PPP1R14C', 'P10253_GAA', 'Q06481_APLP2', 'P12830_CDH1', 'P08603_CFH', 'Q16134_ETFDH', 'O43684_BUB3', 'P37840_SNCA', 'P11182_DBT', 'P40121_CAPG', 'Q96DG6_CMBL', 'P01903_HLA.DRA', 'P59044_NLRP6', 'P61626_LYZ', 'Q9H910_JPT2', 'Q96SI9_STRBP', 'Q96RS6_NUDCD1', 'P51858_HDGF', 'Q15276_RABEP1', 'Q14746_COG2', 'P52655_GTF2A1', 'P00488_F13A1', 'Q9P0M6_MACROH2A2', 'O95810_CAVIN2', 'P04439_HLA.A', 'P10768_ESD', 'P21399_ACO1', 'O75955_FLOT1', 'P37802_TAGLN2', 'P62979_RPS27A', 'P52907_CAPZA1', 'Q08380_LGALS3BP', 'Q12841_FSTL1', 'Q96EE3_SEH1L', 'P00325_ADH1B', 'P02787_TF', 'P83876_TXNL4A', 'P20700_LMNB1', 'Q12906_ILF3', 'Q8IW45_NAXD', 'P80723_BASP1', 'P43251_BTD', 'Q92777_SYN2', 'Q9NUP9_LIN7C', 'P50226_SULT1A2', 'P22748_CA4', 'P07602_PSAP', 'O75173_ADAMTS4', 'P02671_FGA', 'Q8IYI6_EXOC8', 'P43243_MATR3', 'O43521_BCL2L11', 'P30711_GSTT1', 'P61956_SUMO2', 'Q9UKK3_PARP4', 'Q96FZ7_CHMP6', 'Q16822_PCK2', 'P02790_HPX', 'Q9BXN1_ASPN', 'P06576_ATP5F1B', 'Q92522_H1.10', 'P05155_SERPING1', 'P15088_CPA3', 'P05388_RPLP0', 'P07942_LAMB1', 'Q96I24_FUBP3', 'Q99523_SORT1', 'Q9Y3B2_EXOSC1', 'O15439_ABCC4', 'P53582_METAP1', 'P01876_IGHA1', 'P54289_CACNA2D1', 'Q5HYI8_RABL3', 'Q8NE62_CHDH', 'Q00796_SORD')
top_features_ = data.frame(top_features_[1:95])
top_features_$gene = as.character(lapply(top_features_$top_features_, function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[2]
}))
top_features_$uni = as.character(lapply(top_features_$top_features_, function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[1]
}))
top_features_$gene = gsub('.', '-', top_features_$gene, fixed = T)
entr = bitr(top_features_$gene, fromType = 'SYMBOL', toType = 'UNIPROT', OrgDb = 'org.Hs.eg.db')

rownames(top_features_) = top_features_$gene

RPPA = read.csv("Z:/members/wangyingrui/TTTD/other cohort/TCGA-RPPA-pancan-clean.txt", sep = '\t')
EBPlus_dat = read.csv("Z:/members/wangyingrui/TTTD/other cohort/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", sep = '\t')
cliinfo = read.csv("Z:/members/wangyingrui/TTTD/other cohort/clinical_PANCAN_patient_with_followup.tsv", sep = '\t')
unique(cliinfo$acronym)
cliinfo = subset(cliinfo, acronym == "THCA")
dim(cliinfo)
rownames(cliinfo) = cliinfo$bcr_patient_barcode

cols = data.frame(sample = colnames(EBPlus_dat))
cols$id = as.character(lapply(cols$sample, function(x){
  x = strsplit(x, '.', fixed = T)[[1]]
  x = paste0(x[1],'-',x[2],'-', x[3])
  x
}))
cols = cols[cols$id %in% cliinfo$bcr_patient_barcode, ]
cols = cbind(cols, cliinfo[cols$id, ])
write.xlsx(cols, "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/cli_info.xlsx")

dim(EBPlus_dat) # 20531 11070
EBPlus_dat = EBPlus_dat[c("gene_id", cols$sample)]
dim(EBPlus_dat) # 20531 573
EBPlus_dat[1:3,1:3]
EBPlus_dat$gene = as.character(lapply(EBPlus_dat$gene_id, function(x){
  x = strsplit(x, '|', fixed = T)[[1]]
  x[1]
}))
EBPlus_dat$geneId = as.character(lapply(EBPlus_dat$gene_id, function(x){
  x = strsplit(x, '|', fixed = T)[[1]]
  as.numeric(x[2])
}))
EBPlus_dat[1:3,573:575]
EBPlus_dat = EBPlus_dat[c(574:575, 2:573)]
write.xlsx(EBPlus_dat, "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2_TCGA-THCA.xlsx")

intersect(EBPlus_dat$gene, top_features_$gene)
setdiff(top_features_$gene, EBPlus_dat$gene)
entr1 = bitr(EBPlus_dat$gene, fromType = 'SYMBOL', toType = 'UNIPROT', OrgDb = 'org.Hs.eg.db')
entr1 = entr1[entr1$UNIPROT %in% entr$UNIPROT, ]
entr1 = bitr(entr1$UNIPROT, fromType = 'UNIPROT', toType = 'SYMBOL' , OrgDb = 'org.Hs.eg.db')

setdiff(top_features_$gene, entr1$SYMBOL) # "ABRACL"    "APOA4"     "JPT2"      "MACROH2A2" "CAVIN2"    "NAXD"      "ATP5F1B"   "H1-10"
interPros = unique(entr1$SYMBOL)
# [1] "A2M"      "ABCC4"    "ABLIM1"   "ACO1"     "ADAMTS4"  "ADH1B"    "ALDH4A1"  "APCS"     "APLP2"    "ASPN"     "BASP1"    "BCL2L11"  "BTD"      "BUB3"    
# [15] "CA4"      "CAPG"     "CAPZA1"   "CD81"     "CDH1"     "CFH"      "CHMP6"    "CLN5"     "CMBL"     "COG2"     "CPA3"     "CUTA"     "DBT"      "DPP7"    
# [29] "ESD"      "ETFDH"    "EXOC8"    "EXOSC1"   "F11R"     "F13A1"    "FGA"      "FLOT1"    "FSTL1"    "FTH1"     "FUBP3"    "GAA"      "GAL3ST3"  "GSTT1"   
# [43] "GTF2A1"   "HDGF"     "HLA-A"    "HLA-DRA"  "HPX"      "HTRA1"    "ILF3"     "IVD"      "KDSR"     "LAMB1"    "LDHD"     "LGALS3BP" "LIN7C"    "LMNB1"   
# [57] "LYZ"      "MATR3"    "MCFD2"    "METAP1"   "NLRP6"    "NUDCD1"   "PARP4"    "PCK2"     "PPP1R14C" "PSAP"     "RABEP1"   "RPLP0"    "RPS27A"   "S100A6"  
# [71] "SAA1"     "SEH1L"    "SERPINA3" "SERPING1" "SLC43A3"  "SNCA"     "SORT1"    "SSR3"     "STRBP"    "SULT1A2"  "SUMO2"    "SYN2"     "TAGLN2"   "TF"      
# [85] "TPPP"     "TXNL4A"   "VTN" 

dat =subset(EBPlus_dat, gene %in% unique(entr1$SYMBOL))
rownames(dat) = dat$gene
dat = dat[c( 3:ncol(dat))]
dat = data.frame(t(dat))

dat1 = data.frame(row.names = rownames(dat))
dat1[c("ABRACL", "APOA4", "JPT2", "MACROH2A2", "CAVIN2", "NAXD", "ATP5F1B", "H1-10")] = 0

dat = cbind(dat1, dat)
dat$Age = cols$age_at_initial_pathologic_diagnosis
colnames(dat)[1:95] = top_features_[gsub('.', '-',colnames(dat)[1:95], fixed = T),  'top_features_.1.95.']
write.xlsx(dat, "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlusPlus_THCA_95feat.xlsx", rowNames  =T)


dat = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlusPlus_THCA_95feat_pred.xlsx", rowNames  =T )
dat$BCRstatus = ifelse(cols$vital_status=='Alive', 1, 0)
dat$time = ifelse(dat$BCRstatus == 1, cols$days_to_last_followup , cols$days_to_death)
dat[] <- lapply(dat, as.numeric)
table( dat$BCRstatus)
dat$risk = dat$Pred - dat$Age
head(dat)
dat$riskzscore = (dat$risk-mean(dat$risk))/sd(dat$risk)

multicox <- coxph(Surv(time = time, event = BCRstatus) ~ risk , data = dat  )
multicox


sur.cut_test <- surv_cutpoint(dat , time= 'time', event = 'BCRstatus' , variables = 'riskzscore' )
sur.cut_test$cutpoint
sur.cat_test <- surv_categorize(sur.cut_test)
sur.cat_test$riskscore = dat$risk
coxph(Surv(time = time, event = BCRstatus ) ~ risk , data = sur.cat_test )

# multicox <- coxph(Surv(time = time, event = BCRstatus) ~  . , data = dat  )
# multicox
# datcoef = dat
# for(i in names(multicox$coefficients)){
#   datcoef[i] = datcoef[i] * multicox$coefficients[ i]
# }
# datcoef$risk = rowSums(datcoef[names(multicox$coefficients)], na.rm = T)
# 
# sur.cut_test <- surv_cutpoint(datcoef , time= 'time', event = 'BCRstatus' , variables = 'risk' )
# sur.cat_test <- surv_categorize(sur.cut_test)
# sur.cat_test$riskscore = datcoef$risk
# rownames(sur.cat_test) = rownames(datcoef)
# sur.cat_test$riskscore = datCoef$risk
#write.xlsx(sur.cat_test, 'D:/chh/2023workProject/prottalk/code/shenshi_tech/featurevalidation/TCGA_risk.xlsx', rowNames= T)
#sur.cat_test$risk = 'risk'
head(sur.cat_test)
coxph(Surv(time = time, event = BCRstatus ) ~ risk , data = sur.cat_test )

# multicoxsur <- coxph(Surv(time = time, event = BCRstatus) ~  risk , data = datcoef  )
# summary(multicoxsur)

# 
sur.cat_test1 <- surv_cutpoint(dat , time= 'time', event = 'BCRstatus' , variables = 'risk' )
sur.cat_test1$cutpoint
sur.cat_test1 <- surv_categorize(sur.cat_test1)
sur.cat_test1$risk = factor(sur.cat_test1$risk, c('low', 'high'))
head(sur.cat_test1)
sur.cat_test1$riskscore = dat$risk
unicox <- coxph(Surv(time = time, event = BCRstatus ) ~ riskscore , data = sur.cat_test1 )
unicox
ggforest(unicox, data = sur.cat_test1,
         cpositions = c(0.01,0.1,0.3),
         fontsize=1,
         noDigits=2)
# 
# sur.cat_test1 = dat
# sur.cat_test1$risk = ifelse(sur.cat_test1$risk>0, 'high', 'low')
# sur.cat_test1$risk = factor(sur.cat_test1$risk, c('high', 'low'))
fit = survfit(Surv(time , BCRstatus ) ~ risk, data = sur.cat_test )
fit
sumfit = summary(fit )
sumfit
# sumfit$surv
# fitdata <- data.frame( time = fit$time, nrisk = fit$n.risk, nevent= fit$n.event, survival = fit$surv, lower = fit$lower, upper = fit$upper ) 
pdf( paste0(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlusPlus_THCA_95feat_pred.pdf")) )
p = ggsurvplot(fit, pval = T, conf.int = TRUE, 
           #risk.table = TRUE, # Add risk table
           #censor = T,# Add censor table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", ggtheme = theme_classic(), palette = c("#E7B800", "#2E9FDF"), 
           risk.table.fontsize = 5,
           cumcensor = TRUE, #加号数
           # font.x= 15
           #, title = "TCGA"
           # surv.median.line = "hv", # Specify median survival
)
print(p)
dev.off()
dev.off()
dev.off()
dev.off()

temp = dat
head(temp)
temp$BCRstatus = ifelse(temp$BCRstatus==0, 'No', 'Yes')
ggplot(temp, aes(x = BCRstatus, y = riskzscore, color= BCRstatus )) +
  #geom_point(   )+ # stroke = 0: 设置边框宽度为 0，以确保没有边框显示
  geom_boxplot(lwd = 1) + geom_jitter( size=1, alpha = 0.2)+
  labs(y = 'Age Gap')+
  scale_color_d3() +
  theme_classic()+
  geom_signif(comparisons=list(c('No', 'Yes')), step_increase = 0.3, test="wilcox.test",  map_signif_level=T   # 标签样式F为数字，T为*号
  ) + theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


### TPR
delinfo = subset(machine_info, Thyroid_ID %in% TTTDTPR_sur$Thyroid_ID)
dim(delinfo)
dim(TTTDTPR_sur)
for(i in 1:nrow(delinfo)){
  sid = delinfo[i, 'sample']
  tid = delinfo[i, 'Thyroid_ID']
  temp = subset(TTTDTPR_sur, Thyroid_ID == tid)
  #print(c(i, nrow(temp)))
  delinfo[i, 'time'] = temp$Follow_up_time
  delinfo[i, 'BCR'] = temp$Recurrence
}
delinfo[c('Age', 'Pred')] = xgb_results[rownames(delinfo), c('Age', 'Pred')]
delinfo$Age_gap = delinfo$Pred-delinfo$Age
dim(delinfo)#[1] 87 21
dat = delinfo[c('Age_gap', 'time', 'BCR')]
head(dat)
colnames(dat) = c('risk', 'time', 'BCR')
dat[] = lapply(dat, as.numeric )
TPRdat = dat

sur.cut_test <- surv_cutpoint(dat , time= 'time', event = 'BCR' , variables = 'risk' )
sur.cat_test <- surv_categorize(sur.cut_test)
sur.cat_test$riskscore = dat$risk
head(sur.cat_test)


multicoxsur <- coxph(Surv(time = time, event = BCR) ~  risk , data = dat  )
summary(multicoxsur)


sur.cat_test1 = sur.cat_test
sur.cat_test1$risk = factor(sur.cat_test1$risk, c('low', 'high'))
unicox <- coxph(Surv(time = time, event = BCR ) ~ risk , data = sur.cat_test1 )
unicox
ggforest(unicox, data = sur.cat_test, # main = "Hazard ratio",
         cpositions = c(0.01,0.1,0.3),
         fontsize=1,
         noDigits=3)

# sur.cat_test1 = dat
# sur.cat_test1$risk = ifelse(sur.cat_test1$risk>0, 'high', 'low')
sur.cat_test1$risk = factor(sur.cat_test1$risk, c('high', 'low'))
fit = survfit(Surv(time , BCR ) ~ risk, data = sur.cat_test1 )
ggsurvplot(fit, pval = T, conf.int = TRUE, risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", ggtheme = theme_bw()#, title = "TCGA"
           # surv.median.line = "hv", # Specify median survival
)
head(sur.cat_test)
temp = sur.cat_test
temp$BCR = ifelse(temp$BCR==0, 'No', 'Yes')
ggplot(temp, aes(x = BCR, y = riskscore, color= BCR )) +
  #geom_point(   )+ # stroke = 0: 设置边框宽度为 0，以确保没有边框显示
  geom_boxplot(lwd = 1)+geom_jitter( size=1, alpha = 0.1)+
  scale_color_d3() +labs(y = 'Age Gap')+
  theme_classic()+geom_signif(comparisons=list(c('No', 'Yes')), #vecs
                              step_increase = 0.3,
                              test="wilcox.test",  # "t 检验，比较两组（参数）" = "t.test" "wilcox.test"
                              map_signif_level=T   # 标签样式F为数字，T为*号
  ) +theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

c('#46AA6E', "#FF7F0E" , '#9C7822')
