library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(reshape2)
library(ggrepel)
BiocManager::install("biomaRt")
library(biomaRt)

stdscale = function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

coefficients1$Uniprot = sapply(rownames(coefficients1), function(x){
  strsplit(x,'_')[[1]][1]
})
coefficients1$Gene = sapply(rownames(coefficients1), function(x){
  strsplit(x,'_')[[1]][2]
})
coefficients1 = cbind(coefficients1[grepl('qt',colnames(coefficients1))], coefficients1[!grepl('qt',colnames(coefficients1))])
coefficients1[coefficients1$Gene=='Tkt',]



setwd("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid")


mouseinfo = read.csv("Z:/members/wangyingrui/TTTD/Mouse/20231123_AgingPDM_Thyroid_sample_information.CSV" )
mousedat = read.csv("Z:/members/wangyingrui/TTTD/Mouse/Aging_combineMatrix_20230609.csv", row.names = 1 )
mouseinfoA1 = mouseinfo[mouseinfo$Sample_name %in% colnames(mousedat), ]
mouseinfoA1
data.frame(table(mouseinfoA1$Month))

mousedat[1:3,1:3]
mouseinfo[mouseinfo$Group == 'A1', 'Sample_name']
rownames(mouseinfo) = mouseinfo$Sample_name
mousedat = mousedat[colnames(mousedat)[grepl('b23', colnames(mousedat))]]
mousedat$uni = as.character(lapply(rownames(mousedat), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[1]
}))
mousedat$gene = as.character(lapply(rownames(mousedat), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[2]
}))


human <- useMart('ensembl', dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl', dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
listAttributes(human)
m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                values = mousedat$gene , mart = mouse,
                attributesL = c("hgnc_symbol", "chromosome_name","start_position"), martL = human,uniqueRows = T)

m2h.g = data.frame(m2h.g)
m2h.g = subset(m2h.g, HGNC.symbol!='')
head(m2h.g)
m2h.g1 = m2h.g[m2h.g$HGNC.symbol %in% top_features_$gene, ]
dim(m2h.g1)
for(i in 1:nrow(m2h.g1)){
  g = m2h.g1[i, 'MGI.symbol']
  m2h.g1[i,'unigene'] = rownames(subset(mousedat, gene == g ))
}
m2h.g1 =  m2h.g1[ m2h.g1$unigene %in% colnames(mouse_log2_min0.8_zscore), ]
top_features_[!(top_features_$gene %in% m2h.g1$HGNC.symbol), 1 ]
mouse_ML = mouse_log2_min0.8_zscore[c('Age', m2h.g1$unigene)]

mouse_ML$Q03734_Serpina3m = rowMeans(mouse_log2_min0.8_zscore[m2h.g1[m2h.g1$HGNC.symbol == 'SERPINA3', 5]], na.rm = T)

library(dplyr)
m2h.g1 = m2h.g1 %>% distinct(HGNC.symbol, .keep_all = TRUE)
# rownames(m2h.g1)  = m2h.g1$MGI.symbol
mouse_ML = mouse_ML[c('Age', m2h.g1$unigene)]

for(i in 1:nrow(m2h.g1)){
  g = m2h.g1[i, 'HGNC.symbol']
  m2h.g1[i,'human_unigene'] = top_features_[top_features_$gene == g, 1]
}
colnames(mouse_ML)[2:ncol(mouse_ML)] = m2h.g1$human_unigene
top_features_[!(top_features_$gene %in% m2h.g1$HGNC.symbol), 1 ]
mouse_ML[top_features_[!(top_features_$gene %in% m2h.g1$HGNC.symbol), 1 ]] = 0
mouse_ML[1:2,1:3]
write.xlsx(mouse_ML, "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_MLzscore_95feat.xlsx", rowNames = T)

predictResult = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_MLzscore_95feat_pred.xlsx", rowNames = T)
predictResult$Age_gap = predictResult$Pred-predictResult$Age
tempcor = predictResult[c('Age', 'Pred')]
tempcor$type = 'mouse'
pearsonr = cor(tempcor[1:2], method = 'pearson')
spearmanr = cor(tempcor[1:2], method = 'spearman')
tempcor$Age

ggplot(tempcor, aes(Age, Pred, color = type))+
  geom_smooth(aes(fill = type), alpha = 0.15, method = 'lm' )+ # se = F
  geom_point()+
  labs(y = 'predict',x = 'real')+
  #labs(title = unique(tempcor$TrainTest))+
  theme_classic()+ 
  annotate('text',x = 10, y= 55, label = paste0('Pearson r=', round(pearsonr[2], 3) ))+  # 
  annotate('text',x = 10, y= 58, label = paste0('Spearman r=', round(spearmanr[2], 3) ) )+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15),
        legend.position = "none")


mousedat_log2 = log2(mousedat[1:25])
narow = rowSums(is.na(mousedat_log2))
mousedat_log2 = mousedat_log2[narow!=ncol(mousedat_log2),]

mouse_log2_min0.8 = mousedat_log2
mouse_log2_min0.8[is.na(mouse_log2_min0.8)] = min(mouse_log2_min0.8, na.rm = T) * 0.8
mouse_log2_min0.8[1:20,1:2]
mouse_log2_min0.8 = data.frame(t(mouse_log2_min0.8))

dim(mouse_log2_min0.8)
mouse_log2_min0.8$Age = as.numeric(mouseinfo[rownames(mouse_log2_min0.8), 'Month'])
dim(mouse_log2_min0.8)
mouse_log2_min0.8 = cbind(mouse_log2_min0.8['Age'], mouse_log2_min0.8[,1:6685])
mouse_log2_min0.8[1:10,1:2]
mouse_log2_min0.8_zscore = mouse_log2_min0.8
mouse_log2_min0.8_zscore[,2:ncol(mouse_log2_min0.8_zscore)] = apply(mouse_log2_min0.8_zscore[,2:ncol(mouse_log2_min0.8_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})
dim(mouse_log2_min0.8_zscore)
mouse_log2_min0.8_zscore[1:5,1:4]
############ DEswan 1 window  ############
type1 = 'PDTC'
# dat = mouse_log2_min0.8[rownames(delinfo[delinfo$Classificaiton_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
# dat = mouse_log2_min0.8[rownames(delinfo[delinfo$Histopathology_type == type1 & delinfo$TrainTest == 'Discovery', ]),]
glm_overlap =list()
for(type1 in c('')){#  "PTC" "FTC" "ATC" "MTC" "PDTC", 'Borderline'
  #dat = mouse_log2_min0.8_zscore
  dat = EBPlus_dat_zscore
  #dat[] <- lapply(dat, as.numeric)
  cor_result=cor(dat[,1], dat[, 3:ncol(dat)], method = 'spearman')
  cor_result[1:10]
  cor_cut = 0.1
  temp1 = colnames(cor_result)[abs(cor_result)>cor_cut]
  # dat1 = dat[,5+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>cor_cut]))]
  dat1 = dat[,2+c(which(colnames(cor_result) %in% colnames(cor_result)[abs(cor_result)>cor_cut]))]
  dat1[1:3,1:3]
  #covariates = dat[,c(2,4,5)]
  covariates = dat[2]
  head(covariates)
  qt = dat[,1] # Age
  pvalues  <-  NULL
  coefficients  <-  NULL
  qt <- factor(qt)
  for(i in 1:ncol(dat1)){
    deswan.formula  <-  "dat1[, i] ~ qt"
    # deswan.formula  <-  paste(c("dat1[, i] ~ qt", paste("covariates$",colnames(covariates), collapse  =  " + ", sep = "")), collapse = " + ", sep = "")
    # deswan.formula = "dat1[, i] ~ qt + covariates$Gender " # + covariates$SampleType + covariates$Hospital
    # deswan.formula  <-  "dat1[, i] ~ qt + covariates$Gender"
    test.glm  <- try(glm.fit <- glm(as.formula(deswan.formula), family  =  gaussian), silent=TRUE)
    if(class(test.glm)[1] !=  "try-error"){
      print(c(i, colnames(dat1)[i]))
      glm.res <- car::Anova(glm.fit,  type = "2")
      pvalues <- rbind(pvalues, data.frame(variable = colnames(dat1)[i], factor = rownames(glm.res), pvalue=glm.res$`Pr(>Chisq)`, stringsAsFactors  =  F))
      coefficients <- rbind(coefficients, data.frame(variable  =  colnames(dat1)[i] , factor = names(coefficients(glm.fit)), 
               coefficient=coefficients(glm.fit), stringsAsFactors  =  F))
    }
  }
  pvalues = data.frame(pvalues)
  pvalues_p = data.frame(acast(pvalues[c("variable", "factor", "pvalue")], variable~factor))
  pvalues_padjust = data.frame(apply(pvalues_p, 2, function(x){p.adjust(x, method = 'BH')}))
  pvalues_p$p_adjust = p.adjust(pvalues_p$qt, method = 'BH')
  pvalues_p$gender_p_adjust = p.adjust(pvalues_p$covariates.Gender, method = 'BH')
  # temp = subset(pvalues_padjust, qt < 0.05 )
  # if(nrow(temp)>0){
  #   glm_overlap[[type1]] = rownames(temp)
  # }
  head(pvalues, 10)
  coefficients = data.frame(coefficients)
  head(coefficients)
  coefficients1 = acast(coefficients, variable~factor)
  coefficients1 = data.frame(coefficients1)
  head(coefficients1)
  coefficients1$qt_coef = rowMeans(coefficients1[grepl('qt', colnames(coefficients1))], na.rm = T)
  
  coefficients1 = cbind(coefficients1, pvalues_p[rownames(coefficients1),])
  rownames(coefficients1) = gsub('.','-', rownames(coefficients1), fixed = T)
  coefficients1$Uniprot = sapply(rownames(coefficients1), function(x){
    strsplit(x,'_')[[1]][1]
  })
  coefficients1$Gene = sapply(rownames(coefficients1), function(x){
    strsplit(x,'_')[[1]][2]
  })
  
  
  # coefficients1 = coefficients1[c('qt', colnames(coefficients1)[!grepl('qt', colnames(coefficients1))]) ]
  # # coefficients(glm.fit), summary(glm.fit)$coefficients, glm.fit$coefficients
  # # pval_coefbind(pvalues, coefficients1, type1)
  # pval_coefbind(pvalues, coefficients1, type1, )
}
pvalues_p[rownames(coefficients_TCGA),]

uni_gene = mousedat[26:27]
uni_gene[uni_gene$uni=='P43275',]

human <- useMart('ensembl', dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl', dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                values = uni_gene$gene, mart = mouse,
                attributesL = c("hgnc_symbol","chromosome_name","start_position"), martL = human,uniqueRows = T)

head(m2h.g)
m2h.g = data.frame(m2h.g)
write.xlsx(m2h.g, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/m2h_gene.xlsx', rowNames = F)
m2h.g1 = subset(m2h.g, HGNC.symbol!='')
library(dplyr)
m2h.g1 = m2h.g1 %>% distinct(MGI.symbol, .keep_all = TRUE)
rownames(m2h.g1)  = m2h.g1$MGI.symbol


coefficients1$mouse_gene = m2h.g1[coefficients1$Gene, c( 'MGI.symbol')]
coefficients1$human_gene = m2h.g1[coefficients1$Gene, c('HGNC.symbol')]

# write.xlsx(coefficients1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/Mousezscore_allWindows.xlsx', rowNames = T)
# write.xlsx(coefficients1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/Mousezscore_corcut0_allWindows.xlsx', rowNames = T)

coefficients1 = read.xlsx( 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/Mousezscore_corcut0_allWindows.xlsx', rowNames = T)

temp = subset(coefficients1, p_adjust<0.05)
temp$Sig = ifelse(is.na(temp$qt_coef), "None", ifelse(temp$qt_coef>0 & temp$p_adjust<0.05, "Up",
                                                            ifelse( temp$qt_coef< 0& temp$p_adjust<0.05, "Down", "None")))
table(temp$Sig)
ggplot(temp, #[temp$p_adjust>10e-70, ], 
       aes(x = qt_coef, y = -log10(p_adjust), colour=Sig)) +#  , size = -log10(qt_padjust)
  geom_point(alpha=0.7, size=2 ) +
  scale_color_manual(values=c("#546de5","#ff4757"))+# , "#d2dae2"
  geom_text_repel( data = temp[temp$p_adjust<0.05& temp$p_adjust>10e-70 & abs(temp$qt_coef)>0, ],
                   aes(x = qt_coef, y = -log10(p_adjust), label = human_gene),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
                   segment.color = "black", show.legend = FALSE )+
  labs(x="coefficient", y="-log10 (P-adjust)", title = 'Mouse' )+
  theme_classic()+ scale_y_continuous(n.breaks = 10) + #ylim(0, 60)+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

####
human_n = read.xlsx("D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_N_glm_ageContinue.xlsx", rowNames= T)
human_n$Gene = sapply(rownames(human_n), function(x){
  x = strsplit(x,'_')[[1]][2]
  x = gsub('.','-',x, fixed = T)
  x
})


human_n = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEunnormalZscore_Discovery_N_glm_ageContinue.xlsx', rowNames = T)
#human_n[pro_lists,]
human_n$Gene = sapply(rownames(human_n), function(x){
  x = strsplit(x,'_')[[1]][2]
  x = gsub('.','-',x, fixed = T)
  x
})
### monkey ####
monkeydat = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/monkey1-s2.0-S0092867424009140-mmc2.xlsx", sheet = 'gene expression trends' )
unique(monkeydat$X3)
monkeydat=monkeydat[monkeydat$X4 == 'Thyroid gland' & monkeydat$X3 == "protein_coding",]
monkeydat = subset(monkeydat, !is.na( X2) & monkeydat$X5 %in% c("Cluster U", "Cluster D" ))

unique(monkeydat$X5)
monkeyUp = subset(monkeydat, monkeydat$X5 == "Cluster U" & !is.na(monkeydat$X2))
monkeyDown = subset(monkeydat, monkeydat$X5 == "Cluster D"& !is.na(monkeydat$X2))
humanUp = subset(human_n, qt_coef > 0 & qt_padjust<0.05 & !is.na(human_n$Gene))
humanDown = subset(human_n, qt_coef< 0 & qt_padjust<0.05& !is.na(human_n$Gene))
mouseUp = subset(coefficients1, p_adjust<0.05 &  qt_coef>0 & !is.na( human_gene))
mouseDown = subset(coefficients1, p_adjust<0.05 &  qt_coef<0 & !is.na(human_gene))
temp = subset(coefficients1, p_adjust<0.05)
temp = subset(coefficients1, !is.na(coefficients1$human_gene))

intersect(mouseUp$human_gene, humanUp$Gene)
intersect(mouseDown$human_gene, humanDown$Gene)

machine_info=read.xlsx('D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', rowNames = T)
machine_info$Histopathology_type = gsub('PTMC', 'PTC', machine_info$Histopathology_type)
delinfo = subset(machine_info, Tissue_type=='thyroid' & Histopathology_type!='uncertain' & SampleType!='FNA_in_vivo')
delinfo = delinfo[!grepl('ool', rownames(delinfo)),]
delinfo = delinfo[!grepl('mouse', rownames(delinfo)),]
delinfo = delinfo[!grepl('qc', rownames(delinfo)),]
delinfo = subset(delinfo, Histopathology_type!='MEC' )
delinfo = subset(delinfo, qe480!='480' )
QE_deldat_IQR_min0.8 = read.xlsx('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss_naImpute_matrix.xlsx', rowNames = T)
QE_deldat_IQR_min0.8 = QE_deldat_IQR_min0.8[intersect(rownames(QE_deldat_IQR_min0.8), rownames(delinfo)),]

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
QE_deldat_IQR_min0.8_zscore[1:3,1:6]
QE_deldat_IQR_min0.8_zscore$class = machine_info[rownames(QE_deldat_IQR_min0.8_zscore), 'Classificaiton_type']


up = intersect(mouseUp$human_gene, humanUp$Gene)
down = intersect(mouseDown$human_gene, humanDown$Gene)

humanUp$unigene = rownames(humanUp)
rownames(humanUp) = humanUp$Gene

humanDown$unigene = rownames(humanDown)
humanDown = subset(humanDown, !is.na( Gene))
rownames(humanDown) = humanDown$Gene

mouseUp = subset(mouseUp, !is.na( human_gene))
mouseUp$unigene = rownames(mouseUp)

mouseDown = subset(mouseDown, !is.na( human_gene))
mouseDown$unigene = rownames(mouseDown)
mouseDown = subset(mouseDown, mouseDown$human_gene %in% humanDown$Gene)

pros = down
point_plot = QE_deldat_IQR_min0.8_zscore[QE_deldat_IQR_min0.8_zscore$class == 'N', c('Age', humanDown[pros, 'unigene'])]
colnames(point_plot)[2:ncol(point_plot)] = pros
point_plot = melt(point_plot, id.vars = 'Age')
point_plot$species = 'Human'

point_plot = mouse_log2_min0.8_zscore[ c('Age', gsub('-','.', mouseUp$unigene, fixed = T))]
mouseUp[mouseUp$human_gene == 'KLK2', 'unigene']
point_plot$P00755_Klk1b1 = rowMeans(point_plot[mouseUp[mouseUp$human_gene == 'KLK2', 'unigene']], na.rm = T)
mouseUp = mouseUp%>% distinct( human_gene, .keep_all = T)
rownames(mouseUp) = mouseUp$human_gene
mouseUp[mouseUp$human_gene == 'KLK2', 'unigene']
point_plot = point_plot[c('Age', mouseUp[pros, 'unigene'])]
colnames(point_plot)[2:ncol(point_plot)] = pros
point_plot = melt(point_plot, id.vars = 'Age')
point_plot$species = 'mouse'

point_plot = mouse_log2_min0.8_zscore[ c('Age', gsub('-','.', mouseDown$unigene, fixed = T))]
mouseDown[mouseDown$human_gene %in% c('RBM12B', 'SERPINA1', 'SERPINA3'),] #  'unigene'
point_plot$P07759_Serpina3k = rowMeans(point_plot[mouseDown[mouseDown$human_gene == 'SERPINA3', 'unigene']], na.rm = T)
mouseDown = mouseDown%>% distinct( human_gene, .keep_all = T)
rownames(mouseDown) = mouseDown$human_gene
point_plot = point_plot[c('Age', mouseDown[pros, 'unigene'])]
colnames(point_plot)[2:ncol(point_plot)] = pros
point_plot = melt(point_plot, id.vars = 'Age')
point_plot$species = 'mouse'


head(point_plot)
colnames(point_plot ) = c('Age', "Protein", 'Intensity', "species")
ggplot(point_plot, aes(Age, Intensity, color = Protein))+
  geom_smooth(aes(fill = Protein),method = "loess", se = F,  alpha= 0.15, linewidth= 0.01 )+#method = "glm", aes( fill = Protein)
  theme_classic()+
  # scale_color_d3()+
  # scale_fill_d3()+
  labs(y = 'log2 std Expressions', title = unique(point_plot$species))+
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        text = element_text(size = 13))


library(UpSetR)
library(ggVennDiagram)
input = list(mouse =temp[temp$p_adjust<0.05, 'human_gene'] , monkey = monkeydat$X2, human = human_n[human_n$qt_padjust<0.05, 'Gene']  )

input = list(mouse = mouseUp[!is.na(mouseUp$human_gene), 'human_gene'], monkey = monkeyUp[!is.na(monkeyUp$X2), 'X2'], human = humanUp[!is.na(humanUp$Gene), 'Gene'] )
input = list(mouse = mouseDown[!is.na(mouseDown$human_gene), 'human_gene'], monkey = monkeyDown[!is.na(monkeyDown$X2), 'X2'], human = humanDown[!is.na(humanDown$Gene), 'Gene'] )
ggVennDiagram(input, #category.names = c("A","B","C"),
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              edge_lty = "dashed", 
              edge_size = 1) + 
  # labs(title = 'up')+
  #theme(legend.position = 'top')+labs(title = names(x)[1])+
  scale_fill_gradient(low="white",high = "white",name = "count")#library(ggplot2)

library(VennDiagram)
df_inter <- get.venn.partitions(input)
for (i in 1:nrow(df_inter)) df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ', ')
# df_inter[-c(5, 6)]
df_inter = data.frame(df_inter)
############### sankey ##############
pathways=read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid2_monkey/QEunnormalZscore_Discovery_metascape_result.xlsx", rowNames = T)
colnames(pathways)
pathway = c()
gene = c()
for(i in 1:nrow(pathways)){
  temp = c(strsplit(pathways[i, 10], ';')[[1]], strsplit(pathways[i, 15], '; ')[[1]], strsplit(pathways[i, 16], '; ')[[1]] ) #
  pathway = c(pathway, temp)
  gene = c(gene, rep(pathways[i, 6], length(temp)))
}
pathway = data.frame(gene, pathway )
pathway = subset(pathway, !is.na(pathway))
unique(pathway$gene)
head(pathway)
library(ggalluvial)
library(dplyr)
# summarize the data and count the frequencies
frequencies <- pathway %>%
  count(gene, pathway) #%>% arrange(gene, desc(n))
ggplot(data = frequencies[!grepl('GO:', frequencies$pathway),], aes(axis1 = gene, axis2 = pathway, # y = n  Third variable on the X-axis
                      )) +
  geom_alluvium(aes(fill = gene)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=4) +
  # scale_fill_igv() +
  scale_fill_d3(palette =  "category20")+ 
  theme_void()+theme(text = element_text(size = 5))

#### transicr ####
EBPlus_dat = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2_TCGA-THCA.xlsx")
EBPlus_dat[1:3,1:3]
rownames(EBPlus_dat) = paste0('Gene', EBPlus_dat$geneId)
EBPlus_dat = EBPlus_dat[3:ncol(EBPlus_dat)]

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

dim(EBPlus_dat)
EBPlus_dat = data.frame(t(EBPlus_dat))

EBPlus_dat = cbind( cols$age_at_initial_pathologic_diagnosis, EBPlus_dat)
EBPlus_dat_zscore = EBPlus_dat
EBPlus_dat_zscore[,2:ncol(EBPlus_dat_zscore)] = apply(EBPlus_dat_zscore[,2:ncol(EBPlus_dat_zscore)], 2, function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
})
dim(EBPlus_dat_zscore)

EBPlus_dat_zscore[] <- lapply(EBPlus_dat_zscore, as.numeric)
EBPlus_dat_zscore = cbind( EBPlus_dat_zscore)
EBPlus_dat_zscore[1:3,1:3]
EBPlus_dat_zscore = EBPlus_dat_zscore[c(2,1,3:ncol(EBPlus_dat_zscore))]
colnames(EBPlus_dat_zscore)[1:2] = c('Age', 'Gender')
EBPlus_dat_zscore$Age = as.numeric(EBPlus_dat_zscore$Age)
EBPlus_dat_zscore[c( 3:ncol(EBPlus_dat_zscore))]<- lapply(EBPlus_dat_zscore[c( 3:ncol(EBPlus_dat_zscore))], as.numeric)

entr = bitr(as.numeric(gsub('Gene', '', rownames(coefficients_TCGA1))), fromType = 'ENTREZID', toType = c(  'SYMBOL'), OrgDb = 'org.Hs.eg.db')
rownames(entr) = paste0('Gene', entr$ENTREZID)
coefficients_TCGA1[c(  'SYMBOL')] = entr[rownames(coefficients_TCGA1), c( 'SYMBOL')]

write.xlsx(coefficients_TCGA1 , 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/EBPlus_zscore_allWindows.xlsx', rowNames = T)


### 
MouseCoef = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/Mousezscore_allWindows.xlsx", rowNames = T)
MouseCoef = MouseCoef[9:ncol(MouseCoef)]
MouseCoef1 = subset(MouseCoef, !is.na(human_gene) & p_adjust<0.05)
# MouseCoef1 = MouseCoef1[MouseCoef1$human_gene %in% human_n$Gene, ]
MouseCoef1 = MouseCoef1[order(MouseCoef1[,3], decreasing = F),]
MouseCoef1 = MouseCoef1 %>% distinct( human_gene, .keep_all = TRUE)
rownames(MouseCoef1) = MouseCoef1$human_gene

human_n1 = subset(human_n, !is.na(Gene) & qt_padjust<0.05)
#human_n1 = human_n1[human_n1$Gene %in% MouseCoef1$human_gene, ]
human_n1$unigene = rownames(human_n1)
rownames(human_n1) = human_n1$Gene

inter = intersect(human_n1$Gene, MouseCoef1$human_gene)

temp = human_n1[inter, c( 'qt_padjust', 'qt_coef')]
temp[inter,  c('mouse_q', 'mouse_coef')] = MouseCoef1[inter, c('p_adjust', 'qt_coef' )]
temp$Gene = rownames(temp)
temp = temp[temp$qt_padjust<0.05& temp$mouse_q<0.05,]
temp$sig = 'No'
temp[temp$qt_coef<0 & temp$mouse_coef< -0.1, 'sig'] = 'Yes'
temp[temp$qt_coef>0 & temp$mouse_coef>0.1, 'sig'] = 'Yes'

ggplot(temp, #[temp$p_adjust>10e-70, ], 
       aes(x = qt_coef, y = mouse_coef, color = sig )) +#  , size = -log10(qt_padjust)
  geom_point(alpha=0.7, size=1 ) +
  #geom_text(data=temp[temp$sig == 'Yes', ], aes(x=qt_coef, y = mouse_coef ,label=Gene),size=2,vjust=-0.5,position = position_dodge(0.3))+
  geom_text_repel(data=temp[temp$sig == 'Yes', ], aes(x=qt_coef, y = mouse_coef ,label=Gene),  size = 3,  color = "black", # 不固定颜色使其根据Sig修改颜色
                  segment.color = "black", show.legend = FALSE )+
  scale_color_manual(values=c(  "#d2dae2","#ff4757"))+# , "#d2dae2"
  # geom_text_repel( data = temp[temp$p_adjust<0.05& temp$p_adjust>10e-70 & abs(temp$qt_coef)0, ],
  #                  aes(x = qt_coef, y = -log10(p_adjust), label = human_gene),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
  #                  segment.color = "black", show.legend = FALSE )+
  labs(x="human coefficient", y="Mouse coefficient", title = '' )+
  theme_classic()+ #scale_y_continuous(n.breaks = 10) + #ylim(0, 60)+
  geom_hline(yintercept = 0.1, linetype = "dashed")+
  geom_hline(yintercept = -0.1, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
# pheatmap
library(pheatmap)
heatmapmouse = mousedat_log2[paste0(MouseCoef1$Uniprot, '_', MouseCoef1$Gene),]
heatmapmouse = data.frame(t(heatmapmouse))
Age = mouseinfo[rownames(heatmapmouse), 'Month']
heatmapmouse = cbind(Age, heatmapmouse)
heatmapmouse[1:3,1:3]
colnames(heatmapmouse)[2:ncol(heatmapmouse)] = MouseCoef1$human_gene
heatmapHuman = QE_deldat_IQR_min0.8[intersect(rownames(machine_info[machine_info$Classificaiton_type =='N', ]), 
                                              rownames(QE_deldat_IQR_min0.8)),]

heatmapHuman = heatmapHuman[c('Age', human_n1$unigene)]
heatmapHuman[1:3,1:3]
colnames(heatmapHuman)[2:ncol(heatmapHuman)] = human_n1$Gene
heatmapHuman = heatmapHuman[order(heatmapHuman[,1]),]
heatmapmouse = heatmapmouse[order(heatmapmouse[,1]),]
dat = heatmapHuman
dat[rownames(heatmapmouse), colnames(heatmapmouse)] = heatmapmouse

dat[1:3,1:4]

pheatmap(dat[2:ncol(dat)], scale = 'column',
         breaks = c(seq(-1, 1, length=1000)),
         color = colorRampPalette(c("#06C2CE", "black","yellow" ))(1000), 
         #color = colorRampPalette(c("blue", "white","red" ))(1000), 
         show_colnames = T,
         show_rownames = F,
         fontsize_col = 2,
         cellwidth= 3, 
         cellheight= 0.5,
         cluster_rows = F,
         cluster_cols = F)

### mouse anova ####
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


mousedat = read.csv("Z:/members/wangyingrui/TTTD/Mouse/PDM_combineMatrix_20230609.csv", row.names = 1 )
rownames(mousedat) = gsub(';$', '', mousedat$X)
mousedat[1:5,1:5]
mouseinfo = read.csv("Z:/members/wangyingrui/TTTD/Mouse/20231123_AgingPDM_Thyroid_sample_information.CSV" )
rownames(mouseinfo) = mouseinfo$Sample_name
mouseinfo = mouseinfo[mouseinfo$Sample_name%in% colnames(mousedat), ]
dim(mouseinfo)
# temp = subset(mouseinfo, Sex == "Male" & Diet=='M')
# freq = data.frame(table(temp$Month))
# freq$sex = 'Male'
# freq$Diet = 
freq = c()
for(g in c("H", 'M', "L")){
  # g = 'M'
  for(sex in c("Male" , "Female")){
    temp = subset(mouseinfo, Sex == sex & Diet==g)
    if(nrow(temp)==0){next}
    freq1 = data.frame(table(temp$Month))
    freq1$sex = sex
    freq1$Diet = g
    freq = rbind(freq, freq1)
  }
}
freq
# mousedat = mousedat[colnames(mousedat)[grepl('b23', colnames(mousedat))]]
mousedat_log2 = log2(mousedat[2:ncol(mousedat)] )

unique(mouseinfo$Group) #  "O_MA" "M_FE" "M_MA"
yma = mouseinfo[mouseinfo$Group == "Y_MA", ]
mfe = mouseinfo[mouseinfo$Group == "M_FE", ]
mma = mouseinfo[mouseinfo$Group == "M_MA", ]
oma = mouseinfo[mouseinfo$Group == "O_MA", ]

# anova_out = anova_func(mousedat[rownames(yma)], yma$Diet)
# anova_out = anova_func(mousedat[rownames(mfe)], mfe$Diet)
# anova_out = anova_func(mousedat[rownames(mma)], mma$Diet)
anova_out = anova_func(mousedat[rownames(oma)], oma$Diet)
p = nrow(subset(anova_out, Pvalue<0.05 ))
fdr= nrow(subset(anova_out, FDR<0.05 ))
anova_out$uni = as.character(lapply(rownames(anova_out), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[1]
}))
anova_out$gene = as.character(lapply(rownames(anova_out), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  gsub(';$', '', x[2])
}))

# human <- useMart('ensembl', dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# mouse <- useMart('ensembl', dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# m2h.g <- getLDS(values = anova_out$gene, # will be changed
#                 attributes = c("mgi_symbol"), # from type
#                 attributesL = c("hgnc_symbol","chromosome_name","start_position"), # to type
#                 mart = mouse, # from base
#                 martL = human, # to base
#                 filters = "mgi_symbol", uniqueRows = T)
# m2h.g = data.frame(m2h.g)
# m2h.g = m2h.g %>% distinct( MGI.symbol, .keep_all = TRUE)
# rownames(m2h.g) = m2h.g$MGI.symbol

anova_out$human_gene = m2h.g[anova_out$gene, 'HGNC.symbol']
write.xlsx(anova_out, paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_oma_p",p,'_fdr',fdr,"_anova.xlsx"), rowNames = T)

files = list.files("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/")
files = files[grepl('_anova.xlsx',files) & !grepl('log2',files)]
files
temp = subset(coefficients1, p_adjust<0.05 & !is.na( human_gene))
# vennlist = list()

df1 = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_yma_anova.xlsx") , rowNames = T)
df2 = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_oma_p767_fdr33_anova.xlsx") , rowNames = T)

df1 = subset(df1, FDR<0.05 & !is.na( human_gene))
df2 = subset(df2, FDR<0.05 & !is.na( human_gene))
intersect(df1$human_gene, df2$human_gene)

pros = c()
label = c()
f = files[4]#[5] #
anova_out = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/", f) , rowNames = T)
anova_out = subset(anova_out, FDR<0.05 & !is.na( gene))
n = intersect(temp[temp$qt_coef>0 , 'mouse_gene'], anova_out$gene)
pros = c(pros, n)
label = c(label, rep('up2', length(n)))
n = intersect(temp[temp$qt_coef<0, 'mouse_gene'], anova_out$gene)
pros = c(pros, n)
label = c(label, rep('dn2', length(n)))
pros_info = data.frame( label, mousegene = pros, row.names = pros)
pros_info

heatmapdat = mousedat[mouseinfo$Sample_name]
heatmapdat$gene = as.character(lapply(rownames(heatmapdat), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  gsub(';$', '', x[2])
}))
heatmapdat = heatmapdat[heatmapdat$gene %in% pros_info$mousegene,]
human <- useMart('ensembl', dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl', dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
m2h.g <- getLDS(values = heatmapdat$gene, # will be changed
                attributes = c("mgi_symbol"), # from type
                attributesL = c("hgnc_symbol","chromosome_name","start_position"), # to type
                mart = mouse, # from base
                martL = human, # to base
                filters = "mgi_symbol", uniqueRows = T)
m2h.g = data.frame(m2h.g)
m2h.g = m2h.g[m2h.g$MGI.symbol %in% heatmapdat$gene,]

m2h.g = m2h.g %>% distinct( MGI.symbol, .keep_all = TRUE)
rownames(m2h.g) = m2h.g$MGI.symbol
heatmapdat$humangene = m2h.g[heatmapdat$gene,  'HGNC.symbol']
rownames(heatmapdat) = heatmapdat$gene
heatmapdat = heatmapdat[pros_info$mousegene, ]

pros_info$humangene = m2h.g[pros_info$mousegene, 'HGNC.symbol']

heatmapInfo1 = data.frame(Group = mouseinfo$Group, sample = mouseinfo$Sample_name)
heatmapInfo1$Diet = mouseinfo$Diet

heatmapInfo1$Diet = factor(heatmapInfo1$Diet, levels = c("L", "M", "H" ))
heatmapInfo1 <- heatmapInfo1[order(heatmapInfo1$Diet), ]

heatmapInfo1$Group = factor(heatmapInfo1$Group, levels = c("Y_MA", "M_MA", "M_FE", "O_MA" ))
heatmapInfo1 <- heatmapInfo1[order(heatmapInfo1$Group), ]
heatmapInfo1 = data.frame(Group = heatmapInfo1$Group, Diet = heatmapInfo1$Diet, row.names = heatmapInfo1$sample)

pros_info1 = data.frame(label = pros_info$label, row.names = pros_info$mousegene)

pros_info1$label = gsub('1', 'yma', pros_info1$label)
pros_info1$label = gsub('2', 'oma', pros_info1$label)

ComplexHeatmap::pheatmap(heatmapdat[rownames(pros_info1) , rownames(heatmapInfo1)],  scale = 'row', #rownames(heatmapInfo)
         # color = colorRampPalette(c("blue", "white","red" ))(1000),
         color = colorRampPalette(c("#00B0F0", "white","#E74809" ))(1000), 
        breaks = c(seq(-3, 3, length=100)),
         annotation_col = heatmapInfo1, # 样本信息
         annotation_row = pros_info1,
         row_split = pros_info1$label, # gaps_row = c( 14),
         # annotation_colors = ann_colors , 
         labels_row = pros_info$humangene, 
         #angle_col = 45,
         cluster_cols = F, 
         cluster_rows = F, 
         #clustering_method = 'median',
         cellwidth = 4,
         cellheight = 4,#2, # 
         fontsize_col = 1.5,
         fontsize_row = 2,
         show_rownames = T,
         show_colnames = T, name = ' ',
         main = ' ' )


for(f in files[2:5]){
  f = files[5]
  vennlist = list()
  anova_out = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/", f) , rowNames = T)
  anova_out = subset(anova_out, FDR<0.05 & !is.na( human_gene))
  vennlist[[f]] = anova_out$human_gene
  # n = intersect(temp[temp$qt_coef>0 , 'human_gene'], anova_out$human_gene)
  # if(length(n)>0) vennlist[[paste0( 'up')]] = n
  # n = intersect(temp[temp$qt_coef<0, 'human_gene'], anova_out$human_gene)
  # if(length(n)>0) vennlist[[paste0( 'dn')]] = n
  vennlist[[paste0( 'up')]] = temp[temp$qt_coef>0 , 'human_gene']
  vennlist[[paste0( 'dn')]] = temp[temp$qt_coef<0 , 'human_gene']

}

temp = data.frame(x = c('Y-MA', 'M-MA', 'Y-FE', 'O-MA' ), y = c(252, 1,0,33))
temp$x = factor(temp$x, levels = c('Y-MA', 'M-MA', 'Y-FE', 'O-MA' ))
ggplot(temp, aes(x, y, fill = x))+
  geom_bar(stat = 'identity', width  =0.5)+
  geom_text(aes(label= y),vjust = -0.5)+
  theme_classic()+ labs(x = '', y = '')+
  theme(text = element_text(size = 16))+
  scale_fill_jama()#_aaas()# igv()

# narow = rowSums(is.na(mousedat_log2))
# mousedat_log2 = mousedat_log2[narow!=ncol(mousedat_log2),]
# 
# mouse_log2_min0.8 = mousedat_log2
# mouse_log2_min0.8[is.na(mouse_log2_min0.8)] = min(mouse_log2_min0.8, na.rm = T) * 0.8
# mouse_log2_min0.8[1:20,1:2]
# mouse_log2_min0.8 = data.frame(t(mouse_log2_min0.8))
# 
# dim(mouse_log2_min0.8)
# mouse_log2_min0.8$Age = as.numeric(mouseinfo[rownames(mouse_log2_min0.8), 'Month'])
# dim(mouse_log2_min0.8)
# mouse_log2_min0.8 = cbind(mouse_log2_min0.8['Age'], mouse_log2_min0.8[,1:6685])
# mouse_log2_min0.8[1:10,1:2]
# mouse_log2_min0.8_zscore = mouse_log2_min0.8
# mouse_log2_min0.8_zscore[,2:ncol(mouse_log2_min0.8_zscore)] = apply(mouse_log2_min0.8_zscore[,2:ncol(mouse_log2_min0.8_zscore)], 2, function(x){
#   (x-mean(x, na.rm = T))/sd(x, na.rm = T)
# })
# dim(mouse_log2_min0.8_zscore)

###### mouse ML model ######
# mouseml = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/mouse_ML/mouse_zscore_xgb_shap_predict20250214.xlsx", rowNames = T)
mouseml = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/mouse_ML/mouse_zscore_rf_shap_predict20250214.xlsx", rowNames = T)

mouseml = cbind(mouseml, mouseinfo[rownames(mouseml), c("Month","Group","Diet","Sex")])
mouseml$Month
colnames(mouseinfo)
mouseinfo[rownames(mouseml), 'Month']
mouseinfo[rownames(mouseml), c("Month","Group","Diet","Sex" )]
mouseml = mouseml[mouseml$Sex == 'Male' & mouseml$Diet == 'M',]
dim(mouseml)

2:19
for(i in 2:19){
  spcor = cor(mouseml[c(1,i)], method = 'spearman')
  print(c(colnames(mouseml)[i], round(spcor[1,2], 3)))
}


mouseml = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/mouse_ML/mouse_zscore_xgb_shap_predict20250214.xlsx", rowNames = T)
# mouseml = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/mouse_ML/mouse_zscore_rf_shap_predict20250214.xlsx", rowNames = T)
mouseml = cbind(mouseml, mouseinfo[rownames(mouseml), c("Month","Group","Diet","Sex")])
for(n in c(15, 40)){# seq(10,190,5)
  pred = paste0('feat', n)# 55 70 75
  mouseml$Age_Gap = mouseml[,pred] - mouseml$Month
  temp = mouseml[c('Sex', 'Diet', 'Age_Gap')]
  #head(temp)
  temp$Diet =  factor(temp$Diet, levels = c('L', 'M','H'))
  p = ggplot(temp[temp$Sex == 'Male', ], aes(Diet, Age_Gap, color = Diet))+
    geom_boxplot(lwd = 1, width = 0.8)+
    #geom_jitter(alpha = 0.7)+
    geom_point( position = position_jitterdodge(), alpha = 0.5)+
    theme_classic() +
    #geom_smooth(aes(  fill = variable), se = T, method = 'loess', alpha = 0.1)+ #  alpha = 0.15
    labs(title = pred )+
    
    scale_color_d3( )+
    #scale_fill_d3(palette ="category20")+
    theme(axis.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 12),
          text = element_text(size = 15))
  print(p)
}

pred = paste0('feat', n)# 55 70 75
pred = 'feat40'
mouseml$Age_Gap = mouseml[,pred] - mouseml$Month
temp = mouseml[c( 'Group', 'Sex', 'Diet', 'Age_Gap')]
head(temp)
temp$Diet =  factor(temp$Diet, levels = c('L', 'M','H'))
temp$Group = factor(temp$Group, levels = c("Y_MA" , "M_MA", "O_MA","M_FE"))
ggplot(data = temp[temp$Sex == 'Female', ] ,aes(Group , Age_Gap, color = Diet) )+
  geom_boxplot(  lwd = 1, width = 0.5)+
  labs(title = pred )+
  theme_classic() +
  scale_color_d3()+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))
  

for (j in c("Y_MA" , "M_MA", "O_MA")){
  L = temp[temp$Group == j & temp$Diet == 'L' & temp$Sex == 'Male', 'Age_Gap']
  M = temp[temp$Group == j & temp$Diet == 'M'& temp$Sex == 'Male', 'Age_Gap']
  H = temp[temp$Group == j & temp$Diet == 'H'& temp$Sex == 'Male', 'Age_Gap']
  ml = wilcox.test(L, M)
  hl= wilcox.test(L, H)
  mh= wilcox.test(H, M)
  print(c(j, ml$p.value, hl$p.value, mh$p.value))
}
j = "M_FE"
L = temp[temp$Group == j & temp$Diet == 'L' & temp$Sex == 'Female', 'Age_Gap']
M = temp[temp$Group == j & temp$Diet == 'M'& temp$Sex == 'Female', 'Age_Gap']
H = temp[temp$Group == j & temp$Diet == 'H'& temp$Sex == 'Female', 'Age_Gap']
ml = wilcox.test(L, M)
hl= wilcox.test(L, H)
mh= wilcox.test(H, M)
print(c(j, ml$p.value, hl$p.value, mh$p.value))



mousedat = read.csv("Z:/members/wangyingrui/TTTD/Mouse/PDM_combineMatrix_20230609.csv", row.names = 1 )
rownames(mousedat) = gsub(';$', '', mousedat$X)
mousedat[1:3,1:3]
mouseinfo = read.csv("Z:/members/wangyingrui/TTTD/Mouse/20231123_AgingPDM_Thyroid_sample_information.CSV" )
rownames(mouseinfo) = mouseinfo$Sample_name
mouseinfo = mouseinfo[mouseinfo$Sample_name%in% colnames(mousedat), ]
# mousedat = mousedat[colnames(mousedat)[grepl('b23', colnames(mousedat))]]
# mousedat_log2 = log2(mousedat[2:ncol(mousedat)] )

unique(mouseinfo$Group) #  "O_MA" "M_FE" "M_MA"
ss = mouseinfo[mouseinfo$Group %in% c("Y_MA", "M_MA", "O_MA")  & mouseinfo$Diet == 'M', ]

anova_out = anova_func(mousedat[rownames(ss)], ss$Group)
p = nrow(subset(anova_out, Pvalue<0.05 ))
fdr= nrow(subset(anova_out, FDR<0.05 ))
anova_out$uni = as.character(lapply(rownames(anova_out), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  x[1]
}))
anova_out$gene = as.character(lapply(rownames(anova_out), function(x){
  x = strsplit(x, '_', fixed = T)[[1]]
  gsub(';$', '', x[2])
}))

m2h.g <- getLDS(values = anova_out$gene, # will be changed
                attributes = c("mgi_symbol"), # from type
                attributesL = c("hgnc_symbol","chromosome_name","start_position"), # to type
                mart = mouse, # from base
                martL = human, # to base
                filters = "mgi_symbol", uniqueRows = T)
m2h.g = data.frame(m2h.g)
m2h.g = m2h.g %>% distinct( MGI.symbol, .keep_all = TRUE)
rownames(m2h.g) = m2h.g$MGI.symbol

anova_out$human_gene = m2h.g[anova_out$gene, 'HGNC.symbol']
write.xlsx(anova_out, paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_DietM_p",p,'_fdr',fdr,"_anova.xlsx"), rowNames = T)

mfuzzinput = data.frame(
Y = rowMeans(subset(anova_out, Pvalue<0.05, select = rownames(ss[ss$Group=="Y_MA", ]) ), na.rm = T),
M = rowMeans(subset(anova_out, Pvalue<0.05, select = rownames(ss[ss$Group=="M_MA", ]) ), na.rm = T),
O = rowMeans(subset(anova_out, Pvalue<0.05, select = rownames(ss[ss$Group=="O_MA", ]) ), na.rm = T)
)


# library(Biobase)
# library(BiocGenerics)
# library(parallel)
library(RColorBrewer)
library(Mfuzz)
mycol <- c("cyan","yellow","orangered")
mycolor <- colorRampPalette(mycol)(100)

rownames(mfuzzinput) = rownames(subset(anova_out, Pvalue<0.05 ))
mat <- as.matrix(mfuzzinput)

#创建用于Mfuzz的对象；
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
#Features Samples
#9345 4
#过滤缺失值超过25%的基因;
dt.r <- filter.NA(dt, thres=0.25)
#以均值的方式填充缺失值；
dt.f <- fill.NA(dt.r, mode="mean")
#过滤低表达或表达量变化不大的基因；
#由于是差异基因，这里不做过滤；
tmp <- filter.std(dt.f,min.std=0)

#对数据进行标准化；
dt.s <- standardise(tmp)
#查看标准化后的数据；
df.s <- dt.s@assayData$exprs
head(df.s)

#模糊聚类
m1 <- mestimate(dt.s)#使用mestimate函数估计m值；
m1
set.seed(007)
cl <- mfuzz(dt.s,c=9,m=m1)
mfuzz.plot(dt.s, cl,
           mfrow=c(3,3),
           new.window= FALSE,
           time.labels=colnames(dt.s),
           colo = mycolor)

# 查看每个蛋白所属cluster
protein_cluster <- data.frame(cl$cluster)
table(protein_cluster$cl.cluster)
protein_cluster <- cbind(mat[rownames(protein_cluster), ], protein_cluster)
protein_cluster = data.frame(protein_cluster)
write.xlsx(protein_cluster, paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_DietM_mfuzzCluster.xlsx"), rowNames = T)
cluste4 = rownames(protein_cluster[protein_cluster$cl.cluster == 4,])
cluste9 = rownames(protein_cluster[protein_cluster$cl.cluster == 9,])

cluste4 = anova_out[cluste4, ss$Sample_name]
cluste9 = anova_out[cluste9, ss$Sample_name]

ss = mouseinfo[mouseinfo$Group %in% c("Y_MA", "M_MA", "O_MA") , ]

dat = mousedat[rownames(protein_cluster[protein_cluster$cl.cluster == 4,]), ss$Sample_name ]
#dat = dat[intersect( mouseinfo[mouseinfo$Group %in% c("Y_MA", "O_MA", "M_MA"), 'Sample_name'], colnames(dat))]
dat$proteins = rownames(dat)
dat = melt(dat)
dat$Group = ss[dat$variable,  'Group']
dat$Diet = ss[dat$variable,  'Diet']
head(dat)
unique(dat$Group)
unique(dat$Diet)
dat$Group = factor(dat$Group, levels = c("Y_MA" , "M_MA", "O_MA"))
dat$Diet = factor(dat$Diet, levels = c( "L", "M", "H"))

ggplot(dat, aes(Group, log2(value), color = Group))+
  geom_boxplot(lwd = 1, width = 0.8)+
  geom_jitter(alpha = 0.7)+
  #geom_smooth(aes(  fill = variable), se = T, method = 'loess', alpha = 0.1)+ #  alpha = 0.15
  labs(y = 'Log2 Expressions' )+
  theme_classic()+
  scale_color_d3( )+
  #scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))+
  scale_y_continuous(n.breaks = 8)+
  stat_summary(fun.data = function(x) data.frame(y=13 , label = paste("Median=", round(median(x), 3))), geom="text", size = 4)

head(dat)
dat %>% na.omit() %>% 
  group_by(Group, Diet) %>% 
  summarise(mean_value=mean(log2(value), na.rm = T), med_value=median(log2(value), na.rm = T), sd_value=sd(log2(value), na.rm = T)) -> dat1
dat1

dat1 = data.frame(dat1)
dat1$Group1 = paste0(dat1$Group, dat1$Diet)

dat$Group1 = paste0(dat$Group, dat$Diet)
dat$Group1 = factor(dat$Group1, levels = c( "Y_MAL", "Y_MAM", "Y_MAH","M_MAL",  "M_MAM", "M_MAH", "O_MAL", "O_MAM", "O_MAH"))

unique(dat$Group) # Y_MA O_MA M_MA
ggplot()+
  geom_boxplot( aes(Group1, log2(value), color = Diet), data = dat , lwd = 1, width = 0.5)+
  geom_line(aes(Group1, med_value, group = Diet, color = Diet) ,data = dat1, size = 1.1 )+
  geom_ribbon(aes(ymin=med_value-sd_value, ymax=med_value+sd_value, x=Group1, group = Diet , fill = Diet), data = dat1, alpha = 0.1)+
  stat_summary( aes(Group1, log2(value), color = Diet), data = dat ,fun.data = function(x) data.frame(y=min(x) , label = paste("Median=", round(median(x), 3))), geom="text", size = 4)+
  labs(y = 'Log2 Expressions' )+scale_color_d3( )+scale_fill_d3( )+
  theme_classic()+scale_y_continuous(n.breaks = 8)+
  #scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15)) +
  geom_signif(aes(Group1, log2(value), color = Diet), data = dat ,comparisons= list(c('M_MAL', 'M_MAM'), c('M_MAL', 'M_MAH'), c('M_MAM', 'M_MAH'),
                                c('O_MAL', 'O_MAM'), c('O_MAL', 'O_MAH'), c('O_MAM', 'O_MAH'),
                                c('Y_MAL', 'Y_MAM'), c('Y_MAL', 'Y_MAH'), c('Y_MAM', 'Y_MAH') ), #vecs 
              step_increase = 0.3,
              test="t.test",  # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
              map_signif_level=T   # 标签样式F为数字，T为*号
  )
vecs = list()
for (j in unique(dat$Group)){
  vec = combn(unique(dat[dat$Group == j, 'Group1']), 2)
  for(i in seq(1,length(vec), 2)){
    # print((i+ 1) /2)
    vecs[[(i+ 1) /2]] = vec[i:(i+1)]
  }
}

for (j in unique(dat$Group)){
  L = dat[dat$Group == j & dat$Diet == 'L', 'value']
  M = dat[dat$Group == j & dat$Diet == 'M', 'value']
  H = dat[dat$Group == j & dat$Diet == 'H', 'value']
  ml = wilcox.test(L, M)
  hl= wilcox.test(L, H)
  mh= wilcox.test(H, M)
  print(c(j, ml$p.value, hl$p.value, mh$p.value))
}

### exp intensity####
# "A1"
mousedat = read.csv("Z:/members/wangyingrui/TTTD/Mouse/Aging_combineMatrix_20230609.csv", row.names = 1 )
mousedat[1:3,1:3]
rownames(mousedat) = gsub(';$', '', rownames(mousedat))
mousedat[1:3,1:3]
mouseinfo = read.csv("Z:/members/wangyingrui/TTTD/Mouse/20231123_AgingPDM_Thyroid_sample_information.CSV" )
rownames(mouseinfo) = mouseinfo$Sample_name
mouseinfo = mouseinfo[mouseinfo$Sample_name%in% colnames(mousedat), ]

ss = mouseinfo[mouseinfo$Group %in% c("A1") , ]
ss

mousepro = rownames(coefficients1[coefficients1$Gene %in% c('Tkt', 'Acly', 'Dtnbp1', 'Fasn', 'Pdg', 'Thrsp', 'Slc25a1'),])

dat = data.frame( mousedat[mousepro , ss$Sample_name ] )  
dat$proteins = mousepro
dat = melt(dat)
dat$Age = ss[dat$variable, 'Month']

ggplot(dat, aes(Age,  value ))+
  geom_smooth( alpha = 0.15, method = 'loess' )+ # se = F
  geom_point( )+theme_classic()+labs(y = 'Log2 Expressions', title = unique(dat$proteins) ) 
for(pro in mousepro){
  p = ggplot(dat[dat$proteins==pro, ], aes(Age, log2(value) ))+
    geom_smooth( alpha = 0.15, method = 'loess' )+ # se = F
    geom_point()+theme_classic()+
    labs(y = 'Log2 Expressions', title = pro )+
    scale_color_d3( )+scale_fill_d3(palette ="category20")+
    theme(axis.text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 12),
          text = element_text(size = 15)) +
    scale_y_continuous(n.breaks = 5)
  ggsave(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/mouse_A1_log2_", pro, '.pdf'), p)
}


# "Y_MA", "M_MA", "O_MA", "M_FE"
mousedat = read.csv("Z:/members/wangyingrui/TTTD/Mouse/PDM_combineMatrix_20230609.csv", row.names = 1 )
rownames(mousedat) = gsub(';$', '', mousedat$X)
mousedat[1:3,1:3]
mouseinfo = read.csv("Z:/members/wangyingrui/TTTD/Mouse/20231123_AgingPDM_Thyroid_sample_information.CSV" )
rownames(mouseinfo) = mouseinfo$Sample_name
mouseinfo = mouseinfo[mouseinfo$Sample_name%in% colnames(mousedat), ]


mousepro = c('P40142_Tkt') # 'Q62264_Thrsp'
ss = mouseinfo[mouseinfo$Group %in% c("Y_MA", "M_MA", "O_MA", "M_FE") , ]

dat = data.frame(value = as.numeric(mousedat['P40142_Tkt', ss$Sample_name ]), Age = ss$Month)
#dat = dat[intersect( mouseinfo[mouseinfo$Group %in% c("Y_MA", "O_MA", "M_MA"), 'Sample_name'], colnames(dat))]
dat$proteins = 'P40142_Tkt'
#dat = melt(dat)
dat$Group = ss$Group
dat$Diet = ss$Diet
head(dat)
unique(dat$Group)
unique(dat$Diet)
dat$Group = factor(dat$Group, levels = c("Y_MA" , "M_MA", 'M_FE', "O_MA"))
dat$Diet = factor(dat$Diet, levels = c( "L", "M", "H"))

head(dat)

ggplot(dat, aes(Age,  value ))+
  geom_smooth( alpha = 0.15, method = 'loess' )+ # se = F
  geom_point(aes( color = Group))

ggplot(dat, aes(Age, log2(value), color = Group))+
  geom_smooth(aes(fill = Group), alpha = 0.15, method = 'lm' )+ # se = F
  geom_point()+labs(y = 'Log2 Expressions', title = unique(dat$proteins) )+
  theme_classic()+
  scale_color_d3( )+scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))+
  scale_y_continuous(n.breaks = 8)

ggplot(dat, aes(Group, log2(value), color = Group))+
  geom_boxplot(lwd = 1, width = 0.8)+
  geom_jitter(alpha = 0.7)+
  #geom_smooth(aes(  fill = variable), se = T, method = 'loess', alpha = 0.1)+ #  alpha = 0.15
  labs(y = 'Log2 Expressions', title = unique(dat$proteins) )+
  theme_classic()+
  scale_color_d3( )+
  #scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))+
  scale_y_continuous(n.breaks = 8)+
  stat_summary(fun.data = function(x) data.frame(y=24 , label = paste("Median=", round(median(x), 3))), geom="text", size = 4)

head(dat)
dat %>% na.omit() %>% 
  group_by(Group, Diet) %>% 
  summarise(mean_value=mean(log2(value), na.rm = T), med_value=median(log2(value), na.rm = T), sd_value=sd(log2(value), na.rm = T)) -> dat1
dat1

dat1 = data.frame(dat1)
dat1$Group1 = paste0(dat1$Group, dat1$Diet)

dat$Group1 = paste0(dat$Group, dat$Diet)
unique(dat$Group1)
dat$Group1 = factor(dat$Group1, levels = c( "Y_MAL", "Y_MAM", "Y_MAH","M_MAL",  "M_MAM", "M_MAH",  "M_FEL", "M_FEM", "M_FEH", "O_MAL", "O_MAM", "O_MAH"))

unique(dat$Group) # Y_MA O_MA M_MA
head(dat)
ggplot()+
  geom_boxplot( aes(Group1, log2(value), color = Diet), data = dat , lwd = 1, width = 0.5)+
  geom_line(aes(Group1, med_value, group = Diet, color = Diet) ,data = dat1, size = 1.1 )+
  # geom_line(aes(Group1, med_value, group = Group, color = Group) ,data = dat1, size = 1.1 )+
  geom_ribbon(aes(ymin=med_value-sd_value, ymax=med_value+sd_value, x=Group1, group = Diet , fill = Diet), data = dat1, alpha = 0.1)+
  # geom_ribbon(aes(ymin=med_value-sd_value, ymax=med_value+sd_value, x=Group1, group = Group , fill = Group), data = dat1, alpha = 0.1)+
  #stat_summary( aes(Group1, log2(value), color = Diet), data = dat ,fun.data = function(x) data.frame(y=min(x) , label = paste("Median=", round(median(x), 3))), geom="text", size = 4)+
  labs(y = 'Log2 Expressions', x= '', title = unique(dat$proteins) )+
  #scale_color_d3( )+
  scale_color_manual(values = c('#1f77b4', '#ff7f0e',  '#d62728' ))+
  scale_fill_d3( )+
  theme_classic()+scale_y_continuous(n.breaks = 8)+
  #scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15),
        axis.text.x  = element_text(size = 10)) #+
  annotate('text', x = 'M_MAL', y = 24.2, label = 'anova p *')+
  annotate('text', x = 'O_MAL', y = 24.2, label = 'anova p *****')+
  #annotate('text', x = 'M_FEL', y = 24, label = '*')+
  annotate('text', x = 'Y_MAL', y = 24.2, label = 'anova p *')
  # geom_signif(comparisons = list(c('M_MAL', 'M_MAM' , 'M_MAH'),
  #                                c('O_MAL', 'O_MAM', 'O_MAH'),
  #                                c('M_FEL', 'M_FEM' , 'M_FEH'),
  #                                c('Y_MAL','Y_MAM', 'Y_MAH')),
  #            annotations = c('*', "*****", ' ', "*" ))
  # stat_compare_means(aes(Group1, log2(value), color = Group), data = dat, method = "anova")
  #stat_compare_means(aes(Group1, log2(value), color = Diet, group = Group), method = "anova", data = dat )
  geom_signif(aes(Group1, log2(value), color = Group), data = dat ,comparisons= list(c('M_MAL', 'M_MAM' , 'M_MAH'),
                                                                                    c('O_MAL', 'O_MAM', 'O_MAH'),
                                                                                    c('M_FEL', 'M_FEM' , 'M_FEH'),
                                                                                     c('Y_MAL','Y_MAM', 'Y_MAH') ), #vecs 
              step_increase = 0.3,
              test="aov()",  # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
              map_signif_level=T   # 标签样式F为数字，T为*号
  )

  
for(g in unique(dat$Group))  {
  anova_result <- aov( value  ~ Diet , data = dat[dat$Group == g, ] )
  y = summary(anova_result)
  print(c(g, y[[1]][,5][1]))
}

library(ggpubr)
table(dat$Group1)
head(dat)
vecs = list()
for (j in unique(dat$Group1)){
  vec = combn(unique(dat[dat$Group == j, 'Group1']), 2)
  for(i in seq(1,length(vec), 2)){
    # print((i+ 1) /2)
    vecs[[(i+ 1) /2]] = vec[i:(i+1)]
  }
}

####### 5 细胞药物验证 20241203 ########
ttest = function(exp, a, b, c, d){
  Pvalue = c()
  FC = c()
  logFC = c()
  for(i in 1:nrow(exp)){
    # un-impute
    y = try(t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ])),silent=FALSE)
    if('try-error' %in% class(y))
    {
      Pvalue = c(Pvalue, NA)
    }else{
      y = t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]))
      Pvalue = c(Pvalue, y$p.value)
    }
  }
  pre = rowMeans( exp[ , a:b] , na.rm = TRUE)
  post = rowMeans( exp[ , c:d ] , na.rm = TRUE)
  FC = pre/post
  logFC = log2(pre/post)
  differ = pre-post
  
  FDR=p.adjust(Pvalue, "BH")
  out<-cbind(exp, Pvalue, FDR, FC, logFC, differ)
  out
}
ttest1_log2 = function(exp, a, b, c, d){
  Pvalue = c()
  FC = c()
  logFC = c()
  for(i in 1:nrow(exp)){
    # un-impute
    y = try(t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ])),silent=FALSE)
    if('try-error' %in% class(y))
    {
      Pvalue = c(Pvalue, NA)
    }else{
      y = t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]))
      Pvalue = c(Pvalue, y$p.value)
    }
  }
  pre = rowMeans( 2**exp[ , a:b] , na.rm = TRUE)
  post = rowMeans( 2**exp[ , c:d ] , na.rm = TRUE)
  FC = pre/post
  logFC = log2(pre/post)
  differ = pre-post
  
  FDR=p.adjust(Pvalue, "BH")
  out<-cbind(exp, Pvalue, FDR, FC, logFC, differ)
  out
}
ttest_paired = function(exp, a, b, c, d){
  Pvalue = c()
  FC = c()
  logFC = c()
  for(i in 1:nrow(exp)){
    # un-impute
    y = try(t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]), paired = T),silent=FALSE)
    if('try-error' %in% class(y))
    {
      Pvalue = c(Pvalue, NA)
    }else{
      y = t.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]), paired = T)
      Pvalue = c(Pvalue, y$p.value)
    }
  }
  pre = rowMeans( exp[ , a:b] , na.rm = TRUE)
  post = rowMeans( exp[ , c:d ] , na.rm = TRUE)
  FC = pre/post
  logFC = log2(pre/post)
  differ = pre-post
  
  FDR=p.adjust(Pvalue, "BH")
  out<-cbind(exp, Pvalue, FDR, FC, logFC, differ)
  out
}
wilcoxtest_pair = function(exp, a, b, c, d){
  Pvalue = c()
  FC = c()
  logFC = c()
  for(i in 1:nrow(exp)){
    # un-impute
    y = try(wilcox.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]), paired = T),silent=FALSE)
    if('try-error' %in% class(y))
    {
      Pvalue = c(Pvalue, NA)
    }else{
      y = wilcox.test(as.numeric(exp[i, a:b]), as.numeric(exp[i, c:d ]), paired = T)
      Pvalue = c(Pvalue, y$p.value)
    }
  }
  pre = rowMeans( exp[ , a:b] , na.rm = TRUE)
  post = rowMeans( exp[ , c:d ] , na.rm = TRUE)
  FC = pre/post
  logFC = log2(pre/post)
  #differ = pre-post
  
  FDR=p.adjust(Pvalue, "BH")
  out<-cbind(exp, Pvalue, FDR, FC, logFC)#, differ
  out
}

setwd("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5")
list.files("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5")
ptv3_Hsamp = read.csv("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/20241127_ptv3_Hsamp_CAL-62_224samp_9142prot_matrix_quantile.csv")
ptv3_Hsamp_sampinfo= read.csv("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/20241127_ptv3_Hsamp_CAL-62_224samp_9142prot_sampinfo.csv")
All_DEPs = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/All_DEPs_20241204.xlsx")
dim(ptv3_Hsamp)
ptv3_Hsamp[1:3,1:6]
rownames(ptv3_Hsamp) = ptv3_Hsamp$samp_ID
ptv3_Hsamp["H4950",1:6]
ncol(ptv3_Hsamp)
ptv3_Hsamp = ptv3_Hsamp[1:9144]
ptv3_Hsamp_ratio = ptv3_Hsamp
for(nr in 1:nrow(ptv3_Hsamp)){
  sid = rownames(ptv3_Hsamp)[nr]
  temp = ptv3_Hsamp_sampinfo[ptv3_Hsamp_sampinfo$samp_ID == sid, ]
  if (temp$pert_id=='control')next
  
  controlid = temp$control
  ptv3_Hsamp_ratio[nr, 3:9144] = ptv3_Hsamp_ratio[nr, 3:9144]/ptv3_Hsamp_ratio[controlid, 3:9144]
}
ptv3_Hsamp_ratio[1:3,1:6]
ptv3_Hsamp_sampinfo1 = ptv3_Hsamp_sampinfo[c(18, 23, 24, 31  )]
unique(ptv3_Hsamp_sampinfo1$PRISM1st_label_total)

nonRes = subset(ptv3_Hsamp_sampinfo, PRISM1st_label_total == "non-responsive" & samp_ID %in% ptv3_Hsamp$samp_ID)
sensitive = subset(ptv3_Hsamp_sampinfo, PRISM1st_label_total == "sensitive"& samp_ID %in% ptv3_Hsamp$samp_ID)

rownames(ptv3_Hsamp_sampinfo) = ptv3_Hsamp_sampinfo$samp_ID
intersect(colnames(ptv3_Hsamp_ratio), top_features_$top_features_.1.95.)
colnames(ptv3_Hsamp_ratio)[grepl('.', colnames(ptv3_Hsamp_ratio), fixed = T)]
write.xlsx(ptv3_Hsamp_ratio[c(nonRes$samp_ID, sensitive$samp_ID), intersect(colnames(ptv3_Hsamp_ratio), top_features_$top_features_.1.95.)], "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/PTV_95feat.xlsx", rowNames = T)


PTV3_pred_95feat = read.xlsx("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/PTV_95feat.xlsx_pred.xlsx", rowNames = T)

t.test(PTV3_pred_95feat[nonRes$samp_ID, 'Pred'], PTV3_pred_95feat[sensitive$samp_ID, 'Pred'])
wilcox.test(PTV3_pred_95feat[nonRes$samp_ID, 'Pred'], PTV3_pred_95feat[sensitive$samp_ID, 'Pred'])

dat = data.frame(t(ptv3_Hsamp_ratio[c(sensitive$samp_ID, nonRes$samp_ID), 3:ncol(ptv3_Hsamp_ratio)]))
dat[] = lapply(dat, as.numeric)
ttestout = ttest(dat, 1, nrow(sensitive), 1+nrow(sensitive), ncol(dat))
# dim(ttestout[ttestout$FDR<0.05 & abs(ttestout$logFC)>log2(1.5),])
ttestout1 = subset(ttestout, FDR<0.05 & abs( logFC)>log2(1.2) )
ttestout[1:3,1:3]
ttestout$pre = rowMeans( ttestout[sensitive$samp_ID] , na.rm = TRUE)
ttestout$post = rowMeans( ttestout[nonRes$samp_ID] , na.rm = TRUE)
ttestout$sub = ttestout$pre-ttestout$post
ttestout$symbol = as.character(lapply(rownames(ttestout), function(x){
  strsplit(x, '_')[[1]][2]
}))
write.xlsx(ttestout, "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptv3_sensitive-nonres_ttest.xlsx", rowNames = T)
ttestout$sig
ttestout$Sig = ifelse(is.na(ttestout$FDR), "None", ifelse(ttestout$logFC> log2(1.2) & ttestout$FDR<0.05, "Up",
                                                          ifelse( ttestout$logFC< -log2(1.2) & ttestout$FDR<0.05, "Down", "None")))
table(ttestout$Sig)
ttestout1 = subset(ttestout, FDR<0.05 & abs( logFC)>log2(1.2) )
ggplot(ttestout[abs(ttestout$logFC) <2 , ],
       aes(x = logFC, y = -log10(FDR), color = Sig )) +#  , size = -log10(qt_padjust)
  geom_point(alpha=0.7, size=1 )+
  scale_color_manual(values=c( "#546de5", "#d2dae2", "#ff4757"))+# , "#d2dae2"
  geom_text_repel( data = ttestout1,
                   aes(x=logFC, y = -log10(FDR) ,label=symbol),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
                   segment.color = "black", show.legend = FALSE )+
  labs(x="log2(Fold change)", y="-log10(p-adjust)", title = '' )+
  theme_classic()+ scale_y_continuous(n.breaks = 10) + #ylim(0, 60)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = log2(1.2), linetype = "dashed")+
  geom_vline(xintercept = -log2(1.2), linetype = "dashed")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSVA)
CP <- msigdbr(species = "Homo sapiens", category = "C2" ) #%>% dplyr::select(gs_subcat, gs_name, entrez_gene, gene_symbol, human_gene_symbol )
CP = subset(CP, gs_subcat!="CGP")
CP <- msigdbr(species = "Homo sapiens", category = "H" )
head(CP)
unique(CP$gs_subcat)

#dim(ttestout)
entriz = CP[ c("gs_name","gene_symbol")] %>% as.data.frame() # CP$human_gene_symbol %in% ttestout$symbol,
entriz <- split(entriz$gene_symbol, entriz$gs_name)
# gsva1
dat = data.frame(t(ptv3_Hsamp_ratio[c(sensitive$samp_ID, nonRes$samp_ID), 3:ncol(ptv3_Hsamp_ratio)]))
dat[] = lapply(dat, as.numeric)
dim(dat) # 9143  192
dat$symbol = as.character(lapply(rownames(dat), function(x){
  x = strsplit(x, '_')[[1]][2]
  x = gsub('.', '-',x,fixed = T)
  x
}))
dim(dat)

dat$symbol[grepl('NA$',dat$symbol,fixed = T)]

dat = subset(dat,  symbol !='NA' & !is.na(symbol))
dat1 = aggregate(dat, by = list(dat$symbol), FUN = function(x){  mean(x, na.rm = T)} )
dim(dat )
dat1$symbol = dat1
rownames(dat1) = dat1$Group.1
dat1 = dat1[2:193]

gsvaPar1 = gsva(expr = as.matrix(dat1) , 
     gset.idx.list = entriz,
     kcdf = "Gaussian" , # ("Gaussian", "Poisson", "none")
     method = "gsva",  # ("gsva", "ssgsea", "zscore", "plage")
     min.sz = 2,
     ssgsea.norm=T,
     verbose=F, 
     parallel.sz = parallel::detectCores())
gsvaPar1 = data.frame(gsvaPar1)
# write.xlsx(gsvaPar1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/gsva_C2-CP.xlsx', rowNames = T)
write.xlsx(gsvaPar1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/ptv3ratio_gsva_hallmark.xlsx', rowNames = T)
dim(gsvaPar1)

gsvaPar1_ttest = ttest(gsvaPar1[c(sensitive$samp_ID, nonRes$samp_ID)], 1, nrow(sensitive), 1+nrow(sensitive), ncol(gsvaPar1))
gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC')]

write.xlsx(gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC', "pre","post","sub")], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/valid5/gsva_C2-CP_ratio_ttest.xlsx', rowNames = T)

dim(subset(gsvaPar1_ttest, Pvalue<0.05) )#& abs( logFC)>log2(1.5) 
dim(subset(gsvaPar1_ttest, FDR<0.05))

gsvaPar1_ttest1 = subset(gsvaPar1_ttest, FDR<0.05 )
dim(gsvaPar1_ttest1)
gsvaPar1_ttest1$pre = rowMeans( gsvaPar1_ttest1[sensitive$samp_ID] , na.rm = TRUE)
gsvaPar1_ttest1$post = rowMeans( gsvaPar1_ttest1[nonRes$samp_ID] , na.rm = TRUE)
gsvaPar1_ttest1$sub = gsvaPar1_ttest1$pre-gsvaPar1_ttest1$post
dim(gsvaPar1_ttest1)
gsvaPar1_ttest1 = gsvaPar1_ttest1[order(gsvaPar1_ttest1[,'sub']),]
# gsvaPar1_ttest1= gsvaPar1_ttest1[c(1:15, 210:224),]
gsvaPar1_ttest1= gsvaPar1_ttest1[c(1:15, 131:145), ]
dim(gsvaPar1_ttest1)
colnames(ptv3_Hsamp_sampinfo)

heatmapInfo = subset(ptv3_Hsamp_sampinfo, PRISM1st_label_total %in% c("non-responsive", "sensitive"))
heatmapInfo = heatmapInfo[c(1,2, 8, 12,13, 18, 23, 24, 31  )]
rownames(heatmapInfo) = heatmapInfo$samp_ID
# heatmapInfo = heatmapInfo[colnames(gsvaPar1_ttest1)[1:192],]
colnames(heatmapInfo)
heatmapInfo = heatmapInfo[order(heatmapInfo[,5]),]
heatmapInfo  = subset(heatmapInfo , !is.na(heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO))
heatmapInfo1 = data.frame(PRISM = heatmapInfo$PRISM1st_label_total, CCK8 = heatmapInfo$CCK8_cell_viability_before,
                          PRISM1st_logFC = heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO,
                          pert_time = paste0(heatmapInfo$pert_time,'h'),
                          row.names = rownames(heatmapInfo))

rownames(gsvaPar1_ttest1) = tolower(rownames(gsvaPar1_ttest1))


pert_time <- c("#007694", "#7DDFFF")
names(pert_time) <- c(  "24h", "6h")
ann_colors <- list(pert_time = pert_time, PRISM1st_logFC =c( "#366092","#C4D6E9") )
ann_colors


#heatmapInfo1 = subset(heatmapInfo1, !is.na(heatmapInfo1$PRISM1st_logFC))
library(pheatmap)
pheatmap(gsvaPar1_ttest1[rownames(heatmapInfo1) ], # scale = 'row',rownames(heatmapInfo)
         # color = colorRampPalette(c("blue", "white","red" ))(1000),
         color = colorRampPalette(c("#00B0F0", "white","#E74809" ))(1000), 
         #breaks = c(seq(-3, 3, length=100)),
         annotation_col = heatmapInfo1, # 样本信息
         annotation_colors = ann_colors , 
         labels_col = heatmapInfo$drugname, 
         angle_col = 45,
         cluster_cols = F, 
         cluster_rows = T, 
         #clustering_method = 'median',
         cellwidth = 1,
         cellheight = 8,#2, # 
         fontsize_col = 1.5,
         fontsize_row = 7,
         show_rownames = T,
         show_colnames = T,
         name = ' ')


colnames(ttestout)

#### 20241223_throid_diann_RF_report.pg_matrix
ptv3_info = read.csv('Z:/members/wangyingrui/TTTD/drug/20241127_ptv3_Hsamp_CAL-62_224samp_9142prot_sampinfo.csv')
rownames(ptv3_info) = ptv3_info$samp_ID

ptv3 = read.csv("Z:/members/wangyingrui/TTTD/drug/20241223_throid_diann_RF_report.pg_matrix.tsv", sep = '\t')
ptv3unigene = ptv3[c(1,4)]
rownames(ptv3unigene) = ptv3unigene$Protein.Group

ptv3[1:3, 1:6]
rownames(ptv3) = ptv3$Protein.Group
ptv3 = ptv3[6:ncol(ptv3)]
ptv3[1:3, 1:6]
cols = c()
for(x in colnames(ptv3)){
  if(grepl('poo', x)) {
    cols = c(cols, x)
    next
  }
  x = gsub('_rep', 'rep', x)
  x = strsplit(gsub('.mzML', '', x, fixed = T),'_')[[1]]
  n = length(x)
  x = x[n]
  if (x %in% cols){
    print(x)
    cols = c(cols, paste0(x,'rep'))
  }else{
    cols = c(cols, x)
  }
}
colnames(ptv3) = cols
dim(ptv3)
ptv3 = ptv3[colSums(!is.na(ptv3))>2000]
dim(ptv3)
ptv3_min08 = ptv3
ptv3_min08[is.na(ptv3_min08)] = min(ptv3, na.rm = T)*0.8 # 851.88

reps = colnames(ptv3)[grepl('rep', colnames(ptv3))]
reps = c(reps, gsub('rep','', reps))
reps = intersect(sort( reps), colnames(ptv3))
reps = ptv3[reps]

pool = colnames(ptv3)[grepl('poo', colnames(ptv3))]
pool = ptv3[pool]
dim(pool)
narow = rowSums(!is.na(pool))
pool = pool[narow>0, ]
dim(pool)

poolcor = cor(pool, method = "spearman", use = 'pairwise.complete.obs')
poolcor
diag(poolcor) = NA
poolCor = c()
for(i in 1:ncol(poolcor)){
  for (j in i:nrow(poolcor)) {
    poolCor = c(poolCor, poolcor[i, j])
    if(j == nrow(poolcor))break
  }
  if(i == nrow(poolcor))break
}
poolCor = poolCor[!is.na(poolCor)]
poolCor = data.frame(poolCor)

ggplot(poolCor, aes( x = 'Cor', poolCor))+
  geom_violin( fill = "#00B0F0")+
  geom_boxplot(width = 0.1 )+
  labs(x = '', y='Cor')+
  annotate('text', x = "Cor", y = 0.95, label = paste0('Median=',round(median(poolCor$poolCor, na.rm = T), 3)))+
  theme_classic()+
  ylim(0.8, 1)+
  theme(text = element_text(size = 16))

set.seed(2024)
dfpca = prcomp(t(ptv3_min08[colnames(pool)]) )
dfpca = data.frame(dfpca$x)
dfpca[1:3,1:3]
dfpca$label = rownames(dfpca)
ggplot(dfpca, aes( PC1, PC2))+
  geom_point( )+
  geom_text_repel(aes(label=label), color = 'black')+
  theme_classic()+
  theme(text = element_text(size = 16))


repCor = c()
for(i in seq(1, ncol(reps), 2) ){
  temp = reps[c(i, i+1)]
  repcor = cor(temp, method = "spearman", use = 'pairwise.complete.obs')
  repcor
  repCor = rbind(repCor, c(repcor[1,2], colnames(repcor)[1]))
}
repCor = data.frame(repCor)
repCor$X1 = as.numeric(repCor$X1)

ggplot(repCor, aes( x = 'Cor', X1))+
  geom_violin( fill = "#00B0F0")+
  geom_boxplot(width = 0.1 )+
  labs(x = '', y='Cor')+
  annotate('text', x = "Cor", y = 0.98, label = paste0('Median=',round(median(repCor$X1, na.rm = T), 3)))+
  theme_classic()+
  ylim(0.9, 1)+
  theme(text = element_text(size = 16))


setdiff(gsub('rep', '', colnames(ptv3)), ptv3_info$samp_ID)

ptv3_Hsamp_ratio = data.frame(t(ptv3))
for(nr in 1:nrow(ptv3_Hsamp_ratio)){
  sid = rownames(ptv3_Hsamp_ratio)[nr]
  if(grepl('poo', sid)) next
  temp = ptv3_info[ptv3_info$samp_ID == gsub('rep', '', sid), ]
  if(nrow(temp)==0)next
  if (temp$pert_id=='control')next
  controlid = temp$control
  ptv3_Hsamp_ratio[nr, ] = ptv3_Hsamp_ratio[nr, ]/ptv3_Hsamp_ratio[controlid,  ]
}
ptv3_Hsamp_ratio[1:3,1:6]
ptv3_Hsamp_ratio = ptv3_Hsamp_ratio[!(rownames(ptv3_Hsamp_ratio) %in% c(colnames(pool), "H4992", "H5103", "H4993", "H4993rep")), ]

sensitive = subset(ptv3_info,  PRISM1st_label_total=="sensitive" &  samp_ID %in% rownames(ptv3_Hsamp_ratio) )
nonRes = subset(ptv3_info,  PRISM1st_label_total=="non-responsive"& samp_ID %in% rownames(ptv3_Hsamp_ratio))

ttestout = ttest(t(ptv3_Hsamp_ratio[c(sensitive$samp_ID, nonRes$samp_ID), ]), 1, nrow(sensitive), 1+nrow(sensitive), nrow(nonRes)+nrow(sensitive) )
ttestout[1:2, 1:2]
ttestout = data.frame(ttestout)
ttestout$symbol = ptv3unigene[rownames(ttestout),  'Genes']
write.xlsx(ttestout , "D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptv320241223_throid_sensitive-nonres_ttest.xlsx", rowNames = T)
nrow(subset(ttestout, FDR<0.05))
nrow(subset(ttestout, Pvalue<0.05))

ttestout$Sig = ifelse(is.na(ttestout$FDR), "None", ifelse(ttestout$logFC> log2(1.2) & ttestout$FDR<0.05, "Up",
                                                          ifelse( ttestout$logFC< -log2(1.2) & ttestout$FDR<0.05, "Down", "None")))
table(ttestout$Sig)
ttestout1 = subset(ttestout, FDR<0.05 & abs( logFC)>log2(1.2) )
ggplot(subset(ttestout, logFC <1  & logFC > -1.5 ),
       aes(x = logFC, y = -log10(FDR), color = Sig )) +#  , size = -log10(qt_padjust)
  geom_point(alpha=0.7, size=1 )+
  scale_color_manual(values=c( "#546de5", "#d2dae2", "#ff4757"))+# , "#d2dae2"
  geom_text_repel( data = ttestout1,
                   aes(x=logFC, y = -log10(FDR) ,label=symbol),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
                   segment.color = "black", show.legend = FALSE )+
  labs(x="log2(Fold change)", y="-log10(p-adjust)", title = '' )+
  theme_classic()+ scale_y_continuous(n.breaks = 10) + #ylim(0, 60)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = log2(1.2), linetype = "dashed")+
  geom_vline(xintercept = -log2(1.2), linetype = "dashed")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


ptv3_Hsamp_ratio_min08 = ptv3_Hsamp_ratio
ptv3_Hsamp_ratio_min08[is.na(ptv3_Hsamp_ratio_min08)] = min(ptv3_Hsamp_ratio_min08, na.rm = T)*0.8
ptv3_Hsamp_ratio_min08[1:3,1:6]


set.seed(2024)
dfpca = prcomp(ptv3_Hsamp_ratio_min08)
dfpca = data.frame(dfpca$x)
dfpca$label = ptv3_info[gsub('rep', '', rownames(dfpca)) , 'instrument']
dfpca = dfpca[c('PC1', 'PC2', 'label')]
unique(dfpca$label)
dfpca[is.na(dfpca$label),]

library(umap)
library(Rtsne)
set.seed(2024)
dftsne <- Rtsne(ptv3_Hsamp_ratio_min08, dims = 2, perplexity = 10,
                   partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
df1 = data.frame(dftsne$Y) 
rownames(df1) = rownames(ptv3_Hsamp_ratio_min08)
df1$label = ptv3_info[gsub('rep', '', rownames(df1)) , 'instrument']
head(df1)

# ggplot(dfpca, aes( PC1, PC2, color = label ))+
ggplot(df1, aes( X1, X2, color = label ))+
  geom_point( )+
  #geom_text_repel(aes(label=label) )+
  theme_classic()+
  scale_color_d3()+
  theme(text = element_text(size = 16))


library(sva)
ptv3_Hsamp_ratio_min08[1:3,1:6]
dim(ptv3_Hsamp_ratio_min08)

ptv3_info$time = gsub('O', '', ptv3_info$batch)
ptv3_info$time = gsub('CAA', '', ptv3_info$time)
set.seed(2024)
ptv3_Hsamp_ratio_min08combat = ComBat(t(ptv3_Hsamp_ratio_min08), batch = ptv3_info[gsub('rep', '', rownames(ptv3_Hsamp_ratio_min08)) , 'time'] )
ptv3_Hsamp_ratio_min08combat[1:3,1:6]

ptv3_min08[1:3,1:3]
ptv3_min08combat = ptv3_min08[!grepl('mzML', colnames(ptv3_min08))]
ptv3_min08combat = ptv3_min08combat[!(colnames(ptv3_min08combat) %in% c( "H4992", "H5103", "H4993", "H4993rep"))  ]
colnames(ptv3_min08combat)
ptv3_info[gsub('rep', '', colnames(ptv3_min08combat)) , c('pert_id', "time")]
ptv3_min08combat  = ComBat(ptv3_min08combat, batch = ptv3_info[gsub('rep', '', colnames(ptv3_min08combat)) , 'time'] )
ptv3_min08combat[1:3, 1:3]
ptv3_min08combat = data.frame(ptv3_min08combat)
write.xlsx(ptv3_min08combat, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptvmin0.8combat.xlsx', rowNames = T)
ptv3_min08combat_log2 = log2(ptv3_min08combat)
# ptv3_Hsamp_ratio_min08combat = ComBat(ptv3_Hsamp_ratio_min08combat, batch = ptv3_info[gsub('rep', '', rownames(df1)) , 'instrument'] )
# ptv3_Hsamp_ratio_min08combat[1:3,1:6]

set.seed(2024)
dftsne <- Rtsne( t(ptv3_Hsamp_ratio_min08combat), #  t(ptv3_Hsamp_ratio_min08combat)  ptv3_Hsamp_ratio_min08
                dims = 2, perplexity = 10,
                partial_pca=TRUE, verbose = T , check_duplicates = FALSE)
df1 = data.frame(dftsne$Y) 
rownames(df1) = rownames(ptv3_Hsamp_ratio_min08)
df1$label = ptv3_info[gsub('rep', '', rownames(df1)) , 'batch']
df1$instrument = ptv3_info[gsub('rep', '', rownames(df1)) , 'instrument']
df1$PRISM1st = ptv3_info[gsub('rep', '', rownames(df1)) , 'PRISM1st_label']
head(df1)

ggplot(df1, aes( X1, X2, color = instrument ))+
  geom_point( )+
  #geom_text_repel(aes(label=instrument) )+
  #labs(title = 'combat')+
  theme_classic()+
  scale_color_d3(palette ="category20")+
  theme(text = element_text(size = 16))


# ggplot(dfpca, aes( PC1, PC2, color = label ))+
ggplot(df1, aes( X1, X2, color = label ))+
  geom_point( )+
  geom_text_repel(aes(label=instrument) )+
  labs(title = 'combat')+
  theme_classic()+
  scale_color_d3(palette ="category20")+
  theme(text = element_text(size = 16))


library(GSVA)
library(msigdbr)
CP <- msigdbr(species = "Homo sapiens", category = "C2" ) #%>% dplyr::select(gs_subcat, gs_name, entrez_gene, gene_symbol, human_gene_symbol )
CP = subset(CP, gs_subcat!="CGP")
# CP <- msigdbr(species = "Homo sapiens", category = "H" )
head(CP)
unique(CP$gs_subcat)

#dim(ttestout)
entriz = CP[ c("gs_name","gene_symbol")] %>% as.data.frame() # CP$human_gene_symbol %in% ttestout$symbol,
entriz <- split(entriz$gene_symbol, entriz$gs_name)
dat = data.frame(ptv3_Hsamp_ratio_min08combat)
dat = dat[ c(sensitive$samp_ID, nonRes$samp_ID)] 
#dat[] = lapply(dat, as.numeric)
dim(dat) # 9143  192
dat$symbol = ptv3unigene[rownames(dat), 'Genes']
dim(dat)

dat$symbol[grepl('NA$',dat$symbol,fixed = T)]

dat = subset(dat,  symbol !='NA' & !is.na(symbol) & symbol !='' )
dat1 = aggregate(dat, by = list(dat$symbol), FUN = function(x){  mean(x, na.rm = T)} )
dim(dat )
rownames(dat1) = dat1$Group.1
dat1[1:3, 185:194]
dat1 = dat1[2: 193]
# gsva2
gsvaPar1 = gsva(expr = as.matrix(dat1) , 
                gset.idx.list = entriz,
                kcdf = "Gaussian" , # ("Gaussian", "Poisson", "none")
                method = "gsva",  # ("gsva", "ssgsea", "zscore", "plage")
                min.sz = 2,
                ssgsea.norm=T,
                verbose=F, 
                parallel.sz = parallel::detectCores())
gsvaPar1 = data.frame(gsvaPar1)
write.xlsx(gsvaPar1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptv3_ratio_combatTime_gsva.xlsx', rowNames = T)

dim(gsvaPar1)
gsvaPar1_ttest = ttest(gsvaPar1[c(sensitive$samp_ID, nonRes$samp_ID)], 1, nrow(sensitive), 1+nrow(sensitive), ncol(gsvaPar1))
head(gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC')])
dim(subset(gsvaPar1_ttest, FDR<0.05))

write.xlsx(gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC', 'differ')], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptv3_ratio_combatTime_gsva_ttest.xlsx', rowNames = T)

gsvaPar1_ttest1 = subset(gsvaPar1_ttest, Pvalue<0.05)
dim(gsvaPar1_ttest1)
gsvaPar1_ttest1 = gsvaPar1_ttest1[order(gsvaPar1_ttest1[, 'differ']),]
rownames(gsvaPar1_ttest1) = tolower(rownames(gsvaPar1_ttest1))

heatmapInfo1 = data.frame(PRISM = ptv3_info[c(sensitive$samp_ID, nonRes$samp_ID), 'PRISM1st_label_total'],
                          CCK8 =  ptv3_info[c(sensitive$samp_ID, nonRes$samp_ID), 'CCK8_cell_viability_before'],
                          PRISM1st_logFC = ptv3_info[c(sensitive$samp_ID, nonRes$samp_ID), 'PRISM1st_logfold_change_values_relative_to_DMSO'],
                          pert_time = ptv3_info[c(sensitive$samp_ID, nonRes$samp_ID), 'pert_time'],
                          row.names = c(sensitive$samp_ID, nonRes$samp_ID) )
heatmapInfo1$pert_time = paste0(heatmapInfo1$pert_time, 'h')

write.xlsx(heatmapInfo1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptv3_ratio_combatTime_gsva_ttest_pval_differ20_heatmapinfo.xlsx', rowNames = T)

pert_time <- c("#007694", "#7DDFFF")
names(pert_time) <- c(  "24h", "6h")
ann_colors <- list(pert_time = pert_time, PRISM1st_logFC =c( "#366092", "#C4D6E9") )
ann_colors

dim(gsvaPar1_ttest1)
pheatmap(gsvaPar1_ttest1[c(1:20, 234:253) ,1:187],
         color = colorRampPalette(c("#00B0F0", "white", "#E74809"))(1000), 
         annotation_col = heatmapInfo1,
         annotation_colors = ann_colors,
         fontsize = 10,
         angle_col = 45,
         fontsize_row = 7,
         fontsize_col = 1,
         cellwidth = 1.5,
         cellheight = 8,
         labels_col = ptv3_info[c(sensitive$samp_ID, nonRes$samp_ID), 'PRISM1_match'],
         cluster_rows = F,
         cluster_cols = F)

###### 20250103 All_DEPs_20241204, ptv3_Hsamp #####
ptv3_Hsamp = read.csv('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/20241127_ptv3_Hsamp_CAL-62_224samp_9142prot_sampinfo.csv')
alldep = read.xlsx('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204.xlsx')
unis = c('P53985')
for(i in 1:nrow(ptv3_Hsamp)){
  #d = ptv3_Hsamp$drugname
  #uni = 
  unis = c(unique(unis), strsplit(ptv3_Hsamp[i, 'targetv2'],';')[[1]])
}
unis = unique(unis[!is.na(unis)])
uni_drug = data.frame(uniprot = unis,row.names = unis)
colnames(alldep)
head(alldep)
alldep[nrow(alldep)+1,] = NA

for(cols in colnames(alldep)[1:5]){
  # uni = "Q9P2P6"
  # cols="NBM_l_overlap"
  unis = alldep[,cols]
  ds = c()
  for (uni in unis){
    if(is.na(uni)){
      alldep[paste0(cols,'_drugs')] = c(ds, rep(NA, nrow(alldep)-length(ds)))
      break
    }
    temp = subset(ptv3_Hsamp, grepl(uni, ptv3_Hsamp$targetv2) & !is.na(ptv3_Hsamp$drugname)) #  & PRISM1st_label_total =="sensitive"
    if(nrow(temp)==0){
      d = NA
    }else{
      d = paste0(unique(temp$drugname),collapse = ';')
    }
    ds = c(ds, gsub(';NA','',d) )
  }
  
}
library(clusterProfiler)
library(org.Hs.eg.db)

for(cols in colnames(alldep)[6:9]){
  # uni = "Q9P2P6"
  # cols = "Hmonkey"
  alldep[cols] = gsub(' ', '', alldep[,cols])
  genes = gsub(' ', '', alldep[,cols])
  ds = c()
  g_u = bitr(gsub(' ','',genes), fromType = 'SYMBOL' , toType = c('UNIPROT'), OrgDb = 'org.Hs.eg.db')
  for (g in genes){
    if(is.na(g)){
      alldep[paste0(cols,'_drugs')] = c(ds, rep(NA, nrow(alldep)-length(ds)))
      break
    }
    temp_gu = subset(g_u, SYMBOL == g)
    ds1 = c()
    for(uni in temp_gu$UNIPROT){
      temp = subset(ptv3_Hsamp, grepl(uni, ptv3_Hsamp$targetv2) & !is.na(ptv3_Hsamp$drugname))#& PRISM1st_label_total =="sensitive"
      if(nrow(temp)==0){
        d = NA
      }else{
        print(c(cols, g,uni))
        # break
        d = paste0(unique(temp$drugname),collapse = ';')
      }
      ds1 = c(ds1, d)
    }
    ds1 = paste0(ds1,collapse =';')
    ds = c(ds, gsub(';NA','',ds1))
  }
}
write.xlsx(alldep[sort(colnames(alldep))], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204_mapPtv3.xlsx')
#write.xlsx(alldep[sort(colnames(alldep))], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204_mapPtv3sensitive1.xlsx')

alldep = read.xlsx('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204_mapPtv3sensitive1.xlsx')
alldep = alldep[seq(2,ncol(alldep),2)]
alldep$id = rownames(alldep)
sensitiveDrug = melt(alldep, id.vars = 'id')
sensitiveDrug = unique(sensitiveDrug$value)[3:14]
temp=c()
for(i in sensitiveDrug){
  temp=c(temp, strsplit(i,';')[[1]])
}
sensitiveDrug = unique(temp)
sensitiveDrug

dat = data.frame(ptv3_min08)
dat = dat[c(sensitive[sensitive$drugname %in% sensitiveDrug, 'samp_ID'], unique(ptv3_info$control))] 

for(i in sensitiveDrug){
  #print(i)
  sid = sensitive[sensitive$drugname == i, 'samp_ID']
  controlid= ptv3_info[sid,  'control']
  print(c(sid, controlid))
  dat_dsmo = ttest_paired(dat[c(sid, controlid)], 1,2,3,4)
  dat_dsmo$symbol = ptv3unigene[rownames(dat_dsmo),'Genes']

  targetv2 = ptv3_info[sid,  'targetv2']
  print(c(i, unique(target)))
  targetv2 = strsplit(unique(targetv2), ';')[[1]]
  targetv2 = intersect(targetv2, rownames(dat_dsmo))
  if(length(targetv2)>0){
    dat_dsmo[targetv2, 'ptv_targetv2'] = 'ptv_targetv2'
  }
  #dat_dsmo = subset(dat_dsmo, dat_dsmo$Pvalue<0.05)
  pnum = nrow(subset(dat_dsmo, Pvalue <0.05))
  qnum = nrow(subset(dat_dsmo, FDR <0.05))
  write.xlsx(dat_dsmo[c(5:8, 10:ncol(dat_dsmo)) ], paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptvmin0.8_depSensitiveDrug_", i,'_pairedTtest_pnum',pnum,'_qnum',qnum,'.xlsx'), rowNames = T)
  
}

alldep = read.xlsx('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204_mapPtv3sensitive1.xlsx')
files = list.files('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/')
files = files[grepl('ptvmin0.8_', files)]
files = files[grepl('pnum', files)]
files
for(f in files){
  f = "ptvmin0.8_depSensitiveDrug_Ispinesib_pairedTtest_pnum295_qnum0.xlsx"
  d = strsplit(f,'_')[[1]][3]
  dat_dsmo = read.xlsx(paste0('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/',f), rowNames = T)
  if(ncol(dat_dsmo)<6)next
  # print(f)
  unis = subset(dat_dsmo , !is.na(ptv_targetv2) )
  unis = rownames(unis)
  temp = alldep[apply(alldep, 1, function(x) any(grepl(d, x))), ]
  ind = which(apply(temp, 2, function(col) any(grepl(d, col))))
  
  d_TTTDdeps = c()
  for (i in ind){
    if(i<9){
      temp1 = temp[c(i-1,i)]
      temp1 = temp1[!is.na(temp1[2]), 1]
      entriz = bitr(gsub(' ','',temp1), fromType = 'SYMBOL' , toType = c('UNIPROT'), OrgDb = 'org.Hs.eg.db')
      if (nrow(entriz)==0)next
      temp1 = unique(entriz$UNIPROT)
      d_TTTDdeps = c(d_TTTDdeps, temp1)
      next
    }
    temp1 = temp[c(i-1,i)]
    temp1 = temp1[!is.na(temp1[2]), 1]
    d_TTTDdeps = c(d_TTTDdeps, temp1)
  }
  d_TTTDdeps = unique(d_TTTDdeps)
  #print(c(f, d_TTTDdeps))
  overlapPro = intersect(unis, d_TTTDdeps)
  if(length(overlapPro)==0)next
  print(c(d, overlapPro))
  
  dat_dsmo$Sig = ifelse(is.na(dat_dsmo$Pvalue), "None", ifelse(dat_dsmo$logFC> log2(1.2) & dat_dsmo$Pvalue<0.05, "Up",
                                                            ifelse( dat_dsmo$logFC< -log2(1.2) & dat_dsmo$Pvalue<0.05, "Down", "None")))
  #table(dat_dsmo$Sig)
  print(dat_dsmo[overlapPro, ])
  p = ggplot(dat_dsmo,
         aes(x = logFC, y = -log10(Pvalue), color = Sig )) +#  , size = -log10(qt_padjust)
    geom_point(alpha=0.7, size=1 )+
    scale_color_manual(values=c( "#546de5", "#d2dae2", "#ff4757"))+# , "#d2dae2"
    geom_text_repel( data = dat_dsmo[overlapPro, ],
                     aes(x=logFC, y = -log10(Pvalue) ,label=symbol),  size = 3,  # color = "black", # 不固定颜色使其根据Sig修改颜色
                     segment.color = "black", show.legend = FALSE )+
    labs(x="log2(Fold change)", y="-log10(pvalue)", title = '' )+
    theme_classic()+ 
    scale_y_continuous(n.breaks = 10) + #ylim(0, 60)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    geom_vline(xintercept = log2(1.2), linetype = "dashed")+
    geom_vline(xintercept = -log2(1.2), linetype = "dashed")+
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position="right",
          axis.text = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 15, color = "black"),
          axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_text(size = 15, color = "black"),
          #legend.title = element_blank(),
          legend.text = element_text(size = 15))
  ggsave(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/", gsub('xlsx', 'pdf', f)), p)
}
all_dsmo = data.frame(symbol= ptv3unigene[rownames(dat),'Genes'], row.names = rownames(dat))
for(i in sensitiveDrug){
  #i="Luminespib" 
  sid = sensitive[sensitive$drugname == i & sensitive$pert_time == 24, 'samp_ID'] # 24h
  controlid= ptv3_info[sid, 'control']
  print(c(sid, controlid))
  dat_dsmo = dat[c(sid, controlid)]
  # print(ncol(dat_dsmo))
  # dat_dsmo[paste0(i,'_FC')] = dat[,1]/dat[,2]
  # dat_dsmo$equal = dat_dsmo[,1]==dat_dsmo[,2]
  pre = rowMeans(dat_dsmo[ sid],na.rm = T)
  post = rowMeans(dat_dsmo[controlid],na.rm = T)
  dat_dsmo[paste0(i,'_FC')] = dat_dsmo[,1]/dat_dsmo[,2] #pre/post
  #dat_dsmo$symbol = ptv3unigene[rownames(dat_dsmo),'Genes']

  targetv2 = ptv3_info[sid,  'targetv2']
  print(c(i, unique(target)))
  targetv2 = strsplit(unique(targetv2), ';')[[1]]
  targetv2 = intersect(targetv2, rownames(dat_dsmo))
  if(length(targetv2)>0){
    dat_dsmo[targetv2, paste0(i,'_ptvTargetv2')] = 'ptv_targetv2'
  }
  colnames(dat_dsmo)[3:ncol(dat_dsmo)]
  all_dsmo = cbind(all_dsmo, dat_dsmo[3:ncol(dat_dsmo)])
  #
}
write.xlsx(all_dsmo , paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptvmin0.8_depSensitiveDrug_24hsampleRatio.xlsx"), rowNames = T)

all_dsmo1 = all_dsmo
rownames(all_dsmo1) = paste0(rownames(all_dsmo1),'_',all_dsmo1$symbol)
all_dsmo1 = all_dsmo1[seq(2, 27,2)]
for(i in 1:nrow(all_dsmo1)){
  # n = unique(as.numeric(all_dsmo1[i, ]))
  # if (length(n)==1){all_dsmo1[i, ] = NA}
}
all_dsmo1$meanFC = rowMeans(all_dsmo1, na.rm = T)

# all_dsmo1$label = ifelse( all_dsmo1$meanFC>2, 'up', ifelse( all_dsmo1$meanFC<0.5, 'down', 'None'))
all_dsmo1$up = apply(all_dsmo1[1:13], 1, function(x) sum(x > 2))
all_dsmo1$dn = apply(all_dsmo1[1:13], 1, function(x) sum(x < 0.5))
all_dsmo1$label = 'None'
all_dsmo1$label = ifelse( all_dsmo1$up>=7, 'up', ifelse( all_dsmo1$dn>=7, 'down', 'None'))
# all_dsmo[c("meanFC", "up", "dn", "label")] = all_dsmo1[c("meanFC", "up", "dn", "label")]
# write.xlsx(all_dsmo , paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptvmin0.8_depSensitiveDrug_24hsampleRatio.xlsx"), rowNames = T)

dim(all_dsmo1)
ggplot(all_dsmo , aes(dn/13, up/13, color = label))+
  geom_point()+ 
  theme_classic()+
  geom_text_repel(data = all_dsmo[all_dsmo$label !='None', ], aes(label = symbol), segment.color = "black", show.legend = FALSE)+
  scale_color_manual(values=c( "#546de5", "black", "#ff4757")) # "#d2dae2"

# 20250120
sids = sensitive[sensitive$drugname %in% sensitiveDrug , ]
sids$samp_ID
unique(sids$control)
ids = c(sids$samp_ID, unique(sids$control))
temp = ptv3[ids]
temp[1:3, 1:3]
write.xlsx(temp, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptv_13sensitivedrugs_4DMSO_mat.xlsx', rowNames = T)


sid = sensitive[sensitive$drugname %in% sensitiveDrug & sensitive$pert_time == 24, 'samp_ID'] # 24h
controlid= unique(sensitive[sid, 'control'])
dat_combat = data.frame(ptv3_min08combat) 
dat_combat[dat_combat<0] = NA
dat_combat[is.na(dat_combat)] = min(dat_combat, na.rm = T)*0.8
dat_combat = log2(dat_combat[c(sid, controlid)])
ttestout1 = ttest1_log2(dat_combat, 1, length(sid) , 1+length(sid), length(sid) +length(controlid))

# ttestout1 = ttest(dat_combat, 1, length(sid) , 1+length(sid), length(sid) +length(controlid))
ttestout1$Sig = ifelse(is.na(ttestout1$Pvalue), "None", ifelse(ttestout1$logFC> log2(2) & ttestout1$FDR<0.05, "Up",
                                                             ifelse( ttestout1$logFC< -log2(2) & ttestout1$FDR<0.05, "Down", "None")))
colnames(ttestout1)
ttestout1$symbol = ptv3unigene[rownames(ttestout1), 'Genes']
table(ttestout1$Sig)
write.xlsx(ttestout1[c('symbol',"Pvalue", "FDR", "FC", "logFC",'differ', "Sig")], paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptvmin0.8combat_NAchange_log2_depSensitiveDrug_24hsample_FC2_ttest.xlsx"), rowNames = T)

temp = dat_combat#[c('O15066', 'O60333', 'P0CAP2', 'Q58FG1', 'Q8NFQ8'),]
temp$pro = rownames(temp)
temp = melt(temp, id.vars = 'pro')
temp = data.frame(temp)
temp$type =c(rep("sensitive" ,65), rep('DMSO', 10))
head(temp)
temp$type = factor(temp$type, levels = c("sensitive" , 'DMSO'))
ggplot(temp, aes(type, value, color = type))+
  geom_boxplot(lwd = 1, width = 0.8)+
  geom_jitter(alpha = 0.7)+
  theme_classic()+
  #geom_smooth(aes(  fill = variable), se = T, method = 'loess', alpha = 0.1)+ #  alpha = 0.15
  #labs(y = 'Log2 Expressions' )+
  scale_color_d3( )+
  #scale_fill_d3(palette ="category20")+
  theme(axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        text = element_text(size = 15))+
  scale_y_continuous(n.breaks = 8)+
  stat_summary(fun.data = function(x) data.frame(y=13 , label = paste("Median=", round(median(x), 3))), geom="text", size = 4)

# ttestout1 = read.xlsx(paste0("D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/ptvmin0.8combat_depSensitiveDrug_24hsampleRatio_FC2_ttest.xlsx"), rowNames = T)

library(GSVA)
library(msigdbr)
CP <- msigdbr(species = "Homo sapiens", category = "C2" ) #%>% dplyr::select(gs_subcat, gs_name, entrez_gene, gene_symbol, human_gene_symbol )
# CP = subset(CP, gs_subcat!="CGP")
CP <- msigdbr(species = "Homo sapiens", category = "H" )
head(CP)
unique(CP$gs_subcat)

#dim(ttestout)
entriz = CP[ c("gs_name","gene_symbol")] %>% as.data.frame() # CP$human_gene_symbol %in% ttestout$symbol,
entriz <- split(entriz$gene_symbol, entriz$gs_name)

sigpro = rownames(subset(ttestout1, Sig != 'None'))
dat = dat_combat[ c(sid, controlid)]# sigpro,
datgsva = dat[c(sid, controlid)]
#dat[] = lapply(dat, as.numeric)
dim(dat) # 9143  192
datgsva$symbol = ptv3unigene[rownames(datgsva), 'Genes']
dim(datgsva)

datgsva$symbol[grepl('NA$',datgsva$symbol,fixed = T)]

datgsva = subset(datgsva,  symbol !='NA' & !is.na(symbol) & symbol !='' )
datgsva1 = aggregate(datgsva, by = list(datgsva$symbol), FUN = function(x){  mean(x, na.rm = T)} )
dim(datgsva )
rownames(datgsva1) = datgsva1$Group.1

datgsva1 = datgsva1[2:( ncol(datgsva1)-1)]
# gsva2
gsvaPar1 = gsva(expr = as.matrix(datgsva1) , 
                gset.idx.list = entriz,
                kcdf = "Gaussian" , # ("Gaussian", "Poisson", "none")
                method = "gsva",  # ("gsva", "ssgsea", "zscore", "plage")
                min.sz = 2,
                ssgsea.norm=T,
                verbose=F, 
                parallel.sz = parallel::detectCores())
gsvaPar1 = data.frame(gsvaPar1)
#write.xlsx(gsvaPar1, 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptvmin0.8combat_NAchange_depSensitiveDrug_24hsampleRatio_gsva.xlsx', rowNames = T)

dim(gsvaPar1)
gsvaPar1_ttest = ttest(gsvaPar1[c(sid, controlid)], 1, length(sid), 1+length(sid), ncol(gsvaPar1))
colnames(gsvaPar1_ttest)
head(gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC',"differ")])
dim(subset(gsvaPar1_ttest, Pvalue<0.05))

write.xlsx(gsvaPar1_ttest[c('Pvalue', 'FDR', 'FC', 'logFC', 'differ')], 'D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/ptvmin0.8combat_NAchange_log2_depSensitiveDrug_24hsample_allpro_gsva_ttest.xlsx', rowNames = T)

gsvaPar1_ttest1 = gsvaPar1_ttest[gsvaPar1_ttest$FDR<0.05,]#
dim(gsvaPar1_ttest1)
gsvaPar1_ttest1 = gsvaPar1_ttest1[order(gsvaPar1_ttest1[,'differ']),]
heatmapInfo = ptv3_info[colnames(gsvaPar1), ]
heatmapInfo = heatmapInfo[c(1,2, 8, 12,13, 18, 23, 24, 31  )]
rownames(heatmapInfo) = heatmapInfo$samp_ID
# heatmapInfo = heatmapInfo[colnames(gsvaPar1_ttest1)[1:192],]
colnames(heatmapInfo)
heatmapInfo = heatmapInfo[order(heatmapInfo[,5]),]
# heatmapInfo  = subset(heatmapInfo , !is.na(heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO))
heatmapInfo1 = data.frame(
   PRISM = heatmapInfo$PRISM1st_label_total, CCK8 = heatmapInfo$CCK8_cell_viability_before,
                          PRISM1st_logFC = heatmapInfo$PRISM1st_logfold_change_values_relative_to_DMSO,
                          #pert_time = paste0(heatmapInfo$pert_time,'h'),
                          row.names = rownames(heatmapInfo))

pert_time <- c("#007694", "#7DDFFF")
names(pert_time) <- c(  "24h", "6h")
ann_colors <- list(pert_time = pert_time, PRISM1st_logFC =c( "#366092","#C4D6E9") )
ann_colors

#heatmapInfo1 = subset(heatmapInfo1, !is.na(heatmapInfo1$PRISM1st_logFC))
library(pheatmap)
rownames(gsvaPar1_ttest1) = tolower(rownames(gsvaPar1_ttest1))
pheatmap(gsvaPar1_ttest1[c(sid, controlid)] , # scale = 'row',rownames(heatmapInfo)
         # color = colorRampPalette(c("blue", "white","red" ))(1000),
         color = colorRampPalette(c("#00B0F0", "white","#E74809" ))(1000), 
         #breaks = c(seq(-3, 3, length=100)),
         annotation_col = heatmapInfo1, # 样本信息
         # annotation_colors = ann_colors , 
         labels_col = heatmapInfo$drugname, 
         # angle_col = 45,
         cluster_cols = F,  # cluster_rows = T, 
         #clustering_method = 'median',
         cellwidth = 10, cellheight = 7,#2, # 
         fontsize_col = 8.5, fontsize_row = 6,
         show_rownames = T, show_colnames = T,
         name = ' ')
####### 13sensitiveDrugs_targetv2_sankey #####
alldep = read.xlsx('D:/chh/2023workProject/20240821TTTD/survival_cox/20241122valid/20250103/All_DEPs_20241204_mapPtv3sensitive1.xlsx')
library(ggalluvial)
library(dplyr)
gene = c()
sends = c()
cols = c()
# for(i in sensitiveDrug){
#   sid = sensitive[sensitive$drugname == i, 'samp_ID']
#   targetv2 = ptv3_info[sid,  'targetv2']
#   print(c(i, unique(target)))
#   targetv2 = strsplit(unique(targetv2), ';')[[1]]
#   gene = c(gene, targetv2)
#   sends = c(sends, rep(i, length(targetv2)))
# }

for(i in seq(1,ncol(alldep), 2)){
  col_cur = colnames(alldep)[i+1]
  temp = alldep[!is.na(alldep[col_cur]) & alldep[col_cur]!='',]
  if(nrow(temp)==0)next
  #temp = temp[i:(i+1)]
  if (!grepl(';',temp[ i+1])){
    gene = c(gene, temp[,i])
    sends = c(sends, temp[,i+1] )
    cols = c(cols, rep(col_cur, nrow(temp)))
  }else{
    for(j in 1:nrow(temp)){
      ds = strsplit(temp[j,i+1],';')[[1]]
      for(d in ds){
        gene = c(gene, temp[j,i])
        sends = c(sends, d )
        cols = c(cols, col_cur)
      }
    }
  }
  
}
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
rownames(top_features_) = top_features_$uni
for(i in top_features_$uni){
  #i ='P29372'
  temp = subset(ptv3_Hsamp_sampinfo , grepl(i,ptv3_Hsamp_sampinfo$targetv2) & PRISM1st_label_total=="sensitive")
  if(nrow(temp)==0){next}
  print(i)
  # gene = c(gene, top_features_[i, 'gene'])
  # sends = c(sends, temp$drugname)
  # cols = c(cols, col_cur)
}

sensitive_target = data.frame(gene = gene, Drug = sends, Belong = cols)
head(sensitive_target)
#     gene       Drug
# 1 O14558 Luminespib
# 2 P04792 Luminespib
ptv3unigene[sensitive_target$gene, ]
entriz = bitr(unique(sensitive_target$gene), fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
entriz = subset(entriz, !(UNIPROT %in% ptv3unigene$Protein.Group))
entriz = entriz %>% distinct( UNIPROT, .keep_all = TRUE)
rownames(entriz) = entriz$UNIPROT
colnames(entriz) = colnames(ptv3unigene)
ptv3unigene = rbind(ptv3unigene, entriz)
ptv3unigene['Q14568', 'Genes']='HSP90AA2P'
ptv3unigene['Q58FG0', 'Genes']=  'HSP90AA5P'
ptv3unigene['B7ZC32', 'Genes']= 'KIF28P'
# sensitive_target$symbol = ptv3unigene[sensitive_target$gene, 'Genes']
sensitive_target[ 'symbol'] = ptv3unigene[sensitive_target$gene, 'Genes']
sensitive_target$symbol[1:3] = sensitive_target$gene[1:3]


# summarize the data and count the frequencies
frequencies <- sensitive_target %>%
  count(symbol, Drug) #%>% arrange(gene, desc(n))
ggplot(data = sensitive_target , aes(axis1 = Belong, axis2 = symbol, axis3 = Drug# y = n  Third variable on the X-axis
)) +
  geom_alluvium(aes(fill = Drug)) +#  用来绘制流动的线条（表示分类变量之间的流动）
  geom_stratum(aes( )) + #geom_stratum 用来绘制流动的区块（表示每个类别的大小）
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=4) +
  # scale_fill_igv() +
  scale_fill_d3(palette =  "category20")+ 
  theme_void()+
  theme(text = element_text(size = 13))


ptv3[1:3,1:3]
dat = data.frame(ptv3)
dat = dat[c(sensitive[sensitive$drugname %in% sensitiveDrug, 'samp_ID'], unique(ptv3_info$control))] 
#dat[] = lapply(dat, as.numeric)
dim(dat) # 9143  192
dat$symbol = ptv3unigene[rownames(dat), 'Genes']
dim(dat)

dat$symbol[grepl('NA$',dat$symbol,fixed = T)]

dat = subset(dat,  symbol !='NA' & !is.na(symbol) & symbol !='' )
dat1 = aggregate(dat, by = list(dat$symbol), FUN = function(x){  mean(x, na.rm = T)} )
dim(dat )
rownames(dat1) = dat1$Group.1
#dat1[1:3, 185:194]
dim(dat1)
dat1 = dat1[2: 31]
gsvaPar1 = gsva(expr = as.matrix(dat1) , 
                gset.idx.list = entriz,
                kcdf = "Gaussian" , # ("Gaussian", "Poisson", "none")
                method = "gsva",  # ("gsva", "ssgsea", "zscore", "plage")
                min.sz = 2,
                ssgsea.norm=T,
                verbose=F, 
                parallel.sz = parallel::detectCores())
gsvaPar1 = data.frame(gsvaPar1)


##### end #####