import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] =  'Arial' 
# matplotlib.rcParams['font.serif'] = 'Times New Roman'

# https://scikit-learn.org/stable/api/sklearn.ensemble.html
from sklearn import svm, datasets
from sklearn import metrics
from sklearn.metrics import auc, RocCurveDisplay, roc_auc_score, roc_curve , accuracy_score, precision_score, f1_score, recall_score,r2_score, mean_squared_error
from sklearn.model_selection import StratifiedKFold, KFold, cross_val_predict, cross_val_score, GridSearchCV
import xgboost as xgb
from sklearn.svm import SVC
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split
# from boruta import BorutaPy
from sklearn.metrics import  mean_absolute_error, mean_squared_error, r2_score
import scipy.stats as stats
import random
import pickle
from sklearn.preprocessing import MinMaxScaler, StandardScaler
mmscale = MinMaxScaler()
stdscale = StandardScaler()


machine_info = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QC/machine_info.xlsx', index_col=0)
print(machine_info['Histopathology_type'].unique())
machine_info.loc[machine_info['Histopathology_type'] == 'PTMC', ['Histopathology_type']] = 'PTC'


QE_deldat_IQR_minImpute = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.9miss_naImpute_matrix.xlsx', index_col=0)

coefficients1 = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_N_glm_ageContinue.xlsx', index_col=0)
coefficientsM = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_M_glm_ageContinue.xlsx', index_col=0)
coefficientsB = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_B_glm_ageContinue.xlsx', index_col=0)
coefficientsBorderline = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/glm_spearman_cor0.1/10yeas_In_one/QEzscore_Discovery_Borderline_glm_ageContinue.xlsx', index_col=0)
padjust = coefficients1[coefficients1['qt_padjust']<0.05].index.values.tolist()
padjustM = coefficientsM[coefficientsM['qt_padjust']<0.05].index.values.tolist()
padjustB = coefficientsB[coefficientsB['qt_padjust']<0.05].index.values.tolist()
padjustBorderline = coefficientsBorderline[coefficientsBorderline['qt_padjust']<0.05].index.values.tolist()

NBpadjust = list(set(padjust + padjustB))

miss60 = pd.read_excel('D:/chh/2023workProject/20240821TTTD/QE/thyroid_QE_than2000_0.6miss_pro.xlsx')
# miss60['miss60'].values.tolist()
proteins = list(set(padjust) & set (miss60['miss60'].values ))
print(len(padjust), len(proteins))


dat480_deldat_min08 = pd.read_excel('D:/chh/2023workProject/20240821TTTD/thyroid_480_than2000_0.9miss_naImpute_matrix.xlsx', index_col=0)
info480 = dat480_deldat_min08.iloc[:,:5]
matrix480 = dat480_deldat_min08.iloc[:,5:]
featzscore = stdscale.fit_transform(matrix480)
matrix480 = pd.DataFrame(featzscore, columns=matrix480.columns.values, index=matrix480.index.values)
matrix480


matrix = QE_deldat_IQR_minImpute.iloc[:,5:]
info = QE_deldat_IQR_minImpute.iloc[:,:5]
featzscore = stdscale.fit_transform(matrix)
matrix = pd.DataFrame(featzscore, columns=matrix.columns.values, index=matrix.index.values)
matrix['Gender'] = pd.Categorical(info['Gender']).codes
matrix['SampleType'] = pd.Categorical(info['SampleType']).codes
samples = set(machine_info[(machine_info['Classificaiton_type'].isin(['N'])) & (machine_info['TrainTest'] == 'Discovery') & (~pd.isna(machine_info['Gender'] ))].index.values  ) & set(QE_deldat_IQR_minImpute.index.values)
samples = list(samples)
print(samples[:10] )

# matrix['Gender'] = pd.Categorical(info['Gender']).codes
# matrix['Hospital'] = pd.Categorical(info['Hospital']).codes
# matrix['SampleType'] = pd.Categorical(info['SampleType']).codes

trainDat = matrix.loc[samples][padjust]  # padjust+[] NtestProsadjust 'Hospital', padjust, 'SampleType','Gender'
trainY = info.loc[samples]['Age']
print(trainDat.shape)


samplesT1 = set(machine_info[(machine_info['Classificaiton_type'].isin(['N', 'B']) ) & (machine_info['TrainTest'] == 'Test1') & (~pd.isna(machine_info['Gender'] ))].index.values  ) & set(QE_deldat_IQR_minImpute.index.values)
samplesT1 = list(samplesT1)
print(len(samplesT1), samplesT1[:10] )

testDat1 = matrix.loc[samplesT1]#[padjust ]#+['SampleType', 'Hospital',	'Gender']
testY1 = info.loc[samplesT1]['Age']

########### feat select
fold_num = 3
param_grid = { 'max_depth': [None,1 , 2, 3, 4 ],'min_samples_leaf': [1, 2, 4, 6], 'min_samples_split': [2, 4, 6, 8], "n_estimators":[ 500, 1000]} # 20, 30, 50, 100, 700, 800,
# neg_mean_squared_error, r2, neg_root_mean_squared_error, neg_mean_absolute_error
grid_search = GridSearchCV(RandomForestRegressor(random_state = 2024), param_grid, cv = fold_num, scoring='r2', n_jobs = -1 ) # class_weight='balanced', 
grid_search.fit(trainDat, trainY)
params = grid_search.best_params_
print( "RandomForestRegressor", fold_num,  "Best params:{}".format(grid_search.best_params_), "Best score on trainset:{:.3f}".format(grid_search.best_score_))
n_estimators = params["n_estimators"]
min_samples_leaf = params["min_samples_leaf"]
min_samples_split = params["min_samples_split"]
max_depth = params["max_depth"]
modelselect = RandomForestRegressor(n_estimators = n_estimators, min_samples_leaf = min_samples_leaf, min_samples_split = min_samples_split , max_depth = max_depth , random_state = 2024)
# score = model.score(testDat1, testY1)
# print('test1', score)

#### shap value
import shap
#model.fit(trainDat, trainY)
explainer1 = shap.TreeExplainer(modelselect)
shap_values1 = explainer1.shap_values(trainDat )
shap.summary_plot(shap_values1, trainDat , show=False, max_display=10)# , max_display=30
# plt.savefig("D:/chh/2023workProject/20240821TTTD/ML/QE_N"+ str(fold_num) + "fold_top13_shapvalues20241021.pdf")

mean_abs_shap_values = np.mean(np.abs(shap_values1), axis=0)
feats = trainDat.iloc[:, mean_abs_shap_values > 0].columns.values.tolist()
shapfeat = {}
for i in feats:
    index = trainDat.columns.values.tolist().index(i)
    shapfeat[i] = mean_abs_shap_values[index]
shapfeat = sorted(shapfeat.items(),  key=lambda x: x[1], reverse=True)
shapfeat = [i[0] for i in shapfeat ]
print(len(shapfeat), shapfeat[:100])
print(mean_abs_shap_values.shape, trainDat.shape)


### train
models = [ "xgb" ] # "RandomForestRegressor","DecisionTreeClassifier", "KNeighborsClassifier", "SVR",
fold_num = 3
# featimportce
#top_features_ = rf_feature_select# ['SampleType', 'Hospital',	'Gender']+ [ 'Hospital', 'Gender']+
top_features_ = shapfeat# ['SampleType', 'Hospital',	'Gender']+ [ 'Hospital', 'Gender']+

featNum = 100
ml_train = trainDat
print(top_features_[:featNum],'\n')
print(featNum,'Features')
message = "{}\n\n".format(top_features_[:featNum])
# plt.figure()

for featNum in range(105,201, 5): # [45]
    print(featNum,'Features')
    m = "xgb"
# for m in models:
    ml_train = trainDat.reset_index(drop=True)
    y_train = trainY.reset_index(drop=True)

    if m == "RandomForestRegressor":
        # param_grid = { 'max_depth': [None, 2, 3 ],'min_samples_leaf': [2, 4], 'min_samples_split': [2, 4, 8], "n_estimators":[20, 30, 50, 100, 500, 1000]}
        param_grid = { 'max_depth': [None,1 , 2, 3, 4 ],'min_samples_leaf': [1, 2, 4, 6], 'min_samples_split': [2, 4, 6, 8], "n_estimators":[ 500, 700, 1000]} # 20, 30, 50, 100,
        grid_search = GridSearchCV(RandomForestRegressor(random_state = 2024), param_grid, cv = fold_num, scoring='neg_mean_absolute_error', n_jobs = -1) 
        grid_search.fit(trainDat[top_features_[:featNum]], trainY )
        params = grid_search.best_params_
        print("RandomForestRegressor", "Best params:{}".format(grid_search.best_params_), "Best score on trainset:{:.3f}".format(grid_search.best_score_))
        message = message + "RandomForestRegressor Best params:{}, Best score on trainset:{:.3f} \n".format(grid_search.best_params_, grid_search.best_score_)
        ntree = params["n_estimators"]
        max_depth = params["max_depth"]
        min_samples_leaf = params["min_samples_leaf"]
        min_samples_split = params["min_samples_split"]
        model = RandomForestRegressor(n_estimators = ntree, min_samples_leaf = min_samples_leaf, min_samples_split = min_samples_split , max_depth = max_depth , random_state = 2024, n_jobs = -1)
          
    elif m=="SVR":
        param_grid = {'kernel': ['rbf'], 'gamma': [5e-3, 1e-3, 5e-4, 1e-4], 'C': [1, 10, 100, 1000]}
        # param_grid = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [1, 10, 100, 1000]},
        #             {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]
        grid_search = GridSearchCV(SVR(random_state=2024), param_grid, cv = fold_num, scoring='neg_mean_absolute_error', n_jobs = -1) 
        grid_search.fit(ml_train[top_features_[:featNum]], y_train )
        params = grid_search.best_params_
        print("SVR Best params:{}".format(grid_search.best_params_), "Best score on trainset:{:.3f}".format(grid_search.best_score_))
        message = message + "SVR Best params:{} , Best score on trainset:{:.3f} \n".format(grid_search.best_params_, grid_search.best_score_)
        kernel = params["kernel"]
        gamma = params["gamma"]
        C = params["C"]
        model = SVR(C = C, kernel = kernel, gamma = gamma, random_state = 2024)

    elif m=="xgb":
        param_grid = {'max_depth': [1,2,3,4,5], 'gamma': [5e-3, 1e-3, 5e-4, 1e-4, 0.3, 0.5, 0.7, 0.9, 1], 'n_estimators': [500, 1000], # 20, 50, 100, 200,
        'learning_rate':[0.1, 0.05, 0.03 , 0.01,0.001] , 'subsample':[0.3, 0.4, 0.5, 0.7, 0.9], 'reg_alpha' :[1, 2, 4, 8, 10]
        } # 'reg_lambda' :[0.25, 0.2, 0.26 ], 
        grid_search = GridSearchCV(xgb.XGBRegressor( objective='reg:squarederror', random_state=2024), param_grid, cv = fold_num, scoring='neg_mean_absolute_error', n_jobs = -1) # objective ='reg:linear', neg_mean_absolute_error
        grid_search.fit(trainDat[top_features_[:featNum]], y_train )
        params = grid_search.best_params_
        print("xgb Best params:{}".format(grid_search.best_params_), "Best score on trainset:{:.3f}".format(grid_search.best_score_))
        message = message + "xgb Best params:{} , Best score on trainset:{:.3f} \n".format(grid_search.best_params_, grid_search.best_score_)
        max_depth = params["max_depth"]
        gamma = params["gamma"]
        n_estimators = params["n_estimators"]
        learning_rate = params["learning_rate"]
        # reg_alpha = params['reg_alpha']
        # subsample = params["subsample"]
        model = xgb.XGBRegressor(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators, # subsample = subsample,reg_alpha = reg_alpha, 
        objective='reg:squarederror', gamma=gamma, random_state=2024)

    X1 = ml_train[top_features_[:featNum]].reset_index(drop=True)
    kfold = KFold(n_splits = fold_num , shuffle=True, random_state = 2024 )# 设置 K-Fold 参数
    MAE = []# 初始化用于存储每次交叉验证结果的列表
    MSE = []
    RMSE = []
    r2 = []
    probabilitys = []
    y_labels = []

    for i, (train_index, test_index) in enumerate(kfold.split(X1 )):
        X_train, X_val = ml_train[top_features_[:featNum]].loc[train_index], ml_train[top_features_[:featNum]].loc[test_index]
        ml_y_train, y_val = y_train.loc[train_index], y_train.loc[test_index]
        if len( y_train.loc[train_index].unique())<2 or len( y_train.loc[test_index].unique())<2:
            print("A little difference in lable, {} fold validation can't be done.".format(fold_num))
            message = message + "A little difference in lable, {} fold validation can't be done. \n".format(fold_num)
        model.fit(X_train, ml_y_train)# 训练模型
        
        y_pred = model.predict(X_val)

        MAE.append(mean_absolute_error(y_val, y_pred))
        MSE.append(mean_squared_error(y_val, y_pred))
        r2.append(r2_score(y_val, y_pred))
        # np.sqrt(metrics.mean_squared_error(y_test, y_pred)))

    mean_mae = np.mean(MAE)
    mean_mse = np.mean(MSE)
    mean_r2 = np.mean(r2)
    print('Train: mae={}, mse={}, r2={}'.format(round(mean_mae,3), round(mean_mse,3), round(mean_r2,3) ))

    test_pred = model.predict(testDat1[top_features_[:featNum]])
    mae1 = mean_absolute_error(testY1.values, test_pred)
    mse1 = mean_squared_error(testY1.values, test_pred)
    r2_1 = r2_score(testY1.values, test_pred)
    r2weight_1 = r2_score(testY1.values, test_pred,multioutput='variance_weighted')
    corR = stats.spearmanr(testY1.values, test_pred)[0]
    corPearsonR = stats.pearsonr(testY1.values, test_pred)[0]
    print('Test1: mae={}, mse={}, r2={}, pearsonr={}, spearmanr={}'.format(round(mae1,3), round(mse1,3), round(r2_1,3), round(corPearsonR,3), round(corR,3)))
    print(test_pred.tolist(),'\n' )

    # test_pred = model.predict(testDat2[top_features_[:featNum]])
    # mae2 = mean_absolute_error(testY2.values, test_pred)
    # mse2 = mean_squared_error(testY2.values, test_pred)
    # r2_2 = r2_score(testY2.values, test_pred)
    # r2weight_2 = r2_score(testY2.values, test_pred,multioutput='variance_weighted')
    # corR = stats.spearmanr(testY2.values, test_pred)[0]
    # print('Test2: mae={}, mse={}, r2={}, r2weight={}, spearmanr={}'.format(mae2, mse2, r2_2, r2weight_2, corR))
    # print(test_pred.tolist(),'\n')

    message = message + 'Test1: mae1={}, mse={}, r2={}, spearmanr={}\n\n'.format(mae1, mse1, r2_1, corR)
    # with open('D:/chh/2023workProject/20240821TTTD/ML/QE_Nzscore_zscore_xgb_shap_'+str(featNum)+'features20241021.pkl', 'wb') as m:
    with open('D:/chh/2023workProject/20240821TTTD/ML/QE_NBzscore_zscore_xgb_shap_'+str(featNum)+'features20241024.pkl', 'wb') as m:
        pickle.dump(model, m)


# QEzscore RFselect XGBmae shap
top_features_ = ['Q8NI22_MCFD2', 'Q96A11_GAL3ST3', 'P60033_CD81', 'P26440_IVD', 'P02743_APCS', 'O14639_ABLIM1', 'Q9P1F3_ABRACL', 'O94811_TPPP', 'Q9UHL4_DPP7', 'Q8NBI5_SLC43A3', 'P02794_FTH1', 'P30038_ALDH4A1', 
'Q86WU2_LDHD', 'P01011_SERPINA3', 'O75503_CLN5', 'P06727_APOA4', 'Q9UNL2_SSR3', 'O60888_CUTA', 'P0DJI8_SAA1', 'Q06136_KDSR', 'Q9Y624_F11R', 'P04004_VTN', 'Q92743_HTRA1', 'P06703_S100A6', 'P01023_A2M', 
'Q8TAE6_PPP1R14C', 'P10253_GAA', 'Q06481_APLP2', 'P12830_CDH1', 'P08603_CFH', 'Q16134_ETFDH', 'O43684_BUB3', 'P37840_SNCA', 'P11182_DBT', 'P40121_CAPG', 'Q96DG6_CMBL', 'P01903_HLA.DRA', 'P59044_NLRP6', 
'P61626_LYZ', 'Q9H910_JPT2', 'Q96SI9_STRBP', 'Q96RS6_NUDCD1', 'P51858_HDGF', 'Q15276_RABEP1', 'Q14746_COG2', 'P52655_GTF2A1', 'P00488_F13A1', 'Q9P0M6_MACROH2A2', 'O95810_CAVIN2', 'P04439_HLA.A', 'P10768_ESD', 
'P21399_ACO1', 'O75955_FLOT1', 'P37802_TAGLN2', 'P62979_RPS27A', 'P52907_CAPZA1', 'Q08380_LGALS3BP', 'Q12841_FSTL1', 'Q96EE3_SEH1L', 'P00325_ADH1B', 'P02787_TF', 'P83876_TXNL4A', 'P20700_LMNB1', 'Q12906_ILF3',
'Q8IW45_NAXD', 'P80723_BASP1', 'P43251_BTD', 'Q92777_SYN2', 'Q9NUP9_LIN7C', 'P50226_SULT1A2', 'P22748_CA4', 'P07602_PSAP', 'O75173_ADAMTS4', 'P02671_FGA', 'Q8IYI6_EXOC8', 'P43243_MATR3', 'O43521_BCL2L11',
'P30711_GSTT1', 'P61956_SUMO2', 'Q9UKK3_PARP4', 'Q96FZ7_CHMP6', 'Q16822_PCK2', 'P02790_HPX', 'Q9BXN1_ASPN', 'P06576_ATP5F1B', 'Q92522_H1.10', 'P05155_SERPING1', 'P15088_CPA3', 'P05388_RPLP0', 'P07942_LAMB1',
'Q96I24_FUBP3', 'Q99523_SORT1', 'Q9Y3B2_EXOSC1', 'O15439_ABCC4', 'P53582_METAP1', 'P01876_IGHA1', 'P54289_CACNA2D1', 'Q5HYI8_RABL3', 'Q8NE62_CHDH', 'Q00796_SORD']
model = pickle.load(open( "D:/chh/2023workProject/20240821TTTD/ML/matrixzscore_Nzscore/QE_Nzscore_zscore_xgb_shap_95features20241021.pkl", 'rb'))
predfeatnum = 95
model
########## predict
predictresult = info
predictresult['Pred'] = [np.nan] * 6149

# samplesT1 = set(machine_info[(machine_info['Classificaiton_type'].isin(['N', 'B', 'M' 'Borderline'])) & (machine_info['TrainTest'].isin(['Discovery', 'Test1', "Test3"])) & (~pd.isna(machine_info['Gender'] ))].index.values ) & set(QE_deldat_IQR_minImpute.index.values)
#  & (~pd.isna(machine_info['Gender'] ))
samplesT1 = set(machine_info[ (machine_info['TrainTest'].isin(['Test1','Discovery', 'Test2',"Test3"]))].index.values  ) & set(QE_deldat_IQR_minImpute.index.values)
samplesT1 = set(machine_info.index.values  ) & set(QE_deldat_IQR_minImpute.index.values)

samplesT1 = list(samplesT1)
print(len(samplesT1), samplesT1[:10] )

testDatB = matrix.loc[samplesT1 ]#[padjust ]#+['SampleType', 'Hospital',	'Gender']
testYB = info.loc[samplesT1 ]['Age']

test_pred = model.predict(testDatB[top_features_[:predfeatnum]])
mae1 = mean_absolute_error(testYB.values, test_pred)
mse1 = mean_squared_error(testYB.values, test_pred)
r2_1 = r2_score(testYB.values, test_pred)
r2weight_1 = r2_score(testYB.values, test_pred,multioutput='variance_weighted')
corR = stats.spearmanr(testYB.values, test_pred)[0]
corPearsonR = stats.pearsonr(testYB.values, test_pred)[0]
# print(test_pred.tolist())
predictresult.loc[samplesT1 ,'Pred'] = test_pred.tolist()

print('Test1: mae={}, mse={}, r2={}, pearsonr={}, spearmanr={}'.format(round(mae1,3), round(mse1,3), round(r2_1,3), round(corPearsonR,3), round(corR,3)))
# print(testYB.tolist(),'\n' )
# print(test_pred.tolist(),'\n' )
predictresult.iloc[:, :6].to_excel( 'D:/chh/2023workProject/20240821TTTD/ML/QE_xgbmae_shap_'+str(predfeatnum)+'feat_predict20241031.xlsx' )

### QETest3
QETest3 = pd.read_excel("D:/chh/2023workProject/20240821TTTD/QE/QETest3_log2min0.8.xlsx", index_col=0)
trainDatzscore = (QETest3 - QETest3.mean()) / QETest3.std()
na_columns = trainDatzscore.columns[trainDatzscore.isna().any()].tolist()
trainDatzscore = trainDatzscore.drop(columns=na_columns)
predictresult = machine_info.loc[trainDatzscore.index][['Age', 'SampleType', 'DataSet',	'Hospital',	'Gender', 'TrainTest']]
predictresult['Pred'] = [np.nan] * predictresult.shape[0]
# & (~pd.isna(machine_info['Gender'] ))
samplesT1 = set(machine_info[(machine_info['TrainTest'].isin([ 'Test3'])) ].index.values ) & set(trainDatzscore.index.values)
samplesT1 = list(samplesT1)
print(len(samplesT1), samplesT1[:10] )

testDatB = trainDatzscore.loc[samplesT1 ]#[padjust ]#+['SampleType', 'Hospital',	'Gender']
testYB = predictresult.loc[samplesT1 ]['Age']

test_pred = model.predict(testDatB[top_features_[:predfeatnum]])
mae1 = mean_absolute_error(testYB.values, test_pred)
mse1 = mean_squared_error(testYB.values, test_pred)
r2_1 = r2_score(testYB.values, test_pred)
r2weight_1 = r2_score(testYB.values, test_pred,multioutput='variance_weighted')
corR = stats.spearmanr(testYB.values, test_pred)[0]
corPearsonR = stats.pearsonr(testYB.values, test_pred)[0]
# print(test_pred.tolist())
predictresult.loc[samplesT1 ,'Pred'] = test_pred.tolist()
predictresult.to_excel( 'D:/chh/2023workProject/20240821TTTD/ML/QEtest3_xgbmae_shap_'+str(predfeatnum)+'feat_predict20241021.xlsx' )

### 480
predictresult = machine_info.loc[matrix480.index][['Age', 'SampleType', 'DataSet', 'Hospital', 'Gender', 'TrainTest']]
print(predictresult.columns.values)
print(predictresult['TrainTest'].unique())
predictresult
matrix480[['Q9P1F3_ABRACL', 'Q8NBI5_SLC43A3', 'P59044_NLRP6', 'Q9P0M6_MACROH2A2', 'Q12841_FSTL1', 'Q92777_SYN2', 'Q9NUP9_LIN7C', 'Q9UKK3_PARP4', 'Q9Y3B2_EXOSC1']] = 0
# matrix480
predictresult['Pred'] = [np.nan] * 2975
# & (~pd.isna(machine_info['Gender'] ))
samplesT1 = set(machine_info[(machine_info['TrainTest'].isin(['Discovery','Test1', 'Test2', 'Test3'])) ].index.values ) & set(matrix480.index.values)
samplesT1 = list(samplesT1)
print(len(samplesT1), samplesT1[:10] )

testDatB = matrix480.loc[samplesT1 ]#[padjust ]#+['SampleType', 'Hospital',	'Gender']
testYB = predictresult.loc[samplesT1 ]['Age']

test_pred = model.predict(testDatB[top_features_[:predfeatnum]])
mae1 = mean_absolute_error(testYB.values, test_pred)
mse1 = mean_squared_error(testYB.values, test_pred)
r2_1 = r2_score(testYB.values, test_pred)
r2weight_1 = r2_score(testYB.values, test_pred,multioutput='variance_weighted')
corR = stats.spearmanr(testYB.values, test_pred)[0]
corPearsonR = stats.pearsonr(testYB.values, test_pred)[0]
# print(test_pred.tolist())
predictresult.loc[samplesT1 ,'Pred'] = test_pred.tolist() # N B M Borderline

print('Test1: mae={}, mse={}, r2={}, pearsonr={}, spearmanr={}'.format(round(mae1,3), round(mse1,3), round(r2_1,3), round(corPearsonR,3), round(corR,3)))
# print(testYB.tolist(),'\n' )
# print(test_pred.tolist(),'\n' )
predictresult.to_excel( 'D:/chh/2023workProject/20240821TTTD/ML/480_xgbmae_shap_'+str(predfeatnum)+'feat_predict20241021.xlsx' )


dat480Test3 = pd.read_excel("D:/chh/2023workProject/20240821TTTD/dat480/dat480Test3_log2min0.8.xlsx", index_col=0)
trainDatzscore = (dat480Test3 - dat480Test3.mean()) / dat480Test3.std()
na_columns = trainDatzscore.columns[trainDatzscore.isna().any()].tolist()
trainDatzscore = trainDatzscore.drop(columns=na_columns)
trainDatzscore['P59044_NLRP6'] = 0
predictresult = machine_info.loc[trainDatzscore.index][['Age', 'SampleType',	'DataSet',	'Hospital',	'Gender', 'TrainTest']]
predictresult['Pred'] = [np.nan] * predictresult.shape[0]
# & (~pd.isna(machine_info['Gender'] ))
samplesT1 = set(machine_info[(machine_info['TrainTest'].isin([ 'Test3'])) ].index.values ) & set(trainDatzscore.index.values)
samplesT1 = list(samplesT1)
print(len(samplesT1), samplesT1[:10] )

testDatB = trainDatzscore.loc[samplesT1 ]#[padjust ]#+['SampleType', 'Hospital',	'Gender']
testYB = predictresult.loc[samplesT1 ]['Age']

test_pred = model.predict(testDatB[top_features_[:predfeatnum]])
mae1 = mean_absolute_error(testYB.values, test_pred)
mse1 = mean_squared_error(testYB.values, test_pred)
r2_1 = r2_score(testYB.values, test_pred)
r2weight_1 = r2_score(testYB.values, test_pred,multioutput='variance_weighted')
corR = stats.spearmanr(testYB.values, test_pred)[0]
corPearsonR = stats.pearsonr(testYB.values, test_pred)[0]
# print(test_pred.tolist())
predictresult.loc[samplesT1 ,'Pred'] = test_pred.tolist()
predictresult.to_excel( 'D:/chh/2023workProject/20240821TTTD/ML/480test3_xgbmae_shap_'+str(predfeatnum)+'feat_predict20241021.xlsx' )