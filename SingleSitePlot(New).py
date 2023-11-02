# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 14:38:29 2022

@author: Hasan
"""
import pandas as pd
import time
from datetime import datetime
from functools import reduce
import sys
sys.path.insert(1, 'C:/Users/Hasan/HKUST/Haibin SU - group - covid19/SingleSiteAnalysis')
import ImportantFunc as Imp
from operator import itemgetter
start_time = time.process_time()
print("Start =", datetime.now().strftime("%H:%M:%S"))
  
def getList(dict):
    return list(map(itemgetter(0), dict.items()))
def mut_info(MutList, dfTypeAA):
    MutEle = []
    for Mut in MutList:
        for ele in Mut.split(";"): MutEle.append(ele)
    MutEle = sorted(MutEle)
    M = Imp.count_dups(MutEle)
    df = pd.DataFrame({'Mutation': M[0],'Count': M[1]})
    df['Ori'] = [x[0] for x in df['Mutation']]
    df['MutAA'] = [x[-1] for x in df['Mutation']]
    df['Pos'] = df['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    df['WL-Count'] = df['MutAA'].astype(str) + "(" + df['Count'].astype(str) + ")"
    df['Sites'] = df['Ori'].astype(str) + df['Pos'].astype(str)
    df['Count'].astype(int), df['Pos'].astype(int)
    df = pd.merge(df, dfTypeAA, on ='MutAA',how ='left')#Merging the information about AA Type
    df.sort_values(['Pos','Count'],inplace = True, ascending=[True, False])
    df = df.reset_index(drop=True)
    return df
def WL(df,oriAA):
    df["space"] = ""
    R = df["Sites"].tolist()
    N_AA = Imp.count_dups( R ) #N.AA
    N_Count = Imp.sumsimilar( R, df["Count"].tolist() )#N.Count
    WL_AA= Imp.concat_dups( R, df["MutAA"].tolist() ) #WL_AA
    WL_Count = Imp.concat_dups( R, df["WL-Count"].tolist() ) #WL_count    
    df3 = pd.DataFrame({'Sites': oriAA})
    df3['Pos'] = df3['Sites'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    df2 = pd.DataFrame({'Sites': N_AA[0],'N_AA': N_AA[1]})
    df30 = pd.DataFrame({'Sites': N_Count[0],'N_Count': N_Count[1]})
    df31 = pd.DataFrame({'Sites': WL_AA[0],'WL_AA': WL_AA[1]})
    df32 = pd.DataFrame({'Sites': WL_Count[0],'WL_Count': WL_Count[1]})
    dfs = [df3,df2,df30,df31,df32]
    df3 = reduce(lambda left,right:pd.merge(left,right,on ='Sites',how ='outer'), dfs).fillna(0)
    d = {}
    AAType = ["H_Aliphatic","H_Aromatic","Polar","N_Polar","Positive","Negative","Del","Special"]
    for i in AAType:
        d[i] = df.filter(['Sites','Count','MutType'], axis=1).loc[df['MutType'] == i]
        if len(d[i]) == 0:
            d[i] = pd.DataFrame({'Sites':oriAA,i:len(oriAA)*[0]})
            continue
        d[i] = d[i].drop('MutType', axis=1)
        d[i] = d[i].groupby(['Sites']).sum().reset_index(level=['Sites'])
        d[i] = d[i].reset_index(drop=True).rename({'Count':i},axis = 1)
    
    for i0 in getList(d):
        df3 = pd.merge(df3,d[i0],on ='Sites',how ='left').fillna(0)
    
    columns = ['N_AA','WL_AA','N_Count','WL_Count']
    columns.extend(AAType)
    d0 = {}
    for j in columns:
        d0[j] = df3[['Sites', j]].copy()
    return df,df3,d0
###############################################################################################
MonInput = "0304"
dfInput = pd.read_excel('AllUnique_'+MonInput+'-(Corrected).xlsx').dropna(subset=['mutation info']).reset_index(drop=True)
dfTypeAA = pd.read_excel('OriRefSpike.xlsx', sheet_name ='TypeAA', usecols = 'A,B')
dfVarRef = pd.read_excel('OriRefSpike.xlsx', sheet_name ='VarRef', header = 0,usecols = [0,1])
print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's\n')

ref0 = dfVarRef["SeqMSA"].tolist()[0]
oriAA = [str(ref0[m]+str(m+1)) for m in range(len(ref0))]

dfSites = pd.DataFrame({'Sites': oriAA})
D0 = ['N_AA', 'WL_AA', 'N_Count', 'WL_Count', 'H_Aliphatic', 'H_Aromatic', 'Polar', 'N_Polar', 'Positive', 'Negative', 'Del', 'Special']
Df = {}
for Di in D0:
    Df[Di] = dfSites.copy(deep=True)

DV = ['mutation number','Var']
DevVar = {}
for Di2 in DV:
    DevVar[Di2] = pd.DataFrame({Di2:[]})
###############################################################################################
N = 38 # The Number of months
N2 = []
###############################################################################################
for z in range(N,N+1):
    dfInpMon = dfInput.loc[(dfInput['MonthIndex'] <= z)].reset_index(drop=True) #&(dfInput['VarIndex'] == y)
    MutList = dfInpMon['mutation info'].tolist()
    print('iteration -{0} Total Seq = {1}'.format(z,len(MutList)) ) # Variant = {2}  ,y
    if MutList == []: continue
    N2.append(z)
    DVList = {}
    for k in DV:
        DVList[k] = dfInpMon[k].tolist()
    ############################################################################################### 
    """ Deviation & Variant Counting """
    for Di2 in DV:
        E = Imp.count_dups(DVList[Di2])
        dfDev = pd.DataFrame({Di2: E[0],z: E[1]})
        DevVar[Di2] = pd.merge(DevVar[Di2], dfDev, on =Di2,how ='outer')
    ###############################################################################################
    df = mut_info(MutList,dfTypeAA) #PostProcessingDataFrame
    df,dfFinal,d0 = WL(df,oriAA)
    ############################################################################################### 
    """ Creating dfNAA, dfWL_AA, dfNCount, dfWL_Count """
    for j in D0:
        d0[j] = d0[j].rename(columns={j:z})
        Df[j] = pd.merge(Df[j],d0[j],on='Sites',how ='left')
###############################################################################################
""" Generating Mutation Trend """
DomName = ['NTD', 'RBD', 'SD12', 'FC', 'FP', 'S2']
DomRange = [[1,324],[325,540],[541,650],[65,700],[701,850],[851,1273]]
AAType = ["H_Aliphatic","H_Aromatic","Polar","N_Polar","Positive","Negative","Del","Special"]
Ddom = {}
for x in range(len(DomName)):
    DN = DomName[x]
    DR = DomRange[x]
    Ddom[DN] = pd.DataFrame({'Months':N2})
    DAA = []
    for z0 in AAType:
        Ddom[DN][z0] = Df[z0].iloc[DR[0]-1:DR[1], 1:].sum(axis = 0, skipna = True).tolist()       
""" Print Out Basic Information: dfNAA, dfWL_AA, dfNCount, dfWL_Count """
DFinal = pd.DataFrame({'space0':[]})
for k in D0:
    Df[k]['space'] = ""
    Df[k] = Df[k].rename(columns={'Sites':k})
    DFinal = pd.concat([DFinal,Df[k]],axis=1)
###############################################################################################
""" Print Out Deviation,Variant Counting, and PostProcessing"""
DFinal2 = pd.DataFrame({'space0':[]})
for Di2 in DV:
    DevVar[Di2].sort_values([Di2],inplace = True, ascending=[True])
    DevVar[Di2] = DevVar[Di2].reset_index(drop=True)
    DevVar[Di2] = DevVar[Di2].fillna(0)
    DevVar[Di2]['space'] = ""
    DFinal2 = pd.concat([DFinal2,DevVar[Di2]],axis=1)
###############################################################################################
# FileProcessing = "SingleSites_TE_{}_{}-Months.xlsx".format(MonInput,z)
FileProcessing = "SingleSites_TE_{}_All-Months.xlsx".format(MonInput)
DFinal.drop('space0',inplace=True,axis=1), DFinal2.drop('space0',inplace=True,axis=1)
with pd.ExcelWriter(FileProcessing) as writer:
    DFinal.to_excel(writer, sheet_name='BasicData', index=None)
    DFinal2.to_excel(writer, sheet_name='Dev-Var')
    pd.concat([df,dfFinal],axis=1).to_excel(writer, sheet_name='PostProcessing', index=None)
    for i in DomName:
        Ddom[i].to_excel(writer, sheet_name=i, index=None)
print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's\n')