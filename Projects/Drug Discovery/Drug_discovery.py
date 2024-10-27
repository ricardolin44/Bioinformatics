import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from chembl_webresource_client.new_client import new_client

def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    baseData = np.arange(1,1)
    i=0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row =np.array([desc_MolWt,
                       desc_MolLogP,
                       desc_NumHDonors,
                       desc_NumHAcceptors])
        if (i==0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnName = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    descriptors = pd.DataFrame(data=baseData, columns=columnName)
    return descriptors

def pIC50(input):
    pIC50 = []
    for i in input['standard_value_norm']:
        pIC50.append(-np.log10(i*10**-9))
    input['pIC50'] = pIC50
    del input['standard_value_norm']
    # input.drop('standard_value_norm')
    return input

def norm_value(input):
    norm = []
    for i in input['standard_value']:
        if float(i)>100000000:
            i = 100000000
        norm.append(float(i))
    input['standard_value_norm'] = norm
    del input['standard_value']
    return input

target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)

selected_target = targets.target_chembl_id[6]

activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50')
df = pd.DataFrame.from_dict(res)

df.to_csv('bioactivity_data.csv', index=False)

df2 = df[df.standard_value.notna()]

bioactivity_class = []
for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append('inactive')
    elif float(i) <= 1000:
        bioactivity_class.append('active')
    else:
        bioactivity_class.append('intermediate')

selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[selection]
bioactivity_class_df = pd.DataFrame(bioactivity_class, columns=['Bioactivity'])
'''We can also use pd.Series if it has only 1 list to be concated'''
df4 = pd.concat([df3, bioactivity_class_df], axis=1)

df4.to_csv('bioactivity_preprocessed_data.csv', index=False)

# list1 = []           Not a smart way to do
# list2 = []
# list3 = []  iterate all of the list that you want
# data_tuples = list(zip(list1,list2,list3,bioactivity_class))
# df3 = pd.DataFrame(data_tuples, columns=['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'bioactivity_class'])

df_lipinski = lipinski(df4.canonical_smiles)
df_combined = pd.concat([df4, df_lipinski], axis=1)
df_combined.describe()
df_norm = norm_value(df_combined)
df_final = pIC50(df_norm)
df_final

'''Plotting'''
# from matplotlib import pyplot as plt
# import seaborn as sns
# sns.set(style='ticks')
# plt.figure(figsize=(5.5,5.5))
# sns.scatterplot(x='MW', y='pIC50', data=df_final, hue='Bioactivity', size='pIC50', alpha=0.7, edgecolor='black')
# plt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0)
# plt.xlabel('MW',fontweight='bold', fontsize=14)
# plt.ylabel('pIC50',fontweight='bold', fontsize=14)
# plt.title('Drug Discovery')
# plt.show()

'''Model building'''
# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.feature_selection import  VarianceThreshold

# x = df_final.drop('pIC50',axis=1)   #can't run since we haven't convert it into descriptors
# y = df_final.pIC50
# selection = VarianceThreshold(threshold=(.8*(1-.8)))   #remove low variance variable
# x = selection.fit_transform(x)

# np.random.seed(100)
# x_train, x_test, y_train, y_test = train_test_split(x,y, test_size=0.2)
# model = RandomForestRegressor(n_estimators=100)
# model.fit(x_train,y_train)
# r2 = model.score(x_test,y_test)

# y_pred = model.predict(x_test)