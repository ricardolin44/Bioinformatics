import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt

df = pd.read_csv('data/iplauction2023.csv')
df.dropna(inplace=True)
df_new = df[df['base price (in lacs)']>100]
# sns.set(rc={'figure.figsize':(11.7,8.27)})
# g = sns.histplot(x='name',y='base price (in lacs)', data=df_dropna[:10])
df_new


_, mean, std, *_ =  df['base price (in lacs)'].describe()
print(mean, std)

df['revenue'] = df['base price (in lacs)'].apply(lambda x: x/100)
df.columns
df.info()
df.describe()
df.describe(include=['object','bool'])
df['final price (in lacs)'].value_counts(normalize=True)
df.loc[0:10, 'name':'nationality']
df[df['nationality'].apply(lambda x: x[0]=='I')].head()

df.status.value_counts()

d = {'SOLD':'Bye', 'RETAINED':'Welcome'}
df['status'] = df['status'].map(d)
df.head()

col_to_show = ['base price (in lacs)','final price (in lacs)']
df.groupby(['nationality'])[col_to_show].describe(percentiles=[])

# df.pivot_table([col],[line], aggfunc='mean')

diff = (df['final price (in lacs)'] - df['base price (in lacs)'])
df.insert(loc=len(df.columns), column = 'Price differences', value=diff)
df
df = df.rename(columns={'Price difference' : 'Price sum up'})

