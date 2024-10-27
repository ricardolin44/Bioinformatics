import pandas as pd
people = {
    'First' : ['Ricardo', 'Nicholas'],
    'Last' : ['Lin', 'Lilin'],
    'Email' : ['RicardoLin@gmail.com', 'NicholasLilin@gmail.com']
}
df = pd.DataFrame(people)
df.set_index('Email', inplace=False)
df.rename(columns = {'First':'First_name', 'Last':'Last_name'}, inplace=True)
df.loc[1, ['Last_name', 'Email']] = ['Lin', 'NicholasLin@gmail.com']
print(df)