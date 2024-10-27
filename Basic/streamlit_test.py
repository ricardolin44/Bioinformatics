# from Bio import SeqIO
# from Bio import Entrez

# Entrez.email = 'ricardo.lin100@gmail.com'
# handle = Entrez.efetch(db='nucleotide', id='NM_001300741', rettype='gb')
# nuc_data = SeqIO.read(handle,'genbank')

# nuc_data

# import pyLDAvis.gensim_models

# pyLDAvis.enable_notebook()
# vis = pyLDAvis.gensim_models.prepare(lda, corpus, dictionary)
# pyLDAvis.show(vis, local=False)

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

data = pd.read_csv('bioactivity_data.csv')

st.title('Streamlit Tester')

st.markdown('''
This app is used for testing Streamlit App
* **Python libraries:** pandas, streamlit
* **Data source:** [Google.com](https://www.google.com)
''')

st.sidebar.header('User Input Features')
selected_year = st.sidebar.selectbox('Year', list(reversed(range(1950,2020))))

st.header('Display Stats of the Data(s)')
st.write('Data Dimension: '+ str(data.shape[0]) + ' rows and '+str(data.shape[1]))
st.dataframe(data)

