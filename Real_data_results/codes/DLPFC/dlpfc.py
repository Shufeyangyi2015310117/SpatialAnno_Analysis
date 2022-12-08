#!/usr/bin/env python
# coding: utf-8

# In[1]:


## run paga following https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html#Clustering-and-PAGA
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import anndata


# In[2]:


iter = 1
cutoff=[0.03]


# In[3]:


name_ID12=["151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676"]


# In[4]:


## read raw count
adata = anndata.read_text('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_count.txt', first_column_names = True)


# In[5]:


adata


# In[6]:


df = pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_celltype.txt', header=None)


# In[8]:


df


# In[9]:


Ez_u =pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_Ez.txt', header=None, sep=" ")


# In[10]:


Ez_u


# In[11]:


## assign cell names of raw count to embeddings
Ez_u.index = adata.obs.index


# In[12]:


adata.obsm['Ez_u'] = Ez_u


# In[13]:


adata.obs['clusters'] = df.iloc[:, 0]


# In[14]:


adata.obs['clusters'][0:len(Ez_u)] = df.iloc[0:len(Ez_u), 0]


# In[15]:


adata.obs['clusters']


# In[16]:


adata2 = adata[adata.obs['clusters']!='Unknown', :]
adata = adata2


# In[17]:


adata


# In[18]:


sc.pp.neighbors(adata, n_neighbors=4, use_rep = 'Ez_u')
sc.tl.draw_graph(adata)


# In[19]:


## Let’s use the annotated clusters for PAGA.
sc.tl.paga(adata, groups='clusters')


# In[26]:


sc.pl.paga(adata, threshold=0.3, show=False, save='dlpfc_'+name_ID12[iter]+'_discard_Unknown.pdf', fontsize=24)


# In[ ]:





# In[ ]:




