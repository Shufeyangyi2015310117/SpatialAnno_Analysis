#!/usr/bin/env python
# coding: utf-8

# In[274]:


## run paga following https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html#Clustering-and-PAGA

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import anndata


# In[275]:


iter = 11
cutoff=[0.15]


# In[276]:


name_ID12=["151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676"]


# In[277]:


## read raw count
adata = anndata.read_text('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_count.txt', first_column_names = True)


# In[278]:


adata


# In[279]:


df = pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_byesSpace_celltype.txt', header=None)


# In[280]:


df


# In[281]:


Ez_u =pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_byesSpace_pc.txt', header=None, sep=" ")


# In[282]:


Ez_u


# In[283]:


## assign cell names of raw count to embeddings
Ez_u.index = adata.obs.index


# In[284]:


adata.obsm['Ez_u'] = Ez_u


# In[285]:


adata.obs['clusters'] = df.iloc[:, 0]


# In[286]:


adata.obs['clusters'][0:len(Ez_u)] = df.iloc[0:len(Ez_u), 0]


# In[287]:


adata.obs['clusters']


# In[288]:


adata2 = adata[adata.obs['clusters']!='Unknown', :]
adata = adata2


# In[289]:


adata


# In[290]:


sc.pp.neighbors(adata, n_neighbors=4, use_rep = 'Ez_u')
sc.tl.draw_graph(adata)


# In[291]:


## Let’s use the annotated clusters for PAGA.
sc.tl.paga(adata, groups='clusters')


# In[292]:


sc.pl.paga(adata, threshold=0.6, show=False, save='dlpfc_'+name_ID12[iter]+'_bayeSpace_discard_Unknown.pdf', fontsize=24)


# In[ ]:




