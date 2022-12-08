#!/usr/bin/env python
# coding: utf-8

# In[220]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import anndata


# In[221]:


iter = 11
cutoff=[0.15]


# In[222]:


name_ID12=["151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676"]


# In[223]:


## read raw count
adata = anndata.read_text('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_count.txt', first_column_names = True)


# In[224]:


adata


# In[225]:


df = pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_drsc_celltype.txt', header=None)


# In[226]:


df


# In[227]:


Ez_u =pd.read_csv('./Real_data_results/dataFiles/DLPFC/Brain12cross/dlpfc_'+name_ID12[iter]+'_drsc_Ez.txt', header=None, sep=" ")


# In[228]:


Ez_u


# In[229]:


## assign cell names of raw count to embeddings
Ez_u.index = adata.obs.index


# In[230]:


adata.obsm['Ez_u'] = Ez_u


# In[231]:


adata.obs['clusters'] = df.iloc[:, 0]


# In[232]:


adata.obs['clusters'][0:len(Ez_u)] = df.iloc[0:len(Ez_u), 0]


# In[233]:


adata.obs['clusters']


# In[234]:


adata2 = adata[adata.obs['clusters']!='Unknown', :]
adata = adata2


# In[235]:


adata


# In[236]:


sc.pp.neighbors(adata, n_neighbors=4, use_rep = 'Ez_u')
sc.tl.draw_graph(adata)


# In[237]:


## Let’s use the annotated clusters for PAGA.
sc.tl.paga(adata, groups='clusters')


# In[238]:


sc.pl.paga(adata, threshold=0.07, show=False, save='dlpfc_'+name_ID12[iter]+'_drsc_discard_Unknown.pdf', fontsize=24)


# In[ ]:





# In[ ]:




