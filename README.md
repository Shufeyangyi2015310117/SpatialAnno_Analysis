# Analysis Code

We leave emply the folders containing inputs and intermediate outputs due to the size limit. The whole files can be downloaded from [zenodo](https://zenodo.org/record/7408584#.Y5Bbpy2d7RY)

## Simulation 
The simulated dataset examples were in ./Simulation/simulated_datasets folder.

Brief descriptions of simulated scripts (./Simulation/Rcode folder):


**simulated_Scenario1_(5,7,9).R**: Annotation analysis for Scenario I. Test the robustness of SpatialAnno to the erroneous specification of the number of cell/domain types.

**simulated_Scenario2_(0.1,0.2,0.3).R**: Annotation analysis for Scenario II. Evaluate the robustness of SpatialAnno to the marker gene misspecification. 

**simulated_Scenario3_ConCor_(5,7,9,0.1,0.2,0.3).R**: Annotation analysis for Scenario III. Assess the capability of SpatialAnno to utilize high-dimensional non-marker genes in all settings by conditional correlation.

**simulated_Scenario3_ARI_(5,7,9,0.1,0.2,0.3).R**: Annotation analysis for Scenario III. Assess the capability of SpatialAnno to utilize high-dimensional non-marker genes in all settings by ARI.

**simulated_Scenario3_Num_nonmarker_(60,100,500,1000,2000).R**: Annotation analysis for Scenario III. Assess the capability of SpatialAnno to utilize increasing numbers of non-marker genes.

## Real data analysis


Brief descriptions of real data analysis scripts (Real_data_analysis folder):

**dorsolateral_prefrontal_cortex.R**: Annotation analysis for human dorsolateral prefrontal cortex Visium data

**mouse_olfactory_bulb.R**: Annotation analysis for mouse olfactory bulb ST data

**olfactory_bulb.R**: Annotation analysis for mouse olfactory bulb Slide-seqV2 data

**mouse_embryo.R**: Annotation analysis for mouse_embryo seqFISH data



## Real data results 
The real data results were visualized with ggplot2 package (Real_data_results folder).

Brief descriptions of real data viualization  scripts:

**dorsolateral_prefrontal_cortex_viualize.R**: Visualization for human dorsolateral prefrontal cortex Visium data

**mouse_olfactory_bulb_viualize.R**: Visualization for mouse olfactory bulb ST data

**olfactory_bulb_viualize.R**:  Visualization for mouse olfactory bulb Slide-seqV2 data

**mouse_embryo.R**: Visualization for mouse_embryo seqFISH data
