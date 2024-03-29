# Seperation of Hippocampus dorsal and ventral 


#### This project has successfully achieved 
#### the separation of the dorsal and ventral parts 
#### of the subregion from the hippocampus
#### specifically CA1, CA3, DG, and Suniculum. 

## Programming language  
 R

## Data Based

The paper provided two different types of cell-gene matrix data sets, 10xv2 and Ssv4 (SMART-Seq v4 kit). 
The size of the 10xv2 datasets was too big for us to utilize due to our computer ram storage issue. 
Therefore, we planned to choose the data from method Ssv4.
The Seurat package provided in R is the main package we applied for data analysis.
Within the Ssv4 datasets, the paper included the data stored as rda data, datatypes created by R.
Based on the paper the dataset include the part of Hippocampal (HIP) where further separated into multiple clusters.
Since our data is focusing on the hippocampus, we decided to choose the ss.rda data provided in the Ssv4 method. 

## Files Detail

three Rmarkdown files and one Rshiny file were created: 

1. loadWholeCortex.Rmd 
2. Hip_Subregions.Rmd 
3. HIP_SUB_INFO.Rmd 
4. App.R 

In the loadWholeCortex data we created a HIP.s4obj.cluster function with 0.7 resolution and 16 pcs based on the elbow plots results.
By the default setting in Dimplot, we can then check how the cells cluster in umap.
Clusters need to be annotated since the data we re-cluster will only show numbers.
Therefore, we test out several known marker genes; CA1 ventral: Dcn, CA1 dorsal: Wfs1.
Yet, we were facing a problem while processing the annotation steps.

The paper mentions that HIP (hippocampal) includes many other regions such as retrohippocampal (RHP), lateral, and medial entorhinal cortex para-subiculum (PAR), post-subiculum (POST), Pre-subiculum (PRE), subiculum (SUB), pro subiculum (Pros), and a few remaining small regions. Since we only care about the hippocampus, we create another program called Hip_subregions.Rmd that filters out any other clusters other than CA1, CA2, CA3, DG, and Sub-Pros.  We had successfully re-cluster and annotate the dorsal and ventral part for CA1, CA3, and DG by 16 Pcs and 0.4 resolution.

By observing the cluster information, we find that CA2 merges to CA3 while Sub-Pros merges into CA1. We suspect that is because the small number of cells show expression, only 19 out of 4437 cells show CA2 expression, while 33 cells express in SUB-Pros. By comparing the known marker gene for CA1 dorsal, CA1 ventral, CA3 dorsal and CA3 ventral we were able to rename our ident where clusters 3 = CA1 dorsal, 1 = CA1 ventral, 5 = CA3 dorsal, 6 = CA3 ventral. As for the DG part, 5 clusters show expressions (0, 7, 2 ,4, 8). We merge the cluster 0 and 7 as DG dorsal, and 2 4 8 as DG ventral. To find the marker genes with differential expression requires the function called FindAllMarkers. All the results generated by differential expression have been double-checked by the sagitta in situ Allen brain development brain before we rename the clusters. 

Other than HIP, hippocampal, the original dataset also contains information on the subiculum.  We are curious whether we can also find the dorsal and ventral part of the subiculum if we merge the hippocampal part and subiculum part. Therefore, we created another program: HIP_Sub_Info.Rmd. First, we need to filter out any other clusters other than CA1, CA2, CA3, DG and Sub-Pros like we did in the HIP_Subregions.Rmd. We then merge this dataset with the dataset that only contains hippocampal information. By going through the pipeline with 16 PCs and 0.5 resolution, we observed that the subiculum was able to separate into the dorsal and ventral part as well.  By checking the marker genes selected by differential expression like we did in the previous project and validate with Allen brain development brain. Clusters were successfully annotated as 1 = CA1v, 3 = CA1d, 8 = CA3d, 9 = CA3v, 5 = SuBd, 6 = Subv , 0 and 7 = “DGv”, 2, 4, 10 = DGv. Although the dorsal DG part didn’t come out as well as other clusters. We did find a gene Tdo2 that specifically specifies the image of Tdo2 from Allen Bain where we can prove our annotation.


