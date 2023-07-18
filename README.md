# Epigenetic_Project

In this study, we have analysed epigenomic mark H3K27Ac mark distribution pre epigenetic drug treatment in two different cell lines i.e. HCT116 and RH4.
Based on change in gene expression after epigenetic drug treatment, we identified 1000 up and 1000 downregulated genes. For these genes, we analyse the pre-treatment H3K27Ac marks across 21 genomic bins. These bins were created using Promoter (10 bins), TSS and Gene Body (10 bins).
Next, we developed machine learning model for predicting given the pre-treatment transcriptome and epigenomic profile of a sample, the extent to which one can predict locus-specific changes in gene expression upon treatment with HDACi.

We observed that in two cell lines (HCT116 treated with Largazole at 8 doses and RH4 treated with Entinostat at 1µM) where the appropriate data (pre-treatment transcriptome and epigenome as well as post-treatment transcriptome) is publicly available, our model could distinguish the post-treatment up versus downregulated genes with high accuracy (up to ROC of 0.89). Furthermore, a model trained on one cell line is applicable to another cell line suggesting that such a model can be applied to a novel context. 

Here we present a first assessment of the predictability of genome-wide transcriptomic changes upon treatment with HDACi. Lack of appropriate omics data from clinical trials of epigenetic drugs currently hampers the assessment of applicability of our approach to predicting the response to epigenetic drugs in cancer patients.
