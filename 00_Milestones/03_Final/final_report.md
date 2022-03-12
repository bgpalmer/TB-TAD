## Introduction

DNA is a helical three dimensional structure that encodes life through the transcription of deoxynucleotides into messenger RNA (mRNA), which is then used to create new biological structures (ie. translation into proteins). This is essentially the central dogma of biology. Regulatory elements upstream control when DNA should be transcribed, which means that they can control what inputs are provided for downstream physiological processes. Because of this, understanding the causes of upregulation and downregulation of gene transcription and the mechanisms behind them area an important part of understanding biology. 

In the cell, DNA can take on many shapes because of it's lifecycle and environment. During metaphase in the cell cycle, DNA is wound up tightly around simple protein structure called histones, which are then wound up into nucleosomes until the final product is reached: the aligned, diploid chromosomes attached at the centromere. While not replicating, DNA is more of a tangled web within the nucleus (strictly speaking for humans), which allows for it to be opened and transcribed. 

[TODO: need to flow into next paragraph better]

Thanks to on going efforts in next generation sequencing, new assays and techniques have developed recently that fall under the umbrella of 3 Dimensional Chromosome Conformation Capture (3C), which identify DNA coordinates that are near each other. This occurs by ligating (connecting) sequences together while they are within the cell with the goal of sequencing these connected DNA together. One particular 3C technology that has gained traction is Hic data, which allows for the use of next generation sequencing to get high fidelity reads of DNA that are close in space. 

Thanks to these technologies, researchers have learned DNA tends to have a programmed tangled structure that is cell type specific [TODO: source?]. Indeed, there are areas that tend to be more open that allow for transcription to occur and areas that are closed. These areas are known as topologically associating domains (TADs), and have shown to be the basin of biological processes. 

More research has occurred to understand TADs and what characterizes their boundaries. Conformation capture technology shows more contacts been DNA coordinates with a TAD than between coordinates that are in separate TADs, and therefore the boundaries implicitly define what relationship can occur. Thus, the relationships of genes (units of heredity) tend to be contained within tads as demarcated by the boundary [TODO: fix this sentence]. Transcription binding sites and start sites have shown to be more dense around these boundaries as well. Methylation of DNA also tends to be dense around these regions. There are various technologies that exist to capture these data and can provide more information about where TADs exist.

The current issue with 3C technology is that it is expensive. New efforts to identify TADs have been to train machine learning algorithms to predict TADs using epigenetic features that have also shown to be dense around these domains. A paper that recently came out showed success using a supervised boosted gradient machine algorithm to identify tads using various cell lines. The paper made use of preexisting data that is out in public databases, thus making it a viable option to reproduce their steps.

Our goal is to replicate work done in this field using pytorch and to see if we are able to modify the algorithm to increase accuracy. In doing so, we attempt to contribute to TAD prediction. 

