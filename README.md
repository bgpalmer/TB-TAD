# TB-TAD

## Authors

* Tyler Collins
* Brian Palmer

## Project goals and descriptions 

Our goal is to create an unsupervised learning algorithm to detect topologically associating domains (TADs) in chromosomes. TADs are regions with localized highly interacting units of heredity known as genes and have been extensively researched over the last 5 years. Furthermore, sub units in these regions tend to interact more frequently with each other than neighboring regions, thus demarcating boundary lines. Data is generated using chromatin conformation capture (3C) and is stored as matrixes of interaction counts. With this data, the spatial organization of chromosomes can be investigated.

There are a number of published tools already. [ClusterTAD](https://github.com/BDM-Lab/ClusterTAD) has Java and MATLAB implementations to detected TADs through an unsupervised method. [PredTAD](https://github.com/jchyr-sbmi/PredTAD) classifies regions as boundary or non boundary regions by looking at intra-chromosomal and epigenetic data, such as where methylation or protein interactions occur. 

## Member Roles



## Resources



## Reservations

This will be unsupervised model and will be difficult to validate that TAD clusters were identified correctly. We will try to reproduce findings from other published models to determine if we are on the right course. This should also prevent the need to quality control the data, and will give us an idea how to normalize the inputs for biases if not already completed. At a minimum, we hope to identify clusters of TADs, even if we are not as successful as other published models.

## Relationship to your background

Both authors are currently enrolled in the University of Oregon's Bioinformatics and Genomics Masters Program. Hi-C is widely used to detect the 3D dimensional structure of chromatin and new data is expected to continue to be produced. We hope that this algorithm will be useful to researchers studying the conformation of chromatin. 

## To compile for submission

```bash
pandoc README.md -f markdown -t pdf -o tyler_brian_ML_proposal.pdf --metadata title="TB-TAD" --metadata author="Tyler Collins, Brian Palmer"
```