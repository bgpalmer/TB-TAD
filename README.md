# TB-TAD

## Authors

* Tyler Collins
* Brian Palmer

## 1. Project goals and descriptions 

Our goal is to create an unsupervised learning algorithm to detect topologically associating domains (TADs) in chromosomes. TADs are regions with localized highly interacting units of heredity known as genes and have been extensively researched over the last 5 years. Furthermore, sub units in these regions tend to interact more frequently with each other than neighboring regions, thus demarcating boundary lines. Data is generated using chromatin conformation capture (3C) and is stored as matrixes of interaction counts. With this data, the spatial organization of chromosomes can be investigated.

There are a number of published tools already. [ClusterTAD](https://github.com/BDM-Lab/ClusterTAD) has Java and MATLAB implementations to detected TADs through an unsupervised method. [PredTAD](https://github.com/jchyr-sbmi/PredTAD) classifies regions as boundary or non boundary regions by looking at intra-chromosomal and epigenetic data, such as where methylation or protein interactions occur. Using epigenetic markers also helps to expand available inputs (add data from multiple experiments) to train the model, instead of relying entirely on interaction counts found from one single experiment. We think the latter paper will be a good place to start; they use a gradient boosting machine to classify regions as boundary or non boundary areas, which demarcates where TADs exist. In their results, they successfully identify TADs and compete with other models that use Bayesian Additive Regression Trees and Position Specific Linear Models. Additionally, they are able to determine which features are most important for determining TADs in cancer data.

Our goal is to reimplement this approach and compare our results against other published models. Our metrics of success will include finding TAD boundaries and identifying top features to determine boundaries. We will also use Random Forests to classify chromosomal regions as boundary or non-boundary of TADs. We will use data from GEO and ENCODE databases to train our model, which will include contact matrices and epigentic markers. We then hope to use actual data collected by our mentor in a separate project to see how well it predicts TADs. 

## 2. Member Roles

Both members of the team are familiar with the conformation capture technologies and it's output. We will together focus on how we want to shape our model so that we can compete with other published models. 

## 3. Resources

We will rely on GEO and ENCODE databases for data to train and test our model. 

Below enumerates where we will access our training data.

### 3.1 Hic Contact Maps

Here are two counts matrices from two cell lines.  

```bash

# Breast Cancer Hi-C, Hg19, PredTAD

# Readmes
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_OVERALL_README.rtf
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary_README.rtf

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary_intrachromosomal_contact_matrices.tar.gz

# Contains normal breast epithelial cell line data and cancerous counterparts
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66733/suppl/GSE66733_Hi-C_MCF7_MCF10A_processed_HiCfiles.tar.gz

```

An additional matrix is found here that is mapped to a more recent version of the human genome (hg38). This file shows contacts for a colon cancer cell.

```bash

# Colon Cancer Hi-C, Hg38
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE172097&format=file

```

### 3.2 Methylation data

Methylation data for the above cell lines can be found on [ENCODE](https://www.encodeproject.org/). We will download one of the data files when we have a sense which format is correct.

### 3.3 Annotation

Since this work is inspired by PredTAD and we are using similar, we will look to use their annotation scripts found in their [repository] (https://github.com/jchyr-sbmi/PredTAD). They use an r script to create input files for their models to incorporate the data we have listed in sections 3.1 and 3.2.

### 3.3 Additional Resources

If time and availability permits, we will also try to test our model on data created by our mentor in a separate project. Using Random Forests could be memory intensive unless we use Gradient Boosting. In the event that we need more memory, we have access to Talapas, a super computer hosted by the University of Oregon, that permits tasks to allocate memory upwards of approx. 150 GB. 

## 4. Reservations

This will be unsupervised model and will be difficult to validate that TAD clusters were identified correctly. We will try to reproduce findings from other published models to determine if we are on the right course. This should also prevent the need to quality control the data, and will give us an idea how to normalize the inputs for biases if not already completed. At a minimum, we hope to identify clusters of TADs, even if we are not as successful as other published models. If we need to reduce the scope of the project, we may decide to focus on gene interaction maps and forgo epigenetic inputs. This could also make formatting the input less complex, which will be a new challenge for both authors. 

## 5. Relationship to your background

Both authors are currently enrolled in the University of Oregon's Bioinformatics and Genomics Masters Program. Hi-C is widely used to detect the 3D dimensional structure of chromatin and new data is expected to continue to be produced. We hope that this algorithm will be useful to researchers studying the conformation of chromatin. 

## Internal: To compile for submission

```bash
pandoc README.md -f markdown -t pdf -o tyler_brian_ML_proposal.pdf --metadata title="TB-TAD" --metadata author="Tyler Collins, Brian Palmer"
```