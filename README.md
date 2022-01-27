# TB-TAD

## Collaborators

* [Tyler Collins](https://github.com/tcollins12)
* [Brian Palmer](https://github.com/bgpalmer)

## Project goals and descriptions 

Our goal is to create an unsupervised learning algorithm to detect topologically associating domains (TADs) in chromosomes. TADs are regions with localized highly interacting units of heredity known as genes and have been extensively researched over the last 5 years. Furthermore, sub units in these regions tend to interact more frequently with each other than neighboring regions, thus demarcating boundary lines. Data is generated using chromatin conformation capture (3C) and is stored as matrixes of interaction counts. With this data, the spatial organization of chromosomes can be investigated.

There are a number of published tools already. [ClusterTAD](https://github.com/BDM-Lab/ClusterTAD) has Java and MATLAB implementations to detected TADs through an unsupervised method. [PredTAD](https://github.com/jchyr-sbmi/PredTAD) classifies regions as boundary or non boundary regions by looking at intra-chromosomal and epigenetic data, such as where methylation or protein interactions occur. Using epigenetic markers also helps to expand available inputs (add data from multiple experiments) to train the model, instead of relying entirely on interaction counts found from one single experiment. We think the latter paper will be a good place to start; they use a gradient boosting machine to classify regions as boundary or non boundary areas, which demarcates where TADs exist. In their results, they successfully identify TADs and compete with other models that use Bayesian Additive Regression Trees and Position Specific Linear Models. Additionally, they are able to determine which features are most important for determining TADs in cancer data.

Our goal is to reimplement this approach and compare our results against other published models. Our metrics of success will include finding TAD boundaries and identifying top features to determine boundaries. We will also use Random Forests to classify chromosomal regions as boundary or non-boundary of TADs. We will use data from GEO and ENCODE databases to train our model, which will include contact matrices and epigentic markers. We then hope to use actual data collected by our mentor in a separate project to see how well it predicts TADs. 