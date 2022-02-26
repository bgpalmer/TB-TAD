# Preprocessing the data for consumption

`preprocess.r` takes data and processes it into a standard file that can be given to the machine learning step. The preprocessing is largely followed from the preprocessing script found [here](https://github.com/jchyr-sbmi/PredTAD/blob/master/Codes/1_gen_pre_info_08132020.r).


## R packages 

* `Lubridate`
* `devtools`
* `BiocManager`
* `BiocManager::install("FDb.InfiniumMethylation.hg19")`
* `tidyverse`
* `plyr`
* `BiocManager::install("rCGH")
* `install.packages('data.table')`
* `tictoc`
