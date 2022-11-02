# ME-Class3
ME-Class trains a model to use differential DNA methylation to predict expression changes in one or more pairs of samples. ME-Class3 is the latest updated version of ME-Class. This version can integrate methylation information from gene promoter and enhancer regions. It can also integrate data from both 5mC and 5hmC.  

## Requirements
Tested on python 3.10.6. Will definitely not work with python < 3.6.
 * Python Packages
	 *  dataclasses
	 *  matplotlib
     *  numpy
     *  pandas
     *  scipy
     *  scikit-learn
	 *  seaborn

## Installation
    Clone and install with pip:
      git clone https://github.com/edwardslab-wustl/ME-Class3.git
      cd ME-Class3
      pip install .

## Docker information
A dockerized version of ME-Class3 is available at: 

## Usage
    meclass3 <subcommand> <arguments>
    Use meclass3 -h to see full list of subcommands.
    Use meclass3 <subcommand> -h to see help for individual subcommands.

### Basic workflow
  1. **preprocess** - Takes input methylation data for each sample and computes differential methylation for each specified pair.
  1. **overlap_genes_enhancers** (opt) - overlaps gene annotations and a bed file of enhancer locations to link genes with specific enhancers
  1. **interp** - Use preprocessed methylation data with gene or enhancer annotation data to compute methylation signatures
  1. **plot_interp** (opt) - plot methylation signatures
  1. **merge_data** - merge methylation signature from interp step along with expression data to make input data files for the classifier
  1. **classify** - using merged data as input, train model to predict expression labels from methylation signatures
  1. **plot_performance** - use output from one or more classify runs to make ROC and Acc vs Reject Rate plots
  1. **pull_gene_list** - pull list of accurately predicted genes from classify step
  1. **cluster** (opt) - unsupervised clustering of classify output to determine methylation signatures that yield accurate predictions

### Input files

### Example

## License information
ME-Class3 is distributed under the GPL-3.0 Licence. 


## Parameters

 
