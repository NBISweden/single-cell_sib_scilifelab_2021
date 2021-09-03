# <img border="0" src="../logos/deep_learning.png" width="40" height="40" style="vertical-align:middle;"> Deep Learning
***

<br/>

- __*Panagiotis Papasaikas*__, Friedrich Miescher Institute, SIB/Basel, 🇨🇭 Switzerland
- __*Geert van Geest*__, Interfaculty Bioinformatics Unit, SIB/UniBern, 🇨🇭 Switzerland

<br/>

### Pre-course


### Background

Deep learning approaches have found their way on various aspects of the analysis of single cell transcriptomic data. These include deep generative models for dimensionality reduction and visualization,  imputation, out-of-sample inference, batch correction, integration and  exploration of data across different modalities. 

In this section of the course we aim to provide an overview of the field and the different out-of-the box solutions and to take a closer look at the theoretical concepts behind those tools. 

In addition we will provide practical examples for the development, training and evaluation of custom architectures using  ML libraries and apply them to the analysis of provided scRNA-seq datasets. 


#### Software requirements 

We will be using an AWS cloud environment, where we have pre-installed the required software and downloaded the data.

You will be able to access an R-studio and/or a Jupyter notebook environment in order to follow the course examples and to work on the assigned projects. 

The AWS environment will also provide access to GPUs required for the swifter training of the models that will be used in the course.


#### Knowledge requirements

For the practical sections of the course we will be using the TensorFlow and Keras high-level API for neural network development, through the R interface.
In case you are new to the use of Keras you can familiarize yourself with the environment using the gentle tutorials found here:

* Getting started in Deep Learning with keras  ([tutorial](https://machinelearningmastery.com/tensorflow-tutorial-deep-learning-with-tf-keras/))

* Basic ML with RKeras ([tutorial](https://tensorflow.rstudio.com/tutorials/beginners/basic-ml/))


Although we will make an effort to cover material required to follow the course, we assume that students are familiar with basic concepts of Machine Learning and Artificial Neural Networks such as:

To review some of these concepts you can refer to some of this introductory material:

* Introductory material on ANNs and Deep Learning:

([Deep Learning TDS](https://towardsdatascience.com/simply-deep-learning-an-effortless-introduction-45591a1c4abb))

([Deep Learning OREILLY](https://www.oreilly.com/library/view/deep-learning/9781491924570/ch01.html))

([ANNs TDS](https://towardsdatascience.com/basic-concepts-of-neural-networks-1a18a7aa2bd2))

([ANN Sci2lab](https://sci2lab.github.io/ml_tutorial/neural_network/))

* Short introductory video lectures (~5'-10' each) by Andrew Ng and Stanford University (highly recommended). 
This is a full intro course to Machine Learning.
For the purposes of this course focus on lectures 8.1-8.7, 9.1-9.7, 10.3-10.7, 17.1-17.4

([Andrew Ng Lectures](https://www.youtube.com/playlist?list=PLLssT5z_DsK-h9vYZkQkYNWcItqhlRJLN))



### Schedule

For the complete course schedule please refer to:

([Course Schedule](https://nbisweden.github.io/single-cell_sib_scilifelab_2021/schedule.html))


<br/>


### 1. Walkthrough examples 


In this part of the course the student will follow prepared examples in order to familiarize themselves with the prinicples of model development and the Keras environment, and basic concepts for model development for single cell RNAseq data.

#### 1.1. Creation of a simple MLP and vanilla Autoencoder model in Keras


#### 1.2. Building a Variational Autoencoder model for single cell data

* 1.2.1 Data pre-processing and transformations

* 1.2.2 Model specification and architectural choices

* 1.2.3 Loss functions

* 1.2.4 Training callbacks and monitoring

* 1.2.5 Model evaluation and diagnostics

* 1.2.6 Performing inference with a trained model

<br/>



### 2. Group project:

In this section students will be split into two groups of four people and will be assigned analysis tasks on a new dataset using a generative VAE model for single cell RNAseq similar to the one used in the walkthrough examples. 

A guiding script file will also be provided in order to usher you along the different project tasks. 
However you should feel free to step outside the script in case you want to experiment or in case you see a more ingenuous way to carry-out a task!


For the group projects we will use data from the [TabulaMuris](https://tabula-muris.ds.czbiohub.org/) compendium of single-cell mouse transcriptome data. The provided dataset is pre-filtered in terms of genes and cells so no additional pre-processing is required.

It contains data from two **different technologies** (attribute *study*: droplet/facs):

- *microfluidic droplet-based* 3’-end counting: provides a survey of thousands of cells per organ at relatively low coverage and

- *FACS-based* full length transcript analysis: provides higher sensitivity and coverage but at the expense of lower cell counts.

The provided cells are also annotated among others in terms of **animal sex** (attribute *mouse_sex*: F/M) and **tissue of origin** (attribute *cell.class*). There are 5 Tissues of origin:

*Bladder, Limb_Muscle, Lung, Mammary_Gland and Thymus.*


#### 2.1 
The participants will be asked to build one single **VAE model** for the entire dataset.
First:

Group 1 remove all *female bladder cells from the facs data* 

Group 2 remove all *mammary-gland data also from the facs data*. 

These will be used as a holdout test-set for out-of-sample inference using your model later on. 
Finally randomly split 80/20 your remaining data in training  (80%) and validation sets (20%).



#### 2.2
Experiment briefly with **different architectural decisions and hyperparameters**:
Number of layers, layer width, dropout values, beta-parameter for the loss. 
How do these decision impact convergence time, accuracy of the imputed values, number of imputed drop-outs?

In this step we  strongly suggest sub-sampling your dataset to at most 10000 cells in order to speed-up training.


#### 2.3 
Once you have settled on one model architecture, **train your model** using your dataset  using the 80-20 split between training and validation until full convergence. Show your training history. 
From this point on make certain your whole group works using the same fully trained model and corresponding set of weights.


#### 2.4 Evaluate your model.

* 2.4.1 Compare library sizes and mean gene expression values in the input and the output (denoised) data.

* 2.4.2 Show the mean variance trend in the input and outpu (denoised) data.



#### 2.5 Using and interpreting the latent space

* 2.5.1 Use you model's latent space to produce a tSNE visualization of your data. How does it compare to the tSNE visualization produced using as input PCs?

* 2.5.2 Color you projection according to the cell values of different latent nodes.

* 2.5.3 Which are the top genes affected in the denoised output when perturbing the values of individual latent nodes? Do they make sense in terms of which cell-types have higher loadings for those particular latent nodes?

 
#### 2.6 Inference task 1: Dropout-imputation and outlier correction (denoising).


* 2.6.1 Group 1:
Corrupt a random subset of your input data by droupouts (zeroes). Impute the corrupted values using your model. How do they compare between FACS and Droplet data? How do they compare when you repeat the same task in the validation or hold-out datasets? Produce scatterplots to compare initial observed values vs corrupt values and initial observed values vs denoised values.

* 2.6.2 Group 2:
Corrupt a random subset of your input data using large-magnitude outliers. Select genes with non-zero observed counts and force a 2-4x fold change in the observed counts in either direction. How do they compare between FACS and Droplet data? How do they compare when you repeat the same task in the validation or hold-out datasets? Produce scatterplots to compare initial observed values vs corrupt values and initial observed values vs denoised values.


####  2.7 Inference task 2: Batch correction using latent arithmetic.
In this task you will try to correct for the effects of batch resulting from the two different technologies.

* 2.7.1 Group 1: Move the facs data to droplet space by adding to each of their latent representations a technology-delta latent vector. Produce a tSNE representation for the complete dataset using these corrected latent representations. Color cells by technology. Compare to the uncorrected (latent-based) tSNE representation. 

 
* 2.7.1 Group 2:Move the droplet data to facs space by adding to each of their latent representations a technology-delta latent vector. Produce a tSNE representation for the complete dataset using these corrected latent representations. Color cells by technology. Compare to the uncorrected (latent-based) tSNE representation. 

#### 2.8 Inference task 3: Out-of-sample prediction using latent arithmetic.
In this final task you will be asked to use your model to infer location characteristics for the distributions of cells held-out during model construction. Note that there are a few way to go about each prediction task.


* 2.8.1 Group 1: Predict the mean gene expression of *female bladder cells of facs data* by:

- Starting from male bladder facs cells and adding a gender-delta latent vector, also calculated from facs data.
- Starting from female bladder cells of droplet data and adding a technology-delta latent vector. 


* 2.8.1 Group 2: Predict the mean gene expression of *mammary-gland of facs data* by:

- Starting from droplet mammary-gland cells and adding a technology-delta latent vector. 

- Starting from (non-mammary-gland) facs cells and adding mammary-gland delta latent vector calculated using the droplet data. 

*VERY IMPORTANT* for all inference tasks in 2.7 and 2.8: Remember for all calculations of mean latent vectors to use *balanced* samples in terms of their composition for uninteresting (known) sources of variance. For example when calculating a technology-delta latent vector, the two sets of cells used for the calculation of the average facs and average droplet latent vectors should be balanced in terms of tissue AND mouse gender.   




#### 2.9 Optional bonus task: Fine tuning.

Fine-tuning refers to the strategy of starting from a model trained on a large corpus of data and then adjusting the model weights for a new typically smaller (out-of-distribution) dataset. 
This strategy allows e.g for model-training when insufficient data do not allow for full *ab-initio* model training in the smaller dataset.
A simple way to achieve this is to use a very small learning rate for weight adjustment in the fine-tuning step. 

* Starting from your fully trained VAE fine-tune a model specific for the (held-out) mammary-gland facs data. Perform basic evaluation tasks and a latent based tSNE projection (as in sections 2.4  and 2.5.1) to appraise this facs mammary-gland specific model.  


