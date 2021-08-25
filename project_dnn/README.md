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

* Introductory material on ANNs and Deep Learning
([Deep Learning TDS](https://towardsdatascience.com/simply-deep-learning-an-effortless-introduction-45591a1c4abb))
([Deep Learning OREILLY](https://www.oreilly.com/library/view/deep-learning/9781491924570/ch01.html))
([ANNs TDS](https://towardsdatascience.com/basic-concepts-of-neural-networks-1a18a7aa2bd2))
([ANN Sci2lab](https://sci2lab.github.io/ml_tutorial/neural_network/))

* Short introductory video lectures by Andrew Ng (highly recommended). 
This is a full intro course to Machine Learning.
For the purposes of this course focus on lectures 8.1-8.7, 9.1-9.7, 17.1-17.4
([Andrew Ng Lectures](https://www.youtube.com/playlist?list=PLLssT5z_DsK-h9vYZkQkYNWcItqhlRJLN))



### Tentative schedule

For the complete course schedule please refer to:

([Course Schedule](https://nbisweden.github.io/single-cell_sib_scilifelab_2021/schedule.html))


<br/>

### Milestone 1:

1. Walkthrough examples 
In this part of the course the student will follow prepared examples in order to familiarize themselves with the prinicples of model development and the Keras environment, and basic concepts for model development for single cell RNAseq data.

1.1. Creation of a simple MLP and vanilla Autoencoder model in Keras

1.2. Building a Variational Autoencoder model for single cell data

1.2.1 Data preprocessing and transformations

1.2.2 Model specification and arhitectural choices

1.2.3 Loss functions

1.2.4 Training callbacks and monitoring

1.2.5 Model evaluation and diagnostics

1.2.6 Performing inference with a trained model

<br/>

### Milestone 2:

2. Student project

In this section students will be split into two groups and will be assigned analysis tasks on a new dataset using generative models for single cell RNAseq similar to the ones used in the walkthrough examples

