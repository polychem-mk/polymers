## Polymers dataset. Analysis and classification

In this project, we explore the **polymers** dataset from *kaggle.com*
and build a classification model to predict the type of polymers. The data contains Morgan fingerprints for plastics, peptides, and oligosaccharides. 

### Overview
 
Polymers are organic macromolecules made up of repeating  units - monomers.
In chemistry, polymers can be classified in several ways: depending on certain properties, the method of synthesis, the nature of the monomers, etc.
There are a large number of known monomers, and combining different monomers can give us a potentially infinite number of  unique  polymers.   
<img src="images/image.png" align="right" width="25%" height="25%"/>  
Morgan fingerprint represents the chemical structure of a molecule as a vector with values of *0* or *1*. These values indicate the pretense (*1*) or absence (*0*) of a predefined chemical substructure within a circular radius.  And this allows any molecules to be encoded using a finite number of simple predictors.  

The goal of this project is to build a classification model based on Morgan fingerprint data, and it consists of two main parts.  
• The first part is to find out what types of polymers are present in this dataset.  
Although there are labels for three types of polymers - plastics, peptides and 
oligosaccharides - these are general names for polymers and can include very different groups.   

For this purpose, we use the rcdk package to look at the chemical structures of 
the polymers represented in this dataset, and apply unsupervised learning algorithms 
to find polymer types that are not obvious from the chemical structures,
and add polymer subtypes to the existing labels.

• The second part is to find the best model that can be used to classify polymer types.
For this purpose, we will use supervised machine learning methods, where the target 
variable will contain the polymer types that were identified in the first part. 
In this project, we compare three models: K-nearest neighbors, classification tree, and random forest.   

### Data

The dataset is available at  
<https://www.kaggle.com/datasets/victorsabanzagil/polymers/data>   

In this project, the original *polymers_dataset.csv* file was saved as
*polymers_data.rds* and added to the repository.

### Packages
<pre>
R version 4.4.3 (2025-02-28)  

parallelDist_0.2.6    ranger_0.17.0     rpart_4.1.24      caret_7.0-1 
lattice_0.22-6        stringr_1.5.1     rcdk_3.8.1        rcdklibs_2.9
rJava_1.0-11          gt_0.11.1         patchwork_1.3.0   dendextend_1.19.0
ggplot2_3.5.1         tidyr_1.3.1       dplyr_1.1.4.9000  here_1.0.1 
rmarkdown_2.29        knitr_1.49 
</pre>
  
### Files  
<pre>
├── LICENSE 
├── polymers.Rproj 
├── README.md 
├── index.html                             # HTML report
├── data 
│   ├── processed 
│   │   └── polymers_data_processed.rds    # file with processed data
│   └── raw 
│       └── polymers_data.rds              # file with data
├── images 
│   └── image.png                          # Morgan fingerprint image
└── scripts 
    ├── bibliography.bibtex                # bibliography for HTML report
    ├── citation_style.csl                 # citation styles for HTML report
    ├── polymers_analysis.R                # R file with data analysis
    ├── polymers_report.Rmd                # Rmd file to render HTML report
    ├── styles.css                         # styles for HTML report
    ├── functions 
    │   ├── print_gt_table.R 
    │   └── print_molecules.R 
    └── models 
        ├── hclust_model.rds 
        ├── knn_model.rds 
        ├── knn_predicted.rds 
        ├── pca_model_sdev.rds 
        ├── PCs_10.rds 
        ├── ranger_final.rds 
        ├── ranger_model.rds 
        ├── rpart_model.rds 
        └── tune_mtry.rds 
 
</pre>
  
  
The *citation_style.csl* file was downloaded from this repository:
[github.com/citation-style-language](<https://github.com/citation-style-language/styles/blob/master/american-chemical-society.csl>)  
and follows the American Chemical Society citation style.  
  
  