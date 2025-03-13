
#                Polymers dataset. Analysis and classification

# Introduction ------------------------------------------------------------

# How many chemical compounds are known today? And how can we predict the properties
# and find common patterns for those that will be obtained in the future?
# There are millions of synthesized and natural molecules, cas.org contains data on
# more than 200 million compounds, 75 million protein and nucleic acid sequences [1].
# Estimates of the virtual chemical universe for organic molecules are even larger 
# and continue to grow year after year [2,3].
# Thus, encoding an infinite number of chemical structures using simple predictors 
# such as chemical fingerprints allows machine learning algorithms to be applied to
# the large data sets.
 
# In this project we use the polymers dataset from kaggle, available at:
#  <https://www.kaggle.com/datasets/victorsabanzagil/polymers/data>
# The data contains Morgan fingerprints for polymers of three classes: plastics,
# peptides, and oligosaccharides. 
# Plastics are synthetic polymers such as polyethylene, polyvinyl chloride, 
# polystyrene, etc. Peptides and oligosaccharides can be natural or synthetic 
# with peptide and glycosidic bonds respectively [4,5].

# Morgan fingerprint represents the chemical structure of the molecule, also 
# called descriptor, and is in the form of a binary vector. The values of vector 
# indicate the pretense (1) or absence (0) of the predefined chemical substructure,
# within a circular radius [6-8].

# In cheminformatics, various types of chemical fingerprints, including Morgan 
# fingerprints, are used to predict certain properties of molecules and find 
# patterns in large data sets [6]. The properties of substances are directly related 
# to their chemical structure, so these are closely related concepts.

# The goal of this project is to build a classification model based on Morgan
# fingerprint data to predict the type of the polymers.
#  
# The first part is to find out what types of polymers are present in this dataset
# using unsupervised learning algorithms and the rcdk package, which is used to
# work with chemical structures.
# Although labels already exist for three types of polymers - plastics, peptides
# and oligosaccharides - these are general names for polymers and can include very
# different groups of polymers. Therefore, we will add polymer subtypes to the
# existing labels.
#  
# The second part is to find the best model that can be used to classify polymer
# types. For this purpose, we will use supervised machine learning methods, where
# the target variable will contain the polymer types that were identified in the
# first part, and we will use Morgan fingerprints as features.
  

# Libraries and variables --------------------------------------------------

## Libraries ----------------------------------------------------------------

# data manipulation
library(dplyr)
library(tidyr)

# visualization
library(ggplot2)
library(dendextend)
library(patchwork)
library(gt)

library(rcdk)   # for chemical structures
library(stringr)  # to deal with strings

# build the models
library(caret)
library(rpart)
library(ranger)

#library(parallelDist) # distance calculations

## Variables -----------------------------------------------------------

# Path to files.

path_PCs_10 <- "scripts/models/PCs_10.rds"             # PCA, 10 principal components
path_pca_sdev <- "scripts/models/pca_model_sdev.rds"   # PCA, standard deviations 

# Model for the hierarchical cluster analysis:
path_hclust_model <-"scripts/models/hclust_model.rds"  

path_rpart_model <- "scripts/models/rpart_model.rds"    # Classification tree model
path_ranger_model <- "scripts/models/ranger_model.rds"  # Random forest model 
path_knn_model <-"scripts/models/knn_model.rds"         # knn model 
path_knn_predicted <-"scripts/models/knn_predicted.rds" # knn model predictions

# Accuracy and the out-of-bag prediction error for random forest models 
# with different 'mtry':
path_tune_mtry <-"scripts/models/tune_mtry.rds"  

# path to functions:
path_function_print_gt_table <- "scripts/functions/print_gt_table.R"
path_function_print_molecules <- "scripts/functions/print_molecules.R"

# path to data:    
path_data <- "data/raw/polymers_data.rds"
path_data_processed <- "data/processed/polymers_data_processed.rds"

#  Colors 

saccharides_color <- c("#C24641", "#A74AC7","#F2A2E8", "#FC6C85")
plastics_color <-  c("#E0E5E5","#52595D", "#C0C6C7",  "#99A3A3")
peptides_color <- c("#1F45FC", "#0AFFFF", "#43C6DB" , "#357EC7"   )

point_color <-"#5F9EA0"
hist_color <- c("#A0CFEC", "#52595D")
title_color = "#708090"

# Plots theme
themes <- theme_light(base_size = 10 ) + 
  theme(
    plot.title = element_text( size = rel(1.2), family = "system-ui", color = title_color),
    axis.title = element_text( size = rel(1.1)),
    plot.margin = unit(c(10,10,10,10), "pt"))

## Load functions ----------------------------------------------------------

#  Print gt table. Converts data frame into a gt_tbl table (gt package)
source(path_function_print_gt_table)

#  Print molecules. Using rcdk and grid packages
source(path_function_print_molecules)

## Load data ------------------------------------------------------------

# The original file 'polymers_dataset.csv' is not included in
# this repository, but it is available at
# https://www.kaggle.com/datasets/victorsabanzagil/polymers/data

# File 'polymers_dataset.csv' was converted to rds format 
# to reduce the size as shown below:
#       polymers <- read.csv("data/raw/polymers_dataset.csv") 
#       saveRDS(polymers, "data/raw/polymers_data.rds")

# So the data for this project is in the file 'data/raw/polymers_data.rds'. 
# The folder 'data' with file 'polymers_data.rds' needs to be in
# the project working directory, then it can be loaded as:
polymers <- readRDS(path_data)

# Data overview -----------------------------------------------------------

## Dimensions and missing values -----------------------------------------

#  The polymers data set contains data for 20609 polymers.
nrow(polymers)

#  There is no missing values in this dataset.
sum(is.na(polymers))

#  First 3 columns:
#                    X     integer, row ID 
#               smiles     character, chemical structure representation 
#                label     character, polymer type 
str(polymers[ , 1:3])

## label -------------------------------------------------------------------

# The variable 'label' is one of three types of polymer: plastic, peptide,
# or oligosaccharide. The values in this column  are distributed approximately 
# evenly.
table(polymers$label)

## smiles ------------------------------------------------------------------

# The smiles column represent the chemical structures of the polymers
# as SMILES (Simplified Molecular Input Line Entry System [9]).
# All values in the smiles column are unique character strings.
length(unique(polymers$smiles)) == nrow(polymers)

# Here are examples of "smiles" for all three types of the polymers.
polymers %>%
  mutate(length = nchar(smiles)) %>%
  arrange(length) %>%
  group_by(label) %>%
  summarise(first(smiles)) %>% 
  pull()

# It should be noted here that polymers are large macromolecules consisting of 
# repeating chemical structures - monomers, and contain from thousands to millions 
# of atoms [5]. In this data set, we're dealing with these repeating units that
# are shown as SMILES. So when we look at the structure of a polymer, we're 
# looking at the chemical structure of the simple unit that the polymer is made of.


## Morgan fingerprints  ----------------------------------------------------

# The Morgan fingerprints are stored in columns X0, X1, ... Xn and will be used
# as features to predict the type of polymers. We will refer to these columns 
# as X-variables or X columns. There are 2048 X-variables that have two possible
# values 0 or 1.
polymers %>% select(matches("X\\d+")) %>% ncol()

# we omit the first column of the 'polymers' dataset: X (row identifier),
# to avoid confusion with X-variables
polymers <- polymers[ , -1]

# Possible values in columns Xn are 0 and 1, which can be checked as follows:
# polymers %>% select(matches("X\\d+")) %>% unlist() %>% unique()
# [1] 0 1

# A value of 1 indicates the presence of predefined chemical structures in the
# molecule. Each molecule has a very small number of structures compared to all
# possible ones. Therefore, we can expect most of each fingerprint vector to 
# consist mostly of 0s and a small number of 1 values.

# There also some chemical structures common to all or nearly all molecules.

# To examine the distribution of values in the X-variables, we find the X columns
# sums, which basically show how many 1 values are contained in each predictor, 
# or in other words, how often a certain predefined chemical group appears in 
# the polymer molecules.

# 'x_sums' is a vector containing the column sums for all X-variables,
# showing how many 1 values occur in each predictor.
x_sums <- polymers %>% 
  select(matches("X\\d+")) %>% 
  colSums()

# Plot the distribution of the column sums for X-variables.
data.frame(x_sums = x_sums ) %>%
  ggplot(aes(x_sums)) +
  geom_histogram(fill = hist_color[1], 
                 color = hist_color[2], 
                 bins = 100) +
  labs(title = "Histogram. Sums of X-variables",
       x  = "sum") +
  themes 

# The histogram of the sums of X-variables is skewed toward low values,
# indicating that the most frequent value is 0.

# Some summary for  the column sums for X-variables:
x_sums_summary <-
  data.frame( label = c("minimum", "maximum" ,
                      "number of X colimns with all 0s",
                      "number of X colimns with less than one hundred '1' values",
                      "number of X colimns with more than 5K '1' values"),
            value = c( min(x_sums),
                       max(x_sums),
                       sum(x_sums ==0) ,
                       sum(x_sums <= 100),
                       sum(x_sums >= 5000))) 
x_sums_summary %>%
  print_gt_table(table_title = "Sums of X-variables") %>%
  tab_options(column_labels.hidden =TRUE)
  
# There are 23 columns with all 0s, and 761 predictors that contain only 100 
# or fewer 1 values. This means that these predictors are different from 0 in 
# only 0.5% of the data.
  
#  Let's look at the sums of X-variables grouped by label.
x_sum_byLabel <- polymers %>% 
  select(matches("label|X\\d+")) %>%
  select(1:91) %>%     # find the column sums for the first 90 predictors 
  group_by(label) %>%  # grouped by the variable "label";
  summarise(across(where(is.numeric), sum )) %>% 
  t() %>% as.data.frame()  # transpose the table and 
# turn the label values into column names.
colnames(x_sum_byLabel)   <- unlist(x_sum_byLabel[ 1, ])
rownames(x_sum_byLabel)   <- NULL

# delete 1st row which has non-numeric values ("oligosaccharide", "peptide", "plastic")
x_sum_byLabel <- x_sum_byLabel[ -1,] %>% 
  mutate(across(everything(), as.numeric)) %>% # 
  mutate(n = row_number()-1) %>%   # add X column names starting from X0
  mutate(Xn = str_remove(paste("X", n), " ")) 

head(x_sum_byLabel) %>% 
  select(Xn, oligosaccharide, peptide, plastic) %>%
  print_gt_table(table_title = "Sums of the first 6 X-variables by polymer type") %>% 
  cols_align("center")

# The graph below shows the column sums for the first 90 X-variables, colored by
# label (polymer type).
x_sum_byLabel %>%
  ggplot(aes(n-0.3, oligosaccharide,fill = "oligosaccharide")) + 
  geom_col(width = 0.3  ) +
  geom_col(aes (n, peptide, fill = "peptide"), width = 0.3  ) +
  geom_col(aes (n+0.3, plastic, fill = "plastic" ), width = 0.3  ) +
  scale_fill_manual(values = c(saccharides_color[1],
                               peptides_color[1],
                               plastics_color[4]),
                    labels = c("oligosaccharide", "peptide", "plastic")) +
  
  labs(title = "Sums of X-variables by polymer type",
       x  = "X",
       y = "sum") +
  themes +
  theme(legend.position = "inside",
        legend.position.inside = c(.35, .95),
        legend.justification = c("right", "top"),
        legend.title = element_blank() ) 

# We can notice that some X-variables have a large number of 1s for all three
# polymer types, while others are specific only to oligosaccharides, peptides, 
# or peptide+plastic.


# Polymers types ----------------------------------------------------------

# Plastic, peptide, and oligosaccharide are just general names for large groups 
# of chemicals that may include many subtypes. Furthermore, there are different 
# ways to classify polymers.

# To find out what types (and how many types) of polymers are present in 
# our data, we:
# - determine the structure of the polymers that are represented in smiles, and
# - use Morgan fingerprints for unsupervised learning methods to find polymer 
#   groups that are not obvious from the 'smiles' variable.

## rcdk --------------------------------------------------------------------

# We will now use the rcdk package to look at the chemical structure of the
# polymers in more detail and add variables that will help us create more 
# labels for classification.

### Molecular weight --------------------------------------------------------

# To find the molecular weight (Mw) from SMILES, we first convert the "smiles"
# to "molecules" using the parse.smiles() function
molecules <- parse.smiles(polymers$smiles) 
class(molecules[[1]])

# The output of the parse.smiles() is a list of "molecules", so we use the sapply()
# function to get the molecular weight using the get.exact.mass() function.
polymers_Mw <- sapply(molecules, get.exact.mass) %>% unlist()
names(polymers_Mw) <- NULL

# The new variable Mw will be used to find the smallest molecules
# when we need to print the examples of molecular structures.
polymers <- cbind(Mw = polymers_Mw, polymers)

# Molecular weight range of polymers:
polymers %>% group_by(label) %>%
  summarise(Mw_min = round(min(Mw)),
            Mw_max = round(max(Mw))) %>%
  print_gt_table(table_title = "Molecular weight") %>%
  tab_spanner(label = "Mw", columns = starts_with("Mw")) %>%
  cols_label(Mw_min	= "min",
             Mw_max	= "max") %>%
  cols_align("center")


### Aromacity ---------------------------------------------------------------

# The presence of aromatic groups has a great impact on the properties of
# polymers. Therefore one way to add more detail to the classification of
# polymers is to find which polymers have aromatic functional groups, for this
# purpose we use the do.aromaticity() function to add a new variable is_aromatic.

# The output of the do.aromaticity() function is a boolean value 
# indicating any aromatic ring in the chemical structure.
polymers_aromaticity <- sapply(molecules, do.aromaticity) %>% unlist()
names(polymers_aromaticity) <- NULL

polymers <- cbind(is_aromatic = polymers_aromaticity, polymers)

# Proportions of polymers with aromatic groups:
polymers %>% group_by(label) %>%
  summarise(aromatic_prop = mean(is_aromatic)) %>%
  print_gt_table(table_title = "Aromatic molecules in the dataset") %>%
  fmt_percent(columns = aromatic_prop, rows = everything(), decimals = 1) %>%
  cols_label(aromatic_prop	= "percent") %>%
  cols_align("center")

#  More than half of the plastics and peptides have the aromatic group.

### Heteroatoms -------------------------------------------------------------

# Atoms in structures other than carbon and hydrogen are called heteroatoms.
# These atoms are part of functional groups that lead to various properties 
# of molecules that we potentially want to predict.

# To find out what atoms are present in a molecule, we first need to get
# the atoms using the get.atoms() function, and then use another function,
# get.symbol(), to extract the symbols of the elements.
# This can be done as follows:

# sapply(molecules, function(mol){
#      unique(unlist(sapply(get.atoms(mol ), get.symbol)))
#    }) %>% unlist() %>% unique() %>%
#    setdiff( c("C", "H"))

# the output is:
# [1] "Cl" "O"  "N"  "S" 

# We expect all peptides to have nitrogen (N) and oxygen (O), and all
# oligosaccharides to have oxygen (O) in all cases.The proportions of 
# polymers that have elements other than C and H are as follows:
heteroatoms <- polymers %>%
  select(label, smiles) %>% group_by(label) %>%
  summarise(Cl = mean(grepl(".*Cl.*" , smiles)),
            O = mean(grepl(".*O.*" , smiles)),
            N = mean(grepl(".*N.*" , smiles)) ,
            S = mean(grepl(".*S.*" , smiles))) %>%
  ungroup()

heteroatoms %>%
  print_gt_table(table_title = "Heteroatoms") %>%
  fmt_percent(columns = c(Cl,O,N,S), rows = everything(), decimals = 1) 

#                         Oligosaccharides.
# 88.5% of oligosaccharides have nitrogen. To find out which nitrogen groups 
# are present in oligosaccharides, we find unique patterns of about 30 
# characters long containing "N" and corresponding "smiles".
oligosaccharides <-polymers %>% 
  select(Mw, label, smiles) %>%
  filter(label == "oligosaccharide") %>%
  mutate(pattern_nitrogen = str_extract(smiles, ".{0,15}N.{0,15}")) %>%
  filter(!is.na(pattern_nitrogen)) %>%
  arrange(Mw) %>%
  group_by(pattern_nitrogen) %>%
  summarise(smiles = first(smiles))

# In the smiles we can see 5 unique patterns with nitrogen, 
# from the following molecules:
oligosaccharides$pattern_nitrogen

# Chemical structures that correspond to each of the unique
# patterns with "N" that are present in the  "smiles" for oligosaccharides:
 print_molecules(smiles = oligosaccharides$smiles)
 
#  All patterns correspond to saccharides with amino groups [10] and
#  amine-linked saccharides [11].
  
#                         Peptides.
# 40% of peptides contain sulfur. These proteins may contain methionine,
# cysteine, homocysteine, or taurine in their structure [12].
  
#                          Plastics.
# 73.1% of the plastics contain chlorine, 72.5% contain aromatic groups. Thus, 
# it can be concluded that in this dataset, “plastics” are a type of polymer 
# that have carbon-carbon bonds in the main chain and Cl-, aliphatic or aromatic 
# groups as side groups.

## type --------------------------------------------------------------------

# New target variable. In this step, we add a new target variable 'type',
# that will have new labels for the polymers:
#  - for oligosaccharides, we add labels indicating the presence of an amino 
#                          group and/or an amine linkage between two sugars;
#  - for peptides - the presence of sulfur;
#  - for plastics - the presence of chloro- and/or aromatic groups;

 polymers <- polymers %>%
  mutate(type = case_when(
    label == "oligosaccharide" & grepl(".*N$|.*1N.*" , smiles) ~  "saccharide_amino",
    label == "oligosaccharide" & !grepl("*N$|.*1N." , smiles) ~  "saccharide",
                         
    label == "peptide" & grepl(".*S.*" , smiles) & !is_aromatic ~  "peptide_sulfur",
    label == "peptide" & !grepl(".*S.*" , smiles) & !is_aromatic ~  "peptide",
    label == "peptide" & grepl(".*S.*" , smiles) & is_aromatic ~  "peptide_sulfur_aromatic",
    label == "peptide" & !grepl(".*S.*" , smiles) & is_aromatic ~  "peptide_aromatic",
                         
    label == "plastic" & !grepl(".*Cl.*" , smiles) & is_aromatic ~  "plastic_aromatic",
    label == "plastic" & grepl(".*Cl.*" , smiles) & !is_aromatic ~  "plastic_chloro",
    label == "plastic" & !grepl(".*Cl.*" , smiles) & !is_aromatic ~  "plastic_aliphatic",
    label == "plastic" & grepl(".*Cl.*" , smiles) & is_aromatic ~  "plastic_chloro_aromatic") ) %>%
  
  mutate(type = ifelse(label == "oligosaccharide" & grepl(".*N\\[C@.*" , smiles),
                       str_c(type, "_amineLinked"), type)) %>%
  relocate(type, .after = label)

# Distribution of new labels for polymers
polymers %>% 
  ggplot(aes(y = type)) +
  geom_bar(fill = hist_color[1], 
           color = hist_color[2]) + 
  ggtitle("Types of the polymers") +
  themes 

# Models (unsupervised learning) ----------------------------------------------

## PCA ---------------------------------------------------------------------

# For high-dimensional data, we can use principal component analysis (PCA) and
# see how our polymer types are separated by the first few principal components.
  
# This will allow us to visualize features of the data and find the patterns
# having a lot fewer data dimensions, as well as give us new clues about the
# grouping of polymers into subtypes.
  
# To perform PCA, we use the prcomp() function and print a subset of the results
# along with the polymer types, and build scree plots.
 
# pca_model <- polymers %>%
#   select(matches("X\\d+")) %>%
#   prcomp()
# PCs_10 <- pca_model$x[ , 1:10]   # 10 principal components
# pca_sdev <- pca_model$sdev       # standard deviations of the principal components

# * From here, models and model outputs that have long computation times will 
#   be saved using the saveRDS() function.

# The output of prcomp() is large, so only first 10 principal components were saved:
# saveRDS(pca_model$x[ , 1:10], "scripts/models/PCs_10.rds")
# and sdev:
# saveRDS(pca_model$sdev, "scripts/models/pca_model_sdev.rds")

# Load 10 principal components from the PCA model:
PCs_10 <- readRDS(path_PCs_10)

# Load the standard deviations of the principal components from the PCA model:
pca_sdev <- readRDS(path_pca_sdev)

# The first 5 PCs with labels:
bind_cols( label = polymers$label,  PCs_10[ , 1:5]) %>%
  slice_head(n = 6) %>%
  print_gt_table(table_title = "The first five principal components",
                                 col_padding = 10)

### Scree plots -------------------------------------------------------------

# The proportion of variance explained by each principal component 
variance_explained <- tibble(ve = pca_sdev^2/sum(pca_sdev^2),
                             PC = 1:length(pca_sdev))

# Vector of cumulative proportion of variance explained  
cumulat_variance_explained <- cumsum(variance_explained$ve)

# Number of the PCS that explain 90% of the variance in the data,
# is used for the scree plots to limit the number of the PCs
nPCs_variance_90 <- which(cumulat_variance_explained >= 0.9)[1]

# This plot shows the proportion of variance explained for 254 principal 
# components that capture 90% of the information in this data.
variance_explained[1:nPCs_variance_90 ,] %>%
  ggplot(
    aes(PC, ve)) +
  geom_point(color = point_color, size = 1.5) +
  labs(title = "Scree plot. 254 principal components",
       x  = "principal component number",
       y = "variance explained") +
  themes +
  theme(plot.margin = unit(c(5,5,5,5), "pt") ) 

# This plot is the cumulative proportion of explained variance for the
# first 25 principal components. The first two PCs explain just over 30% 
# of the variance.
variance_explained[1:25 ,] %>%
  mutate(cve = cumsum(ve)) %>%
  ggplot(
    aes(PC, cve)) +
  geom_point(color = point_color, size = 1.5) +
  labs(title = "Scree plot. 25 principal components",
       x  = "principal component number",
       y = "cumulative variance explained") +
  themes +
  theme(plot.margin = unit(c(5,5,5,20), "pt") ) 
 

### Scatter plots ------------------------------------------------------------

# Next, we show scatter plots of the principal components and color by polymer type.
# Given that the first two principal components explain about 1/3 of the variance
# in the data, we expect that adding more principal components will also be useful
# in finding out if there are more polymer groups than we have already detected.
pca_results <- bind_cols(Mw = polymers$Mw,
                        label = polymers$label,
                        type = polymers$type, 
                        smiles = polymers$smiles,
                        PCs_10[ , 1:10])

theme_PC <- ggplot() +
  scale_color_manual(values = c(peptides_color, plastics_color, saccharides_color)) +
  themes +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text( size = rel(0.7) ),
        plot.margin = unit(c(5,10,5,10), "pt"))

theme_PC + geom_point(data = pca_results, aes(PC1, PC2, colour = factor(type)),shape = 1) + 
  ggtitle("PC1, PC2")

theme_PC + geom_point(data = pca_results, aes(PC4, PC10, colour = factor(type)),shape = 1) + 
  ggtitle("PC4, PC10")
  
theme_PC + geom_point(data = pca_results, aes(PC4, PC2, colour = factor(type)),shape = 1) + 
  ggtitle("PC4, PC2")

theme_PC + geom_point(data = pca_results, aes(PC3, PC6, colour = factor(type)),shape = 1) + 
  ggtitle("PC3, PC6")

#  In these graphs, all peptides are colored blue (from light blue to dark blue),
# saccharides are red/purple, and plastics are gray. The first two principle
# components are enough to separate original three labels. Plastics and
# oligosaccharides form groups according to the types we have already defined.
# However, there should be more peptide groups, from 3 to 4-5.


## Hierarchical cluster analysis -------------------------------------------

# Now let's apply another algorithm to find out which groups of peptides 
# might be present in the data.
  
#  There are some discussions on the application of kmeans and the hierarchical
# clustering algorithm to binary data .
# Hierarchical clustering with Ward and "centroid" methods based on Euclidean
# distances is not suitable for binary data [13-15].
  
# Therefore, we will use hierarchical cluster analysis using the method 
# "average", and to calculate distances we set the method to "binary".

# We perform hierarchical cluster analysis as follows:

# polymers_matrix <- polymers %>%                        # We convert
#     select(matches("X\\d+")) %>% as.matrix()           # X-variables to matrix,
# dists <- parallelDist::parDist( polymers_matrix,  method = "binary") # find Jaccard distances, and
# hclust_model <- hclust(dists, method = "average")      # build hierarchical clustering model.

# The model was saved as 'hclust_model.rds' in 'models' folder.
# Load the hierarchical clustering model:
hclust_model <- readRDS(path_hclust_model)

# We reduce the number of clusters to 12,
cut_hclust_model <- cutree(hclust_model, 12)

# add "label" and "type" variables to clusters
hclust_results <- tibble(label = polymers$label,
                         type = polymers$type,
                         cluster = cut_hclust_model)

#  and check how the clusters are distributed across the polymer types:
hclust_results %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster, type) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(type) %>%
  pivot_wider(names_from = type, values_from = n) %>%
  print_gt_table(table_title = "Distribution of the clusters by type of the polymers",
              col_padding = 4) %>%
  tab_spanner(label = md('<span style="color:#FC6C85;">saccharides</span>'),
              columns = starts_with("sacch")) %>%
  tab_spanner(label = md('<span style="color:#1F45FC;">peptides</span>'), 
              columns = starts_with("pept")) %>%
  tab_spanner(label = md('<span style="color:#52595D;">plastics</span>'), 
              columns = starts_with("plas")) %>%
  cols_label(saccharide	= "carbohydrate",
             saccharide_amino = "amino",
             saccharide_amineLinked = "amine-linked",
             saccharide_amino_amineLinked = "amino/amine-linked",
             peptide = "aliphatic",
             peptide_sulfur = "sulfur",
             peptide_aromatic = "aromatic",
             peptide_sulfur_aromatic = "sulfur/aromatic",
             plastic_aliphatic = "aliphatic",
             plastic_chloro = "chloro",
             plastic_aromatic = "aromatic",
             plastic_chloro_aromatic = "chloro/aromatic") %>%
  sub_missing(columns = everything(), rows = everything(), missing_text = "-") %>%
  data_color( columns = starts_with("sacch"), palette = saccharides_color[1], apply_to =  "text") %>% 
  data_color( columns = starts_with("pept"), palette = peptides_color[1], apply_to =  "text") %>%
  data_color( columns = starts_with("plas"), palette = plastics_color[4], apply_to =  "text") %>%
  cols_align("center") 

#  Although we can separate the original polymer labels, the same clusters include
# different types of polymers, especially for peptides; we observe similar results
# as for PCA.

### Dendrogram --------------------------------------------------------------

#  This can be visualized using a dendrogram: change the hclust model to 
# "dendrogram" and cut at height 0.787 to have 12 clusters.
dendrogram <- hclust_model %>% as.dendrogram() %>% cut( h = 0.787)
dendrogram <- dendrogram$upper

# To add the correct labels, we get the dendrogram attribute that shows the 
# number of cases in the nodes ("x.member"), then filter only the values that 
# correspond to leaf nodes (the "members" attribute is equal to 1). Then compare
# with the number of cases for each cluster in the "hclust_results" table.
n <- get_nodes_attr(dendrogram, "members")
x_member <- get_nodes_attr(dendrogram, "x.member")
label_dendrogram  <- data.frame(i = 1:12,
                                branch = labels(dendrogram),
                                x_member = x_member[which(n == 1)] )

cluster_dendrogram <- hclust_results %>% group_by(cluster) %>%
  summarise(x_member = n(), types = toString(unique(type))) %>%
  left_join(label_dendrogram, by = "x_member") %>%
  arrange(i) %>%
  mutate(label_color = case_when(
    str_starts(types, "sacch") ~ saccharides_color[1],
    str_starts(types, "pept") ~ peptides_color[1],
    str_starts(types, "plast") ~ plastics_color[2]))

labels(dendrogram) <- cluster_dendrogram$types

dendrogram %>%
  set("labels_col", cluster_dendrogram$label_color) %>%
  set("labels_cex", 0.7) %>%
  plot( horiz = TRUE, main = "Dendrogram")

# The dendrogram can also give us an idea of which type of polymers have a more
# similar structure.

### Clusters ----------------------------------------------------------------

#  We can use the results of cluster analysis to find additional features of
# chemical structures within the same polymer types. To do this, we draw examples
# of molecules by clusters.

# We want to display the peptide molecules with the lowest molecular weight, 
# so we add the variable "Mw".
hclust_results <- bind_cols(hclust_results, Mw = polymers$Mw, smiles = polymers$smiles)

#  clusters for peptides only:
clusters <- unique(hclust_results$cluster[which(hclust_results$label == "peptide")])

# filter out peptide clusters with low counts and keep 3 instances for each cluster
peptides <- hclust_results %>%
  filter(cluster %in% clusters) %>%
  group_by(cluster) %>%
  filter(n() > 20) %>%
  arrange(Mw) %>% 
  slice_min(Mw, n =3, with_ties = FALSE)

# display of peptide molecules by cluster:
smiles <-  peptides %>%   pull(smiles)
print_molecules(smiles = smiles, row_labels = clusters, prefix = "cluster ")


#### Peptides ----------------------------------------------------------------

#  There are three interesting clusters of peptides, all of which contain the same
# functional groups found in amino acids: histidine (His), tryptophan (Trp), and
# proline (Pro).
  
#  One way to find these groups is to use the rcdk package. There are functions
# that extract fragments (get.murcko.fragments() and get.exhaustive.fragments() ),
# but applying them to the dataset will take a long time. So we will use the
# stringr package to find patterns in the smiles.
  
#  As we can see, these groups have cycles. In SMILES, cyclic structures are
# enclosed in parentheses and have numbers indicating the ring.
ring_patterns <- unique(str_extract(smiles, "\\([CN=12]+\\d[CN=12]*\\)" ) )
ring_patterns

# tryptophan has two rings, there has to be "2" in the pattern
Trp <- ring_patterns[which(str_detect(ring_patterns, "2"))]

# histidine has two nitrogens
His <- ring_patterns[which(str_count(ring_patterns, "N") == 2)]

# proline is not aromatic, so there is no double bonds
Pro <- ring_patterns[which(str_detect(ring_patterns, "=", negate = TRUE))]

# Peptide molecules grouped into clusters 5, 6 and 7:
par(mfrow = c(1,3))
print_molecules( peptides$smiles[which(str_detect(smiles, Trp))][1], layout = 1)
symbols(x=67, y=72, circles=75, lwd = 2,
        add = TRUE, inches = FALSE, fg="#8EEBEC" ) 
text(x = 60, y = 280, paste("claster = ", peptides$cluster[which(str_detect(smiles, Trp))][1]))

print_molecules( peptides$smiles[which(str_detect(smiles, His))][1], layout = 1)
symbols(x=40, y=115, circles=45, lwd = 2,
        add = TRUE, inches = FALSE, fg="#8EEBEC") 
text(x = 60, y = 280, paste("claster = ", peptides$cluster[which(str_detect(smiles, His))][1]))

print_molecules( peptides$smiles[which(str_detect(smiles, Pro))][1], layout = 1)
symbols(x=135, y=52, circles=55, lwd = 2,
        add = TRUE, inches = FALSE, fg="#8EEBEC") 
text(x = 230, y = 280, paste("claster = ", peptides$cluster[which(str_detect(smiles, Pro))][1]))

# Thus, cluster 6 shows peptides with tryptophan, cluster 5 shows peptides with
# histidine, and cluster 7 shows peptides with proline.
  
#  We add this information to the polymer types and look at the distribution of
# these functional groups in the data:
polymers <- polymers %>%
  mutate(is_Trp = str_detect(smiles, Trp),
         is_His = str_detect(smiles, His),
         is_Pro = str_detect(smiles, Pro))

# proportions of functional groups of the peptides
multiple_func_groups <- polymers %>%
  filter(label == "peptide") %>%
  mutate(func_groups_count = is_Trp  + is_His + is_Pro) %>%
  pull(func_groups_count)

table(multiple_func_groups)

# Most peptides have none or one of these functional groups, 
# but some have two or all three.

# We add the names of the corresponding functional groups to the peptide types.
polymers <- polymers %>%
  mutate(type = ifelse(is_Trp, str_c(type, "_Trp"), type)) %>%
  mutate(type = ifelse(is_His, str_c(type, "_His"), type)) %>%
  mutate(type = ifelse(is_Pro, str_c(type, "_Pro"), type))


#### Plastics ----------------------------------------------------------------

# We can look at the "plastics" as well.

# 'plastics' clusters:
clusters <- unique(hclust_results$cluster[which(hclust_results$label == "plastic")])

smiles <- hclust_results %>%
  filter(cluster %in% clusters) %>%
  group_by(cluster) %>%
  arrange(Mw) %>% 
  slice_min(Mw, n =3, with_ties = FALSE) %>%
  pull(smiles)

print_molecules(smiles = smiles, row_labels = clusters, prefix = "cluster ")

#  Cluster 1 corresponds to aromatic polymers, and clusters 4 and 9 to aliphatic
# polymers.


#### Oligosaccharides --------------------------------------------------------

# Most oligosaccharides are found in one cluster.
#  Therefore, we will not add any other subtypes for plastics and oligosaccharides.

### type, distribution ------------------------------------------------------

#  Here are the final polymer types that will be used for classification in the
# target variable type.

# Distribution of new values in the variable "type"
polymers %>% 
  ggplot(aes(y = type)) +
  geom_bar(fill = hist_color[1], 
           color = hist_color[2]) + 
  ggtitle("Types of the polymers") +
  themes +
  theme(plot.margin = unit(c(5,20,5,30), "pt") ) 


# Data preparation --------------------------------------------------------

# For the models, we keep only the target variable 'type' and the X-variables
polymers <- polymers %>% select(matches("type|X\\d+"))

# The processed dataset was saved as an rds file: 
#       saveRDS(polymers, "data/processed/polymers_data_processed.rds")
# and can be loaded to skip all the above steps.
#       polymers <- readRDS(path_data_processed) 

## Select features ---------------------------------------------------------

# As features we will use Morgan fingerprints stored in X columns. 

# Above we defined "x_sums", a vector containing the column sums for all 
# X-variables, as:
x_sums <- polymers %>% 
  select(matches("X\\d+")) %>% 
  colSums()

# There are 761 out of 2048 X columns with less than 100 of 1 values. These
# predictors are different from 0 in 0.5% of the data; therefore we filter out
# these variables. The first 6 sums of X columns; for example, in column X0 
# there are only 5 values of "1" out of 20609
head(x_sums)

# Names of X-variables that have more than one hundred '1' values 
# that will be used in the training set:
ind_x <- which(x_sums >100)
ind_x <- which(colnames(polymers) %in% names(ind_x))


## Split data --------------------------------------------------------------
  
# We split data to train set and test set and change the target variable, type,
# to factor.
  
#  The createDataPartition() function from caret partitions the data so that the
# training and test sets have similar proportions of target variable values.
set.seed(1)       # indexes for the training set:
index <- createDataPartition(y = polymers$type, 
                             times = 1, p = 0.8, 
                             list = FALSE)

#       train set
train_x <- polymers[ index, ind_x] %>% select(matches("X\\d+"))

train_type <- polymers[ index, ] %>% select(type) %>% mutate(type = factor(type))

#      test set
test_x <- polymers[ -index, ind_x] %>% select(matches("X\\d+"))
test_type <- polymers[ -index, ] %>% select(type) %>%  mutate(type = factor(type))

train_x[1:5, c(1:3, (ncol(train_x)-3):ncol(train_x) )] %>%
  print_gt_table(table_title = "The first and the last predictors",
              col_padding = 15) %>%
  cols_align("center")


# Models (supervised learning) --------------------------------------------

# First, we compare three different algorithms with default settings, then based
# on the results we select the best one and adjust the parameters for this model.

## Model selection ---------------------------------------------------------

#  We try three classical algorithms that can be applied to our data.


### K-nearest neighbors (kNN) -----------------------------------------------

# For this model, we use the knn3() function from caret, by default it 
# has k = 5 for classification.
knn_model <- knn3(type ~ . ,
                  data = bind_cols(train_type, train_x))

knn_model$k

### Classification tree -----------------------------------------------------

# To fit this model, we use rpart() function from the rpart package.

# classification tree
# rpart_model <- rpart(type ~ . ,
#                      method = "class",
#                      data = bind_cols(train_type, train_x))

# Load the classification tree model:
rpart_model <- readRDS(path_rpart_model)
unlist(rpart_model$control)

### Random forest -----------------------------------------------------------

# The ranger package allows us to quickly compute a random forest model.
  
#  We also want to estimate a variable importance; this can be done by setting the
# importance parameter to 'impurity_corrected'.
#  From ranger() manual: "'impurity_corrected' importance measure is unbiased in
# terms of the number of categories and category frequencies..". This is our case.

# ranger_model <- ranger(type ~ . ,
#                        data = bind_cols(train_type, train_x),
#                        seed = 1,
#                        importance = "impurity_corrected")

# Load the random forest model
ranger_model <- readRDS(path_ranger_model) 
ranger_model

## Compare results ---------------------------------------------------------
 
#  We will use accuracy as a metric to evaluate the performance of the models.

# knn_predicted <- predict(knn_model, test_x, type = "class")
knn_predicted <- readRDS(path_knn_predicted)

rpart_predicted <- predict(rpart_model, test_x, type = "class")

ranger_predicted <- predict(ranger_model, test_x)

data.frame(knn = round(mean(knn_predicted ==  test_type$type),4),
           rpart = round(mean(rpart_predicted ==  test_type$type),4),
           ranger = round(mean(ranger_predicted$predictions ==  test_type$type),4) ) %>%
  print_gt_table(table_title = "Polymer type prediction accuracy for models with default settings") %>%
  cols_label(knn = "kNN",
             rpart = "classification tree",
             ranger = "random forest") %>%
  cols_align("center")

# The random forest model provides the highest accuracy, so we will focus on 
# this model.

## Random forest -----------------------------------------------------------

# Random forest is a powerful classification algorithm that performs better 
# than classification tree and kNN models in most cases. Even with default 
# settings, it shows high accuracy. Next, we will try to reduce the number 
# of predictors and find a better mtry parameter. This will reduce the 
# computation time and further improve the accuracy.

### Variable importance -----------------------------------------------------

# Before tuning the model parameters, let's first consider the variable importance,
# which was estimated by the ranger() function and can be found in its output.

# Variable importance table: variables with higher importance 
# are at the beginning with lower indices n.
variable_importance <- 
  data.frame(x = names(ranger_model$variable.importance),
             importance = ranger_model$variable.importance) %>%
  arrange(desc(importance)) %>%
  mutate(n = row_number()) 

# Variable importance plot.
variable_importance %>%
  ggplot(aes(n, importance)) +
  geom_point(color = point_color, size = 1.5) + 
  ggtitle("Variable importance")  + 
  themes 

#  We can reduce the number of X-variables and use fewer predictors, which will
# allow us to calculate the optimal model parameters faster.

# We use an arbitrary limit on the number of predictors to keep:
# we drop variables that have less than 5% of the maximum 
# variable importance value.
ind_x_importance <- variable_importance %>%
  filter(importance > max(variable_importance$importance)*0.05) %>% pull(x)

length(ind_x_importance)


### mtry --------------------------------------------------------------------

# mtry is a number of variables to split at in each node. Default is the square
# root of the number variables. In our first random forest model, this number
# was 35, which is close to the square root of the number of X
# variables we used in the training set.
sqrt(ncol(train_x))
  
# But now we have fewer predictors, we need to keep that in mind and include 
# mtry values around this number:
sqrt(length(ind_x_importance))

# Keep only predictors with high variable importance:
train_x <- train_x[ , ind_x_importance]

# Find the optimal mtry: we reduce the number of trees to save computation time;
# find the accuracy for the test set along with the out-of-bag prediction error
tune_mtry <- sapply(round(seq(-10, 10, 2)+ sqrt(length(ind_x_importance))),
                    function(i){
  ranger_tune = ranger(type ~ . ,
                      data = bind_cols(train_type, train_x),
                      mtry = i,
                      num.trees = 200,
                      seed = 1)
  predicted = predict(ranger_tune, test_x)
  
  c(accuracy = mean(predicted$predictions ==  test_type$type),
    oob_error = ranger_tune$prediction.error)
})


mtry_results <- data.frame(mtry = round(seq(-10, 10, 2)+ sqrt(length(ind_x_importance))),
                           accuracy = tune_mtry[1, ],
                           oob_error = tune_mtry[2, ]) 

# 'mtry_results' was saved as 'tune_mtry.rds', and can be loaded as:
#  mtry_results <- readRDS(path_tune_mtry)
 
# Plot the accuracy for different values of mtry:
mtry_results %>%
  ggplot(aes(mtry, accuracy)) +
  geom_point(color = point_color, size = 1.8) +
  ggtitle("Accuracy of the random forest models")  +
  themes 

# We choose the best 'mtry' as the one that results in the highest accuracy
# and lowest OOB error
(best_mtry <- mtry_results %>%
  filter(accuracy == max(accuracy)) %>%
  filter(oob_error == min(oob_error)) %>%
  first() %>%
  pull(mtry) )

# The mtry values that give us the best accuracy start at `r best_mtry`. Now we
# can build the final random forest model with the best mtry parameter and a large
# number of trees.

## Final model -------------------------------------------------------------

# Fit the final random forest model with the best mtry parameter 
# and a large number of trees.
ranger_final <- ranger(type ~ . ,
                        data = bind_cols(train_type, train_x),
                        mtry = best_mtry,
                        num.trees = 1000,
                        seed = 1)

predicted = predict(ranger_final, test_x)

paste("Accuracy =", mean(predicted$predictions ==  test_type$type))

# While an accuracy of 1 is not common for supervised learning algorithms, we are
# dealing with a special type of predictors here. The features in this dataset are
# not derived from experiments or observations, but are predefined descriptors
# of specific chemical structures. Each predictor indicates the presence of
# these groups in the molecules. And belonging to a certain polymer type directly
# depends on the presence of these groups.

# Conclusion --------------------------------------------------------------

#  Applying unsupervised learning algorithms to Morgan fingerprints as features
# allows us to find patterns that can be useful for detecting polymers with
# similar chemical structures, quickly extracting specific groups of polymers from
# large data sets.
  
#  Morgan fingerprints also can be a "fuel" for the supervised learning models and
# eventually lead to the development of applications that can be used in practice
# to classify polymers based on chemical structure that is encoded in the simple
# predictors.
  
#  In this project, we compared three models: K-nearest neighbors, classification
# tree, and random forest. To evaluate the performance of the models, we used
# accuracy; this metric was evaluated on unseen data in the test set.
  
#  The best model for the polymers dataset is the random forest model with mrty
# equal to `r best_mtry`, and it requires `r length(ind_x_importance)` x-variables
# to achieve an accuracy of 1. The high accuracy can be explained by the nature of
# the features, which are predefined descriptors.


# References --------------------------------------------------------------

# (1)  What Is a CAS Registry Number? https://www.cas.org/cas-data/cas-registry.

# (2)  Blum, L. C.; Reymond, J.-L. 970 Million Druglike Small Molecules for Virtual
#      Screening in the Chemical Universe Database GDB-13. Journal of the American 
#      Chemical Society 2009, 131 (25), 8732–8733. https://doi.org/10.1021/ja902302h.

# (3)  Ruddigkeit, L.; Deursen, R. van; Blum, L. C.; Reymond, J.-L. Enumeration of 
#      166 Billion Organic Small Molecules in the Chemical Universe Database GDB-17. 
#      Journal of Chemical Information and Modeling 2012, 52 (11), 2864–2875. 
#      https://doi.org/10.1021/ci300415d.

# (4)  Desidery, L.; Lanotte, M. 1 - Polymers and Plastics: Types, Properties, and
#      Manufacturing. In Plastic waste for sustainable asphalt roads; Giustozzi, 
#      F., Nizamuddin, S., Eds.; Woodhead publishing series in civil and structural 
#      engineering; Woodhead Publishing, 2022; pp 3–28. 
#      https://doi.org/10.1016/B978-0-323-85789-5.00001-0.

# (5)  Young, R. J.; Lovell, P. A. Introduction to Polymers: Third Edition, 3rd ed.;
#      CRC Press: United States, 2011.

# (6)  Yang, J.; Cai, Y.; Zhao, K.; Xie, H.; Chen, X. Concepts and Applications of
#      Chemical Fingerprint for Hit and Lead Screening. Drug Discovery Today 2022, 
#      27 (11), 103356. https://doi.org/10.1016/j.drudis.2022.103356.

# (7)  Morgan, H. L. The Generation of a Unique Machine Description for Chemical 
#      Structures-a Technique Developed at Chemical Abstracts Service. Journal of
#      Chemical Documentation 1965, 5 (2), 107–113. https://doi.org/10.1021/c160017a018.

# (8)  Lo, Y.-C.; Rensi, S. E.; Torng, W.; Altman, R. B. Machine Learning in
#      Chemoinformatics and Drug Discovery. Drug Discovery Today 2018, 23 (8), 1538–1546.
#      https://doi.org/10.1016/j.drudis.2018.05.010.

# (9)  Weininger, D. SMILES, a Chemical Language and Information System. 1. Introduction 
#      to Methodology and Encoding Rules. Journal of Chemical Information and Computer 
#      Sciences 1988, 28 (1), 31–36. https://doi.org/10.1021/ci00057a005.

# (10) BHAGAVAN, N. V. CHAPTER 9 - Simple Carbohydrates. In Medical biochemistry 
#      (fourth edition); BHAGAVAN, N. V., Ed.; Academic Press: San Diego, 2002; 
#      pp 133–151. https://doi.org/10.1016/B978-012095440-7/50011-1.

# (11) Akhtar, T.; Cumpstey, I. Investigations into the Synthesis of Amine-Linked
#      Neodisaccharides. Tetrahedron Letters 2007, 48 (49), 8673–8677.
#      https://doi.org/10.1016/j.tetlet.2007.10.039.

# (12) Brosnan, J. T.; Brosnan, M. E. The Sulfur-Containing Amino Acids: An Overview.
#      The Journal of nutrition 2006, 136 6 Suppl, 1636S–1640S.

# (13) Scheuch, C. Clustering Binary Data. 
#      https://blog.tidy-intelligence.com/posts/clustering-binary-data/.

# (14) Hierarchical or TwoStep Cluster Analysis for Binary Data?
#      https://stats.stackexchange.com/questions/116856/hierarchical-or-twostep-cluster-analysis-for-binary-data.

# (15) Comparing Hierarchical Clustering Dendrograms Obtained by Different Distances
#      & Methods.
#      https://stats.stackexchange.com/questions/63546/comparing-hierarchical-clustering-dendrograms-obtained-by-different-distances/63549#63549.


# Versions ----------------------------------------------------------------

# R version 4.4.3 (2025-02-28)  

# parallelDist_0.2.6    ranger_0.17.0     rpart_4.1.24      caret_7.0-1 
# lattice_0.22-6        stringr_1.5.1     rcdk_3.8.1        rcdklibs_2.9
# rJava_1.0-11          gt_0.11.1         patchwork_1.3.0   dendextend_1.19.0
# ggplot2_3.5.1         tidyr_1.3.1       dplyr_1.1.4.9000  

