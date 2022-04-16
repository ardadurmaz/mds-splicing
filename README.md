# mds-splicing

- This repository stores the scripts required for preprocessing and integrative modeling of gene expression, mutation and percent-spliced in (psi) data of MDS cases.
- The data subfolder holds sample datasets in the form of observations (rows) and features (columns) which can be used as a template for data input.
- "latent_model_vae.py" is the script for the latent model based on variational inference and structured similar to CCA. Code is tested and run on Tensorflow 2.0.
- "plot_consensus.py" script is used to generate unsupervised clustering across different number of clusters (K) using Gaussian Mixture Models (GMM) and keeps track of BIC.
- In order to process the workflow with new datasets, "latent_model_vae.py" followed by "plot_consensus.py" scripts can be run with appropriate input data structures and parameters.
