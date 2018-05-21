
GENAVi (Gene Expression Normalization Analysis and Visualization) is an rshiny web application that provides a GUI based platform for the analysis of gene expression data. 

GENAVi combines several R packages commonly used for normalizing, clustering, visualizing, and performing differential expression analysis (DEA) on RNA-seq data.

Within our application we have also included an RNA-seq data set on a panel of 20 cell lines commonly used for the study of breast and ovarian cancer.
This dataset can serve as an introduction to the various functions of GENAVi.
In addition to being able to query and analyze the dataset we provided, users can also upload their own gene expression count matrix and apply all of the functions within GENAVi to analyze their own data.

GENAVi is hosted on the junkdnalab rshiny server which can be accessed in any browser with internet connection through this link: https://junkdnalab.shinyapps.io/GENAVi/

Additionally, GENAVi can be run on a users local machine or server by entering the command: shiny::runGitHub("alpreyes/GENAVi") in an RStudio window. 
With this method, all necessary R packages will be automatically installed, all source code used to build GENAVi will be downloaded directly from the GitHub repository, and the application will be hosted and run locally.

Lastly, the user can also run GENAVi on a local machine or server via a docker image through this link: https://hub.docker.com/r/cedarscompbio/genavi/

We recommend running GENAVi on a local machine or server when uploading a large dataset that would require a more powerful server than the junkdnalab rshiny server for analysis.
