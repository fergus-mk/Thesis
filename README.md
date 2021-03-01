# Thesis

This is a repository for my Earth Science MSc thesis. The thesis is an attempt to the quality of tagged
flood defences in OpenStreetMap for the purpose of hydrological modelling.

_I will include a fuller descritpion at a later stage, I will also upload neater versions of the code if I have time_

__Known Issues/ Improvements Recommended__

If someone found this topic interesting and wanted to take it further I think it needs some improvments:

- Code needs to be made neater/ more in line with PEP 8 (this is my first project and I need to improve this)
- More data points should be used for the machine learning part of the project(e.g. total 1,000 or more)
- The GFP250M (flood plain) data euclidean distance should be recalculated with a higher resolution
- Additional input datasets (such as one inorporating urban area/ a drainage basin dataset should be used)
- The scikit-learn component should be improved with a. better bootstrapping, and b. hyperparameter optimisation 


The files uploaded so far are:

- The initial query to access flood defences from Geofabrik - GeofabrikQuery.ipynb
- Calculating summary statistics for classified flood defence points - SummaryStats.py
- Using these stats to classify the entire datset - ClassPointsFinal.ipynb

_Thanks for visiting this repository_

Fergus
