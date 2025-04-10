Started 18 March 2025  
By Victor

#### mastree folder: 
- mastree_subset_fagus.csv: subset of the full `MASTREE+` dataset, used in Journé et al. 2024 ( https://doi.org/10.1038/s41477-024-01651-w), and downloaded from https://osf.io/s2cd4 (accessed 18 Mar 2025)  
  From README file in OSF repo:
  > Last version of data was available here: https://github.com/JJFoest/MASTREEplus/tree/main/Data (MASTREEplus_2022-02-03_V1.csv)
  > Data have been restricted here to > 14 years of observation, we removed, "flower", "index", "pollen" data and kept only continous data (VarType == "C"). We kept observations after 1951
  > All columns names are reported in table 1 of https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16130  
  > \+ new columns names here:  
  > sitenewname: here, we created a new column name based on (Alpha_Number, "\_",Site_number, "\_", Species_code), to avoid any redundancy in site names.  
n: new length time series based on the observation grouped by sitenewname. (ie number of observation for each sitenewname factor group)  
Date2: conversion of Date (we used here strptime() function from R to convert Date from character)

### ebms folder
- Ebs stand for "English Beech Masting Survey"
- Come from here: https://osf.io/dqwzn (paper: https://www.pnas.org/doi/10.1073/pnas.2423181122)
- Note: this is NOT seed trap data
The ground below each tree was searched for seeds for 7 min, and the collected seeds were later sorted and counted in the laboratory