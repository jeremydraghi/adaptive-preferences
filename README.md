# adaptive-preferences
Code for "Asymmetric evolvability leads to specialization without trade-offs"

This repository contains the basic code used for all simulations: an R script named "revised_parallel_model.R" and C++ functions in the file "pickParents.cpp". The basic model is that a numerical experiment is designed and executed in parallel from the R script, while the processor-heavy tasks of choosing parents for each generation is conducted with faster C++ functions.

While most of the experimental design is intuitively automated, there are several features which rely on manual code changes. In particular, selection of parameter values for each replicate in lines 64-71 would need to be modified to account for changes in which parameters are varied across a set of replicates. Also, the function pickParentsFreePref() is not controled by a varible, but can be commented in/out as an option. 
