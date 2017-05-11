# BioGeoBears-phylogeography
Scripts for using the BioGeoBears for phylogeographic model selection and inference. 

Collapse_tree_tips.R: BioGeoBears assumes that 
your phylogeny is a phylogeny of species/monophyletic populations, 
rather than a tree of specimen samples. If you have a larger tree of specimens or
haplotypes, it can be reduced to a tree of suppored clades/species/populations by 
pruning with this script.

BioGeoBears_Pmac_test.R: example script for testing phylogeography. Full example available in folder: Pmac_Pruned_tree_12pops_example

References:
Matzke NJ (2012) Founder-event speciation in BioGeoBEARS package dramatically improves likelihoods and alters parameter inference in Dispersal-Extinction-Cladogenesis (DEC) analyses. Frontiers of Biogeography 4(suppl. 1), 210.

Matzke NJ (2013) BioGeoBEARS: BioGeography with Bayesian (and Likelihood) Evolutionary
Analysis in R Scripts, University of California, Berkeley.

Matzke NJ (2014) Model Selection in Historical Biogeography Reveals that Founder-Event Speciation Is a Crucial Process in Island Clades. Systematic Biology 63, 951-970.
