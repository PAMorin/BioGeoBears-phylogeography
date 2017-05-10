# BioGeoBears-phylogeography
Scripts for using the BioGeoBears for phylogeographic model selection and inference. 

Collapse_tree_tips.R: BioGeoBears assumes that 
your phylogeny is a phylogeny of species/monophyletic populations, 
rather than a tree of specimen samples. If you have a larger tree of specimens or
haplotypes, it can be reduced to a tree of suppored clades/species/populations by 
pruning with this script.

BioGeoBears_Pmac_test.R: example script for testing phylogeography. Full example available in folder: Pmac_Pruned_tree_12pops_example