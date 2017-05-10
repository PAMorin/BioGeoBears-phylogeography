#######################################################
# How (and whether) to collapse tips
# http://phylo.wikidot.com/example-biogeobears-scripts#pruning_a_tree
#######################################################
# People with trees that mix between-species and within-species variation 
# often need to reduce their tree to have just a tree of just species 
# or just monophyletic populations.
#
# This is because BioGeoBEARS (and Lagrange etc.) assume that 
# your phylogeny is a phylogeny of species/monophyletic populations, 
# rather than a tree of specimen samples.
#
# Options include:
# 
# 1. Re-do the analysis with just one OTU per species/monophyletic population 
#    (OTU = "operational taxonomic unit" = tip in a tree)
# 
# 2. Do a full gene tree/species tree analysis, e.g. with starBEAST
# 
# 3. Do a somewhat crude pruning in R, after the fact.
#
# Please note that, officially, it is a somewhat bad idea to estimated a dated 
# tree in BEAST with a mixture of within-species and between species OTUs (this 
# is the precursor to #3).  
#
# This is because the standard BEAST analysis will assume that the tree process is 
# EITHER phylogenetic (when you use a pure-birth (Yule) *OR* birth-death tree prior) 
# or population genetic (when you use a coalescent prior).  When your data mixes 
# these, there is the potential that some part of the tree dating is getting 
# screwed up by having an incorrect prior.
# 
# In my experience, it doesn't seem to matter much, at least for cases where the 
# data are a decent mix of within- and between-population data, there are a number
# of calibration points, and a phylogenetic tree prior is being used. 
# 
# I have seen bad dating problems, though, if there are few calibrations, and the data 
# are mostly within one coalescing population, and a birth-death prior is used.  This 
# pushed the divergence times way too deep and spread out the nodes way too much, 
# compared to the quite clocklike pattern (all the tips about the same height) and 
# coalescent pattern (a large clade of sequences, with very recent divergence times)
# in the undated tree.
# 
# Anyway, reviewers may complain if you use strategy #3. It's up to you to make a 
# good decision about these strategies and decide how to defend it in your paper.
# 
# All that said, here's the code for pruning trees.
###############################################################################
rm(list = ls())
options(stringsAsFactors = FALSE)

library(ape)    # for read.tree, drop.tip, etc.
library(gdata)     # for trim

library(BioGeoBEARS)    # for prune_specimens_to_species etc.

workdir <- "/Users/phil.morin/Documents/Mol\ Ecol\ Lab/\ Projects/\ Sperm\ whales/Pmac_mitogenomics/Phylogeog_BioGeoBears_results"
wd = np(workdir)
#PM changed to folder in the BioGeoBears folder; change to appropriate 
#folder in project directory by copying path from Termina window.
setwd(wd)
# Double-check your working directory with getwd()
getwd()

# REQUIRED patches for this to work
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")

# Here is a newick string (you could also load a tree from a newick file with read.tree)
#newick_string = "(((sp1:0.03610559932,sp2:0.03610559932)1.0000:0.1776629349,((((sp3:0.008442378445,sp4:0.008442378445)1.0000:0.09122692306,sp5:0.09966930151)0.8380:0.05350909945,sp6:0.153178401)0.7060:0.04118371017,((sp20:0.1646354697,(((((sp21:0.05541959218,sp22:0.05541959218)0.7220:0.02482956203,sp23:0.08024915421)0.3250:0.01090832874,(sp24:0.07252867627,sp25:0.07252867627)0.5000:0.01862880668)0.3060:0.01016242573,sp26:0.1013199087)0.8940:0.05355067998,(sp27:0.1465197179,(sp28:0.1343470201,(sp29:0.1300709381,(((sp30:0.02406821269,sp31:0.02406821269)1.0000:0.03706382207,sp32:0.06113203476)1.0000:0.06390204876,((((sp33:0.02740065156,sp34:0.02740065156)0.6250:0.007508656628,sp35:0.03490930819)0.8450:0.02363732707,sp36:0.05854663526)1.0000:0.06308298802,(sp37:0.09852184783,(sp38:0.01910350891,(sp39:0.007259454956,sp40:0.007259454956)0.8450:0.01184405395)1.0000:0.07941833892)0.8610:0.02310777545)0.3450:0.003404460232)0.1360:0.005036854609)0.1050:0.004276082013)0.2300:0.0121726978)0.1700:0.008350870719)0.4160:0.009764880993)0.7460:0.02738230786,(sp19:0.1322081512,(sp18:0.1108740943,(sp17:0.09751406895,(sp16:0.07347283163,(sp15:0.05565075583,(sp14:0.04956539576,(sp13:0.04331045159,(sp12:0.02278042935,(sp11:0.01471299318,((sp9:0.006538713623,sp10:0.006538713623)0.7890:0.007881519684,(sp7:0.006528757196,sp8:0.006528757196)0.8350:0.007891476111)0.1960:0.0002927598732)0.6500:0.008067436169)1.0000:0.02053002224)0.4060:0.00625494417)0.2690:0.006085360071)0.4080:0.0178220758)0.9690:0.02404123732)0.8090:0.0133600254)0.8600:0.02133405685)0.9990:0.05980962632)0.2810:0.002344333608)0.5880:0.01940642309):0.7862314658,(sp41:0.05865850646,(sp42:0.01278826063,sp43:0.01278826063)1.0000:0.04587024583):0.9413414935);"

tr = read.tree(file="All_Pmac_Unique_hap_seqs_strict_100M_rename.newick")
tr
# Look at the tip labels
tr$tip.label
#write.csv(tr$tip.label, "old_tipnames.csv")
# open file and delete first column, first row, to just leave the tip labels

Pmac_old_tipnames <- read.csv("old_tipnames.csv", header = FALSE)
# Make a list of the tip labels
# (you could also load this from a text or Excel file)

# Make a list of the NEW tip labels. This must be in 
# the SAME ORDER as the tips in "old_tipnames"
# (it doesn't matter if the tips are in the same 
#  order as the tree)
#
#

Pmac_new_tipnames <- read.csv("New_tipnames2.csv", header = FALSE)

# Make a data.frame
OTUs = Pmac_old_tipnames
species = Pmac_new_tipnames
tmpmat = cbind(OTUs, species)
df = as.data.frame(tmpmat, stringsAsFactors=FALSE)
names(df) = c("OTUs", "species")

# Look at the beginning of the table
head(df)

library(gdata) # for trim

# We'll make a PDF so you can see all the steps
pdffn = "pruning_PDF.pdf"
pdf(pdffn, height=6, width=6)

result = prune_specimens_to_species(original_tr=tr, xls=df, group_name="default", 
                          titletxt="", areas_abbr=NULL, plot_intermediate=TRUE)

# stop PDF and then open it
dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

# Here's the new tree
pruned_tree = result$tr
# We can write it to a file
write.tree(pruned_tree, file="pruned_tree2.newick")

# Look at the Newick file
moref("pruned_tree2.newick")

# Look at the tip labels
pruned_tree$tip.label
write.csv(pruned_tree$tip.label, "pruned_tree_tip_labels2.csv")
# open file and delete first column, first row, to just leave the tip labels



################
# didn't continue to do the geography below

#######################################################
# If you need to merge geography, you can do that also
#######################################################

# These also have to be in the same order
# (this is all easier to do in Excel with 
#  gdata::read.xls -- but this is just 
#  showing the concept)
region = c("A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A|B",
           "C|D",
           "E",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "A",
           "F",
           "C",
           "A",
           "A",
           "A")

# Make a data.frame
OTUs = old_tipnames
species = new_tipnames
tmpmat = cbind(OTUs, species, region)
df = as.data.frame(tmpmat, stringsAsFactors=FALSE)
names(df) = c("OTUs", "species", "region")

# Look at the beginning of the table:
head(df)

# We'll make a PDF so you can see all the steps
pdffn = "pruning_PDF2.pdf"
pdf(pdffn, height=6, width=6)

result2 = prune_specimens_to_species(original_tr=tr, xls=df, group_name="default", titletxt="", areas_abbr=NULL, plot_intermediate=TRUE)

# stop PDF and then open it
dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

# Here's the new tree
pruned_tree = result2$tr

# We can write it to a file
write.tree(pruned_tree, file="pruned_tree.newick")

# Look at the Newick file
moref("pruned_tree.newick")

# This time we have a geography object also
tipranges = result2$tipranges

# Write the tipranges to a text file
lgdata_fn = "pruned_geog.data"
save_tipranges_to_LagrangePHYLIP(tipranges=tipranges, lgdata_fn=lgdata_fn)

# Look at the text file
moref(lgdata_fn)