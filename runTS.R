#-----------------------------------------------------------------------------#
# Parameters for the analysis, uncomment and edit them if you prefer to leave
# it automated for future uses.


#taxon ids to sequence names
#idsFile <- "inputs/test1"
idsFile <- "data/validation/taxonIDs_2_sequenceIDs.txt"

#fasta files
#multifasta <- "inputs/multifasta"
multifasta <- "data/validation/sequences_with_taxonIDs.fasta"

#NCBI taxonomy files
taxondir <- "taxdump/"

#number of sequences to get
m <- 100

#method to use (either 'diversity' or 'balance')
method <- "diversity"

#wheter to randomize 
randomize <- "yes"

#allow TS to repeat IDs if needed? ('no' is better to get a higher diversity)
replacement <- "no"

ignoreIDs <- NULL
requireIDs <- NULL
ignoreNonLeafID <- NULL
outFile <- "output.fasta"

#-----------------------------------------------------------------------------#

source("bin/TaxonSampling/TaxonSampling.R")

#library("CHNOSZ")
#library("ape")

nodes <- suppressMessages(getnodes(taxondir))
countIDs <- TS_TaxonomyData(idsFile, nodes)

nodes <- Simplify_Nodes(nodes, countIDs)
outputIDs <- TS_Algorithm(1, m, nodes, countIDs, method, randomize,
                          replacement, ignoreIDs, requireIDs,
                          ignoreNonLeafID)

WriteFasta(idsFile, multifasta, outputIDs, outFile)
