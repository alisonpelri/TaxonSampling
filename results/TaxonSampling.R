
# install.packages("CHNOSZ")
library("CHNOSZ")

idsFile <- "../KOMODO2/projects/ast/tax_ids_non_redundant.txt"
taxondir <- "../KOMODO2/taxdump"
nodes <- suppressMessages(getnodes(taxondir))

# -=-=-=- Taxonomy Sampling (TS) -=-=-=-


# TODO: adapt ignoreIDs (or create a new parameter) for when two taxa are 
#       parent/child and we want to ignore only the parent.
# TODO: test ignoreNonLeafID parameter


TS_TaxonomyData <- function(idsFile, nodes) {
  # Evaluates the taxonomy structure from the input IDs and returns the count
  # of each taxon ID above zero in the taxonomy graph.
  #
  # Args:
  #   idsFile: (char/int) either a path to a file with the input taxon IDs, or
  #                       a vector with the IDs themselves, in either integer
  #                       or character format.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure. Created by getnodes() from the CHNOSZ
  #                       package.
  # Returns:
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  
  # Get ids to search
  if (!is.null(idsFile) & isTRUE(file.exists(as.character(idsFile)))) {
    ids <- read.table(idsFile, sep = "\t", header = FALSE, comment.char = "",
                      colClasses = "integer", strip.white = TRUE)
    ids <- ids[, 1]
  } else {
    ids <- as.integer(idsFile)  # assume it's a test.name or back.name, for now
  }
  
  # Sanity check: filter IDs that aren't part of NCBI notation.
  if (!all(is.element(ids, nodes$id))) {
    cat("Warning: the following inputs are not part of NCBI taxonomy IDs",
        "and will be ignored.\n",
        ids[!is.element(ids, nodes$id)], "\n")
    ids <- ids[is.element(ids, nodes$id)]
  }
  
  # Sanity check: unique ID inputs, remove duplicates.
  if (any(duplicated(ids))) {
    cat("Warning: some ids are repeated, using only one instance of each.\n",
        ids[duplicated(ids)], "\n")
    ids <- unique(ids)
  }
  
  # Counting how often each taxonomy ID occurs in the dataset
  countIDs <- rep.int(0, nrow(nodes))
  names(countIDs) <- nodes$id
  nodes_parent <- nodes$parent
  names(nodes_parent) <- nodes$id
  searchIDs <- ids
  countIDs[as.character(searchIDs)] <- 1
  while (length(searchIDs) > 0) {
    searchIDs <- nodes_parent[as.character(searchIDs)]
    parentage <- table(searchIDs)
    countIDs[names(parentage)] <- countIDs[names(parentage)] + parentage
    searchIDs <- searchIDs[searchIDs != 1]
  }
  
  countIDs <- countIDs[countIDs > 0]
  return(countIDs)
}



TS_Algorithm <- function(taxon, m, nodes, countIDs, replacement = "no",
                         method = "balanced", strict = "no",
                         requireIDs = NULL) {
  # An algorithm that receives a group of Taxonomy IDs and the size m of the
  # sample to obtain from them. Returns a vector with a maximized taxonomy 
  # diversity. Assumes that every input id is unique.
  # If replacement = "yes", then it may return repeated IDs if it increases
  # the taxonomy diversity when strict = "yes".
  #
  # Args:
  #   taxon: (char/integer) Taxon from which to start sampling children taxa.
  #   m: (integer) size of the sample to generate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure. Created by getnodes() from the CHNOSZ
  #                       package.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   replacement: (char) whether the algorithm allows to repeat IDs in order
  #                       to maximize taxonomy diversity and to reach m IDs
  #                       in the output.
  #   method: (char) whether the algorithm will choose IDs randomly or 
  #                  maintaining a balanced allocation (m_i differing by no
  #                  more than 1 if the maximum possible value wasn't reached).
  #   strict: (char) whether it favors balanced taxa representation 
  #                  (strict = "yes") or maximized taxa representation ("no").
  #   requireIDs: (char) IDs that must appear in the output.
  # Returns:
  #   outputIDs: (vector) vector of IDs with maximized taxonomy diversity.
  
  # Sanity check
  if (m <= 0) {
    Print("Error: m less or equal than zero during recursion.\n")
    exit(0)
  }
  
  # First step: Find the sub-taxa (children nodes) of the current taxon
  taxon <- as.integer(taxon)
  children <- nodes$id[nodes$parent == taxon & nodes$id != taxon]
  children <- intersect(children, names(countIDs))
  
  # Condition to end recursion
  if (length(children) == 0) {
    if (replacement == "no") {
      return(as.character(taxon))
    } else {
      return(rep(as.character(taxon), m))
    }
  }
  
  childrenCount <- countIDs[as.character(children)]
  
  # For cases when one taxon isn't a leaf node, but is an input ID that should
  # be available to sample along its children.
  # Example of this case is an input with Homo sapiens (9606), 
  # H. sapiens neanderthalensis (63221) and H. sapiens ssp. denisova (741158).
  if (sum(childrenCount) < countIDs[as.character(taxon)]) {
    childrenCount <- c(childrenCount, countIDs[as.character(taxon)])
    childrenCount[as.character(taxon)] <- 1
  }
  
  m_i <- rep.int(0, length(childrenCount))
  names(m_i) <- names(childrenCount)
  
  # Allocating m_i to the children taxa at last.
  # Relaxed mode:
  # Balanced: allocate ensuring that any two taxa differ by at most 1 if
  #           m_i < n_i (or m_i < childrenCount), but allow more if one taxon
  #           has m_i == n_i.
  # Randomized: if we don't have the m fully distributed over m_i (children),
  #             choose a random child from childrenCount that still has taxa
  #             available to choose.
  # Strict mode:
  # Balanced: ensure two taxa allocation (m_i) differ at most by 1.
  # Randomized: fully random child sampling, uniform distribution among taxa.
  if (method == "balanced" & strict == "no") {
    while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m) {
      child <- names(childrenCount[childrenCount > m_i])
      m_i[child] <- m_i[child] + 1
      m <- m - length(child)
    }
    child <- sample(names(childrenCount[childrenCount > m_i]), m)
    m_i[child] <- m_i[child] + 1
  } else if (method == "randomized" & strict == "no") {
    while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
      child <- sample(names(childrenCount[childrenCount > m_i]), 1)
      m_i[child] <- m_i[child] + 1
      m <- m - 1
    }
  } else if (method == "balanced" & strict == "yes") {
    m_i <- m_i + floor(m / length(names(childrenCount)))
    sampledChildren <- sample(names(childrenCount), m - sum(m_i))
    m_i[sampledChildren] <- m_i[sampledChildren] + 1
  } else if (method == "randomized" & strict == "yes") {
    m_i <- table(sample(names(childrenCount), m, replace = TRUE))
  }
  
  
  # Require specific IDs (and their parents/ancestors) to be present.
  # If any required ID has a lower m_i allocated than needed (informed by
  # requireIDs), reallocate from another ID in m_i.
  if (!is.null(requireIDs)) {
    req <- is.element(names(m_i), names(requireIDs))
    toAdd <- m_i[req] < requireIDs[names(m_i[req])]
    if (method == "balanced") {
      while (any(toAdd)) {
        subtract <- m_i[setdiff(names(m_i), names(toAdd))]
        subtract <- sample(names(subtract[subtract == max(subtract)]), 1)
        m_i[subtract] <- m_i[subtract] - 1
        add <- sample(names(toAdd[toAdd == TRUE]), 1)
        m_i[add] <- m_i[add] + 1
        toAdd <- m_i[req] < requireIDs[names(m_i[req])]
      }
    } else if (method == "randomized") {
      while (any(toAdd)) {
        subtract <- m_i[setdiff(names(m_i), names(toAdd))]
        subtract <- sample(names(subtract[subtract > 0]), 1)
        m_i[subtract] <- m_i[subtract] - 1
        add <- sample(names(toAdd[toAdd == TRUE]), 1)
        m_i[add] <- m_i[add] + 1
        toAdd <- m_i[req] < requireIDs[names(m_i[req])]
      }
    }
  }
  
  outputIDs <- character(0)
  for (id in names(m_i)) {
    if (m_i[id] == 0) {
      next
    } else if (id == as.character(taxon)) {
      outputIDs <- c(outputIDs, id)
    } else {
      outputIDs <- c(outputIDs,
                     TS_Algorithm(id, m_i[id], nodes, countIDs, replacement,
                                  method, strict, requireIDs))
    }
  }
  
  return(outputIDs)
}


Wrapper_TS_Algorithm <- function(taxon, m, nodes, countIDs, replacement = "no",
                                 method = "balanced", strict = "no", 
                                 ignoreIDs = NULL, requireIDs = NULL,
                                 ignoreNonLeafID = NULL) {
  # Wrapper for the TS algorithm. Provides support for performance, for
  # ignoring specific taxon IDs, and for requiring specific taxon IDs.
  #
  # Args:
  #   taxon: (char/integer) Taxon from which to start sampling children taxa.
  #   m: (integer) size of the sample to generate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure. Created by getnodes() from the CHNOSZ
  #                       package.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   replacement: (char) whether the algorithm allows to repeat IDs in order
  #                       to maximize taxonomy diversity and to reach m IDs
  #                       in the output.
  #   method: (char) whether the algorithm will choose IDs randomly or 
  #                  maintaining a balanced allocation (m_i differing by no
  #                  more than 1 if the maximum possible value wasn't reached).
  #   strict: (char) whether it favors balanced taxa representation 
  #                  (strict = "yes") or maximized taxa representation ("no").
  #   ignoreIDs: (char) IDs that mustn't appear in the output.
  #   requireIDs: (char) IDs that must appear in the output.
  #   ignoreNonLeafID: (char) (testing) a non-leaf ID to ignore; won't apply
  #                           to leaf nodes and won't exclude its children from
  #                           the analysis, unlike ignoreIDs.
  # Returns:
  #   outputIDs: (vector) vector of IDs with maximized taxonomy diversity.
  
  # Reduce the node information to he necessary only, reduces search time.
  nodes <- nodes[is.element(nodes$id, names(countIDs)), 1:2]
  
  
  if (!is.null(ignoreIDs)) {
    ignoreIDs <- as.integer(ignoreIDs)
    
    # Sanity check: filter IDs that aren't part of NCBI notation.
    if (!all(is.element(ignoreIDs, nodes$id))) {
      cat("Warning: the following inputs are not part of NCBI taxonomy IDs",
          "and will be ignored.\n",
          ignoreIDs[!is.element(ignoreIDs, nodes$id)], "\n")
      ignoreIDs <- ignoreIDs[is.element(ignoreIDs, nodes$id)]
    }
    
    # Sanity check: unique ID inputs, remove duplicates.
    if (any(duplicated(ignoreIDs))) {
      cat("Warning: some required IDs are repeated,",
          "using only one instance of each.\n",
          ignoreIDs[duplicated(ignoreIDs)], "\n")
      ignoreIDs <- unique(ignoreIDs)
    }
    
    # Make sure none of the ignoreIDs is a parent/ancestor of another.
    # It eases the processing of the following steps.
    for (id in ignoreIDs) {
      if (any(is.element(ignoreIDs[ignoreIDs != id],
                         allparents(id, nodes = nodes)))) {
        ignoreIDs <- ignoreIDs[ignoreIDs != id]
      }
    }
    
    # Subtract the presence of the ignoreIDs from the countIDs.
    for (id in ignoreIDs) {
      parentIDs <- as.character(allparents(id, nodes = nodes))
      countIDs[parentIDs] <- countIDs[parentIDs] - countIDs[as.character(id)]
    }
    countIDs <- countIDs[countIDs > 0]
    
    # Pruning the children of removed nodes. 
    orphans <- nodes$id[!is.element(nodes$parent, names(countIDs)) &
                         is.element(nodes$id, names(countIDs))]
    while (length(orphans) > 0) {
      countIDs <- countIDs[setdiff(names(countIDs), orphans)]
      orphans <- nodes$id[!is.element(nodes$parent, names(countIDs)) &
                           is.element(nodes$id, names(countIDs))]
    }
  }
  
  if (!is.null(ignoreNonLeafID)) {
    ignoreNonLeafID <- as.integer(ignoreNonLeafID)
    
    # Sanity check: filter IDs that aren't part of NCBI notation.
    if (!all(is.element(ignoreNonLeafID, nodes$id))) {
      cat("Warning: the following inputs are not part of NCBI taxonomy IDs",
          "and will be ignored.\n",
          ignoreNonLeafID[!is.element(ignoreNonLeafID, nodes$id)], "\n")
      ignoreNonLeafID <- ignoreNonLeafID[is.element(ignoreNonLeafID, nodes$id)]
    }
    
    # Sanity check: unique ID inputs, remove duplicates.
    if (any(duplicated(ignoreNonLeafID))) {
      cat("Warning: some required IDs are repeated,",
          "using only one instance of each.\n",
          ignoreNonLeafID[duplicated(ignoreNonLeafID)], "\n")
      ignoreNonLeafID <- unique(ignoreNonLeafID)
    }
    
    # Make sure every ignoreNonLeafID is not a leaf node
    # (i.e. not just a consequence of its children).
    for (id in ignoreNonLeafID) {
      children <- nodes$id[nodes$parent == id & nodes$id != id]
      children <- intersect(children, names(countIDs))
      childrenSum <- sum(countIDs[as.character(children)])
      if (childrenSum > 0 & childrenSum < countIDs[as.character(id)]) {
        countIDs[as.character(id)] <- childrenSum
      }
    }
  }
  
  if (!is.null(requireIDs)) {
    requireIDs <- as.integer(requireIDs)
    
    # Sanity check: only IDs that are present in the input.
    if (!all(is.element(requireIDs, nodes$id))) {
      cat("Warning: the following required IDs are not part of your input",
          "and will be ignored.\n",
          requireIDs[!is.element(requireIDs, nodes$id)], "\n")
      requireIDs <- requireIDs[is.element(requireIDs, nodes$id)]
    }
    
    # Sanity check: remove IDs that has any ignoreIDs as its parent.
    if (!all(is.element(requireIDs, names(countIDs)))) {
      cat("Warning: the following required IDs are children or part of your",
          "ignored IDs and will be ignored.\n",
          requireIDs[!is.element(requireIDs, names(countIDs))], "\n")
      requireIDs <- requireIDs[is.element(requireIDs, names(countIDs))]
    }
    
    # Sanity check: unique ID inputs, remove duplicates.
    if (any(duplicated(requireIDs))) {
      cat("Warning: some required IDs are repeated,",
          "using only one instance of each.\n",
          requireIDs[duplicated(requireIDs)], "\n")
      requireIDs <- unique(requireIDs)
    }
    
    if (length(requireIDs) > 0) {
      requireIDs <- TS_TaxonomyData(requireIDs, nodes)
    } else {
      requireIDs <- NULL
    }
  }
  
  # Reduce, once again, the node information to he necessary only,
  # reduces search time.
  nodes <- nodes[is.element(nodes$id, names(countIDs)), 1:2]
  

  # Ensure m <= number of valid ids.
  if (!is.null(idsFile) & isTRUE(file.exists(as.character(idsFile)))) {
    ids <- read.table(idsFile, sep = "\t", header = FALSE, comment.char = "",
                      colClasses = "integer", strip.white = TRUE)
    ids <- ids[, 1]
  } else {
    ids <- as.integer(idsFile)  # assume it's a test.name or back.name, for now
  }
  if (m > length(intersect(ids, names(countIDs)))) {
    m <- length(intersect(ids, names(countIDs)))
  }



  # Call the TS algorithm itself.
  outputIDs <- TS_Algorithm(taxon, m, nodes, countIDs, replacement,
                            method, strict, requireIDs)

  return(outputIDs)
}


RandomSampling <- function(idsFile, nodes, n = 100) {
  # Checks valid taxonomy IDs and samples them under an uniform distribution.
  #
  # Args:
  #   idsFile: (char/int) either a path to a file with the input taxon IDs, or
  #                       a vector with the IDs themselves, in either integer
  #                       or character format.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure. Created by getnodes() from the CHNOSZ
  #                       package.
  #   n: (int) number of IDs to sample.
  # Returns:
  #   randomIDs: (int) sampled IDs under an uniform distribution.

  # Get ids to search
  if (!is.null(idsFile) & isTRUE(file.exists(as.character(idsFile)))) {
    ids <- read.table(idsFile, sep = "\t", header = FALSE, comment.char = "",
                      colClasses = "integer", strip.white = TRUE)
    ids <- ids[, 1]
  } else {
    ids <- as.integer(idsFile)  # assume it's a test.name or back.name, for now
  }
  
  # Sanity check: filter IDs that aren't part of NCBI notation.
  if (!all(is.element(ids, nodes$id))) {
    cat("Warning: the following inputs are not part of NCBI taxonomy IDs",
        "and will be ignored.\n",
        ids[!is.element(ids, nodes$id)], "\n")
    ids <- ids[is.element(ids, nodes$id)]
  }
  
  # Sanity check: unique ID inputs, remove duplicates.
  if (any(duplicated(ids))) {
    cat("Warning: some ids are repeated, using only one instance of each.\n",
        ids[duplicated(ids)], "\n")
    ids <- unique(ids)
  }

  # Random sampling
  randomIDs <- sample(ids, n)
  return(randomIDs)

}



Evaluate_TS <- function(outputIDs, nodes, countIDs) {
  # Collect metrics about the taxon diversity of a sample (outputIDs).
  #
  # Args:
  #   outputIDs: (vector) vector of IDs to evaluate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure. Created by getnodes() from the CHNOSZ
  #                       package.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   level: (int) which rank level (from NCBI's taxonomy structure) to
  #                evaluate outoutIDs.
  # Returns:
  #   taxonDiversity: (int) proportion of sampled taxa compared to how many 
  #                         were possible in each taxonomy level.  
  
  selectedIDs <- TS_TaxonomyData(outputIDs, nodes) 
#  taxonDiversity <- numeric(0)
  listTaxon <- list()
  
  
  # Find the sub-taxa (children nodes) of the current taxon
  children <- 1
  while (length(children) > 0) {
    taxon <- as.integer(children)
    children <- nodes$id[is.element(nodes$parent, taxon) & 
                         !is.element(nodes$id, taxon)]

    children <- intersect(children, names(countIDs))
    selectedChildren <- intersect(children, names(selectedIDs))

    listTaxon <- c(listTaxon, list(selectedChildren))

#    # Metric: representation/diversity
#    taxonDiversity <- c(taxonDiversity,
#                        length(selectedChildren)/min(length(outputIDs),
#                                                     length(children)))
  }
  
  
  # Metric: balance
  #selectedIDs[children]

#  return(taxonDiversity)
  return(listTaxon) 
}



