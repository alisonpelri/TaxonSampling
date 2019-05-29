

n <- 100
x <- c(50, 100, 150, 200, 250, 300, 350, 400)

source("TaxonSampling.R")
library("ggplot2")

# For parallel processing. for a serial run, do "cores <- 1"
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
cores <- 1
if (cores > 1) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}


#comb <- function(x, ...) {
#  lapply(seq_along(x),
#    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
#}
#
#oper <- foreach(i=1:10, .combine='comb', .multicombine=TRUE,
#                .init=list(list(), list())) %dopar% {
#  list(i+2, i+3)
#}


#foreach (i = 1:n, .export = c("Wrapper_TS_Algorithm", "RandomSampling",
#                              "Evaluate_TS", "TS_TaxonomyData",
#                              "TS_Algorithm")) %dopar% {


outputDiversity <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                        "4" = numeric(0), "5" = numeric(0))
for (i in 1:n) {
  listOutput <- list()
  listOutput$"50" <- Wrapper_TS_Algorithm(1, 50, nodes, countIDs)
  listOutput$"100" <- Wrapper_TS_Algorithm(1, 100, nodes, countIDs)
  listOutput$"150" <- Wrapper_TS_Algorithm(1, 150, nodes, countIDs)
  listOutput$"200" <- Wrapper_TS_Algorithm(1, 200, nodes, countIDs)
  listOutput$"250" <- Wrapper_TS_Algorithm(1, 250, nodes, countIDs)
  listOutput$"300" <- Wrapper_TS_Algorithm(1, 300, nodes, countIDs)
  listOutput$"350" <- Wrapper_TS_Algorithm(1, 350, nodes, countIDs)
  listOutput$"400" <- Wrapper_TS_Algorithm(1, 400, nodes, countIDs)
  
  evalOutput <- list()
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  evalOutput$"150" <- Evaluate_TS(listOutput$"150", nodes, countIDs)
  evalOutput$"200" <- Evaluate_TS(listOutput$"200", nodes, countIDs)
  evalOutput$"250" <- Evaluate_TS(listOutput$"250", nodes, countIDs)
  evalOutput$"300" <- Evaluate_TS(listOutput$"300", nodes, countIDs)
  evalOutput$"350" <- Evaluate_TS(listOutput$"350", nodes, countIDs)
  evalOutput$"400" <- Evaluate_TS(listOutput$"400", nodes, countIDs)
  
  
  for (level in 3:5) {
    diversity <- numeric(0)
    for (number in names(evalOutput)) {
      diversity <- c(diversity, length(evalOutput[[number]][[level]]))
    }
    outputDiversity[[level]] <- rbind(outputDiversity[[level]], diversity)
  }
}


randomDiversity <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                        "4" = numeric(0), "5" = numeric(0))
for (i in 1:n) {
  listRandom <- list()
  listRandom$"50" <- RandomSampling(idsFile, nodes, 50)
  listRandom$"100" <- RandomSampling(idsFile, nodes, 100)
  listRandom$"150" <- RandomSampling(idsFile, nodes, 150)
  listRandom$"200" <- RandomSampling(idsFile, nodes, 200)
  listRandom$"250" <- RandomSampling(idsFile, nodes, 250)
  listRandom$"300" <- RandomSampling(idsFile, nodes, 300)
  listRandom$"350" <- RandomSampling(idsFile, nodes, 350)
  listRandom$"400" <- RandomSampling(idsFile, nodes, 400)
  
  evalRandom <- list()
  evalRandom$"50" <- Evaluate_TS(listRandom$"50", nodes, countIDs)
  evalRandom$"100" <- Evaluate_TS(listRandom$"100", nodes, countIDs)
  evalRandom$"150" <- Evaluate_TS(listRandom$"150", nodes, countIDs)
  evalRandom$"200" <- Evaluate_TS(listRandom$"200", nodes, countIDs)
  evalRandom$"250" <- Evaluate_TS(listRandom$"250", nodes, countIDs)
  evalRandom$"300" <- Evaluate_TS(listRandom$"300", nodes, countIDs)
  evalRandom$"350" <- Evaluate_TS(listRandom$"350", nodes, countIDs)
  evalRandom$"400" <- Evaluate_TS(listRandom$"400", nodes, countIDs)

  for (level in 3:5) {
    diversity <- numeric(0)
    for (number in names(evalRandom)) {
      diversity <- c(diversity, length(evalRandom[[number]][[level]]))
    }
    randomDiversity[[level]] <- rbind(randomDiversity[[level]], diversity)
  }
}


save.image("Eval.RData")

confidence <- .995  # 99% = .995, 95% = .975
outputMeans <- list()
randomMeans <- list()
outputCI <- list()
randomCI <- list()
for (level in 3:5) {
  outputMeans[[level]] <- colMeans(outputDiversity[[level]])  
  randomMeans[[level]] <- colMeans(randomDiversity[[level]])
  outputCI[[level]] <- apply(outputDiversity[[level]], 2, sd)
  randomCI[[level]] <- apply(randomDiversity[[level]], 2, sd)
  outputCI[[level]] <- outputSD[[level]]/sqrt(n)
  randomCI[[level]] <- randomSD[[level]]/sqrt(n)
  outputCI[[level]] <- qt(confidence, df = n-1) * outputSD[[level]]
  randomCI[[level]] <- qt(confidence, df = n-1) * randomSD[[level]]
}


for (level in 3:5) {
  df <- data.frame(x, outputMeans = outputMeans[[level]], 
                      randomMeans = randomMeans[[level]])
  
  imageName <- paste0("level", level, ".png")
  png(imageName)
  print(ggplot(df, aes(x)) + 
          geom_point(aes(y=outputMeans, colour="TS")) +
          geom_line(aes(y=outputMeans, colour="TS")) + 
          geom_errorbar(aes(ymin=outputMeans - outputCI[[level]],
                            ymax=outputMeans + outputCI[[level]],
                            colour = "TS"), width=1) +
          geom_point(aes(y=randomMeans, colour="RS")) +
          geom_line(aes(y=randomMeans, colour="RS")) +
          geom_errorbar(aes(ymin=randomMeans - randomCI[[level]],
                            ymax=randomMeans + randomCI[[level]],
                            colour = "RS"), width=1) +
          xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
          scale_color_manual("Method",
                             values = c("TS" = "red", "RS" = "blue")))
  dev.off()
}


