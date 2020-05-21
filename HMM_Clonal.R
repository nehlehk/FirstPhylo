library("HMM")
library(dplyr)
library(tidyr)

data = scan("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/12.txt", character(), quote = "")
data = noquote(data)
data = replace_na(data,0)
data = as.numeric(data)
n = length(data)
avg = mean(data)
std = sd(data)




States = c("Clonal", "non_Clonal")
Symbols = 0:9
transProbs = matrix(c(0.9, 0.10, 0.10, 0.90), c(length(States),length(States)), byrow = TRUE)
startProbs = c(0.97 , 0.03)
emissionProbs = matrix(c(rnorm(n = n , mean = avg , sd = std), rnorm(n = n , mean = avg , sd = std)) ,nrow = 2, byrow = TRUE)

hmm = initHMM(States, Symbols, startProbs =startProbs , transProbs = transProbs, emissionProbs = emissionProbs)

vit = viterbi(hmm, data)