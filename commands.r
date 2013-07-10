dyn.load("rnaeval_r.so")

RT             <- 0.0019872370936902486 * (37 + 273.15) # RT in kcal / mol at 37C
structures     <- read.csv("test_10nt_all_str.txt")
structures$mfe <- apply(
  structures, 
  1, 
  function(row) .C("energy_of_struct_r", sequence = colnames(structures), structure = row[1], energy = as.double(1))$energy
)
energies           <- structures[, 2]
key.indices        <- as.list(match(c(0, min(energies)), energies))
names(key.indices) <- c("empty", "mfe")

row.rate.matrix <- matrix(nrow = nrow(structures), ncol = nrow(structures), byrow = TRUE)

for(i in 1:nrow(row.rate.matrix)) {
  for(j in 1:ncol(row.rate.matrix)) {
    if (i != j) {
      row.rate.matrix[i, j] <- min(1, exp(-(energies[j] - energies[i]) / RT))
    }
  }

  row.rate.matrix[i, i] <- -sum(Filter(function(x) !is.na(x), row.rate.matrix[i,]))
}

col.rate.matrix                 <- t(row.rate.matrix)
eigen.rate.data                 <- eigen(col.rate.matrix)
eigen.rate.data$inverse.vectors <- solve(eigen.rate.data$vectors)


prob.at.time <- function(eigen.rate.data, structure.index, time, key.indices) {
  Reduce(function(sum.j, j) {
    sum.j + (Reduce(function(sum.k, k) {
      # This ifelse statement can get pulled out as a pretty significant optimization 
      #   since it will only be non-zero when key.indices$empty == k
      sum.k + (
        eigen.rate.data$vectors[structure.index, j] * 
        eigen.rate.data$inverse.vectors[j, k]       * 
        ifelse(key.indices$empty == k, 1, 0)
      )
    }, 1:ncol(eigen.rate.data$vectors), 0) * exp(eigen.rate.data$values[j] * time))
  }, 1:ncol(eigen.rate.data$vectors), 0)
}

prob.at.time(eigen.rate.data, key.indices$empty, 0, key.indices)