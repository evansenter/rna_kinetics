prob.at.time <- function(eigen.rate.data, structure.index, time) {
  Reduce(function(sum.j, j) {
    # This code is hard-wired to only consider the kinetics for folding from the empty structure, to save an order of complexity.
    sum.j + (
      eigen.rate.data$vectors[structure.index, j] * 
      eigen.rate.data$inverse.vectors[j, eigen.rate.data$key.indices$empty] * 
      exp(eigen.rate.data$values[j] * time)
    )
  }, 1:ncol(eigen.rate.data$vectors), 0)
}