argv <- commandArgs(TRUE)

if (length(argv) != 2) {
  cat("./Rscript serializer.r INPUT_STRUCTURES OUTPUT_FILENAME\n")
  q("no")
}

if (is.na(file.info(argv[1])$size)) {
  cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
  q("no")
}

if (!is.na(file.info(argv[2])$size)) {
  cat(paste(argv[2], "is a file that already exists, cowering out.\n"))
  q("no")
}

dyn.load("rnaeval_r.so")

RT             <- 0.0019872370936902486 * (37 + 273.15) # RT in kcal / mol at 37C
structures     <- read.csv(argv[1])
structures$mfe <- apply(
  structures, 
  1, 
  function(row) .C("energy_of_struct_r", sequence = colnames(structures), structure = row[1], energy = as.double(1))$energy
)
energies           <- structures[, 2]
key.indices        <- as.list(match(c(0, min(energies)), energies))
names(key.indices) <- c("empty", "mfe")

col.rate.matrix <- matrix(nrow = nrow(structures), ncol = nrow(structures), byrow = TRUE)

for(i in 1:nrow(col.rate.matrix)) {
  for(j in 1:ncol(col.rate.matrix)) {
    if (i != j) {
      col.rate.matrix[j, i] <- min(1, exp(-(energies[j] - energies[i]) / RT))
    }
  }

  col.rate.matrix[i, i] <- -sum(Filter(function(x) !is.na(x), col.rate.matrix[, i]))
}

eigen.rate.data                 <- eigen(col.rate.matrix)
eigen.rate.data$inverse.vectors <- solve(eigen.rate.data$vectors)
eigen.rate.data$key.indices     <- key.indices

saveRDS(eigen.rate.data, argv[2])