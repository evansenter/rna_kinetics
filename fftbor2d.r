# argv <- commandArgs(TRUE)

# if (length(argv) != 2) {
#   cat("./Rscript serializer.r INPUT_STRUCTURES OUTPUT_FILENAME\n")
#   q("no")
# }

# if (is.na(file.info(argv[1])$size)) {
#   cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
#   q("no")
# }

# if (!is.na(file.info(argv[2])$size)) {
#   cat(paste(argv[2], "is a file that already exists, cowering out.\n"))
#   q("no")
# }

# dyn.load("rnaeval_r.so")

# RT             <- 0.0019872370936902486 * (37 + 273.15) # RT in kcal / mol at 37C
# structures     <- read.csv(argv[1])
# structures$mfe <- apply(
#   structures, 
#   1, 
#   function(row) .C("energy_of_struct_r", sequence = colnames(structures), structure = row[1], energy = as.double(1))$energy
# )
# energies           <- structures[, 2]
# key.indices        <- as.list(match(c(0, min(energies)), energies))
# names(key.indices) <- c("empty", "mfe")

# col.rate.matrix <- matrix(nrow = nrow(structures), ncol = nrow(structures), byrow = TRUE)

# for(i in 1:nrow(col.rate.matrix)) {
#   for(j in 1:ncol(col.rate.matrix)) {
#     if (i != j) {
#       col.rate.matrix[j, i] <- min(1, exp(-(energies[j] - energies[i]) / RT))
#     }
#   }

#   col.rate.matrix[i, i] <- -sum(Filter(function(x) !is.na(x), col.rate.matrix[, i]))
# }

# eigen.rate.data                 <- eigen(col.rate.matrix)
# eigen.rate.data$inverse.vectors <- solve(eigen.rate.data$vectors)
# eigen.rate.data$key.indices     <- key.indices

# saveRDS(eigen.rate.data, argv[2])

################################################################################################################################ 
# Can make a wrapper for calling FFTbor2D directly from R using .C
# The (main) problem right now seems to be in the encoding

library(Matrix)
fftbor.data   <- "test_6.txt"
seq.length    <- as.numeric(system(paste0("ruby -r vienna_rna -e 'p(RNA.from_fasta(\"", fftbor.data, "\").seq.length)'"), intern = T))
xition.ncol   <- seq.length + 1
fftbor.file   <- tempfile()
write(system(paste("FFTbor2D -S", fftbor.data), intern = T), fftbor.file)
fftbor.data   <- read.delim(fftbor.file, col.names = c("i", "j", "p", "ensemble"))
fftbor.data   <- within(fftbor.data, ij <- i * xition.ncol + j)
shift.moves   <- function(x) {
  x + as.numeric(lapply(Map(function(x) x, Map(function(i) { (xition.ncol + 1i) * 1i ^ i }, 0:3)), function(x) Re(x) + Im(x)))
}
valid.moves   <- function(x) { Filter(function(y) {
  sum((c(floor(y / xition.ncol), y %% xition.ncol) - c(floor(x / xition.ncol), x %% xition.ncol)) ^ 2) == 2
}, shift.moves(x)) }

transition.list            <- as.data.frame(do.call(rbind, apply(fftbor.data[c("ij", "p")], 1, function(row) do.call(rbind, lapply(c(valid.moves(row[1]), row[1]), function(move) c(row[1], move))))))
row.names(transition.list) <- NULL

transition.move.prob <- function(from, to) {
  ifelse(
    nrow(fftbor.data[fftbor.data$ij == to,]) == 0, 
    0, 
    min(1, fftbor.data[fftbor.data$ij == to,]$p / fftbor.data[fftbor.data$ij == from,]$p) / length(valid.moves(from))
  )
}
transition.stay.prob <- function(from) {
  1 - Reduce(function(sum, x) sum + transition.move.prob(from, x), valid.moves(from), 0)
}

transition.list <- cbind(transition.list, apply(transition.list, 1, function(row) ifelse(row[1] == row[2], transition.stay.prob(row[1]), transition.move.prob(row[1], row[2]))))
names(transition.list) <- c("from", "to", "p")

if (nrow(transition.list[transition.list$from == xition.ncol ^ 2 - 1 & transition.list$to == xition.ncol ^ 2 - 1,]) == 0) {
  transition.list <- rbind(transition.list, c(xition.ncol ^ 2 - 1, xition.ncol ^ 2 - 1, 0))
}

transition.matrix <- sparseMatrix(i = transition.list$from, j = transition.list$to, x = transition.list$p, index1 = F)

# > example.matrix <- matrix(sapply(runif(100), function(x) ifelse(x < 0.8, 1e-10 * rnorm(1), x)), nrow = 10)
# > summary(example.matrix)
#        V1                V2               V3               V4        
#  Min.   :0.00000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#  1st Qu.:0.00000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#  Median :0.00000   Median :0.0000   Median :0.0000   Median :0.0000  
#  Mean   :0.09437   Mean   :0.3668   Mean   :0.2740   Mean   :0.2665  
#  3rd Qu.:0.00000   3rd Qu.:0.8693   3rd Qu.:0.6576   3rd Qu.:0.6164  
#  Max.   :0.94369   Max.   :0.9970   Max.   :0.9327   Max.   :0.9925  
#        V5               V6                   V7                V8        
#  Min.   :0.0000   Min.   :-8.578e-11   Min.   :0.00000   Min.   :0.0000  
#  1st Qu.:0.0000   1st Qu.:-2.125e-11   1st Qu.:0.00000   1st Qu.:0.0000  
#  Median :0.0000   Median : 1.166e-11   Median :0.00000   Median :0.0000  
#  Mean   :0.1802   Mean   : 9.039e-12   Mean   :0.09636   Mean   :0.1696  
#  3rd Qu.:0.0000   3rd Qu.: 3.765e-11   3rd Qu.:0.00000   3rd Qu.:0.0000  
#  Max.   :0.9259   Max.   : 1.012e-10   Max.   :0.96358   Max.   :0.8483  
#        V9              V10         
#  Min.   :0.0000   Min.   :0.00000  
#  1st Qu.:0.0000   1st Qu.:0.00000  
#  Median :0.0000   Median :0.00000  
#  Mean   :0.3467   Mean   :0.08685  
#  3rd Qu.:0.8384   3rd Qu.:0.00000  
#  Max.   :0.9492   Max.   :0.86849  
# > solve(example.matrix)
#                [,1]          [,2]          [,3]          [,4]          [,5]
#  [1,]  4.823554e-11  4.163170e-01 -1.859031e-10 -3.755312e+00  1.059669e+00
#  [2,]  1.417369e-10  1.611454e-01 -9.639641e-11 -1.017220e+00  1.145411e-10
#  [3,] -1.593782e-10 -1.211507e-01  1.760503e-10  1.905317e+00 -1.314745e-10
#  [4,] -1.927770e-10 -3.344075e-01 -1.177074e-10  2.110926e+00 -2.265881e-10
#  [5,]  1.407325e-10  5.915764e-05 -8.856924e-11 -3.734304e-04  7.900745e-13
#  [6,] -1.328237e+00 -3.869238e+09 -1.325439e+00  2.442432e+10 -1.753286e+00
#  [7,] -2.049295e-10 -2.560033e-01  1.037795e+00  1.616005e+00 -9.750630e-11
#  [8,] -1.277228e-10 -8.582732e-01  2.203684e-10 -2.758019e+00  4.225506e-10
#  [9,]  2.113216e-10  1.299855e+00  8.277366e-11 -8.985461e-01 -4.071271e-11
# [10,]  1.151419e+00 -1.398628e+00 -9.679265e-11  1.676792e+00 -1.847376e-10
#                [,6]          [,7]          [,8]          [,9]         [,10]
#  [1,] -2.473826e+00  2.475546e+00 -3.291821e-10  3.629006e+00 -7.392144e-02
#  [2,] -9.575531e-01  9.582187e-01 -2.255642e-10  2.070024e+00 -2.861304e-02
#  [3,]  7.198981e-01 -7.203985e-01  4.339562e-10 -1.774911e+00  2.151157e-02
#  [4,]  1.987106e+00 -1.988487e+00  3.292711e-10 -2.113845e+00  1.066920e+00
#  [5,] -3.515262e-04  3.517706e-04  1.141616e+00  3.739470e-04 -1.050399e-05
#  [6,]  2.299167e+10 -2.300765e+10  3.990532e+00 -2.445810e+10  6.870237e+08
#  [7,]  1.521215e+00 -1.522272e+00  2.831500e-10 -1.618239e+00  4.545606e-02
#  [8,] -1.417469e+00  2.598047e+00 -1.930054e-10  1.526173e+00 -1.053696e+00
#  [9,] -8.458402e-01  8.464282e-01 -1.830985e-10  8.997885e-01 -2.527490e-02
# [10,]  1.578437e+00 -1.579534e+00 -1.217021e+00 -1.679110e+00  4.716592e-02