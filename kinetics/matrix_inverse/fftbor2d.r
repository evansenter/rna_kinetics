# Finish implementing the triangle inequality

argv <- commandArgs(TRUE)

if (!is.element(length(argv), 1:2)) {
  cat("./Rscript fftbor2d.r INPUT_FA_FILE [INPUT_FFTBOR2D_RUN]\n")
  q("no")
}

if (is.na(file.info(argv[1])$size)) {
  cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
  q("no")
}

if (length(argv) == 2 & is.na(file.info(argv[2])$size)) {
  cat(paste(argv[2], "doesn't appear to exist, cowering out.\n"))
  q("no")
}

library(Matrix)
library(gtools)
require(corpcor)
library(rbenchmark)

unpruned.transition.matrix <- function(fa.input, fftbor2d.output) {
  seq.length    <- as.numeric(system(paste0("ruby -r vienna_rna -e 'p(RNA.from_fasta(\"", fa.input, "\").seq.length)'"), intern = T))
  xition.ncol   <- seq.length + 1
  bp.dist       <- as.numeric(system(paste0(
    "ruby -r vienna_rna -e 'p ViennaRna::Rna.bp_distance(*File.read(\"",
    fa.input,
    "\").chomp.split(/\\n/)[-2..-1])'"), intern = T
  ))

  if (missing(fftbor2d.output)) {
    fftbor2d.output <- tempfile()
    write(system(paste("FFTbor2D -S -E /Users/evansenter/Source/fftbor/rna_turner_2.1.2.par", fa.input), intern = T), fftbor2d.output)
  }

  fftbor.data   <- read.delim(fftbor2d.output, header = F, col.names = c("i", "j", "p", "ensemble"))
  if (nrow(fftbor.data[fftbor.data$i == bp.dist & fftbor.data$j == 0,]) == 0) {
    fftbor.data <- rbind(fftbor.data, c(bp.dist, 0, 0, 0))
  }
  
  two.d.to.rmoi <- function(i, j) {
    i * xition.ncol + j
  }
  
  rmoi.to.two.d <- function(ij) {
    c(floor(ij / xition.ncol), ij %% xition.ncol)
  }

  fftbor.data <- within(fftbor.data, ij <- two.d.to.rmoi(i, j))
  fftbor.data <- fftbor.data[order(fftbor.data$ij),]
  

  fftbor.matrix <- sparseMatrix(
    i      = fftbor.data$i, 
    j      = fftbor.data$j, 
    x      = fftbor.data$p, 
    index1 = F
  )

  shift.moves <- function(x) {
    x + as.numeric(lapply(Map(function(x) x, Map(function(i) { (xition.ncol + 1i) * 1i ^ i }, 0:3)), function(x) Re(x) + Im(x)))
  }
  valid.moves <- function(x) { Filter(function(y) {
    y >= 0 & 
    sum(rmoi.to.two.d(y)) >= bp.dist &
    sum((rmoi.to.two.d(y) - rmoi.to.two.d(x)) ^ 2) == 2
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

  transition.list  <- transition.list[transition.list$p > 0,]
  uniq.transitions <- unique(c(transition.list$from, transition.list$to))
  mapping          <- cbind(uniq.transitions, order(uniq.transitions))
  index.of         <- function(unconsolidated) { mapping[mapping[,1] == unconsolidated][2] }
  unindex.of       <- function(consolidated) { mapping[mapping[,2] == consolidated][1] }

  matrix.data <- {}

  # MAPPING RETURNS A 1-INDEXED LIST OF MAPPINGS FROM THE 0-INDEXED TRANSITION MATRIX
  matrix.data$unpruned.transition.matrix <- sparseMatrix(
    i      = sapply(transition.list$from, index.of), 
    j      = sapply(transition.list$to, index.of), 
    x      = transition.list$p
  )

  matrix.data$index.mapper   <- index.of
  matrix.data$index.unmapper <- unindex.of
  matrix.data$bp.dist        <- bp.dist
  matrix.data$xition.ncol    <- xition.ncol
  matrix.data$fftbor.matrix  <- fftbor.matrix

  matrix.data
}

transition.matrix.for.inversion <- function(matrix.data) {
  pruned.matrix <- matrix.data$unpruned.transition.matrix[
    -matrix.data$index.mapper(matrix.data$xition.ncol * matrix.data$bp.dist), 
    -matrix.data$index.mapper(matrix.data$xition.ncol * matrix.data$bp.dist) 
  ]

  matrix.data$inversion.matrix <- diag(nrow(pruned.matrix)) - pruned.matrix

  matrix.data
}

visualize.transition.matrix.from.fa.file <- function(fa.input) {
  matrix.data <- transition.matrix.for.inversion(unpruned.transition.matrix(fa.input))

  quartz("Heatmap", 6, 6)
  image(
    x    = 1:dim(matrix.data)[[1]], 
    y    = 1:dim(matrix.data)[[2]], 
    z    = as.matrix(matrix.data), 
    col  = rev(gray(0:64 / 64)),
    xlab = "Remapped column index",
    ylab = "Remapped row index"
  )
  title(paste("Transition matrix for", basename(fa.input)))
}

visualize.fftbor2d.matrix.from.fa.file <- function(fa.input) {
  matrix.data <- unpruned.transition.matrix(fa.input)$fftbor.matrix

  quartz("Heatmap", 6, 6)
  image(
    x    = 1:dim(matrix.data)[[1]], 
    y    = 1:dim(matrix.data)[[2]], 
    z    = as.matrix(matrix.data),
    col  = rev(gray(0:64 / 64)),
    xlab = "Column index",
    ylab = "Row index"
  )
  title(paste("FFTbor2D matrix for", basename(fa.input)))
}

mfpt.from.fa.using.fftbor2d <- function(fa.input, fftbor2d.output) {
  # This check happens by definition of transition.stay.prob
  # if (!all(sapply(apply(transition.matrix, 1, sum), function(x) { x %in% 0:1 }))) {
  #   warning("One or more row-sums are not 0 or 1:")
  #   print(apply(transition.matrix, 1, sum)[apply(transition.matrix, 1, sum) > 0])
  # }
  # 
  # Hurts the runtime so I'm removing this check.
  # rev.utriangle.matrix <- upperTriangle(transition.matrix[, rev(seq_len(ncol(transition.matrix)))], diag = T)
  # if (nrow(transition.list[transition.list$p > 0,]) != length(subset(rev.utriangle.matrix, rev.utriangle.matrix > 0))) {
  #   warning("The y-mirrored transition matrix is not upper triangular (only values in the upper-left quadrant are permitted)")
  # }
  # 
  # This check needs to toggle odd/even parity on the bp.dist parity.
  # if (
  #   length(subset(transition.list$from, odd(transition.list$from))) != length(transition.list$from) | 
  #   length(subset(transition.list$to,   odd(transition.list$to)))   != length(transition.list$to)
  # ) {
  #   warning("The matrix doesn't satisfy the parity condition.")
  # }

  matrix.data      <- transition.matrix.for.inversion(unpruned.transition.matrix(fa.input, fftbor2d.output))

  pseudo.mfpt.list <- pseudoinverse(matrix.data$inversion.matrix) %*% as.matrix(rep(1, nrow(matrix.data$inversion.matrix)))
  print(pseudo.mfpt.list[matrix.data$index.mapper(matrix.data$bp.dist)])
  # 
  # mfpt.list        <- solve(matrix.data$inversion.matrix) %*% as.matrix(rep(1, nrow(matrix.data$inversion.matrix)))
  # mfpt.to.mfe      <- mfpt.list[matrix.data$index.mapper(matrix.data$bp.dist)]
  # print(mfpt.to.mfe)
  # 
  # mfpt.to.mfe
}

mfpt.from.fa.using.fftbor2d(argv[1])

# benchmark(
#   mfpt.from.fa.using.fftbor2d(argv[1]),
#   replications = 1,
#   columns      = c("elapsed")
# )