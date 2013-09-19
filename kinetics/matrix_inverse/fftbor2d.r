argv <- commandArgs(TRUE)

if (!is.element(length(argv), c(1, 3))) {
  cat("./Rscript fftbor2d.r INPUT_FA_FILE\n")
  cat("./Rscript fftbor2d.r SEQUENCE STR_1 STR_2\n")
  q("no")
}

if (length(argv) == 1 & is.na(file.info(argv[1])$size)) {
  cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
  q("no")
}

# Easiest way to homogenize the interface to the actual code, if data was provided at the command line rather
# than via a file, just make a tempfile and put it in there.
if (length(argv) == 1) {
  fa.input <- argv[1]
} else if (length(argv) == 3) {
  fa.input <- tempfile()
  write(do.call(paste, c(as.list(argv), sep = "\n")), fa.input)
}

# Matrix gives us sparseMatrix, corpcor gives us pseudoinverse.
require(corpcor, quietly = T)

mfpt.from.fa.using.fftbor2d <- function(fa.input) {
  fa.data <- suppressWarnings(readLines(fa.input))
  fa.data <- fa.data[seq(length(fa.data) - 2, length(fa.data))]
  
  seq.length    <- nchar(fa.data[1])
  xition.ncol   <- seq.length + 1

  fftbor2d.output <- tempfile()
  write(system(paste("FFTbor2D -E ~/bin/rna_turner_2.1.2.par", fa.input), intern = T), fftbor2d.output)
  
  bp.dist       <- read.csv(text = readLines(fftbor2d.output, 4)[4], header = F)[1, 1]
  fftbor.data   <- read.delim(
    text = do.call(paste, c(as.list(readLines(fftbor2d.output)[-1:-5]), sep = "\n")), 
    header = F, 
    col.names = c("i", "j", "p", "ensemble")
  )
  
  if (bp.dist <= 0) {
    cat("BP distance between the two structures (bp.dist) looks weird. Quitting.\n")
    q("no")
  }
  
  if (nrow(fftbor.data) <= 0) {
    cat("FFTbor data in fftbor.data looks weird. Quitting.\n")
    q("no")
  }
  
  if (nrow(fftbor.data[fftbor.data$i == bp.dist & fftbor.data$j == 0,]) == 0) {
    cat("p(bp.dist, 0) = 0\nINFINITY\n")
    q("no")
  }
  
  if (nrow(fftbor.data[fftbor.data$i == 0 & fftbor.data$j == bp.dist,]) == 0) {
    print("p(0, bp.dist) = 0\nINFINITY\n")
    q("no")
  }
  
  two.d.to.rmoi <- function(i, j) {
    i * xition.ncol + j
  }
  
  rmoi.to.two.d <- function(ij) {
    c(floor(ij / xition.ncol), ij %% xition.ncol)
  }

  fftbor.data <- within(fftbor.data, ij <- two.d.to.rmoi(i, j))
  fftbor.data <- fftbor.data[order(fftbor.data$ij),]
  move.size   <- nrow(fftbor.data) - 1
  
  valid.moves <- function(x) { fftbor.data$ij[fftbor.data$ij != x] }

  transition.move.prob <- function(from, to) {
    min(1, (fftbor.data[fftbor.data$ij == to,]$p / fftbor.data[fftbor.data$ij == from,]$p)) / move.size
  }

  transition.list            <- expand.grid(fftbor.data$ij, fftbor.data$ij)
  transition.list            <- cbind(transition.list, rep(NA, nrow(transition.list)))
  names(transition.list)     <- c("from", "to", "p")
  
  transition.list[transition.list$from != transition.list$to,]$p <- apply(transition.list[transition.list$from != transition.list$to,], 1, function(row) transition.move.prob(row[1], row[2]))
  
  transition.list[transition.list$from == transition.list$to,]$p <- sapply(fftbor.data$ij, function(index) 1 - sum(transition.list[transition.list$from == index & transition.list$to != index,]$p))
  
  mapping    <- cbind(fftbor.data$ij, 1:nrow(fftbor.data))
  index.of   <- function(unconsolidated) { mapping[mapping[,1] == unconsolidated][2] }
  unindex.of <- function(consolidated) { mapping[mapping[,2] == consolidated][1] }

  # --------------------------------------------------------------------------------------
  # index.of RETURNS A 1-INDEXED LIST OF MAPPINGS FROM THE 0-INDEXED TRANSITION DATA FRAME
  # --------------------------------------------------------------------------------------
  pruned.matrix <- matrix(data = transition.list$p, nrow = nrow(fftbor.data), ncol = nrow(fftbor.data))[-index.of(xition.ncol * bp.dist), -index.of(xition.ncol * bp.dist)]
  
  inversion.matrix <- diag(nrow(pruned.matrix)) - pruned.matrix
  pseudo.mfpt.list <- pseudoinverse(inversion.matrix) %*% as.matrix(rep(1, nrow(inversion.matrix)))
  mfpt.list        <- solve(inversion.matrix) %*% as.matrix(rep(1, nrow(inversion.matrix)))
  
  cat(pseudo.mfpt.list[index.of(bp.dist)])
  cat("\n")
  cat(mfpt.list[index.of(bp.dist)])
  cat("\n")
}

mfpt.from.fa.using.fftbor2d(fa.input)
