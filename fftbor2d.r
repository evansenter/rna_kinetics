# argv <- commandArgs(TRUE)

# if (length(argv) != 1) {
#   cat("./Rscript fftbor2d.r INPUT_FA_FILE\n")
#   q("no")
# }

# if (is.na(file.info(argv[1])$size)) {
#   cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
#   q("no")
# }

library(Matrix)
library(gdata)

mfpt.from.fa.using.fftbor2d <- function(fftbor.input) {
  seq.length    <- as.numeric(system(paste0("ruby -r vienna_rna -e 'p(RNA.from_fasta(\"", fftbor.input, "\").seq.length)'"), intern = T))
  bp.dist       <- as.numeric(system(paste0(
    "ruby -r vienna_rna -e 'p ViennaRna::Rna.bp_distance(*File.read(\"",
    fftbor.input,
    "\").chomp.split(/\n/)[-2..-1])'"), intern = T
  ))
  xition.ncol   <- seq.length + 1
  fftbor.file   <- tempfile()
  write(system(paste("FFTbor2D -S -P 9", fftbor.input), intern = T), fftbor.file)
  fftbor.data   <- read.delim(fftbor.file, header = F, col.names = c("i", "j", "p", "ensemble"))
  fftbor.data   <- within(fftbor.data, ij <- i * xition.ncol + j)
  shift.moves   <- function(x) {
    x + as.numeric(lapply(Map(function(x) x, Map(function(i) { (xition.ncol + 1i) * 1i ^ i }, 0:3)), function(x) Re(x) + Im(x)))
  }
  valid.moves   <- function(x) { Filter(function(y) {
    y >= 0 & sum((c(floor(y / xition.ncol), y %% xition.ncol) - c(floor(x / xition.ncol), x %% xition.ncol)) ^ 2) == 2
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

  transition.matrix <- sparseMatrix(i = transition.list$from, j = transition.list$to, x = transition.list$p, dims = c(xition.ncol ^ 2 - 1, xition.ncol ^ 2 - 1), index1 = F)
  pruned.matrix     <- transition.matrix[-(xition.ncol * bp.dist + 1), -(xition.ncol * bp.dist + 1)]

  if (!all(sapply(apply(transition.matrix, 1, sum), function(x) { x %in% 0:1 }))) {
    warning("One or more row-sums are not 0 or 1")
  }

  rev.utriangle.matrix <- upperTriangle(transition.matrix[, rev(seq_len(ncol(transition.matrix)))], diag = T)
  if (nrow(transition.list[transition.list$p > 0,]) != length(subset(rev.utriangle.matrix, rev.utriangle.matrix > 0))) {
    warning("The y-mirrored transition matrix is not upper triangular (only values in the upper-left quadrant are permitted)")
  }

  mfpt.list <- solve(diag(nrow(pruned.matrix)) - pruned.matrix) %*% as.matrix(rep(1, nrow(pruned.matrix)))
  mfpt.list[bp.dist + 1]
}

synthetic.analysis <- sapply(sort(list.files(path = "synthetic_seq_files", pattern = "seq_*", full.names = T)), function(file) mfpt.from.fa.using.fftbor2d(file))
