argv <- commandArgs(TRUE)

if (length(argv) != 1) {
  cat("./Rscript fftbor2d.r INPUT_FA_FILE\n")
  q("no")
}

if (is.na(file.info(argv[1])$size)) {
  cat(paste(argv[1], "doesn't appear to exist, cowering out.\n"))
  q("no")
}

suppressMessages(library(Matrix))
suppressMessages(require(corpcor))

mfpt.from.fa.using.fftbor2d <- function(fa.input) {
  rt            <- 0.0019872370936902486 * 310.15 # kcal / mol / K
  seq.length    <- as.numeric(system(paste0("ruby -r vienna_rna -e 'p(RNA.from_fasta(\"", fa.input, "\").seq.length)'"), intern = T))
  xition.ncol   <- seq.length + 1
  bp.dist       <- as.numeric(system(paste0(
    "ruby -r vienna_rna -e 'p ViennaRna::Rna.bp_distance(*File.read(\"",
    fa.input,
    "\").chomp.split(/\\n/)[-2..-1])'"), intern = T
  ))

  fftbor2d.output <- tempfile()
  write(system(paste("FFTbor2D -S -E ~/bin/rna_turner_2.1.2.par", fa.input), intern = T), fftbor2d.output)

  fftbor.data   <- read.delim(fftbor2d.output, header = F, col.names = c("i", "j", "p", "ensemble"))
  
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

  fftbor.matrix <- sparseMatrix(
    i      = fftbor.data$i, 
    j      = fftbor.data$j, 
    x      = fftbor.data$p, 
    index1 = F
  )
  
  valid.moves <- function(x) { fftbor.data$ij[fftbor.data$ij != x] }

  transition.list <- expand.grid(fftbor.data$ij, fftbor.data$ij)
  row.names(transition.list) <- NULL

  transition.move.prob <- function(from, to) {
    min(1, (fftbor.data[fftbor.data$ij == to,]$p / fftbor.data[fftbor.data$ij == from,]$p)) / length(valid.moves(from))
  }
  transition.stay.prob <- function(from) {
    1 - Reduce(function(sum, x) sum + transition.move.prob(from, x), valid.moves(from), 0)
  }

  transition.list <- cbind(transition.list, apply(transition.list, 1, function(row) ifelse(row[1] == row[2], transition.stay.prob(row[1]), transition.move.prob(row[1], row[2]))))
  names(transition.list) <- c("from", "to", "p")
  
  uniq.transitions <- unique(c(transition.list$from, transition.list$to))
  mapping          <- cbind(uniq.transitions, order(uniq.transitions))
  index.of         <- function(unconsolidated) { mapping[mapping[,1] == unconsolidated][2] }
  unindex.of       <- function(consolidated) { mapping[mapping[,2] == consolidated][1] }

  # ----------------------------------------------------------------------------------
  # index.of RETURNS A 1-INDEXED LIST OF MAPPINGS FROM THE 0-INDEXED TRANSITION MATRIX
  # ----------------------------------------------------------------------------------
  unpruned.transition.matrix <- sparseMatrix(
    i      = sapply(transition.list$from, index.of), 
    j      = sapply(transition.list$to, index.of), 
    x      = transition.list$p
  )
  
  pruned.matrix <- unpruned.transition.matrix[
    -index.of(xition.ncol * bp.dist), 
    -index.of(xition.ncol * bp.dist) 
  ]
  
  inversion.matrix <- diag(nrow(pruned.matrix)) - pruned.matrix
  pseudo.mfpt.list <- pseudoinverse(inversion.matrix) %*% as.matrix(rep(1, nrow(inversion.matrix)))
  mfpt.list        <- solve(inversion.matrix) %*% as.matrix(rep(1, nrow(inversion.matrix)))
  
  cat(pseudo.mfpt.list[index.of(bp.dist)])
  cat("\n")
  cat(mfpt.list[index.of(bp.dist)])
  cat("\n")
}

mfpt.from.fa.using.fftbor2d(argv[1])
