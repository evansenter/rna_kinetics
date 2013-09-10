# If we change temperature here, it needs to be changed in the call to RNAeval as well.

argv <- commandArgs(TRUE)

if (!is.element(length(argv), 3:4)) {
  cat("./Rscript eigen_grapher.r EIGEN_BIN_FILE OUTPUT_FILENAME TIME_RANGE (TITLE)\n")
  cat("    TIME_RANGE: Must be an R expression indicating the log_{10}(time) points of interest, ie. seq(-4, 1, 1e-1)\n")
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

log.time.range <- eval(parse(text = argv[3]))
if (!is.numeric(log.time.range)) {
  cat(paste(argv[2], "is not a valid time range, cowering out.\n"))
  q("no")
}

eigen.rate.data <- readRDS(argv[1])
graph.title     <- ifelse(!is.na(argv[4]), argv[4], paste(argv[1], argv[3]))

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

graph.kinetics <- function(rate.data, log.time.range, graph.title) {
  time.range <- Map(function(i) { 10 ^ i }, log.time.range)
  plot(
    log.time.range, 
    Map(function(i) { prob.at.time(eigen.rate.data, eigen.rate.data$key.indices$empty, i) }, time.range), 
    main = graph.title, xlab = bquote(log[10] ~ "(time)"), ylab = "Probability", type = "l", col = "red"
  )

  points(
    log.time.range, 
    Map(function(i) { prob.at.time(eigen.rate.data, eigen.rate.data$key.indices$mfe, i) }, time.range), 
    type = "l", col = "blue"
  )

  legend(
    "topright", 
    c("Empty", "MFE"), 
    lty = c(1, 1), col = c("red", "blue")
  )
}

pdf(argv[2], 6, 6)
graph.kinetics(eigen.rate.data, log.time.range, graph.title)
