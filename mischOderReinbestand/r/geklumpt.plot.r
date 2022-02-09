#Rscript --vanilla pic1.plot.r PDFNAME=../pic/geklumpt.pdf

. <- grep("^PDFNAME=", commandArgs(), value=TRUE)
PDFNAME <- if(length(.) > 0) sub("^PDFNAME=", "", .[1]) else "../pic/geklumpt.pdf"

library(deldir)

set.seed(0)
n <- 100
k <- n/4
l <- 50
x <- matrix(runif(2*n, 0, l), n)

a <- kmeans(x, k)
x <- x + (a$centers[a$cluster,] - x) / 2
a <- deldir(x, rw=c(0,l,0,l))

pdf(PDFNAME)
par(mar=c(0,0,0,0))
plot(a, wlines="tess", showrect = TRUE, showpoints = FALSE)
symbols(x, circles=pmax(.1, sqrt(a$summary$del.area)/20), bg=1, inches=FALSE, add=TRUE)
dev.off()
