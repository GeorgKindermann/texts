#Rscript --vanilla pic1.plot.r PDFNAME=../pic/regelmaessig.pdf

. <- grep("^PDFNAME=", commandArgs(), value=TRUE)
PDFNAME <- if(length(.) > 0) sub("^PDFNAME=", "", .[1]) else "../pic/regelmaessig.pdf"

library(deldir)

set.seed(0)
n <- 400
k <- 100
l <- 50
x <- matrix(runif(2*n, 0, l), n)

a <- kmeans(x, k)
x <- a$centers

a <- deldir(x, rw=c(0,l,0,l))

pdf(PDFNAME)
par(mar=c(0,0,0,0))
#plot(tile.list(a), close = TRUE, showpoints = FALSE, asp=1, axes=FALSE)
plot(a, wlines="tess", showrect = TRUE, showpoints = FALSE)
symbols(x, circles=pmax(.1, sqrt(a$summary$del.area)/20), bg=1, inches=FALSE, add=TRUE)
dev.off()
