#Rscript --vanilla pic1.plot.r PDFNAME=../pic/zufaellig.pdf

. <- grep("^PDFNAME=", commandArgs(), value=TRUE)
PDFNAME <- if(length(.) > 0) sub("^PDFNAME=", "", .[1]) else "../pic/zufaellig.pdf"

library(deldir)

set.seed(0)
n <- 100
l <- 50
x <- data.frame(x=runif(n, 0, l), y=runif(n, 0, l))
a <- deldir(x$x, x$y, rw=c(0,l,0,l))
x$d <- sqrt(a$summary$del.area)/20

pdf(PDFNAME)
par(mar=c(0,0,0,0))
plot(a, wlines="tess", showrect = TRUE, showpoints = FALSE)
symbols(x$x, x$y, circles=pmax(.1, x$d), bg=1, inches=FALSE, add=TRUE)
dev.off()
