tmp <- commandArgs(trailingOnly=TRUE)

x <- read.table(paste0(tmp, "/border.txt"))
x <- x[-nrow(x),]
i <- duplicated(x) & !is.na(x[,1])
. <- diff(i)
. <- c(0, .) == 1 | c(., 0) == -1
x$fix <- do.call(paste, x) %in% do.call(paste, x[.,])
x[i & !x$fix,1:2] <- NA
. <- is.na(x[,1])
x <- x[!(c(.[-1], FALSE) & .),]

. <- lapply(split(x, cumsum(is.na(x[,1]))), \(x) {
    . <- x[!is.na(x[,1]),]
    if(nrow(.) > 1) .
    else NULL
})
. <- .[lengths(.) > 0]							      

P <- do.call(rbind, lapply(., \(x) {
    rbind(NA, round(aggregate(x, list(cumsum(c(0, sqrt(diff(x[,1])^2 + diff(x[,2])^2)) %/% 200), cumsum(c(0, abs(diff(x$fix))))), mean)[c("V1", "V2")]))
	     }))[-1,]

ext <- matrix(c(177000, 5013027, 815612, 5452628), 2)

pdf(paste0(tmp, "/border.pdf"), diff(ext[1,])/1e5, diff(ext[2,])/1e5)
par(mar=rep(0, 4),mai=rep(0, 4), xpd = NA)
plot(P, type="l", col="#FF00007F", xlim=ext[1,], ylim=ext[2,], axes=FALSE, frame.plot=FALSE, ann=FALSE, asp=1)
points(ext[1,], ext[2,], pch=".", col="#01000001")
dev.off()
