#Rscript --vanilla SpTab.tab.r TEXTABNAME=../rtab/SpTab.tex

. <- grep("^TEXTABNAME=", commandArgs(), value=TRUE)
TEXTABNAME <- if(length(.) > 0) sub("^TEXTABNAME=", "", .[1]) else "./rtab/SpTab.tex"

x <- read.csv("./tabs/SpTab.csv")

FILE=TEXTABNAME

rT <- ifelse(is.na(x$Tl), "", paste0("\\ratingC{", 0, "}{", x$Tl, "}"))
rP <- ifelse(is.na(x$Pl) | is.na(x$Ph), "", paste0("\\ratingC{", 100-x$Pl, "}{", 100-x$Ph, "}"))
rL <- ifelse(is.na(x$Ll), "", paste0("\\ratingC{", 0, "}{", x$Ll, "}"))
rR <- ifelse(is.na(x$Rl) | is.na(x$Rh), "", paste0("\\ratingC{", x$Rl, "}{", x$Rh, "}"))
rN <- ifelse(is.na(x$Nl), "", paste0("\\ratingC{", 0, "}{", x$Nl, "}"))
rK <- ifelse(is.na(x$Kl), "", paste0("\\ratingC{", 0, "}{", x$Kl, "}"))
rNFix <- ifelse(is.na(x$NFix), "", paste0("\\jaNein{", x$NFix, "}"))
NL <- rep(r"(\\)", nrow(x))
NL[c(rep(FALSE, 4), TRUE)] <- r"(\\[.3em])"
NL[length(NL)] <- ""

rSalz <- ifelse(is.na(x$Salztolerant), "", paste0("\\jaNein{", x$Salztolerant, "}"))
rSMet <- ifelse(is.na(x$Schwermetalltolerant), "", paste0("\\jaNein{", x$Schwermetalltolerant, "}"))
rWind <- ifelse(is.na(x$Windbestaeubt) | x$Windbestaeubt!=1, "", paste0("\\jaNein{1}"))
rGift <- ifelse(is.na(x$giftig) | x$giftig!=1, "", paste0("\\jaNein{1}"))

rAlter <- ifelse(is.na(x$maxAlter), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$maxAlter/4), "}"))
rWurzel <- ifelse(is.na(x$Wurzeltiefe), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Wurzeltiefe/3), "}"))
rSamAlter <- ifelse(is.na(x$Samenueberdauerung), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Samenueberdauerung), "}"))
rHoehe <- ifelse(is.na(x$Baumhoehe), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Baumhoehe/6*10), "}"))

rKonk <- ifelse(is.na(x$KonkStrat), "", paste0("\\ratingC{", 0, "}{", c("Konk" = 100, "Stress" = 50, "Ruderal" = 0), "}"))

#cat(paste0(paste(x$co, x$sp, x$de, x$ba, rT, rK, rP, rL, rR, rN, rNFix, sep=" & "), NL), file=FILE, sep="\n")

cat(paste0(paste(x$co, x$sp, x$de, x$ba, rHoehe, rAlter, rKonk, rWurzel, rT, rK, rP, rL, rR, rN, rNFix, rSalz, rSMet, rWind, rSamAlter, rGift, sep=" & "), NL), file=FILE, sep="\n")


