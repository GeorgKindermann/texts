#Rscript --vanilla SpTab.tab.r TEXTABNAME=../rtab/SpTab.tex

. <- grep("^TEXTABNAME=", commandArgs(), value=TRUE)
TEXTABNAME <- if(length(.) > 0) sub("^TEXTABNAME=", "", .[1]) else "./rtab/SpTab.tex"

x <- read.csv("./tabs/SpTab.csv", comment.char = "#")
x[] <- lapply(x, function(x) if(is.character(x)) type.convert(gsub(",", ".", x), as.is=TRUE) else x)

FILE=TEXTABNAME

rTyp <- paste0("\\typ{", x$Typ, "}{", x$Immergruen, "}")

rT <- ifelse(is.na(x$Tl), "", paste0("\\ratingC{", 0, "}{", x$Tl, "}"))
rP <- ifelse(is.na(x$Pl) | is.na(x$Ph), "", paste0("\\ratingC{", 100-x$Pl, "}{", 100-x$Ph, "}"))
rL <- ifelse(is.na(x$Ll), "", paste0("\\ratingC{", 0, "}{", x$Ll, "}"))
rR <- ifelse(is.na(x$Rl) | is.na(x$Rh), "", paste0("\\ratingC{", x$Rl, "}{", x$Rh, "}"))
rN <- ifelse(is.na(x$Nl), "", paste0("\\ratingC{", 0, "}{", x$Nl, "}"))
rK <- ifelse(is.na(x$Kl), "", paste0("\\ratingC{", 0, "}{", x$Kl, "}"))
NL <- rep(r"(\\)", nrow(x))
NL[c(rep(FALSE, 4), TRUE)] <- r"(\\[.3em])"
NL[length(NL)] <- ""

#rSalz <- ifelse(is.na(x$Salztolerant), "", paste0("\\jaNein{", x$Salztolerant, "}"))
#rSMet <- ifelse(is.na(x$Schwermetalltolerant), "", paste0("\\jaNein{", x$Schwermetalltolerant, "}"))
#rWind <- ifelse(is.na(x$Windbestaeubt) | x$Windbestaeubt!=1, "", paste0("\\jaNein{1}"))
#rNFix <- ifelse(is.na(x$NFix), "", paste0("\\jaNein{", x$NFix, "}"))
rSalzMetWindN <- paste0("\\jaNeinQ{", x$Salztolerant, "}{", x$Schwermetalltolerant, "}{", x$Windbestaeubt, "}{" , x$NFix, "}")


#rGift <- ifelse(is.na(x$giftig) | x$giftig!=1, "", paste0("\\jaNein{1}"))
#rEssbar <- ifelse(is.na(x$essbar) | x$essbar!=1, "", paste0("\\jaNein{1}"))
rGiftigEssbar <- paste0("\\jaNeinP{", x$giftig, "}{", x$essbar, "}")

rAlter <- ifelse(is.na(x$maxAlter), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$maxAlter/4), "}"))
rWurzel <- ifelse(is.na(x$Wurzeltiefe), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Wurzeltiefe/3), "}"))
rSamAlter <- ifelse(is.na(x$Samenueberdauerung), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Samenueberdauerung), "}"))
rHoehe <- ifelse(is.na(x$Baumhoehe), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Baumhoehe/6*10), "}"))

rKonk <- ifelse(is.na(x$KonkStrat), "", paste0("\\ratingC{", 0, "}{", c("Konk" = 100, "Stress" = 50, "Ruderal" = 0), "}"))

rH15 <- ifelse(is.na(x$HoeheAlter15), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$HoeheAlter15/16*100), "}"))
rZiel <- ifelse(is.na(x$Planbarkeit), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Planbarkeit), "}"))
rStreu <- ifelse(is.na(x$Streuabbau), "", paste0("\\ratingC{", 0, "}{", pmin(100, pmax(0, 5 - x$Streuabbau)*20), "}"))
rAusschlag <- ifelse(is.na(x$Stockausschlag), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Stockausschlag), "}"))
rWurzelbrut <- ifelse(is.na(x$Wurzelbrut), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Wurzelbrut), "}"))
rDichte <- ifelse(is.na(x$Holzdichte), "", paste0("\\ratingC{", 0, "}{", pmin(100, pmax(0, x$Holzdichte - 350)/4.5 ), "}"))
rDauer <- ifelse(is.na(x$DauerhaftigkeitDesHolzes), "", paste0("\\ratingC{", 0, "}{", pmin(100, pmax(0, 5-x$DauerhaftigkeitDesHolzes)*25), "}"))

rMaeuse <- ifelse(is.na(x$Maeuse), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Maeuse), "}"))
rVerbiss <- ifelse(is.na(x$Verbiss), "", paste0("\\ratingC{", 0, "}{", pmin(100,100-(x$Verbiss-1)*100/3), "}"))
rFegenSchlagen <- ifelse(is.na(x$FegenSchlagen), "", paste0("\\ratingC{", 0, "}{", pmin(100,100-(x$FegenSchlagen)*100/3), "}"))
rSchaelen <- ifelse(is.na(x$Schaelen), "", paste0("\\ratingC{", 0, "}{", pmin(100,100-(x$Schaelen-1)*100/3), "}"))

rWindwurf <- ifelse(is.na(x$Windwurf), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Windwurf), "}"))
rSchneebruch <- ifelse(is.na(x$Schneebruch), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Schneebruch), "}"))
rPilze <- ifelse(is.na(x$Pilze), "", paste0("\\ratingC{", 0, "}{", pmin(100,x$Pilze), "}"))

#Typ
#Immergrün

#cat(paste0(paste(x$co, x$sp, x$de, x$ba, rT, rK, rP, rL, rR, rN, rNFix, sep=" & "), NL), file=FILE, sep="\n")

cat(paste0(paste(x$co, x$sp, x$de, x$ba, rTyp, rHoehe, rH15, rAlter, rKonk, rWurzel, rT, rK, rP, rL, rR, rN, rSalzMetWindN, rSamAlter, rGiftigEssbar, rZiel, rStreu, rAusschlag, rWurzelbrut, rDichte, rDauer, rMaeuse, rVerbiss, rFegenSchlagen, rSchaelen, rWindwurf, rSchneebruch, rPilze, sep=" & "), NL), file=FILE, sep="\n")


