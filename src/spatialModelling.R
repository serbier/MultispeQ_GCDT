library(FactoMineR)
library(factoextra)
library(dplyr)
library(missMDA)
library(readxl)
library(asrml)
library(mrbean)
library("corrplot")

setwd("~/Documents/CIAT/GCDT_F56/MultiSpeQ") 

blups <- read_excel('./data/blupsGCDT20-01_parents_imputed_recoded_multicollinear.xlsx', sheet = "Direct_measures")
blups <- blups %>% 
  select(any_of(c('Genotype', "YdPl-M4INV2", "LA-M1_INV2","LA-M2_INV2","FLPLO-M2_INV2","MSWE-M4INV2", "STWPL-M4INV2", "PHIPLO-M4INV2")))

targetCols <- c('f23g','LTD','Ambient.Humidity','Ambient.Temperature','Leaf.Angle','LEF','Light.Intensity..PAR.','NPQt','Phi2','PhiNO','PhiNPQ','PS1.Active.Centers','PS1.Open.Centers','PS1.Over.Reduced.Centers','PS1.Oxidized.Centers','Thickness','sitio.de.cosecha','Ambient.Pressure','Leaf.Temperature','SPAD_650','Device.ID','dayM','phyInx','Carga','row','col','family','Genotype','das','sampling','m_thick', 'Genotype')

gcdtDf <- read.delim("./proccesed_data/20-01_GCDT_MultispecQ_filtered.csv", stringsAsFactors=T, sep = ',', na.strings =c('NA','.',''))
gcdtDf <- gcdtDf %>% 
  mutate(LTD = Leaf.Temperature - Ambient.Temperature) %>% 
  select(any_of(targetCols))


# one sample --------------------------------------------------------------
sample <- '34_morning_INV-M3-2'
s1 <- gcdtDf %>% 
  filter(sampling == sample)

s1$f23g <- as.factor(s1$f23g)


filledflDf <- fill.asreml(s1, rows = 'row', ranges = 'col')
filledflDf$row <- as.factor(filledflDf$row)
filledflDf$col <- as.factor(filledflDf$col)
filledflDf$Device.ID <- as.factor(filledflDf$Device.ID)

whole.m1 <- asreml(SPAD_650 ~1,
                   random =~f23g+Device.ID,
                   residual = ~ar1(row):ar1(col),
                   data = filledflDf[order(filledflDf$row, filledflDf$col),],
                   na.action = na.method(y = c("include", "omit", "fail"), x = c("include", "include","omit")))




dir.create(file.path('./results/', 'spatial_modelling'), showWarnings = FALSE)
dir.create(file.path('./results/spatial_modelling', sample), showWarnings = FALSE)
dir.create(file.path('./results/spatial_modelling', sample,'blups'), showWarnings = FALSE)
dir.create(file.path('./results/spatial_modelling', sample,'residuals'), showWarnings = FALSE)
dir.create(file.path('./results/spatial_modelling', sample,'outs'), showWarnings = FALSE)
dir.create(file.path('./results/spatial_modelling', sample,'variograms'), showWarnings = FALSE)

filledflDf <- fill.asreml(s1, rows = 'row', ranges = 'col', by = 'sitio.de.cosecha')





fill.asreml <- function (x, rows = NULL, ranges = NULL, by, extra) 
{
  if(is.null(rows)|rows=="") return()
  if(is.null(ranges)|ranges=="") return()
  
  checo <- which(colnames(x) %in% c(rows, ranges))
  if (length(checo) != 2) {
    stop("Please double check the rows and ranges argument that you provided. We did not find such columns.\n")
  }
  roro <- x[, which(colnames(x) == rows)]
  raro <- x[, which(colnames(x) == ranges)]
  if (!is.numeric(roro) | !is.numeric(raro)) {
    stop("Please make sure that the columns; ", rows, " and ", 
         ranges, " are both numeric.\n", call. = FALSE)
  }
  if (!missing(extra)) {
    extras = TRUE
  }
  else {
    extras = FALSE
  }
  fac <- which(unlist(lapply(x, class)) == "factor")
  if (length(fac) > 0) {
    for (u in 1:length(fac)) {
      fac2 <- fac[u]
      x[, fac2] <- as.character(x[, fac2])
    }
  }
  if (!missing(by)) {
    y <- split(x, x[, by])
    xnew <- lapply(y, function(x, extra, extras) {
      x <- x[order(x[, rows], x[, ranges]), ]
      roro <- x[, which(colnames(x) == rows)]
      raro <- x[, which(colnames(x) == ranges)]
      bybo <- na.omit(unique(x[, by]))
      ro1 <- seq(min(roro, na.rm = TRUE), max(roro, na.rm = TRUE))
      ra1 <- seq(min(raro, na.rm = TRUE), max(raro, na.rm = TRUE))
      needed <- expand.grid(ro1, ra1)
      colnames(needed) <- c(rows, ranges)
      needed <- needed[order(needed[, rows], needed[, 
                                                    ranges]), ]
      head(needed)
      x <- x[order(x[, rows], x[, ranges]), ]
      head(x)
      dis <- dim(x)
      dis2 <- dim(needed)
      newf <- data.frame(matrix(NA, dis2[1], dis[2]))
      colnames(newf) <- colnames(x)
      head(newf)
      sto <- rownames(x)
      rownames(x) <- paste(x[, rows], x[, ranges], sep = ".")
      rownames(needed) <- rownames(newf) <- paste(needed[, 
                                                         rows], needed[, ranges], sep = ".")
      newf[rownames(x), ] <- x
      newf[, c(rows, ranges)] <- needed
      head(newf)
      rownames(newf) <- NULL
      newf[, by] <- bybo
      if (extras) {
        for (u in 1:length(extra)) {
          pextra <- extra[u]
          pox <- table(newf[, c(rows, ranges, pextra)])
          leve <- dim(pox)[3]
          levelnames <- na.omit(unique(newf[, pextra]))
          for (o in 1:leve) {
            init1 <- which(apply(as.matrix(pox[, , o]), 
                                 1, sum) > 0)
            stend1 <- c(init1[1], init1[length(init1)])
            init2 <- which(apply(as.matrix(pox[, , o]), 
                                 2, sum) > 0)
            stend2 <- c(init2[1], init2[length(init2)])
            kk <- which(newf[, rows] >= stend1[1] & 
                          newf[, rows] <= stend1[2] & newf[, ranges] >= 
                          stend2[1] & newf[, ranges] <= stend2[2])
            newf[kk, pextra] <- levelnames[o]
          }
        }
      }
      return(newf)
    }, extra = extra, extras = extras)
    xnew <- do.call(rbind, xnew)
  }
  else {
    roro <- x[, which(colnames(x) == rows)]
    raro <- x[, which(colnames(x) == ranges)]
    cat("Argument 'by' not provided. Single field assumed.\n")
    ro1 <- seq(min(roro, na.rm = TRUE), max(roro, na.rm = TRUE))
    ra1 <- seq(min(raro, na.rm = TRUE), max(raro, na.rm = TRUE))
    needed <- expand.grid(ro1, ra1)
    colnames(needed) <- c(rows, ranges)
    needed <- needed[order(needed[, rows], needed[, ranges]), 
    ]
    head(needed)
    x <- x[order(x[, rows], x[, ranges]), ]
    head(x)
    dis <- dim(x)
    dis2 <- dim(needed)
    newf <- data.frame(matrix(NA, dis2[1], dis[2]))
    colnames(newf) <- colnames(x)
    head(newf)
    sto <- rownames(x)
    rownames(x) <- paste(x[, rows], x[, ranges], sep = ".")
    rownames(needed) <- rownames(newf) <- paste(needed[, 
                                                       rows], needed[, ranges], sep = ".")
    newf[rownames(x), ] <- x
    newf[, c(rows, ranges)] <- needed
    head(newf)
    rownames(newf) <- NULL
    if (extras) {
      for (u in 1:length(extra)) {
        pextra <- extra[u]
        pox <- table(newf[, c(rows, ranges, pextra)])
        leve <- dim(pox)[3]
        levelnames <- na.omit(unique(newf[, pextra]))
        for (o in 1:leve) {
          init1 <- which(apply(as.matrix(pox[, , o]), 
                               1, sum) > 0)
          stend1 <- c(init1[1], init1[length(init1)])
          init2 <- which(apply(as.matrix(pox[, , o]), 
                               2, sum) > 0)
          stend2 <- c(init2[1], init2[length(init2)])
          kk <- which(newf[, rows] >= stend1[1] & newf[, 
                                                       rows] <= stend1[2] & newf[, ranges] >= stend2[1] & 
                        newf[, ranges] <= stend2[2])
          newf[kk, pextra] <- levelnames[o]
        }
      }
    }
    xnew <- newf
  }
  return(xnew)
}
