library(FactoMineR)
library(factoextra)
library(dplyr)
library(missMDA)
library(readxl)
library("corrplot")
setwd("~/Documents/CIAT/GCDT_F56/MultiSpeQ") 
blups <- read_excel('./data/blupsGCDT20-01_parents_imputed_recoded_multicollinear.xlsx', sheet = "Direct_measures")
blups <- blups %>% 
  select(any_of(c('Genotype', "YdPl-M4INV2", "LA-M1_INV2","LA-M2_INV2","FLPLO-M2_INV2","MSWE-M4INV2", "STWPL-M4INV2", "PHIPLO-M4INV2")))

targetCols <- c('LTD','Ambient.Humidity','Ambient.Temperature','Leaf.Angle','LEF','Light.Intensity..PAR.','NPQt','Phi2','PhiNO','PhiNPQ','PS1.Active.Centers','PS1.Open.Centers','PS1.Over.Reduced.Centers','PS1.Oxidized.Centers','Thickness','sitio.de.cosecha','Ambient.Pressure','Leaf.Temperature','SPAD_650','Device.ID','dayM','phyInx','Carga','row','col','family','Genotype','das','sampling','m_thick', 'Genotype')
gcdtDf <- read.delim("./proccesed_data/20-01_GCDT_MultispecQ_filtered.csv", stringsAsFactors=T, sep = ',', na.strings =c('NA','.',''))
gcdtDf <- gcdtDf %>% 
  mutate(LTD = Leaf.Temperature - Ambient.Temperature) %>% 
  select(any_of(targetCols))



# test for one sampling ---------------------------------------------------

#results PCA
dir.create(file.path('./results', 'pca'), showWarnings = FALSE)
for (sample in unique(gcdtDf$sampling)){
  dir.create(file.path('./results/pca', sample), showWarnings = FALSE)
  #get sampling df
  s1 <- gcdtDf %>% 
    filter(sampling == sample)
  
  #filter by missing data by column  
  s1md <- s1[, which(colMeans(!is.na(s1)) > 0.8)]
  
  #merge with direct measures
  s1merg <- merge(x = s1md, y = blups, by = "Genotype", all.x=TRUE)
  
  #remove str Cols
  strCols <- c('Genotype','Leaf.Angle','Device.ID','Thickness','sitio.de.cosecha','dayM','Carga','row','col','family','das','sampling')
  pcaD <- s1merg %>% 
    select(!any_of(strCols))
  
  #impute by pca method
  res.comp <- imputePCA(pcaD, ncp = 5)
  
  supVars <- match(c("Ambient.Humidity","Ambient.Temperature","Light.Intensity..PAR.","Ambient.Pressure","Leaf.Temperature","YdPl-M4INV2",
          "LA-M1_INV2", "LA-M2_INV2", "FLPLO-M2_INV2", "MSWE-M4INV2", "STWPL-M4INV2", "PHIPLO-M4INV2"),colnames(res.comp$completeObs))
  
  #PCA
  res.pca = PCA(res.comp$completeObs,ind.sup = c(3,4,5),quanti.sup = supVars, graph = F)
  #informative Plots
  var <- get_pca_var(res.pca)
  eig.val <- as.data.frame(get_eigenvalue(res.pca))
  #write eigenvalues
  write.table(eig.val,file=paste0('./results/pca/',sample,"/",sample,"_eig_val.tsv"), quote = F, sep='\t')
  # save indv plot by device ID.
  pdf(file=paste0('./results/pca/',sample,"/",sample,"_plots.pdf"))
  indPlot <- fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             habillage = s1merg[-c(3,4,5),]$Device.ID, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Device ID",
             title=sample)
  
  biplotG <- fviz_pca_var(res.pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
               repel = TRUE, # Avoid text overlapping
               title = sample )
  corPlot <- corrplot(var$cos2,title=sample, is.corr=FALSE)
  print(indPlot)
  print(biplotG)
  print(corPlot)
  dev.off()

}














