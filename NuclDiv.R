NuclDiv <- function(fastas){
  archivos <- fastas
  sec <- list()
  for(i in 1:length(archivos)){
    sec.i <- readDNAStringSet(archivos[i])
    sec.i <- unlist(sec.i)
    sec[[i]] <- sec.i
  }
  par <- combn(1:length(sec), m = 2, simplify = T)
  nucl.div <- NULL
  for(c in 1:ncol(par)){
    secs <- DNAStringSet(list(sec[[par[1,c]]], sec[[par[2,c]]]))
    ali <- AlignSeqs(secs)
    ali.b <- as.DNAbin(ali)
    nucl.div[c] <- pegas::nuc.div(ali.b)
  }
  
  boxplot.default(nucl.div)
  NucleotideDiversity <<- nucl.div
}