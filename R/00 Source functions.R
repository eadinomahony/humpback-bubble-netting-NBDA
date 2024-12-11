#########################################################
#########################################################
# Source all function files

# Set working directory to folder that this R file is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################################################
#########################################################

lf <- list.files()
funks <- grep("function",lf) ; funks
lfunk <- lf[funks] ; lfunk
lfunk <- lfunk[-grep("00",lfunk)] ; lfunk

i=1
for(i in 1:length(lfunk)){
  lfi <- lfunk[i] ; print(lfi)
  source(lfi)
}

#########################################################
#########################################################
