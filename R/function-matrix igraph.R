#########################################################
#########################################################
# Functions
# Eric Keen, v. June 2020

#########################################################
#########################################################
# Test data for function development

testdata <- FALSE
if(testdata){

  # Set working directory to folder that this R file is in
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

  dyads <- read.csv("../data/dyads-1.csv",stringsAsFactor=FALSE) ; dyads
  dyads <- ai(dyads)
  var <- "sri"


}

#########################################################
#########################################################

dyads2igraph <- function(dyads,vari="sri"){
  # Convert dyads data i-graph object matrix (also useful for NBDA analyses)

  #vari <- 'sri'

  data <- dyads
  Aid <- unique(as.character(data$A)) ; Aid
  Bid <- unique(as.character(data$B)) ; Bid
  all.ids <- c(Aid,Bid) ; all.ids
  all.ids <- unique(all.ids) ; all.ids

  ABs <- paste0(data$A," ",data$B) ; ABs
  BAs <- paste0(data$B," ",data$A) ; BAs

  MAT <- matrix(data = 0, nrow = length(all.ids), ncol = length(all.ids), byrow = FALSE,
                dimnames = list(all.ids,all.ids))
  head(MAT)

  i=1
  for(i in 1:nrow(data)){
    print(paste(i," out of ",nrow(data)))

    A <- as.character(data$A[i]) ; A
    B <- as.character(data$B[i]) ; B
    X <- data[i,which(names(data)==vari)] ; X

    # A row first
    row.match <- which(row.names(MAT)==A) ; row.match
    col.match <- which(colnames(MAT)==B) ; col.match
    MAT[row.match,col.match] <- as.numeric(X)

    # B row second
    row.match <- which(row.names(MAT)==B) ; row.match
    col.match <- which(colnames(MAT)==A) ; col.match
    MAT[row.match,col.match] <- as.numeric(X)

  }

  MAT[1:6,1:6]
  #hist(MAT,ylim=c(0,100))

  library(igraph)
  g=graph.adjacency(MAT,mode="undirected",weighted=TRUE)
  return(list(g=g,mx=MAT))
}

#########################################################
#########################################################

#dyads <- fast.dyads(hw)
#dyads <- ai(dyads)
#var <- "sri"

dyad2matrix <- function(dyads,var){

  head(dyads)
  data <- dyads
  X <- data[,which(names(data)==var)] ; X

  MAT <- matrix(data = X,
                nrow = length(unique(data$A)),
                ncol = length(unique(data$A)),
                byrow = FALSE,
                dimnames = list(unique(data$a),unique(data$b)))
  MAT[1:6,1:6]

  library(igraph)
  g=graph.adjacency(MAT,mode="undirected",weighted=TRUE)
  return(list(g=g,mx=MAT))

}

#mat <- dyads2igraph(dyads)$mx
#mat[1:6,1:6]


#########################################################
#########################################################
# igraph formatting functions

##############################
# Get values to format Vertices (dots)

## This function essentially extracts values of a specific variable for each 
## node in the graph g, based on matching IDs between the graph nodes and the 
## ID column in the ILV dataframe. It returns these values in a vector, where 
## each element corresponds to a node in the graph. Overall, this function is 
## useful for extracting specific information from a dataframe (ILV) based on 
## the nodes in a graph.

Vvalues <- function(g,ILV,var){
  cname <- which(names(ILV)==var) # finds the index of the column in the ILV dataframe that matches the given variable name.
  Ns <- vector() ; Ns
  i=1
  for(i in 1:length(V(g)$name)){
    idi <- V(g)$name[i] ; idi
    row.match <- which(ILV$id==idi)[1] ; row.match
    ni <- ILV[row.match,cname] ; ni
    Ns[i] <- ni
  }
  Ns
  return(Ns)
}

##############################
# Format edge weight
efunk <- function(drow){
  result <- 1 ;  if(drow$sri.sim.p < 0.05){result <- 5} ; return(result)
}

# Get a value based on edge data

Ecolor <- function(g,dyads,efunk){
  newi <- rep(NA,times=length(E(g))) ; newi
  i=2018
  for(i in 1:length(E(g))){
    E(g)[i]
    indivs <- V(g)[inc(i)]$name ; indivs
    row.match <- which(indivs[1]==dyads$A & indivs[2]==dyads$B) ; row.match
    if(length(row.match)>0){
      newi[i] <- efunk(dyads[row.match,])
    }else{
      print(i)
    }
  }
  newi
  return(newi)
}

##############################
# Color vertex by quantile of ILV value

color.by.quantile <- function(var,varcol=NULL,na.col="grey",legvals){
  var
  if(is.null(varcol)){
    #varcol <- colorRampPalette(rev(wes_palette("Zissou1")))
    varcol <- colorRampPalette(wes_palette("Zissou1"))
  }
  varcol <- varcol(100)
  qs <- quantile(var,na.rm=TRUE,seq(0,1,length=100)) ; qs
  qcol <- rep(na.col,times=length(var))
  for(i in 1:length(var)){
    vari <- var[i]
    if(!is.na(vari)){
      qi <- which.min(abs(var[i] - qs)) ; qi
      qcol[i] <- varcol[qi]
    }
  }
  qcol

  # Legend
  #legvals <- c(0,25,50,75,100) ; legvals
  legcol <- c()
  i=1
  for(i in 1:length(legvals)){
    qi <- which.min(abs(legvals[i] - qs)) ; qi
    legcol[i] <- varcol[qi]
  }
  legvals
  legcol
  return(
    list(varcols=qcol,legval=legvals,legcol=legcol)
  )
}
