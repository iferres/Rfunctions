#pangStats - pangenome statistics
#' @name pangStats
#' @title Compute Pangenome Statistics
#' @description Takes a panmatrix, a coregenome alignment and information about
#' recombination and produces a dataframe with computed fluidity, nucleotide
#' diversity, accesory genome distance, and shared recombination blocks between
#' random sampled genomes. For each row of the \code{data.frame}, the values 
#' are calculated ...
#' 
#' Details on input files below.
#' @param panmat Path to 'gene_presence_absence.csv' file returned by Roary.
#' @param coreali Path to coregenome aligment in fasta format. Roary's 
#' 'core_gene_alignment.aln' is _NOT_ recommended. The output of 
#' \strong{progressiveMauve} is highly suggested.
#' @param gubbgff Path to Gubbins gff output file.
#' @param ssize The size of the subsampled population.
#' @param snum The number of taken subsamples.
#' @param cd A factor between 1 and 0. The fraction of genomes a gene must be
#' in to be considered as a core-gene (notation comes from Roary). This 
#' parameter is important for identifying non-core genes, i.e. accesory genes
#' when calculating the distances between genenomes.
#' @param dist Distance method. See [ade4::dist.binary]. Default: 1, which
#' corresponds to Jaccard index.
#' @details All the inputs should have been produced from the same set of 
#' original genomes, otherwise it will fail.
#' @value A \code{data.frame}.

pangStats <- function (
  
  panmat,
  coreali,
  gubbgff,
  ssize=10,
  snum=500,
  accs=.95,
  dist=1
  
) {
  
  # Load dependencies, install if any is missing
  pkgs <- c('seqinr', 'pegas', 'ade4')
  rq <- sapply(pkgs, require, character.only=TRUE)
  if (!all(rq)){
    print('Installing required packages...')
    install.packages(pkgs[!rq])
    sapply(pkgs[!rq], require, character.only=TRUE)
  }
  
  if(accs>=1 | accs<=0){
    stop('accs should be a factor between 1 and 0.')
  }
  
  # Load shared recombination blocks function
  #  Returns an incidence matrix with the number of shared recombination
  #  blocks between each pair of genomes.
  #  Genome names in both gubbins ('taxa' field) and the panmat must be
  #  the same.
  sharedRec <- function(gubbgff, panmat){
    
    #Creates empthy output matrix
    cln <- strsplit(readLines(panmat, n = 1L), ',')[[1]]
    nams <- gsub('\"','',cln[15:length(cln)])
    norgs <- length(nams)
    m <- as.data.frame(matrix(0L, nrow = norgs, ncol = norgs))
    dimnames(m) <- rep(list(nams), 2L)
    
    #Read gff
    gff <- readLines(gubbgff)
    dat <- sapply(strsplit(gff[-c(1,2)],'\t'),'[',9)
    taxa <- sapply(strsplit(dat, ';'), '[', 3)
    taxa <- strsplit(gsub('taxa=\"|\"','',taxa), ' +')
    taxa <- sapply(taxa, function(x){
      if (x[1]=='') x[-1] else x
    })
    
    #Check
    fctr <- unique(unlist(taxa))
    if(!all(fctr%in%nams)){
      stop("Names in gubbins gff file ('taxa' field) doesn't match with panmatrix names.")
    }
    
    #Compute
    for (i in 1:length(taxa)){
      m[taxa[[i]], taxa[[i]]] <- m[taxa[[i]], taxa[[i]]] + 1L
    }
    diag(m) <- 0L
    
    return(m)
  }
  
  
  
  # Panmatrix
  
  print('[1/5] Processing pan-matrix...')
  
  panmatrix <- read.table(panmat,header=T,sep=',',check.names=F) 
  rnames <- as.vector(panmatrix$Gene) 
  
  panmatrix <- panmatrix[,15:dim(panmatrix)[2]]
  panmatrix <- apply(panmatrix,2,nchar)
  
  panmatrix[which(panmatrix>0)] <- 1L
  rownames(panmatrix) <- rnames
  
  gnum <- dim(panmatrix)[1]
  pnam <- rownames(panmatrix)
  pcom <- combn(pnam,2)
  dimy <- dim(pcom)[2]
  
  if (gnum<=ssize) {
    stop('Sample size must be smaller than number of genomes')
  }
  
  # Core genome
  
  print('[2/5] Processing core alignment...')
  
  fas <- read.alignment(coreali,format='fasta')
  dna <- as.DNAbin(fas)
  nams <- gsub('>','',system(paste('grep ">"',coreali),intern=T))
  rownames(dna) <- nams
  
  # Jaccard distance (Accsesory genome)
  
  print('[3/5] Calculating accesory genome distances...')
  acs <- which(rowSums(panmatrix) < round(ncol(panmatrix)*cd))
  pacs <- t(panmatrix[acs, ])
  dd <- as.matrix(ade4::dist.binary(pacs, method = dist))
  
  # Gff
  
  print('[4/5] Processing Gubbins gff file...')
  
  m <- sharedRec(gubbgff, panmat)
  
  print('[5/5] Calculating...')
  
  # Data analysis
  
  seqn      <- seq(1,snum)
  
  rfluidity <- rep(NA, snum)
  rcoregene <- rep(NA, snum)
  raccsdist <- rep(NA, snum)
  rnumrecom <- rep(NA, snum)
  
  pb <- txtProgressBar(min = 0, max = snum, style = 3)
  for (i in seqn) {
    
    ii <- sample(1:dimy,ssize)
    yy <- pcom[,ii]
    
    # Fluidity
    
    rfluidity[i] <- mean(vapply(1:ssize, function(y){
      p1 <- panmatrix[which(pnam==yy[1,y]), ]
      p2 <- panmatrix[which(pnam==yy[2,y]), ]
      (sum(p1>0L & p2==0L) + sum(p1==0L & p2>0L)) / (sum(p1) + sum(p2))
    }, FUN.VALUE = 1))
    
    # Core genome diversity
    
    snams <- unique(as.vector(yy))
    nmat <- which(nams%in%snams)
    dna2 <- dna[nmat,]
    rcoregene[i] <- nuc.div(dna2)
    
    # Accesory genome Jaccard distance
    
    raccsdist[i] <- mean(vapply(1:ncol(yy), function(y){
      dd[ yy[1,y] , yy[2,y] ]
    }, FUN.VALUE = 1))
    
    
    # Sharing of recombination blocks
    
    rnumrecom[i] <- mean(vapply(1:ncol(yy), function(y){
      m[ yy[1,y] , yy[2,y] ]
      }, FUN.VALUE = 1L))
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  result <- list(rfluidity,rcoregene, raccsdist, rnumrecom)
  names(result) <- c('fluidity','nucdiv', 'accsdist', 'numrec')
  
  result <- as.data.frame(result)
  
  return(result)
}


