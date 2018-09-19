#Builds a 'coregenome' using a xmfa (as returned by progressiveMauve) 
# and ordering the LCBs according to the order in reference. Returns a FASTA
# file with the 'core' alignment.
#Selects those LCBs which have sequence from the reference, and assumes these 
# LCBs are part of the coregenome.
#It is guaranteed that the complete reference genome will appear in the final
# alignment in the corresponding order, not the others. 
#If LCBs do not contain some of the genomes (except the reference, which 
# always should be present), this function adds gaps to fill these spaces.
#NOTE: progressiveMauve returns LBCs with many gaps in some cases, and also 
# this function, as mentioned above, adds gaps in cases where the LCB don't 
# contain all the taxa (nco parameter). Don't freak out if you see many gaps in
# the final alignment. Understand what this function does in first place. If
# you want to build a phylogeny, you may find appropiate to trim the alignment 
# first.
#Suitable to pass output to gubbins (recombination detection), in theory. If
# you are considering this, then you shouldn't trim the alignment because you
# would be creating too many artifitial junctions. This script tries to avoid 
# this issue by using the order of a reference, but is inevitable that some of
# the other genomes will be artificialy glued together, specially if the taxa
# have very plastic and variable genomes. This is an open problem, as long as
# I know. Another thing to consider is to extract complete LCBs (blocks with 
# all taxa represented) and run gubbins over separated fasta files, each one 
# with a different LCB.
#Tips to choose reference: try to select a closed genome, with few 'n's, and 
# the largest one. If more than one chromosome available, recomended to subset
# individual closed chromosomes first, and run progressiveMauve using this
# subset and the rest of genomes instead of runing all together. After 
# progressiveMauve, run this script over the xmfa output.

#NOTE2: I'm thinking that there may be an issue when concatenates the LCBs, 
# because at it is now, is not considering the orientation of the sequences 
# (information that is given in the xmfa file). I'm not sure if in the final
# alignment I'm pasting together blocks in the correct order but incorrect
# orientation. I put this here as a personal reminder to check this out later. 
# Any thoughts are welcome.

#Usage: 
	# xmfa: path to xmfa file
	# dout: output directory name. Creates it of doesnt exists.
	# ref: A reference to use as guide. Dont use the complete or relative
	#  path, just the basename of the reference.
	#  Example: If reference is at ~/hello/world/reference.fasta, then just
	#  use 'reference.fasta'.
	# nco: the number of genomes in the xmfa alignment.
buildCoreGenome <- function(xmfa,
                            dout,
                            ref, 
			    nco){ 
  
  if (missing(nco)){
	  stop('You must provide the number of genomes present in the alignment')
  }
  
  if (missing(ref)){
    stop('ref is empthy, please select a reference')
  }
  
  if (missing(dout)){
    dout <- 'out/'
  }
  
  if(!dir.exists(dout)){
    dir.create(dout)
  }
  
  ncom <- (nco * 2L) + 2L
  rl <- readLines(xmfa)
  
  na <- rl[2*1:((ncom/2)-1)]
  idx <- 1:length(na)
  na <- sapply(lapply(strsplit(na, '/'), rev), '[', 1)
  idx <- data.frame(Index=idx, Basename=na)
  #ref on idx?
  if(!ref%in%idx$Basename){
    stop('ref is not in file')
  }
  
  rl <- rl[-(1L:ncom)]
  
  #index of chunks (blocks) of sequences in xmfa
  eq <- which(rl=='=')
  vp <- vapply(1:length(eq), function(x){
    
    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)
    
  }, FUN.VALUE = c(1L,1L))
  
  
  #which indexes are on each chunk
  whix <- apply(vp, 2, function(x){
    rr <- rl[x[1]:x[2]]
    gp <- grep('> ', rr, fixed = TRUE, value = TRUE)
    as.integer(gsub('^[>] |[:].+','', gp))
  })
  
  #which chunks have reference on it
  refix <- idx$Index[which(idx$Basename==ref)]
  whre <- which(vapply(whix, function(x){refix%in%x}, FUN.VALUE = NA))
  
  
  #Identify start and end of reference to sort them
    #First identify lines in rl where those chunks are located
  vpref <- vp[, whre]
  
    #Second, identify line on each chunk and subset start and end
  pos <- apply(vpref, 2, function(x){
    rr <- rl[x[1]:x[2]]
    gp <- grep(paste0('> ', refix, ':'), rr, fixed = TRUE, value = TRUE)
    gs <- gsub('> \\d{1,2}:| ([+]|[-]) .+', '', gp)
  })
  
  pos <- lapply(strsplit(pos, '-'), as.integer)
  pos <- as.data.frame(do.call(rbind, pos))
  
    #Third, order rows by start
  or <- order(pos$V1)
  
  #Order chunks
  orch <- vpref[, or]
  
  #Extract sequences from ordered chunks, write to temporary fasta
   #create files
  tmps <- lapply(1:nco, function(x){
    # txt <- paste0('>', idx$Basename[which(idx$Index==x)])
    con <- paste0(dout,'/tmp.seq', x, '.fasta')
    file(con, open = 'w')
  })
   #headers
  lapply(1:nco, function(x) {
    txt <- paste0('>', idx$Basename[which(idx$Index==x)])
    writeLines(txt, con = tmps[[x]], sep = '\n')
  })
  apply(orch, 2, function(x){
    #chunk to sequence list
    sq <- chunk2seq(rl[x[1]:x[2]]) #function defined below
    #which not present
    whno <- which(!(1:nco)%in%as.integer(names(sq)))
    if(length(whno>0)){
      sql <- nchar(sq[[1]])
      #fill
      sq[as.character(whno)] <- strrep('-', sql)
      #order
      sq <- sq[as.character(1:nco)]
    }
    #write temporary files with concatenated lcbs per genome
    sapply(1:nco, function(x){
      writeLines(sq[[x]], con = tmps[[x]], sep = '')
    })
    
  })
  #enter
  lapply(tmps, function(x) {writeLines('\n', con = x)})
  #close
  lapply(tmps, close)
  
  #Append temp files and output
  tmps <- paste0(dout,'tmp.seq', 1:nco, '.fasta')
  out <- paste0(dout, 'Ref_',ref,'_alignment.fasta')
  file.create(out)
  for (x in tmps){
   file.append(out, x) 
  }
  
  #Return
  out
}


###############################################################################



#Taken and adapted from seqinr package.
#Takes a string vector (result of >readLines(*aFastaFile*) ) and transforms it
# into a sequence list.
chunk2seq <- function(lines){
  
  ind <- which(substr(lines, 1L, 1L) == ">")
  nseq <- length(ind)
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                       collapse = ""))
  
  nomseq <- sapply(ind, function(i){gsub('^[>] |[:].+', '', lines[i])})
  
  sequences <- as.list(tolower(sequences))
  
  names(sequences) <- nomseq
  return(sequences)
  
}



