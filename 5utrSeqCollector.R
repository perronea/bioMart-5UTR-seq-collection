library(biomaRt)
library(seqinr)
#library(foreach) #not sure how to best implement MC functionality
#library(doMC)
#registerDoMC()

#set working directory
#system("mkdir /Users/andersperrone/workDir/genomeNetSeqs")
#system("mkdir /Users/andersperrone/workDir/genomeNetSeqs/5utrSeq")
#system("mkdir /Users/andersperrone/workDir/genomeNetSeqs/randomUpstream")
setwd("/Users/andersperrone/workDir/genomeNetSeqs")

#list all species in the ensembl mart
speciesList <- listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))[,1]
#list all species already created
collectedSpecies <-system("ls \ 5utrSeq", intern = T)

#Creates a data frame to summarize the # of sequences that have been collected
summary_df <- data.frame("speciesName"=character(),"all"=numeric(),"filtered"=numeric(),"seqs"=numeric(),"up"=numeric(), stringsAsFactors = FALSE)

for (species in setdiff(speciesList, collectedSpecies)){
  print(paste("Current species", species))
  dir.create(paste("5utrSeq/", species, sep="")) #creates species dir in 5utr
  dir.create(paste("randomUpstream/", species, sep="")) #creates species dir in rU
  #get the species mart
  currSpeciesMart <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset=species)
  currSpeciesGenes <- getBM(
    attributes=c("ensembl_gene_id","gene_biotype"), 
    mart=currSpeciesMart,
    checkFilters=FALSE)
  print(paste("# all 5'UTRs: ", length(currSpeciesGenes[,1])))

  #filter for protein_coding genes
  filteredGenes <- currSpeciesGenes[currSpeciesGenes[,2]=="protein_coding",]
  geneIds <- filteredGenes[,1]
  print(paste("# protein coding 5'UTRs: ", length(geneIds)))

  print("Collecting sequences in blocks of 1000...")
  #break up sequence queries into blocks of 1000
  genes1000 <- length(geneIds) %/% 1000
  genesRem <- length(geneIds) %% 1000
  #create data frame to hold all ids and sequences
  seq_df <- data.frame("gene_id"=character(), "utr"=character(), stringsAsFactors = FALSE)
  if(genes1000>0) {
    for(geneBlock in 1:genes1000) {
      print(geneBlock)
      #try-catch for stochastic getBM error caused possible by server overload
      tryResult <- try(biomaRt::getSequence(
        id=geneIds[((geneBlock-1)*1000+1):(geneBlock*1000)],
        type='ensembl_gene_id',
        seqType='5utr',
        mart=currSpeciesMart
      ))
      numTries <- 0
      #while there's an error print error and try again
      while(inherits(tryResult, "try-error") & numTries<100) {
        tryResult <- try(biomaRt::getSequence(id=geneIds[((geneBlock-1)*1000+1):(geneBlock*1000)],type='ensembl_gene_id',seqType='5utr',mart=currSpeciesMart))
        numTries <- numTries+1
        print(paste("Number tries ", as.character(numTries), sep=''))	
      }
      if(!(inherits(tryResult, "try-error"))) {
        utrSeqs <- tryResult
        #check each utr for length and add it to the data frame
        for(i in 1:1000) {
          if(nchar(utrSeqs[i,1])<21 | nchar(utrSeqs[i,1])>600) next()
          new_utr <- c(utrSeqs[i,2], utrSeqs[i,1])
          #seq_df <- rbind(seq_df, setNames(as.list(new_utr), names(seq_df)))
          seq_df[nrow(seq_df)+1,] <- new_utr
        }
      }
    }
  }
  tryResult <- try(biomaRt::getSequence(
    id =geneIds[(geneBlock*1000+1):(geneBlock*1000+genesRem)],
    type='ensembl_gene_id',
    seqType='5utr',
    mart=currSpeciesMart
  ))
  numTries <- 0
  while(inherits(tryResult, "try-error") & numTries<100) {
    tryResult <- try(biomaRt::getSequence(id =geneIds[(geneBlock*1000+1):(geneBlock*1000+genesRem)],type='ensembl_gene_id',seqType='5utr',mart=currSpeciesMart))
    numTries <- numTries+1
    print(paste("Number tries ", as.character(numTries), sep=''))	
  }
  if(!(inherits(tryResult, "try-error"))) {
    remUtrSeqs <- tryResult
    #sometimes the remaining ids do not get sequences back so remUtrSeqs will be empty, but if it's not...
    if(length(remUtrSeqs[,1])!=0) {
      for(i in 1:length(remUtrSeqs[,1])) {
        if(nchar(remUtrSeqs[i,1])<21 | nchar(utrSeqs[i,1])>600) next()
        new_utr <- c(remUtrSeqs[i,2], remUtrSeqs[i,1])
        seq_df[nrow(seq_df)+1,] <- new_utr
      }
      print(paste(length(seq_df$gene_id), "valid 5'UTRs"))
    }
  }
  
  
  #####Time to check the upstream sequences######
  
  
  utr1000 <- length(seq_df$gene_id) %/% 1000
  utrRem <- length(seq_df$gene_id) %% 1000
  utrBlock <- 0
  #add a column to the dataframe to record the upstream sequence
  upstream <- vector(,length(seq_df$gene_id))
  seq_df <- cbind(seq_df, upstream)
  print("Checking upstream region...")
  if(utr1000>0) {
    for(utrBlock in 1:utr1000) {
      print(utrBlock)
      tryResult <- try(biomaRt::getSequence(
        id=geneIds[((utrBlock-1)*1000+1):(utrBlock*1000)],
        type='ensembl_gene_id',
        seqType='coding_gene_flank',
        upstream=5000,
        mart=currSpeciesMart
      ))
      numTries <- 0
      while(inherits(tryResult, "try-error") & numTries<100) {
        tryResult <- try(biomaRt::getSequence(id=geneIds[((utrBlock-1)*1000+1):(utrBlock*1000)],type='ensembl_gene_id',seqType='coding_gene_flank',upstream=5000,mart=currSpeciesMart))
        numTries <- numTries+1
        print(paste("Number tries ", as.character(numTries), sep=''))	
      }
      if(!(inherits(tryResult, "try-error"))) {
        upstreamSeqs <- tryResult
        for(i in 1:1000) {
          #if the upstream seq is incomplete or the gene_id is not already in the df skip it
          if(nchar(upstreamSeqs[i,1])<5000 | !(upstreamSeqs[i,2] %in% seq_df$gene_id)) next()
          #looks for a random sequence of length==5utr that is 8x 5utr lengths away + random length in range of 5utr length
          thisUtr <- seq_df$utr[seq_df$gene_id==upstreamSeqs[i,2]]
          randSeqStart <- nchar(thisUtr)*8+sample(nchar(thisUtr),1)
          randSeqEnd <- randSeqStart+nchar(thisUtr)
          randSeq <- substr(upstreamSeqs[1,"coding_gene_flank"], randSeqStart, randSeqEnd)
          seq_df$upstream[seq_df$gene_id == upstreamSeqs[i,2]] <- randSeq
        }
      }
    }
  }
  print(utrBlock+1)
  tryResult <- try(biomaRt::getSequence(
    id =geneIds[(utrBlock*1000+1):(utrBlock*1000+utrRem)],
    type='ensembl_gene_id',
    seqType='coding_gene_flank',
    upstream=5000,
    mart=currSpeciesMart
  ))
  numTries <- 0
  while(inherits(tryResult, "try-error") & numTries<0) {
    tryResult <- try(biomaRt::getSequence(id =geneIds[(utrBlock*1000+1):(utrBlock*1000+utrRem)],type='ensembl_gene_id',seqType='coding_gene_flank',upstream=5000,mart=currSpeciesMart))
    numTries <- numTries+1
    print(paste("Number tries ", as.character(numTries), sep=''))	
  }
  if(!(inherits(tryResult, "try-error"))) {
    remUpstreamSeqs <- tryResult
    if(length(remUpstreamSeqs[,1])!=0) {
      for(i in 1:length(remUpstreamSeqs[,1])) {
        if(nchar(remUpstreamSeqs[i,1])<5000 | !(remUpstreamSeqs[i,2] %in% seq_df$gene_id)) next() 
        thisUtr <- seq_df$utr[seq_df$gene_id==remUpstreamSeqs[i,2]]
        randSeqStart <- nchar(thisUtr)*8+sample(nchar(thisUtr),1)
        randSeqEnd <- randSeqStart+nchar(thisUtr)
        randSeq <- substr(remUpstreamSeqs[1,"coding_gene_flank"], randSeqStart, randSeqEnd)
        seq_df$upstream[seq_df$gene_id == remUpstreamSeqs[i,2]] <- randSeq
      }
    }
  }
  lengthBefore <- length(seq_df$gene_id)

  print(paste("# seqs before checking upstream:", length(seq_df$gene_id)))
  #take out all utrs and gene_ids that do not have a valid upstream sequence
  seq_df <- seq_df[seq_df$upstream!=FALSE,]
  print(paste("# seqs after checking upstream:", length(seq_df$gene_id)))
  
  #add relevant information to summary data frame
  sum <- c(species, length(currSpeciesGenes[,1]), length(geneIds), lengthBefore, length(seq_df$gene_id))
  summary_df[nrow(summary_df)+1,] <- sum
  write.table(summary_df, file = "summary.csv", sep=",")
  
  #write seqs from the data frame to each species file in fasta forma
  for(i in 1:length(seq_df$gene_id)) {
    fileNameUtr=paste("5utrSeq/",species,"/",seq_df$gene_id[i], "_5utr_",species,".fasta",sep='')
    write.fasta(sequences= s2c(seq_df$utr[i]),names= "", file.out=fileNameUtr, open="w")
    fileNameUp=paste("randomUpstream/",species,"/",seq_df$gene_id[i], "_upstream_",species,".fasta",sep='')
    write.fasta(sequences= s2c(seq_df$utr[i]),names= "", file.out=fileNameUp, open="w")
  }
  print("done")
}
print("5'UTR retrieval complete")


