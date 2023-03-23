contentAnalysis <- function(genomeFastaFile, fChrOK, fPos, fMut){
  
  vRes = c(GTfreqPolar = 0, GTfreqNonpolar = 0, minusFreq = 0, nbMut = 0)
  
  # File of chromosomes (sizes) to use in the analysis (= excluding Mito/Chloroplast)
  tcOK = read.table(fChrOK, h=T, as.is=T, sep="\t")
  tcOK$start = 1
  tcOK$end = tcOK$size
  tcOK$strand = "*"
  gChrOK = makeGRangesFromDataFrame(tcOK)
  
  
  # Positions of RNAPIII-transcribed regions
  tp = read.table(fPos, h=T, as.is=T, sep="\t")
  
  tpr = tp[which(is.element(tp$chromosome, tcOK$chromosome)==T) , ]
  gpIII = reduce(makeGRangesFromDataFrame(tp, keep.extra.columns = T))
  
  
  # Loading the genome
  genome = readDNAStringSet(genomeFastaFile)
  
  # Converting chromosome IDs into the short version (For example: II dna:chromosome chromosome:R64-1-1:II:1:813184:1 REF --> II)
  for(i in (1:length(names(genome)))){
    fullID = names(genome)[i]
    chrID = unlist(strsplit(fullID, " "))[1]
    names(genome)[i] = chrID
  }
  
  # Runs only for B. subtilis
  if(length(names(genome)) > 1 & names(genome)[2] == "bacillus_subtillis_3610_plasmid") {
    names(genome)[1] = "seq71907"
    
    tcOK$chromosome[1] = "seq71907"
    gChrOK = makeGRangesFromDataFrame(tcOK)
  }
  
  # Fixing genome names for R. toruloides
  if(sp == "Rhodosporidium_toruloides") {
    for(i in 1:length(genome)) {
      tmp.ls = strsplit(names(genome)[i], "\\|")
      tmp.ls = unlist(tmp.ls)
      names(genome)[i] = tmp.ls[2]
    }
  }
  
  # Fixing genome names for A. thaliana
  if(sp == "Arabidopsis_thaliana") {
    names(genome) = c(1:5)
  }
  
  # Working with the mutations
  tm = read.table(fMut, h=T, as.is=T, sep="\t")
  
  # Runs only for B. subtilis
  if(length(names(genome)) > 1 & names(genome)[2] == "bacillus_subtillis_3610_plasmid") {
    tm$chromosome = "seq71907"
  }
  
  tmr = tm[which(is.element(tm$chromosome, tcOK$chromosome)==T) , ]
  tm = tmr
  tm$start = tm$gPos
  tm$end = tm$start
  tm$strand = "*"
  gm = makeGRangesFromDataFrame(tm, keep.extra.columns = TRUE)
  
  co = countOverlaps(gm, gpIII, ignore.strand = T)
  
  gmIII = gm[which(co>0)]
  gmRest = gm[which(co==0)]
  
  # Extracting the sequence from RNAPIII-transcribed regions
  seqIII = getSeq(genome, gpIII)
  # Computing the GC content
  alfIII = alphabetFrequency(seqIII, collapse = T)
  gcIII = sum(alfIII[c("G", "C")]) / sum(alfIII[c("A", "C", "G", "T")])
  gcIII
  
  gRest = setdiff(gChrOK, gpIII, ignore.strand = T)
  
  cat("GC content for tRNA genes: ", gcIII, "\n", sep = "")
  
  # Computing the GC content for the rest of the genome:
  seqRest = getSeq(genome, gRest)
  
  #vChr = unique(seqnames(gRest))
  #for(chr in vChr){
  #  seqLength = lengths(genome[chr])
  #  gg = gRest[which(seqnames(gRest)==chr)]
  #  ends = end(gg)
  #  starts = start(gg)
  #  maxPos = max(ends, starts)
  #  cat("Testing if the probem is with chromosome: ", chr, " (seqLength:", seqLength, " vs maxPos:", maxPos, ") problem?", maxPos>=seqLength, "\n")
  #}
  
  alfRest = alphabetFrequency(seqRest, collapse = T)
  gcRest = sum(alfRest[c("G", "C")]) / sum(alfRest[c("A", "C", "G", "T")])
  gcRest
  
  cat("GC content for rest of genome: ", gcRest, "\n", sep = "")
  
  tm = read.table(fMut, h=T, as.is=T, sep="\t")
  
  # Runs only for B. subtilis
  if(length(names(genome)) > 1 & names(genome)[2] == "bacillus_subtillis_3610_plasmid") {
    tm$chromosome = "seq71907"
  }
  
  tmr = tm[which(is.element(tm$chromosome, tcOK$chromosome)==T) , ]
  tm = tmr
  tm$start = tm$gPos
  tm$end = tm$start
  tm$strand = "*"
  gm = makeGRangesFromDataFrame(tm, keep.extra.columns = TRUE)
  
  co = countOverlaps(gm, gpIII, ignore.strand = T)
  
  gmIII = gm[which(co>0)]
  gmRest = gm[which(co==0)]
  
  #############################################################
  # Simple expectation (= without correcting for GC content)
  
  sizeOK = sum(lengths(gChrOK))
  sizeRest = sum(lengths(gRest))
  sizeIII = sum(lengths(gpIII))
  
  # Working only with base-substitutions for now:
  nbBS_total = length(which(tm$type=="BASE_SUB"))
  nbBS_rest = length(which(gmRest$type=="BASE_SUB"))
  nbBS_III = length(which(gmIII$type=="BASE_SUB"))
  
  nbBS_expected_III = sizeIII * (nbBS_rest / sizeRest)
  nbBS_expected_III
  
  #vRes["nbBS_expected_noGC"] = nbBS_expected_III
  #vRes["nbBS_observed"] = nbBS_III
  #vRes["basesCovered"] = sizeIII
  
  nbBS_III
  cat("Base-substitutions in RNAPIII-transcribed regions: ", nbBS_III, " observed vs ", nbBS_expected_III, " expected\n")
  
  vRes["nbMut"] = nbBS_III
  
  #############################################################
  # More complex version (= correcting for GC content)
  vBases = c("A", "T", "C", "G")
  nbBS_expected_III_total = 0
  
  nbBS_rest_bases = rep(0, length(vBases))
  nbBS_III_bases = rep(0, length(vBases))
  nbBS_expected_III_bases = rep(0, length(vBases))
  names(nbBS_rest_bases) = vBases
  names(nbBS_III_bases) = vBases
  names(nbBS_expected_III_bases) = vBases
  
  
  for(base in vBases){
    nbBS_III_bases[base] = length(which(gmIII$type=="BASE_SUB" & gmIII$from==base))
    nbBS_rest_bases[base] = length(which(gmRest$type=="BASE_SUB" & gmRest$from==base))
    nbBS_expected_III_bases[base] = alfIII[base] * (nbBS_rest_bases[base] / alfRest[base])
    cat(base, ": ", nbBS_III_bases[base], " observed vs ", nbBS_expected_III_bases[base], " expected.\n")
  }
  
  
  # Checking errors to make sure they're actually there
  if(length(gmIII) >= 1) {
    gmIII.tab = as.data.frame(gmIII)
    gmIII.tab = gmIII.tab[,c("seqnames", "gPos", "from", "to")]
    gmIII.tab$strand = "*"
    gmIII.tab$strand = factor(gmIII.tab$strand, levels = c("+", "-", "*"))
    
    # Getting only the SNPs
    gmIII.tab = gmIII.tab[which(nchar(gmIII.tab$from) == 1 & nchar(gmIII.tab$to) == 1),]
    
    gpIII.tab = as.data.frame(gpIII)
    
    # Getting the strand for each mutation
    for(i in 1:nrow(gmIII.tab)) {
      chr = gmIII.tab$seqnames[i]
      pos = gmIII.tab$gPos[i]
      
      tmp.tab = gpIII.tab[which(gpIII.tab$seqnames == chr),]
      str = tmp.tab$strand[which(tmp.tab$start <= pos & tmp.tab$end >= pos)]
      
      gmIII.tab$strand[i] = str
    }
    
    gmIII.tab$refBase = "X"
    gmIII.tab$mutBase = "X"
    for(i in 1:nrow(gmIII.tab)) {
      if(gmIII.tab$strand[i] == "+") {
        gmIII.tab$refBase[i] = gmIII.tab$from[i]
        gmIII.tab$mutBase[i] = gmIII.tab$to[i]
      }
      
      else {
        tmpFrom = gmIII.tab$from[i]
        tmpTo = gmIII.tab$to[i]
        
        if(tmpFrom == "A") {
          gmIII.tab$refBase[i] = "T"
        }
        if(tmpTo == "A") {
          gmIII.tab$mutBase[i] = "T"
        }
        
        if(tmpFrom == "T") {
          gmIII.tab$refBase[i] = "A"
        }
        if(tmpTo == "T") {
          gmIII.tab$mutBase[i] = "A"
        }
        
        if(tmpFrom == "G") {
          gmIII.tab$refBase[i] = "C"
        }
        if(tmpTo == "G") {
          gmIII.tab$mutBase[i] = "C"
        }
      
        if(tmpFrom == "C") {
          gmIII.tab$refBase[i] = "G"
        }
        if(tmpTo == "C") {
          gmIII.tab$mutBase[i] = "G"
        }
      }
    }
  
    # Number of G->T errors on the + and - strands
    nbGTPolar = nrow(gmIII.tab[which(gmIII.tab$refBase == "G" & gmIII.tab$mutBase == "T"),])
    nbGTNonpolar = nrow(gmIII.tab[which(gmIII.tab$from == "G" & gmIII.tab$to == "T"),])
    
    # Number of mutations that were found on the - strand
    nbMinus = nrow(gmIII.tab[which(gmIII.tab$strand == "-"),])
    
    # Getting frequencies
    GTfreqPolar = nbGTPolar/nrow(gmIII.tab)
    GTfreqNonpolar = nbGTNonpolar/nrow(gmIII.tab)
    minusFreq = nbMinus/nrow(gmIII.tab)
    
    # Printing frequencies
    cat("G->T frequency in tRNAs (polarized): ", GTfreqPolar, '\n')
    cat("G->T frequency in tRNAs (not polarized): ", GTfreqNonpolar, '\n')
    cat("-strand frequency in tRNAs: ", minusFreq, '\n')
    
    # Storing frequencies in vRes
    vRes["GTfreqNonpolar"] = GTfreqNonpolar
    vRes["GTfreqPolar"] = GTfreqPolar
    vRes["minusFreq"] = minusFreq
  }
  
  # Example for how to look at one specific base in the reference genome
  # subseq(genome, start = 1, end = 1)
  
  # Returning vRes from function
  vRes
}
