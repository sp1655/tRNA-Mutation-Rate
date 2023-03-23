RNAPIII_mutations_analysis <- function(genomeFastaFile, fChrOK, fPos, fMut){
  
  # Variable to be returned at end of function
  vRes = c(nbBS_observed = 0, nbBS_expected_GC = 0, nbBS_expected_noGC = 0, basesCovered = 0, 
           AT_GC_intr = 0, GC_AT_intr = 0, AT_TA_intr = 0, AT_CG_intr = 0, GC_TA_intr = 0, GC_CG_intr = 0,
           AT_GC_rest = 0, GC_AT_rest = 0, AT_TA_rest = 0, AT_CG_rest = 0, GC_TA_rest = 0, GC_CG_rest = 0,
           AT_intr = 0, AG_intr = 0, AC_intr = 0, TA_intr = 0, TG_intr = 0, TC_intr = 0, 
           GA_intr = 0, GT_intr = 0, GC_intr = 0, CA_intr = 0, CT_intr = 0, CG_intr = 0,
           AT_rest = 0, AG_rest = 0, AC_rest = 0, TA_rest = 0, TG_rest = 0, TC_rest = 0, 
           GA_rest = 0, GT_rest = 0, GC_rest = 0, CA_rest = 0, CT_rest = 0, CG_rest = 0,
           aFreq = 0, tFreq = 0, gFreq = 0, cFreq = 0, aFreqRest = 0, tFreqRest = 0, gFreqRest = 0, cFreqRest = 0)
  
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
  
  # Extracting the sequence from RNAPIII-transcribed regions
  seqIII = getSeq(genome, gpIII)
  # Computing the GC content
  alfIII = alphabetFrequency(seqIII, collapse = T)
  gcIII = sum(alfIII[c("G", "C")]) / sum(alfIII[c("A", "C", "G", "T")])
  gcIII
  
  gRest = setdiff(gChrOK, gpIII, ignore.strand = T)
  
  cat("GC content for intr genes: ", gcIII, "\n", sep = "")
  
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
  
  j = 0
  if(length(gmIII) > 0) {
    gmIII = gmIII[which(nchar(gmIII$from) == 1 & nchar(gmIII$to) == 1),]
    
    lsRepl = c("A", "T", "G", "C")
    names(lsRepl) = c("T", "A", "C", "G")
    
    j = 0
    
    # This takes a long time to run - parallelize? would need to write separate function
    # Assigning strand to mutations
    for(i in 1:length(gmIII)) {
      chr = as.character(gmIII@seqnames[i])
      pos = gmIII$gPos[i]
      
      tmp = as.data.frame(gpIII[which(gpIII@seqnames == chr),])
      str = tmp$strand[which(tmp$start <= pos & tmp$end >= pos)]
      
      # Assigning strand of gene match to mutation
      if(length(str) == 1) {
        gmIII@strand[i] = str
      }
      
      # Assigning strand of first match to mutations with multiple gene matches
      else {
        str = str[1]
        gmIII@strand[i] = str
        #cat("Error: Mutation at position ", pos, " in chromosome ", chr, " matches to multiple genes. Recording as first 
        #  gene listed in gmIII.\n", sep = "")
          
        j = j + 1  
      }
    
      # Changing the to and from columns to the appropriate match for -strand mutations
      if(str == "-") {
        gmIII$from[i] = lsRepl[gmIII$from[i]]
        gmIII$to[i] = lsRepl[gmIII$to[i]]
      }
    }
  }
  
  cat("Number of mutations matched to multiple genes: ", j, "\n", sep = "")
  
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
  
  vRes["nbBS_expected_noGC"] = nbBS_expected_III
  vRes["nbBS_observed"] = nbBS_III
  vRes["basesCovered"] = sizeIII
  
  nbBS_III
  cat("Base-substitutions in RNAPIII-transcribed regions: ", nbBS_III, " observed vs ", nbBS_expected_III, " expected\n")
  
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
  
  ###################################################################################################################
  # Finding numbers for each specific transition/transversion
  
  #vNewBases <- c("A", "G")
  vBPS = c("A>G", "T>C", "G>A", "C>T", "A>T", "T>A", "A>C", "T>G", "G>T", "C>A", "G>C", "C>G")
  
  nbBS_III_BPS = rep(0, length(vBPS))
  nbBS_rest_BPS = rep(0, length(vBPS))
  nbBS_expected_III_BPS = rep(0, length(vBPS))
  
  names(nbBS_III_BPS) = vBPS
  names(nbBS_rest_BPS) = vBPS
  names(nbBS_expected_III_BPS) = vBPS
  
  for(newBase in vBases) {
    for(base in vBases) {
      if(base != newBase) {
        #cat(newBase, " to ", base, "\n")
        tmpLabel = paste(newBase, ">", base, sep = "")
        #cat("label: ", tmpLabel, "\n", sep = "")
        
        nbBS_III_BPS[tmpLabel] = length(which(gmIII$type == "BASE_SUB" & gmIII$from == newBase & gmIII$to == base))
        nbBS_rest_BPS[tmpLabel] = length(which(gmRest$type == "BASE_SUB" & gmRest$from == newBase & gmRest$to == base))
        nbBS_expected_III_BPS[tmpLabel] = alfIII[newBase] * (nbBS_rest_BPS[tmpLabel]/alfRest[newBase])
        #cat("BS: ", newBase, " to ", base, " - ", tmp1, " observed vs ", tmpExp, " expected\n", sep = "")
      }
    }
  }
  
  
  vBPS_full = c("A:T>G:C", "G:C>A:T", "A:T>T:A", "A:T>C:G", "G:C>T:A", "G:C>C:G")
  
  nbBS_rest_BPS_full = rep(0, length(vBPS_full))
  nbBS_III_BPS_full = rep(0, length(vBPS_full))
  nbBS_expected_III_BPS_full = rep(0, length(vBPS_full))
  
  names(nbBS_III_BPS_full) = vBPS_full
  names(nbBS_expected_III_BPS_full) = vBPS_full
  
  # Combining appropriate mutation types to one value (full mutation type)
  nbBS_rest_BPS_full["A:T>G:C"] = nbBS_rest_BPS["A>G"] + nbBS_rest_BPS["T>C"]
  nbBS_rest_BPS_full["G:C>A:T"] = nbBS_rest_BPS["G>A"] + nbBS_rest_BPS["C>T"]
  nbBS_rest_BPS_full["A:T>T:A"] = nbBS_rest_BPS["A>T"] + nbBS_rest_BPS["T>A"]
  nbBS_rest_BPS_full["A:T>C:G"] = nbBS_rest_BPS["A>C"] + nbBS_rest_BPS["T>G"]
  nbBS_rest_BPS_full["G:C>T:A"] = nbBS_rest_BPS["G>T"] + nbBS_rest_BPS["C>A"]
  nbBS_rest_BPS_full["G:C>C:G"] = nbBS_rest_BPS["G>C"] + nbBS_rest_BPS["C>G"]
  
  nbBS_III_BPS_full["A:T>G:C"] = nbBS_III_BPS["A>G"] + nbBS_III_BPS["T>C"]
  nbBS_III_BPS_full["G:C>A:T"] = nbBS_III_BPS["G>A"] + nbBS_III_BPS["C>T"]
  nbBS_III_BPS_full["A:T>T:A"] = nbBS_III_BPS["A>T"] + nbBS_III_BPS["T>A"]
  nbBS_III_BPS_full["A:T>C:G"] = nbBS_III_BPS["A>C"] + nbBS_III_BPS["T>G"]
  nbBS_III_BPS_full["G:C>T:A"] = nbBS_III_BPS["G>T"] + nbBS_III_BPS["C>A"]
  nbBS_III_BPS_full["G:C>C:G"] = nbBS_III_BPS["G>C"] + nbBS_III_BPS["C>G"]

  nbBS_expected_III_BPS_full["A:T>G:C"] = nbBS_expected_III_BPS["A>G"] + nbBS_expected_III_BPS["T>C"]
  nbBS_expected_III_BPS_full["G:C>A:T"] = nbBS_expected_III_BPS["G>A"] + nbBS_expected_III_BPS["C>T"]
  nbBS_expected_III_BPS_full["A:T>T:A"] = nbBS_expected_III_BPS["A>T"] + nbBS_expected_III_BPS["T>A"]
  nbBS_expected_III_BPS_full["A:T>C:G"] = nbBS_expected_III_BPS["A>C"] + nbBS_expected_III_BPS["T>G"]
  nbBS_expected_III_BPS_full["G:C>T:A"] = nbBS_expected_III_BPS["G>T"] + nbBS_expected_III_BPS["C>A"]
  nbBS_expected_III_BPS_full["G:C>C:G"] = nbBS_expected_III_BPS["G>C"] + nbBS_expected_III_BPS["C>G"]
  
  # Printing values
  for(bp in vBPS_full) {
    cat(bp, ": ", nbBS_III_BPS_full[bp], " observed vs ", nbBS_expected_III_BPS_full[bp], " expected\n", sep = "")
  }
  
  # Filling in vRes with the appropriate numbers for each mutation type
  vRes["AT_intr"] = nbBS_III_BPS["A>T"]
  vRes["AG_intr"] = nbBS_III_BPS["A>G"]
  vRes["AC_intr"] = nbBS_III_BPS["A>C"]
  vRes["TA_intr"] = nbBS_III_BPS["T>A"]
  vRes["TG_intr"] = nbBS_III_BPS["T>G"]
  vRes["TC_intr"] = nbBS_III_BPS["T>C"]
  vRes["GA_intr"] = nbBS_III_BPS["G>A"]
  vRes["GT_intr"] = nbBS_III_BPS["G>T"]
  vRes["GC_intr"] = nbBS_III_BPS["G>C"]
  vRes["CA_intr"] = nbBS_III_BPS["C>A"]
  vRes["CT_intr"] = nbBS_III_BPS["C>T"]
  vRes["CG_intr"] = nbBS_III_BPS["C>G"]
  
  vRes["AT_rest"] = nbBS_rest_BPS["A>T"]
  vRes["AG_rest"] = nbBS_rest_BPS["A>G"]
  vRes["AC_rest"] = nbBS_rest_BPS["A>C"]
  vRes["TA_rest"] = nbBS_rest_BPS["T>A"]
  vRes["TG_rest"] = nbBS_rest_BPS["T>G"]
  vRes["TC_rest"] = nbBS_rest_BPS["T>C"]
  vRes["GA_rest"] = nbBS_rest_BPS["G>A"]
  vRes["GT_rest"] = nbBS_rest_BPS["G>T"]
  vRes["GC_rest"] = nbBS_rest_BPS["G>C"]
  vRes["CA_rest"] = nbBS_rest_BPS["C>A"]
  vRes["CT_rest"] = nbBS_rest_BPS["C>T"]
  vRes["CG_rest"] = nbBS_rest_BPS["C>G"]
  
  vRes["AT_GC_intr"] = nbBS_III_BPS_full["A:T>G:C"]
  vRes["GC_AT_intr"] = nbBS_III_BPS_full["G:C>A:T"]
  vRes["AT_TA_intr"] = nbBS_III_BPS_full["A:T>T:A"]
  vRes["AT_CG_intr"] = nbBS_III_BPS_full["A:T>C:G"]
  vRes["GC_TA_intr"] = nbBS_III_BPS_full["G:C>T:A"]
  vRes["GC_CG_intr"] = nbBS_III_BPS_full["G:C>C:G"]
  
  vRes["AT_GC_rest"] = nbBS_rest_BPS_full["A:T>G:C"]
  vRes["GC_AT_rest"] = nbBS_rest_BPS_full["G:C>A:T"]
  vRes["AT_TA_rest"] = nbBS_rest_BPS_full["A:T>T:A"]
  vRes["AT_CG_rest"] = nbBS_rest_BPS_full["A:T>C:G"]
  vRes["GC_TA_rest"] = nbBS_rest_BPS_full["G:C>T:A"]
  vRes["GC_CG_rest"] = nbBS_rest_BPS_full["G:C>C:G"]
  
  ##########################################################################################################################################
  
  nbBS_expected_III_corrected = sum(nbBS_expected_III_bases)
  cat("After correcting for base composition: ", sum(nbBS_III_bases), " observed vs ", sum(nbBS_expected_III_bases), " expected.\n")
  
  vRes["nbBS_expected_GC"] = sum(nbBS_expected_III_bases)
  
  ##########################################################################################################################################
  # Getting nucleotide frequencies in genes of interest (intr or protein-coding genes)
  gpIII.tab = as.data.frame(gpIII)
  
  aFreq = 0
  tFreq = 0
  gFreq = 0
  cFreq = 0
  
  for(i in 1:nrow(gpIII.tab)) {
    chr = as.character(gpIII.tab$seqnames[i])
    start = gpIII.tab$start[i]
    end = gpIII.tab$end[i]
    
    tmpseq = subseq(genome[chr], start = start, end = end)
    
    strand = gpIII.tab$strand[i]
    
    if(strand == "+") {
      aFreq = aFreq + letterFrequency(tmpseq, letters = c("A"))
      tFreq = tFreq + letterFrequency(tmpseq, letters = c("T"))
      gFreq = gFreq + letterFrequency(tmpseq, letters = c("G"))
      cFreq = cFreq + letterFrequency(tmpseq, letters = c("C"))
    }
    
    if(strand == "-") {
      aFreq = aFreq + letterFrequency(tmpseq, letters = c("T"))
      tFreq = tFreq + letterFrequency(tmpseq, letters = c("A"))
      gFreq = gFreq + letterFrequency(tmpseq, letters = c("C"))
      cFreq = cFreq + letterFrequency(tmpseq, letters = c("G"))
    }
  }
  
  vRes["aFreq"] = aFreq
  vRes["tFreq"] = tFreq
  vRes["gFreq"] = gFreq
  vRes["cFreq"] = cFreq
  
  # Getting nucleotide frequencies in the rest of the genome
  aFreqRest = sum(letterFrequency(genome, letters = c("A"))) - aFreq
  tFreqRest = sum(letterFrequency(genome, letters = c("T"))) - tFreq
  gFreqRest = sum(letterFrequency(genome, letters = c("G"))) - gFreq
  cFreqRest = sum(letterFrequency(genome, letters = c("C"))) - cFreq
  
  vRes["aFreqRest"] = aFreqRest
  vRes["tFreqRest"] = tFreqRest
  vRes["gFreqRest"] = gFreqRest
  vRes["cFreqRest"] = cFreqRest
  
  # Returning vRes
  vRes
}
