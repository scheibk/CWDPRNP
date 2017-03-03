############PRNP Functions##############
####Checks install packages####
CheckCRANPackages <- function(packages, repos="http://cran.r-project.org",
                                depend=c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")) {
  for(p in packages) {
    pckg=try( require(p, character.only=TRUE) )
    if(!pckg) {
      install.packages(p, repos=repos, dependencies=depend)
    }
    require(p, character.only=TRUE) 
  }
} 

CheckBioconductorPackages <- function(packages, repos="https://bioconductor.org/biocLite.R",
                              depend=c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")) {
  for(p in packages) {
    pckg=try( require(p, character.only=TRUE) )
    if(!pckg) {
      biocLite(p, repos=repos, dependencies=depend)
    }
    require(p, character.only=TRUE) 
  }
} 
####Necessary Background Functions####
list.DNAStringSet<-function(x){
  DNAStringSet(unlist(x, use.names=T),use.names=T)
}

list.contig<-function(x){
  DNAString(as.character(ConsensusSequence(x, noConsensusChar = "S")))
}

align.method<-function(x, method){
  if(method=="DECIPHER"){
    align<-DECIPHER::AlignSeqs(x)
    return(align)
  }
  if(method=="ClustalW"){
    align<-msa(x, method=c("ClustalW"))
    return(align)
  }
  if(method=="ClustalOmega"){
    align<-msa(x, method=c("ClustalOmega"))
    return(align)
  }
  if(method=="Muscle"){
    align<-msa(x, method=c("Muscle"))
  return(align)
  }
}
########

####Load Reference Sequence####
RefSeq<-function(x){
  if(x=="CECA"){    
    elk<-read.GenBank("AF016228")
    write.dna(elk, file="CECAprion.fasta", format="fasta", nbcol=10, colsep="")
    elk.prion<-readDNAStringSet("CECAprion.fasta")
    elk.prion<-unlist(elk.prion)
    elk.prion<-as.character(elk.prion)
    elk.prion<-paste(elk.prion, sep="", collapse="")
    elk.prion<-DNAString(elk.prion)
    return(elk.prion)}
  
  if(x=="ODVI"){
    wtd<-read.GenBank("AY275712")
    write.dna(wtd, file="ODVIprion.fasta", format="fasta",nbcol=10, colsep="")
    wtd.prion<-readDNAStringSet("ODVIprion.fasta")
    wtd.prion<-unlist(wtd.prion)
    wtd.prion<-as.character(wtd.prion)
    wtd.prion<-paste(wtd.prion, sep="", collapse="")
    wtd.prion<-DNAString(wtd.prion)
    return(wtd.prion)
  }
  if(x=="ODHE"){
    mule<-read.GenBank("U25965")
    write.dna(mule, file="ODHEprion.fasta", format="fasta",nbcol=10, colsep="")
    mule.prion<-readDNAStringSet("ODHEprion.fasta")
    mule.prion<-unlist(mule.prion)
    mule.prion<-as.character(mule.prion)
    mule.prion<-paste(mule.prion, sep="", collapse="")
    mule.prion<-DNAString(mule.prion)
    return(mule.prion)
  }
  if(x=="RATA"){
    rein<-read.GenBank("AY639093")
    write.dna(rein, file="RATAprion.fasta", format="fasta", nbcol=10, colsep="")
    rein.prion<-readDNAStringSet("RATAprion.fasta")
    rein.prion<-unlist(rein.prion)
    rein.prion<-as.character(rein.prion)
    rein.prion<-paste(rein.prion, sep="", collapse="")
    rein.prion<-DNAString(rein.prion)
    return(rein.prion)
  }
  if(x==y){
    ref<-read.GenBank(y)
    write.dna(ref, file="refprion.fasta", format="fasta",nbcol=10, colsep="")
    ref.prion<-readDNAStringSet("refprion.fasta")
    ref.prion<-unlist(ref.prion)
    ref.prion<-as.character(ref.prion)
    ref.prion<-paste(ref.prion, sep="", collapse="")
    ref.prion<-DNAString(ref.prion)
    return(ref.prion)
  }
}

####

#########Call Bases
CallBases<-function(flab, rlab, label){
  temp<-list.files(pattern=flab)
  temp2<-list.files(pattern=rlab)
  sub<-gsub(label, "", temp2)
  names(temp)<-names(temp2)<-sub
  forward<-lapply(temp, readsangerseq)
  forward<-lapply(forward, makeBaseCalls)
  reverse<-lapply(temp2, readsangerseq)
  reverse<-lapply(reverse, makeBaseCalls)
  list<-list(forward,reverse)
  return(list)
}
####

###############Write Chromatograms
Chrom<-function(forward,reverse){
  key<-names(forward)
  pdf1<-paste(key, "F.pdf", sep="")
  for(i in 1:length(forward)){
    chromatogram(forward[[i]],showcalls = c("both"), filename=pdf1[[i]], showhets=T)
  }
  pdf2<-paste(key, "R.pdf", sep="")
  for(i in 1:length(reverse)){
    chromatogram(reverse[[i]], showcalls = c("both"), filename=pdf2[[i]], showhets=T)
  }
}

##################


####Alignments####
Align<-function(ref, forward, reverse, method){
  ref<-cbind(c("Ref"), as.character(ref))
  key<-names(forward)
  key1<-c(key, "Ref")
  
  f1<-list()
  for(i in 1:length(forward)){
    f1[[i]]<-as.character(unlist(forward[[i]]@primarySeq))
  }
  f1<-cbind(key,as.matrix(unlist(f1)))
  colnames(f1)<-c("ID", "Sequence")
  f1<-rbind(f1, as.character(ref))
  f1prim<-list()
  for(i in 1:length(f1[,2])){
    f1prim[[i]]<-DNAString(f1[i,2])
  }
  names(f1prim)<-f1[,1]
  f1align<-list.DNAStringSet(f1prim)
  f1align<-align.method(f1align, method)
  f1trim<-as.matrix(f1align)
  f1trim<-f1trim[,-which(f1trim["Ref",]=="-")]
  f1trim<-as.data.frame(f1trim)
  f1trim<-f1trim[match(key1, rownames(f1trim)),]
  f1trim<-data.frame(lapply(f1trim, as.character), stringsAsFactors=F)
  
  f2<-list()
  for(i in 1:length(forward)){
    f2[[i]]<-as.character(unlist(forward[[i]]@secondarySeq))
  }
  f2<-cbind(key, as.matrix(unlist(f2)))
  colnames(f2)<-c("ID","Sequence")
  f2<-rbind(f2, as.character(ref))
  f2prim<-list()
  for(i in 1:length(f2[,2])){
    f2prim[[i]]<-DNAString(f2[i,2])
  }
  names(f2prim)<-f2[,1]
  f2align<-list.DNAStringSet(f2prim)
  f2align<-align.method(f2align, method)
  f2trim<-as.matrix(f2align)
  f2trim<-f2trim[,-which(f2trim["Ref",]=="-")]
  f2trim<-as.data.frame(f2trim)
  f2trim<-f2trim[match(key1, rownames(f2trim)),]
  f2trim<-data.frame(lapply(f2trim, as.character), stringsAsFactors=F)
  
  comp1<-list()
  for(i in 1:length(reverse)){
    comp1[[i]]<-reverse[[i]]@primarySeq
  }
  comp1<-DNAStringSet(comp1)
  comp1<-reverseComplement(comp1)
  comp1<-cbind(key, as.character(comp1))
  r1<-rbind(comp1, as.character(ref))
  colnames(r1)<-c("ID","Sequence")
  r1prim<-list()
  for(i in 1:length(r1[,2])){
    r1prim[[i]]<-DNAString(r1[i,2])
  }
  names(r1prim)<-r1[,1]
  r1align<-list.DNAStringSet(r1prim)
  r1align<-align.method(r1align, method)
  r1trim<-as.matrix(r1align)
  r1trim<-r1trim[,-which(r1trim["Ref",]=="-")]
  r1trim<-as.data.frame(r1trim)
  r1trim<-r1trim[match(key1, rownames(r1trim)),]
  r1trim<-data.frame(lapply(r1trim, as.character), stringsAsFactors=F)
  
  comp2<-list()
  for(i in 1:length(reverse)){
    comp2[[i]]<-reverse[[i]]@secondarySeq
  }
  comp2<-DNAStringSet(comp2)
  comp2<-reverseComplement(comp2)
  comp2<-cbind(key, as.character(comp2))
  r2<-rbind(comp2, as.character(ref))
  colnames(r2)<-c("ID","Sequence")
  r2prim<-list()
  for(i in 1:length(r2[,2])){
    r2prim[[i]]<-DNAString(r2[i,2])
  }
  names(r2prim)<-r2[,1]
  r2align<-list.DNAStringSet(r2prim)
  r2align<-align.method(r2align, method)
  r2trim<-as.matrix(r2align)
  r2trim<-r2trim[,-which(r2trim["Ref",]=="-")]
  r2trim<-as.data.frame(r2trim)
  r2trim<-r2trim[match(key1, rownames(r2trim)),]
  r2trim<-data.frame(lapply(r2trim, as.character), stringsAsFactors=F)
  rownames(f1trim)<-rownames(f2trim)<-rownames(r1trim)<-rownames(r2trim)<-key1
  list<-list(f1trim, f2trim, r1trim, r2trim)
  return(list)
}
##################


####Align Raw####
Poly<-function(poly, align, ref){
  if(poly==F){
    fallele1<-align[[1]]
    fallele2<-align[[2]]
    rallele1<-align[[3]]
    rallele2<-align[[4]]
    key<-rownames(fallele1)
    a1lab<-paste(key, "_a1", sep="")
    a2lab<-paste(key, "_a2", sep="")
    rownames(fallele1)<-a1lab
    rownames(fallele2)<-a2lab
    rownames(rallele1)<-a1lab
    rownames(rallele2)<-a2lab
    newf<-rbind(fallele1, fallele2)
    newr<-rbind(rallele1, rallele2)
    newf<-newf[!rownames(newf) %in% c("Ref_a1", "Ref_a2"),]
    newr<-newr[!rownames(newr) %in% c("Ref_a1", "Ref_a2"),]
    newf<-newf[order(row.names(newf)),]
    newr<-newr[order(row.names(newr)),]
    ref1<-t(as.matrix(strsplit(as.character(ref), split="")[[1]]))
    rownames(ref1)<-c("Ref")
    newf<-rbind(newf, ref1)
    newr<-rbind(newr, ref1)
    colnames(newf)<-colnames(newr)<-c(1:ncol(newr))
    new<-list(newf, newr)
    return(new)
  }
  if(poly==T){
    fallele1<-align[[1]]
    fallele2<-align[[2]]
    rallele1<-align[[3]]
    rallele2<-align[[4]]
    key<-rownames(fallele1)
    a1lab<-paste(key, "_a1", sep="")
    a2lab<-paste(key, "_a2", sep="")
    rownames(fallele1)<-a1lab
    rownames(fallele2)<-a2lab
    rownames(rallele1)<-a1lab
    rownames(rallele2)<-a2lab
    newf<-rbind(fallele1, fallele2)
    newr<-rbind(rallele1, rallele2)
    newf<-newf[!rownames(newf) %in% c("Ref_a1", "Ref_a2"),]
    newr<-newr[!rownames(newr) %in% c("Ref_a1", "Ref_a2"),]
    newf<-newf[order(row.names(newf)),]
    newr<-newr[order(row.names(newr)),]
    ref1<-t(as.matrix(strsplit(as.character(ref), split="")[[1]]))
    rownames(ref1)<-c("Ref")
    newf<-rbind(newf, ref1)
    newr<-rbind(newr, ref1)
    colnames(newf)<-colnames(newr)<-c(1:ncol(newr))
    new1<-matrix(0,nrow=1,ncol=ncol(newf))
    for(i in 1:ncol(newf)){
      new1[i]<-length(unique(newf[,i]))
    }
    new2<-matrix(0,nrow=1,ncol=ncol(newr))
    for(i in 1:ncol(newr)){
      new2[i]<-length(unique(newr[,i]))
    }
    rownames(new1)<-rownames(new2)<-c("Poly")
    colnames(new1)<-colnames(new2)<-colnames(newf)
    newf<-rbind(newf,new1)
    newr<-rbind(newr,new2)
    new<-list(newf,newr)
    return(new)
  }
}
######################################

####ExtractSNPs####
ExtractSNPs<-function(align){
  alignf<-align[[1]]
  alignr<-align[[2]]
  alignf<-as.data.frame(t(alignf))
  Polyf<-as.numeric(alignf$Poly)
  alignf<-alignf[,-c(ncol(alignf))]
  alignf$Poly<-Polyf
  alignf<-alignf[alignf$Poly>1,]
  alignf<-t(alignf)
  alignr<-as.data.frame(t(alignr))
  Polyr<-as.numeric(alignr$Poly)
  alignr<-alignr[,-c(ncol(alignr))]
  alignr$Poly<-Polyr
  alignr<-alignr[alignr$Poly>1,]
  alignr<-t(alignr)
  polyf<-as.numeric(colnames(alignf))
  polyr<-as.numeric(colnames(alignr))
  polyf1<-which(polyf %in% polyr)
  polyr1<-which(polyr %in% polyf)
  polyf<-alignf[,polyf1]
  polyr<-alignr[,polyr1]
  ref<-t(as.matrix(polyr["Ref",]))
  rownames(ref)<-c("Ref")
  polyf<-polyf[!rownames(polyf) %in% c("Ref", "Poly"),]
  polyr<-polyr[!rownames(polyr) %in% c("Ref", "Poly"),]
  a1lab<-paste(rownames(polyf), "F", sep="")
  a2lab<-paste(rownames(polyr), "R", sep="")
  key<-c(a1lab,a2lab)
  polyfull<-rbind(polyf,polyr)
  rownames(polyfull)<-key
  polyfull<-polyfull[order(rownames(polyfull)),]
  polyfull<-rbind(polyfull,ref)
  SNPs<-list(polyfull,polyf, polyr)
  return(SNPs)
}
################################

####Nucleotide Table####
NucTab<-function(x, align, nuc){
  if(x=="CECA"){
    PRNPlabels<-c(394)
    key<-rownames(align[[1]])
    align2<-list()
    for(i in 1:length(PRNPlabels)){
      align2[[i]]<-cbind(align[[1]][,PRNPlabels[i]], align[[2]][,PRNPlabels[i]], align[[3]][,PRNPlabels[i]], align[[4]][,PRNPlabels[i]])
    }
    nuc394f<-apply(align2[[1]][,c(1,2)], 1, paste, collapse="/")
    nuc394r<-apply(align2[[1]][,c(3,4)], 1, paste, collapse="/")
    nuc394geno<-matrix(0, nrow=NROW(nuc394f), ncol=1, dimnames=list(c(rownames(align2[[1]])), c("Genotype")))
    nuc394geno[,1]<-ifelse((nuc394f=="A/A")&(nuc394r=="A/A"),c("A/A"),nuc394geno)
    nuc394geno[,1]<-ifelse(((nuc394f=="A/T")|(nuc394f=="T/A"))&((nuc394r=="A/T")|(nuc394r=="T/A")), c("A/T"),nuc394geno)
    nuc394geno[,1]<-ifelse((nuc394f=="T/T")&(nuc394r=="T/T"), c("T/T"), nuc394geno)
    nuc394geno[,1]<-ifelse(nuc394geno==0, c("ERROR"), nuc394geno)
    colnames(nuc394geno)<-c("nuc394")
    rownames(nuc394geno)<-key
    return(nuc394geno)
  }
  if(x=="ODVI"){
    PRNPlabels<-c(285, 286, 347)
    key<-rownames(align[[1]])
    align2<-list()
    for(i in 1:length(PRNPlabels)){
      align2[[i]]<-cbind(align[[1]][,PRNPlabels[i]], align[[2]][,PRNPlabels[i]], align[[3]][,PRNPlabels[i]], align[[4]][,PRNPlabels[i]])
    }
    nuc285<-align2[[1]]
    nuc285f<-apply(nuc285[,c(1,2)], 1, paste, collapse="/")
    nuc285r<-apply(nuc285[,c(3,4)], 1, paste, collapse="/")
    nuc285geno<-matrix(0, nrow=NROW(nuc285f), ncol=1, dimnames=list(c(rownames(nuc285)), c("Genotype")))
    nuc285geno[,1]<-ifelse((nuc285f=="A/A")&(nuc285r=="A/A"),c("A/A"), nuc285geno)
    nuc285geno[,1]<-ifelse(((nuc285f=="A/C")|(nuc285f=="C/A"))&((nuc285r=="A/C")|(nuc285r=="C/A")), c("A/C"),nuc285geno)
    nuc285geno[,1]<-ifelse((nuc285f=="C/C")&(nuc285r=="C/C"), c("C/C"), nuc285geno)
    nuc285geno[,1]<-ifelse(nuc285geno==0, c("ERROR"), nuc285geno)
    nuc286<-align2[[2]]
    nuc286f<-apply(nuc286[,c(1,2)], 1, paste, collapse="/")
    nuc286r<-apply(nuc286[,c(3,4)], 1, paste, collapse="/")
    nuc286geno<-matrix(0, nrow=NROW(nuc286f), ncol=1, dimnames=list(c(rownames(nuc286)), c("Genotype")))
    nuc286geno[,1]<-ifelse((nuc286f=="G/G")&(nuc286r=="G/G"),c("G/G"), nuc286geno)
    nuc286geno[,1]<-ifelse(((nuc286f=="G/A")|(nuc286f=="A/G"))&((nuc286r=="G/A")|(nuc286r=="A/G")), c("G/A"),nuc286geno)
    nuc286geno[,1]<-ifelse((nuc286f=="A/A")&(nuc286r=="A/A"), c("A/A"), nuc286geno)
    nuc286geno[,1]<-ifelse(nuc286geno==0, c("ERROR"), nuc286geno)
    nuc347<-align2[[3]]
    nuc347f<-apply(nuc347[,c(1,2)], 1, paste, collapse="/")
    nuc347r<-apply(nuc347[,c(3,4)], 1, paste, collapse="/")
    nuc347geno<-matrix(0, nrow=NROW(nuc347f), ncol=1, dimnames=list(c(rownames(nuc347)), c("Genotype")))
    nuc347geno[,1]<-ifelse((nuc347f=="C/C")&(nuc347r=="C/C"),c("C/C"), nuc347geno)
    nuc347geno[,1]<-ifelse(((nuc347f=="C/G")|(nuc347f=="G/C"))&((nuc347r=="C/G")|(nuc347r=="G/C")), c("C/G"),nuc347geno)
    nuc347geno[,1]<-ifelse((nuc347f=="G/G")&(nuc347r=="G/G"), c("G/G"), nuc347geno)
    nuc347geno[,1]<-ifelse(nuc347geno==0, c("ERROR"), nuc347geno)
    ODVIgeno<-cbind(nuc285geno, nuc286geno, nuc347geno)
    colnames(ODVIgeno)<-c("nuc285", "nuc286", "nuc347")
    rownames(ODVIgeno)<-key
    return(ODVIgeno)
  }
  if(x=="ODHE"){
    PRNPlabels<-c(674)
    key<-rownames(align[[1]])
    align2<-list()
    for(i in 1:length(PRNPlabels)){
      align2[[i]]<-cbind(align[[1]][,PRNPlabels[i]], align[[2]][,PRNPlabels[i]], align[[3]][,PRNPlabels[i]], align[[4]][,PRNPlabels[i]])
    }
    nuc674f<-apply(align2[[1]][,c(1,2)], 1, paste, collapse="/")
    nuc674r<-apply(align2[[1]][,c(3,4)], 1, paste, collapse="/")
    nuc674geno<-matrix(0, nrow=NROW(nuc674f), ncol=1, dimnames=list(c(rownames(align2[[1]])), c("Genotype")))
    nuc674geno[,1]<-ifelse((nuc674f=="C/C")&(nuc674r=="C/C"),c("C/C"),nuc674geno)
    nuc674geno[,1]<-ifelse(((nuc674f=="C/T")|(nuc674f=="T/C"))&((nuc674r=="C/T")|(nuc674r=="T/C")), c("C/T"),nuc674geno)
    nuc674geno[,1]<-ifelse((nuc674f=="T/T")&(nuc674r=="T/T"), c("T/T"), nuc674geno)
    nuc674geno[,1]<-ifelse(nuc674geno==0, c("ERROR"), nuc674geno)
    colnames(nuc674geno)<-c("nuc674")
    rownames(nuc674geno)<-key
    return(nuc674geno)
  }
  if(x=="other"){
    PRNPlabels<-nuc
    key<-rownames(align[[1]])
    align2<-list()
    for(i in 1:length(PRNPlabels)){
      align2[[i]]<-cbind(align[[1]][,PRNPlabels[i]], align[[2]][,PRNPlabels[i]], align[[3]][,PRNPlabels[i]], align[[4]][,PRNPlabels[i]])
    }
    
    for(i in 1:length(align2)){
      align2[[i]]<-ifelse(align2[[i]]!="A"&align2[[i]]!="G"&align2[[i]]!="C"&align2[[i]]!="T","-",align2[[i]])
    }
    
    NucMat<-c(1,5,6,7,"-",5,2,8,9,"-",6,8,3,10,"-",7,9,10,4,"-","-","-","-","-","-")
    NucMat<-matrix(NucMat, nrow=5, ncol=5)
    rownames(NucMat)<-colnames(NucMat)<-c("A","T","C","G","-")
    GenList<-list()
    for(i in 1:length(align2)){
      GenList[[i]]<-matrix(0, nrow=NROW(align2[[1]]), ncol=2, dimnames=list(c(rownames(align2)), c("GenF","GenR")))  
    }
    error<-list()
    for(i in 1:length(align2)){
      error[[i]]<-matrix(0, nrow=NROW(align2[[1]]), ncol=1)
    }
    for(i in 1:length(GenList)){
      for(j in 1:nrow(GenList[[1]])){
        GenList[[i]][j,1]<-NucMat[align2[[i]][j,1],align2[[i]][j,2]]
        GenList[[i]][j,2]<-NucMat[align2[[i]][j,3],align2[[i]][j,4]]
        error[[i]][j,1]<-ifelse(GenList[[i]][j,1]==GenList[[i]][j,2],GenList[[i]][j,1], "Error")
        error[[i]][j,1]<-ifelse(GenList[[i]][j,1]=="-","Error",error[[i]][j,1])
        GenList[[i]]<-cbind(GenList[[i]],error[[i]])
        GenList[[i]]<-cbind(GenList[[i]][,c(1,2)],GenList[[i]][,ncol(GenList[[i]])])
      }
    }
    list<-list()
    for(i in 1:length(align2)){
      list[[i]]<-matrix("-", nrow=NROW(align2[[i]]), ncol=1, dimnames=list(c(rownames(align2[[i]])), c("Genotype")))
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="1", "A/A", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="2", "T/T", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="3", "C/C", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="4", "G/G", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="5", "A/T", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="6", "A/C", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="7", "A/G", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="8", "T/C", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="9", "T/G", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="10", "C/G", list[[i]][,1])
      list[[i]][,1]<-ifelse(GenList[[i]][,3]=="Error", "Error", list[[i]][,1])
    }
    names(list)<-PRNPlabels
    NucTab<-matrix(0, nrow=NROW(list[[1]]), ncol=length(list), dimnames=list(c(rownames(list[[1]])), c(PRNPlabels)))
    for(i in 1:length(list)){
      for(j in 1:nrow(list[[1]])){
        NucTab[j,i]<-list[[i]][j,1]
      }
    }
    rownames(NucTab)<-key
    return(NucTab)
  }
}


####################################

####Codon Table########
CodTab<-function(x, Tab, align){
  if(x=="CECA"){
    NucTab<-Tab
    CECAcodon<-matrix(0, nrow=NROW(NucTab), ncol=1, dimnames=list(c(rownames(NucTab)), c("cod132")))
    CECAcodon[,1]<-ifelse(NucTab[,c("nuc394")]=="A/A", "M/M", CECAcodon[,1])
    CECAcodon[,1]<-ifelse(NucTab[,c("nuc394")]=="A/T", "M/L", CECAcodon[,1])
    CECAcodon[,1]<-ifelse(NucTab[,c("nuc394")]=="T/T", "L/L", CECAcodon[,1])
    CECAcodon[,1]<-ifelse(NucTab[,c("nuc394")]=="ERROR", "ERROR", CECAcodon[,1])
    return(CECAcodon)
  }
  if(x=="ODVI"){
    NucTab<-Tab
    ODVIcodon<-matrix(0, nrow=NROW(NucTab), ncol=3, dimnames=list(c(rownames(NucTab)), c("cod95", "cod96", "cod116")))
    ODVIcodon[,1]<-ifelse(NucTab[,c("nuc285")]=="A/A", "Q/Q", ODVIcodon[,1])
    ODVIcodon[,1]<-ifelse(NucTab[,c("nuc285")]=="A/C", "Q/H", ODVIcodon[,1])
    ODVIcodon[,1]<-ifelse(NucTab[,c("nuc285")]=="C/C", "H/H", ODVIcodon[,1])
    ODVIcodon[,1]<-ifelse(NucTab[,c("nuc285")]=="ERROR", "ERROR", ODVIcodon[,1])
    ODVIcodon[,2]<-ifelse(NucTab[,c("nuc286")]=="G/G", "G/G", ODVIcodon[,2])
    ODVIcodon[,2]<-ifelse(NucTab[,c("nuc286")]=="G/A", "G/S", ODVIcodon[,2])
    ODVIcodon[,2]<-ifelse(NucTab[,c("nuc286")]=="A/A", "S/S", ODVIcodon[,2])
    ODVIcodon[,2]<-ifelse(NucTab[,c("nuc286")]=="ERROR", "ERROR", ODVIcodon[,2])
    ODVIcodon[,3]<-ifelse(NucTab[,c("nuc347")]=="C/C", "A/A", ODVIcodon[,3])
    ODVIcodon[,3]<-ifelse(NucTab[,c("nuc347")]=="C/G", "A/G", ODVIcodon[,3])
    ODVIcodon[,3]<-ifelse(NucTab[,c("nuc347")]=="G/G", "G/G", ODVIcodon[,3])
    ODVIcodon[,3]<-ifelse(NucTab[,c("nuc347")]=="ERROR", "ERROR", ODVIcodon[,3])
    return(ODVIcodon)
  }
  if(x=="ODHE"){
    NucTab<-Tab
    ODHEcodon<-matrix(0, nrow=NROW(NucTab), ncol=1, dimnames=list(c(rownames(NucTab)), c("cod225")))
    ODHEcodon[,1]<-ifelse(NucTab[,c("nuc674")]=="C/C", "S/S", ODHEcodon[,1])
    ODHEcodon[,1]<-ifelse(NucTab[,c("nuc674")]=="C/T", "S/F", ODHEcodon[,1])
    ODHEcodon[,1]<-ifelse(NucTab[,c("nuc674")]=="T/T", "F/F", ODHEcodon[,1])
    ODHEcodon[,1]<-ifelse(NucTab[,c("nuc674")]=="ERROR", "ERROR", ODHEcodon[,1])
    return(ODHEcodon)
  }
  if(x=="other"){
    AATab<-align
    colnames(AATab[[1]])<-colnames(AATab[[2]])<-colnames(AATab[[3]])<-colnames(AATab[[4]])<-c(1:ncol(AATab[[1]]))
    PRNPlabels<-colnames(Tab)
    Nuc<-as.numeric(PRNPlabels)
    Cod<-vector("integer", length(Nuc))
    for(i in 1:length(Cod)){
      Cod[i]<-Nuc[i]/3
    }
    Pos<-vector("integer", length(Cod))
    Pos<-ifelse(Cod-trunc(Cod)>0 & Cod-trunc(Cod)<0.34,1,Pos)
    Pos<-ifelse(Cod-trunc(Cod)>0.34 & Cod-trunc(Cod)<0.67,2,Pos)
    Pos<-ifelse(Cod-trunc(Cod)==0,3, Pos)
    Cod<-ceiling(Cod)
    CodList<-list()
    for(i in 1:length(Nuc)){
      CodList[[i]]<-vector("integer",3)
    }
    for(i in 1:length(Nuc)){
      if(Pos[i]==3){
        CodList[[i]][1]<-Nuc[i]-2
        CodList[[i]][2]<-Nuc[i]-1
        CodList[[i]][3]<-Nuc[i]
      }
      if(Pos[i]==2){
        CodList[[i]][1]<-Nuc[i]-1
        CodList[[i]][2]<-Nuc[i]
        CodList[[i]][3]<-Nuc[i]+1
      }
      if(Pos[i]==1){
        CodList[[i]][1]<-Nuc[i]
        CodList[[i]][2]<-Nuc[i]+1
        CodList[[i]][3]<-Nuc[i]+2
      }
    }
    NucList<-list()
    for(i in 1:length(align)){
      NucList[[i]]<-cbind(align[[i]][,unlist(CodList)])
      colnames(NucList[[i]])<-unlist(CodList)
    }
    Nuc1<-list()
    for(i in 1:length(align)){
      Nuc1[[i]]<-align[[i]][,Nuc]
      Nuc1[[i]]<-as.matrix(Nuc1[[i]])
      colnames(Nuc1[[i]])<-Nuc
    }
    NucList<-list(Cod, Nuc1, NucList)
    
    Key<-NucList[[1]]
    codons<-NucList[[3]]
    for1<-codons[[1]]
    for1<-as.matrix(for1)
    for1<-ifelse(for1!="A"&for1!="G"&for1!="C"&for1!="T","-",for1)
    er1<-for1[for1[,1]=="-",]
    er2<-for1[for1[,2]=="-",]
    er3<-for1[for1[,3]=="-",]
    er<-c(rownames(er1),rownames(er2),rownames(er3))
    er<-unique(er)
    for1<-for1[!rownames(for1)%in%er,]
    DNA1<-list()
    for(i in 1:nrow(for1)){
      DNA1[i]<-paste(for1[i,], collapse="")
    }
    names(DNA1)<-rownames(for1)
    for(i in 1:length(DNA1)){
      DNA1[[i]]<-DNAString(DNA1[[i]])
    }
    AA1<-list()
    for(i in 1:length(DNA1)){
      AA1[[i]]<-Biostrings::translate(DNA1[[i]])
      AA1[[i]]<-as.character(AA1[[i]])
      AA1[[i]]<-strsplit(AA1[[i]], split="")
    }
    AA1<-t(matrix(as.vector(unlist(AA1)),nrow=length(Cod),ncol=length(DNA1), dimnames=list(round(Cod),c(names(DNA1)))))
    
    for2<-codons[[2]]
    for2<-as.matrix(for2)
    for2<-ifelse(for2!="A"&for2!="G"&for2!="C"&for2!="T","-",for2)
    er1<-for2[for2[,1]=="-",]
    er2<-for2[for2[,2]=="-",]
    er3<-for2[for2[,3]=="-",]
    er<-c(rownames(er1),rownames(er2),rownames(er3))
    er<-unique(er)
    for2<-for2[!rownames(for2)%in%er,]
    DNA2<-list()
    for(i in 1:nrow(for2)){
      DNA2[i]<-paste(for2[i,], collapse="")
    }
    names(DNA2)<-rownames(for2)
    for(i in 1:length(DNA2)){
      DNA2[[i]]<-DNAString(DNA2[[i]])
    }
    AA2<-list()
    for(i in 1:length(DNA2)){
      AA2[[i]]<-Biostrings::translate(DNA2[[i]])
      AA2[[i]]<-as.character(AA2[[i]])
      AA2[[i]]<-strsplit(AA2[[i]], split="")
    }
    AA2<-t(matrix(as.vector(unlist(AA2)),nrow=length(Cod),ncol=length(DNA2), dimnames=list(round(Cod),c(names(DNA2)))))
    
    rev1<-codons[[3]]
    rev1<-as.matrix(rev1)
    rev1<-ifelse(rev1!="A"&rev1!="G"&rev1!="C"&rev1!="T","-",rev1)
    er1<-rev1[rev1[,1]=="-",]
    er2<-rev1[rev1[,2]=="-",]
    er3<-rev1[rev1[,3]=="-",]
    er<-c(rownames(er1),rownames(er2),rownames(er3))
    er<-unique(er)
    rev1<-rev1[!rownames(rev1)%in%er,]
    DNA3<-list()
    for(i in 1:nrow(rev1)){
      DNA3[i]<-paste(rev1[i,], collapse="")
    }
    names(DNA3)<-rownames(rev1)
    for(i in 1:length(DNA3)){
      DNA3[[i]]<-DNAString(DNA3[[i]])
    }
    AA3<-list()
    for(i in 1:length(DNA3)){
      AA3[[i]]<-Biostrings::translate(DNA3[[i]])
      AA3[[i]]<-as.character(AA3[[i]])
      AA3[[i]]<-strsplit(AA3[[i]], split="")
    }
    AA3<-t(matrix(as.vector(unlist(AA3)),nrow=length(Cod),ncol=length(DNA3), dimnames=list(round(Cod),c(names(DNA3)))))
    
    rev2<-codons[[4]]
    rev2<-as.matrix(rev2)
    rev2<-ifelse(rev2!="A"&rev2!="G"&rev2!="C"&rev2!="T","-",rev2)
    er1<-rev2[rev2[,1]=="-",]
    er2<-rev2[rev2[,2]=="-",]
    er3<-rev2[rev2[,3]=="-",]
    er<-c(rownames(er1),rownames(er2),rownames(er3))
    er<-unique(er)
    rev2<-rev2[!rownames(rev2)%in%er,]
    DNA4<-list()
    for(i in 1:nrow(rev2)){
      DNA4[i]<-paste(rev2[i,], collapse="")
    }
    names(DNA4)<-rownames(rev2)
    for(i in 1:length(DNA4)){
      DNA4[[i]]<-DNAString(DNA4[[i]])
    }
    AA4<-list()
    for(i in 1:length(DNA4)){
      AA4[[i]]<-Biostrings::translate(DNA4[[i]])
      AA4[[i]]<-as.character(AA4[[i]])
      AA4[[i]]<-strsplit(AA4[[i]], split="")
    }
    AA4<-t(matrix(as.vector(unlist(AA4)),nrow=length(Cod),ncol=length(DNA4), dimnames=list(round(Cod),c(names(DNA4)))))
    AAlist<-list(AA1,AA2,AA3,AA4)
    n1<-rownames(AAlist[[1]])
    n2<-rownames(AAlist[[2]])
    n3<-rownames(AAlist[[3]])
    n4<-rownames(AAlist[[4]])
    n<-n1[n1%in%n2]
    n<-n[n%in%n3]
    n<-n[n%in%n4]
    for(i in 1:length(AAlist)){
      AAlist[[i]]<-AAlist[[i]][n,]
    }
    PRNPlabels<-as.character(Cod)
    AATab<-AAlist
    for(i in 1:length(AATab)){
      AATab[[i]]<-as.matrix(AATab[[i]])
      colnames(AATab[[i]])<-PRNPlabels
    }
    AAlist<-list()
    for(i in 1:length(PRNPlabels)){
      AAlist[[i]]<-cbind(AATab[[1]][,PRNPlabels[i]],AATab[[2]][,PRNPlabels[i]], AATab [[3]][,PRNPlabels[i]], AATab[[4]][,PRNPlabels[i]])
    }
    GenList<-list()  
    for(i in 1:ncol(AATab[[1]])){
      GenList[[i]]<-matrix(0, nrow=NROW(AATab[[1]]), ncol=4)
    }
    for(i in 1:length(GenList)){
      for(j in 1:ncol(GenList[[1]])){
        GenList[[i]][,j]<-AATab[[j]][,i]
      }
    }
    GenList1<-list()
    for(i in 1:length(GenList)){
      GenList1[[i]]<-matrix(0, nrow=NROW(GenList[[1]]), ncol=2)
    }
    for(i in 1:length(GenList1)){
      for(j in 1:nrow(GenList[[1]])){
        GenList1[[i]][j,1]<-paste(min(GenList[[i]][j,1:2]), max(GenList[[i]][j,1:2]), sep="/")
        GenList1[[i]][j,2]<-paste(min(GenList[[i]][j,3:4]), max(GenList[[i]][j,3:4]), sep="/")
      }
    }
    AAlist1<-list()
    for(i in 1:length(GenList1)){
      AAlist1[[i]]<-matrix(0, nrow=NROW(GenList1[[1]]), ncol=1)
    }
    for(i in 1:length(AAlist1)){
      for(j in 1:nrow(GenList1[[1]])){
        AAlist1[[i]][,1]<-ifelse(GenList1[[i]][,1]==GenList1[[i]][,2], GenList1[[i]][,1], AAlist1[[i]])
        AAlist1[[i]][,1]<-ifelse(GenList1[[i]][,1]!=GenList1[[i]][,2], c("Error"), AAlist1[[i]])
      }
    }
    AATab1<-matrix(0, nrow=NROW(GenList[[1]]), ncol=length(AAlist1), dimnames=list(c(rownames(AATab[[1]])),PRNPlabels))
    for(i in 1:nrow(AATab1)){
      for(j in 1:ncol(AATab1)){
        AATab1[i,j]<-AAlist1[[j]][i]
      }
    }
    return(AATab1)
  }
}

##################################

####Prop Tables####
PropTab<-function(Pop, x, CodTab, PopLabels){
  if(Pop==F){
    if(x=="CECA"){
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-c("cod132")
      Cod132<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("M/M", "M/L", "L/L", "Total")))
      Cod132[,c("M/M")]<-sum(CodTab[,c("cod132")]==c("M/M"))
      Cod132[,c("M/L")]<-sum(CodTab[,c("cod132")]==c("M/L"))
      Cod132[,c("L/L")]<-sum(CodTab[,c("cod132")]==c("L/L"))
      Cod132[,c("Total")]<-sum(Cod132)
      prop132<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("M/M", "M/L", "L/L")))
      prop132[,c("M/M")]<-Cod132[,c("M/M")]/Cod132[,c("Total")]
      prop132[,c("M/L")]<-Cod132[,c("M/L")]/Cod132[,c("Total")]
      prop132[,c("L/L")]<-Cod132[,c("L/L")]/Cod132[,c("Total")]
      Allele132<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("M", "L", "Total")))
      Allele132[,c("M")]<-(2*sum(CodTab[,c("cod132")]==c("M/M")))+sum(CodTab[,c("cod132")]==c("M/L"))
      Allele132[,c("L")]<-(2*sum(CodTab[,c("cod132")]==c("L/L")))+sum(CodTab[,c("cod132")]==c("M/L"))
      Allele132[,c("Total")]<-sum(Allele132)
      aprop132<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("M","L")))
      aprop132[,c("M")]<-Allele132[,c("M")]/Allele132[,c("Total")]
      aprop132[,c("L")]<-Allele132[,c("L")]/Allele132[,c("Total")]
      GenProp<-matrix(0, nrow=3,ncol=1, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod132")))
      GenProp[,c("cod132")]<-t(prop132)
      AlleleProp<-matrix(0, nrow=2, ncol=1, dimnames=list(c("Major", "Minor"), c("cod132")))
      AlleleProp[,c("cod132")]<-t(aprop132)
      GenProp<-list(GenProp, AlleleProp, prop132, aprop132)
      names(GenProp)<-c("Genotype","Allele", "Geno132", "Allele132")
      return(GenProp)
    }
    if(x=="ODVI"){
      CodTab<-CodTab[!rownames(CodTab)%in% c("Ref"),]
      Cod95<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("Q/Q", "Q/H", "H/H", "Total")))
      Cod95[,c("Q/Q")]<-sum(CodTab[,c("cod95")]==c("Q/Q"))
      Cod95[,c("Q/H")]<-sum(CodTab[,c("cod95")]==c("Q/H"))
      Cod95[,c("H/H")]<-sum(CodTab[,c("cod95")]==c("H/H"))
      Cod95[,c("Total")]<-sum(Cod95)
      prop95<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("Q/Q", "Q/H", "H/H")))
      prop95[,c("Q/Q")]<-Cod95[,c("Q/Q")]/Cod95[,c("Total")]
      prop95[,c("Q/H")]<-Cod95[,c("Q/H")]/Cod95[,c("Total")]
      prop95[,c("H/H")]<-Cod95[,c("H/H")]/Cod95[,c("Total")]
      Allele95<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("Q", "H","Total")))
      Allele95[,c("Q")]<-(2*sum(CodTab[,c("cod95")]==c("Q/Q")))+sum(CodTab[,c("cod95")]==c("Q/H"))
      Allele95[,c("H")]<-(2*sum(CodTab[,c("cod95")]==c("H/H")))+sum(CodTab[,c("cod95")]==c("Q/H"))
      Allele95[,c("Total")]<-sum(Allele95)
      aprop95<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("Q","H")))
      aprop95[,c("Q")]<-Allele95[,c("Q")]/Allele95[,c("Total")]
      aprop95[,c("H")]<-Allele95[,c("H")]/Allele95[,c("Total")]
      Cod96<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("G/G", "G/S", "S/S", "Total")))
      Cod96[,c("G/G")]<-sum(CodTab[,c("cod96")]==c("G/G"))
      Cod96[,c("G/S")]<-sum(CodTab[,c("cod96")]==c("G/S"))
      Cod96[,c("S/S")]<-sum(CodTab[,c("cod96")]==c("S/S"))
      Cod96[,c("Total")]<-sum(Cod96)
      prop96<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("G/G", "G/S", "S/S")))
      prop96[,c("G/G")]<-Cod96[,c("G/G")]/Cod96[,c("Total")]
      prop96[,c("G/S")]<-Cod96[,c("G/S")]/Cod96[,c("Total")]
      prop96[,c("S/S")]<-Cod96[,c("S/S")]/Cod96[,c("Total")]
      Allele96<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("G", "S","Total")))
      Allele96[,c("G")]<-(2*sum(CodTab[,c("cod96")]==c("G/G")))+sum(CodTab[,c("cod96")]==c("G/S"))
      Allele96[,c("S")]<-(2*sum(CodTab[,c("cod96")]==c("S/S")))+sum(CodTab[,c("cod96")]==c("G/S"))
      Allele96[,c("Total")]<-sum(Allele96)
      aprop96<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("G","S")))
      aprop96[,c("G")]<-Allele96[,c("G")]/Allele96[,c("Total")]
      aprop96[,c("S")]<-Allele96[,c("S")]/Allele96[,c("Total")]
      Cod116<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("A/A", "A/G", "G/G", "Total")))
      Cod116[,c("A/A")]<-sum(CodTab[,c("cod116")]==c("A/A"))
      Cod116[,c("A/G")]<-sum(CodTab[,c("cod116")]==c("A/G"))
      Cod116[,c("G/G")]<-sum(CodTab[,c("cod116")]==c("G/G"))
      Cod116[,c("Total")]<-sum(Cod116)
      prop116<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("A/A", "A/G", "G/G")))
      prop116[,c("A/A")]<-Cod116[,c("A/A")]/Cod116[,c("Total")]
      prop116[,c("A/G")]<-Cod116[,c("A/G")]/Cod116[,c("Total")]
      prop116[,c("G/G")]<-Cod116[,c("G/G")]/Cod116[,c("Total")]
      Allele116<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("A", "G","Total")))
      Allele116[,c("A")]<-(2*sum(CodTab[,c("cod116")]==c("A/A")))+sum(CodTab[,c("cod116")]==c("G/G"))
      Allele116[,c("G")]<-(2*sum(CodTab[,c("cod116")]==c("G/G")))+sum(CodTab[,c("cod116")]==c("A/G"))
      Allele116[,c("Total")]<-sum(Allele116)
      aprop116<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("A","G")))
      aprop116[,c("A")]<-Allele116[,c("A")]/Allele116[,c("Total")]
      aprop116[,c("G")]<-Allele116[,c("G")]/Allele116[,c("Total")]
      GenProp<-matrix(0, nrow=3,ncol=3, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod95","cod96","cod116")))
      GenProp[,c("cod95")]<-t(prop95)
      GenProp[,c("cod96")]<-t(prop96)
      GenProp[,c("cod116")]<-t(prop116)
      AlleleProp<-matrix(0, nrow=2, ncol=3, dimnames=list(c("Major", "Minor"), c("cod95", "cod96", "cod116")))
      AlleleProp[,c("cod95")]<-t(aprop95)
      AlleleProp[,c("cod96")]<-t(aprop96)
      AlleleProp[,c("cod116")]<-t(aprop116)
      GenProp<-list(GenProp, AlleleProp, prop95, prop96, prop116, aprop95, aprop96, aprop116)
      names(GenProp)<-c("Genotype","Allele", "Geno95", "Geno96", "Geno116", "Allele95", "Allele96", "Allele116")
      return(GenProp)
    }
    if(x=="ODHE"){
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-c("cod225")
      Cod225<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("S/S", "S/F", "F/F", "Total")))
      Cod225[,c("S/S")]<-sum(CodTab[,c("cod225")]==c("S/S"))
      Cod225[,c("S/F")]<-sum(CodTab[,c("cod225")]==c("S/F"))
      Cod225[,c("F/F")]<-sum(CodTab[,c("cod225")]==c("F/F"))
      Cod225[,c("Total")]<-sum(Cod225)
      prop225<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("S/S", "S/F", "F/F")))
      prop225[,c("S/S")]<-Cod225[,c("S/S")]/Cod225[,c("Total")]
      prop225[,c("S/F")]<-Cod225[,c("S/F")]/Cod225[,c("Total")]
      prop225[,c("F/F")]<-Cod225[,c("F/F")]/Cod225[,c("Total")]
      Allele225<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("S", "F", "Total")))
      Allele225[,c("S")]<-(2*sum(CodTab[,c("cod225")]==c("S/S")))+sum(CodTab[,c("cod225")]==c("S/F"))
      Allele225[,c("F")]<-(2*sum(CodTab[,c("cod225")]==c("F/F")))+sum(CodTab[,c("cod225")]==c("F/F"))
      aprop225<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("S","F")))
      aprop225[,c("S")]<-Allele225[,c("S")]/Allele225[,c("Total")]
      aprop225[,c("F")]<-Allele225[,c("F")]/Allele225[,c("Total")]
      GenProp<-matrix(0, nrow=3,ncol=1, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod225")))
      GenProp[,c("cod225")]<-t(prop225)
      AlleleProp<-matrix(0, nrow=2, ncol=1, dimnames=list(c("Major", "Minor"), c("cod225")))
      AlleleProp[,c("cod225")]<-t(aprop225)
      GenProp<-list(GenProp, AlleleProp, prop225, aprop225)
      names(GenProp)<-c("Genotype","Allele", "Geno225", "Allele225")
      return(GenProp)
    }
    if(x=="other"){
      CodTab<-CodTab
      col<-colnames(CodTab)
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-col
      colnames(CodTab)<-paste("cod", colnames(CodTab), sep="")
      codname<-colnames(CodTab)
      CodTab[CodTab=="Error"]<-NA
      er<-as.matrix(CodTab[rowSums(is.na(CodTab))>0,])
      er<-rownames(er)
      CodTab<-CodTab[!rownames(CodTab)%in% er,]
      CodTab<-as.matrix(CodTab)
      labels<-list()
      for(i in 1:ncol(CodTab)){
        labels[[i]]<-unique(na.omit(CodTab[,i]))
        labels[[i]]<-strsplit(labels[[i]], split="/")
        labels[[i]]<-as.vector(unlist(labels[[i]]))
        labels[[i]]<-unique(na.omit(labels[[i]]))
        labels[[i]]<-sort(labels[[i]])
      }
      names(labels)<-codname
      n1<-vector()
      for(i in 1:length(labels)){
        n1[[i]]<-length(labels[[i]])
      }
      if(n1<2){
        labels1<-subset(labels, n1<2)
        for(i in 1:length(labels1)){
          labels1[[i]]<-c(paste(labels1[[i]][1], labels1[[i]][1], sep="/"), paste(labels1[[i]][1], c("-"), sep="/"), paste(c("-"),c("-"), sep="/"))
        }
        labels2<-labels1
      }
      if(n1==2){
        labels2<-subset(labels, n1==2)
        for(i in 1:length(labels2)){
          labels2[[i]]<-c(paste(labels2[[i]][1], labels2[[i]][1], sep="/"), paste(labels2[[i]][1], labels2[[i]][2], sep="/"), paste(labels2[[i]][2], labels2[[i]][2], sep="/"))
        }
        labels2<-labels2
      }
      if(n1==2&n1<2){
        labels1<-subset(labels, n1<2)
        for(i in 1:length(labels1)){
          labels1[[i]]<-c(paste(labels1[[i]][1], labels1[[i]][1], sep="/"), paste(labels1[[i]][1], c("-"), sep="/"), paste(c("-"),c("-"), sep="/"))
        }
        labels2<-subset(labels, n1==2)
        for(i in 1:length(labels2)){
          labels2[[i]]<-c(paste(labels2[[i]][1], labels2[[i]][1], sep="/"), paste(labels2[[i]][1], labels2[[i]][2], sep="/"), paste(labels2[[i]][2], labels2[[i]][2], sep="/"))
        }
        labels2<-append(labels1,labels2)
      }
      GenSum<-list()
      for(i in 1:length(labels2)){
        GenSum[[i]]<-matrix(0, nrow=1,ncol=4)
        rownames(GenSum[[i]])<-c("Count")
        colnames(GenSum[[i]])<-c(labels2[[i]][1],labels2[[i]][2],labels2[[i]][3], "Total")
      }
      names(GenSum)<-names(labels2)
      CodTab<-as.matrix(CodTab)
      for(i in 1:length(GenSum)){
        GenSum[[i]][1,1]<-sum(CodTab[,i]==labels2[[i]][1])
        GenSum[[i]][1,2]<-sum(CodTab[,i]==labels2[[i]][2])
        GenSum[[i]][1,3]<-sum(CodTab[,i]==labels2[[i]][3])
        GenSum[[i]][1,4]<-sum(GenSum[[i]][1,1:3])
      }
      GenProp<-matrix(0,3,length(GenSum), dimnames=list(c("Homozygote A", "Heterozygote", "Homozygote B"), names(GenSum)))
      for(i in 1:length(GenSum)){
        GenProp[1,i]<-(GenSum[[i]][1,1]/GenSum[[i]][1,4])
        GenProp[2,i]<-(GenSum[[i]][1,2]/GenSum[[i]][1,4])
        GenProp[3,i]<-(GenSum[[i]][1,3]/GenSum[[i]][1,4])
      }
      AlleleProp<-matrix(0, 2, ncol(GenProp), dimnames=list(c("Allele A", "Allele B"),colnames(GenProp)))
      for(i in 1:length(GenSum)){
        AlleleProp[1,i]<-((2*GenSum[[i]][1,1])+(GenSum[[i]][1,2]))/(2*GenSum[[i]][1,4])
        AlleleProp[2,i]<-((2*GenSum[[i]][1,3])+(GenSum[[i]][1,2]))/(2*GenSum[[i]][1,4])
      }
      labels3<-list()
      for(i in 1:length(labels2)){
        labels3[[i]]<-labels2[[i]][2] 
        labels3[[i]]<-strsplit(labels3[[i]], split="/")
      }
      GProp<-list()
      for(i in 1:length(GenSum)){
        GProp[[i]]<-t(GenProp[,i])
        colnames(GProp[[i]])<-unlist(labels2[[i]])
        rownames(GProp[[i]])<-c("Prop")
      }
      AProp<-list()
      for(i in 1:length(GenSum)){
        AProp[[i]]<-t(AlleleProp[,i])
        colnames(AProp[[i]])<-unlist(labels3[[i]])
        rownames(AProp[[i]])<-c("Prop")
      }
      GenoLab<-paste("Geno", col, sep="")
      AlleleLab<-paste("Allele", col, sep="")
      GenProp<-c(list(GenProp), list(AlleleProp), GProp, AProp)
      names(GenProp)<-c("Genotype", "Allele", GenoLab, AlleleLab)
      return(GenProp)
    }
  }
  if(Pop==T){
    if(x=="CECA"){
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-c("cod132")
      CodTab<-cbind(PopLabels, CodTab)
      PopLabels<-unique(CodTab[,1])
      PopList<-list()
      for(i in 1:length(unique(CodTab[,1]))){
        PopList[[i]]<-CodTab[CodTab[,1]==PopLabels[i],]
      }
      
      Cod132<-list()
      for(i in 1:length(PopList)){
        Cod132[[i]]<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("M/M", "M/L", "L/L", "Total")))
        Cod132[[i]][,c("M/M")]<-sum(PopList[[i]][,c("cod132")]==c("M/M"))
        Cod132[[i]][,c("M/L")]<-sum(PopList[[i]][,c("cod132")]==c("M/L"))
        Cod132[[i]][,c("L/L")]<-sum(PopList[[i]][,c("cod132")]==c("L/L"))
        Cod132[[i]][,c("Total")]<-sum(Cod132[[i]])
      }
      prop132<-list()
      for(i in 1:length(PopList)){
        prop132[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("M/M", "M/L", "L/L")))
        prop132[[i]][,c("M/M")]<-Cod132[[i]][,c("M/M")]/Cod132[[i]][,c("Total")]
        prop132[[i]][,c("M/L")]<-Cod132[[i]][,c("M/L")]/Cod132[[i]][,c("Total")]
        prop132[[i]][,c("L/L")]<-Cod132[[i]][,c("L/L")]/Cod132[[i]][,c("Total")]
      }
      Allele132<-list()
      for(i in 1:length(PopList)){
        Allele132[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("M", "L", "Total")))
        Allele132[[i]][,c("M")]<-(2*sum(PopList[[i]][,c("cod132")]==c("M/M")))+sum(PopList[[i]][,c("cod132")]==c("M/L"))
        Allele132[[i]][,c("L")]<-(2*sum(PopList[[i]][,c("cod132")]==c("L/L")))+sum(PopList[[i]][,c("cod132")]==c("M/L"))
        Allele132[[i]][,c("Total")]<-sum(Allele132[[i]])
      }
      aprop132<-list()
      for(i in 1:length(Allele132)){
        aprop132[[i]]<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("M","L")))
        aprop132[[i]][,c("M")]<-Allele132[[i]][,c("M")]/Allele132[[i]][,c("Total")]
        aprop132[[i]][,c("L")]<-Allele132[[i]][,c("L")]/Allele132[[i]][,c("Total")]
      }
      GenProp<-list()
      for(i in 1:length(prop132)){    
        GenProp[[i]]<-matrix(0, nrow=3,ncol=1, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod132")))
        GenProp[[i]][,c("cod132")]<-t(prop132[[i]])
      }
      AlleleProp<-list()
      for(i in 1:length(aprop132)){
        AlleleProp[[i]]<-matrix(0, nrow=2, ncol=1, dimnames=list(c("Major", "Minor"), c("cod132")))
        AlleleProp[[i]][,c("cod132")]<-t(aprop132[[i]])
      }
      GenAll<-list()
      for(i in 1:ncol(GenProp[[1]])){
        GenAll[[i]]<-matrix(0,length(GenProp),3, dimnames=list(c(PopLabels),c("Homozygote Major", "Heterozygote", "Homozygote Minor"))) 
      }
      for(i in 1:length(GenProp)){
        for(j in 1:ncol(GenProp[[1]])){
          GenAll[[j]][i,]<-GenProp[[i]][,j]
        }
      }
      AlleleAll<-list()
      for(i in 1:ncol(AlleleProp[[1]])){
        AlleleAll[[i]]<-matrix(0,length(AlleleProp),2, dimnames=list(c(PopLabels), c("Major","Minor"))) 
      }
      for(i in 1:length(AlleleProp)){
        for(j in 1:ncol(AlleleProp[[1]])){
          AlleleAll[[j]][i,]<-AlleleProp[[i]][,j]
        }
      }
      names1<-paste("Genotype",PopLabels, sep="_")
      names2<-paste("Allele", PopLabels, sep="_")
      names3<-paste("Geno132", PopLabels, sep="_")
      names4<-paste("Allele132", PopLabels, sep="_")
      names5<-paste("Genotype", "All", sep="_")
      names6<-paste("Allele", "All", sep="_")
      names<-c(names1, names2, names3, names4, names5, names6)
      GenProp<-c(GenProp, AlleleProp, prop132, aprop132, GenAll,AlleleAll)
      names(GenProp)<-names
      return(GenProp)
    }
    if(x=="ODVI"){
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-c("cod95", "cod96", "cod116")
      CodTab<-cbind(PopLabels, CodTab)
      PopLabels<-unique(CodTab[,1])
      PopList<-list()
      for(i in 1:length(unique(CodTab[,1]))){
        PopList[[i]]<-CodTab[CodTab[,1]==PopLabels[i],]
      }
      Cod95<-list()
      for(i in 1:length(PopList)){
        Cod95[[i]]<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("Q/Q", "Q/H", "H/H", "Total")))
        Cod95[[i]][,c("Q/Q")]<-sum(PopList[[i]][,c("cod95")]==c("Q/Q"))
        Cod95[[i]][,c("Q/H")]<-sum(PopList[[i]][,c("cod95")]==c("Q/H"))
        Cod95[[i]][,c("H/H")]<-sum(PopList[[i]][,c("cod95")]==c("H/H"))
        Cod95[[i]][,c("Total")]<-sum(Cod95[[i]])
      }
      prop95<-list()
      for(i in 1:length(Cod95)){
        prop95[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("Q/Q", "Q/H", "H/H")))
        prop95[[i]][,c("Q/Q")]<-Cod95[[i]][,c("Q/Q")]/Cod95[[i]][,c("Total")]
        prop95[[i]][,c("Q/H")]<-Cod95[[i]][,c("Q/H")]/Cod95[[i]][,c("Total")]
        prop95[[i]][,c("H/H")]<-Cod95[[i]][,c("H/H")]/Cod95[[i]][,c("Total")]
      }
      Allele95<-list()
      for(i in 1:length(PopList)){
        Allele95[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("Q", "H","Total")))
        Allele95[[i]][,c("Q")]<-(2*sum(PopList[[i]][,c("cod95")]==c("Q/Q")))+sum(PopList[[i]][,c("cod95")]==c("Q/H"))
        Allele95[[i]][,c("H")]<-(2*sum(PopList[[i]][,c("cod95")]==c("H/H")))+sum(PopList[[i]][,c("cod95")]==c("Q/H"))
        Allele95[[i]][,c("Total")]<-sum(Allele95[[i]])
      }
      aprop95<-list()
      for(i in 1:length(Allele95)){
        aprop95[[i]]<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("Q","H")))
        aprop95[[i]][,c("Q")]<-Allele95[[i]][,c("Q")]/Allele95[[i]][,c("Total")]
        aprop95[[i]][,c("H")]<-Allele95[[i]][,c("H")]/Allele95[[i]][,c("Total")]
      }
      Cod96<-list()
      for(i in 1:length(PopList)){
        Cod96[[i]]<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("G/G", "G/S", "S/S", "Total")))
        Cod96[[i]][,c("G/G")]<-sum(PopList[[i]][,c("cod96")]==c("G/G"))
        Cod96[[i]][,c("G/S")]<-sum(PopList[[i]][,c("cod96")]==c("G/S"))
        Cod96[[i]][,c("S/S")]<-sum(PopList[[i]][,c("cod96")]==c("S/S"))
        Cod96[[i]][,c("Total")]<-sum(Cod96[[i]])
      }
      prop96<-list()
      for(i in 1:length(Cod96)){
        prop96[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("G/G", "G/S", "S/S")))
        prop96[[i]][,c("G/G")]<-Cod96[[i]][,c("G/G")]/Cod96[[i]][,c("Total")]
        prop96[[i]][,c("G/S")]<-Cod96[[i]][,c("G/S")]/Cod96[[i]][,c("Total")]
        prop96[[i]][,c("S/S")]<-Cod96[[i]][,c("S/S")]/Cod96[[i]][,c("Total")]
      }
      Allele96<-list()
      for(i in 1:length(PopList)){
        Allele96[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("G", "S","Total")))
        Allele96[[i]][,c("G")]<-(2*sum(PopList[[i]][,c("cod96")]==c("G/G")))+sum(PopList[[i]][,c("cod96")]==c("G/S"))
        Allele96[[i]][,c("S")]<-(2*sum(PopList[[i]][,c("cod96")]==c("S/S")))+sum(PopList[[i]][,c("cod96")]==c("G/S"))
        Allele96[[i]][,c("Total")]<-sum(Allele96[[i]])
      }
      aprop96<-list()
      for(i in 1:length(Allele96)){
        aprop96[[i]]<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("G","S")))
        aprop96[[i]][,c("G")]<-Allele96[[i]][,c("G")]/Allele96[[i]][,c("Total")]
        aprop96[[i]][,c("S")]<-Allele96[[i]][,c("S")]/Allele96[[i]][,c("Total")]
      }
      Cod116<-list()
      for(i in 1:length(PopList)){
        Cod116[[i]]<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("A/A", "A/G", "G/G", "Total")))
        Cod116[[i]][,c("A/A")]<-sum(PopList[[i]][,c("cod116")]==c("A/A"))
        Cod116[[i]][,c("A/G")]<-sum(PopList[[i]][,c("cod116")]==c("A/G"))
        Cod116[[i]][,c("G/G")]<-sum(PopList[[i]][,c("cod116")]==c("G/G"))
        Cod116[[i]][,c("Total")]<-sum(Cod116[[i]])
      }
      prop116<-list()
      for(i in 1:length(Cod116)){
        prop116[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("A/A", "A/G", "G/G")))
        prop116[[i]][,c("A/A")]<-Cod116[[i]][,c("A/A")]/Cod116[[i]][,c("Total")]
        prop116[[i]][,c("A/G")]<-Cod116[[i]][,c("A/G")]/Cod116[[i]][,c("Total")]
        prop116[[i]][,c("G/G")]<-Cod116[[i]][,c("G/G")]/Cod116[[i]][,c("Total")]
      }
      Allele116<-list()
      for(i in 1:length(PopList)){
        Allele116[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("A", "G","Total")))
        Allele116[[i]][,c("A")]<-(2*sum(PopList[[i]][,c("cod116")]==c("A/A")))+sum(PopList[[i]][,c("cod116")]==c("G/G"))
        Allele116[[i]][,c("G")]<-(2*sum(PopList[[i]][,c("cod116")]==c("G/G")))+sum(PopList[[i]][,c("cod116")]==c("A/G"))
        Allele116[[i]][,c("Total")]<-sum(Allele116[[i]])
      }
      aprop116<-list()
      for(i in 1:length(Allele116)){
        aprop116[[i]]<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("A","G")))
        aprop116[[i]][,c("A")]<-Allele116[[i]][,c("A")]/Allele116[[i]][,c("Total")]
        aprop116[[i]][,c("G")]<-Allele116[[i]][,c("G")]/Allele116[[i]][,c("Total")]
      }
      GenProp<-list()
      for(i in 1:length(prop95)){    
        GenProp[[i]]<-matrix(0, nrow=3,ncol=3, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod95", "cod96", "cod116")))
        GenProp[[i]][,c("cod95")]<-t(prop95[[i]])
        GenProp[[i]][,c("cod96")]<-t(prop96[[i]])
        GenProp[[i]][,c("cod116")]<-t(prop116[[i]])
      }
      AlleleProp<-list()
      for(i in 1:length(aprop95)){
        AlleleProp[[i]]<-matrix(0, nrow=2, ncol=3, dimnames=list(c("Major", "Minor"), c("cod95", "cod96", "cod116")))
        AlleleProp[[i]][,c("cod95")]<-t(aprop95[[i]])
        AlleleProp[[i]][,c("cod96")]<-t(aprop96[[i]])
        AlleleProp[[i]][,c("cod116")]<-t(aprop116[[i]])
      }
      names1<-paste("Genotype",PopLabels, sep="_")
      names2<-paste("Allele", PopLabels, sep="_")
      names3<-paste("Geno95", PopLabels, sep="_")
      names4<-paste("Geno96", PopLabels, sep="_")
      names5<-paste("Geno116", PopLabels, sep="_")
      names6<-paste("Allele95", PopLabels, sep="_")
      names7<-paste("Allele96", PopLabels, sep="_")
      names8<-paste("Allele116", PopLabels, sep="_")   
      names<-c(names1, names2, names3, names4, names5, names6, names7, names8)
      GenProp<-c(GenProp, AlleleProp, prop95, prop96, prop116, aprop95, aprop96, aprop116)
      names(GenProp)<-names
      return(GenProp)
    }
    if(x=="ODHE"){
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-c("cod225")
      CodTab<-cbind(PopLabels, CodTab)
      PopLabels<-unique(CodTab[,1])
      PopList<-list()
      for(i in 1:length(unique(CodTab[,1]))){
        PopList[[i]]<-CodTab[CodTab[,1]==PopLabels[i],]
      }
      Cod225<-list()
      for(i in 1:length(PopList)){
        Cod225[[i]]<-matrix(0, nrow=1, ncol=4, dimnames=list(c("Count"), c("S/S", "S/F", "F/F", "Total")))
        Cod225[[i]][,c("S/S")]<-sum(PopList[[i]][,c("cod225")]==c("S/S"))
        Cod225[[i]][,c("S/F")]<-sum(PopList[[i]][,c("cod225")]==c("S/F"))
        Cod225[[i]][,c("F/F")]<-sum(PopList[[i]][,c("cod225")]==c("F/F"))
        Cod225[[i]][,c("Total")]<-sum(Cod225[[i]])
      }
      prop225<-list()
      for(i in 1:length(Cod225)){
        prop225[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Prop"), c("S/S", "S/F", "F/F")))
        prop225[[i]][,c("S/S")]<-Cod225[[i]][,c("S/S")]/Cod225[[i]][,c("Total")]
        prop225[[i]][,c("S/F")]<-Cod225[[i]][,c("S/F")]/Cod225[[i]][,c("Total")]
        prop225[[i]][,c("F/F")]<-Cod225[[i]][,c("F/F")]/Cod225[[i]][,c("Total")]
      }
      Allele225<-list()
      for(i in 1:length(PopList)){
        Allele225[[i]]<-matrix(0, nrow=1, ncol=3, dimnames=list(c("Count"), c("S", "F", "Total")))
        Allele225[[i]][,c("S")]<-(2*sum(PopList[[i]][,c("cod225")]==c("S/S")))+sum(PopList[[i]][,c("cod225")]==c("S/F"))
        Allele225[[i]][,c("F")]<-(2*sum(PopList[[i]][,c("cod225")]==c("F/F")))+sum(PopList[[i]][,c("cod225")]==c("F/F"))
      }
      aprop225<-list()
      for(i in 1:length(Allele225)){
        aprop225[[i]]<-matrix(0,nrow=1,ncol=2, dimnames=list(c("Prop"), c("S","F")))
        aprop225[[i]][,c("S")]<-Allele225[[i]][,c("S")]/Allele225[[i]][,c("Total")]
        aprop225[[i]][,c("F")]<-Allele225[[i]][,c("F")]/Allele225[[i]][,c("Total")]
      }
      GenProp<-list()
      for(i in 1:length(prop225)){    
        GenProp[[i]]<-matrix(0, nrow=3,ncol=1, dimnames=list(c("Homozygote Major", "Heterozygote", "Homozygote Minor"), c("cod225")))
        GenProp[[i]][,c("cod225")]<-t(prop225[[i]])
      }
      AlleleProp<-list()
      for(i in 1:length(aprop225)){
        AlleleProp[[i]]<-matrix(0, nrow=2, ncol=1, dimnames=list(c("Major", "Minor"), c("cod225")))
        AlleleProp[[i]][,c("cod225")]<-t(aprop225[[i]])
      }
      names1<-paste("Genotype",PopLabels, sep="_")
      names2<-paste("Allele", PopLabels, sep="_")
      names3<-paste("Geno225", PopLabels, sep="_")
      names4<-paste("Allele225", PopLabels, sep="_")    
      names<-c(names1, names2, names3, names4)
      GenProp<-c(GenProp, AlleleProp, prop225, aprop225)
      names(GenProp)<-names
      return(GenProp)
    }
    if(x=="other"){
      CodTab<-CodTab
      col<-colnames(CodTab)
      CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
      colnames(CodTab)<-col
      CodTab<-cbind(PopLabels, CodTab)
      CodTab[CodTab=="Error"]<-NA
      PopLabels<-unique(CodTab[,1])
      PopList<-list()
      for(i in 1:length(unique(CodTab[,1]))){
        PopList[[i]]<-CodTab[CodTab[,1]==PopLabels[i],,drop=F]
      }
      names(PopList)<-PopLabels
      for(i in 1:length(PopList)){
        PopList[[i]]<-as.matrix(PopList[[i]][,-c(1)])
      }
      for(i in 1:length(PopList)){
        colnames(PopList[[i]])<-col
      }
      CodTab1<-CodTab
      for(i in 1:length(PopList)){
        colnames(PopList[[i]])<-paste("cod", colnames(PopList[[i]]), sep="")
      }
      labels<-list()
      for(i in 1:length(PopList)){
        labels[[i]]<-list()
      }
      for(j in 1:length(labels)){
        for(i in 1:ncol(PopList[[i]])){
          labels[[j]][[i]]<-unique(na.omit(PopList[[j]][,i]))
          labels[[j]][[i]]<-strsplit(labels[[j]][[i]], split="/")
          labels[[j]][[i]]<-as.vector(unlist(labels[[j]][[i]]))
          labels[[j]][[i]]<-unique(na.omit(labels[[j]][[i]]))
          labels[[j]][[i]]<-sort(labels[[j]][[i]])
        }
      }
      for(i in 1:length(labels)){
        names(labels[[i]])<-colnames(PopList[[i]])
      }
      names(labels)<-PopLabels
      n1<-list()
      for(i in 1:length(labels)){
        n1[[i]]<-list()
      }
      for(i in 1:length(labels)){
        for(j in 1:length(labels[[i]])){
          n1[[i]][[j]]<-length(labels[[i]][[j]]) 
        }
      }
      labels1<-list()
      labels2<-list()
      for(i in 1:length(n1)){
        labels1[[i]]<-list()
        labels2[[i]]<-list()
      }
      for(i in 1:length(labels)){
        for(j in 1:length(labels[[i]])){
          labels1[[i]]<-subset(labels[[i]], n1[[i]]<2)
          labels2[[i]]<-subset(labels[[i]], n1[[i]]==2)
        }
      }
      names(labels1)<-names(labels2)<-PopLabels
      container1<-list()
      for(i in 1:length(labels1)){
        if(length(labels1[[i]])>0){
          container1[[i]]<-labels1[[i]]
        }
      }
      container2<-list()
      for(i in 1:length(labels2)){
        if(length(labels2[[i]])>0){
          container2[[i]]<-labels2[[i]]
        }
      }
      labels1<-container1
      labels2<-container2
      if(length(labels1)>0){
      for(i in 1:length(labels1)){
        labels1[[i]]<-c(paste(labels1[[i]][1], labels1[[i]][1], sep="/"), paste(labels1[[i]][1], c("-"), sep="/"), paste(c("-"),c("-"), sep="/"))
      }
      }
      if(length(labels1)<1){
        labels1<-labels1
      }
      if(length(labels2)>0){
      for(i in 1:length(labels2)){
        labels2[[i]]<-c(paste(labels2[[i]][[1]][1], labels2[[i]][[1]][1], sep="/"), paste(labels2[[i]][[1]][1], labels2[[i]][[1]][2], sep="/"), paste(labels2[[i]][[1]][2], labels2[[i]][[1]][2], sep="/"))                                                                  
      }
      }
      if(length(labels2)<1){
        labels2<-labels2
      }
      labels2<-append(labels1, labels2)
      names(labels2)<-PopLabels
      labels3<-list()
      for(i in 1:length(PopLabels)){
        labels3[[i]]<-subset(labels2, names(labels2)==c(PopLabels[[i]]))
      }
      names(labels3)<-PopLabels
      list<-list()
      for(i in 1:length(labels3[[i]])){
        list[[i]]<-matrix(0, nrow=1,ncol=4)
        rownames(list[[i]])<-c("Count")
        colnames(list[[i]])<-c(labels3[[i]][[i]][1], labels3[[i]][[i]][2], labels3[[i]][[i]][3], "Total")
      }
      names(list)<-names(labels[[1]])
      GenSum<-list()
      for(i in 1:length(labels3)){
        GenSum[[i]]<-list
      }
      names(GenSum)<-PopLabels
      
      for(i in 1:length(GenSum)){
        for(j in 1:length(GenSum[[i]])){
          GenSum[[i]][[j]][1,1]<-sum(PopList[[i]][,j]==labels3[[i]][[j]][1])
          GenSum[[i]][[j]][1,2]<-sum(PopList[[i]][,j]==labels3[[i]][[j]][2])
          GenSum[[i]][[j]][1,3]<-sum(PopList[[i]][,j]==labels3[[i]][[j]][3])
          GenSum[[i]][[j]][1,4]<-sum(GenSum[[i]][[j]][1,1:3])
        }
      }
      GenProp<-list()
      for(i in 1:length(labels3)){
        GenProp[[i]]<-matrix(0,3,length(GenSum[[i]]))
        colnames(GenProp[[i]])<-names(GenSum[[i]])
        rownames(GenProp[[i]])<-c("Homozygote A", "Heterozygote", "Homozygote B")
      }
      for(i in 1:length(GenProp)){
        for(j in 1:length(GenSum[[i]])){
          GenProp[[i]][1,j]<-(GenSum[[i]][[j]][1,1]/GenSum[[i]][[j]][1,4])
          GenProp[[i]][2,j]<-(GenSum[[i]][[j]][1,2]/GenSum[[i]][[j]][1,4])
          GenProp[[i]][3,j]<-(GenSum[[i]][[j]][1,3]/GenSum[[i]][[j]][1,4])
        }
      }
      AlleleProp<-list()
      for(i in 1:length(GenProp)){
        AlleleProp[[i]]<-matrix(0,2, ncol(GenProp[[i]]))
        colnames(AlleleProp[[i]])<-colnames(GenProp[[i]])
        rownames(AlleleProp[[i]])<-c("Allele A", "Allele B")
      }
      for(i in 1:length(AlleleProp)){
        for(j in 1:ncol(GenProp[[i]])){
          AlleleProp[[i]][1,j]<-((2*GenSum[[i]][[j]][1,1])+(GenSum[[i]][[j]][1,2]))/(2*GenSum[[i]][[j]][1,4])
          AlleleProp[[i]][2,j]<-((2*GenSum[[i]][[j]][1,3])+(GenSum[[i]][[j]][1,2]))/(2*GenSum[[i]][[j]][1,4])
        }
      }
      GProp<-list()
      for(i in 1:length(GenSum)){
        GProp[[i]]<-list()
      }
      for(i in 1:length(GenSum)){
        for(j in 1:length(GenSum[[1]])){
          GProp[[i]][[j]]<-t(GenProp[[i]][,j])
          colnames(GProp[[i]][[j]])<-labels3[[i]][[j]]
          rownames(GProp[[i]][[j]])<-c("Prop")
        }
      }
      labels4<-list()
      for(i in 1:length(labels3)){
        labels4[[i]]<-list()
      }
      for(i in 1:length(labels3)){
        for(j in 1:length(labels3[[1]])){
          labels4[[i]][[j]]<-labels3[[i]][[j]][[2]]
          labels4[[i]][[j]]<-strsplit(labels4[[i]][[j]], "/")
        }
      }
      names(labels4)<-names(labels3)
      
      AProp<-list()
      for(i in 1:length(GenSum)){
        AProp[[i]]<-list()
      }
      for(i in 1:length(GenSum)){
        for(j in 1:length(GenSum[[1]])){
          AProp[[i]][[j]]<-t(AlleleProp[[i]][,j])
          colnames(AProp[[i]][[j]])<-unlist(labels4[[i]][[j]])
          rownames(AProp[[i]][[j]])<-c("Prop")
        }
      }
      Genotype_All<-matrix(0,length(PopLabels),nrow(GenProp[[1]]))
      for(i in 1:length(GenProp)){
        Genotype_All[i,]<-t(GenProp[[i]])
      }
      
      Allele_All<-matrix(0, length(PopLabels), nrow(AlleleProp[[1]]))
      for(i in 1:length(AlleleProp)){
        Allele_All[i,]<-t(AlleleProp[[i]])
      }
      rownames(Allele_All)<-rownames(Genotype_All)<-PopLabels
      codlab<-vector()
      colnew<-vector()
      for(i in 1:length(col)*3){
        colnew<-c(rep(col[1],3))
      }
      zygote<-c("Homozygote A", "Heterozygote", "Homozygote B")
      for(i in 1:length(colnew)){
        codlab[i]<-paste(colnew[i],zygote[i], sep="_")
      }
      all<-c("Allele A", "Allele B")
      colnew2<-vector()
      for(i in 1:length(col)*2){
        colnew2<-c(rep(col[1],2))
      }
      codlab2<-vector()
      for(i in 1:length(colnew2)){
        codlab2[i]<-paste(colnew2[i], all[i], sep="_")
      }
      colnames(Genotype_All)<-codlab
      colnames(Allele_All)<-codlab2
      
      list1<-c(GenProp,AlleleProp, unlist(GProp, recursive=F), unlist(AProp, recursive=F), list(Genotype_All),list(Allele_All))
      GenoLab<-paste("Genotype", PopLabels, sep="_")
      AlleleLab<-paste("Allele", PopLabels,sep="_")

      GLab<-paste("Geno",col, sep="")
      GLab1<-paste(GLab, PopLabels,sep="_")
      ALab<-paste("Allele",col, sep="")
      ALab1<-paste(ALab, PopLabels, sep="_")
      names<-c(GenoLab, AlleleLab, GLab1, ALab1, "Genotype_All","Allele_All")
      names(list1)<-names
      return(list1)
    }
  }
}


#################################

###Write Individual Genotypes###
write.Tabs<-function(x, NucTab, CodTab, filename){
  if(x=="CECA"){
    NucTab<-as.matrix(NucTab[!rownames(NucTab)%in% c("Ref"),])
    colnames(NucTab)<-c("nuc394")
    CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
    colnames(CodTab)<-c("cod132")
    Ind<-list(NucTab, CodTab)
    names(Ind)<-c("SNPs", "AminoAcids")
    write.xlsx(Ind, file=filename, rowNames=T)
  }
  if(x=="ODVI"){
    NucTab<-as.matrix(NucTab[!rownames(NucTab) %in% c("Ref"),])
    colnames(NucTab)<-c("nuc285","nuc286","nuc347")
    CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
    colnames(CodTab)<-c("cod95","cod96","cod116")
    Ind<-list(NucTab, CodTab)
    names(Ind)<-c("SNPs", "AminoAcids")
    write.xlsx(Ind, file=filename, rowNames=T)
  }
  if(x=="other"){
    nuclab<-colnames(NucTab)
    NucTab<-as.matrix(NucTab[!rownames(NucTab)%in% c("Ref"),])
    colnames(NucTab)<-nuclab
    codlab<-colnames(CodTab)
    CodTab<-as.matrix(CodTab[!rownames(CodTab)%in% c("Ref"),])
    colnames(CodTab)<-codlab
    Ind<-list(NucTab, CodTab)
    names(Ind)<-c("SNPs","AminoAcids")
    write.xlsx(Ind, file=filename, rowNames=T, overwrite=T)
  }
}
#################################

###Write Prop Table###
write.PropTab<-function(PropTab, filename){
  write.xlsx(PropTab, file=filename, rowNames=T)
}
##################

####SpatialInterp####
IDW<-function(method,Sus,PropTab,Centroid,X,Y,power,shapefile){
  if(method=="Genotype"){
    if(Sus=="More"){
      PropTab<-PropTab
      More<-PropTab[,1]
      Centroid<-Centroid
      More<-cbind(More, Centroid)
      Count<-as.matrix(table(PopLabels))
      colnames(Count)<-c("N")
      shape<-readOGR(".", layer=shapefile)
      extent<-extent(shape)
      ext<-as.matrix(extent)
      loc<-as.data.frame(More[,2:3])
      colnames(loc)<-c("x","y")
      coordinates(loc)<-~x+y
      grd<-expand.grid(x=seq(from=ext[1,1], to=ext[1,2],by=X), y=seq(from=ext[2,1], to=ext[2,2],by=Y))
      coordinates(grd)<-~x+y
      gridded(grd)<-TRUE
      moreidw<-idw(More[,1]~1, locations=loc, newdata=grd, idp=power)
      idwraster<-raster(moreidw)
      cr<-crop(idwraster, extent(shape), snap="out")
      fr<-rasterize(shape, cr)
      lr<-raster::mask(x=cr, mask=fr)
      return(lr)
    }
    if(Sus=="Less"){
      PropTab<-PropTab
      Less<-PropTab[,2]
      Centroid<-Centroid
      Less<-cbind(Less, Centroid)
      Count<-as.matrix(table(PopLabels))
      colnames(Count)<-c("N")
      shape<-readOGR(".", layer=shapefile)
      extent<-extent(shape)
      ext<-as.matrix(extent)
      loc<-as.data.frame(Less[,2:3])
      colnames(loc)<-c("x","y")
      coordinates(loc)<-~x+y
      grd<-expand.grid(x=seq(from=ext[1,1], to=ext[1,2],by=X), y=seq(from=ext[2,1], to=ext[2,2],by=Y))
      coordinates(grd)<-~x+y
      gridded(grd)<-TRUE
      lessidw<-idw(Less[,1]~1, locations=loc, newdata=grd, idp=power)
      lessraster<-raster(lessidw)
      lesscr<-crop(lessraster, extent, snap="out")
      lessfr<-rasterize(shape, lesscr)
      lesslr<-raster::mask(x=lesscr, mask=lessfr)
      return(lesslr)
    }
  }
}
#######################################

####KDE####
KDE<-function(method,Sus,locus,genotypes,shapefile,Codons,XYdat){
  if(method=="Genotype"){
    if(Sus=="More"){
      poly<-readOGR(".", layer=shapefile)
      Codons<-Codons
      inddata<-as.matrix(Codons[!rownames(Codons)%in% c("Ref"),])
      colnames(inddata)<-colnames(Codons)
      genotypes<-genotypes
      More<-as.matrix(inddata[inddata[,locus]==genotypes[1],])
      More<-as.data.frame(More)
      XYdat<-XYdat
      XYdat<-XYdat[rownames(XYdat)%in%rownames(More),]
      XYdat<-as.data.frame(XYdat)
      colnames(XYdat)<-c("x","y")
      bbox2<-extent(poly)
      bbox2<-as.matrix(bbox2, nrow=2, ncol=2)
      bbox2<-t(bbox2)
      bbox<-matrix(0,5,2)
      bbox[1,]<-c(bbox2[1,1], bbox2[1,2])
      bbox[2,]<-c(bbox2[1,1], bbox2[2,2])
      bbox[3,]<-c(bbox2[2,1], bbox2[2,2])
      bbox[4,]<-c(bbox2[2,1], bbox2[1,2])
      bbox[5,]<-c(bbox2[1,1], bbox2[1,2])
      mse.MM<-mse2d(as.points(XYpts),poly=bbox,nsmse=100,100000)
      bw<-mse.MM$h[which.min(mse.MM$mse)]
      plot(mse.MM$h[5:100], mse.MM$mse[5:100], type="l")
      kde<-kernel2d(as.points(XYdat),poly=bbox,h0=bw, nx=1000, ny=1000)
      kderaster<-raster(kde)
      cr<-crop(kderaster, extent(poly), snap="out")
      fr<-rasterize(poly, cr)
      lr<-raster::mask(x=cr, mask=fr)
      return(lr)
    }
    if(Sus=="Less"){
      poly<-readOGR(".", layer=shapefile)
      Codons<-Codons
      inddata<-as.matrix(Codons[!rownames(Codons)%in% c("Ref"),]) 
      colnames(inddata)<-colnames(Codons)
      genotypes<-genotypes
      Less<-as.matrix(inddata[inddata[,locus]==genotypes[2],])
      Less2<-as.matrix(inddata[inddata[,locus]==genotypes[3],])
      Less<-rbind(Less, Less2)
      Less<-as.data.frame(Less)
      XYdat<-XYdat
      XYdat<-XYdat[rownames(XYdat)%in%rownames(Less),]
      XYdat<-as.data.frame(XYdat)
      colnames(XYdat)<-c("x","y")
      bbox2<-extent(poly)
      bbox2<-as.matrix(bbox2, nrow=2, ncol=2)
      bbox2<-t(bbox2)
      bbox<-matrix(0,5,2)
      bbox[1,]<-c(bbox2[1,1], bbox2[1,2])
      bbox[2,]<-c(bbox2[1,1], bbox2[2,2])
      bbox[3,]<-c(bbox2[2,1], bbox2[2,2])
      bbox[4,]<-c(bbox2[2,1], bbox2[1,2])
      bbox[5,]<-c(bbox2[1,1], bbox2[1,2])
      mse.MM<-mse2d(as.points(XYdat),poly=bbox,nsmse=100,100000)
      bw<-mse.MM$h[which.min(mse.MM$mse)]
      plot(mse.MM$h[5:100], mse.MM$mse[5:100], type="l")
      kde<-kernel2d(as.points(XYdat),poly=bbox,h0=bw, nx=1000, ny=1000)
      kderaster<-raster(kde)
      cr<-crop(kderaster, extent(poly), snap="out")
      fr<-rasterize(poly, cr)
      lr<-raster::mask(x=cr, mask=fr)
      return(lr)
    }
    }
}
########################################

####floating pie graph####
Pie<-function(shapefile, PropTab, XYdat, radius){
  poly<-readOGR(".", layer=shapefile)
  plot(poly)
  for(i in 1:nrow(PropTab)){
    floating.pie(XYdat[i,"x"], XYdat[i,"y"],
                 c(PropTab[i,1], PropTab[i,2])+PropTab[i,3], radius=sum(PropTab[i,])*radius, col=c("red", "blue"))
  }
  legend(x="topright", legend=c("More Susceptible", "Less Susceptible"), pch=c(19,19),col=c("red","blue"))
}
#########################################