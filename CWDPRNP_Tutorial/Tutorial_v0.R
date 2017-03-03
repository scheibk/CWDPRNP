#############################Load Package################################
###Description:
##Loads package and associated dependencies.  

###Arguments:
##There are no arguments to call

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial")
source("CWDPRNP_ver0.R")
CheckCRANPackages(packages=c("seqinr", "ape","openxlsx","sp","splancs","rgdal","maptools","gstat", "raster","plotrix"))
CheckBioconductorPackages(packages=c("DECIPHER", "sangerseqR","Biostrings","msa"))

#Note, this package also requires Rtools if running on Windows (https://cran.r-project.org/bin/windows/Rtools/)

#Note:  Make sure to set set all working directories to your specified file extensions.  The extensions provided here are for example.

#Note:  Folders have results from tutorial already available.  R can overwrite these files as long as these are not open in another program.  Rerunning this tutorial should not produce an error related to file access.

#####################################################

###################Module 1: Alignment and Variation#####################

####Call Reference Sequence####
###Description:  
##Calls reference sequence from GenBank.  Default commands utilize 771 bp fragments that span from start codon to stop codon.  User can specify accession number of other GenBank files if desired.
###

###Arguments:
##x:  A character string specifying species of interest or accession number (see details below) 
###

###Details:
##Cervus canadensis (x==c("CECA"))
##Odocoileus virginianus (x==c("ODVI"))
##Odocoileus hemionus (x==c("ODHE"))
##Other accession number (x==c("Accession Number")) -> ex. (x==c("AY639093")) for Rangifer tarandus
###

###Values:
##Returns an R object of class DNAString. 
###

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/ab1_Data")
ref<-RefSeq(x=c("CECA"))

##########################################

####Make Base Calls####
###Description:  
##Uploads all .ab1 files from working directory and performs base calls.  Function is meant to work with two files from each sample, one forward read and one reverse read.  It is suggested to label each file in the following format:
##Forward Files:  "SampleID_F.ab1"
##Reverse Files:  "SampleID_R.ab1"
###

###Arguments:
##flab:  a character string identifying the common label of all forward reads
##rlab:  a character string identifying the common label of all reverse reads
##label:  a character string identifying the point at which sample ID can be subset from rest of filename ("*text.ab1")
###

###Values:
##Returns a list of length 2.  The first element is a list containing all forward reads as sangerseq objects.  The second element is a list containing all reverse reads as sangerseq objects.  Note that you must assign each list to its own object in order to proceed.
###

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/ab1_Data")
myfiles<-CallBases(flab=c("*F_223F.ab1"), rlab=c("*R_224R.ab1"), label=c("\\R_224R.*"))
forward<-myfiles[[1]]
reverse<-myfiles[[2]]

###########################################

####Write Chromatograms####

###Description:  
##Function will produce a .pdf file for each element of the called list.
###

###Arguments:
##files:  A list of sangerseq objects (see: "Make Base Calls")
##ref:  A DNAString object containing the reference sequence (see:  "Call Reference Sequence")
###

###Values:
##Produces a .pdf file in the specified working directory for each element the list specified by argument 'files'.  The name of each file will be the name of each individual sangerseq object in the specified list. 
###

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Chromatogram")
Chrom(forward, reverse)

############################################

####Perform Sequence Alignment####

###Description:
##This function aligns and clips all forward and reverse reads to the extent of the reference sequence.  
###

###Arguments:
##ref:  A DNAString object containing the reference sequence (see:  "Call Reference Sequence")
##forward:  A list of forward sequences of class sangerseq
##reverse:  A list of reverse sequences of class sangerseq
##method:  Indicates which alignment algorithm to call ("DECIPHER","ClustalW", "ClustalOmega","Muscle")
###

###Details:
##Alignments are performed for primary (allele 1) and secondary (allele 2) calls for forward and reverse sequences.  Reverse sequences are transformed into their forward complement sequence prior to alignment.
###

###Values:
##Produces as list of length 4.  Each element contains a matrix of aligned sequences.  List[[1]]=forward primary, list[[2]]=forward secondary, list[[3]]=reverse primary, list[[4]]=reverse secondary.
###

###Example:
align<-Align(ref=ref, forward=forward, reverse=reverse, method="ClustalW")

#Note:  We have noticed misalignments or errors with several methods, particularly with large dataset.  We have found ClustalW to work in the vast majority of cases; however, this method is the slowest.  If you experience any issues (e.g. misalignments) with other methods, we would suggest going back to ClustalW.  It is up to the user to determine the tradeoff between speed and accuracy.
#############################################

####Identify Polymorphic Sites####
###Description:
##This function is used to identify potentially polymorphic sites.  This is useful when a researcher is interested in identifying polymorphic sites outside of SNPs previously associated with disease susceptibility
###

###Arguments:
##poly:  A logical argument (True, False) that identifies whether a row should be added indicating level of polymorphism
##align:  A list object containing aligned sequences
##ref:  A DNAString object containing reference sequence
###

###Values:
##Produces a list of length 2.  Each element contains a matrix of aligned sequences.  List[[1]]=forward primary and secondary calls, list[[2]]=reverse primary and secondary calls.  If poly=T, then the last row will indicate the level of polymorphism (1,...,n; depending on the number of potential alleles).
###

###Example:
polysites<-Poly(poly=T, align=align, ref=ref)

#######################################################

####Extract Polymorphic Sites####
###Description:
##This function extracts polymorphic sites from the results of the full alignment.
###

###Arguments:
##align:  A list object produced from the poly function
###

###Details:
##This function produces three tables:  (1) Extracted sites with forward and reverse calls combined, (2) Extracted sites for forward calls only, and (3)  Extracted sites for reverse calls only.  The table of interest, in most cases, is list[[1]].  This table is used to guide the researcher in choosing sites (column names) that are polymorphic.  Note that sequencing errors, labeled using some non-nucleic acid indicator will flag a site as polymorphic.  It is up to the researcher to determine which sites should be used for subsequent functions.  The row indicating the level of polymorphism in the Poly() function is also a helpful guide.  It is highly unlikely anything with a level of polymorphism >2 is a true SNP;however, genotyping errors can cause polymorphic sites to flag with a level of polymorphism >3.  

###Example
SNPs<-ExtractSNPs(polysites)

#################################################


#############Module 2:  Summarize Prion Gene Variation###################
#####################Option 1:  Default SNPs#############################
###The code below provides extracts genotypes and genotype frequencies for disease-associated loci###

####Nucleotide Table####
###Description:  
##Function to extract genotypes for SNPs associated with reduced CWD susceptibility.
###

###Arguments:
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##align:  A list object containing aligned sequences
###

###Details:
##Subsets known polymorphisms that have been shown to influence disease susceptibility (see Robinson et al. 2012; Haley and Hoover 2015).  
###

###Values:
##Returns a matrix of individual genotypes at each disease-associated locus.
###

###Example:
Nucleotides<-NucTab(x=c("CECA"), align=align)

###########################################

####Codon Table####
###Description:
##Function to extract amino acid polymorphisms associated with reduced CWD susceptibility.
###

###Arguments:
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##NucTab:  A matrix of SNP genotypes for each sampled individual (see: "Nucleotide Table")
###

###Values:
##Returns a matrix of individual amino acid types at each disease-associated locus.
###

###Example:
Codons<-CodTab(x=c("CECA"), Nucleotides)

####Write Individual Results####
###Description:
##Function to write an .xlsx document containing a table of SNPs and amino acid polymorphisms, referenced by individual.
###

###Arguments:
##x:  A character string specifying species of interest (see:  "Call Reference Sequence")
##NucTab:  A matrix of SNP genotypes for each sampled individual (see: "Nucleotide Table")
##CodTab:  A matrix of amino acid types (see:  "Codon Table")
##filename:  A character string containing the desired excel file name
###

###Values:
##Produces a .xlsx file with two worksheets.  Worksheet 1 (SNPs) contains the NucTab matrix.  Worksheet 2 (AminoAcids) contains the CodTab matrix.
###

###Example
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Results")
write.Tabs(x="CECA", Nucleotides, Codons, filename="elkdefaulttest.xlsx")

###############################################

####Genotype/Allele Frequency Tables####
###Description:  
##A function to calculate the genotype and allele frequencies for each disease-associated locus.
###

###Arguments:  
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##CodTab:  A matrix of amino acid types for each sampled individual (see:  "Codon Table")
###

###Values:
##Produces a list object containing several matrices (varies by species).  This list contains a matrix summarizing genotype frequencies across all loci, a matrix summarizing genotype frequencies for each locus, a matrix summarizing allele frequencies across all loci, a matrix summarizing allele frequencies for each locus.  
###

###Example:
##Without Population Labels
Prop<-PropTab(Pop=F, x=c("CECA"),Codons)

##With Population Labels
#Load Population labels and XY coordinates#
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Spatial_Data")
Popxy<-read.csv("LocalityInfo.csv")
Popxy<-Popxy[order((Popxy[,1])),]
labels<-Popxy[,1]
Popxy<-Popxy[,-1]
PopLabels<-as.character(Popxy[,1])

Prop<-PropTab(Pop=T, x=c("CECA"), Codons, PopLabels=PopLabels)

###############################################

####Write Proportion Tables to Excel####
###Description:
##Writes an excel file (.xlsx) containing allele and genotype frequencies
###

###Arguments:
##PropTab:  A list item containing allele and genotype frequency tables
##filename:  A character string containing the desired excel file name
###

###Values:
##Produces a .xlsx file with a worksheet for every element of PropTab
###

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Results")
write.PropTab(Prop, filename="ElkPropTab.xlsx")
###############################################

####################Option 2:  User-defined SNPs#########################
####Nucleotide Table####
###Description:  
##Function to extract genotypes for SNPs associated with reduced CWD susceptibility.
###

###Arguments:
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##align:  A list object containing aligned sequences
##nuc:  A vector of SNPs
###

###Details:
##Subsets alignment items by user-defined SNPs.  
###

###Values:
##Returns a matrix of individual genotypes.
###

###Example:
Nucleotides<-NucTab(x=c("other"), align=align, nuc=c(63))

###########################################

####Codon Table####
###Description:
##Function to extract amino acid polymorphisms for user-defined SNPs.
###

###Arguments:
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##Tab:  A matrix of SNP genotypes for each sampled individual (see: "Nucleotide Table")
##align:  Item containing alignment matrices ("See Align function")
###

###Values:
##Returns a matrix of individual amino acid types for user-defined loci.
###

###Example:
Codons<-CodTab(x=c("other"), Nucleotides, align=align)

####Write Individual Results####
###Description:
##Function to write an .xlsx document containing a table of SNPs and amino acid polymorphisms, referenced by individual.
###

###Arguments:
##x:  A character string specifying species of interest (see:  "Call Reference Sequence")
##NucTab:  A matrix of SNP genotypes for each sampled individual (see: "Nucleotide Table")
##CodTab:  A matrix of amino acid types (see:  "Codon Table")
##filename:  A character string containing the desired excel file name
###

###Values:
##Produces a .xlsx file with two worksheets.  Worksheet 1 (SNPs) contains the NucTab matrix.  Worksheet 2 (AminoAcids) contains the CodTab matrix.
###

###Example
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Results")
write.Tabs(x="other", Nucleotides, Codons, filename="elkusertest.xlsx")

###############################################

####Genotype/Allele Frequency Tables####
###Description:  
##A function to calculate the genotype and allele frequencies for each disease-associated locus.

###IMPORTANT NOTE:  This function has a known bug that only allows a user to pass one user-defined SNP at a time.  For example, if x="other" you can only run SNP 63 or SNP 394, not both at the same time.  Make sure to truncate NucTab and CodTab functions to only pull one SNP until problem is resolved.  
###
###

###Arguments:  
##x:  A character string specifying species of interest (see: "Call Reference Sequence")
##CodTab:  A matrix of amino acid types for each sampled individual (see:  "Codon Table")
###

###Values:
##Produces a list object containing several matrices (varies by species).  This list contains a matrix summarizing genotype frequencies across all loci, a matrix summarizing genotype frequencies for each locus, a matrix summarizing allele frequencies across all loci, a matrix summarizing allele frequencies for each locus.  
###

###Example:
##Without Population Labels
Prop<-PropTab(Pop=F, x=c("other"),Codons)

##With Population Labels
#Load Population labels and XY coordinates#
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Spatial_Data")
Popxy<-read.csv("LocalityInfo.csv")
Popxy<-Popxy[order((Popxy[,1])),]
labels<-Popxy[,1]
Popxy<-Popxy[,-1]
PopLabels<-as.character(Popxy[,1])

Prop<-PropTab(Pop=T, x=c("other"), Codons, PopLabels=PopLabels)

###############################################

####Write Proportion Tables to Excel####
###Description:
##Writes an excel file (.xlsx) containing allele and genotype frequencies
###

###Arguments:
##PropTab:  A list item containing allele and genotype frequency tables
##filename:  A character string containing the desired excel file name
###

###Values:
##Produces a .xlsx file with a worksheet for every element of PropTab
###

###Example:
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Results")
write.PropTab(Prop, filename="ElkUserTab.xlsx")
###############################################

####################Module 3:  Spatial Analyses##########################
####Note:  All following examples refer to analysis on default SNPs (Option 1).  It is possible to pass these functions over user-defined SNPs as well.
####

####IMPORTANT NOTE:  Matrix with geographic coordinates MUST be in the same spatial projection as the defined shapefile BEFORE running the following functions.
####

####Rerun default examples with PropTab(Pop=T,...).  This will reproduce example figures from publication.  Note that the kernel density map may be slightly different, as jittering was used to produce random noise in individual locations:

Nucleotides<-NucTab(x=c("CECA"), align=align)
Codons<-CodTab(x=c("CECA"), Nucleotides)

#Population Labels#
setwd("~/Side Projects/CWDPRNP/CWDPRNP_Tutorial/Spatial_Data")
Popxy<-read.csv("LocalityInfo.csv")
Popxy<-Popxy[order((Popxy[,1])),]
labels<-Popxy[,1]
Popxy<-Popxy[,-1]
PopLabels<-as.character(Popxy[,1])

Prop<-PropTab(Pop=T, x=c("CECA"), Codons, PopLabels=PopLabels)
####

####Individual Analysis - Kernel Density Estimation####
###Description:
##Produces a raster item containing a kernel density map for the defined level of susceptibility
###

###Arguments:
##Method:  Method must be set to =c("Genotype").  Future updates will provide additional functionality.
##Sus:  A character string indicating whether to map more susceptible genotypes (Sus="More") or less susceptible genotypes (Sus="Less").
##locus:  A character string indicating codon (e.g. "cod132")
##genotypes:  Possible genotypes for locus (e.g. c("M/M","M/L","L/L"))
##shapefile:  A character string indicating the name of a GIS shapefile in the working directory.  This should represent the sampling area.
##Codons:  A matrix of amino acid types for each sampled individual (see:  "Codon Table").
##XYdat:  A matrix of individual locations.  Must match order of matrix indicated in Codons argument.  
###

###Values:
##Produces a raster item containing a kernel density map for the defined level of susceptibility.  
###

###Example:
#Jitter XY coordinates for KDE example#
Xj<-jitter(Popxy[,2],factor=6)
Yj<-jitter(Popxy[,3],factor=6)
XYj<-cbind(Xj, Yj)
Popxy<-cbind(as.numeric(XYj[,1]),as.numeric(XYj[,2]))
colnames(Popxy)<-c("x","y")
rownames(Popxy)<-labels

poly<-readOGR(".", layer="ElkSelection2")
plot(poly)
points(Popxy)
kdetest<-KDE(method="Genotype",Sus="Less",locus="cod132",genotypes=c("M/M","M/L","L/L"), shapefile="ElkSelection2",Codons=Codons, XYdat=Popxy[,1:2])
plot(poly, axes=T)
plot(kdetest, add=T, alpha=0.75)

##########################################

####Population Analysis - Spatial Interpolation####
###Description:
##Produces a raster item containing interpolated population-level proportions of indicated genotype.
###

###Arguments:
##Method:  Method must be set to =c("Genotype").  Future updates will provide additional functionality.
##Sus:  A character string indicating whether to map more susceptible genotypes (Sus="More") or less susceptible genotypes (Sus="Less").
##PropTab:  A table containing the population-level genotype frequencies.  This corresponds to the table named "Genotype_All" in Prop item  
##Centroid:  A table containing geographic coordinates for population centroids
##X:  Numerical value determining cell-size for spatial grid on x-axis.
##Y:  Numerical value determining cell-size for spatial grid on y-axis.
##power:  Numerical value determining shape of inverse distance weighting function for spatial interpolation.
##shapefile:  A character string indicating the name of a GIS shapefile in the working directory.  This should represent the sampling area.
##
###

###Values:
##Produces a raster item containing interpolated population-level proportions of indicated genotype.  Note a trade-off between computing power and grid size parameters.  Larger values are best but take a considerable amount of computing time.  The power parameter is decided by the user.  It is probably best to try multiple values.  
###

###Example:
GenProp<-Prop[["Genotype_All"]]

#Format Centroid file#
Popxy<-read.csv("LocalityInfo.csv")
Popxy<-Popxy[order((Popxy[,1])),]
labels<-Popxy[,1]
Popxy<-Popxy[,-1]
PopLabels<-as.character(Popxy[,1])
Popxy<-unique(Popxy)
Cent<-Popxy
Cent<-cbind(as.numeric(Cent[,2]), as.numeric(Cent[,3]))
colnames(Cent)<-c("x","y")

idwraster<-IDW(method="Genotype", Sus="Less", PropTab=GenProp, Centroid=Cent,X=100, Y=100, power=3, shapefile="ElkSelection2")

#Map Results#
poly<-readOGR(".", layer="ElkSelection2")
plot(poly, axes=T)
plot(idwraster, add=T, alpha=0.75)

#################################

####Population Analysis - Population Frequencies####
###Description:
##Produces a map of population-specific genotype frequencies.
###

###Arguments:
##shapefile:  A character string indicating the name of a GIS shapefile in the working directory.  This should represent the sampling area.
##PropTab:  A table containing the population-level genotype frequencies.  This corresponds to the table named "Genotype_All" in Prop item  
##XYdat:  A table containing geographic coordinates for population centroids
##radius:  A scaling-constant to decide pie chart size (scaled by population size).
##
###

###Values:
##Produces a map of population-specific genotype frequencies.  Pie charts are mapped to centroid.  The radius argument is a scaling constant.  Try different values and find one that best fits study area.  
###

###Example:
Pie(shapefile="ElkSelection2",PropTab=GenProp,XYdat=Cent, radius=10000)

#######################################