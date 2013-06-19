library(GEOquery)
source('heatmap.plus.R') 
library(heatmap.plus)
source(file="Func_List.r")
library(matlab)
library(stringr)
library(plyr)
library("hgu133plus2.db") # Affymetrix Human Genome U133 Plus 2.0 Array annotation data
library("hwgcod.db") # GE/Amersham Codelink Human Whole Genome Bioarray 
library("hgug4112a.db") # 
options(stringsAsFactors=FALSE)
##TODO - get more datasets

## ovarian cancer tissue from ascites-cytology-positive or -negative patients
## log2 GC-RMA signal
#gse39204 <- getGEO("GSE39204") 
#gse39204.d <- exprs(gse39204$GSE39204_series_matrix.txt.gz)
#gse39204.s <- phenoData(gse39204$GSE39204_series_matrix.txt.gz)@data
#
## 38 ovarian cancer cell lines (13 OCCC cell lines and 25 non-OCCC cell lines
## log 2 values of RMA signal intensity
#gse29175 <- getGEO("GSE29175") 
#gse29175.d <- exprs(gse29175$GSE29175_series_matrix.txt.gz)
#gse29175.s <- phenoData(gse29175$GSE29175_series_matrix.txt.gz)@data
#
## 10 clear cell ovarian cancer specimens and 10 normal ovarian surface epithelium
## log2 RMA signal 
#gse29450 <- getGEO("GSE29450") 
#gse29450.d <- exprs(gse29450$GSE29450_series_matrix.txt.gz)
#gse29450.s <- phenoData(gse29450$GSE29450_series_matrix.txt.gz)@data
#
## 58 stage III-IV ovarian cancer patients treated with Carboplatin and Taxol agents
## log2 signal processed by RMA
#gse30161 <- getGEO("GSE30161") 
#gse30161.d <- exprs(gse30161$GSE30161_series_matrix.txt.gz)
#gse30161.s <- phenoData(gse30161$GSE30161_series_matrix.txt.gz)@data
#
###tothill data (Endometriosis and Ovarian cancer - no normals.)
#tothill <- getGEO("GSE9891")
#tothill.d1 <- exprs(tothill$`GSE9891_series_matrix-1.txt.gz`)
#tothill.d2 <- exprs(tothill$`GSE9891_series_matrix-2.txt.gz`)
#all(dimnames(tothill.d1)[[1]]==dimnames(tothill.d2)[[1]])
##combined 
#tothill.d <- cbind(tothill.d1, tothill.d2)
## uses phenoData method on one particular object of the tothill list and 
##access the data slot(a.k.a. attribute)
#tothill.s1 <- phenoData(tothill$`GSE9891_series_matrix-1.txt.gz`)@data
#tothill.s2 <- phenoData(tothill$`GSE9891_series_matrix-2.txt.gz`)@data
#all(names(tothill.s1)==names(tothill.s2))
#tothill.s <- rbind(tothill.s1, tothill.s2) ##This is metadata!
#
###GSE6008 - 103 samples
#gse6008 <- getGEO("GSE6008")
#gse6008.d <- (exprs(gse6008$GSE6008_series_matrix.txt.gz))
#gse6008.s <- phenoData(gse6008$GSE6008_series_matrix.txt.gz)@data
#all(dimnames(gse6008.s)[[1]]==names(gse6008.d))
#
###GSE7305 - 20 samples 10 endo/ovary and 10 endo normaiil 
## platform Affymetrix U133 plus 2.0 array
#gse7305 <- getGEO("GSE7305")
#gse7305.d <- (exprs(gse7305$GSE7305_series_matrix.txt.gz))
#gse7305.s <- phenoData(gse7305$GSE7305_series_matrix.txt.gz)@data
#all(dimnames(gse7305.s)[[1]]==names(gse7305.d))
#
###GSE7307 - list of a variety of tumor/normal samples (need to select only relevant cases
#sample.info.gse7307 <- read.delim(file="/data/htorres/kate/GSE7307_GEO_Sample_Info.csv")
#gse7307.disease <- (subset(sample.info.gse7307, Disease.type!="Normal"))
#gse7307.endo.normal <- (subset(sample.info.gse7307, Tissue.Cell.Line..C.=="endometrium"))
#gse7307.s <- rbind(gse7307.disease, gse7307.endo.normal)
#
#
##GSE37837 # NOT NORMALIZED apparently
#gse37837 <- getGEO('GSE37837')
#gse37837.d <- as.data.frame(exprs(gse37837[[1]]))
#gse37837.s <- phenoData(gse37837$GSE37837_series_matrix.txt.gz)@data
#
#
##GSE5108 
## already normalized (Why kate says otherwise on spreadsheet?)
## 22 samples - 11 Eutopic Endometrium vs 11 Ectopic Endometriosis
## I will consider all strings containing:
## - 'eutopic' as normal controls.
## - 'endometriosis' as cancer.  
## Thus: 11 normals, 11 cancer.
## We may stratify different cancer grades latter.
#gse5108 = getGEO('GSE5108')
#gse5108.d = as.data.frame(exprs(gse5108[[1]]))
#gse5108.s = phenoData(gse5108$GSE5108_series_matrix.txt.gz)@data 
# all(dimnames(gse5108.s)[[1]] == names(gse5108.d))
#
#dataframe <- NULL
#for(i in 1:length(gse7307.s$Sample)){
#	if(is.null(dataframe)){
#		test <- getGEO(as.character(gse7307.s[i,1]))
#		dataframe <- Table(test)
#		dataframe[,2] <- as.numeric(dataframe[,2])
#		dimnames(dataframe)[[2]] <- c("ID_REF",as.character(gse7307.s[i,1]))
#		#print("Entrou IF") #Debug
#	}
#	else{
#		test <- getGEO(as.character(gse7307.s[i,1]))
#		test.d <- Table(test)
#		test.d[,2] <- as.numeric(test.d[,2])
#		dimnames(test.d)[[2]] <- c("ID_REF",as.character(gse7307.s[i,1]))
#		dataframe = merge(dataframe, test.d, by.x="ID_REF", by.y="ID_REF")
#		#print("ELSE") #Debug
#	}
#	#Progresso
#	print(i)
#	#print(files[i])
#}
#gse7307.d <- dataframe
#dimnames(gse7307.d)[[1]] <- gse7307.d[,1]
#gse7307.d <- gse7307.d[,-1]
#all(names(gse7307.d)==gse7307.s$Sample)
#
###GSE6364 - list of a variety of tumor/normal samples (need to select only relevant cases
#gse6364 <- getGEO("GSE6364")
#gse6364.d <- (exprs(gse6364$GSE6364_series_matrix.txt.gz))
#gse6364.s <- phenoData(gse6364$GSE6364_series_matrix.txt.gz)@data
#
###GSE23339 18 endmetrioma(ovarian EMS) vs 20 eutopic endometrium
###platform: illumina human6-v2.0 (GPL6102) - seems normalized.
#gse23339 <- getGEO("GSE23339")
#gse23339.d <- (exprs(gse23339$GSE23339_series_matrix.txt.gz))
#gse23339.s <- phenoData(gse23339$GSE23339_series_matrix.txt.gz)@data
#
#save(
#    gse29175.d,
#    gse29175.s,
#    gse39204.s,
#    gse39204.d,
#    gse29450.d,
#    gse29450.s,
#    gse30161.d, 
#    gse30161.s,
#    gse6364.d, 
#    gse6364.s, 
#    gse7307.d, 
#    gse7307.s, 
#    gse7305.d, 
#    gse7305.s, 
#    tothill.d, 
#    tothill.s, 
#    gse6008.d, 
#    gse6008.s,
#    gse5108.s,
#    gse5108.d, 
#    gse23339.s,
#    gse23339.d,
#    gse37837.s,
#    gse37837.d,
#    file="GEOfiles_associated_withENDO.downloaded.April2013.rda"
#    )
####start from here#####
setwd("/data/htorres/kate")
load(file="GEOfiles_associated_withENDO.downloaded.April2013.rda")
#
### New Src signature data frame
Src.signature = read.table('src_signature.txt', sep='\t') #from Bild et al 2006 suplemental table1
colnames(Src.signature) = c('ProbeID', 'GeneSymbol', 'Description', 'LocusLink', 'FoldChange')
### adjust this criteria
Src.signature$logFC <- log2(Src.signature$FoldChange)
Src.signature$direction <- NA
#Src.signature$direction[Src.signature$FoldChange<1] <- c("DN.BILD") # just to denote this info is not from our analysis (BILD) 
#Src.signature$direction[Src.signature$FoldChange>1] <- c("UP.BILD")
Src.signature$direction[Src.signature$logFC<0] <- c("DN.BILD") # just to denote this info is not from our analysis (BILD) 
Src.signature$direction[Src.signature$logFC>0] <- c("UP.BILD")

# src Gene signature file 
# obtained by conveting affy probe names of each platform with DAVID (http://david.abcc.ncifcrf.gov/)
#david.src = read.table('David.Src.txt', header=T) #maybe I should use the corresponding R/BioC Annotation package for consistency
#david.src$To=toupper(as.character(david.src$To))
#k <- c('PTHLH', 'ZNF12', 'SRC', 'PMS2L3', 'IRS1')
#colnames(david.src) <- c('affyProbe', 'GeneSymbol')
# getting illumina IDs:
# davidIllumina <- read.table('david_affy_illumina.txt', header=T, sep='\t') #submitted david.src$From to DAVID conversion tool and got Illumina IDs 
#indices <- match(david.src$affyProbe, Src.signature$ProbeID)

#update Src.signature with updated gene symbols from the bioconductor bimap
Src.signature$GeneSymbol <- unlist(
    mget(Src.signature$ProbeID, hgu133plus2SYMBOL, ifnotfound=NA), 
    use.names=F)
rename <- with(Src.signature, which(is.na(GeneSymbol)))
Src.signature$GeneSymbol[rename] <- Src.signature$ProbeID[rename]
Src.signature.out <- Src.signature[,c(1:2,7)]
colnames(Src.signature.out)[3] <- 'Bild.direction'
src = Src.signature #some of houtan's code reference this variable name
# used to subset Codelink probeIDs:
tmp <- unlist(
    mget(unique(Src.signature$GeneSymbol), hwgcodALIAS2PROBE, ifnotfound=NA),
    use.names=F)
codelinkSrc <- unique(tmp[!is.na(tmp)]) #Src.signature GE/Amersham probe IDs 
rm(tmp)
#used to subset agilent 4112 probe IDs
ag44112aSrc <- unlist(
    mget(unique(Src.signature$GeneSymbol), hgug4112aALIAS2PROBE, ifnotfound=NA),
    use.names=F)
##get probe specific to SRC signature
gse6008.src.probe <- intersect(src[,1], dimnames(gse6008.d)[[1]])
tothill.src.probe <- intersect(src[,1], dimnames(tothill.d)[[1]])
gse7305.src.probe <- intersect(src[,1], dimnames(gse7305.d)[[1]])
gse7307.src.probe <- intersect(src[,1], dimnames(gse7307.d)[[1]])
gse6364.src.probe <- intersect(src[,1], dimnames(gse6364.d)[[1]])
gse39204.src.probe <- intersect(src[,1], dimnames(gse39204.d)[[1]])
gse29175.src.probe <- intersect(src[,1], dimnames(gse29175.d)[[1]])
gse29450.src.probe <- intersect(src[,1], dimnames(gse29450.d)[[1]])
gse30161.src.probe <- intersect(src[,1], dimnames(gse30161.d)[[1]])

##The probenames should be prefixed by 'GE' but for some reason its not. WARNING! Verify why.
rownames(gse5108.d) <- paste("GE",  rownames(gse5108.d), sep='') 
gse5108.src.probe <- intersect(codelinkSrc, dimnames(gse5108.d)[[1]])
gse37837.src.probe <- intersect(ag44112aSrc, rownames(gse37837.d))

##subset data to only SRC specific and create dataframe. 
tothill.src <- tothill.d[tothill.src.probe,]; tothill.src <- as.data.frame(tothill.src)
gse6008.src <- gse6008.d[gse6008.src.probe,]; gse6008.src <- as.data.frame(gse6008.src)
gse6364.src <- gse6364.d[gse6364.src.probe,]; gse6364.src <- as.data.frame(gse6364.src)
gse7307.src <- gse7307.d[gse7307.src.probe,]; gse7307.src <- as.data.frame(gse7307.src)
gse7305.src <- gse7305.d[gse7305.src.probe,]; gse7305.src <- as.data.frame(gse7305.src)
gse39204.src <- gse39204.d[gse39204.src.probe,]; gse39204.src <- as.data.frame(gse39204.src)
gse29175.src <- gse29175.d[gse29175.src.probe,]; gse29175.src <- as.data.frame(gse29175.src)
gse29450.src <- gse29450.d[gse29450.src.probe,]; gse29450.src <- as.data.frame(gse29450.src)
gse30161.src <- gse30161.d[gse30161.src.probe,]; gse30161.src <- as.data.frame(gse30161.src)
## TODO insert rows of NAs in place of the missing Src signature genes (non-GPL570 platforms)
gse5108.src <- as.data.frame(gse5108.d[gse5108.src.probe,])
gse37837.src <- as.data.frame(gse37837.d[gse37837.src.probe,])

sort.data.frame <- function(x, key, ...) {
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}

###cleanup sample manifest #metadata) # WARNING: I Don't understand this section
#further code expects every manifest dataframe (df.s) to contain a characteristics_ch1a column
tothill.s$characteristics_ch1a <- tothill.s$characteristics_ch1.2
gse6008.nn <- strsplit(as.character(gse6008.s$characteristics_ch1), ": ")
gse6008.nnn <- as.data.frame(as.character(unlist(lapply((gse6008.nn), "[", 2) )))
gse6008.s$characteristics_ch1a <- gse6008.nnn[,1]
### endometrium is normal & endometrium/ovary is diseased
gse7307.s$characteristics_ch1a <- gse7307.s$Tissue.Cell.Line..C. 
gse7305.nn <- strsplit(as.character(gse7305.s$characteristics_ch1), " ")
gse7305.nnn <- as.data.frame(as.character(unlist(lapply((gse7305.nn), "[", 1) )))
gse7305.s$characteristics_ch1a <- gse7305.nnn[,1]
gse5108.s$characteristics_ch1a <- gse5108.s$characteristics_ch1
gse6364.s$characteristics_ch1a <- gse6364.s$characteristics_ch1
gse37837.s$characteristics_ch1a <- gse37837.s$characteristics_ch1.2
gse39204.s <- within(gse39204.s, characteristics_ch1a <- characteristics_ch1.3) #Does the same thing.
gse29175.s$characteristics_ch1a <- gse29175.s$characteristics_ch1.4
gse29450.s$characteristics_ch1a <- gse29450.s$characteristics_ch1.1
gse30161.s$characteristics_ch1a <- gse30161.s$characteristics_ch1.4

#Sort the metadataframes:
tothill.s.sort <- sort.data.frame(tothill.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse6008.s.sort <- sort.data.frame(gse6008.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse6364.s.sort <- sort.data.frame(gse6364.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse7307.s.sort <- sort.data.frame(gse7307.s[,c("Sample","characteristics_ch1a")], 
    key = "characteristics_ch1a")
gse7305.s.sort <- sort.data.frame(gse7305.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse5108.s.sort <- sort.data.frame(gse5108.s[,c("geo_accession","characteristics_ch1")], 
	key = "characteristics_ch1")
gse37837.s.sort <- sort.data.frame(gse37837.s[,c("geo_accession", "characteristics_ch1a")],
    key= "characteristics_ch1a")
gse39204.s.sort <- sort.data.frame(gse39204.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse29175.s.sort <- sort.data.frame(gse29175.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse29450.s.sort <- sort.data.frame(gse29450.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")
gse30161.s.sort <- sort.data.frame(gse30161.s[,c("geo_accession","characteristics_ch1a")], 
	key = "characteristics_ch1a")

## ColSideColors - red is cancer, grey is normal, #CB4B16' #not endometrioid cancer
color.map.gse29175.type <- function(mol.biol) {
    if (str_detect(mol.biol, 'non-OCCC')) '#CB4B16' #different ovarian cancer
    else if (mol.biol == 'histology: OCCC') 'red' #OCCC is considered endometriosis-derived cancer
    else '#FFFFFF' #white
}
color.map.gse39204.type <- function(mol.biol) {
    if ( str_detect(mol.biol, 'histology: clear') | 
         str_detect(mol.biol, 'endometrioid') ) 'red'
    else if ( str_detect(mol.biol, 'serous') | 
              str_detect(mol.biol, 'mucinous') | 
              str_detect(mol.biol, 'carcinosarcoma') )'#CB4B16' #different kind of cancer 
    else '#FFFFFF' #white
    }
color.map.gse39204.type2 <- function(mol.biol) {
    if ( str_detect(mol.biol, 'histology: clear') ) '#FF3300' #lightred 
    else if ( str_detect(mol.biol, 'endometrioid') ) '#FFCC00' #yellowish
    else if ( str_detect(mol.biol, 'serous')  )'#0000FF' 
    else if ( str_detect(mol.biol, 'mucinous') )'#6600FF'
    else if ( str_detect(mol.biol, 'carcinosarcoma') )'#CC00FF' #different kind of cancer 
    else '#FFFFFF' #white
    }
color.map.gse29450.type <- function(mol.biol) {
    if ( str_detect(mol.biol, 'clear cell ovarian cancer') ) 'red'
    else if ( str_detect(mol.biol, 'Normal') ) "#9C897E" #grey
    else '#FFFFFF' #white
    }
color.map.gse30161.type <- function(mol.biol) {
    if ( str_detect(mol.biol, 'Endometrioid') |
         str_detect(mol.biol, 'Clear') ) 'red'
    else if ( str_detect(mol.biol, 'Mucinous') |
              str_detect(mol.biol, 'Serous') ) '#CB4B16' #different kind of cancer
    else '#FFFFFF' # white
    }
color.map.gse37837.type <- function(mol.biol) {
    if ( str_detect(mol.biol, 'ectopic') ) "red"
    else if ( str_detect(mol.biol, 'eutopic') ) "#9C897E" #grey
    else "#FFFFFF" #white
    }
color.map.tothill.type.houtan <- function(mol.biol) {
	if (mol.biol=="Fallopian") "green" #green
	else if (mol.biol=="Ovary") "#red" #red
	else if (mol.biol=="Peritoneum") "black"
	else "#FFFFFF" #white
    } 
color.map.tothill.type <- function(mol.biol) {
    if (str_detect(mol.biol, 'Endo')) 'red'
    else if ( str_detect(mol.biol, 'Ser') ) '#CB4B16' #orangy
    else if ( str_detect(mol.biol, 'Adeno') ) '#B58900' #Yellowy
    else "#FFFFFF"
    }
color.map.gse6008.type <- function(mol.biol) {
	if (mol.biol=="Clear_Cell") "green" #green
	else if (mol.biol=="Endometrioid") "red" #red
	else if (mol.biol=="Mucinous") "black"
	else if (mol.biol=="Serous") "yellow"
	else "#FFFFFF"
    } 
color.map.gse7307.type <- function(mol.biol) {
	if (mol.biol=="myometrium") "green" #green
	else if (mol.biol=="endometrium/ovary") "red" 
	else if (mol.biol=="endometrium") "#9C897E" #grey
	else "#FFFFFF"
    }       
color.map.gse7305.type <- function(mol.biol) {
	if (mol.biol=="Endometrium/Ovary-Disease") "red" 
	else if (mol.biol=="Endometrium-Normal") "#9C897E" #grey
	else "#FFFFFF"
    } 
color.map.gse6364.type <- function(mol.biol) {
	if (mol.biol=="Early Secretory Phase Endometriosis") "red"
	else if (mol.biol=="Proliferative Phase Endometriosis") "red" 
	else if (mol.biol=="Mid Secretory Phase Endometriosis") "red" 
	else if (mol.biol=="Mid Secretory Phase Normal") "#9C897E" 
	else if (mol.biol=="Proliferative Phase Normal") "#9C897E" 
	else if (mol.biol=="Early Secretory Phase Normal") "#9C897E" 
	else "#FFFFFF"
    } 
color.map.gse5108.type <- function(mol.biol) {
    if (str_detect(mol.biol, 'endometriosis')) "red"
    else if (str_detect(mol.biol, 'eutopic endometrium')) "#9C897E" #grey
    else '#FFFFFF'
    }
#rowsidecolors - red is upregulated, green is downregulated
bildColors.colorize <- function(abbrv){ 
    if ( is.na(abbrv) ) '#FFFFFF' # missing
	else if (abbrv == 'UP.BILD') 'red'
	else if (abbrv == 'DN.BILD') 'green'
	else '#FFFFFF' #insignificant cases 
    }    
klColors.colorize = function(abbrv){
	if (abbrv == 'UP.KL') {return('red')}
	else if (abbrv == 'DN.KL') {return('green')}
	else return('#FFFFFF') #insignificant cases
    }
build.side.map = function(df){
	### takes dataframe, merges with the Src signature and creates RowSideColors for our plots
	## df is the full dataframe, with FC and P values
    ## ONLY WORKS FOR affy human genome U133 v2.0
	df$ProbeID = dimnames(df)[[1]]
    Src.signature$ProbeID <- as.character(Src.signature$ProbeID)
    # The following line was messing up the row (probes) ordering
    src.sig <- join(df, Src.signature, by=c("ProbeID"))  
	src.sig$Ldirection = NA # our expression directionality vector
	src.sig$Ldirection[src.sig$FC<0]  <- c("DN.KL") 
	src.sig$Ldirection[src.sig$FC>0] <- c("UP.KL")
	src.sig$Lsignificance[src.sig$p.value<=0.05] <- c("Significant")
	#I want the white color to mark the statistically insignificant cases, so:
	src.sig$Ldirection[which(is.na(src.sig$Lsignificance))] <- 'insignificant'
	bildColors = unlist(lapply(src.sig$direction, bildColors.colorize))
	klColors = unlist(lapply(src.sig$Ldirection, klColors.colorize))
	colors = cbind(bildColors,klColors)
	colnames(colors) = c('Bild et al.', 'Lawrenson et al.')
	return(colors)
}

build.side.map.plus <- function(df, probemap=hgu133plus2GENENAME) {
	### takes dataframe, merges with the david.src signature and creates RowSideColors
	## df: the Bild SRC subsetted GSE dataframe, with samples, FC, mean.T, mean.N and p.value
    ## probe2geneSymbol.db is a biocondutor "bimap" annotation object 
    df$probeID <- rownames(df)
    df$GeneSymbol <- unlist( mget(df$probeID, probemap, ifnotfound=NA), use.names=F )
    #having the geneSymbol, get the gene expression directionality in the Src.signature
    indices <- match(df$GeneSymbol, Src.signature$GeneSymbol)
    df$direction <- Src.signature$direction[indices]
    df$direction[is.na(df$direction)] <- 'unknown' 
    df$Ldirection <- NA
    df$Ldirection[df$FC > 0] <- c('DN.KL')
	df$Ldirection[df$FC>0] <- c("UP.KL")
	df$Lsignificance[df$p.value<=0.05] <- c("Significant")
	#I want the white color to mark the statistically insignificant cases, so:
	df$Ldirection[which(is.na(df$Lsignificance))] <- 'insignificant'
    df$Ldirection[is.na(df$Ldirection)] <- 'unknown' #DAVID makes ZN12 disappear!
	bildColors = unlist(lapply(df$direction, bildColors.colorize))
	klColors = unlist(lapply(df$Ldirection, klColors.colorize))
	colors = cbind(bildColors,klColors)
	colnames(colors) = c('Bild et al.', 'Lawrenson et al.')
	return(colors) 
}

build.top.map = function(df, df.s, df.colormap.function, df.colormap.function2=NA,
        labels=c('group1', 'group2') ){
	#builds the heatmap top rows for helping group visualization
	## df is the full dataframe, with FC and P values
	## df.colormap.function maps metadata with colors
	cols = as.matrix(colnames(df[1:(ncol(df)-4)]))
	if ('Sample' %in% colnames(df.s) ) {
	    df.src.type = apply(cols,
		    1,
		    func.list$vlookup,
		    df.s[,c("Sample","characteristics_ch1a")],
		    "characteristics_ch1a")
	} else { 
        df.src.type = apply(cols,
	    1,
	    func.list$vlookup,df.s[,c("geo_accession","characteristics_ch1a")],
            "characteristics_ch1a")
	}
	df.src.type = as.character(df.src.type)
	df.src.type[is.na(df.src.type)] = c('0')
	cc.df.src.type = unlist(lapply(df.src.type, df.colormap.function))
	ColSideColors = matrix(as.character(c(cc.df.src.type)), 
			nrow = length(cc.df.src.type),
			ncol = 1)
    if (is.na(df.colormap.function2) ) {
	    ColSideColors = cbind(ColSideColors, ColSideColors) #2columns make it thicker 
    } else {
        df.src.type2 = as.character(df.src.type)
        df.src.type2[is.na(df.src.type2)] = c('0')
        cc.df.src.type2 = unlist(lapply(df.src.type2, df.colormap.function2))
        ColSideColors2 = matrix(as.character(c(cc.df.src.type2)), 
                nrow = length(cc.df.src.type2),
                ncol = 1)
	    ColSideColors = cbind(ColSideColors2, ColSideColors) #2columns make it thicker 
        colnames(ColSideColors) <- labels
    }
	return(ColSideColors)
}

removeColumns <- function(df) {
    #removes spurious columns that get in the way of plotting just the samples
	## df is the full dataframe, with FC and P values
    filterout <- c('FC', 'p.value', 'mean.T', 'mean.N')
    # R syntax can become ugly. This is how I acchieve removal by colname
    filterout <-
    -1 * unlist(lapply(filterout, function(x){which(colnames(df) == x)} ) )
    df <- df[,filterout]
    return(df)
}

srcGeneLookup <- function(probes) {
    return(
        as.character(unlist(
            lapply(probes, func.list$vlookup, Src.signature, 'GeneSymbol')
            )
        )
    )
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make volcano plot plus ##TODO: GENESYMBOL MAPPING NOT WORKING
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
volcano.plus<-
function(x,fold.change.col,pv.col, title.plot, cut.line, fold.cut1, 
    fold.cut2, pv.adj.col,ncolors=1, text=NA, angle=-45, probemap=NA) {
    #text <-- a logical vector used to subset x    
    #cut.line <-- p.value cutoff line (not log10 transformed)
    op <- par(no.readonly = TRUE)
    p<-x[,pv.col]
    if (all(is.na(p))==FALSE) {
        M<-x[,fold.change.col]
        upr = p<=cut.line & M >=fold.cut1
        dwr = p<=cut.line & M<=fold.cut2
        q <- x[,pv.adj.col]
        q.f <- x[(upr|dwr),pv.adj.col]
        minq.f <- 0
        maxq.f <- 0
        if (length(q.f) > 0) {
            minq.f <- min(q.f)
            maxq.f <- max(q.f)
            if (is.na(minq.f)) {
                print(q.f)
            }
        }
        library('marray')
        Gcol <- maPalette(low = "#C8FFC8", high = "#006400", k = ncolors)
        Rcol <- maPalette(low = "#FFC8C8", high = "#640000", k = ncolors)
        xlim.r<-(range(x[,fold.change.col]))
        if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
            par(fig=c(0, 0.8, 0, 1), mar=c(4, 4, 4, 1))
        }
        if(abs(xlim.r[1])>abs(xlim.r[2])){
            plot(
                    x[,fold.change.col], #x-axis
                    -1*log10(x[,pv.col]), #y-axis
                    xlim=c(xlim.r[1], abs(xlim.r[1])), #x-axis limits
                    main=title.plot,
                    xlab="Gene Expression\nlog2(fold change)",
                    ylab="-1 * log10 of the Significance",
                    cex=.5,
                    pch=20
            )
            if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
                points(
                    M[upr], -log10(p[upr]), 
                    col=Rcol[(as.integer(((q[upr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],
                    cex=1,
                    pch=20)
                points(M[dwr], -log10(p[dwr]),
                    col=Gcol[(as.integer(((q[dwr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],
                    cex=1,
                    pch=20)
            }   
            else {
                points(M[upr], -log10(p[upr]), col=Rcol[ncolors],cex=1,pch=20)
                points(M[dwr], -log10(p[dwr]), col=Gcol[ncolors],cex=1,pch=20)         
            }

            #points(M[upr], -log10(p[upr]), col="red")
            #points(M[dwr], -log10(p[dwr]), col="green")

            abline(h= -log10(cut.line), lty=3, lwd=1)
            abline(v= fold.cut1, lty=3, lwd=1, col="black")
            abline(v= fold.cut2, lty=3, lwd=1, col="black")
        }else{
            plot(
                    x[,fold.change.col], #x-axis
                    -1*log10(x[,pv.col]), #y-axis
                    xlim=c(-1*xlim.r[2], xlim.r[2]), #x-axis limits
                    main=title.plot,
                    xlab="Gene Expression\nlog2(fold change)",
                    ylab="-1 * log10 of the Significance",
                    cex=.5,
                    pch=20
            )

            if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
                points(M[upr],
                    -log10(p[upr]),
                    col=Rcol[(as.integer(((q[upr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],
                    cex=1,
                    pch=20)
                points(M[dwr], 
                    -log10(p[dwr]), 
                    col=Gcol[(as.integer(((q[dwr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],
                    cex=1,
                    pch=20)
            }
            else {
                points(M[upr], -log10(p[upr]), col=Rcol[ncolors],cex=1,pch=20)
                points(M[dwr], -log10(p[dwr]), col=Gcol[ncolors],cex=1,pch=20)         
            }
            abline(h= -log10(cut.line), lty=3, lwd=1)
            abline(v= fold.cut1, lty=3, lwd=1, col="black")
            abline(v= fold.cut2, lty=3, lwd=1, col="black")
        }
        if( !all(is.na(text)) ){
            probes <- rownames(x)[text]
            #print( c('probes', probes) )
            if (is.na(probemap) | missing(probemap) ) { 
                genelabels <- srcGeneLookup(probes)
            } else {
                genelabels <- unlist(
                    mget(probes, probemap, ifnotfound=NA), 
                    use.names=F)
            }
            #print( c('probemap', probemap) )
            #print( c('genelabels', genelabels) )
            text(M[text],
                -log10(p[text]), 
                labels=genelabels,
                cex=0.5,
                srt=angle,
                offset=0.25,
                pos=4)        
        }
        if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
            ColorLevels <- seq( minq.f, maxq.f, ((maxq.f-minq.f)/(ncolors-1)))
            #print(paste(minq.f,maxq.f))
            #print(ncolors)
            #print(ColorLevels)
            #Stay on same page and set up region and coordinates for legend.
            par(new=TRUE)
            par(fig=c(0.8, 0.9, 0.2, 0.8), mar=c(4, 0, 4, 3))
            #maColorBar(ColorLevels, col=Gcol, horizontal=FALSE, k=0, cex.axis=.8)
            image(1, seq(1,ncolors,1),
                matrix(data=seq(1,ncolors,1), ncol=ncolors, nrow=1),
                col=Gcol,
                xlab="",ylab="",
                axes=FALSE
            )
            axis(4, at = seq(0.5, (ncolors-0.5), 1), 
                    labels = rep('',ncolors),cex.axis=.5,mgp=c(0, .2, 0)
            )
            par(new=TRUE)
            par(fig=c(0.9, 1, 0.2, 0.8), mar=c(4, 0, 4, 3))
            #maColorBar(ColorLevels, col=Rcol, horizontal=FALSE, k=ncolors, cex.axis=.8)
            image(1, seq(1,ncolors,1),
                    matrix(data=seq(1,ncolors,1), ncol=ncolors, nrow=1),
                    col=Rcol,
                    xlab="",ylab="",
                    axes=FALSE
            )
            axis(4, at = seq(0.5, (ncolors-0.5), 1), 
                    labels = signif(ColorLevels,2),cex.axis=.5,mgp=c(0, .2, 0),
            )
        }
    }   
    par(op)
}
plotVolcano <- function(df, 
    filename='volcano.svg', 
    title='Volcano Plot', 
    probemap=NA) {
    #custom volcano plot customized for kate's analysis
    svg(filename=filename)
#    png(filename=filename)
    volcano.plus(df, 
        fold.change.col = "FC", 
        pv.col = "p.value", 
        title.plot= title, 
        cut.line = 0.05, 
        fold.cut1 = 0,
        fold.cut2 = 0,
        pv.adj.col = df$p.value,
        ncolors = 5,
        #text = (df$FC < -0.5 | df$FC > 0.5) & -log10(df$p.value) > -log10(0.05),
        text = -log10(df$p.value) > -log10(0.05),
        angle = 45,
        probemap = probemap
        ) 
    dev.off()
} 
plotHeatmap <-
function(df, df.s, colormap, 
    colormap2=NA, 
    filename='Heatmap.svg', 
    title='Heatmap', 
    probemap=NA, 
    labels=c('group1', 'group2')) {
        svg(filename=filename)
        #colnames I always want removed:
        filterout <- c('FC', 'p.value', 'mean.T', 'mean.N')
        # R syntax can become ugly. This is how I acchieve removal by colname
        filterout <-
        -1 * unlist(lapply(filterout, function(x){which(colnames(df) == x)} ) )
        filtered = df[,filterout]
        if ( is.na(probemap) | missing(probemap) ){
            labRow <- srcGeneLookup(rownames(df))
            RowSideColors <- build.side.map(df)
        } else {
            labRow <- unlist(mget(rownames(df), probemap), use.names=F)
            RowSideColors <- build.side.map.plus(df, probemap)
        }
        hv <- heatmap.plus(
            as.matrix(filtered), # PLOTS ONLY THE SAMPLE COLUMNS
            na.rm=TRUE,
            scale="none",
            #scale="row",
            ColSideColors=build.top.map(df, df.s, colormap, colormap2, labels),
            RowSideColors=RowSideColors, 
            col=jet.colors(75),
            key=FALSE,
            #key=T, 
            symkey=FALSE,
            #symkey=T,
            density.info="none",
            trace="none",
            Rowv=TRUE,
            Colv=NA,
            cexRow=0.6,
            cexCol=0.6,
            keysize=2,
            dendrogram=c("none"),
            main = paste(title,
                dim(filtered)[1],
                "Probes; ",
                dim(filtered)[2],
                "Samples"
                ),
            labRow=labRow,
            margin = c(10,10)
        )
        dev.off()
        return(hv)
    }
plotBildDirectionality <-
function(df, df.s, colormap, filename='Heatmap.svg', title='Heatmap', probemap=NA) {
    #Produces a heatmap sorted according to Bild Directionality
    bildExpDir <- build.side.map(df)[,1] # Bild et al. expression directionality
    up <- which(bildExpDir == 'red')
    down <- which(bildExpDir == 'green')
    reorderRows <- c( up, down ) 
    df <- df[reorderRows,]
    svg(filename=filename)
    #colnames I always want removed:
    filterout <- c('FC', 'p.value', 'mean.T', 'mean.N')
    # R syntax can become ugly. This is how I acchieve removal by colname
    filterout <-
    -1 * unlist(lapply(filterout, function(x){which(colnames(df) == x)} ) )
    filtered = df[,filterout]
    hv <- heatmap.plus(
        as.matrix(filtered), # PLOTS ONLY THE SAMPLE COLUMNS
        na.rm=TRUE,
        scale="none",
        ColSideColors=build.top.map(df, df.s, colormap),
        RowSideColors=build.side.map(df), # most differentially expressed genes
        col=jet.colors(75),
        key=FALSE,
        #key=T, 
        symkey=FALSE,
        #symkey=T,
        density.info="none",
        trace="none",
        Rowv=NA,
        Colv=NA,
        cexRow=0.6,
        cexCol=0.6,
        keysize=2,
        dendrogram=c("none"),
        main = paste(title,
            dim(filtered)[1],
            "Probes; ",
            dim(filtered)[2],
            "Samples"
            ),
        probemap=probemap,
        #labCol=NA
        labRow=srcGeneLookup(rownames(df)),
        margin = c(10,10)
    )
    dev.off()
    return(hv)
}

#TODO: figure out how to insert color legend in plots
### Processing Tothill (GSE9891) data:
#reorder columns accordind to the manifest
dfTothill <- tothill.src[,tothill.s.sort$geo_accession]
#TODO: remove COL1A1? #Houtan says this is an outlier (seems consistent in every dataset!)
dfTothill <- dfTothill[-which(rownames(dfTothill) == "1556499_s_at"),]
#dfTothill <- log2(dfTothill)


endometriosis <- rownames(tothill.s.sort)[
    str_detect(tothill.s.sort$characteristics_ch1a, 'Endo')]
otherovarian <- rownames(tothill.s.sort)[
    str_detect(tothill.s.sort$characteristics_ch1a, 'Ser') | 
    str_detect(tothill.s.sort$characteristics_ch1a, 'Adeno') 
    ]

#endometriosis <- which(str_detect(tothill.s.sort$characteristics_ch1a, 'Endo')) 
#otherovarian <- which(!str_detect(tothill.s.sort$characteristics_ch1a, 'Endo'))
mean.n <- apply(dfTothill[,endometriosis], 1, mean, na.rm=T)
mean.t <- apply(dfTothill[,otherovarian], 1, mean, na.rm=T)
p.value <-
   apply(dfTothill, 1, func.list$studentT, s1=endometriosis, s2=otherovarian)
dfTothill <- cbind.data.frame(dfTothill[,endometriosis], dfTothill[,otherovarian])
dfTothill$mean.N <- mean.n
dfTothill$mean.T <- mean.t
dfTothill$p.value <- p.value
dfTothill$FC <- with(dfTothill, mean.T - mean.N)
#Plotting
plotVolcano(dfTothill,
    title="Volcano Plot, Endometriosis vs Other Ovarian Carcinomas (GSE9891)", 
    filename='gseTothill.volcano.svg'
    )
dfTothill.hv <- plotHeatmap(dfTothill, tothill.s, 
    colormap=color.map.tothill.type,
    title = 'Tothill',
    filename = 'gseTothill.heatmap.svg'
    )
dfTothill.hvb <- plotBildDirectionality(dfTothill, tothill.s, color.map.tothill.type,
    title = 'Tothill Bild Sorted',
    filename = 'gseTothill.bild.sorted.svg'
    )
dfTothill.out<-dfTothill[,c('p.value', 'FC')] 
dfTothill.out$id <- dimnames(dfTothill.out)[[1]] #row

## Processing GSE39204 #already log2-converted
df39204 <- gse39204.src[,gse39204.s.sort$geo_accession]    
df39204 <- df39204[-which(rownames(df39204) == "1556499_s_at"),] #COL1A1
#remove <- rownames(
#    gse39204.s[which(
#        gse39204.s$characteristics_ch1a == "histology: serous/clear" |
#        gse39204.s$characteristics_ch1a == "histology: undifferentiated" |
#        gse39204.s$characteristics_ch1a == "histology: carcinosarcoma"
#        
#        ) ,]
#    )
#df39204 <- df39204[,setdiff(colnames(df39204), remove)]
## subsetting part is specific to each dataset (not easily generalizable)
#endometriosis <- which( 
#    str_detect(gse39204.s.sort$characteristics_ch1a, 'histology: clear') | 
#    str_detect(gse39204.s.sort$characteristics_ch1a,'endometrioid') 
#    )
endometriosis <- rownames(gse39204.s.sort)[
    str_detect(gse39204.s.sort$characteristics_ch1a, 'histology: clear') | 
    str_detect(gse39204.s.sort$characteristics_ch1a,'endometrioid') 
    ]
otherovarian <- rownames(gse39204.s.sort)[
    str_detect(gse39204.s.sort$characteristics_ch1a, 'histology: serous') | 
    str_detect(gse39204.s.sort$characteristics_ch1a, 'mucinous') |
    str_detect(gse39204.s.sort$characteristics_ch1a, 'histology: serous low grade') 
    ]
    
## end of specific subsetting 
mean.n <- apply(df39204[,endometriosis], 1, mean, na.rm=T)
mean.t <- apply(df39204[,otherovarian], 1, mean, na.rm=T)
p.value <-
   apply(df39204, 1, func.list$studentT, s1=endometriosis, s2=otherovarian)
df39204 <- cbind.data.frame(df39204[,endometriosis], df39204[,otherovarian])
df39204$mean.N <- mean.n
df39204$mean.T <- mean.t
df39204$p.value <- p.value
df39204$FC <- with(df39204, mean.T - mean.N)
#Plotting
plotVolcano(df39204,
    title="Volcano Plot, Endometriosis vs Other Ovarian Carcinomas (GSE39204)", 
    filename='gse39204.volcano.svg'
    )
df39204.hv <- plotHeatmap(df39204, gse39204.s, 
    colormap=color.map.gse39204.type,
    colormap2=color.map.gse39204.type2,
    labels=c('subtypes', 'Endo vs other'),
    title = 'GSE39204',
    filename = 'gse39204.heatmap.svg'
    )
df39204.hvb <- plotBildDirectionality(df39204, gse39204.s, color.map.gse39204.type,
    title = 'GSE39204 Bild Sorted',
    filename = 'gse39204.bild.sorted.svg'
    )
df39204.out<-df39204[,c('p.value', 'FC')] 
df39204.out$id <- dimnames(df39204.out)[[1]] #row

## Processing GSE29175 #already log2-converted
df29175 <- gse29175.src[,gse29175.s.sort$geo_accession]    
## subsetting part is specific to each dataset (not easily generalizable)
endometriosis <- which( 
    str_detect(gse29175.s.sort$characteristics_ch1a, 'histology: OCCC')  
    ) 
otherovarian <- which(str_detect(gse29175.s.sort$characteristics_ch1a, 'non-OCCC'))
## end of specific subsetting 
mean.n <- apply(df29175[,endometriosis], 1, mean, na.rm=T)
mean.t <- apply(df29175[,otherovarian], 1, mean, na.rm=T)
p.value <-
   apply(df29175, 1, func.list$studentT, s1=endometriosis, s2=otherovarian)
df29175 <- cbind.data.frame(df29175[,endometriosis], df29175[,otherovarian])
df29175$mean.N <- mean.n
df29175$mean.T <- mean.t
df29175$p.value <- p.value
df29175$FC <- with(df29175, mean.T - mean.N)
#Plotting
plotVolcano(df29175,
    title="Volcano Plot, Endometriosis vs Other Ovarian Carcinomas (GSE29175)", 
    filename='gse29175.volcano.svg'
    )
df29175.hv <- plotHeatmap(df29175, gse29175.s, 
    colormap=color.map.gse29175.type,
    title = '29175',
    filename = 'gse29175.heatmap.svg'
    )
df29175.hvb <- plotBildDirectionality(df29175, gse29175.s, color.map.gse29175.type,
    title = '29175 Bild Sorted',
    filename = 'gse29175.bild.sorted.svg'
    )
df29175.out<-df29175[,c('p.value', 'FC')] 
df29175.out$id <- dimnames(df29175.out)[[1]] #row

## Processing GSE29450 #already log2-converted
df29450 <- gse29450.src[,gse29450.s.sort$geo_accession]    
#removing COL1A1 
df29450 <- df29450[-which(rownames(df29450) == "1556499_s_at"),]
## subsetting part is specific to each dataset (not easily generalizable)
endometriosis <- which( 
    str_detect(gse29450.s.sort$characteristics_ch1a, 'clear cell ovarian cancer')  
    ) 
otherovarian <- which(str_detect(gse29450.s.sort$characteristics_ch1a, 'epithelium'))
## end of specific subsetting 
mean.n <- apply(df29450[,endometriosis], 1, mean, na.rm=T)
mean.t <- apply(df29450[,otherovarian], 1, mean, na.rm=T)
p.value <-
   apply(df29450, 1, func.list$studentT, s1=endometriosis, s2=otherovarian)
df29450 <- cbind.data.frame(df29450[,endometriosis], df29450[,otherovarian])
df29450$mean.N <- mean.n
df29450$mean.T <- mean.t
df29450$p.value <- p.value
df29450$FC <- with(df29450, mean.T - mean.N)
#Plotting
plotVolcano(df29450,
    title="Volcano Plot, Endometriosis vs Other Ovarian Carcinomas (GSE29450)", 
    filename='gse29450.volcano.svg'
    )
df29450.hv <- plotHeatmap(df29450, gse29450.s, 
    colormap=color.map.gse29450.type,
    title = 'gse29450',
    filename = 'gse29450.heatmap.svg'
    )
df29450.hvb <- plotBildDirectionality(df29450, gse29450.s, color.map.gse29450.type,
    title = '29450 Bild Sorted',
    filename = 'gse29450.bild.sorted.svg'
    )
df29450.out<-df29450[,c('p.value', 'FC')] 
df29450.out$id <- dimnames(df29450.out)[[1]] #row

## Processing GSE30161 #already log2-converted
df30161 <- gse30161.src[,gse30161.s.sort$geo_accession]
df30161 <- df30161[-which(rownames(df30161) == "1556499_s_at"),] # outlier
df30161 <- df30161[-which(rownames(df30161) == "224321_at"),] # outlier
## subsetting part is specific to each dataset (not easily generalizable)
endometriosis <- rownames(gse30161.s.sort)[
    str_detect(gse30161.s.sort$characteristics_ch1a, 'Clear') | 
    str_detect(gse30161.s.sort$characteristics_ch1a,'Endometrioid') 
    ]
otherovarian <- rownames(gse30161.s.sort)[
    str_detect(gse30161.s.sort$characteristics_ch1a, 'Serous') | 
    str_detect(gse30161.s.sort$characteristics_ch1a, 'Mucinous') 
    ]
# end of specific subsetting 
mean.n <- apply(df30161[,endometriosis], 1, mean, na.rm=T)
mean.t <- apply(df30161[,otherovarian], 1, mean, na.rm=T)
p.value <-
   apply(df30161, 1, func.list$studentT, s1=endometriosis, s2=otherovarian)
df30161 <- cbind.data.frame(df30161[,endometriosis], df30161[,otherovarian])
df30161$mean.N <- mean.n
df30161$mean.T <- mean.t
df30161$p.value <- p.value
df30161$FC <- with(df30161, mean.T - mean.N)
#Plotting
plotVolcano(df30161,
    title="Volcano Plot, Endometriosis vs Other Ovarian Carcinomas (GSE30161)",
    filename='gse30161.volcano.svg'
    )
df30161.hv <- plotHeatmap(df30161, gse30161.s,
    colormap=color.map.gse30161.type,
    title = 'gse30161',
    filename = 'gse30161.heatmap.svg'
    )
df30161.hvb <- plotBildDirectionality(df30161, gse30161.s, color.map.gse30161.type,
    title = '30161 Bild Sorted',
    filename = 'gse30161.bild.sorted.svg'
    )
df30161.out<-df30161[,c('p.value', 'FC')] 
df30161.out$id <- rownames(df30161.out) #row
    
### Processing GSE7305
df7305 <- as.matrix(gse7305.src[,as.character(gse7305.s.sort[,1])])
#df7305 <- df7305[-4,]  #"1556499_s_at probe" an outlier that skews the coloring! Gene COL1A1
df7305 <- df7305[-which(rownames(df7305) == "1556499_s_at"),] 
df7305 <- as.data.frame(df7305);
df7305 <- log2(df7305);
df7305.n <- apply(df7305[,1:10], 1, mean, na.rm=T); 
df7305$mean.N <- df7305.n
df7305.t <- apply(df7305[,11:20], 1, mean, na.rm=T); 
df7305$mean.T <- df7305.t
pv <- apply(df7305, 1, func.list$studentT, s1=c(1:10),s2=c(11:20))
df7305$p.value <- as.numeric(pv)
df7305$FC <- df7305$mean.T - df7305$mean.N
#### Plotting
plotVolcano(df7305,
    title="Volcano Plot, Endometrium/Ovary Disease vs Endometrium-Normal (GSE7305)", 
    filename='gse7305.volcano.svg'
    )
gse7305.hv <- plotHeatmap(df7305, gse7305.s, 
    colormap=color.map.gse7305.type,
    title = 'GSE7305',
    filename = 'gse7305.heatmap.svg'
    )
gse7305.hvb <- plotBildDirectionality(df7305, gse7305.s, color.map.gse7305.type,
    title = 'GSE7305 Bild Sorted',
    filename = 'gse7305.bild.sorted.svg'
    )
df7305.out<-df7305[,c('p.value', 'FC')] 
df7305.out$id <- dimnames(df7305.out)[[1]] #row

### Preprocessing GSE7307
df7307 <- as.matrix(gse7307.src[,as.character(gse7307.s.sort[,1])])
df7307 <- df7307[-which(rownames(df7307) == "1556499_s_at"),] 
df7307 <- as.data.frame(df7307);
gse7307.norm <- as.character(subset(gse7307.s.sort, 
    characteristics_ch1a=="endometrium")$Sample)
gse7307.tum <- as.character(subset(gse7307.s.sort, 
    characteristics_ch1a=="endometrium/ovary")$Sample)
df7307 <- cbind.data.frame(df7307[gse7307.norm], df7307[gse7307.tum])
df7307 <- log2(df7307) #previously only RMA normalized but not log2 transformed
df7307$mean.N <- apply(df7307[,gse7307.norm], 1, mean, na.rm=T)
df7307$mean.T <- apply(df7307[,gse7307.tum], 1, mean, na.rm=T)
pv.7307 <- apply(df7307, 1, func.list$studentT, s1=c(gse7307.norm),s2=c(gse7307.tum))
df7307$p.value <- as.numeric(pv.7307)
df7307$FC <- df7307$mean.T - df7307$mean.N
### Plotting
plotVolcano(df7307,
    title = "Volcano Plot, Endometrium/Ovary vs Endometrium (GSE7307)",
    filename = 'gse7307.volcano.svg'
    )
gse7307.hv <- plotHeatmap(df7307, gse7307.s, color.map.gse7307.type,
    title = 'GSE7307',
    filename = 'gse7307.heatmap.svg'
    )
gse7307.hvb <- plotBildDirectionality(df7307, gse7307.s, color.map.gse7307.type,
    title = 'GSE7307 Bild Sorted',
    filename = 'gse7307.bild.sorted.svg'
    )
df7307.out<-df7307[,c('p.value', 'FC')] 
df7307.out$id <- dimnames(df7307.out)[[1]] #row

### Preprocessing GSE6364
df6364 <- as.matrix(gse6364.src[,as.character(gse6364.s.sort[,1])])
df6364 <- df6364[-which(rownames(df6364) == "1556499_s_at"),] # outlier
df6364 <- as.data.frame(df6364);
df6364 <- log2(df6364);
#normals
df6364.n = df6364[,as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])]
#tumors
df6364.t = df6364[,as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"])]
df6364.n.mean <- apply(df6364.n, 1, mean, na.rm=T) 
df6364.t.mean <- apply(df6364.t, 1, mean, na.rm=T) 
df6364 = cbind(df6364.n, df6364.t)
df6364$mean.T <- df6364.t.mean
df6364$mean.N <- df6364.n.mean
pv.df6364 <- apply(df6364, 1, func.list$studentT, 
    s1=c(as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])),
    s2=c(as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"]))
)
df6364$p.value <- as.numeric(pv.df6364)
df6364$FC <- df6364$mean.T - df6364$mean.N
### Plotting
plotVolcano(df6364,
    title = "Volcano Plot, Endometriosis vs Normal (GSE6364)",
    filename = 'gse6364.volcano.svg'
    )
gse6364.hv <- plotHeatmap(df6364, gse6364.s, color.map.gse6364.type,
    title = 'GSE6364',
    filename = 'gse6364.heatmap.svg'
    )
gse6364.hvb <- plotBildDirectionality(df6364, gse6364.s, color.map.gse6364.type,
    title = 'GSE6364 Bild Sorted',
    filename = 'gse6364.bild.sorted.svg'
    )
df6364.out<-df6364[,c('p.value', 'FC')] 
df6364.out$id <- dimnames(df6364.out)[[1]] #row

##################################
#  FINAL TABLE                   #
##################################
list.names <- c("df6364.out", "df30161.out", "df29450.out", "df29175.out", "dfTothill.out", "df39204.out", "df7307.out", "df7305.out")
list.objects <- list(df6364.out, df30161.out, df29450.out, df29175.out, dfTothill.out, df39204.out, df7307.out, df7305.out)
names(list.objects) <- list.names
for ( i in list.names ) {
    final.table <- merge(Src.signature.out, list.objects[[i]], by.x='ProbeID', by.y='id', all.x=T)
    write.table(final.table, 
                file=paste(i, '.txt', sep=''), 
                quote=F, row.names=F, sep='\t'
                )
    }
###################################
#### NON-AFFY U133 v2.0 DATASETS  #
###################################

### Processing GSE5108
#df5108 <- as.matrix(gse5108.src[,as.character(gse5108.s.sort[,1])])
#df5108 <- as.data.frame(df5108) 
#ColNormals <- which(str_detect(gse5108.s.sort$characteristics_ch1, 'endometriosis')) 
#ColTumors <- which(str_detect(gse5108.s.sort$characteristics_ch1, 'eutopic endometrium')) 
#df5108 <- df5108[,c(ColNormals, ColTumors)]
#df5108$mean.T <- apply(df5108[,ColTumors], 1, mean, na.rm=T)
#df5108$mean.N <- apply(df5108[,ColNormals], 1, mean, na.rm=T)
#df5108$p.value <- as.numeric(
#    apply(df5108, 1, func.list$studentT, s1=ColNormals, s2=ColTumors)
#    )
#df5108$FC <- df5108$mean.T - df5108$mean.N
##plotting
#plotVolcano(df5108,
#    title='GSE5108 SRC genes: Endometriosis vs eutopic endometrium',
#    filename='gse5108.volcano.svg',
#    probemap=hwgcodSYMBOL
#    )
#gse5108.hv <- plotHeatmap(df5108, gse5108.s, color.map.gse5108.type, 
#    probemap=hwgcodSYMBOL,
#    title='GSE5108',
#    filename='gse5108.heatmap.svg'
#    )
#gse5108.hvb <- plotBildDirectionality(df5108, gse5108.s, color.map.gse5108.type,
#    title = 'GSE5108 Bild Sorted',
#    filename = 'gse5108.bild.sorted.svg',
#    probemap = hwgcodSYMBOL
#    )
##TODO append probenames and GeneSymbols for each df
### Processing GSE37837
#df37837 <- gse37837.src[,gse37837.s.sort$geo_accession] 
#ColNormals <- which(str_detect(gse37837.s.sort$characteristics_ch1a, 'eutopic'))
#ColTumors <- which(str_detect(gse37837.s.sort$characteristics_ch1a, 'ectopic'))
#df37837 <- df37837[,c(ColNormals,ColTumors)]
#df37837$mean.T <- apply(df37837[,ColTumors], 1, mean, na.rm=T)
#df37837$mean.N <- apply(df37837[,ColNormals], 1, mean, na.rm=T)
#df37837$p.value <- as.numeric(
#    apply(df37837, 1, func.list$studentT, s1=ColNormals, s2=ColTumors)
#    )
#df37837$FC <- df37837$mean.T - df37837$mean.N
##plotting
#plotVolcano(df37837,
#    filename='gse37837.volcano.svg',
#    title='GSE37837 SRC genes: Endometriosis vs eutopic endometrium',
#    probemap=hgug4112aSYMBOL
#    )
#gse37837.hv <- plotHeatmap(df37837, gse37837.s, color.map.gse37837.type, 
#    probemap=hgug4112aSYMBOL,
#    title='GSE37837',
#    filename='gse37837.heatmap.svg'
#    )
#endometriosis
#endometriosis
