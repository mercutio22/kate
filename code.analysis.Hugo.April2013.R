library(GEOquery)
source('heatmap.plus.R') 
library(heatmap.plus)
source(file="Func_List.r")
library(matlab)
library(stringr)
library(plyr)
library("hgu133plus2.db") # Affymetrix Human Genome U133 Plus 2.0 Array annotation data
library("hwgcod.db") # GE/Amersham Codelink Human Whole Genome Bioarray 
options(stringsAsFactors=FALSE)

###tothill data
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
#sample.info.gse7307 <- read.delim(file="GSE7307_GEO_Sample_Info.csv")
#gse7307.disease <- (subset(sample.info.gse7307, Disease.type!="Normal"))
#gse7307.endo.normal <- (subset(sample.info.gse7307, Tissue.Cell.Line..C.=="endometrium"))
#gse7307.s <- rbind(gse7307.disease, gse7307.endo.normal)
#
##TODO - get more datasets
##GSE37837
##GSE23339
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
#gse5108.s = phenoData(gse5108$GSE5108_series_matrix.txt.gz)@data #Holy Cow, 11 sample kinds.
# all(dimnames(gse5108.s)[[1]] == names(gse5108.d))
#
#dataframe <- NULL
#for(i in 1:length(gse7307.s$Sample)){
#	
#	if(is.null(dataframe)){
#		test <- getGEO(as.character(gse7307.s[i,1]))
#		dataframe <- Table(test)
#		dataframe[,2] <- as.numeric(dataframe[,2])
#		dimnames(dataframe)[[2]] <- c("ID_REF",as.character(gse7307.s[i,1]))
#		#print("Entrou IF") #Debug
#	}
#	
#	else{
#		test <- getGEO(as.character(gse7307.s[i,1]))
#		test.d <- Table(test)
#		test.d[,2] <- as.numeric(test.d[,2])
#		dimnames(test.d)[[2]] <- c("ID_REF",as.character(gse7307.s[i,1]))
#		dataframe = merge(dataframe, test.d, by.x="ID_REF", by.y="ID_REF")
#		#print("ELSE") #Debug
#	}
#
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
#
#save(
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
david.src = read.table('David.Src.txt', header=T) #maybe I should use the corresponding R/BioC Annotation package for consistency
david.src$To=toupper(as.character(david.src$To))
david.src = david.src[,c(2,1)] #inverting column order because I am going to add probe names of other platforms
colnames(david.src) <- c('GeneSymbol', 'affyProbe')
# getting illumina IDs:
# davidIllumina <- read.table('david_affy_illumina.txt', header=T, sep='\t') #submitted david.src$From to DAVID conversion tool and got Illumina IDs 
#update Src.signature with updated gene symbols
indices = match(david.src$affyProbe, Src.signature$ProbeID)
Src.signature$GeneSymbol[indices] <- david.src$GeneSymbol 
## confirming this with a bimap ###
## unlist(mget(Src.signature$ProbeID, hgu133plus2SYMBOL), use.names=F)
## Mostly agrees.
src = Src.signature #some of houtan's code reference this variable name

# used to subset Codelink probeIDs:
tmp <- unlist(
    mget(david.src$GeneSymbol, revmap(hwgcodSYMBOL), ifnotfound=NA),
    use.names=F
    )
codelinkSrc <- unique(tmp[!is.na(tmp)]) # GE/Amersham probe IDs corresponding to affy U133 v2.0
rm(tmp)

##get probe specific to SRC signature
gse6008.src.probe <- intersect(src[,1], dimnames(gse6008.d)[[1]])
tothill.src.probe <- intersect(src[,1], dimnames(tothill.d)[[1]])
gse7305.src.probe <- intersect(src[,1], dimnames(gse7305.d)[[1]])
gse7307.src.probe <- intersect(src[,1], dimnames(gse7307.d)[[1]])
gse6364.src.probe <- intersect(src[,1], dimnames(gse6364.d)[[1]])
##The probenames should be prefixed by 'GE' but for some reason its not. WARNING! Verify why.
rownames(gse5108.d) <- paste("GE",  rownames(gse5108.d), sep='') 
gse5108.src.probe <- intersect(codelinkSrc, dimnames(gse5108.d)[[1]])
#
##subset data to only SRC specific and create dataframe. 
tothill.src <- tothill.d[tothill.src.probe,]; tothill.src <- as.data.frame(tothill.src)
gse6008.src <- gse6008.d[gse6008.src.probe,]; gse6008.src <- as.data.frame(gse6008.src)
gse6364.src <- gse6364.d[gse6364.src.probe,]; gse6364.src <- as.data.frame(gse6364.src)
gse7307.src <- gse7307.d[gse7307.src.probe,]; gse7307.src <- as.data.frame(gse7307.src)
gse7305.src <- gse7305.d[gse7305.src.probe,]; gse7305.src <- as.data.frame(gse7305.src)
gse5108.src <- as.data.frame(gse5108.d[gse5108.src.probe,])


#WHat does this do? clueless..
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
tothill.nn <- strsplit(as.character(tothill.s$characteristics_ch1), " : ")
tothill.nnn <- as.data.frame(as.character(unlist(lapply((tothill.nn), "[", 2) )))
tothill.s$characteristics_ch1a <- tothill.nnn[,1]

gse6008.nn <- strsplit(as.character(gse6008.s$characteristics_ch1), ": ")
gse6008.nnn <- as.data.frame(as.character(unlist(lapply((gse6008.nn), "[", 2) )))
gse6008.s$characteristics_ch1a <- gse6008.nnn[,1]

#further code expects every manifest dataframe (df.s) to contain a characteristics_ch1 column
gse6364.s$characteristics_ch1a <- gse6364.s$characteristics_ch1

### endometrium is normal & endometrium/ovary is diseased
gse7307.s$characteristics_ch1a <- gse7307.s$Tissue.Cell.Line..C. 

gse7305.nn <- strsplit(as.character(gse7305.s$characteristics_ch1), " ")
gse7305.nnn <- as.data.frame(as.character(unlist(lapply((gse7305.nn), "[", 1) )))
gse7305.s$characteristics_ch1a <- gse7305.nnn[,1]

gse5108.s$characteristics_ch1a <- gse5108.s$characteristics_ch1

## why do we need to do this? 
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

## ColSideColors
color.map.tothill.type <- function(mol.biol) {
	if (mol.biol=="Fallopian") "green" #green
	else if (mol.biol=="Ovary") "#red" #red
	else if (mol.biol=="Peritoneum") "black"
	else "#FFFFFF" #white
} #
color.map.gse6008.type <- function(mol.biol) {
	if (mol.biol=="Clear_Cell") "green" #green
	else if (mol.biol=="Endometrioid") "red" #red
	else if (mol.biol=="Mucinous") "black"
	else if (mol.biol=="Serous") "yellow"
	else "#FFFFFF"
} #
color.map.gse7307.type <- function(mol.biol) {
	if (mol.biol=="myometrium") "green" #green
	else if (mol.biol=="endometrium/ovary") "red" 
	else if (mol.biol=="endometrium") "#9C897E" #grey
	else "#FFFFFF"
} #
color.map.gse7305.type <- function(mol.biol) {
	if (mol.biol=="Endometrium/Ovary-Disease") "red" 
	else if (mol.biol=="Endometrium-Normal") "#9C897E" #grey
	else "#FFFFFF"
} #
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

#rowsidecolors 
bildColors.colorize <- function(abbrv){ 
	if (abbrv == 'UP.BILD') 'red'
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
	### takes dataframe, merges with the david.src signature and creates RowSideColors for our plots
	## df is the Bild SRC subsetted GSE dataframe, with samples, FC, mean.T, mean.N and p.value
    ## probe2geneSymbol.db is a biocondutor "bimap" annotation object such as hgu133plus2GENENAME 
    df$probeID <- rownames(df)
    df$GeneSymbol <- unlist( mget(df$probeID, probemap, ifnotfound=NA), use.names=F )
    #having the geneSymbol, lookup whats the gene expression directionality in the Src.signature
    indices <- match(df$GeneSymbol, Src.signature$GeneSymbol)
    df$direction <- Src.signature$direction[indices]
    df$direction[is.na(df$direction)] <- 'unknown' #DAVID makes ZN12 disappear! DISCUSS with Houtan how we should proceed.
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

build.top.map = function(df, df.s, df.colormap.function){
	#builds the heatmap top rows for helping group visualization
	## df is the full dataframe, with FC and P values
	## df.colormap.function maps metadata with colors
	cols = as.matrix(colnames(df[1:(ncol(df)-4)]))
	if ('Sample' %in% colnames(df.s)){
	df.src.type = apply(cols,
		1,
		func.list$vlookup,
		df.s[,c("Sample","characteristics_ch1a")],
		"characteristics_ch1a")
	}
	else { df.src.type = apply(cols,
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
	ColSideColors = cbind(ColSideColors, ColSideColors) #2columns make it thicker 
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

#build.top.map.plus <- function(df, df.s.sort, df.colormap.function) {
#    df <- removeColumns(df)    
#    colors <- unlist(lapply(df.s.sort$characteristics_ch1a, df.colormap.function)) 
#	ColSideColors = matrix(as.character(c(colors)), 
#			nrow = length(colors),
#			ncol = 1)
#    colors <- cbind(ColSideColors, ColSideColors)
#    return(colors)
#}

#srcGeneLookup <- function(probes) {
#    return(
#        as.character(unlist(
#            lapply(probes, func.list$vlookup, david.src, 'GeneSymbol')
#            )
#        )
#    )
#}
srcGeneLookup <- function(affyP) {
    indices <- unlist(lapply(david.src$affyProbe, match, affyP), use.names=F)
    genes <- (david.src$GeneSymbol[indices])
    genes[is.na(genes)] <- 'unknown'
}   

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make volcano plot plus
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
volcano.plus<-
function(x,fold.change.col,pv.col, title.plot, cut.line, fold.cut1, 
    fold.cut2, pv.adj.col,ncolors=1, text=NA, angle=-45, probemap=hgu133plus2SYMBOL) {
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

            #points(M[upr], -log10(p[upr]), col="red")
            #points(M[dwr], -log10(p[dwr]), col="green")

            abline(h= -log10(cut.line), lty=3, lwd=1)
            abline(v= fold.cut1, lty=3, lwd=1, col="black")
            abline(v= fold.cut2, lty=3, lwd=1, col="black")
        }
        if( !all(is.na(text)) ){
            probes = rownames(df)[text]
            text(M[text],
                -log10(p[text]), 
                #labels=probes,
                #My hack for labeling kate's data according to gene symbols
                labels=srcGeneLookup(probes),
#                labels=unlist(
#                    mget(probes, probemap, ifnotfound=NA), 
#                    use.names=F
#                    ),
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
plotVolcano <- function(df, filename='volcano.svg', 
    title='Volcano Plot', 
    probemap=hgu133plus2SYMBOL) {
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
        angle = 45
        ) 
    dev.off()
} 
plotHeatmap <-
function(df, df.s, colormap, filename='Heatmap.svg', title='Heatmap', probemap=hgu133plus2SYMBOL) {
    svg(filename=filename)
    #colnames I always want removed:
    filterout <- c('FC', 'p.value', 'mean.T', 'mean.N')
    # R syntax can become ugly. This is how I acchieve removal by colname
    filterout <-
    -1 * unlist(lapply(filterout, function(x){which(colnames(df) == x)} ) )
    filtered = df[,filterout]
    labRow=srcGeneLookup(rownames(df)),
    labRow=unlist(mget(rownames(df), probemap), use.names=F),
    hv <- heatmap.plus(
        as.matrix(filtered), # PLOTS ONLY THE SAMPLE COLUMNS
        na.rm=TRUE,
        scale="none",
        #scale="row",
        ColSideColors=build.top.map(df, df.s, colormap),
        RowSideColors=build.side.map.plus(df, probemap), # most differentially expressed genes
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
        #labCol=NA
        #labRow=srcGeneLookup(rownames(df)),
        labRow=unlist(mget(rownames(df), probemap), use.names=F),
        margin = c(10,10)
    )
    dev.off()
    return(hv)
}
plotBildDirectionality <-
function(df, df.s, colormap, filename='Heatmap.svg', title='Heatmap') {
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
        #labCol=NA
        labRow=srcGeneLookup(rownames(df)),
        margin = c(10,10)
    )
    dev.off()
    return(hv)
}

#cc.col.merged <- rbind(cc.col.t, cc.col.gse6008)
#colnames(cc.col.merged) = c("Tothill & GSE6008","Tothill & GSE6008")
#cc.col.merged <- cbind(cc.col.merged, cc.col.merged)
#png(filename = "tothill.src.png", bg="white", res=300, width=3000, height=3000)
#tothill.hv<-heatmap.plus(
#		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
#		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
#		as.matrix(tothill.src[,tothill.s.sort[,1]]), 
#		#temp.cc,
#		na.rm=TRUE,
#		scale="none",
#		#RowSideColor=probe.cc,
#		ColSideColors=cc.col.t,
#		col=jet.colors(75),
#		#key=FALSE,
#		symkey=FALSE,
#		density.info="none",
#		trace="none",
#		Rowv=TRUE,
#		Colv=NA,
#		cexRow=1,
#		cexCol=0.6,
#		keysize=1,
#		dendrogram=c("none"),
#		main = paste("Tothill; SRC only (73 probes) ",
#			dim(tothill.src)[1],
#			"Probes; ",
#			dim(tothill.src)[2],
#			"Samples"
#			),
#		#labCol=NA
#		labRow=NA
#)
#dev.off()
#png(filename = "gse6008.src.png", bg="white", res=300, width=3000, height=3000)
#gse6008.hv<-heatmap.plus(
#		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
#		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
#		as.matrix(gse6008.src[,gse6008.s.sort[,1]]), 
#		#temp.cc,
#		na.rm=TRUE,
#		scale="none",
#		#RowSideColor=probe.cc,
#		ColSideColors=cc.col.gse6008,
#		col=jet.colors(75),
#		#key=FALSE,
#		symkey=FALSE,
#		density.info="none",
#		trace="none",
#		Rowv=TRUE,
#		Colv=NA,
#		cexRow=1,
#		cexCol=0.6,
#		keysize=1,
#		dendrogram=c("none"),
#		main = paste("gse6008; SRC only (73 probes) ",
#			dim(gse6008.src)[1],"Probes; ",
#			dim(gse6008.src)[2],
#			"Samples"
#			),
#		#labCol=NA
#		labRow=NA
#		)
#dev.off()

#this table compares our analysis with Bilds's
#tcomparison <- table(subset(src.sig, label == "Significant")$direction.y, subset(src.sig, label == "Significant")$direction.x)

### Preprocessing GSE7305
tttt <- as.matrix(gse7305.src[,as.character(gse7305.s.sort[,1])])
#tttt <- tttt[-4,]  #"1556499_s_ati probe" an outlier that skews the coloring!
tttt <- tttt[-which(rownames(tttt) == "1556499_s_at"),] 
tttt <- as.data.frame(tttt);
tttt <- log2(tttt);
tttt.n <- apply(tttt[,1:10], 1, mean, na.rm=T); 
tttt$mean.N <- tttt.n
tttt.t <- apply(tttt[,11:20], 1, mean, na.rm=T); 
tttt$mean.T <- tttt.t
pv <- apply(tttt, 1, func.list$studentT, s1=c(1:10),s2=c(11:20))
tttt$p.value <- as.numeric(pv)
tttt$FC <- tttt$mean.T - tttt$mean.N

#### Plotting
plotVolcano(tttt,
    title="Volcano Plot, Endometrium/Ovary Disease vs Endometrium-Normal (GSE7305)", 
    filename='gse7305.volcano.svg'
    #filename='gse7305.volcano.png'
    )
#gse7305.hv <- plotHeatmap(tttt, gse7305.s, color.map.gse7305.type,
#    title = 'GSE7305',
#    filename = 'gse7305.heatmap.svg'
#    )
#gse7305.hvb <- plotBildDirectionality(tttt, gse7305.s, color.map.gse7305.type,
#    title = 'GSE7305 Bild Sorted',
#    filename = 'gse7305.bild.sorted.svg'
#    )
#
#### Preprocessing GSE7307
#tt <- as.matrix(gse7307.src[,as.character(gse7307.s.sort[,1])])
#tt <- tt[-which(rownames(tt) == "1556499_s_at"),] 
#tt <- as.data.frame(tt);
#gse7307.norm <- as.character(subset(gse7307.s.sort, 
#    characteristics_ch1a=="endometrium")$Sample)
#gse7307.tum <- as.character(subset(gse7307.s.sort, 
#    characteristics_ch1a=="endometrium/ovary")$Sample)
#tt.t <- apply(tt[,gse7307.tum], 1, mean, na.rm=T); tt$mean.T <- tt.t
#tt$FC <- tt$mean.T - tt$mean.N
#### Plotting
#plotVolcano(tt,
#    title = "Volcano Plot, Endometrium/Ovary vs Endometrium (GSE7307)",
#    filename = 'gse7307.volcano.svg'
#    )
#gse7307.hv <- plotHeatmap(tt, gse7307.s, color.map.gse7307.type,
#    title = 'GSE7307',
#    filename = 'gse7307.heatmap.svg'
#    )
#gse7307.hvb <- plotBildDirectionality(tt, gse7307.s, color.map.gse7307.type,
#    title = 'GSE7307 Bild Sorted',
#    filename = 'gse7307.bild.sorted.svg'
#    )
#
#### Preprocessing GSE6364
#ttt <- as.matrix(gse6364.src[,as.character(gse6364.s.sort[,1])])
#ttt <- ttt[-which(rownames(ttt) == "1556499_s_at"),] # outlier
#ttt <- as.data.frame(ttt);
#ttt <- log2(ttt);
##normals
#ttt.n = ttt[,as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])]
##tumors
#ttt.t = ttt[,as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"])]
#ttt.n.mean <- apply(ttt.n, 1, mean, na.rm=T) 
#ttt.t.mean <- apply(ttt.t, 1, mean, na.rm=T) 
#ttt = cbind(ttt.n, ttt.t)
#ttt$mean.T <- ttt.t.mean
#ttt$mean.N <- ttt.n.mean
#pv.ttt <- apply(ttt, 1, func.list$studentT, s1=c(as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])),s2=c(as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"])))
#ttt$p.value <- as.numeric(pv.ttt)
#ttt$FC <- ttt$mean.T - ttt$mean.N
#### Plotting
#plotVolcano(ttt,
#    title = "Volcano Plot, Endometriosis vs Normal (GSE6364)",
#    filename = 'gse6364.volcano.svg'
#    )
#gse6364.hv <- plotHeatmap(ttt, gse6364.s, color.map.gse6364.type,
#    title = 'GSE6364',
#    filename = 'gse6364.heatmap.svg'
#    )
#gse6364.hvb <- plotBildDirectionality(ttt, gse6364.s, color.map.gse6364.type,
#    title = 'GSE6364 Bild Sorted',
#    filename = 'gse6364.bild.sorted.svg'
#    )

## Preprocessing GSE5108
df5108 <- as.matrix(gse5108.src[,as.character(gse5108.s.sort[,1])])
df5108 <- as.data.frame(df5108) 
ColNormals <- which(str_detect(gse5108.s.sort$characteristics_ch1, 'endometriosis')) 
ColTumors <- which(str_detect(gse5108.s.sort$characteristics_ch1, 'eutopic endometrium')) 
df5108 <- df5108[,c(ColNormals, ColTumors)]
df5108$mean.T <- apply(df5108[,ColTumors], 1, mean, na.rm=T)
df5108$mean.N <- apply(df5108[,ColNormals], 1, mean, na.rm=T)
df5108$p.value <- as.numeric(
    apply(df5108, 1, func.list$studentT, s1=ColNormals, s2=ColTumors)
    )
df5108$FC <- df5108$mean.T - df5108$mean.N
#plotting
plotVolcano(df5108,
    title='GSE5108 SRC genes: Endometriosis vs eutopic endometrium',
    filename='gse5108.volcano.svg'
#    filename='gse5108.volcano.png'
    )
gse5108.hv <- plotHeatmap(df5108, gse5108.s, color.map.gse5108.type, 
    probemap=hwgcodSYMBOL,
    title='GSE5108',
    filename='gse5108.heatmap.svg'
    )
