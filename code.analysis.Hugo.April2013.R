library(GEOquery)
source('heatmap.plus.R') #Houtan's heatmat function??
library(heatmap.plus)
source(file="Func_List.r")
library(matlab)

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
###GSE7305 - 20 samples 10 endo/ovary and 10 endo normal
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
##GSE3151 - Bild et al. Oncogene Signature Dataset
##GSE5108

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
#all(dimnames(gse6364.s)[[1]]==names(gse6364.d))
#
#
#save(gse6364.d, gse6364.s, gse7307.d, gse7307.s, gse7305.d, gse7305.s, tothill.d, tothill.s, gse6008.d, gse6008.s, file="GEOfiles_associated_withENDO.downloaded.April2013.rda")


###start from here#####
setwd("/data/htorres/kate")
load(file="GEOfiles_associated_withENDO.downloaded.April2013.rda")


## New Src signature data frame
Src.signature = read.table('src_signature.txt', sep='\t') # src_signature.txt from suplemental table1
colnames(Src.signature) = c('ProbeID', 'GeneSymbol', 'Description', 'LocusLink', 'FoldChange')
Src.sig$direction <- NA
Src.signature$direction[Src.signature$FoldChange<=1] <- c("DN.BILD") # just to denote this info is not from our analysis (BILD) 
Src.signature$direction[Src.signature$FoldChange>1] <- c("UP.BILD")

##src signature #Bild et al. 2006
src <- read.delim(file="bild.src.txt") #where did this file come from? Maybe GSE3151
length(intersect(src[,1], dimnames(tothill.d)[[1]]))
length(intersect(src[,1], dimnames(gse6008.d)[[1]]))
length(intersect(src[,1], dimnames(gse7305.d)[[1]]))
length(intersect(src[,1], dimnames(gse7307.d)[[1]]))
length(intersect(src[,1], dimnames(gse6364.d)[[1]]))
#get probe specific to SRC signature
gse6008.src.probe <- intersect(src[,1], dimnames(gse6008.d)[[1]])
tothill.src.probe <- intersect(src[,1], dimnames(tothill.d)[[1]])
gse7305.src.probe <- intersect(src[,1], dimnames(gse7305.d)[[1]])
gse7307.src.probe <- intersect(src[,1], dimnames(gse7307.d)[[1]])
gse6364.src.probe <- intersect(src[,1], dimnames(gse6364.d)[[1]])
#subset data to only SRC specific and create dataframe. 
tothill.src <- tothill.d[tothill.src.probe,]; tothill.src <- as.data.frame(tothill.src)
gse6008.src <- gse6008.d[gse6008.src.probe,]; gse6008.src <- as.data.frame(gse6008.src)
gse6364.src <- gse6364.d[gse6364.src.probe,]; gse6364.src <- as.data.frame(gse6364.src)
gse7307.src <- gse7307.d[gse7307.src.probe,]; gse7307.src <- as.data.frame(gse7307.src)
gse7305.src <- gse7305.d[gse7305.src.probe,]; gse7305.src <- as.data.frame(gse7305.src)

sort.data.frame <- function(x, key, ...) {
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}


##cleanup sample manifest #metadata)
tothill.nn <- strsplit(as.character(tothill.s$characteristics_ch1), " : ")
tothill.nnn <- as.data.frame(as.character(unlist(lapply((tothill.nn), "[", 2) )))
tothill.s$characteristics_ch1a <- tothill.nnn[,1]


gse6008.nn <- strsplit(as.character(gse6008.s$characteristics_ch1), ": ")
gse6008.nnn <- as.data.frame(as.character(unlist(lapply((gse6008.nn), "[", 2) )))
gse6008.s$characteristics_ch1a <- gse6008.nnn[,1]

gse6364.s$characteristics_ch1a <- gse6364.s$characteristics_ch1

## endometrium is normal & endometrium/ovary is diseased
gse7307.s$characteristics_ch1a <- gse7307.s$Tissue.Cell.Line..C. 

gse7305.nn <- strsplit(as.character(gse7305.s$characteristics_ch1), " ")
gse7305.nnn <- as.data.frame(as.character(unlist(lapply((gse7305.nn), "[", 1) )))
gse7305.s$characteristics_ch1a <- gse7305.nnn[,1]


#This repetitive block I will refactor into a function.
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

#tothill.gse6008.src <- cbind(tothill.d[gse6008.src.probe,], gse6008.d[gse6008.src.probe,]); tothill.gse6008.src <- as.data.frame(tothill.gse6008.src) ## combining two dataset on common src probe

#tothill.gse6008.src.sort <- tothill.gse6008.src[,c(tothill.s.sort[,1],gse6008.s.sort[,1])]


## ColSideColors
color.map.tothill.type <- function(mol.biol) {
	if (mol.biol=="Fallopian") "green" #green
	else if (mol.biol=="Ovary") "#701C22" #red
	else if (mol.biol=="Peritoneum") "black"
	else "#FFFFFF" #white
} #
color.map.gse6008.type <- function(mol.biol) {
	if (mol.biol=="Clear_Cell") "green" #green
	else if (mol.biol=="Endometrioid") "#701C22" #red
	else if (mol.biol=="Mucinous") "black"
	else if (mol.biol=="Serous") "yellow"
	else "#FFFFFF"
} #
color.map.gse7307.type <- function(mol.biol) {
	if (mol.biol=="myometrium") "green" #green
	else if (mol.biol=="endometrium/ovary") "#701C22" #red
	else if (mol.biol=="endometrium") "#9C897E" #grey
	else "#FFFFFF"
} #
color.map.gse7305.type <- function(mol.biol) {
	if (mol.biol=="Endometrium/Ovary-Disease") "#701C22" #green
	else if (mol.biol=="Endometrium-Normal") "#9C897E" #grey
	else "#FFFFFF"
} #
color.map.gse6364.type <- function(mol.biol) {
	if (mol.biol=="Early Secretory Phase Endometriosis") "701C22" #red
	else if (mol.biol=="Proliferative Phase Endometriosis") "701C22" #red
	else if (mol.biol=="Mid Secretory Phase Endometriosis") "701C22" #red
	else if (mol.biol=="Mid Secretory Phase Normal") "#9C897E" 
	else if (mol.biol=="Proliferative Phase Normal") "#9C897E" 
	else if (mol.biol=="Early Secretory Phase Normal") "#9C897E" 
	else "#FFFFFF"
} #

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
	df$direction = NA
	df$direction[df$FC<0]  <- c("DN.KL") 
	df$direction[df$FC>=0] <- c("UP.KL")
	df$id = dimnames(df)[[1]]
	df$significance[df$p.value<=0.05] <- c("Significant")
	src.sig <- merge(df, Src.signature, by.x = 'id', by.y = 'ProbeID', all.x = T) 
	#I want the white color to mark the statistically insignificant cases, so:
	src.sig$direction.x[which(is.na(src.sig$significance))] <- 'insignificant'
	bildColors = unlist(lapply(src.sig$direction.y, bildColors.colorize))
	klColors = unlist(lapply(src.sig$direction.x, klColors.colorize))
	colors = cbind(bildColors,klColors)
	colnames(colors) = c('Bild et al.', 'Lawrenson et al.')
	return(colors)
}

## WORKING HERE ##TODO: This sux. Should need fewer parameters
build.top.map = function(df, df.s, df.colormap.function){
	#builds the heatmap top rows for helping group visualization
	## df is the full dataframe, with FC and P values
	## df.colormap.function maps metadata with colors
	cols = as.matrix(colnames(df[1:(ncol(df)-4)]))
	df.src.type = apply(cols,1,func.list$vlookup,df.s[,c("Sample","characteristics_ch1a")],"characteristics_ch1a")
	df.src.type = as.character(df.src.type)
	df.src.type[is.na(df.src.type)] = c('0')
	cc.df.src.type = unlist(lapply(df.src.type, df.colormap.function))
	ColSideColors = matrix(as.character(c(cc.df.src.type)), 
			nrow = length(cc.df.src.type),
			ncol = 1)
	ColSideColors = cbind(ColSideColors, ColSideColors) #a trick to make the top plot part thicker 
	return(ColSideColors)
}


####################################################### THIS BLOCK WILL BE GONE #################
#################################################################################################
#tothill.type
tothill.src.type<-apply(as.matrix(names(tothill.src[,tothill.s.sort[,1]])),
	1,
	func.list$vlookup,
	tothill.s[,c("geo_accession","characteristics_ch1a")],
	"characteristics_ch1a")
tothill.src.type<-as.character(tothill.src.type)
tothill.src.type[is.na(tothill.src.type)]<-c("0")
cc.tothill.src.type<-unlist(lapply(tothill.src.type, color.map.tothill.type))
cc.col.t<-matrix(as.character(c(cc.tothill.src.type)), nrow=285,ncol=1)
colnames(cc.col.t) = c("Tothill")
cc.col.t <- cbind(cc.col.t, cc.col.t)

#gse6008.type
gse6008.src.type<-apply(as.matrix(names(gse6008.src[,gse6008.s.sort[,1]])),1,func.list$vlookup,gse6008.s[,c("geo_accession","characteristics_ch1a")],"characteristics_ch1a")
gse6008.src.type<-as.character(gse6008.src.type)
gse6008.src.type[is.na(gse6008.src.type)]<-c("0")
cc.gse6008.src.type<-unlist(lapply(gse6008.src.type, color.map.gse6008.type))

cc.col.gse6008<-matrix(as.character(c(cc.gse6008.src.type)), nrow=103,ncol=1)
colnames(cc.col.gse6008) = c("GSE6008")
cc.col.gse6008 <- cbind(cc.col.gse6008, cc.col.gse6008)

#gse6364.type ##PROBLEM: the resulting dataframe is 37 by 2 - heatmap expects 41 by 2
gse6364.src.type<-apply(as.matrix(names(gse6364.src[,gse6364.s.sort[,1]])),1,func.list$vlookup,gse6364.s[,c("geo_accession","characteristics_ch1a")],"characteristics_ch1a")
gse6364.src.type<-as.character(gse6364.src.type)
gse6364.src.type[is.na(gse6364.src.type)]<-c("0")
cc.gse6364.src.type<-unlist(lapply(gse6364.src.type, color.map.gse6364.type))

cc.col.gse6364<-matrix(as.character(c(cc.gse6364.src.type)), 
	nrow=37,
	#nrow=41,
	ncol=1)
colnames(cc.col.gse6364) = c("GSE6364")
cc.col.gse6364 <- cbind(cc.col.gse6364, cc.col.gse6364)

#gse7305.type
gse7305.src.type<-apply(
	as.matrix(names(gse7305.src[,gse7305.s.sort[,1]])),
	1,
	func.list$vlookup,
	gse7305.s[,c("geo_accession","characteristics_ch1a")],
	"characteristics_ch1a"
	)
gse7305.src.type<-as.character(gse7305.src.type)
gse7305.src.type[is.na(gse7305.src.type)]<-c("0")
cc.gse7305.src.type<-unlist(lapply(gse7305.src.type, color.map.gse7305.type))

cc.col.gse7305<-matrix(as.character(c(cc.gse7305.src.type)), nrow=20,ncol=1)
colnames(cc.col.gse7305) = c("GSE7305")
cc.col.gse7305 <- cbind(cc.col.gse7305, cc.col.gse7305)

#gse7307.type
gse7307.src.type<-apply(
	as.matrix(names(gse7307.src[,as.character(gse7307.s.sort[,1])])),
	1,
	func.list$vlookup,
	gse7307.s[,c("Sample","characteristics_ch1a")],
	"characteristics_ch1a"
	)
gse7307.src.type<-as.character(gse7307.src.type)
gse7307.src.type[is.na(gse7307.src.type)]<-c("0")
cc.gse7307.src.type<-unlist(lapply(gse7307.src.type, color.map.gse7307.type))

cc.col.gse7307<-matrix(as.character(c(cc.gse7307.src.type)), nrow=134,ncol=1)
colnames(cc.col.gse7307) = c("GSE7307")
cc.col.gse7307 <- cbind(cc.col.gse7307, cc.col.gse7307)
#######################################################################################
################################GONE###################################################


#cc.col.merged <- rbind(cc.col.t, cc.col.gse6008)
#colnames(cc.col.merged) = c("Tothill & GSE6008","Tothill & GSE6008")
#cc.col.merged <- cbind(cc.col.merged, cc.col.merged)

png(filename = "tothill.src.png", bg="white", res=300, width=3000, height=3000)
tothill.hv<-heatmap.plus(
		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
		as.matrix(tothill.src[,tothill.s.sort[,1]]), 
		#temp.cc,
		na.rm=TRUE,
		scale="none",
		#RowSideColor=probe.cc,
		ColSideColors=cc.col.t,
		col=jet.colors(75),
		#key=FALSE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=TRUE,
		Colv=NA,
		cexRow=1,
		cexCol=0.6,
		keysize=1,
		dendrogram=c("none"),
		main = paste("Tothill; SRC only (73 probes) ",
			dim(tothill.src)[1],
			"Probes; ",
			dim(tothill.src)[2],
			"Samples"
			),
		#labCol=NA
		labRow=NA
)
dev.off()
png(filename = "gse6008.src.png", bg="white", res=300, width=3000, height=3000)
gse6008.hv<-heatmap.plus(
		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
		as.matrix(gse6008.src[,gse6008.s.sort[,1]]), 
		#temp.cc,
		na.rm=TRUE,
		scale="none",
		#RowSideColor=probe.cc,
		ColSideColors=cc.col.gse6008,
		col=jet.colors(75),
		#key=FALSE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=TRUE,
		Colv=NA,
		cexRow=1,
		cexCol=0.6,
		keysize=1,
		dendrogram=c("none"),
		main = paste("gse6008; SRC only (73 probes) ",
			dim(gse6008.src)[1],"Probes; ",
			dim(gse6008.src)[2],
			"Samples"
			),
		#labCol=NA
		labRow=NA
		)
dev.off()

## This block is new ####
tttt <- as.matrix(gse7305.src[,as.character(gse7305.s.sort[,1])])
tttt <- tttt[-4,]  #"it refers to the 1556499_s_ati probe" an outlier that skews the coloring!
tttt <- as.data.frame(tttt);
tttt <- log2(tttt);
tttt.n <- apply(tttt[,1:10], 1, mean, na.rm=T); tttt$mean.N <- tttt.n
tttt.t <- apply(tttt[,11:20], 1, mean, na.rm=T); tttt$mean.T <- tttt.t
pv <- apply(tttt, 1, func.list$studentT, s1=c(1:10),s2=c(11:20))
tttt$p.value <- as.numeric(pv)
tttt$FC <- tttt$mean.T - tttt$mean.N
func.list$volcano(tttt, fold.change.col = "FC", pv.col = "p.value", title.plot= "Volcano Plot, Endometrium vs Normal (GSE7305)", cut.line = -log10(0.05), foldcut = 0.5) 
#this table compares our analysis with Bilds's
tcomparison <- table(subset(src.sig, label == "Significant")$direction.y, subset(src.sig, label == "Significant")$direction.x)
png(filename = "gse7305.src.png", bg="white", res=300, width=3000, height=3000)
gse7305.hv<-heatmap.plus(
		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
		#as.matrix(gse7305.src[,gse7305.s.sort[,1]]), 
		as.matrix(tttt[,1:20]), # PLOTS ONLY THE SAMPLE COLUMNS
		#temp.cc,
		na.rm=TRUE,
		scale="none",
		#RowSideColor=probe.cc,
		#ColSideColors=cc.col.gse7305,
		ColSideColors=build.top.map(gse7305.src,
			gse7305.s,
			gse7305.s.sort, 
			label=GSE7305,
			color.map.gse7305.type
			),
		RowSideColors=build.side.map(tttt), # most differentially expressed genes
		col=jet.colors(75),
		#key=FALSE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=TRUE,
		Colv=NA,
		cexRow=1,
		cexCol=0.6,
		keysize=1,
		dendrogram=c("none"),
		main = paste("gse7305; SRC only (73 probes) ",
			dim(gse7305.src)[1],
			"Probes; ",
			dim(gse7305.src)[2],
			"Samples"
			),
		#labCol=NA
		labRow=NA
)
dev.off()
png(filename = "gse7307.src.png", bg="white", res=300, width=3000, height=3000)
tt <- as.matrix(gse7307.src[,as.character(gse7307.s.sort[,1])])
tt <- tt[-4,] # removed "1556499_s_at"
tt <- as.data.frame(tt);
gse7307.norm <- as.character(subset(gse7307.s.sort, characteristics_ch1a=="endometrium")$Sample)
gse7307.tum <- as.character(subset(gse7307.s.sort, characteristics_ch1a=="endometrium/ovary")$Sample)
tt = cbind(tt[,gse7307.norm],tt[,gse7307.tum]) 
tt <- log2(tt);
pv.tt <- apply(tt, 1, func.list$studentT, s1=c(gse7307.norm),s2=c(gse7307.tum))
tt$p.value <- as.numeric(pv.tt)
tt.n <- apply(tt[,gse7307.norm], 1, mean, na.rm=T); tt$mean.N <- tt.n
tt.t <- apply(tt[,gse7307.tum], 1, mean, na.rm=T); tt$mean.T <- tt.t
tt$FC <- tt$mean.T - tt$mean.N
func.list$volcano(tt, fold.change.col = "FC", pv.col = "p.value", title.plot= "Volcano Plot, Endometrium/Ovary vs Normal (GSE7307)", cut.line = -log10(0.05), foldcut = 0.5)
gse7307.hv<-heatmap.plus(
	#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
	#as.matrix(temp[,(as.character(cl.sort[,1]))]),
	#as.matrix(gse7307.src[,as.character(gse7307.s.sort[,1])]), 
	as.matrix(tt[,1:(ncol(tt)-4)]),
	#temp.cc,
	na.rm=TRUE,
	scale="none",
	#RowSideColor=probe.cc,
	ColSideColors=build.top.map(tt, gse7307.s, color.map.gse7307.type),
	RowSideColors=build.side.map(tt),
	col=jet.colors(75),
	#key=FALSE,
	symkey=FALSE,
	density.info="none",
	trace="none",
	Rowv=TRUE,
	Colv=NA,
	cexRow=1,
	cexCol=0.6,
	keysize=1,
	dendrogram=c("none"),
	main = paste("gse7307; SRC only (73 probes) ",
		dim(tt[,1:(ncol(tt)-4)])[1],
		"Probes; ",
		 dim(tt[,1:(ncol(tt)-4)])[2],
		"Samples"),
	#labCol=NA
	labRow=NA
)
dev.off()

png(filename = "gse6364.src.png", bg="white", res=300, width=3000, height=3000)
ttt <- as.matrix(gse6364.src[,as.character(gse6364.s.sort[,1])])
ttt <- ttt[-4,] # outlier
ttt <- as.data.frame(ttt);
#ttt <- log2(ttt);
#### This dataset has a variety of tumors and normals. We need to select only the relevant ones
#ttt.n <- apply(ttt[,as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])], 1, mean, na.rm=T); ttt$mean.N <- ttt.n
#ttt.t <- apply(ttt[,as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"])], 1, mean, na.rm=T); ttt$mean.T <- ttt.t
#pv.ttt <- apply(ttt, 1, func.list$studentT, s1=c(as.character(gse6364.s.sort[c(7:9,19:26,33:37),"geo_accession"])),s2=c(as.character(gse6364.s.sort[c(1:6,10:18,27:32),"geo_accession"])))
#ttt$p.value <- as.numeric(pv.ttt)
#ttt$FC <- ttt$mean.T - ttt$mean.N
#func.list$volcano(ttt, fold.change.col = "FC", pv.col = "p.value", title.plot= "Volcano Plot, Endometrium/Ovary vs Normal (GSE6364)", cut.line = -log10(0.05), foldcut = 0.5)
gse6364.hv<-heatmap.plus(
		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
		#as.matrix(gse6364.src[,gse6364.s.sort[,1]]),
		as.matrix(ttt), 
		#temp.cc,
		na.rm=TRUE,
		scale="none",
		#RowSideColor=probe.cc,
		ColSideColors=cc.col.gse6364,
		col=jet.colors(75),
		#key=FALSE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=TRUE,
		Colv=NA,
		cexRow=1,
		cexCol=0.6,
		keysize=1,
		dendrogram=c("none"),
		main = paste("gse6364; SRC only (73 probes) ",dim(gse6364.src)[1],"Probes; ",dim(gse6364.src)[2],"Samples"),
		#labCol=NA
		labRow=NA
)
dev.off()
png(filename = "differentiated_genes.png", bg="white", res=300, width=3000, height=3000)
diff.genes.hv<-heatmap.plus(
		#as.matrix(temp[hv1[[1]],(as.character(cl.sort[,1]))]),
		#as.matrix(temp[,(as.character(cl.sort[,1]))]),
		#as.matrix(gse6364.src[,gse6364.s.sort[,1]]),
		as.matrix(subset(tt, p.value<=0.05 & FC <= -0.5)[,c(gse7307.norm,gse7307.tum)]), 
		#temp.cc,
		na.rm=TRUE,
		scale="none",
		#RowSideColor=probe.cc,
		#ColSideColors=cc.col.gse7307,
		col=jet.colors(75),
		#key=FALSE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		Rowv=TRUE,
		Colv=NA,
		cexRow=1,
		cexCol=0.6,
		keysize=1,
		dendrogram=c("none"),
		main = paste("gse7307; SRC only (73 probes) ",dim(gse6364.src)[1],"Probes; ",dim(gse6364.src)[2],"Samples"),
		#labCol=NA
		labRow=NA
)
dev.off()
