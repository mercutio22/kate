#func.list files:
func.list<-list()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### plot histogram
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$plot.hist<-function(df,dist.meth,hc.meth){

        plot(
                hclust(
                        dist(
                                t(df),method=dist.meth
                                )
                        ,method=hc.meth),
                labels=names(t(df)),cex=0.75,hang= -3
        )

}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make volcano plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$volcano<-function(x,fold.change.col,pv.col, title.plot, cut.line, foldcut) {
	xlim.r<-(range(x[,fold.change.col]))
	if(abs(xlim.r[1])>abs(xlim.r[2])){
#png(file="volcano_plot.png", bg="white", res=300, width=3000, height=3000)
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(xlim.r[1], abs(xlim.r[1])), #x-axis limits
				main=title.plot,
				xlab="Gene Expression\nlog2(fold change)",
				ylab="-1 * log10 of the Significance",
				type = "n"
		
		)
		foldcut.neg <- -1*foldcut
		abline(v=foldcut, lty=1, lwd=3)
		abline(v=foldcut.neg, lty=1, lwd=3)
		abline(h=cut.line, lty=2, lwd=2)
		
		p <- x[,pv.col];
		M <- x[,fold.change.col];
		points(M, -log10(p), col="black", pch=20, cex = 0.5) 
		upr = -1*log10(p) >= cut.line & M >= foldcut
		points(M[upr], -log10(p)[upr], col="red", pch=20)


		numr = cbind(M[upr],-log10(p)[upr])

		numr = dim(numr)[1]

		upl = -1*log10(p) >= cut.line & M <= foldcut.neg
		points(M[upl], -log10(p)[upl], col="darkgreen", pch=20)

		numl = cbind(M[upl],-log10(p)[upl])

		numl = dim(numl)[1]

		#all.s = M[M!=upl] & M[M!=upr] & -1*log10(p) < cut.line

		#points(M[all.s], -log10(p)[all.s], col="black", pch=20, cex = 0.5)  

     		op <- par(bg="white")
		legend("topright", paste("UP Reg = ",numr, sep=""), pch=20, col="red", inset = .02, cex = 0.5)
		legend("topleft", paste("DN Reg = ",numl, sep=""), pch=20, col="darkgreen", inset = .02, cex = 0.5)
     		par(op)
#dev.off()

	}else{
#png(file="volcano_plot.png", bg="white", res=300, width=3000, height=3000)
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(-1*xlim.r[2], xlim.r[2]), #x-axis limits
				main=title.plot,
				xlab="Gene Expression\nlog2(fold change)",
				ylab="-1 * log10 of the Significance",
				type = "n"
		
		)
		foldcut.neg <- -1*foldcut
		abline(v=foldcut, lty=1, lwd=3)
		abline(v=foldcut.neg, lty=1, lwd=3)
		abline(h=cut.line, lty=2, lwd=2)
		
		p <- x[,pv.col];
		M <- x[,fold.change.col];
		points(M, -log10(p), col="black", pch=20, cex = 0.5) 
		upr = -1*log10(p) >= cut.line & M >= foldcut
		points(M[upr], -log10(p)[upr], col="red", pch=20)

		numr = cbind(M[upr],-log10(p)[upr])

		numr = dim(numr)[1]

		upl = -1*log10(p) >= cut.line & M <= foldcut.neg
		points(M[upl], -log10(p)[upl], col="darkgreen", pch=20)

		numl = cbind(M[upl],-log10(p)[upl])

		numl = dim(numl)[1]

		#all.s = M[M!=upl] & M[M!=upr] & -1*log10(p) < cut.line

		#points(M[all.s], -log10(p)[all.s], col="black", pch=20, cex = 0.5)  

     		op <- par(bg="white")
		legend("topright", paste("UP Reg = ",numr, sep=""), pch=20, col="red", inset = .02, cex = 0.5)
		legend("topleft", paste("DN Reg = ",numl, sep=""), pch=20, col="darkgreen", inset = .02, cex = 0.5)
     		par(op)
#dev.off()

	}
}
func.list$volcano.dnamethylation<-function(x,fold.change.col,pv.col, title.plot, cut.line, foldcut) {
	xlim.r<-(range(x[,fold.change.col]))
	if(abs(xlim.r[1])>abs(xlim.r[2])){
#png(file="volcano_plot.png", bg="white", res=300, width=3000, height=3000)
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(xlim.r[1], abs(xlim.r[1])), #x-axis limits
				main=title.plot,
				xlab="DNA Methylation\nBeta-value difference",
				ylab="-1 * log10 of the Significance",
				type = "n"
		
		)
		foldcut.neg <- -1*foldcut
		abline(v=foldcut, lty=1, lwd=3)
		abline(v=foldcut.neg, lty=1, lwd=3)
		abline(h=cut.line, lty=3, lwd=3)
		
		p <- x[,pv.col];
		M <- x[,fold.change.col];
		points(M, -log10(p), col="black", pch=20, cex = 0.5) 
		upr = -1*log10(p) >= cut.line & M >= foldcut
		points(M[upr], -log10(p)[upr], col="red", pch=20)


		numr = cbind(M[upr],-log10(p)[upr])

		numr = dim(numr)[1]

		upl = -1*log10(p) >= cut.line & M <= foldcut.neg
		points(M[upl], -log10(p)[upl], col="darkblue", pch=20)

		numl = cbind(M[upl],-log10(p)[upl])

		numl = dim(numl)[1]

		#all.s = M[M!=upl] & M[M!=upr] & -1*log10(p) < cut.line

		#points(M[all.s], -log10(p)[all.s], col="black", pch=20, cex = 0.5)  

     		op <- par(bg="white")
		legend("topright", paste("Hyper = ",numr, sep=""), pch=20, col="red", inset = .02, cex = 0.5)
		legend("topleft", paste("Hypo = ",numl, sep=""), pch=20, col="darkblue", inset = .02, cex = 0.5)
     		par(op)
#dev.off()

	}else{
#png(file="volcano_plot.png", bg="white", res=300, width=3000, height=3000)
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(-1*xlim.r[2], xlim.r[2]), #x-axis limits
				main=title.plot,
				xlab="DNA Methylation\nBeta-value difference",
				ylab="-1 * log10 of the Significance",
				type = "n"
		
		)
		foldcut.neg <- -1*foldcut
		abline(v=foldcut, lty=1, lwd=3)
		abline(v=foldcut.neg, lty=1, lwd=3)
		abline(h=cut.line, lty=3, lwd=3)
		
		p <- x[,pv.col];
		M <- x[,fold.change.col];
		points(M, -log10(p), col="black", pch=20, cex = 0.5) 
		upr = -1*log10(p) >= cut.line & M >= foldcut
		points(M[upr], -log10(p)[upr], col="red", pch=20)

		numr = cbind(M[upr],-log10(p)[upr])

		numr = dim(numr)[1]

		upl = -1*log10(p) >= cut.line & M <= foldcut.neg
		points(M[upl], -log10(p)[upl], col="darkblue", pch=20)

		numl = cbind(M[upl],-log10(p)[upl])

		numl = dim(numl)[1]

		#all.s = M[M!=upl] & M[M!=upr] & -1*log10(p) < cut.line

		#points(M[all.s], -log10(p)[all.s], col="black", pch=20, cex = 0.5)  

     		op <- par(bg="white")
		legend("topright", paste("Hyper = ",numr, sep=""), pch=20, col="red", inset = .02, cex = 0.5)
		legend("topleft", paste("Hypo = ",numl, sep=""), pch=20, col="darkblue", inset = .02, cex = 0.5)
     		par(op)
#dev.off()

	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### student t.test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$studentT.VAREQUAL<-function(x,s1,s2) {
    x1 <- x[s1]
    x2 <- x[s2]
    x1 <- as.numeric(x1)
    x2 <- as.numeric(x2)
    t.out <- t.test(x1,x2, alternative="two.sided", na.rm=TRUE, var.equal=TRUE)
    out <- as.numeric(t.out$p.value)
    return(out)
}
func.list$studentT<-function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- t.test(x1,x2, alternative="two.sided", na.rm=TRUE)
	out <- as.numeric(t.out$p.value)
	return(out)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fisher exact test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$fisher.test<-function(x,tot.exp, tot.obs,sub.obs, sub.exp) {
	x1 <- x[sub.obs]
	x2 <- x[sub.exp]
	#x1 <- as.numeric(x1)
	#x2 <- as.numeric(x2)
	t.out <- fisher.test(matrix(c(tot.exp-x2,tot.obs-x1,x2,x1),nrow=2,ncol=2))
	out <- as.numeric(t.out$p.value)
	return(out)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### vlookup method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$vlookup<-function(val, df, col){

               df[df[1] == val, col][1]
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### histogram of a specific gene
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$hist<-function(data,gene){

par(mfcol=c(2,1))
hist(as.numeric(data[gene,gcimp]), xlim=c(0,1), main="Histogram - gCIMP")
hist(as.numeric(data[gene,ngcimps]), xlim=c(0,1), main="Histogram - NONgCIMP")

}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### wilcox exact test (deals with ties)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$wilcox.exact<-function(x,s1,s2) {
    x1 <- x[s1]
    x2 <- x[s2]
    x1 <- as.numeric(x1)
    x2 <- as.numeric(x2)
    t.out <- wilcox.exact(x1,x2, alternative="two.sided", exact=TRUE, na.rm=TRUE)
    out <- as.numeric(t.out$p.value)
    return(out)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### plot correlation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$plot.cor<-function(df){

dat.cor<-cor(df)
image(dat.cor,axes=F)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])
axis(3,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])

}

func.list$plot.cor.hclust<-function(df, dist.method, hclust.method, title.main, col.lab, path){

#t.df <-t(df) #transpose dat
#dat.dist<-dist(t.df,method=dist.method)
#dat.clust<-hclust(dat.dist,method=hclust.method)



dat.cor<-cor(df, use=hclust.method, method=dist.method)
#image(dat.cor,axes=F,col=jet.colors(75))
#axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])
#axis(3,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])
dat.cor1<-format(round(dat.cor, 2))
#x11()
heatmap.2(
		as.matrix(dat.cor),
		na.rm=TRUE,
		scale="none",
		col=jet.colors(75),
		key=TRUE,
		symkey=FALSE,
		density.info="none",
		trace="none",
		cexRow=1,
		cexCol=1,
		dendrogram="both",
		Colv=TRUE,
		Rowv=TRUE,
		#cellnote=dat.cor1,
		#notecex=1,
		#notecol="white",
		keysize=1,
		ColSideColors=col.lab,
		RowSideColors=col.lab,
		xlab = "SAMPLES",
		ylab = "SAMPLES",
		main = paste(title.main, "\n", hclust.method," : ",dist.method)
	)
	png(filename = path, bg="white", res=300, width=3000, height=3000)
	heatmap.2(
			as.matrix(dat.cor),
			na.rm=TRUE,
			scale="none",
			col=jet.colors(75),
			key=TRUE,
			symkey=FALSE,
			density.info="none",
			trace="none",
			cexRow=0.5,
			cexCol=0.5,
			dendrogram="both",
			Colv=TRUE,
			Rowv=TRUE,
			#cellnote=dat.cor1,
			#notecex=1,
			#notecol="white",
			keysize=1,
			ColSideColors=col.lab,
			RowSideColors=col.lab,
			xlab = "SAMPLES",
			ylab = "SAMPLES",
			main = title.main
	)
	dev.off()
#x11()
#plot(
#		dat.clust,
#		labels=names(t(df)), cex=0.75, hang= -3, main = title.main
#	)
}
func.list$consensus.clustering<-function(df, title.main, col.lab, path){
	heatmap.2(
			as.matrix(df),
			na.rm=TRUE,
			scale="none",
			col=jet.colors(75),
			key=TRUE,
			symkey=FALSE,
			density.info="none",
			trace="none",
			cexRow=1,
			cexCol=1,
			dendrogram="none",
			Colv=FALSE,
			Rowv=FALSE,
			#cellnote=dat.cor1,
			#notecex=1,
			#notecol="white",
			keysize=1,
			ColSideColors=col.lab,
			RowSideColors=col.lab,
			xlab = "SAMPLES",
			ylab = "SAMPLES",
			main = title.main
	)
	png(filename = path, bg="white", res=300, width=3000, height=3000)
	heatmap.2(
			as.matrix(df),
			na.rm=TRUE,
			scale="none",
			col=jet.colors(75),

			key=TRUE,
			symkey=FALSE,
			density.info="none",
			trace="none",
			cexRow=1,
			cexCol=1,
			dendrogram="none",
			Colv=FALSE,
			Rowv=FALSE,
			#cellnote=dat.cor1,
			#notecex=1,
			#notecol="white",
			keysize=1,
			ColSideColors=col.lab,
			RowSideColors=col.lab,
			xlab = "SAMPLES",
			ylab = "SAMPLES",
			main = title.main
	)
	dev.off()
}
func.list$bb.test<-function(x,obs,totobs,obsref,totobsref) {
	x1 <- x[obs]
	x3 <- x[obsref]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else{
    x1 <- as.numeric(x1)
    x3 <- as.numeric(x3)
    #need to subtract the obs value from the totobs
    t.out <- binom.test(c(x1,as.numeric(totobs)-x1), p=x3/as.numeric(totobsref), alternative="two.sided")
    out <- as.numeric(t.out$p.value)
    return(out)
	}
}
odd.ratio<-function(x,ob,tot.ob,ex,tot.ex){
	x1 <- x[ob]
	x3 <- x[ex]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x3 <- as.numeric(x3)
		top <- (x1/tot.ob)/(1-(x1/tot.ob))
		bot <- (x3/tot.ex)/(1-(x3/tot.ex))
		out <- as.numeric(top/bot)
		return(out)
	}
}
log.odd.ratio<-function(x,ob,tot.ob,ex,tot.ex){
	x1 <- x[ob]
	x3 <- x[ex]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x3 <- as.numeric(x3)
		top <- (x1/tot.ob)/(1-(x1/tot.ob))
		bot <- (x3/tot.ex)/(1-(x3/tot.ex))
		out <- log10(as.numeric(top/bot))
		return(out)
	}
}
odd.ratio.mut<-function(x,m.pos,m.neg,w.pos,w.neg){
	x1 <- x[m.pos]
	x2 <- x[m.neg]
	x3 <- x[w.pos]
	x4 <- x[w.neg]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x2)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else if(is.na(x4)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x2 <- as.numeric(x2)
		x3 <- as.numeric(x3)
		x4 <- as.numeric(x4)
		top <- (x1/x2)
		bot <- (x3/x4)
		out <- as.numeric(top/bot)
		return(out)
	}
}
log.odd.ratio.mut<-function(x,m.pos,m.neg,w.pos,w.neg){
	x1 <- x[m.pos]
	x2 <- x[m.neg]
	x3 <- x[w.pos]
	x4 <- x[w.neg]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x2)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else if(is.na(x4)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x2 <- as.numeric(x2)
		x3 <- as.numeric(x3)
		x4 <- as.numeric(x4)
		top <- (x1/x2)
		bot <- (x3/x4)
		out <- log10(as.numeric(top/bot))
		return(out)
	}
}
fisher.test.mut<-function(x){
	x1 <- x[1]
	x2 <- x[2]
	x3 <- x[4]
	x4 <- x[5]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x2)){
		return("NA")
	}else if(is.na(x3)){
		return("NA")
	}else if(is.na(x4)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x2 <- as.numeric(x2)
		x3 <- as.numeric(x3)
		x4 <- as.numeric(x4)
		out <- fisher.test(matrix(c(x1,x2,x3,x4),ncol=2,nrow=2,byrow=T))
		out <- out$p.value
		return(out)
	}
}
fisher.test.mut1<-function(x, t, e){
	x1 <- x[3]
	x2 <- x[9]
	if(is.na(x1)){
		return("NA")
	}else if(is.na(x2)){
		return("NA")
	}else{
		x1 <- as.numeric(x1)
		x2 <- as.numeric(x2)
		x3 <- as.numeric(t)
		x4 <- as.numeric(e)
		out <- fisher.test(matrix(c(x1,x2,x3-x1,x4-x2),ncol=2,nrow=2,byrow=T))
		out <- out$p.value
		return(out)
	}
}
func.list$countif<-function(x) { #x = data frame; cols = the columns you want to evaluate; what = the character you want to find, e.g. "NA"
	x1 <- apply(is.na(x),2,sum) #extract values based on cols
	x1.n<- dim(x)[1]
	per <- (x1/(x1.n))*100 #makes percent of the what per number of cols.  use this value to extract out number of NAs
	per <- round(per,2) #round the per so 2 digits decimals
	return(per)
}
func.list$countif.na<-function(x,cols,what) { #x = data frame; cols = the columns you want to evaluate; what = the character you want to find, e.g. "NA"
	x1 <- x[cols] #extract values based on cols
	x1 <- as.numeric(x1)
	out <- length(x1[x1[]==what]) #sums all 'what' within the evaluated cols.
	out <- as.numeric(out)
	per <- (out/length(cols))*100 #makes percent of the what per number of cols.  use this value to extract out number of NAs
	per <- round(per,0) #round the per so no decimals
	return(per)
}

#to sort data.frame matrix.  Use this function.  then use sort.data.frame(dat, key = "LOC")
sort.data.frame <- function(x, key, ...) {
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}


#read from clipboard for Mac (thanks to Ken Knoblauch for this hint
#for PCs, use read.table(file("clipboard"))
# cities <- read.clipboard(header="TRUE")
read.clipboard<-function(header=TRUE,...) {
	MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
	if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,...)
		else read.clipboard<-read.table(file("clipboard"),...) }
	else {
		if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,...)
		else read.clipboard<- read.table(pipe("pbpaste") ,...)}
}

read.clipboard.csv<-function(header=TRUE,sep=',',...) {  #same as read.clipboard(sep=',')
	MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
	if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,...)
		else read.clipboard<-read.table(file("clipboard"),...) }
	else {
		if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,...)
		else read.clipboard<- read.table(pipe("pbpaste") ,...)}
}

#### correlations matrix between groups of samples.  outputs plot with marked significant correlations
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- (cor(x, y, use = "complete.obs", method="spearman"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

	test <- cor.test(x,y, alternative="two.sided")
	# borrowed from printCoefmat
	Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
			cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
			symbols = c("***", "**", "*", ".", " "))

	text(0.5, 0.5, txt, cex = cex * r)
	text(.8, .8, Signif, cex=cex, col=2)
}
#example how to plot 'panel.cor'
#pairs(USJudgeRatings[,c(2:3,6,1,7)], lower.panel=panel.smooth, upper.panel=panel.cor)

plotmatrix.hh<-function (data, mapping = aes(), colour = "black")
{
	grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
	grid <- subset(grid, x != y)
	all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
						xcol <- grid[i, "x"]
						ycol <- grid[i, "y"]
						data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol],
								x = data[, xcol], y = data[, ycol], data)
					}))
	all$xvar <- factor(all$xvar, levels = names(data))
	all$yvar <- factor(all$yvar, levels = names(data))
	densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
						data.frame(xvar = names(data)[i], yvar = names(data)[i],
								x = data[, i])
					}))
	mapping <- defaults(mapping, aes_string(x = "x", y = "y"))
	class(mapping) <- "uneval"
	ggplot(all, mapping) + facet_grid(xvar ~ yvar) + #removed 'scales="free"'
			geom_point(colour = colour) + stat_density(aes(x = x,
							y = ..scaled.. * diff(range(x)) + min(x)), data = densities,
					position = "identity", colour = "grey20", geom = "line")
}
theme_white <- function(colss) {
	theme_update (

			plot.background = theme_blank(),
			panel.background=theme_rect(colour=colss, size=2),
			axis.text.x= theme_text(colour="black", angle=90, hjust=1),
			axis.text.y= theme_text(colour="black", hjust=1),
			axis.title.x =theme_text(colour="black",face="bold"),
			axis.title.y =theme_text(colour="black",face="bold", angle = 90)
	#panel.grid.major = theme_bw()
	)
}


####

#use tabletrend function

print.ordtest=function(l,...)
{
	tmp=matrix(c(l$estimate,l$ASE),nrow=1)
	dimnames(tmp)=list(l$name,c("Estimate","ASE"))
	print(round(tmp,4),...)
}


compADPQ=function(x)
{
	nr=nrow(x)
	nc=ncol(x)
	Aij=matrix(0,nrow=nr,ncol=nc)
	Dij=matrix(0,nrow=nr,ncol=nc)
	for (i in 1:nr)	{
		for (j in 1:nc)	{
			
			Aij[i,j]=sum(x[1:i,1:j])+sum(x[i:nr,j:nc])-sum(x[i,])-sum(x[,j])
			
			Dij[i,j]=sum(x[i:nr,1:j])+sum(x[1:i,j:nc])-sum(x[i,])-sum(x[,j])
		}}
	P=sum(x*Aij)
	Q=sum(x*Dij)
	return(list(Aij=Aij,Dij=Dij,P=P,Q=Q))
}


scores=function(x,MARGIN=1,method="table",...)
{
	# MARGIN
	#	1 - row
	# 	2 - columns
	
	# Methods for ranks are
	#
	# x - default
	# rank
	# ridit
	# modridit
	
	if (method=="table")
	{
		if (is.null(dimnames(x))) return(1:(dim(x)[MARGIN]))
		else {
			options(warn=-1)
			if
			(sum(is.na(as.numeric(dimnames(x)[[MARGIN]])))>0)
			{
				out=(1:(dim(x)[MARGIN]))
			}
			else
			{
				out=(as.numeric(dimnames(x)[[MARGIN]]))
			}
			options(warn=0)
		}
	}
	else	{
		### method is a rank one
		Ndim=dim(x)[MARGIN]
		OTHERMARGIN=3-MARGIN
		
		ranks=c(0,(cumsum(apply(x,MARGIN,sum))))[1:Ndim]+(apply(x,MARGIN,sum)+1)/2
		if (method=="ranks") out=ranks
		if (method=="ridit") out=ranks/(sum(x))
		if (method=="modridit") out=ranks/(sum(x)+1)
	}
	
	return(out)
}


#------------------------------------ FUNCTIONS------------------------------------#
		
		tablegamma=function(x)
{
# Statistic
	tmp=compADPQ(x)
	P=tmp$P
	Q=tmp$Q
	gamma=(P-Q)/(P+Q)
# ASE
	Aij=tmp$Aij
	Dij=tmp$Dij
	tmp1=4/(P+Q)^2
	tmp2=sqrt(sum((Q*Aij - P*Dij)^2 * x))
	gamma.ASE=tmp1*tmp2
# Test
	var0=(4/(P+Q)^2) * (sum(x*(Aij-Dij)^2) - ((P-Q)^2)/sum(x))
	tb=gamma/sqrt(var0)
	p.value=2*(1-pnorm(tb))
# Output
	
	out=list(estimate=gamma,ASE=gamma.ASE,statistic=tb,p.value=p.value,name=
					"Gamma",bornes=c(-1,1))
	class(out)="ordtest"
	return(out)
}


tabletauc=function(x)
{
	tmp=compADPQ(x)
	P=tmp$P
	Q=tmp$Q
	m=min(dim(x))
	n=sum(x)
# statistic
	
	tauc=(m*(P-Q))/(n^2*(m-1))
# ASE
	Aij=tmp$Aij
	Dij=tmp$Dij
	dij=Aij-Dij
	tmp1=2*m/((m-1)*n^2)
	tmp2= sum(x * dij^2) - (P-Q)^2/n
	ASE=tmp1*sqrt(tmp2)
	
# Test
	tb=tauc/ASE
	p.value=2*(1-pnorm(tb))
# Output
	
	out=list(estimate=tauc,ASE=ASE,statistic=tb,p.value=p.value,name="Kendal
					l's tau-c",bornes=c(-1,1))
	class(out)="ordtest"
	return(out)
}

tabletaub=function(x)
{
# Statistic
	tmp=compADPQ(x)
	P=tmp$P
	Q=tmp$Q
	n=sum(x)
	wr=n^2 - sum(apply(x,1,sum)^2)
	wc=n^2 - sum(apply(x,2,sum)^2)
	taub=(P-Q)/sqrt(wr*wc)
# ASE
	Aij=tmp$Aij
	Dij=tmp$Dij
	w=sqrt(wr*wc)
	dij=Aij-Dij
	nidot=apply(x,1,sum)
	ndotj=apply(x,2,sum)
	n=sum(x)
	vij=outer(nidot,ndotj, FUN=function(a,b) return(a*wc+b*wr))
	tmp1=1/(w^2)
	tmp2= sum(x*(2*w*dij + taub*vij)^2)
	tmp3=n^3*taub^2*(wr+wc)^2
	tmp4=sqrt(tmp2-tmp3)
	taub.ASE=tmp1*tmp4
# Test
	var0=4/(wr*wc) * (sum(x*(Aij-Dij)^2) - (P-Q)^2/n)
	tb=taub/sqrt(var0)
	p.value=2*(1-pnorm(tb))
# Output
	
	out=list(estimate=taub,ASE=taub.ASE,statistic=tb,p.value=p.value,name="K
					endall's tau-b",bornes=c(-1,1))
	class(out="ordtest")
	return(out)
}

tablesomersD=function(x,dep=2)
{
	# dep: which dimension stands for the dependant variable
	# 1 - ROWS
	# 2 - COLS
# Statistic
	if (dep==1) x=t(x)
	tmp=compADPQ(x)
	P=tmp$P
	Q=tmp$Q
	n=sum(x)
	wr=n^2 - sum(apply(x,1,sum)^2)
	somers=(P-Q)/wr
# ASE
	Aij=tmp$Aij
	Dij=tmp$Dij
	dij=Aij-Dij
	tmp1=2/wr^2
	tmp2=sum(x*(wr*dij - (P-Q)*(n-apply(x,1,sum)))^2)
	ASE=tmp1*sqrt(tmp2)
# Test
	var0=4/(wr^2) * (sum(x*(Aij-Dij)^2) - (P-Q)^2/n)
	tb=somers/sqrt(var0)
	p.value=2*(1-pnorm(tb))
# Output
	if (dep==1) dir="R|C" else dir= "C|R"
	name=paste("Somer's D",dir)
	
	out=list(estimate=somers,ASE=ASE,statistic=tb,p.value=p.value,name=name,
			bornes=c(-1,1))
	class(out)="ordtest"
	return(out)
	
	
}

#out=table.somersD(data)



tablepearson=function(x,scores.type="table")
{
	
# Statistic
	sR=scores(x,1,scores.type)
	sC=scores(x,2,scores.type)
	n=sum(x)
	Rbar=sum(apply(x,1,sum)*sR)/n
	Cbar=sum(apply(x,2,sum)*sC)/n
	ssr=sum(x*(sR-Rbar)^2)
	ssc=sum(t(x)* (sC-Cbar)^2)
	tmpij=outer(sR,sC,FUN=function(a,b) return((a-Rbar)*(b-Cbar)))
	ssrc= sum(x*tmpij)
	v=ssrc
	w=sqrt(ssr*ssc)
	r=v/w
# ASE
	bij=outer(sR,sC, FUN=function(a,b)return((a-Rbar)^2*ssc +
								(b-Cbar)^2*ssr))
	tmp1=1/w^2
	tmp2=x*(w*tmpij - (bij*v)/(2*w))^2
	tmp3=sum(tmp2)
	ASE=tmp1*sqrt(tmp3)
# Test
	var0= (sum(x*tmpij) - (ssrc^2/n))/ (ssr*ssc)
	tb=r/sqrt(var0)
	p.value=2*(1-pnorm(tb))
# Output
	
	out=list(estimate=r,ASE=ASE,statistic=tb,p.value=p.value,name="Pearson
					Correlation",bornes=c(-1,1))
	class(out)="ordtest"
	return(out)
}

# table.pearson(data)


tablespearman=function(x)
{
	# Details algorithme manuel SAS PROC FREQ page 540
# Statistic
	n=sum(x)
	nr=nrow(x)
	nc=ncol(x)
	tmpd=cbind(expand.grid(1:nr,1:nc))
	ind=rep(1:(nr*nc),as.vector(x))
	tmp=tmpd[ind,]
	rhos=cor(apply(tmp,2,rank))[1,2]
# ASE
	Ri=scores(x,1,"ranks")- n/2
	Ci=scores(x,2,"ranks")- n/2
	sr=apply(x,1,sum)
	sc=apply(x,2,sum)
	F=n^3 - sum(sr^3)
	G=n^3 - sum(sc^3)
	w=(1/12)*sqrt(F*G)
	vij=data
	for (i in 1:nrow(x))
	{
		qi=0
		if (i<nrow(x))
		{
			for (k in i:nrow(x)) qi=qi+sum(x[k,]*Ci)
		}
	}
	for (j in 1:ncol(x))
	{
		qj=0
		if (j<ncol(x))
		{
			for (k in j:ncol(x)) qj=qj+sum(x[,k]*Ri)
		}
		vij[i,j]=n*(Ri[i]*Ci[j] +
					0.5*sum(x[i,]*Ci)+0.5*sum(data[,j]*Ri) +qi+qj)
	}
	
	
	v=sum(data*outer(Ri,Ci))
	wij=-n/(96*w)*outer(sr,sc,FUN=function(a,b) return(a^2*G+b^2*F))
	zij=w*vij-v*wij
	zbar=sum(data*zij)/n
	vard=(1/(n^2*w^4))*sum(x*(zij-zbar)^2)
	ASE=sqrt(vard)
# Test
	vbar=sum(x*vij)/n
	p1=sum(x*(vij-vbar)^2)
	p2=n^2*w^2
	var0=p1/p2
	stat=rhos/sqrt(var0)
	
# Output
	out=list(estimate=rhos,ASE=ASE,name="Spearman
					correlation",bornes=c(-1,1))
	class(out)="ordtest"
	return(out)
}

#tablespearman(data)




tablelambdasym=function(x)
{
# Statistic
	ri = apply(x,1,max)
	r=max(apply(x,2,sum))
	n=sum(x)
	cj=apply(x,2,max)
	c=max(apply(x,1,sum))
	sri=sum(ri)
	w=2*n - r -c
	v=2*n - sri - sum(cj)
	lambda=(w-v)/w
# ASE ...
	
	tmpSi=0
	l=min(which(apply(x,2,sum)==r))
	for (i in 1:length(ri))
	{
		li=min(which(x[i,]==ri[i]))
		if (li==l) tmpSi=tmpSi+x[i,li]
	}
	
	tmpSj=0
	k=min(which(apply(x,1,sum)==c))
	for (j in 1:length(cj))
	{
		kj=min(which(x[,j]==cj[j]))
		if (kj==k) tmpSj=tmpSj+x[kj,j]
	}
	
	rk=max(x[k,])
	cl=max(x[,l])
	tmpx=tmpSi+tmpSj+rk+cl
	y=8*n-w-v-2*tmpx
	
	nkl=x[k,l]
	tmpSij=0
	for (i in 1:nrow(x))
	{
		for (j in 1:ncol(x))
		{
			li=min(which(x[i,]==ri[i]))
			kj=min(which(x[,j]==cj[j]))
			tmpSij=tmpSij+x[kj,li]
		}
	}
	
	
	ASE=(1/(w^2))*sqrt(w*v*y-(2*(w^2)*(n-tmpSij))-2*(v^2)*(n-nkl))
# Output
	
	out=list(estimate=lambda,ASE=ASE,name="Lambda
					Symetric",bornes=c(0,1))
	class(out)="ordtest"
	return(out)
	
}
#tablelambdasym(data)

tablelambdaasym=function(x,transpose=FALSE)
{
# Statistic
	if (transpose==TRUE) x=t(x)
	ri = apply(x,1,max)
	r=max(apply(x,2,sum))
	sri=sum(ri)
	n=sum(x)
	lambda=(sum(ri)-r)/(n-r)
# ASE
	l=min(which(apply(x,2,sum)==r))
	tmp=0
	for (i in 1:length(ri))
	{
		li=min(which(x[i,]==ri[i]))
		if (li==l) tmp=tmp+x[i,li]
	}
	ASE=sqrt(((n-sri)/(n-r)^3) *(sri+r-2*tmp))
# Output
	if (transpose) dir="R|C" else dir= "C|R"
	name=paste("Lambda asymetric",dir)
	out=list(estimate=lambda,ASE=ASE,name=name,bornes=c(0,1))
	class(out)="ordtest"
	return(out)
}




tableUCA=function(x,transpose=TRUE)
{
	if (transpose==TRUE) x=t(x)
	# Statistic
	n=sum(x)
	ni=apply(x,1,sum)
	nj=apply(x,2,sum)
	Hx=-sum((ni/n)*log(ni/n))
	Hy=-sum((nj/n)*log(nj/n))
	Hxy=-sum((x/n)*log(x/n))
	v=Hx+Hy- Hxy
	w=Hy
	U=v/w
	# ASE
	tmp1=1/((n)*(w^2))
	tmpij=0
	for (i in 1:nrow(x))
	{
		for (j in 1:ncol(x))
		{
			tmpij=tmpij+(  x[i,j]*  (
							Hy*log(x[i,j]/ni[i])+(Hx-Hxy)*    log(nj[j]/n)               )^2     )
		}
	}
	ASE=tmp1*sqrt(tmpij)
	# Output
	if (transpose) dir="R|C" else dir= "C|R"
	name=paste("Uncertainty Coefficient",dir)
	out=list(estimate=U,ASE=ASE,name=name,bornes=c(0,1))
	class(out)="ordtest"
	return(out)
}

#tableUCA(data)

tableUCS=function(x)
{
	
	# Statistic
	n=sum(x)
	ni=apply(x,1,sum)
	nj=apply(x,2,sum)
	Hx=-sum((ni/n)*log(ni/n))
	Hy=-sum((nj/n)*log(nj/n))
	Hxy=-sum((x/n)*log(x/n))
	U=(2*(Hx+Hy-Hxy))/(Hx+Hy)
	# ASE
	tmpij=0
	for (i in 1:nrow(x))
	{
		for (j in 1:ncol(x))
		{
			tmpij=tmpij+(  x[i,j]*
						(Hxy*log(ni[i]*nj[j]/n^2) - (Hx+Hy)*log(x[i,j]/n))^2   /(n^2*(Hx+Hy)^4)
						)
		}
	}
	ASE=2*sqrt(tmpij)
	# Output
	name="Uncertainty Coefficient Symetric"
	out=list(estimate=U,ASE=ASE,name=name,bornes=c(0,1))
	class(out)="ordtest"
	return(out)
}





tablelinear=function(x,scores.type="table")
{
	r=tablepearson(x,scores.type)$estimate
	n=sum(x)
	ll=r^2*(n-1)
	out=list(estimate=ll)
	return(out)
}

tablephi=function(x)
{
	if (all.equal(dim(x),c(2,2))==TRUE)
	{
		rtot=apply(x,1,sum)
		ctot=apply(x,2,sum)
		phi= det(x)/sqrt(prod(rtot)*prod(ctot))
	}
	else {
		Qp=chisq.test(x)$statistic
		phi=sqrt(Qp/sum(x))
	}
	names(phi)="phi"
	return(phi=phi)
}


tableCramerV=function(x)
{
	if (all.equal(dim(x),c(2,2))==TRUE)
	{
		cramerV=tablephi(x)
	}
	else
	{
		Qp=tableChisq(x)$estimate
		cramerV=sqrt((Qp/n)/min(dim(x)-1))
	}
	names(cramerV)="Cramer's V"
	return(cramerV)
}


tableChisq=function(x)
{
	nidot=apply(x,1,sum)
	ndotj=apply(x,2,sum)
	n=sum(nidot)
	eij=outer(nidot,ndotj,"*")/n
	R=length(nidot)
	C=length(ndotj)
	dll=(R-1)*(C-1)
	Qp=sum((x-eij)^2/eij)
	p.value=1-pchisq(Qp,dll)
	
	out=list(estimate=Qp,dll=dll,p.value=p.value,dim=c(R,C),name="Pearson's
					Chi-square")
	return(out)
}


tableChisqLR=function(x)
{
# Likelihood ratio Chi-squared test
	nidot=apply(x,1,sum)
	ndotj=apply(x,2,sum)
	n=sum(nidot)
	eij=outer(nidot,ndotj,"*")/n
	R=length(nidot)
	C=length(ndotj)
	dll=(R-1)*(C-1)
	G2=2*sum(x*log(x/eij))
	p.value=1-pchisq(G2,dll)
	
	out=list(estimate=G2,dll=dll,p.value=p.value,dim=c(R,C),name="Likelihood
					ratio Chi-square")
	return(out)
}

tableChisqCA=function(x)
{
	if (all.equal(dim(x),c(2,2))==TRUE)
	{
		nidot=apply(x,1,sum)
		ndotj=apply(x,2,sum)
		n=sum(nidot)
		eij=outer(nidot,ndotj,"*")/n
		R=length(nidot)
		C=length(ndotj)
		dll=(R-1)*(C-1)
		tmp=as.vector(abs(x-eij))
		tmp=pmax(tmp-0.5,0)
		tmp=matrix(tmp,byrow=TRUE,ncol=C)
		Qc=sum(tmp^2/eij)
		p.value=1-pchisq(Qc,dll)
		
		out=list(estimate=Qc,dll=dll,p.value=p.value,dim=c(R,C),name="Continuity
						adjusted Chi-square")
		return(out)
	}
	else
	{ stop("Continuity-adjusted chi-square must be used with
						(2,2)-tables",call.=FALSE) }
}

tableChisqMH=function(x)
{
	n=sum(x)
	G2=(n-1)*(tablepearson(x)$estimate^2)
	dll=1
	p.value=1-pchisq(G2,dll)
	
	out=list(estimate=G2,dll=dll,p.value=p.value,dim=dim(x),name="Mantel-Hae
					nszel Chi-square")
	return(out)
	
}

tableCC=function(x)
{
	Qp=tableChisq(x)$estimate
	n=sum(x)
	P=sqrt(Qp/(Qp+n))
	m=min(dim(x))
	
	out=list(estimate=P,dim=dim(x),bornes=c(0,sqrt((m-1)/m)),name="Contingen
					cy coefficient")
	return(out)
	
}

tabletrend=function(x,transpose=FALSE)
{
	if (any(dim(x)==2))
	{
		if (transpose==TRUE) {
			x=t(x)
		}
		
		if (dim(x)[2]!=2){stop("Cochran-Armitage test for trend must be
							used with a (R,2) table. Use transpose argument",call.=FALSE) }
		
		nidot=apply(x,1,sum)
		n=sum(nidot)
		
		Ri=scores(x,1,"table")
		Rbar=sum(nidot*Ri)/n
		
		s2=sum(nidot*(Ri-Rbar)^2)
		pdot1=sum(x[,1])/n
		T=sum(x[,1]*(Ri-Rbar))/sqrt(pdot1*(1-pdot1)*s2)
		p.value.uni=1-pnorm(abs(T))
		p.value.bi=2*p.value.uni
		
		out=list(estimate=T,dim=dim(x),p.value.uni=p.value.uni,p.value.bi=p.value.bi,name="Cochran-Armitage test for trend")
		return(out)
		
	}
	else {stop("Cochran-Armitage test for trend must be used with a
						(2,C) or a (R,2) table",call.=FALSE) }
}




