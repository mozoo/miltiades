##############################################################################################################################
#                                                                                                                            #
# Miltiades is a simple R script to perform sliding window-based analyses on protein composition.                            #
#                                                                                                                            #
# Copyright (C) 2016 Federico Plazzi                                                                                         #
#                                                                                                                            #
# This program is free software: you can redistribute it and/or modify                                                       #
# it under the terms of the GNU General Public License as published by                                                       #
# the Free Software Foundation, either version 3 of the License, or                                                          #
# (at your option) any later version.                                                                                        #
#                                                                                                                            #
# This program is distributed in the hope that it will be useful,                                                            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                              #
# GNU General Public License for more details.                                                                               #
#                                                                                                                            #
# You should have received a copy of the GNU General Public License                                                          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                                      #
#                                                                                                                            #
##############################################################################################################################

#Miltiades version: 1.0

#Loading seqinr package.
library(seqinr)

#Setting the default name for the FASTA alignment file and pointing to the first sequence.
infile <- "infile"
protein.position <- 1

#Checking for a user-specified alignment file.
for (i in 1:length(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE))) {
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###infile###",infile <- scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1],NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###sequence###",protein.position <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	}

#Reading the alignment file and defining its length (by default, the analysis is on the complete sequence).
alignment <- read.fasta(file=paste("../",infile,sep=""),seqtype="AA",forceDNAtolower=FALSE,strip.desc=TRUE)
protein.sequence <- alignment[protein.position][[1]]
protein.length <- length(protein.sequence)
stop <- protein.length

#Setting defaults.
outfile <- "outfile.pdf"
start <- 1
criteria <- c("charge","hydropathy","structure","chemical","standard area","accessible area","loop propensity")
criterion <- "hydropathy"
wsize <- 20
wstep <- 1
alpha <- 0.05
replicates <- 1000

#Reading user-defined parameters.
for (i in 1:length(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE))) {
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###outfile###",outfile <- scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1],NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###start###",start <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###stop###",stop <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###criterion###",criterion <- scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1],NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###wsize###",wsize <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###wstep###",wstep <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###alpha###",alpha <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	ifelse(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i] == "###replicates###",replicates <- as.numeric(scan("../miltiades.conf",what=character(),sep="\n",quiet=TRUE)[i+1]),NA)
	}

#Aminoacid groupings.
IUPAC.degeneracy <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","B","Z","X","?")
acidic <- c("D","E")
basic <- c("H","K","R")
NC <- c("A","C","F","G","I","L","M","N","P","Q","S","T","V","W","Y")
polar.uncharged <- c("C","G","N","Q","S","T","Y")
hydrophobic.nonpolar <- c("A","F","I","L","M","P","V","W")
external <- c("D","E","H","K","N","Q","R")
internal <- c("F","I","L","M","V")
ambivalent <- c("A","C","G","P","S","T","W","Y")
aliphatic <- c("I","L","V")
aliphatic.small <- c("A","G")
amide <- c("N","Q")
aromatic <- c("F","W","Y")
hydroxyl <- c("S","T")
imino <- "P"
sulfur <- c("C","M")
unknown <- c("X","?")

#Aminoacid areas and loop propensities.
standard.area <- matrix(c(88.1,"G",118.2,"A",129.8,"S",146.1,"C",146.8,"P",152.5,"T",158.7,"D",164.5,"V",165.5,"N",181.0,"I",186.2,"E",193.1,"L",193.2,"Q",202.5,"H",203.3,"M",222.8,"F",225.8,"K",238.8,"Y",256.0,"R",266.2,"W",162.1,"B",189.7,"Z",181.955,"X",181.955,"?"),nrow=24,ncol=2,byrow=TRUE)
accessible.area <- matrix(c(13.9,"C",23.0,"I",23.5,"V",25.2,"G",28.7,"F",29.0,"L",30.5,"M",31.5,"A",41.7,"W",44.2,"S",46.0,"T",46.7,"H",53.7,"P",59.1,"Y",60.9,"D",62.2,"N",72.3,"E",74.0,"Q",93.8,"R",110.3,"K",61.55,"B",73.15,"Z",48.51,"X",48.51,"?"),nrow=24,ncol=2,byrow=TRUE)
loop.propensities <- matrix(c(0.2299,"N",-0.2615,"A",-0.01515,"C",-0.2047,"E",-0.2434,"W",-0.2075,"Y",0.33793,"L",0.5523,"P",-0.1766,"R",0.0089,"T",-0.3862,"V",0.4332,"G",-0.4222,"I",-0.22590,"M",-0.1877,"Q",0.1429,"S",0.2276,"D",-0.2256,"F",-0.0012,"H",-0.100092,"K",0.22875,"B",-0.1962,"Z",-0.0362506,"X",-0.0362506,"?"),nrow=24,ncol=2,byrow=TRUE)

#Computing sequence length without non-aminoacid characters.
fixed.protein.sequence <- character()
for (i in 1:protein.length) {
	if (protein.sequence[i] %in% IUPAC.degeneracy) fixed.protein.sequence <- c(fixed.protein.sequence,protein.sequence[i])
	}
protein.length <- length(fixed.protein.sequence)

#Setting sliding windows centers.
centering <- (((stop-start+1)-2*(wsize%/%2))%%wstep) %/% 2
xs <- start+centering+wsize%/%2
i <- xs + wstep
while (i < stop-wsize%/%2) {
	xs <- c(xs,i)
	i <- i + wstep
	}

#Setting analysis for criterion = "all" (all six criteria are being used).
par(mfrow=c(1,1))
criterion.check <- criterion
if (criterion == "all") {
	par(mfrow=c(2,4))
	criterion <- criteria
	}

#Using charge criterion (acidic | basic).
if ("charge" %in% criterion) {
	count.exp.acidic <- 0
	count.exp.basic <- 0
	count.exp.NC <- 0
	for (i in 1:protein.length) {
		ifelse(fixed.protein.sequence[i] %in% acidic,count.exp.acidic <- count.exp.acidic+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.acidic <- count.exp.acidic+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% basic,count.exp.basic <- count.exp.basic+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.basic <- count.exp.basic+0.15,NA)
		ifelse(fixed.protein.sequence[i] %in% NC,count.exp.NC <- count.exp.NC+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.NC <- count.exp.NC+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.NC <- count.exp.NC+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.NC <- count.exp.NC+0.75,NA)
		}
	freq.exp.acidic <- count.exp.acidic/protein.length
	freq.exp.basic <- count.exp.basic/protein.length
	freq.exp.NC <- count.exp.NC/protein.length
	p.values <- numeric()
	local.acidic <- numeric()
	local.basic <- numeric()
	local.NC <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		freq.obs.acidic <- 0
		freq.obs.basic <- 0
		freq.obs.NC <- 0
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			ifelse(protein.sequence[i] %in% acidic,freq.obs.acidic <- freq.obs.acidic+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.acidic <- freq.obs.acidic+0.5,NA)
			ifelse(protein.sequence[i] == "Z",freq.obs.acidic <- freq.obs.acidic+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.acidic <- freq.obs.acidic+0.1,NA)
			ifelse(protein.sequence[i] %in% basic,freq.obs.basic <- freq.obs.basic+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.basic <- freq.obs.basic+0.15,NA)
			ifelse(protein.sequence[i] %in% NC,freq.obs.NC <- freq.obs.NC+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.NC <- freq.obs.NC+0.5,NA)
			ifelse(protein.sequence[i] == "Z",freq.obs.NC <- freq.obs.NC+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.NC <- freq.obs.NC+0.75,NA)
			}
		if (freq.obs.acidic+freq.obs.basic+freq.obs.NC == 0) {
			p.values <- c(p.values,NA)
			local.acidic <- c(local.acidic,NA)
			local.basic <- c(local.basic,NA)
			local.NC <- c(local.NC,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			chi.square.results <- chisq.test(c(freq.obs.acidic,freq.obs.basic,freq.obs.NC),p=c(freq.exp.acidic,freq.exp.basic,freq.exp.NC),simulate.p.value=TRUE,B=replicates)
			p.values <- c(p.values,chi.square.results$p.value)
			local.acidic <- c(local.acidic,freq.obs.acidic/(2*(wsize%/%2)+1))
			local.basic <- c(local.basic,freq.obs.basic/(2*(wsize%/%2)+1))
			local.NC <- c(local.NC,freq.obs.NC/(2*(wsize%/%2)+1))
			ifelse(chi.square.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(chi.square.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using charge criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],acidic=local.acidic[p.values<alpha],basic=local.basic[p.values<alpha],NC=local.NC[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_charge.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",freq.exp.acidic,freq.exp.basic,freq.exp.NC,sep="\t"),file=paste("../",outfile,"_charge.out",sep=""),append=TRUE)
	}

#Using hydropathy criterion (acidic | basic | polar uncharged | hydrophobic nonpolar)
if ("hydropathy" %in% criterion) {
	count.exp.acidic <- 0
	count.exp.basic <- 0
	count.exp.polar.uncharged <- 0
	count.exp.hydrophobic.nonpolar <- 0
	for (i in 1:protein.length) {
		ifelse(fixed.protein.sequence[i] %in% acidic,count.exp.acidic <- count.exp.acidic+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.acidic <- count.exp.acidic+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% basic,count.exp.basic <- count.exp.basic+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.basic <- count.exp.basic+0.15,NA)
		ifelse(fixed.protein.sequence[i] %in% polar.uncharged,count.exp.polar.uncharged <- count.exp.polar.uncharged+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.polar.uncharged <- count.exp.polar.uncharged+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.polar.uncharged <- count.exp.polar.uncharged+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.polar.uncharged <- count.exp.polar.uncharged+0.35,NA)
		ifelse(fixed.protein.sequence[i] %in% hydrophobic.nonpolar,count.exp.hydrophobic.nonpolar <- count.exp.hydrophobic.nonpolar+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.hydrophobic.nonpolar <- count.exp.hydrophobic.nonpolar+0.4,NA)
		}
	freq.exp.acidic <- count.exp.acidic/protein.length
	freq.exp.basic <- count.exp.basic/protein.length
	freq.exp.polar.uncharged <- count.exp.polar.uncharged/protein.length
	freq.exp.hydrophobic.nonpolar <- count.exp.hydrophobic.nonpolar/protein.length
	p.values <- numeric()
	local.acidic <- numeric()
	local.basic <- numeric()
	local.polar.uncharged <- numeric()
	local.hydrophobic.nonpolar <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		freq.obs.acidic <- 0
		freq.obs.basic <- 0
		freq.obs.polar.uncharged <- 0
		freq.obs.hydrophobic.nonpolar <- 0
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			ifelse(protein.sequence[i] %in% acidic,freq.obs.acidic <- freq.obs.acidic+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.acidic <- freq.obs.acidic+0.5,NA)	
			ifelse(protein.sequence[i] == "Z",freq.obs.acidic <- freq.obs.acidic+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.acidic <- freq.obs.acidic+0.1,NA)
			ifelse(protein.sequence[i] %in% basic,freq.obs.basic <- freq.obs.basic+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.basic <- freq.obs.basic+0.15,NA)
			ifelse(protein.sequence[i] %in% polar.uncharged,freq.obs.polar.uncharged <- freq.obs.polar.uncharged+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.polar.uncharged <- freq.obs.polar.uncharged+0.5,NA)
			ifelse(protein.sequence[i] == "Z",freq.obs.polar.uncharged <- freq.obs.polar.uncharged+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.polar.uncharged <- freq.obs.polar.uncharged+0.35,NA)
			ifelse(protein.sequence[i] %in% hydrophobic.nonpolar,freq.obs.hydrophobic.nonpolar <- freq.obs.hydrophobic.nonpolar+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.hydrophobic.nonpolar <- freq.obs.hydrophobic.nonpolar+0.4,NA)
			}
		if (freq.obs.acidic+freq.obs.basic+freq.obs.polar.uncharged+freq.obs.hydrophobic.nonpolar == 0) {
			p.values <- c(p.values,NA)
			local.acidic <- c(local.acidic,NA)
			local.basic <- c(local.basic,NA)
			local.polar.uncharged <- c(local.polar.uncharged,NA)
			local.hydrophobic.nonpolar <- c(local.hydrophobic.nonpolar,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			chi.square.results <- chisq.test(c(freq.obs.acidic,freq.obs.basic,freq.obs.polar.uncharged,freq.obs.hydrophobic.nonpolar),p=c(freq.exp.acidic,freq.exp.basic,freq.exp.polar.uncharged,freq.exp.hydrophobic.nonpolar),simulate.p.value=TRUE,B=replicates)
			p.values <- c(p.values,chi.square.results$p.value)
			local.acidic <- c(local.acidic,freq.obs.acidic/(2*(wsize%/%2)+1))
			local.basic <- c(local.basic,freq.obs.basic/(2*(wsize%/%2)+1))
			local.polar.uncharged <- c(local.polar.uncharged,freq.obs.polar.uncharged/(2*(wsize%/%2)+1))
			local.hydrophobic.nonpolar <- c(local.hydrophobic.nonpolar,freq.obs.hydrophobic.nonpolar/(2*(wsize%/%2)+1))
			ifelse(chi.square.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(chi.square.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using hydropathy criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],acidic=local.acidic[p.values<alpha],basic=local.basic[p.values<alpha],polar.uncharged=local.polar.uncharged[p.values<alpha],hydrophobic.nonpolar=local.hydrophobic.nonpolar[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_hydropathy.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",freq.exp.acidic,freq.exp.basic,freq.exp.polar.uncharged,freq.exp.hydrophobic.nonpolar,sep="\t"),file=paste("../",outfile,"_hydropathy.out",sep=""),append=TRUE)
	}

#Using structure criterion (external | internal | ambivalent).
if ("structure" %in% criterion) {
	count.exp.external <- 0
	count.exp.internal <- 0
	count.exp.ambivalent <- 0
	for (i in 1:protein.length) {
		ifelse(fixed.protein.sequence[i] %in% external,count.exp.external <- count.exp.external+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.external <- count.exp.external+1,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.external <- count.exp.external+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.external <- count.exp.external+0.35,NA)
		ifelse(fixed.protein.sequence[i] %in% internal,count.exp.internal <- count.exp.internal+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.internal <- count.exp.internal+0.25,NA)
		ifelse(fixed.protein.sequence[i] %in% ambivalent,count.exp.ambivalent <- count.exp.ambivalent+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.ambivalent <- count.exp.ambivalent+0.4,NA)
		}
	freq.exp.external <- count.exp.external/protein.length
	freq.exp.internal <- count.exp.internal/protein.length
	freq.exp.ambivalent <- count.exp.ambivalent/protein.length
	p.values <- numeric()
	local.external <- numeric()
	local.internal <- numeric()
	local.ambivalent <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		freq.obs.external <- 0
		freq.obs.internal <- 0
		freq.obs.ambivalent <- 0
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			ifelse(protein.sequence[i] %in% external,freq.obs.external <- freq.obs.external+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.external <- freq.obs.external+1,NA)
			ifelse(protein.sequence[i] == "Z",freq.obs.external <- freq.obs.external+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.external <- freq.obs.external+0.35,NA)
			ifelse(protein.sequence[i] %in% internal,freq.obs.internal <- freq.obs.internal+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.internal <- freq.obs.internal+0.25,NA)
			ifelse(protein.sequence[i] %in% ambivalent,freq.obs.ambivalent <- freq.obs.ambivalent+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.ambivalent <- freq.obs.ambivalent+0.4,NA)
			}
		if (freq.obs.external+freq.obs.internal+freq.obs.ambivalent == 0) {
			p.values <- c(p.values,NA)
			local.external <- c(local.external,NA)
			local.internal <- c(local.internal,NA)
			local.ambivalent <- c(local.ambivalent,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			chi.square.results <- chisq.test(c(freq.obs.external,freq.obs.internal,freq.obs.ambivalent),p=c(freq.exp.external,freq.exp.internal,freq.exp.ambivalent),simulate.p.value=TRUE,B=replicates)
			p.values <- c(p.values,chi.square.results$p.value)
			local.external <- c(local.external,freq.obs.external/(2*(wsize%/%2)+1))
			local.internal <- c(local.internal,freq.obs.internal/(2*(wsize%/%2)+1))
			local.ambivalent <- c(local.ambivalent,freq.obs.ambivalent/(2*(wsize%/%2)+1))
			ifelse(chi.square.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(chi.square.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using structure criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],external=local.external[p.values<alpha],internal=local.internal[p.values<alpha],ambivalent=local.ambivalent[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_structure.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",freq.exp.external,freq.exp.internal,freq.exp.ambivalent,sep="\t"),file=paste("../",outfile,"_structure.out",sep=""),append=TRUE)
	}

#Using chemical (rasmol) funcionality (acidic | aliphatic | aliphatic (small) | amide | aromatic | basic | hydroxhyl | imino | sulphur).
if ("chemical" %in% criterion) {
	count.exp.acidic <- 0
	count.exp.aliphatic <- 0
	count.exp.aliphatic.small <- 0
	count.exp.amide <- 0
	count.exp.aromatic <- 0
	count.exp.basic <- 0
	count.exp.hydroxyl <- 0
	count.exp.imino <- 0
	count.exp.sulfur <- 0
	for (i in 1:protein.length) {
		ifelse(fixed.protein.sequence[i] %in% acidic,count.exp.acidic <- count.exp.acidic+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.acidic <- count.exp.acidic+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.acidic <- count.exp.acidic+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% aliphatic,count.exp.aliphatic <- count.exp.aliphatic+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.aliphatic <- count.exp.aliphatic+0.15,NA)
		ifelse(fixed.protein.sequence[i] %in% aliphatic.small,count.exp.aliphatic.small <- count.exp.aliphatic.small+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.aliphatic.small <- count.exp.aliphatic.small+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% amide,count.exp.amide <- count.exp.amide+1,NA)
		ifelse(fixed.protein.sequence[i] == "B",count.exp.amide <- count.exp.amide+0.5,NA)
		ifelse(fixed.protein.sequence[i] == "Z",count.exp.amide <- count.exp.amide+0.5,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.amide <- count.exp.amide+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% aromatic,count.exp.aromatic <- count.exp.aromatic+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.aromatic <- count.exp.aromatic+0.15,NA)
		ifelse(fixed.protein.sequence[i] %in% basic,count.exp.basic <- count.exp.basic+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.basic <- count.exp.basic+0.15,NA)
		ifelse(fixed.protein.sequence[i] %in% hydroxyl,count.exp.hydroxyl <- count.exp.hydroxyl+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.hydroxyl <- count.exp.hydroxyl+0.1,NA)
		ifelse(fixed.protein.sequence[i] %in% imino,count.exp.imino <- count.exp.imino+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.imino <- count.exp.imino+0.05,NA)
		ifelse(fixed.protein.sequence[i] %in% sulfur,count.exp.sulfur <- count.exp.sulfur+1,NA)
		ifelse(fixed.protein.sequence[i] %in% unknown,count.exp.sulfur <- count.exp.sulfur+0.1,NA)
		}
	freq.exp.acidic <- count.exp.acidic/protein.length
	freq.exp.aliphatic <- count.exp.aliphatic/protein.length
	freq.exp.aliphatic.small <- count.exp.aliphatic.small/protein.length
	freq.exp.amide <- count.exp.amide/protein.length
	freq.exp.aromatic <- count.exp.aromatic/protein.length
	freq.exp.basic <- count.exp.basic/protein.length
	freq.exp.hydroxyl <- count.exp.hydroxyl/protein.length
	freq.exp.imino <- count.exp.imino/protein.length
	freq.exp.sulfur <- count.exp.sulfur/protein.length
	p.values <- numeric()
	local.acidic <- numeric()
	local.aliphatic <- numeric()
	local.aliphatic.small <- numeric()
	local.amide <- numeric()
	local.aromatic <- numeric()
	local.basic <- numeric()
	local.hydroxyl <- numeric()
	local.imino <- numeric()
	local.sulfur <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		freq.obs.acidic <- 0
		freq.obs.aliphatic <- 0
		freq.obs.aliphatic.small <- 0
		freq.obs.amide <- 0
		freq.obs.aromatic <- 0
		freq.obs.basic <- 0
		freq.obs.hydroxyl <- 0
		freq.obs.imino <- 0
		freq.obs.sulfur <- 0
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			ifelse(protein.sequence[i] %in% acidic,freq.obs.acidic <- freq.obs.acidic+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.acidic <- freq.obs.acidic+0.5,NA)	
			ifelse(protein.sequence[i] == "Z",freq.obs.acidic <- freq.obs.acidic+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.acidic <- freq.obs.acidic+0.1,NA)
			ifelse(protein.sequence[i] %in% aliphatic,freq.obs.aliphatic <- freq.obs.aliphatic+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.aliphatic <- freq.obs.aliphatic+0.15,NA)
			ifelse(protein.sequence[i] %in% aliphatic.small,freq.obs.aliphatic.small <- freq.obs.aliphatic.small+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.aliphatic.small <- freq.obs.aliphatic.small+0.1,NA)
			ifelse(protein.sequence[i] %in% amide,freq.obs.amide <- freq.obs.amide+1,NA)
			ifelse(protein.sequence[i] == "B",freq.obs.amide <- freq.obs.amide+0.5,NA)
			ifelse(protein.sequence[i] == "Z",freq.obs.amide <- freq.obs.amide+0.5,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.amide <- freq.obs.amide+0.1,NA)
			ifelse(protein.sequence[i] %in% aromatic,freq.obs.aromatic <- freq.obs.aromatic+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.aromatic <- freq.obs.aromatic+0.15,NA)
			ifelse(protein.sequence[i] %in% basic,freq.obs.basic <- freq.obs.basic+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.basic <- freq.obs.basic+0.15,NA)
			ifelse(protein.sequence[i] %in% hydroxyl,freq.obs.hydroxyl <- freq.obs.hydroxyl+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.hydroxyl <- freq.obs.hydroxyl+0.1,NA)
			ifelse(protein.sequence[i] %in% imino,freq.obs.imino <- freq.obs.imino+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.imino <- freq.obs.imino+0.05,NA)
			ifelse(protein.sequence[i] %in% sulfur,freq.obs.sulfur <- freq.obs.sulfur+1,NA)
			ifelse(protein.sequence[i] %in% unknown,freq.obs.sulfur <- freq.obs.sulfur+0.1,NA)
			}
		if (freq.obs.acidic+freq.obs.aliphatic+freq.obs.aliphatic.small+freq.obs.amide+freq.obs.aromatic+freq.obs.basic+freq.obs.hydroxyl+freq.obs.imino+freq.obs.sulfur == 0) {
			p.values <- c(p.values,NA)
			local.acidic <- c(local.acidic,NA)
			local.aliphatic <- c(local.aliphatic,NA)
			local.aliphatic.small <- c(local.aliphatic.small,NA)
			local.amide <- c(local.amide,NA)
			local.aromatic <- c(local.aromatic,NA)
			local.basic <- c(local.basic,NA)
			local.hydroxyl <- c(local.hydroxyl,NA)
			local.imino <- c(local.imino,NA)
			local.sulfur <- c(local.sulfur,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			chi.square.results <- chisq.test(c(freq.obs.acidic,freq.obs.aliphatic,freq.obs.aliphatic.small,freq.obs.amide,freq.obs.aromatic,freq.obs.basic,freq.obs.hydroxyl,freq.obs.imino,freq.obs.sulfur),p=c(freq.exp.acidic,freq.exp.aliphatic,freq.exp.aliphatic.small,freq.exp.amide,freq.exp.aromatic,freq.exp.basic,freq.exp.hydroxyl,freq.exp.imino,freq.exp.sulfur),simulate.p.value=TRUE,B=replicates)
			p.values <- c(p.values,chi.square.results$p.value)
			local.acidic <- c(local.acidic,freq.obs.acidic/(2*(wsize%/%2)+1))
			local.aliphatic <- c(local.aliphatic,freq.obs.aliphatic/(2*(wsize%/%2)+1))
			local.aliphatic.small <- c(local.aliphatic.small,freq.obs.aliphatic.small/(2*(wsize%/%2)+1))
			local.amide <- c(local.amide,freq.obs.amide/(2*(wsize%/%2)+1))
			local.aromatic <- c(local.aromatic,freq.obs.aromatic/(2*(wsize%/%2)+1))
			local.basic <- c(local.basic,freq.obs.basic/(2*(wsize%/%2)+1))
			local.hydroxyl <- c(local.hydroxyl,freq.obs.hydroxyl/(2*(wsize%/%2)+1))
			local.imino <- c(local.imino,freq.obs.imino/(2*(wsize%/%2)+1))
			local.sulfur <- c(local.sulfur,freq.obs.sulfur/(2*(wsize%/%2)+1))
			ifelse(chi.square.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(chi.square.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using chemical criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],acidic=local.acidic[p.values<alpha],aliphatic=local.aliphatic[p.values<alpha],aliphatic.small=local.aliphatic.small[p.values<alpha],amide=local.amide[p.values<alpha],aromatic=local.aromatic[p.values<alpha],basic=local.basic[p.values<alpha],hydroxyl=local.hydroxyl[p.values<alpha],imino=local.imino[p.values<alpha],sulfur=local.sulfur[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_chemical.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",freq.exp.acidic,freq.exp.aliphatic,freq.exp.aliphatic.small,freq.exp.amide,freq.exp.aromatic,freq.exp.basic,freq.exp.hydroxyl,freq.exp.imino,freq.exp.sulfur,sep="\t"),file=paste("../",outfile,"_chemical.out",sep=""),append=TRUE)
	}

#Using standard area criterion (as taken from texshade).
if ("standard area" %in% criterion) {
	sum.standard.areas <- 0
	for (i in 1:protein.length) sum.standard.areas <- sum.standard.areas + as.numeric(standard.area[standard.area[,2] == fixed.protein.sequence[i],1])
	mean.standard.area <- sum.standard.areas/protein.length
	p.values <- numeric()
	local.mean.standard.area <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		obs.standard.areas <- numeric()
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			if (protein.sequence[i] %in% IUPAC.degeneracy) {
				obs.standard.areas <- c(obs.standard.areas,as.numeric(standard.area[standard.area[,2] == protein.sequence[i],1]))
				}
			}
		if (length(obs.standard.areas) %in% c(0,1)) {
			p.values <- c(p.values,NA)
			local.mean.standard.area <- c(local.mean.standard.area,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else if (sd(obs.standard.areas) == 0) {
			p.values <- c(p.values,NA)
			local.mean.standard.area <- c(local.mean.standard.area,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			wilcox.test.results <- wilcox.test(obs.standard.areas,mu=mean.standard.area,exact=FALSE)
			p.values <- c(p.values,wilcox.test.results$p.value)
			local.mean.standard.area <- c(local.mean.standard.area,mean(obs.standard.areas))
			ifelse(wilcox.test.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(wilcox.test.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using standard area criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],mean.standard.area=local.mean.standard.area[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_standard_area.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",mean.standard.area,sep="\t"),file=paste("../",outfile,"_standard_area.out",sep=""),append=TRUE)
	}

#Using accessible area criterion (as taken from texshade).
if ("accessible area" %in% criterion) {
	sum.accessible.areas <- 0
	for (i in 1:protein.length) sum.accessible.areas <- sum.accessible.areas + as.numeric(accessible.area[accessible.area[,2] == fixed.protein.sequence[i],1])
	mean.accessible.area <- sum.accessible.areas/protein.length
	p.values <- numeric()
	local.mean.accessible.area <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		obs.accessible.areas <- numeric()
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			if (protein.sequence[i] %in% IUPAC.degeneracy) {
				obs.accessible.areas <- c(obs.accessible.areas,as.numeric(accessible.area[accessible.area[,2] == protein.sequence[i],1]))
				}
			}
		if (length(obs.accessible.areas) %in% c(0,1)) {
			p.values <- c(p.values,NA)
			local.mean.accessible.area <- c(local.mean.accessible.area,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else if (sd(obs.accessible.areas) == 0) {
			p.values <- c(p.values,NA)
			local.mean.accessible.area <- c(local.mean.accessible.area,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			wilcox.test.results <- wilcox.test(obs.accessible.areas,mu=mean.accessible.area,exact=FALSE)
			p.values <- c(p.values,wilcox.test.results$p.value)
			local.mean.accessible.area <- c(local.mean.accessible.area,mean(obs.accessible.areas))
			ifelse(wilcox.test.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(wilcox.test.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using accessible area criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],mean.accessible.area=local.mean.accessible.area[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_accessible_area.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",mean.accessible.area,sep="\t"),file=paste("../",outfile,"_accessible_area.out",sep=""),append=TRUE)
	}

#Using loop propensities criterion (as taken from Via A, Rother K, Tramontano A (2014), "Managing Your Biological Data with Python", Taylor & Francis, Boca Raton, p.81).
if ("loop propensity" %in% criterion) {
	sum.loop.propensities <- 0
	for (i in 1:protein.length) sum.loop.propensities <- sum.loop.propensities + as.numeric(loop.propensities[loop.propensities[,2] == fixed.protein.sequence[i],1])
	mean.loop.propensity <- sum.loop.propensities/protein.length
	p.values <- numeric()
	local.mean.loop.propensity <- numeric()
	plot.cexs <- numeric()
	plot.colors <- character()
	for (w in 1:length(xs)) {
		obs.loop.propensities <- numeric()
		for (i in (xs[w]-wsize%/%2):(xs[w]+wsize%/%2)) {
			if (protein.sequence[i] %in% IUPAC.degeneracy) {
				obs.loop.propensities <- c(obs.loop.propensities,as.numeric(loop.propensities[loop.propensities[,2] == protein.sequence[i],1]))
				}
			}
		if (length(obs.loop.propensities) %in% c(0,1)) {
			p.values <- c(p.values,NA)
			local.mean.loop.propensity <- c(local.mean.loop.propensity,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else if (sd(obs.loop.propensities) == 0) {
			p.values <- c(p.values,NA)
			local.mean.loop.propensity <- c(local.mean.loop.propensity,NA)
			plot.cexs <- c(plot.cexs,1)
			plot.colors <- c(plot.colors,"#000000")
			}
			else {
			wilcox.test.results <- wilcox.test(obs.loop.propensities,mu=mean.loop.propensity,exact=FALSE)
			p.values <- c(p.values,wilcox.test.results$p.value)
			local.mean.loop.propensity <- c(local.mean.loop.propensity,mean(obs.loop.propensities))
			ifelse(wilcox.test.results$p.value<alpha,plot.cexs <- c(plot.cexs,1),plot.cexs <- c(plot.cexs,0.5))
			ifelse(wilcox.test.results$p.value<alpha,plot.colors <- c(plot.colors,"#FF0000"),plot.colors <- c(plot.colors,"#000000"))
			}
		}
	plot(xs,p.values,type="p",xlim=c(start,stop),ylim=c(0,1),main="Using loop propensity criterion",ylab="p-value",xlab=paste("Residue (sliding window size: ",wsize,"; window step: ",wstep,")",sep=""),pch=21,col=plot.colors,bg=plot.colors,cex=plot.cexs)
	abline(h=alpha)
	results.data.frame <- data.frame(start=xs[p.values<alpha]-wsize%/%2,stop=xs[p.values<alpha]+wsize%/%2,p.value=p.values[p.values<alpha],mean.loop.propensity=local.mean.loop.propensity[p.values<alpha])
	write.table(results.data.frame,file=paste("../",outfile,"_loop_propensity.out",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat(paste("Whole protein\t\t",mean.loop.propensity,sep="\t"),file=paste("../",outfile,"_loop_propensity.out",sep=""),append=TRUE)
	}

if (!criterion.check %in% c(criteria,"all")) message("The selected criterion is not available.")

#Saving the output to outfile.
if (criterion.check %in% c(criteria,"all")) dev.copy2pdf(file=paste("../",outfile,".pdf",sep=""))
