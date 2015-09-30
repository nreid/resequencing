#this script takes a gff file and spits out bed files containing 1st, 2nd, 3rd codon positions for the whole genome. 
#it uses what is presumably the gff file uploaded to NCBI, then backtranslated ("liftdown") to the version of the genome I've been using. 
#probably a much better way to do this...

library(dplyr)
library(stringr)
library(magrittr)

read.table("../killifish20130322asm.fai",stringsAsFactors=FALSE)->fai
read.table("../kfish2rae5h_fc17_subupd4.gff3.liftdown.gz",stringsAsFactors=FALSE,sep="\t",quote="")->gff
gff[,8]<-as.numeric(gff[,8])
gffm<- gff[grep("^ID=Funhe[^;]+t1;",gff[,9]),]
cbind(gffm,str_extract(gffm[,9],perl("(?<=ID=)[^;]+")))->gffm
cbind(gffm,str_extract(gffm[,9],"Map:[^,]+"))->gffm
cbind(gffm,str_extract(gffm[,9],"Express:[^,]+"))->gffm
cbind(gffm,str_extract(gffm[,9],"Homology:[^,]+"))->gffm
cbind(gffm,str_extract(gffm[,9],"Intron:[^,]+"))->gffm
cbind(gffm,str_extract(gffm[,9],"Protein:[^;]+"))->gffm
cbind(gffm,str_extract(gffm[,9],"Class:[^,]+"))->gffm
names(table(gffm[,10])[table(gffm[,10])==1])->sub
gffm2<-gffm[gffm[,10]%in%sub,]

filter(gff, V3=="CDS",grepl("t1;",gff[,9])) ->gffc
cbind(gffc,str_extract(gffc[,9],"Funh[^;]+t1"))->gffc
gffc[gffc[,5]-gffc[,4]>10,]->gffc

gffsub<- gff[grep("^Parent=Funhe[^;]+t1;",gff[,9]),]
gffsub<-cbind(gffsub,gsub(";.*","",gffsub[,9]))
colnames(gffsub)[10]<-"V10"
table(gffsub[,10],gffsub[,3])->fcount
cbind(rownames(fcount),fcount)->fcount

read.table("../kfish2asm.ncbiscafids",stringsAsFactors=FALSE)->scaftrans

#scaftrans[-(6429:6430),]->scaftrans


cat("",file="fiveprime.bed")
cat("",file="threeprime.bed")

for(i in 1:dim(fcount)[1]){
	
	if(fcount[i,2]==0|fcount[i,3]==0){print("bailing 1");next}
	
	g1<-filter(gffsub,V10==fcount[i,1])
	
	if(length(unique(g1[,1]))>1){print("bailing 2");next}
	
	if(!(g1[1,7]=="+"|g1[1,7]=="-")){print("bailing 3");next}
	
	g2<-g1[g1[,3]=="CDS",]
	g3<-g1[g1[,3]=="exon",]
	
	cdsbound<-range(c(g2[,4],g2[,5]))
	exopos<-c()
	for(j in 1:length(g3[,1])){
		exopos<-c(exopos,g3[j,4]:g3[j,5])
		}
	
	if(g1[1,7]=="+"){
		
		fiveprime<-exopos[exopos<cdsbound[1]]
		threeprime<-exopos[exopos>cdsbound[2]]
		
		if(length(threeprime)){
			t1<-cbind(g1[1,1],threeprime-1,threeprime)
			write.table(t1,file="threeprime.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")
			}
		if(length(fiveprime)){
			t2<-cbind(g1[1,1],fiveprime-1,fiveprime)
			write.table(t2,file="fiveprime.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")			
			}
			
		}

	if(g1[1,7]=="-"){
		
		threeprime<-exopos[exopos<cdsbound[1]]
		fiveprime<-exopos[exopos>cdsbound[2]]

		if(length(threeprime)){
			t1<-cbind(g1[1,1],threeprime-1,threeprime)
			write.table(t1,file="threeprime.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")			
			}
		if(length(fiveprime)){
			t2<-cbind(g1[1,1],fiveprime-1,fiveprime)
			write.table(t2,file="fiveprime.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")						
			}
			
		}
	
	if((i%%100)==0){
		cat(i," ")
		}
	}

	
	





cat("",file="c1.bed")
cat("",file="c2.bed")
cat("",file="c3.bed")

c1out<-c()
c2out<-c()
c3out<-c()
for(i in 1:length(gffc[,1])){
	if(gffc[i,7]=="+"){
		c1<-seq(from=gffc[i,4]+gffc[i,8],to=gffc[i,5],by=3)
		c1<-cbind(gffc[i,1],c1-1,c1)
		#c1out<-rbind(c1out,c1)
		write.table(c1,file="c1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")
		
		c2<-seq(from=gffc[i,4]+gffc[i,8]+1,to=gffc[i,5],by=3)
		c2<-cbind(gffc[i,1],c2-1,c2)
		#c2out<-rbind(c2out,c2)
		write.table(c2,file="c2.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")
		
		c3<-seq(from=gffc[i,4]+gffc[i,8]+2,to=gffc[i,5],by=3)
		c3<-cbind(gffc[i,1],c3-1,c3)
		#c3out<-rbind(c3out,c3)
		write.table(c3,file="c3.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")
		}
	if(gffc[i,7]=="-"){
		c1<-seq(from=gffc[i,5]-gffc[i,8],to=gffc[i,4],by=-3)
		c1<-rev(c1)
		c1<-cbind(gffc[i,1],c1-1,c1)
		#c1out<-rbind(c1out,c1)
		write.table(c1,file="c1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")

		c2<-seq(from=gffc[i,5]-gffc[i,8]-1,to=gffc[i,4],by=-3)
		c2<-rev(c2)
		c2<-cbind(gffc[i,1],c2-1,c2)
		#c2out<-rbind(c2out,c2)
		write.table(c2,file="c2.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")

		c3<-seq(from=gffc[i,5]-gffc[i,8]-2,to=gffc[i,4],by=-3)
		c3<-rev(c3)
		c3<-cbind(gffc[i,1],c3-1,c3)
		#c3out<-rbind(c3out,c3)
		write.table(c3,file="c3.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep="\t")

		}
	if(i%%1000==0){
		cat(i," ")
		}
	}






