library(hierfstat)

msparse.2pops<-function(x,n1,n2,nsites){
	
	he<-grep("segsites", x)
	ndat<-length(he)
	ind<-he[2]-he[1]-3
	
	out1pi<-vector("list", ndat)
	out1tw<-vector("list", ndat)
	out1tj<-vector("list", ndat)
	out2pi<-vector("list", ndat)
	out2tw<-vector("list", ndat)
	out2tj<-vector("list", ndat)
	outfst<-vector("list", ndat)

	
 	piest<-function(x){
		
		inds<-dim(x)[1]
		sum(colSums(x)*(inds-colSums(x))/choose(inds,2))
		
		}
	
	twest<-function(x){
		
		inds<-dim(x)[1]
		sum(colSums(x)>0)/sum(1/(1:(inds-1)))
		
		}
	
	fstest<-function(x,num1,num2){
		x<-x+10
		pops<-c(rep(1,num1),rep(2,num2))
		x<-cbind(pops,x)
		wc(x,diploid=FALSE)$FST
		}
	
	for(i in 1:ndat){
		
		if(x[he[i]]=="segsites: 0"){
			out1pi[[i]]<-0
			out1tw[[i]]<-0
			out1tj[[i]]<-0
			out2pi[[i]]<-0
			out2tw[[i]]<-0
			out2tj[[i]]<-0
			outfst[[i]]<-0

			}
			else{
		
				outmat<-x[(he[i]+2):(he[i]+1+ind)]
				outmat<-do.call(rbind,strsplit(outmat,""))
				class(outmat)	<-"numeric"
				out1pi[[i]]<-piest(outmat[1:n1,])
				out1tw[[i]]<-twest(outmat[1:n1,])
				out1tj[[i]]<-tajimas(out1pi[[i]],out1tw[[i]],n1)
				out2pi[[i]]<-piest(outmat[(n1+1):(n1+n2),])
				out2tw[[i]]<-twest(outmat[(n1+1):(n1+n2),])
				out2tj[[i]]<-tajimas(out2pi[[i]],out2tw[[i]],n2)
				outfst[[i]]<-fstest(outmat,n1,n2)
					}
		if(ndat%%10==0){print(i)}
		}	
	
	cbind(unlist(out1pi),unlist(out1tw),unlist(out1tj),unlist(out2pi),unlist(out2tw),unlist(out2tj),unlist(outfst))

	}


tajimas<-function(pi,tw,n){
	
	if(any(is.na(c(pi,tw)))){
		return(NA)
		}
	
	nm<-n-1
	a1<-sum(1/1:nm)
	a2<-sum(1/(1:nm)^2)
	b1<-(n+1)/(3*(n-1))
	b2<-(2*((n^2)+n+3))/(9*n*(n-1))
	c1<-b1-(1/a1)
	c2<-b2-((n+2)/(a1*n))+a2/(a1^2)
	e1<-c1/a1
	e2<-c2/((a1^2)+a2)
	
	
	
	S<-tw*a1
	num<-pi-tw
	denom<-sqrt((e1*S)+(e2*S*(S-1)))
	
	tajd<-num/denom
	
	return(tajd)
	
	}
	

nsites<-5000
the<-0.017213686777507266*nsites
rho<-4*(nsites-1)*15000*2*1e-8

scan(pipe(paste("~/bin/msdir_2/ms",
"192 20000", 
"-t ", the, 
"-r", rho, nsites,  
"-I 2 96 96",
"-n 1 0.0918114",
"-n 2 0.0949811",
"-em 0 2 1 4.31936",
"-em 0 1 2 2.11734",
"-ej 0.0005 2 1",
"-en 0.0005 1 0.0735887",
"-en 0.001 1 1",
sep=" ")),what="character",sep="\n")->junk

msparse.2pops(junk,n1=96,n2=96,nsites=5000)->testERKC
testERKC<-cbind(testERKC,testERKC[,4]-testERKC[,1])
testERKC<-cbind(testERKC,testERKC[,6]-testERKC[,3])
testERKC[,8]<-(testERKC[,8]-mean(testERKC[,8]))/sd(testERKC[,8])
testERKC[,9]<-(testERKC[,9]-mean(testERKC[,9]))/sd(testERKC[,9])

write.table(testERKC,"testERKC.txt")
