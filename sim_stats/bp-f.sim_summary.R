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
the<-0.0020053876480760741*nsites
rho<-4*(nsites-1)*15000*2*1e-8

scan(pipe(paste("~/bin/msdir_2/ms",
"192 20000", 
"-t ", the, 
"-r", rho, nsites,  
"-I 2 96 96",
"-n 1 9.99999",
"-n 2 0.485873",
"-em 0 1 2 0.002",
"-em 0 2 1 0.030756",
"-ej 0.01011685 2 1",
"-en 0.01011685 1 7.11875",
"-en 4.735732 1 1",
sep=" ")),what="character",sep="\n")->junk

msparse.2pops(junk,n1=96,n2=96,nsites=5000)->testBPF
testBPF<-cbind(testBPF,testBPF[,4]-testBPF[,1])
testBPF<-cbind(testBPF,testBPF[,6]-testBPF[,3])
testBPF[,8]<-(testBPF[,8]-mean(testBPF[,8]))/sd(testBPF[,8])
testBPF[,9]<-(testBPF[,9]-mean(testBPF[,9]))/sd(testBPF[,9])

write.table(testBPF,"testBPF.txt")
