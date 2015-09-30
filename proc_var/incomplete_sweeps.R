library(mscr)
library(pvclust)
library(cluster)
#library(NMF)

#this is a coarse attempt to look for incomplete sweeps in populations sequenced at very low coverage (0.5x)
#this script tries to split individuals into two clusters based on their genotypes in sliding windows across the genome. 
#it then calculates the average distance within and between clusters and the ratios of within/between

#read in 50kb/10kb sliding window bed file. 
win <- read.table("~/popgen/kfish3/50kb10kb.windows.bed",stringsAsFactors=FALSE)
win[,2] <- as.character(win[,2])
win[,3] <- as.character(win[,3])
win <- as.matrix(win)

#output file name
fname <- "~/popgen/variants/bowfree/ALL1/incomplete.sweeps.50kb10kb.txt"

#initialize output
cat("\\#incomplete sweep scan\n",file=fname)

#populations to analyze
pops <- c("BI","NBH","F-","BP","SH","NYC","KC","ER")

#loop through windows
for(i in 1:length(win[,1])){
	
	dati <- c()
	wini <- win[i,]

	dati <- wini

	#tabix format for window
	tabix<-paste(wini[1],":",wini[2],"-",wini[3],sep="")

	#pull genotypes for region, reformat them. return NA for window if none exist. 
	g0 <- NULL
	try(
		g0 <- pullgenos(tabix)
		)
	if(is.null(g0)){
		dati <- c(dati,rep(NA,7*length(pops)))
		cat(dati,"\n",file=fname,append=TRUE)
		#dat <- rbind(dat, dati)
		next()
		}
	#exclude individuals with too much missing data
	g0 <- g0[,(colSums(is.na(g0))/length(g0[,1]))<.8]
	
	#loop through populations for each window
	for(i in pops){
		p1 <- i
		
		#calc distance matrix between individuals within each population
		cld <- dist(t(g0[,grep(p1,colnames(g0))]),method="manhattan")
		#if any pairwise comparisons lack any comparisons, return NAs
		if(any(is.na(cld))){
			dati <- c(dati,rep(NA,7))
			next()
			}
		#cluster, k=2
		pm1 <- pam(cld,k=2,diss=TRUE)
		
		#organize output
		cl1 <- names(pm1$clustering[pm1$clustering==1])
		cl2 <- names(pm1$clustering[pm1$clustering==2])
		
		len <- c(length(cl1),length(cl2))
	
		m <- c(mean(as.matrix(cld)[cl1,cl1]),mean(as.matrix(cld)[cl2,cl2]))
	
		inter <- mean(as.matrix(cld)[cl1,cl2])
	
		rat <- c(mean(as.matrix(cld)[cl1,cl1])/mean(as.matrix(cld)[cl1,cl2]),mean(as.matrix(cld)[cl2,cl2])/mean(as.matrix(cld)[cl1,cl2]))
	
		ord <- order(rat)
			
		dati <- c(dati,len[ord],m[ord],inter,rat[ord])
	}

	cat(dati,"\n",file=fname,append=TRUE)
	#dat <- rbind(dat, dati)
	
	#scaf start end clus1_length clus2_length mean_clus1 mean_clus2 mean_clus1_2 mclus1/mclus1_2 mclus2/mclus1_2

	}

