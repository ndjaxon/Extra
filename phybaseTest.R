#The goal here is to get an idea of the relationship between lnL and number of taxa (which dictates the tree space within which 
#the probability of seeing a particular tree is calculated). If this relationship is linear (which it seems to be), then if we
#can calculate the slope, then given a particular log likelihood of a subsampled tree (from which we can estimate the 
#y-intercept), we should be able to estimate what the likelihood of the full tree should be, given the likelihood of a subsampled
#tree. 

#To demonstrate this linear relationship, I have simulated a species tree (using ms), simulated a gene tree (30 individuals
#each for 3 populations + outgroup - 91 tips total) under this species tree, and then subsampled this gene tree using the phrapl
#function in different ways. I subsampled 1, 2, 3, 4, 5, 6, 10, 15, and 30 individuals per population (+ outgroup) both with
#replacement (bootstrapping) and without replacement (jackknifing). I subsampled enough replicates to have sampled all individuals
#in the tree once (e.g., subsample of 1 was taken 30 times, 2 was taken 15 times ..., 30 was taken once). Then I imported these
#subsampled trees in to phybase, and for each, calculated the actual likelihood of the species tree, given each gene tree. By,
#summing likelihoods across subsample replicates for a given subsampling regime, I can tell in what sense log likelihoods are
#are additive (e.g., do lnLs from 30 subsamples of 1 add up to the lnL from 1 subsample of 30? The answer is not quite). To
#get a sense of what the relationship actually is, I plotted the number of taxa in a subsample by the mean lnL across replicates 
#(there is a steep negative linear relationship).  By mutliplying the number of taxa by the slope, and adding the intercept, I can
#get the log likelihood. We can use this in phrapl, assuming we can calculate the slope and intercept. This will likely require
#subsampling 2 or three different regimes and caluclating the slope for each model (as the slope slightly changes when the model
#changes.) I in fact recalculated log likelihoods based on a tweaked (false) species tree. The slope is slightly steeper, suggeting
#that this can result in amplified differences in log likelihoods between models when sample sizes are greater, a feature that
#makes intuitive sense. It doesn't look like the actual independence of a subsample makes a difference when it comes to the log
#likelihoods: whether I jackknifed (ensuring the individuals in the overall dataset were only sampled once) or bootstrapped (in
#which individuals could be re-sampled in different subsamples), the relationships were the same.

#Note that trees must be in Nexus format to be imported into phybase. I used phyutility to do this (deleting the taxon name header
#by hand each time, as phybase can't read this header properly).

library(phybase)

#Load tree
speciesTree<-read.tree.string(file="speciesTree.tre") #Read in tree in nexus format
speciesNodeMatrix<-read.tree.nodes(str=speciesTree$tree,name=speciesTree$names)$nodes #get node matrix
speciesRootNode<-rootoftree(speciesNodeMatrix) #declare root node

#Generate single 30 sample (91 total individuals) gene tree under species tree
#Copy and paste this to tre file and use for subsampling
speciesSimTree<-sim.coaltree.sp(rootnode=speciesRootNode,nodematrix=speciesNodeMatrix,
	nspecies=4,seq=c(30,30,30,1),name=speciesTree3$names)$gt	

#Bring in bootstrapped subsampled trees
speciesBootTrees<-list()
for(rep in c(1:6,10,15)){
	inputfile<-paste("/Users/nathanjackson/Desktop/PhyBaseTest/BootstrappedSubsampling_rooted/Subsample",
		rep,"/observed.nex",sep="")
	speciesBootTrees[[length(speciesBootTrees) + 1]]<-read.tree.string(file=inputfile)$tree
}

#Bring in jackknife subsampled trees
speciesJackTrees<-list()
for(rep in c(1:6,10,15)){
	inputfile<-paste("/Users/nathanjackson/Desktop/PhyBaseTest/JackknifedSubsampling_rooted/Subsample",
		rep,"/observed.nex",sep="")
	speciesJackTrees[[length(speciesJackTrees) + 1]]<-read.tree.string(file=inputfile)$tree
}

#Bring in single 30 sample gene tree and add this to the bootstrapped and jackniffed lists
geneTree_30<-read.tree.string(file="JackknifedSubsampling_rooted/observed.nex")
speciesBootTrees[[length(speciesBootTrees) + 1]]<-geneTree_30$tree
speciesJackTrees[[length(speciesJackTrees) + 1]]<-geneTree_30$tree


#Make assignment matrices for the different subsample sizes
#Make empty matrices
speciesSimAssign<-list()
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,4)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,7)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,10)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,13)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,16)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,19)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,31)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,46)
speciesSimAssign[[length(speciesSimAssign) + 1]]<-matrix(0,4,91)

#Fill in matrices to match dataset
for(rep in 1:length(speciesSimAssign)){
	startCol<-1
	endCol<-(ncol(speciesSimAssign[[rep]]) - 1) / 3
	for(rep2 in 1:(nrow(speciesSimAssign[[1]]) - 1)){
		speciesSimAssign[[rep]][rep2,startCol:endCol]<-1
		startCol<-startCol + ((ncol(speciesSimAssign[[rep]]) - 1) / 3)
		endCol<-endCol + ((ncol(speciesSimAssign[[rep]]) - 1) / 3)
	}
	speciesSimAssign[[rep]][4,startCol]<-1
}

#Create individual names for the subsampled datasets (doing this because the "species.name" function gets them out of order)
taxonSubNames<-list()
for(rep in 1:length(speciesSimAssign)){
	taxonSubNames[[length(taxonSubNames) + 1]]<-as.character(c(1:ncol(speciesSimAssign[[rep]])))
}

#Create species names for the species tree
speciesNames<-c("1","2","3","4")


#Calculate lnL of species tree give each gene tree subsample (which includes the whole 30 sample tree)
#bootstraps
lnL_boot<-list()
for(rep in 1:length(speciesBootTrees)){
	lnL_boot_current<-array()
	for(rep2 in 1:length(speciesBootTrees[[rep]])){
		lnL_boot_current<-c(lnL_boot_current,loglikeSP(speciesBootTrees[[rep]][rep2],speciesTree$tree,taxonSubNames[[rep]],speciesNames,
			speciesSimAssign[[rep]]))
	}
	lnL_boot[[length(lnL_boot) + 1]]<-lnL_boot_current[!is.na(lnL_boot_current)]
}

#jackknifes
lnL_jack<-list()
for(rep in 1:length(speciesJackTrees)){
	lnL_jack_current<-array()
	for(rep2 in 1:length(speciesJackTrees[[rep]])){
		lnL_jack_current<-c(lnL_jack_current,loglikeSP(speciesJackTrees[[rep]][rep2],speciesTree$tree,taxonSubNames[[rep]],speciesNames,
			speciesSimAssign[[rep]]))
	}
	lnL_jack[[length(lnL_jack) + 1]]<-lnL_jack_current[!is.na(lnL_jack_current)]
}

#Get sum of each set of subsamples (are they equal?)
#bootstraps
lnL_boot_sum<-array()
for(rep in 1:length(lnL_boot)){
	lnL_boot_sum<-c(lnL_boot_sum,sum(lnL_boot[[rep]]))
}
lnL_boot_sum<-lnL_boot_sum[!is.na(lnL_boot_sum)]

#Get sum of each set of subsamples (are they equal?)
#jackknifes
lnL_jack_sum<-array()
for(rep in 1:length(lnL_jack)){
	lnL_jack_sum<-c(lnL_jack_sum,sum(lnL_jack[[rep]]))
}
lnL_jack_sum<-lnL_jack_sum[!is.na(lnL_jack_sum)]

#Get mean of each set of subsamples
#bootstraps
lnL_boot_mean<-array()
for(rep in 1:length(lnL_boot)){
	lnL_boot_mean<-c(lnL_boot_mean,mean(lnL_boot[[rep]]))
}
lnL_boot_mean<-lnL_boot_mean[!is.na(lnL_boot_mean)]

#Get mean of each set of subsamples
#jackknifes
lnL_jack_mean<-array()
for(rep in 1:length(lnL_jack)){
	lnL_jack_mean<-c(lnL_jack_mean,mean(lnL_jack[[rep]]))
}
lnL_jack_mean<-lnL_jack_mean[!is.na(lnL_jack_mean)]





########Now also create a tweaked species tree, and recalculate likelihoods#####
#Load tree
speciesTree<-read.tree.string(file="speciesTree_wrong.tre") #Read in tree in nexus format
speciesNodeMatrix<-read.tree.nodes(str=speciesTree$tree,name=speciesTree$names)$nodes #get node matrix
speciesRootNode<-rootoftree(speciesNodeMatrix) #declare root node

#Calculate lnL of species tree give each gene tree subsample (which includes the whole 30 sample tree)
#bootstraps
lnL_boot_wrong<-list()
for(rep in 1:length(speciesBootTrees)){
	lnL_boot_wrong_current<-array()
	for(rep2 in 1:length(speciesBootTrees[[rep]])){
		lnL_boot_wrong_current<-c(lnL_boot_wrong_current,loglikeSP(speciesBootTrees[[rep]][rep2],speciesTree$tree,taxonSubNames[[rep]],speciesNames,
			speciesSimAssign[[rep]]))
	}
	lnL_boot_wrong[[length(lnL_boot_wrong) + 1]]<-lnL_boot_wrong_current[!is.na(lnL_boot_wrong_current)]
}

#jackknifes
lnL_jack_wrong<-list()
for(rep in 1:length(speciesJackTrees)){
	lnL_jack_wrong_current<-array()
	for(rep2 in 1:length(speciesJackTrees[[rep]])){
		lnL_jack_wrong_current<-c(lnL_jack_wrong_current,loglikeSP(speciesJackTrees[[rep]][rep2],speciesTree$tree,taxonSubNames[[rep]],speciesNames,
			speciesSimAssign[[rep]]))
	}
	lnL_jack_wrong[[length(lnL_jack_wrong) + 1]]<-lnL_jack_wrong_current[!is.na(lnL_jack_wrong_current)]
}

#Get sum of each set of subsamples (are they equal?)
#bootstraps
lnL_boot_wrong_sum<-array()
for(rep in 1:length(lnL_boot_wrong)){
	lnL_boot_wrong_sum<-c(lnL_boot_wrong_sum,sum(lnL_boot_wrong[[rep]]))
}
lnL_boot_wrong_sum<-lnL_boot_wrong_sum[!is.na(lnL_boot_wrong_sum)]

#Get sum of each set of subsamples (are they equal?)
#jackknifes
lnL_jack_wrong_sum<-array()
for(rep in 1:length(lnL_jack_wrong)){
	lnL_jack_wrong_sum<-c(lnL_jack_wrong_sum,sum(lnL_jack_wrong[[rep]]))
}
lnL_jack_wrong_sum<-lnL_jack_wrong_sum[!is.na(lnL_jack_wrong_sum)]

#Get mean of each set of subsamples
#bootstraps
lnL_boot_wrong_mean<-array()
for(rep in 1:length(lnL_boot_wrong)){
	lnL_boot_wrong_mean<-c(lnL_boot_wrong_mean,mean(lnL_boot_wrong[[rep]]))
}
lnL_boot_wrong_mean<-lnL_boot_wrong_mean[!is.na(lnL_boot_wrong_mean)]

#Get mean of each set of subsamples
#jackknifes
lnL_jack_wrong_mean<-array()
for(rep in 1:length(lnL_jack_wrong)){
	lnL_jack_wrong_mean<-c(lnL_jack_wrong_mean,mean(lnL_jack_wrong[[rep]]))
}
lnL_jack_wrong_mean<-lnL_jack_wrong_mean[!is.na(lnL_jack_wrong_mean)]




###########PLOTS#########
#Plot lnL against the number of taxa in a subsample (or possible trees for given sample size, which gives same answer)
ntax=c(3*c(1:6,10,15,30)+1) #sample size vector
loghow<-(-log(sapply(ntax, howmanytrees))) #log-likelihood of a tree given the total tree space (this gives same relationship
#as the number of taxa in a subsample)
lm(lnL_boot_mean ~ ntax) #slope for bootstraped
lm(lnL_jack_mean ~ ntax) #slope for jackknifed

#plot bootstrapped subsamples versus jackknifed subsamples
plot(lnL_boot_mean ~ ntax,type='n')
points(lnL_boot_mean ~ ntax,col="red2")
points(lnL_jack_mean ~ ntax,col="green4")
abline(lm(lnL_boot_mean ~ ntax),col="red2")
abline(lm(lnL_jack_mean ~ ntax),col="green4")

#plot lnLs from right species tree versus the wrong species tree
plot(c(lnL_boot_mean[1:8],lnL_boot_wrong_mean[9]) ~ ntax,type='n',ylab="mean lnL",las=1)
points(lnL_boot_mean ~ ntax,col="red2")
points(lnL_boot_wrong_mean ~ ntax,col="green4")
abline(lm(lnL_boot_mean ~ ntax),col="red2")
abline(lm(lnL_boot_wrong_mean ~ ntax),col="green4")



####Additional Plots
colors<-c("red","darkorange2","goldenrod","forestgreen","darkturquoise",
	"darkblue","blueviolet") #red,orange,yellow,green,aqua,blue,purple

#plot how slopes change based on 2 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) -2)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+1)] ~ ntax[which.lnL:(which.lnL+1)]),col=colors[i])
	which.lnL<-which.lnL + 1
}	

#plot how slopes change based on 3 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 3)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+2)] ~ ntax[which.lnL:(which.lnL+2)]),col=colors[i])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on 4 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 4)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+3)] ~ ntax[which.lnL:(which.lnL+3)]),col=colors[i])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on sampling 2, with spread of 2
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 3)){
	abline(lm(lnL_boot_mean[c(which.lnL,(which.lnL+2))] ~ ntax[c(which.lnL,(which.lnL+2))]),col=colors[i])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on sampling 2, with spread of 3
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 4)){
	abline(lm(lnL_boot_mean[c(which.lnL,(which.lnL+3))] ~ ntax[c(which.lnL,(which.lnL+3))]),col=colors[i])
	which.lnL<-which.lnL + 1
}



#####Repeat the above plots, but including slopes for an incorrect model
#plot how slopes change based on 2 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
points(lnL_boot_wrong_mean ~ ntax,col=colors[6])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) -2 )){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+1)] ~ ntax[which.lnL:(which.lnL+1)]),col=colors[1])
	abline(lm(lnL_boot_wrong_mean[which.lnL:(which.lnL+1)] ~ ntax[which.lnL:(which.lnL+1)]),col=colors[6])
	which.lnL<-which.lnL + 1
}	

#plot how slopes change based on 3 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
points(lnL_boot_wrong_mean ~ ntax,col=colors[6])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 3)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+2)] ~ ntax[which.lnL:(which.lnL+2)]),col=colors[1])
	abline(lm(lnL_boot_wrong_mean[which.lnL:(which.lnL+2)] ~ ntax[which.lnL:(which.lnL+2)]),col=colors[6])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on 4 step intervals
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
points(lnL_boot_wrong_mean ~ ntax,col=colors[6])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 4)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+3)] ~ ntax[which.lnL:(which.lnL+3)]),col=colors[1])
	abline(lm(lnL_boot_wrong_mean[which.lnL:(which.lnL+3)] ~ ntax[which.lnL:(which.lnL+3)]),col=colors[6])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on sampling 2, with spread of 2
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
points(lnL_boot_wrong_mean ~ ntax,col=colors[6])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 3)){
	abline(lm(lnL_boot_mean[c(which.lnL,(which.lnL+2))] ~ ntax[c(which.lnL,(which.lnL+2))]),col=colors[1])
	abline(lm(lnL_boot_wrong_mean[c(which.lnL,(which.lnL+2))] ~ ntax[c(which.lnL,(which.lnL+2))]),col=colors[6])
	which.lnL<-which.lnL + 1
}

#plot how slopes change based on sampling 2, with spread of 3
plot(lnL_boot_mean ~ ntax,type='n',ylab="mean lnL")
points(lnL_boot_mean ~ ntax,col=colors[1])
points(lnL_boot_wrong_mean ~ ntax,col=colors[6])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) - 4)){
	abline(lm(lnL_boot_mean[c(which.lnL,(which.lnL+3))] ~ ntax[c(which.lnL,(which.lnL+3))]),col=colors[1])
	abline(lm(lnL_boot_wrong_mean[c(which.lnL,(which.lnL+3))] ~ ntax[c(which.lnL,(which.lnL+3))]),col=colors[6])
	which.lnL<-which.lnL + 1
}



####One idea is that instead of calculating the slope based on the mean lnL for a set of points, we should
#calculate the slope for each subsample, and then use these to build a condidence interval, such that we can determine
#how much confidence in the average slope. However, looking at these plots, this will only give us a confidence interval
#around the slope for these sets of points, it doesn't help us get a better sense of how well we can extrapolate the
#slope to larger samples sizes.
colors<-c("red","darkorange2","goldenrod","forestgreen","darkturquoise",
	"darkblue","blueviolet") #red,orange,yellow,green,aqua,blue,purple

#plot how slopes change based on 2 step intervals
plot(lnL_boot_mean ~ ntax,type='n')
points(lnL_boot_mean ~ ntax,col=colors[1])
which.lnL<-1
for(i in 1:(length(lnL_boot_mean) -2)){
	abline(lm(lnL_boot_mean[which.lnL:(which.lnL+1)] ~ ntax[which.lnL:(which.lnL+1)]),col=colors[i])
	which.lnL<-which.lnL + 1
}

#Here, I am overlaying slopes based on each subsample separately (for the outlier slope (green) above)
#These slopes are in black. Their range is, not suprisingly, even greater. Also, note that 
for(j in 1:6){
thispair<-c(lnL_boot[[4]][j],lnL_boot[[5]][j])
	abline(lm(thispair ~ ntax[4:5]),col="black")
}


