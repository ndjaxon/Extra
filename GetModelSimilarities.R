#To rescale similarity score to between 0 and 1
Normalize<-function(x){
	(x - min(x)) / (max(x) - min(x))
}

#Calculate model similarity to the true model
GetModelSimilarityScore<-function(trueModel,currentModel){
	similarity=0
	###COLLAPSE
	#If same collapse history in first time period...
	if(paste(trueModel[1:3],collapse="_") %in% paste(currentModel[1:3],collapse="_")){
		similarity=similarity + 1
	}
	#If there is some collapse...
	if(sum(currentModel[1:3]) > 0){
		similarity=similarity + 1
	}
	#If there is some collapse in second time period...
	if(length(currentModel[4:6][which(!is.na(currentModel[4:6]))]) > 0){	
		if(sum(currentModel[4:6][which(!is.na(currentModel[4:6]))]) > 0){
			similarity=similarity + 1
		}
	}

	###MIGRATION
	#If the same number of migration rates in the first time period...
	if(sum(currentModel[7:12]) == sum(trueModel[7:12])){
		similarity=similarity + 1
	}
	#If the migration matrix matches for the first time period...
	if(paste(trueModel[7:12],collapse="_") %in% paste(currentModel[7:12],collapse="_")){
		similarity=similarity + 1
	}
	#If the same number of migration rates in the second time period...
	if(length(currentModel[13:16][which(!is.na(currentModel[13:16]))]) > 0){
		if(sum(currentModel[13:16][which(!is.na(currentModel[13:16]))]) == 
			sum(trueModel[13:16][which(!is.na(trueModel[13:16]))])){
			similarity=similarity + 1
		}
	}	
	
	return(similarity)
}

#Import phrapl results table
table<-read.table("results.all.txt",header=TRUE,
	stringsAsFactors=FALSE,sep="\t")

#Get similarity scores for each model in the table
similarityVec<-c()
for(i in sequence(nrow(table))){
	if(table$migration[i] == 0){
		trueModel<-table[which(table$models == 21),][1,31:46]
	}else{
		trueModel<-table[which(table$models == 22),][1,31:46]
	}
	currentModel<-table[i,31:46]
	similarityVec<-append(similarityVec,GetModelSimilarityScore(trueModel,currentModel))	
}

#Scale to between zero and one
similarity<-round(Normalize(similarityVec),2)
table<-cbind(table,similarity)



#######Plot similarity across wAIC

#Repete for each nLoci treatment
for(m in sequence(length(unique(table$nLoci)))){



	
	#For now, just do 1 m (locus number) at a time
	thisTable<-table[which(table$nLoci == unique(table$nLoci)[m]),]
	treatments<-paste("t",thisTable$divTime,"_m",thisTable$migration,sep="")
	utreat<-unique(treatments)
	
	#Make 9 panel plot with each migration_divergence treatment combo
	par(mfrow=c(3,3))
	for(k in sequence(length(utreat))){
		currentTreat<-thisTable[which(treatments==utreat[k]),]

		#Make color vector such that the top model is red
		colorVec<-c("red",rep("black",(nrow(currentTreat) - 1)))

		#Plot
		plot(currentTreat$wAIC~currentTreat$similarity,xlab=
			"Similarity to true model",ylab="AIC weight",pch=16,
			axes=FALSE,xlim=c(0,1),ylim=c(0,1),col=colorVec)
		title(main=utreat[k])
		axis(1,at=c(0.0,0.2,0.4,0.6,0.8,1.0),label=
			c("0.0","0.2","0.4","0.6","0.8","1.0"))
		axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),label=
			c("0.0","0.2","0.4","0.6","0.8","1.0"),las=1)
	}




}