#The purpose of this script is to extract migration rate and tau values from 
#the Hey and Pinho 2012 meta-analysis table. I organize these according
#to taxanomic rank and species/population status
#I can print these values on the delimitation index heat map

tab<-read.table("Hey&Pinho_Evolution_2012_Table.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

#Make subsets of data (excluding subspecies)
pops<-tab[which(tab$rank=="populations"),]
species<-tab[which(tab$rank=="species"),]
#birds
pops_bird<-pops[which(pops$group=="bird"),]
species_bird<-species[which(species$group=="bird"),]
#mammals
pops_mammal<-pops[which(pops$group=="mammal"),]
species_mammal<-species[which(species$group=="mammal"),]
#insects
pops_insect<-pops[which(pops$group=="insect"),]
species_insect<-species[which(species$group=="insect"),]
#plants
pops_plant<-pops[which(pops$group=="plant"),]
species_plant<-species[which(species$group=="plant"),]
#fish
pops_fish<-pops[which(pops$group=="fish"),]
species_fish<-species[which(species$group=="fish"),]


#Build new dataframe to place median values
med.values<-data.frame("group"=c("all","bird","mammal","insect",
	"plant","fish"),"n_pops"=c((nrow(pops)*2),(nrow(pops_bird)*2),(nrow(pops_mammal)*2),
	(nrow(pops_insect)*2),(nrow(pops_plant)*2),(nrow(pops_fish)*2)),"n_species"=
	c((nrow(species)*2),(nrow(species_bird)*2),(nrow(species_mammal)*2),(nrow(species_insect)*2),
	(nrow(species_plant)*2),(nrow(species_fish)*2)),"tau_pops_med"=NA,"tau_species_med"=NA,
	"mig_pops_med"=NA,"mig_species_med"=NA)

#Calculate median values
#For tau, just add the two taus to get in 4N
#for 2Nm values, take the average between the two and then divide by 2 to get 4Nm
pops_all_med<-sapply(pops[13:16],median)
species_all_med<-sapply(species[13:16],median)
med.values$tau_pops_med[1]<-pops_all_med[2] + pops_all_med[4]
med.values$tau_species_med[1]<-species_all_med[2] + species_all_med[4]
med.values$mig_pops_med[1]<-mean(c(pops_all_med[1],pops_all_med[3])) / 2
med.values$mig_species_med[1]<-mean(c(species_all_med[1],species_all_med[3])) / 2

pops_bird_med<-sapply(pops_bird[13:16],median)
species_bird_med<-sapply(species_bird[13:16],median)
med.values$tau_pops_med[2]<-pops_bird_med[2] + pops_bird_med[4]
med.values$tau_species_med[2]<-species_bird_med[2] + species_bird_med[4]
med.values$mig_pops_med[2]<-mean(c(pops_bird_med[1],pops_bird_med[3])) / 2
med.values$mig_species_med[2]<-mean(c(species_bird_med[1],species_bird_med[3])) / 2

pops_mammal_med<-sapply(pops_mammal[13:16],median)
species_mammal_med<-sapply(species_mammal[13:16],median)
med.values$tau_pops_med[3]<-pops_mammal_med[2] + pops_mammal_med[4]
med.values$tau_species_med[3]<-species_mammal_med[2] + species_mammal_med[4]
med.values$mig_pops_med[3]<-mean(c(pops_mammal_med[1],pops_mammal_med[3])) / 2
med.values$mig_species_med[3]<-mean(c(species_mammal_med[1],species_mammal_med[3])) / 2

pops_insect_med<-sapply(pops_insect[13:16],median)
species_insect_med<-sapply(species_insect[13:16],median)
med.values$tau_pops_med[4]<-pops_insect_med[2] + pops_insect_med[4]
med.values$tau_species_med[4]<-species_insect_med[2] + species_insect_med[4]
med.values$mig_pops_med[4]<-mean(c(pops_insect_med[1],pops_insect_med[3])) / 2
med.values$mig_species_med[4]<-mean(c(species_insect_med[1],species_insect_med[3])) / 2

pops_plant_med<-sapply(pops_plant[13:16],median)
species_plant_med<-sapply(species_plant[13:16],median)
med.values$tau_pops_med[5]<-pops_plant_med[2] + pops_plant_med[4]
med.values$tau_species_med[5]<-species_plant_med[2] + species_plant_med[4]
med.values$mig_pops_med[5]<-mean(c(pops_plant_med[1],pops_plant_med[3])) / 2
med.values$mig_species_med[5]<-mean(c(species_plant_med[1],species_plant_med[3])) / 2

pops_fish_med<-sapply(pops_fish[13:16],median)
species_fish_med<-sapply(species_fish[13:16],median)
med.values$tau_pops_med[6]<-pops_fish_med[2] + pops_fish_med[4]
med.values$tau_species_med[6]<-species_fish_med[2] + species_fish_med[4]
med.values$mig_pops_med[6]<-mean(c(pops_fish_med[1],pops_fish_med[3])) / 2
med.values$mig_species_med[6]<-mean(c(species_fish_med[1],species_fish_med[3])) / 2

#Print out summary table
write.table(Hey&Pinho_med.values,file="median.values.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)