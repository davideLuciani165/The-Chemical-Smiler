

stringColl<- function(matr)
	{cumRes<- c()
	parent1<- c("", "")
	parent2<- c("[", "]")
	parent3<- c("<", ">")
	parent4<- c("", "")
	parLs<- list(parent1)
	if(ncol(matr)>1)
		{for (h in 2: ncol(matr))
			{#h<- 1; h<- h+1
			if(h==2)
				{parLs[[length(parLs)+1]]<-parent2}
			if(h==3)
				{parLs[[length(parLs)+1]]<-parent3}
			if(h>3)
				{parLs[[length(parLs)+1]]<-parent1}
			}
		}
	for (i in 1:nrow(matr))
		{#i<- 1; i<- i+1
		cumString<- c()
		for (j in 1: ncol(matr))
			{#j<- 1; j<- j+1
			cumString<- paste(cumString, parLs[[j]][1], matr[i,j],  parLs[[j]][2], sep="", collapse="")
			}
		cumRes<- c(cumRes, cumString)
		}
	return(cumRes)
	}

#matr<- matr_n
stringXColl<- function(matr)
	{cumRes<- c()
	parent1<- c("", "")
	parent2<- c("[", "]")
	parent3<- c("<", ">")
	parent4<- c("", "")
	parLs<- list(parent1)
	if(ncol(matr)>1)
		{for (h in 2: ncol(matr))
			{#h<- 1; h<- h+1
			if(h==2)
				{parLs[[length(parLs)+1]]<-parent2}
			if(h==3)
				{parLs[[length(parLs)+1]]<-parent3}
			if(h>3)
				{parLs[[length(parLs)+1]]<-parent1}
			}
		}
	for (i in 1:nrow(matr))
		{#i<- 1; i<- i+1
		cumString<- c()
		for (j in 1: ncol(matr))
			{#j<- 1; j<- j+1
			cumString<- paste(cumString, ifelse(is.na(matr[i,j]), "", parLs[[j]][1]), ifelse(is.na(matr[i,j]), "", matr[i,j]),  ifelse(is.na(matr[i,j]), "", parLs[[j]][2]), sep="", collapse="")
			}
		cumRes<- c(cumRes, cumString)
		}
	return(cumRes)
	}


#options(warn = 2)
## m – item in the first section of the Identifier, but not in the isotopic segment
#emptyToo<- TRUE
#expandAbbrev<- TRUE
decompInChI<- function(inchiData, emptyToo, expandAbbrev)
	{fullStruct<- list()
	if(!is.matrix(inchiData))
		{inchiData<- as.matrix(inchiData)
		colnames(inchiData)<- c("CAS_NO", "InChI")
		}
	for (a in 1: nrow(inchiData))
		{#a<-1 3 1; a<- a+1

		inChIStr<- rep(list(list("")), 2)
		PinChI<- c(rep(list(""),2))
		names(PinChI)<- c("disconnected", "r_iconnected")

		#1) Find the first slash. The slash is preceded by the version and followed by
		#a string (call it S) that contains all other layers of the identifier.
		currInChI<- inchiData[a, "InChI"]
		SinChI <- paste(unlist(strsplit(currInChI, "/", fixed=TRUE))[-1], collapse="/")

		#2) Search for “/r” in S.
		#(P[1] represents the whole identifier or an identifier of a disconnected
		#structure; P[2] if not empty represents an identifier of a “reconnected”
		#structure”.)
		app_P<- unlist(strsplit(SinChI, "/r", fixed=TRUE)) 
		testR<- length(app_P)>1
		if(testR)
			{#If “/r” is found then copy preceding “/r” substring to P[1]
			#and the following “/r” string to P[2] 
			PinChI[[1]]<- app_P[1]
			PinChI[[2]]<- app_P[2]		
			}
		if(!testR)
			{#else copy S to P[1]
			PinChI[[1]]<- SinChI
			}
		#3)Search for “/f” in each non-empty P[i].
		for (i in 1: length(PinChI))
			{#i<- 1; i<- i+1
			QinChI<-  c(rep(list(""),2))
			names(QinChI)<- c("main", "f_ixed")
			if(nchar(PinChI[[i]])>0)
				{app_Q<- unlist(strsplit(PinChI[[i]], "/f", fixed=TRUE))
				testF<- length(app_Q)>1
				if(testF)
					{#If “/f” was found 
					#then copy the preceding string to Q[i][1] 
					QinChI[[1]]<- app_Q[1]
					#and the following string to Q[i][2]
					QinChI[[2]]<- app_Q[2]	
					#QinChI[[2]]<- abbrevSubst(QinChI)
			
					}
				if(!testF)
					{#else copy P[i] to Q[i][1]
					QinChI[[1]]<- PinChI[[i]]
					#(Q[i][1] represents the Main layer; Q[i][2] represents fixed-H layer)
					}
				#4) Search for “/i" in each non-empty Q[i][j]. 
				for (j in 1: length(QinChI))
					{#j<- 1; j<- j+1
					if(nchar(QinChI[[j]])>0)
						{RinChI<-  c(rep(list(""),2))
						names(RinChI)<- c("nonIsotopic", "i_sotopic")
						app_I<- unlist(strsplit(QinChI[[j]], "/i", fixed=TRUE))
						testI<- length(app_I)>1
						if(testI)
							{#If “/i" was found 
							 #then copy the preceding string into R[i][j][1] 
							RinChI[[1]]<- app_I[1]
							#and the following string into R[i][j][2] 
							RinChI[[2]]<- app_I[2]
							}
						if(!testI)
							{#else copy Q[i][j] to R[i][j][1]
							RinChI[[1]]<- QinChI[[j]]
							}
						inChIStr[[i]][[j]]<- RinChI
						names(inChIStr[[i]])<- names(QinChI)[seq(1,j)]
						}
					}
				}
			names(inChIStr)<- names(PinChI)[seq(1,i)]
			}
		if(!emptyToo)
			{if(!is.null(inChIStr$disconnected$main$i_sotopic))
				{if(inChIStr$disconnected$main$i_sotopic=="")
					{inChIStr$disconnected$main$i_sotopic<- NULL}
				}
			if(!is.null(inChIStr$disconnected$f_ixed$i_sotopic))
				{if(inChIStr$disconnected$f_ixed$i_sotopic=="")
					{inChIStr$disconnected$f_ixed$i_sotopic<- NULL}
				}
			if(!is.null(inChIStr$disconnected$f_ixed$i_sotopic))
				{if(inChIStr$disconnected$f_ixed$i_sotopic=="")
					{inChIStr$disconnected$f_ixed$i_sotopic<- NULL}
				}
			if(length(inChIStr$riconnected[[1]])>1)
				{for (m in 1:length(inChIStr$r_iconnected[[1]]))
					{if(!is.null(inChIStr$r_iconnected[[1]][[m]]))
						{if(inChIStr$r_iconnected[[1]][[m]]=="")
							{inChIStr$r_iconnected[[1]][[m]]<- NULL}
						}
					}
				}
			if(length(inChIStr$r_iconnected[[1]])==1)
				{if(!is.null(inChIStr$r_iconnected[[1]]))
					{if(length(inChIStr$r_iconnected[[1]])==0 | inChIStr$r_iconnected[[1]]=="")
						{inChIStr$r_iconnected[[1]]<- NULL}
					}			
				}
			if(is.null(inChIStr$r_iconnected)|length(inChIStr$r_iconnected)==0)
				{inChIStr$r_iconnected <- NULL		
				}
			if(!is.null(inChIStr$r_iconnected$main$i_sotopic))
				{if(inChIStr$r_iconnected$main$i_sotopic=="")
					{inChIStr$r_iconnected$main$i_sotopic<- NULL}
				}
			}
		if(expandAbbrev)
			{fullStruct[[1+length(fullStruct)]]<-  abbrevSubst(inChIStr) }
		if(!expandAbbrev)
			{fullStruct[[1+length(fullStruct)]]<-  inChIStr }		
		}
 
	names(fullStruct)<- inchiData[, "CAS_NO"]
	return(fullStruct)
	}

#extraData<- manInchiData[notMissInchi, c("ID", "idS","compNum","ElCount","MolWeight","Charge","saltClass","complexClass","metalClass","identCateg","eraseComp","flagOrg","WarningFromNormal","normFlag")]
#colnames(extraData)

#selData<- auxData[keepInx,]
keepUniqueComp<- function(selData, vars=c("ID", "identCateg"))
	{if(!is.matrix(selData))
		{selData <- t(as.matrix(selData))}
	uniID<- unique(selData[,vars[1]])
	keepRec<- c()
	for (z in 1: length(uniID))
		{#z<- 1; z- z+1
		selID<- which(selData[,vars[1]]==uniID[z])
		keepRec<- c(keepRec, selID[!duplicated(selData[selID,vars[2]])])
		}
	return(keepRec)
	}

# selData<- extAuxData[O_keepInx,]
#vars=c("ID", "MolWeight")

keepLargestComp<- function(selData, vars=c("ID", "MolWeight"), orgType)
	{if(!is.matrix(selData))
		{selData <- t(as.matrix(selData))}
	uniID<- unique(selData[,vars[1]])
	keepRec<- c()
	for (z in 1: length(uniID))
		{#z<- 1; z- z+1
		if(orgType=="all")
			{selID<- which(selData[,vars[1]]==uniID[z])}
		if(orgType!="all")
			{ selID<- which(selData[,vars[1]]==uniID[z] & selData[,"flagOrg"]==1)}
		if(length(selID)==0)
			{selID<- which(selData[,vars[1]]==uniID[z])}
			
		molWeights<-as.numeric(gsub(" ", "", selData[selID,vars[2]]))
			
		keepRec<- c(keepRec, selID[which(molWeights==max(molWeights)[1])])
		#keepRec<- c(keepRec, selID[which(selData[selID,vars[2]]==max(selData[selID,vars[2]])[1])])
		}
	return(keepRec)
	}

#La funzione keepLargestOrgComp seleziona la componente piu grande solo fra quelle organiche
keepLargestOrgComp<- function(selData, vars=c("ID", "MolWeight"))
	{if(!is.matrix(selData))
		{selData <- t(as.matrix(selData))}
	uniID<- unique(selData[,vars[1]])
	keepRec<- c()
	for (z in 1: length(uniID))
		{#z<- 1; z- z+1
		selID<- which(selData[,vars[1]]==uniID[z] & selData[,"flagOrg"]==1)
		if(length(selID)==0)
			{selID<- which(selData[,vars[1]]==uniID[z])}
		keepRec<- c(keepRec, selID[which(selData[selID,vars[2]]==max(selData[selID,vars[2]])[1])])
		}
	return(keepRec)
	}

#xData<- auxData
#compVar=""
#xData[,"flagOrg"]
#which(xData[,"flagOrg"]== -1)
restructComp<- function(xData, vars=c("ID", "metalClass", "flagOrg", "WarningFromNormal"), compVar="compNum")
	{if(!is.matrix(xData))
		{xData <- t(as.matrix(xData))}
	uniID<- unique(xData[,vars[1]])
	keepRec<- c()
	xData[, "WarningFromNormal"]<- as.character(xData[, "WarningFromNormal"])
	for (z in 1: length(uniID))
		{#z<- 1; z- z+1
		selID<- which(xData[,vars[1]]==uniID[z])
		for (v in 2: length(vars))
			{#v<- 1; v<- v+1
			currVar<- vars[v]
			if(currVar=="metalClass")
				{metGroup<- as.integer(xData[selID, "metalClass"])
				testMet<- sum(metGroup)>0
				if(testMet)
					{xData[selID, "metalClass"]<- 1}
				}
			if(currVar=="flagOrg")
				{organGroup<- xData[selID, "flagOrg"]
				orgComp<- which(organGroup==1)
				inorgComp<- which(organGroup==0 & xData[selID, "metalClass"]==0)
				unkonwnComp<- which(organGroup== -1)
				testUnknown<- (length(unkonwnComp)>0 & ( length(orgComp)==0 | length(inorgComp)==0 ))
				testMixed<- (length(orgComp)>0 & length(inorgComp)>0)
				testOnlyOrgan<- (!testUnknown & !testMixed & length(orgComp)>0)
				testOnlyInorgan<- (!testUnknown & !testMixed & length(inorgComp)>0)
				if(testMixed)
					{xData[selID, "flagOrg"]<- 2}
				if(testOnlyOrgan)
					{xData[selID, "flagOrg"]<- 1}
				if(testOnlyInorgan)
					{xData[selID, "flagOrg"]<- 0}
				if(testUnknown)
					{xData[selID, "flagOrg"]<- -1}
				}
			if(currVar=="WarningFromNormal")
				{warnGroup<- xData[selID, "WarningFromNormal"]
				notRepeatWarn<- unique(warnGroup)
				notRepeatWarn<- notRepeatWarn[which(notRepeatWarn!="")]
				xData[selID, "WarningFromNormal"]<- rep(paste(notRepeatWarn, collapse="/"), length(selID))
				}
			}
		}
	if(compVar=="compNum")
		{firCompInx<- which(xData[,"compNum"]==1)
		restructData<- xData[firCompInx,]}
	else {restructData<- xData}
	if(!is.data.frame(restructData))
	{restructData <- as.data.frame(restructData)}
	return(restructData)
	}

#options(warn=2)
#colnames(restructData)
#restructData$saltClass[500:855]
#restructData$metalClass[500:855]
#table(restructData$metalClass)
#restructData$flagOrg[500:855]
#table(restructData$flagOrg)
#inchiData<- finchiData
##>> which(is.na(inchiData[,"ID"]))
#colnames(extraData)
 
rewriteInChI<- function(inchiData, auxData=extraData, #missData=missDat,
				tautom=c(0, 1, 2)[1], 
				stereoChem=FALSE, 
				neutralize=c("no", "yes")[2], 
				atomCounterIones=c("no", "yes")[1],
				largestComponent=c("no", "yes")[1],
				orgType=c("all", "onlyOrganic", "alsoOrgInorg")[3],
				noMetallic=c("no", "yes")[2],
				reduce=FALSE, standard=FALSE, extendCAS=TRUE, bestMatch=TRUE,
				mainDir="E:/Data/Environment/projects/StandardSmiles/R") 
	{repackData<- c()
	#keepInx<- 0

	if(!is.null(auxData))
		{numNotIdS<- which(colnames(inchiData)!="idS")
		originalFeature<- rep("", nrow(auxData))
		restrFullDat<- restructComp(auxData, vars=c("ID", "metalClass", "flagOrg", "WarningFromNormal","InChI_SingComp", "Charge"), compVar="")
		for (b in 1:length(originalFeature))
			{#b<- 1; b<- b+1
			feaSalt<- ifelse(restrFullDat$saltClass[b]==1, "salt", "neutral")
			if(restrFullDat$metalClass[b]==1)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==1, "organometallic", ifelse(restrFullDat$flagOrg[b]==-1,"unkonwn-metallic", "metallic"))}
			if(restrFullDat$metalClass[b]==0)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==2, "organic-inorganic", ifelse(restrFullDat$flagOrg[b]==1, "organic",  ifelse(restrFullDat$flagOrg[b]==-1,"unknown", "inorganic")))}
			originalFeature[b]<- paste(feaSalt, feaOrgMetal, sep=", ")
			}
		extAuxData<- cbind(auxData, originalFeature)
		colnames(extAuxData)<- c(colnames(auxData), "originalFeature")
		if(neutralize=="yes")
			{keep1Inx <- which(extAuxData[,"normFlag"]!=0)
			## O_keep1Inx <- which(extAuxData[,"normFlag"]==0)
			## appkeep1Inx<- which(extAuxData[keep1Inx, "ID"] %in% extAuxData[O_keep1Inx, "ID"])
			## O_keep1Inx<- O_keep1Inx[appkeep1Inx]

			keep2Inx   <- which(extAuxData[keep1Inx,"eraseComp"]==0|extAuxData[keep1Inx,"saltClass"]!=0)
			## O_keep2Inx   <- which(extAuxData[O_keep1Inx,"eraseComp"]==0|extAuxData[O_keep1Inx,"saltClass"]!=0)
			if(orgType=="onlyOrganic")
				{keep3Inx<- which(extAuxData[keep1Inx,"flagOrg"]==1)
				## O_keep3Inx<- which(extAuxData[O_keep1Inx,"flagOrg"]==1)
				}
			if(orgType=="alsoOrgInorg")
				{exclFromKeep3<- which(restrFullDat[keep1Inx,"flagOrg"]==0)
				## O_exclFromKeep3<- which(restrFullDat[O_keep1Inx,"flagOrg"]==0)
				keep3Inx<- setdiff(seq(1, nrow(extAuxData[keep1Inx,])), exclFromKeep3)  
				## O_keep3Inx<- setdiff(seq(1, nrow(extAuxData[O_keep1Inx,])), O_exclFromKeep3) 
				}
			if(orgType=="all")
				{keep3Inx<-  seq(1, nrow(extAuxData[keep1Inx,]))
				## O_keep3Inx<-  seq(1, nrow(extAuxData[O_keep1Inx,]))
				}
			if(noMetallic=="yes")
				{keep4Inx<- which(extAuxData[keep1Inx,"metalClass"]!=1)
				## O_keep4Inx<- which(extAuxData[O_keep1Inx,"metalClass"]!=1)
				}
			if(noMetallic=="no")
				{keep4Inx<-  seq(1, nrow(extAuxData[keep1Inx,]))
				## O_keep4Inx<-  seq(1, nrow(extAuxData[O_keep1Inx,]))
				}
			if(atomCounterIones=="yes") #saltClass
				{keep5Inx<- which(extAuxData[keep1Inx ,"complexClass"]!=-1|extAuxData[keep1Inx ,"saltClass"]!=0)
				## O_keep5Inx<- which(extAuxData[O_keep1Inx ,"complexClass"]!=-1|extAuxData[O_keep1Inx ,"saltClass"]!=0)
				}
			if(atomCounterIones=="no")
				{keep5Inx<-  seq(1, nrow(extAuxData[keep1Inx,]))
				## O_keep5Inx<-  seq(1, nrow(extAuxData[O_keep1Inx,]))
				}
			keepInx<- intersect(intersect(intersect(intersect(keep1Inx, keep2Inx), keep3Inx), keep4Inx), keep5Inx)
			## O_keepInx<- intersect(intersect(intersect(intersect(O_keep1Inx, O_keep2Inx), O_keep3Inx), O_keep4Inx), O_keep5Inx)

			keepInx<-   intersect(intersect(intersect(keep2Inx, keep3Inx), keep4Inx), keep5Inx)
			## O_keepInx<- intersect(intersect(intersect(O_keep2Inx, O_keep3Inx), O_keep4Inx), O_keep5Inx)

#which(extAuxData$InChI_SingComp[O_keepInx]=="")
			if(largestComponent=="yes") #saltClass
				{keepLargInx<- keepLargestComp(extAuxData[keepInx,], vars=c("ID", "MolWeight"), orgType)
				 ## O_keepLargInx<- keepLargestComp(extAuxData[O_keepInx,], vars=c("ID", "MolWeight"), orgType)
				 }
#which(extAuxData$InChI_SingComp[O_keepLargInx]=="")
			if(largestComponent=="no")
				{appextAuxData<-extAuxData[keepInx,]
				if(!is.matrix(appextAuxData))
					{appextAuxData<-t(as.matrix(appextAuxData))
					}
				keepLargInx<-  seq(1, nrow(appextAuxData))
				}
			defKeepInx<- keepInx[keepLargInx]			
			## O_defKeepInx<- O_keepInx[O_keepLargInx]
#which(extAuxData$InChI_SingComp[O_defKeepInx]=="")
## length(defKeepInx)
## length(O_defKeepInx)
## length(keepInx)
## length(O_keepInx)
## length(defKeepInx[appDefKeep1Inx])
## appDefKeep1Inx<- which(extAuxData[O_defKeepInx, "ID"] %in% extAuxData[defKeepInx, "ID"] )
## length(intersect(extAuxData[O_defKeepInx, "ID"], extAuxData[defKeepInx, "ID"]))
			sinchiData<- inchiData[defKeepInx,]
			## O_sinchiData<- inchiData[O_defKeepInx,]
 
			selData<- extAuxData[defKeepInx,]
			## O_selData<- extAuxData[O_defKeepInx,]
#which(O_selData$InChI_SingComp=="")
			keep3Inx   <- keepUniqueComp(selData, vars=c("ID", "identCateg"))
			## O_keep3Inx <- keepUniqueComp(O_selData, vars=c("ID", "identCateg"))
			sAuxData<- selData[keep3Inx,]
			## O_sAuxData<- O_selData[O_keep3Inx,]
#which(O_sAuxData$InChI_SingComp=="")
			#sinchiData<- sinchiData[keep3Inx, -numIdS]
			sinchiData<- sinchiData[keep3Inx, ]
			## O_sinchiData<- O_sinchiData[O_keep3Inx, ]
			}
		#Si evita l'effetto del nodo Normalizer sulla carica nelle sostanze che non sono sali
		#Non solo: lo si evita anche quando il nodo Normalizer ha fallito

		if(neutralize=="yes")
			{neutralInx<- which(sAuxData[,"originalFeature"]=="neutral, organic"
			 			 |sAuxData[,"originalFeature"]=="neutral, metallic"
			 			 |sAuxData[,"originalFeature"]=="neutral, inorganic"
			 			 |sAuxData[,"originalFeature"]=="neutral, organic-inorganic")
##NB: quando si tratta di failedNormInx o failedCharge allora devi produrre un warning

			failedNormInx<- which(sAuxData[,"Charge"]!= 0)
			sAuxData[failedNormInx,"WarningFromNormal"]<- paste("notRecharged ", sAuxData[failedNormInx,"WarningFromNormal"], sep="")
			allInx<- union(neutralInx, failedNormInx)
			allInx<- allInx[order(allInx)]
			keepOrigInx<- which(extAuxData[,"normFlag"]==0)
			cumOrigInx<- c()
			cumOrigBackInx<- c()
			for (y in 1: length(sinchiData[allInx,"idS"]))
				{#y<- 1; y<- y+1
				setSearch<- sinchiData[allInx[y],"idS"]
				if(length(setSearch)>0)
					{cumOrigInx<- c(cumOrigInx, which(inchiData[keepOrigInx, "idS"]==setSearch))
					cumOrigBackInx<- c(cumOrigBackInx, which(sinchiData[allInx, "idS"]==setSearch) )
					}
				}

			#orig_neutralInx<- which(extAuxData[,"originalFeature"]=="neutral, organic"
			#			 |extAuxData[,"originalFeature"]=="neutral, metallic"
			#			 |extAuxData[,"originalFeature"]=="neutral, inorganic"
			#			 |extAuxData[,"originalFeature"]=="neutral, organic-inorganic")

 			#inxFromNotNorm <- which(inchiData[which(extAuxData[,"normFlag"]==0),"idS"] %in% sinchiData[neutralInx,"idS"])

			#sinchiData[neutralInx,"InChI"]<- O_sinchiData[neutralInx, "InChI"]
			sinchiData[allInx[cumOrigBackInx],"InChI"]<- extAuxData[keepOrigInx[cumOrigInx], "InChI_SingComp"]
#100-85-6[66854]<pubChem>
## appTTT<- gsub("100-85-6[66854]<pubChem>", "", sinchiData[allInx,"ID"], fixed=TRUE)
## foundInxProblem<- which(nchar(sinchiData[allInx,"ID"])> nchar(appTTT))
## sinchiData[allInx[foundInxProblem],"ID"]
## O_sinchiData[allInx[foundInxProblem],"ID"]

## sinchiData[allInx[foundInxProblem],"InChI"]
## O_sAuxData[allInx[foundInxProblem],"InChI_SingComp"]

## length( O_sinchiData[, "ID"])
## length( sinchiData[ , "ID"] )

## length(unique( O_sinchiData[, "ID"]))
## length(unique( sinchiData[ , "ID"] ))

## which(sinchiData[ , "ID"] %in% O_sinchiData[, "ID"])
##  O_sinchiData[, "ID"]
## length( O_sAuxData[, "ID"])
## length(  sAuxData[ , "ID"] )

## appTTT2<- gsub("980-71-2", "", O_sinchiData[neutralInx,"ID"], fixed=TRUE)
## foundInxProblem2<- which(nchar(O_sinchiData[neutralInx,"ID"])> nchar(appTTT2))
## O_sinchiData[neutralInx[foundInxProblem2],"ID"]

			#sinchiData[neutralInx,"InChI"]<- inchiData[inxFromNotNorm,"InChI"]
			#sinchiData<- sinchiData[, -numIdS]
			#sinchiData<- sinchiData[, c("ID", "InChI")]
			##sinchiData<- sinchiData[,numNotIdS,with=FALSE]
			sinchiData<- sinchiData[,numNotIdS]
			
			if(is.matrix(sinchiData))
				{sinchiData<-as.matrix(sinchiData)
				colnames(sinchiData)<- c("CAS_NO", "InChI")
				}
			
			}
		if(neutralize=="no")
			{keep1Inx<- which(extAuxData[,"normFlag"]==0)
			notNeutralData<- extAuxData[keep1Inx,]
			restrDat<- restructComp(notNeutralData, vars=c("ID", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal"), compVar="")
			firCompInx<- which(notNeutralData[, "compNum"]==1)
			sinchiData<- inchiData[keep1Inx[firCompInx],]
			restrDat<- restrDat[firCompInx,]
			if(orgType=="onlyOrganic")
				{keep2Inx<- which(restrDat[,"flagOrg"]==1)}
			if(orgType=="alsoOrgInorg")
				{keep2Inx<- which(restrDat[,"flagOrg"]==2|restrDat[,"flagOrg"]==1)}
			if(orgType=="all")
				{keep2Inx<-  seq(1, nrow(restrDat))}
			if(noMetallic=="yes")
				{keep3Inx<- which(restrDat[,"metalClass"]!=1)}
			if(noMetallic=="no")
				{keep3Inx<-  seq(1, nrow(restrDat))} 
			keepInx<-  intersect(keep2Inx, keep3Inx)
			sAuxData<- restrDat[keepInx,]
			#sinchiData<- sinchiData[keepInx,-numIdS]
			#sinchiData<- sinchiData[keepInx,c("ID", "InChI")]
			##sinchiData<- sinchiData[,numNotIdS,with=FALSE]
			sinchiData<- sinchiData[,numNotIdS]
			
			if(is.matrix(sinchiData))
				{sinchiData<-as.matrix(sinchiData)
				colnames(sinchiData)<- c("CAS_NO", "InChI")
				}
			
			}
		}
	if(length(sinchiData)[1]>0)
		{
## appTTT<- gsub("980-71-2", "", sinchiData[,"ID"], fixed=TRUE)
## foundInxProblem<- which(nchar(sinchiData[,"ID"])> nchar(appTTT))
## sinchiData[foundInxProblem,"ID"]

#>> sinchiData[,"ID"]
#>> which(is.na(sinchiData[,"ID"]))
#>> sinchiData[1:40,]
#>> colnames(sinchiData)
#>> which(sinchiData[,"InChI"]=="")
#>> which(is.na(sinchiData[,"InChI"]))

		inchiDataLs<- decompInChI(sinchiData, TRUE, FALSE)
		wAuxData<- sAuxData[,c("ElCount","MolWeight","saltClass","metalClass", "flagOrg","originalFeature", "WarningFromNormal","normFlag")]
#RICONTROLLA DA repackInChI perchè restituisce solo un elemento attivo di 80 (il campione totale era di 100)
		resLs<- repackInChI(inchiDataLs, wAuxData, #missData, 
							tautom, stereoChem, reduce, standard, extendCAS, bestMatch)
		repackData<- resLs[[1]]

## colnames(repackData)
## colnames(resLs[[2]])
## which(repackData[,"InChI"]=="")
## appTTT<- gsub("980-71-2", "", repackData[,"CAS_NO"], fixed=TRUE)
## foundInxProblem<- which(nchar(repackData[,"CAS_NO"])> nchar(appTTT))
## repackData[foundInxProblem,"CAS_NO"]

## appTTT<- gsub("980-71-2", "", resLs[[2]][,"ID"], fixed=TRUE)
## foundInxProblem<- which(nchar(resLs[[2]][,"ID"])> nchar(appTTT))
## resLs[[2]][foundInxProblem,]

#colnames(sinchiData)

## appTTT<- gsub("980-71-2", "", sinchiData[,"ID"], fixed=TRUE)
## foundInxProblem<- which(nchar(sinchiData[,"ID"])> nchar(appTTT))
## sinchiData[foundInxProblem,"ID"]


## appTTT<- gsub("980-71-2", "", names(inchiDataLs), fixed=TRUE)
## foundInxProblem<- which(nchar(names(inchiDataLs))> nchar(appTTT))
## names(inchiDataLs)[foundInxProblem]

		if(length(resLs)>1)
			{write.table(resLs[[2]],file=file.path(mainDir,"data/duplData.csv"),sep=";",na="",row.names=F) }
		if(!stereoChem)
			{if(!is.matrix(repackData)){repackData<- t(as.matrix(repackData))}
			appStereo<- repackData[,"WarningFromNormal"]
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)			
			appStereo<- gsub("STEREO_FORCED_BAD/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_ERROR", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR", "", appStereo, fixed = FALSE)
			repackData[,"WarningFromNormal"]<- appStereo
			}
		}
	#write options set
	numRec<- nrow(repackData)
	tautomer<- rep(tautom, numRec)
	stereoChemistry<- rep(stereoChem,numRec)
	appOption<- cbind(tautomer, stereoChemistry)
	neutralization<- rep(neutralize,numRec)
	appOption<- cbind(appOption, neutralization)
	if(neutralize=="no")
		{noSingleAtomIones<- rep("-",numRec)
		onlyHeavierComponent<- rep("-", numRec)
		}
	if(neutralize=="yes")
		{noSingleAtomIones<- rep(atomCounterIones,numRec)
		onlyHeavierComponent<- rep(largestComponent, numRec)
		}
	appOption<- cbind(appOption, noSingleAtomIones)	
	appOption<- cbind(appOption, onlyHeavierComponent)
	selectForOrganicity<- rep(orgType, numRec)
	appOption<- cbind(appOption, selectForOrganicity)
	noMetallicComponents<- rep(noMetallic, numRec)
	optionSet<- cbind(appOption, noMetallicComponents)
	colnames(optionSet)<- c("tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")
	repackData<- cbind(repackData, optionSet)
	return(repackData)
	}

#diagData<- wAuxData
#colnames(diagData)
 
####NEUTRALIZED
rewNeuInChI<- function(inchiData, auxData=extraData, #missData=missDat,
				tautom=c(0, 1, 2)[1], 
				stereoChem=FALSE,
				largestComponent=c("no", "yes", "largestOrganic")[3],
				orgType=c("all", "onlyOrganic", "alsoOrgInorg")[3],
				noMetallic=c("no", "yes")[2],
				reduce=FALSE, standard=FALSE, extendCAS=TRUE, bestMatch=TRUE,
				mainDir="E:/Data/Environment/projects/StandardSmiles/R") 
	{repackData<- c()
	if(!stereoChem)
		{fullVarNames <- c("CAS_NO", "InChI", "ElCount", "MolWeight", "saltClass", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal", "normFlag","notHInchi", "tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")}
	if(stereoChem)
		{fullVarNames <- c("CAS_NO", "InChI", "StereoDiagn", "ElCount", "MolWeight", "saltClass", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal", "normFlag","notHInchi", "tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")}
	
	atomCounterIones<- c("no", "yes")[2]
	if(!is.null(auxData))
		{numNotIdS<- which(colnames(inchiData)!="idS")
		originalFeature<- rep("", nrow(auxData))
		restrFullDat<- restructComp(auxData, vars=c("ID", "metalClass", "flagOrg", "WarningFromNormal","InChI_SingComp", "Charge"), compVar="")
		for (b in 1:length(originalFeature))
			{#b<- 1; b<- b+1
			feaSalt<- ifelse(restrFullDat$saltClass[b]==1, "salt", "neutral")
			if(restrFullDat$metalClass[b]==1)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==1, "organometallic", ifelse(restrFullDat$flagOrg[b]==-1,"unkonwn-metallic", "metallic"))}
			if(restrFullDat$metalClass[b]==0)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==2, "organic-inorganic", ifelse(restrFullDat$flagOrg[b]==1, "organic",  ifelse(restrFullDat$flagOrg[b]==-1,"unknown", "inorganic")))}
			originalFeature[b]<- paste(feaSalt, feaOrgMetal, sep=", ")
			}	

			
 		extAuxData<- cbind(auxData, originalFeature)
		colnames(extAuxData)<- c(colnames(auxData), "originalFeature")
		collectExcluded<- cbind(as.character(extAuxData[,"ID"]), rep("", nrow(extAuxData)))
		colnames(collectExcluded)<- c("ID", "Info")
        	sinchiData<- inchiData
		colOfSinchi<- c("CAS_NO", "InChI", "idS")
		colnames(sinchiData)<- colOfSinchi
		colOfAuxData<- colnames(extAuxData)
		keep1Inx<- seq(1, nrow(extAuxData))
		if(noMetallic=="yes")
			{keepMetCmpInx<- which(extAuxData[,"metalClass"]==1)
			metMol<- unique(extAuxData[keepMetCmpInx,"ID"])
			#names(metMol)<- NULL
			complKeepNotMetInx<-c()
			if(length(metMol)>0)
				for (j in 1: length(metMol))
					{#j<- 1; j<- j+1
						{complKeepNotMetInx<- c(complKeepNotMetInx, which(extAuxData[,"ID"]==metMol[j]))}
					}
			keep1Inx<- setdiff(keep1Inx,complKeepNotMetInx)
			if(length(keep1Inx)>0)
				{extAuxData<- extAuxData[keep1Inx, ]	
				sinchiData<- sinchiData[keep1Inx, ]
				if(!is.matrix(sinchiData))
					{sinchiData<-t(as.matrix(sinchiData))
					extAuxData<-t(as.matrix(extAuxData))
					colnames(sinchiData)<- colOfSinchi
					colnames(extAuxData)<- colOfAuxData
					}
				}
			if(length(keep1Inx)==0)
				{extAuxData<- extAuxData[keep1Inx, ]	
				sinchiData<- sinchiData[keep1Inx, ]
				}
			}			
		testEmptyDat<- nrow(extAuxData)>0	
     	if(orgType=="onlyOrganic" & testEmptyDat)
			{keep2Inx<- which(extAuxData[,"flagOrg"]==1)
			if(length(keep2Inx)>0)
				{extAuxData<- extAuxData[keep2Inx, ]	
			    sinchiData<- sinchiData[keep2Inx, ]
				if(!is.matrix(sinchiData))
					{sinchiData<-t(as.matrix(sinchiData))
					extAuxData<-t(as.matrix(extAuxData))
					}
				}

			if(length(keep2Inx)==0)
				{extAuxData<- extAuxData[keep2Inx, ]	
				sinchiData<- sinchiData[keep2Inx, ]
				}
			}
		if(orgType=="alsoOrgInorg" & testEmptyDat)
			{keepInorgCmpInx<- which(extAuxData[,"flagOrg"]!=1)
			inorgMol<- unique(extAuxData[keepInorgCmpInx,"ID"])
			complKeepInorgInx<-c()
			if(length(inorgMol)>0)
				for (j in 1: length(inorgMol))
					{#j<- 1; j<- j+1
						{TestInorgComp<- prod((extAuxData[which(extAuxData[,"ID"]==inorgMol[j]),"flagOrg"]==0)*1)==1
						if(TestInorgComp)
							{complKeepInorgInx<- c(complKeepInorgInx, which(extAuxData[,"ID"]==inorgMol[j]))}
						}
					}
			complKeepInorgInx<- unique(complKeepInorgInx)
			keep2Inx<-  setdiff(seq(1, nrow(extAuxData)),complKeepInorgInx)
			if(length(keep2Inx)>0)
				{extAuxData<- extAuxData[keep2Inx, ]	
				sinchiData<- sinchiData[keep2Inx, ]
				if(!is.matrix(sinchiData))
					{sinchiData<-t(as.matrix(sinchiData))
					extAuxData<-t(as.matrix(extAuxData))
					}	
				}
			if(length(keep2Inx)==0)
				{extAuxData<- extAuxData[keep2Inx, ]	
				sinchiData<- sinchiData[keep2Inx, ]
				}
			}
		testEmptyDat<- nrow(extAuxData)>0	
		#basic neutralization options	
		if(testEmptyDat)
			{keep3aInx   <- which(extAuxData[,"eraseComp"]==0)
			keep3bInx   <- keepUniqueComp(extAuxData, vars=c("ID", "identCateg"))
			keep3Inx<- intersect(keep3aInx, keep3bInx)
			keep3Inx<- keep3Inx[order(keep3Inx)]
			complKeep3Inx<- setdiff(seq(1, nrow(extAuxData)),keep3Inx)
			if(length(keep3Inx)>0)
				{extAuxData<- extAuxData[keep3Inx, ]	
				sinchiData<- sinchiData[keep3Inx, ]
				if(!is.matrix(sinchiData))
					{sinchiData<-t(as.matrix(sinchiData))
					extAuxData<-t(as.matrix(extAuxData))
					colnames(sinchiData)<- colOfSinchi
                	colnames(extAuxData)<- colOfAuxData
					}
				}
			if(length(keep3Inx)==0)
				{extAuxData<- extAuxData[keep3Inx, ]	
				sinchiData<- sinchiData[keep3Inx, ]
				}
			}
		testEmptyDat<- nrow(extAuxData)>0	
			
		if(largestComponent=="yes" & testEmptyDat) 
			{keep4Inx<- keepLargestComp(extAuxData, vars=c("ID", "MolWeight"), "all")
			if(length(keep4Inx)>0)
				{extAuxData<- extAuxData[keep4Inx, ]	
				sinchiData<- sinchiData[keep4Inx, ]			
				if(!is.matrix(sinchiData))
					{sinchiData<-t(as.matrix(sinchiData))
					extAuxData<-t(as.matrix(extAuxData))
					colnames(sinchiData)<- colOfSinchi
					colnames(extAuxData)<- colOfAuxData
					}
				}
			if(length(keep4Inx)==0)
				{extAuxData<- extAuxData[keep4Inx, ]	
				sinchiData<- sinchiData[keep4Inx, ]
				}
			}			
		sAuxData<- extAuxData
	
		}
		
	failMol <- setdiff(inchiData[,"ID"],sinchiData[,"CAS_NO"]) 
	emtMat<- matrix("",nrow=length(failMol),ncol=length(fullVarNames)) 
	colnames(emtMat)<-fullVarNames 
	emtMat[,"CAS_NO"] <- failMol 
	emtMat[,"InChI"] <- "Excluded by options" 
	
	testEmtSinchi <- nrow(sinchiData)==0 #
	if(!testEmtSinchi)
		{inchiDataLs<- decompInChI(sinchiData, TRUE, FALSE)
		wAuxData<- sAuxData[,c("ElCount","MolWeight","saltClass","metalClass", "flagOrg","originalFeature", "WarningFromNormal","normFlag")]
		resLs<- repackInChI(inchiDataLs, wAuxData, #missData, 
							tautom, stereoChem, reduce, standard, extendCAS, bestMatch)
		#repackData<- resLs[[1]]
		repackData<- resLs
		#if(length(resLs)>1)
		#	{write.table(resLs[[2]],file=file.path(mainDir,"data/duplData.csv"),sep=";",na="",row.names=F) }
		if(!stereoChem)
			{if(!is.matrix(repackData)){repackData<- t(as.matrix(repackData))}
			appStereo<- repackData[,"WarningFromNormal"]
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)			
			appStereo<- gsub("STEREO_FORCED_BAD/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_ERROR", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR", "", appStereo, fixed = FALSE)
			repackData[,"WarningFromNormal"]<- appStereo
			}
		#write options set
		numRec<- nrow(repackData)
		tautomer<- rep(tautom, numRec)
		stereoChemistry<- rep(stereoChem,numRec)
		appOption<- cbind(tautomer, stereoChemistry)
		neutralization<- rep("yes",numRec)
		appOption<- cbind(appOption, neutralization)
		noSingleAtomIones<- rep(atomCounterIones,numRec)
		onlyHeavierComponent<- rep(largestComponent, numRec)
		appOption<- cbind(appOption, noSingleAtomIones)	
		appOption<- cbind(appOption, onlyHeavierComponent)
		selectForOrganicity<- rep(orgType, numRec)
		appOption<- cbind(appOption, selectForOrganicity)
		noMetallicComponents<- rep(noMetallic, numRec)
		optionSet<- cbind(appOption, noMetallicComponents)
		colnames(optionSet)<- c("tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")
		repackData<- cbind(repackData, optionSet)
		repackExtData<- rbind(repackData, emtMat)
		} 
	if(testEmtSinchi) 
		{repackExtData<- emtMat 
		} 
		
	return(repackExtData)
	}



####NOTNEUTRALIZED
rewNotNeuInChI<- function(inchiData, auxData=extraData, #missData=missDat,
				tautom=c(0, 1, 2)[1], 
				stereoChem=FALSE, 
				largestComponent=c("no", "yes", "largestOrganic")[1],
				orgType=c("all", "onlyOrganic", "alsoOrgInorg")[3],
				noMetallic=c("no", "yes")[2],
				reduce=FALSE, standard=FALSE, extendCAS=TRUE, bestMatch=TRUE,
				mainDir="E:/Data/Environment/projects/StandardSmiles/R") 
	{repackData<- c()
	atomCounterIones<- c("no", "yes")[1]
	
	if(!stereoChem)
		{fullVarNames <- c("CAS_NO", "InChI", "ElCount", "MolWeight", "saltClass", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal", "normFlag", "notHInchi","tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")}
	if(stereoChem)
		{fullVarNames <- c("CAS_NO", "InChI", "StereoDiagn", "ElCount", "MolWeight", "saltClass", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal", "normFlag","notHInchi", "tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")}
		
	if(!is.null(auxData))
		{numNotIdS<- which(colnames(inchiData)!="idS")
		originalFeature<- rep("", nrow(auxData))
		restrFullDat<- restructComp(auxData, vars=c("ID", "metalClass", "flagOrg", "WarningFromNormal","InChI_SingComp", "Charge"), compVar="")
		
		for (b in 1:length(originalFeature))
			{#b<- 1; b<- b+1
			feaSalt<- ifelse(restrFullDat$saltClass[b]==1, "salt", "neutral")
			if(restrFullDat$metalClass[b]==1)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==1|restrFullDat$flagOrg[b]==2, "organometallic", ifelse(restrFullDat$flagOrg[b]==-1,"unkonwn-metallic", "metallic"))}
			if(restrFullDat$metalClass[b]==0)
				{feaOrgMetal<- ifelse(restrFullDat$flagOrg[b]==2, "organic-inorganic", ifelse(restrFullDat$flagOrg[b]==1, "organic",  ifelse(restrFullDat$flagOrg[b]==-1,"unknown", "inorganic")))}
			originalFeature[b]<- paste(feaSalt, feaOrgMetal, sep=", ")
			}
		extAuxData<- cbind(auxData, originalFeature)
		colnames(extAuxData)<- c(colnames(auxData), "originalFeature")
		notNeutralData<- extAuxData
		restrDat<- restructComp(notNeutralData, vars=c("ID", "metalClass", "flagOrg", "originalFeature", "WarningFromNormal"), compVar="")
		firCompInx<- which(notNeutralData[, "compNum"]==1)
		sinchiData<- inchiData[firCompInx,]
		if(!is.matrix(sinchiData))
			{sinchiData<-t(as.matrix(sinchiData))}
		restrDat<- restrDat[firCompInx,]
		if(!is.matrix(restrDat))
			{restrDat<-as.matrix(restrDat)}

		collectExcluded<- cbind(as.character(restrDat[,"ID"]), rep("", nrow(restrDat)))
		colnames(collectExcluded)<- c("ID", "Info")

		if(orgType=="onlyOrganic")
			{keep2Inx<- which(restrDat[,"flagOrg"]==1)
			complKeep2Inx<- setdiff(seq(1, nrow(restrDat)),keep2Inx)
			if(length(complKeep2Inx)>0)
				{collectExcluded[complKeep2Inx,"Info"]<- rep("notOrganic", length(complKeep2Inx)) }		
			}
		if(orgType=="alsoOrgInorg")
			{keep2Inx<- which(restrDat[,"flagOrg"]==2|restrDat[,"flagOrg"]==1)
			complKeep2Inx<-  setdiff(seq(1, nrow(restrDat)),keep2Inx)
			if(length(complKeep2Inx)>0)
				{collectExcluded[complKeep2Inx,"Info"]<- rep("fullyInorganic",length(complKeep2Inx)) }		
			}
		if(orgType=="all")
			{keep2Inx<-  seq(1, nrow(restrDat))}
		if(noMetallic=="yes")
			{keep3Inx<- which(restrDat[,"metalClass"]!=1)
			complKeep3Inx<-  setdiff(seq(1, nrow(restrDat)),keep3Inx)
			if(length(complKeep3Inx)>0)
				{collectExcluded[complKeep3Inx,"Info"]<- paste(collectExcluded[complKeep3Inx,"Info"], rep("metallic",length(complKeep3Inx)), sep="_") }
			}
		if(noMetallic=="no")
			{keep3Inx<-  seq(1, nrow(restrDat))} 
		keepInx<-  intersect(keep2Inx, keep3Inx)
		sAuxData<- restrDat[keepInx,]
		
		if(!is.matrix(sAuxData))
			{sAuxData<-t(as.matrix(sAuxData))}
		
		sinchiData<- sinchiData[keepInx,numNotIdS]
		
		if(is.matrix(sinchiData))
			{sinchiData<-as.matrix(sinchiData)
			colnames(sinchiData)<- c("CAS_NO", "InChI")
			}
		
		if(!is.matrix(sinchiData))
			{sinchiData<-t(as.matrix(sinchiData))
			colnames(sinchiData)<- c("CAS_NO", "InChI")
			}
		}
		
	failMol <- setdiff(inchiData[,"ID"],sinchiData[,"CAS_NO"])
	emtMat<- matrix("",nrow=length(failMol),ncol=length(fullVarNames))
	colnames(emtMat)<-fullVarNames
	emtMat[,"CAS_NO"] <- failMol
	
	testEmtSinchi <- length(sinchiData)[1]==0
	
	if(!testEmtSinchi)
		{inchiDataLs<- decompInChI(sinchiData, TRUE, FALSE)
		wAuxData<- sAuxData[,c("ElCount","MolWeight","saltClass","metalClass", "flagOrg","originalFeature", "WarningFromNormal","normFlag")]
		#if(!is.matrix(wAuxData))
		#{wAuxData<-t(as.matrix(wAuxData))
		#}
		resLs<- repackInChI(inchiDataLs, wAuxData, #missData, 
							tautom, stereoChem, reduce, standard, extendCAS, bestMatch)
		#repackData<- resLs[[1]]
		repackData<- resLs
		#if(length(resLs)>1)
		#	{write.table(resLs[[2]],file=file.path(mainDir,"data/duplData.csv"),sep=";",na="",row.names=F) }
		if(!stereoChem)
			{if(!is.matrix(repackData)){repackData<- t(as.matrix(repackData))}
			appStereo<- repackData[,"WarningFromNormal"]
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("DUBIOUS_STEREO_REMOVED", "", appStereo, fixed = FALSE)			
			appStereo<- gsub("STEREO_FORCED_BAD/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_FORCED_BAD", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_TRANSFORMED", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR/", "", appStereo, fixed = FALSE)
			appStereo<- gsub("/STEREO_ERROR", "", appStereo, fixed = FALSE)
			appStereo<- gsub("STEREO_ERROR", "", appStereo, fixed = FALSE)
			repackData[,"WarningFromNormal"]<- appStereo
			}
		#write options set
		numRec<- nrow(repackData)
		tautomer<- rep(tautom, numRec)
		stereoChemistry<- rep(stereoChem,numRec)
		appOption<- cbind(tautomer, stereoChemistry)
		neutralization<- rep("no",numRec)
		appOption<- cbind(appOption, neutralization)
		noSingleAtomIones<- rep("-",numRec)
		onlyHeavierComponent<- rep("-", numRec)
		appOption<- cbind(appOption, noSingleAtomIones)	
		appOption<- cbind(appOption, onlyHeavierComponent)
		selectForOrganicity<- rep(orgType, numRec)
		appOption<- cbind(appOption, selectForOrganicity)
		noMetallicComponents<- rep(noMetallic, numRec)
		optionSet<- cbind(appOption, noMetallicComponents)
		colnames(optionSet)<- c("tautomer", "stereoChemistry", "neutralization", "noSingleAtomIones", "onlyHeavierComponent", "selectForOrganicity", "noMetallicComponents")
		repackData<- cbind(repackData, optionSet)
		colSet<- colnames(repackData)

		collResExcluded<- collectExcluded[which(collectExcluded[,"Info"]!=""), ]
		if(!is.matrix(collResExcluded)){collResExcluded<- t(as.matrix(collResExcluded))}
		uniExclMolSet<- unique(collResExcluded[,"ID"])
		repackSet<- unique(repackData[,"CAS_NO"])
		if(length(uniExclMolSet)>1)
			{for (x in 1: length(uniExclMolSet))
				{#x<- 1; x<- x+1
				setUni<- uniExclMolSet[x]
				dblInx<- which(collResExcluded[,"ID"]==setUni)
				appString<- ""
				collResExcluded[dblInx[1],"Info"]<- ifelse(substr(collResExcluded[dblInx[1],"Info"], 1, 1)=="_", substr(collResExcluded[dblInx[1],"Info"], 2, nchar(collResExcluded[dblInx[1],"Info"])), collResExcluded[dblInx[1],"Info"])
				if(length(dblInx)>1)
					{for (y in 2:length(dblInx))
						{#y<- 2; y<- y+1
						collResExcluded[dblInx[y],"ID"]<- ""
						modString<- ifelse(substr(collResExcluded[dblInx[y],"Info"], 1, 1)=="_", substr(collResExcluded[dblInx[y],"Info"], 2, nchar(collResExcluded[dblInx[y],"Info"])), collResExcluded[dblInx[y],"Info"])
						appString<- paste(unique(c(appString, modString)), sep="_")
						}
					if(length(which(appString==""))>0) {appString<- appString[-which(appString=="")]}
					collResExcluded[dblInx[1],"Info"]<- paste(unique(c(collResExcluded[dblInx[1],"Info"],appString)), collapse="_")
					testRepack<- length(which(repackSet==setUni))>0
					if(testRepack){collResExcluded[dblInx[1],"ID"]<- ""}
					}
				}
			}
		appendExclud<- matrix("", nrow=length(which(collResExcluded[,"ID"]!="")), ncol=length(colSet))
		colnames(appendExclud)<- colSet
		appendExclud[,"CAS_NO"]<- collResExcluded[which(collResExcluded[,"ID"]!=""), "ID"]
		appendExclud[,"InChI"]<- collResExcluded[which(collResExcluded[,"ID"]!=""), "Info"]
		repackExtData<- rbind(repackData,appendExclud)
		}
	
	if(testEmtSinchi) 
		{repackExtData<- emtMat} 

	if(!testEmtSinchi) 
		{repackExtData<- rbind(repackExtData, emtMat)
		} 
		
	return(repackExtData)
	}


repackInChI<- function(inchiDataLs, diagData, #missData, 
					   tautom, stereoChem, reduce, standard, extendCAS, bestMatch) 
	{#VER/(chqpbtms/i(hbtms)/f(hqbtms/i(btms)o))/r(chqpbtms/i(hbtms)/f(hqbtms/i(btms)o))
	resLs<- list()
	ID_names<- names(inchiDataLs)
	if(!is.matrix(diagData))
		{diagData<-t(as.matrix(diagData))}
	diagVars<- colnames(diagData) 
	fullColNames<- c("CAS_NO", "InChI", "StereoDiagn")
	if(stereoChem)
		{allVars<- union(fullColNames, diagVars)
		repackData<- matrix("", nrow=0, ncol=length(allVars))
		colnames(repackData)<- allVars
		}
	if(!stereoChem)
		{allVars<- union(fullColNames[-3], diagVars)
		repackData<- matrix("", nrow=0, ncol=length(allVars))
		colnames(repackData)<- allVars
		}
	addedRecord<- rep("", ncol(repackData))
	cumInchi <- c()
	for (h in 1: length(inchiDataLs))#
		{#h<- 1; h<- h+1
		repackData<- rbind(repackData, addedRecord)
		repackData[h,"CAS_NO"] <- ID_names[h]
		if(standard) 
			{repackInchi<- "InChI=1S"
			repackInchi2<- "InChI=1S"}
		if(!standard) 
			{repackInchi<- "InChI=1"
			repackInchi2<- "InChI=1"
			}
		nameOfLayer<- "/"
		currPack<- inchiDataLs[[h]]
		namesOfPack<- names(currPack)
		diagnSterRec<- "0"
		diagnStereo<-  "0"
		for (j in 1: length(currPack))
			{#j<- 1; j<- j+1
			appNam<- unlist(strsplit(namesOfPack[j], "_", fixed=TRUE))
			testSplit<- length(appNam)>1 
			if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }			
			subPack<- currPack[[j]]
			namesOfSubPack<- names(subPack)
			for (k in 1: length(subPack))
				{#k<- 1; k<- k+1
				if(!is.null(namesOfSubPack))
					{appNam<- unlist(strsplit(namesOfSubPack[k], "_", fixed=TRUE))
					testSplit<- length(appNam)>1
					}
				if(is.null(namesOfSubPack)){testSplit<- FALSE}
				if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }
				sub2Pack<- subPack[[k]]
				namesOfSub2Pack<- names(sub2Pack)
				for (z in 1: length(sub2Pack))
					{#z<- 1; z<- z+1
					if(!is.null(namesOfSub2Pack))
						{appNam<- unlist(strsplit(namesOfSub2Pack[z], "_", fixed=TRUE))
						testSplit<- length(appNam)>1
						}
					if(is.null(namesOfSub2Pack)){testSplit<- FALSE}
					if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }
					sub3Pack<- sub2Pack[[z]]
					if(sub3Pack!="")
						{repackInchi2<- paste(repackInchi2, nameOfLayer, sub3Pack, collapse="", sep="")
						if(tautom==0)
							{repackInchi<- paste(repackInchi, nameOfLayer, sub3Pack, collapse="", sep="") }
						if(tautom==1)
							{if(nameOfLayer!="/f")
								{repackInchi<- paste(repackInchi, nameOfLayer, sub3Pack, collapse="", sep="") 
								}							
							}
						if(tautom==2)
							{if(nameOfLayer!="/f")
								{appFragms<- unlist(strsplit(sub3Pack, "/", fixed=TRUE))
								appNams<- substr(appFragms, 1,1)
								locOfNotLayer<- which(appNams!="h" & appNams!="")
								sub4Pack<- paste(appFragms[locOfNotLayer], "/", collapse="", sep="")
								repackInchi<- paste(repackInchi, nameOfLayer, sub4Pack, collapse="", sep="")
								}
							}
						if(stereoChem)
							{appFragms<- unlist(strsplit(repackInchi, "/", fixed=TRUE))
							appNams<- substr(appFragms, 1,1)
							
							locOfStereoLayers<- which(appNams=="b" | appNams=="t" | appNams=="m" | appNams=="s" )
							if(length(locOfStereoLayers)>0)
								{for (d in 1: length(locOfStereoLayers))
									{#d<- 1; d<- d+1
									uncerStereo<- appFragms[locOfStereoLayers[d]]
									appUn1<- unlist(strsplit(uncerStereo, "?", fixed=TRUE))
									appUn2<- unlist(strsplit(uncerStereo, "u", fixed=TRUE))
									testUn1<- (length(appUn1)>1 | (nchar(appUn1[1])<nchar(uncerStereo)))
									testUn2<- (length(appUn2)>1 | (nchar(appUn2[1])<nchar(uncerStereo)))
									diagnStereo<- ifelse(!testUn1 & !testUn2, "0", ifelse(testUn1, "?", "u"))
									}
								}
							if(diagnSterRec=="0")
								{diagnSterRec<- diagnStereo}
							if(diagnSterRec=="?" & diagnStereo=="u")
								{diagnSterRec<- "?u"}
							if(diagnSterRec=="u" & diagnStereo=="?")
								{diagnSterRec<- "?u"}
							}
						if(!stereoChem)
							{appFragms<- unlist(strsplit(repackInchi, "/", fixed=TRUE))
							appNams<- substr(appFragms, 1,1)
							locOfNotStereoLayers<- which(appNams!="b" & appNams!="t" & appNams!="m" & appNams!="s" )
							repackInchi<- paste(appFragms[locOfNotStereoLayers], "/", collapse="", sep="")							
							appFragms2<- unlist(strsplit(repackInchi2, "/", fixed=TRUE))
							appNams2<- substr(appFragms2, 1,1)
							locOfNotStereoLayers2<- which(appNams2!="b" & appNams2!="t" & appNams2!="m" & appNams2!="s" & appNams2!="")
							repackInchi2<- paste(appFragms2[locOfNotStereoLayers2],  "/",  collapse="", sep="")
							}								
						}
					}				
				}
			}
		repackData[h,"InChI"] <- repackInchi2
		##repackData[h,diagVars] <- as.matrix(diagData[h,diagVars,with=FALSE])
		
		if(!is.matrix(diagData))
			{diagData<-t(as.matrix(diagData))
			repackData[h,diagVars] <- as.matrix(diagData[h,diagVars])
			}
		
		if(is.matrix(diagData))
			{repackData[h,diagVars] <- as.matrix(diagData[h,diagVars])
			}	
		cumInchi <- c(cumInchi, repackInchi)
		if(stereoChem) 
			{appFragms3<- unlist(strsplit(repackInchi2, "/", fixed=TRUE))
			appNams3<- substr(appFragms3, 1,1)
			locOfNotStereoLayers3<- which(appNams3=="b" | appNams3=="t" | appNams3=="m" | appNams3=="s")
			if(length(locOfNotStereoLayers3)==0 & diagnSterRec==0)
				{diagnSterRec<- -1}
			repackData[h,"StereoDiagn"] <- diagnSterRec
			}
		}
	rownames(repackData)<- NULL

	#if (tautom==2) 
		#{UnicumInchi<-unique(cumInchi)
		#duplInchi<-c()
		#uniDuplLog<-c()
		#uniDuplInx<-seq(1, length(cumInchi))
		#for (q in 1: length(UnicumInchi))
			#{#q<-1;  q<-q+1#; q<-30
			#appDuplInx<-which(cumInchi==UnicumInchi[q])
			#testDupl<-length(appDuplInx)>1
			#if (testDupl)
				#{duplInx<-appDuplInx[-1]
				#uniDuplInx[duplInx]<-rep(appDuplInx[1],length(duplInx))
				#}
			#}
		#cUniDuplInx<-unique(uniDuplInx)
		
		#for (z in 1: length(cUniDuplInx))
			#{#z<-1;  z<-z+1#;
			#currInx<-cUniDuplInx[z]
			#dplInx<-which(uniDuplInx==currInx)
			#repackData[dplInx,"InChI"]<-rep(repackData[dplInx[1],"InChI"],length(dplInx))
			#}
		#}

	#if (tautom!=2) 
		#{cUniDuplInx<-seq(1, length(cumInchi))
		#}
		
	notHydrInchi<-as.matrix(cumInchi)

	#if(extendCAS)
		#{monlyCAS<- ID_names
		#unOnlyCAS<- unique(ID_names)		
		#elimInx<- setdiff(seq(1,nrow(repackData)), cUniDuplInx)
		#elimCAS<- repackData[elimInx, c("CAS_NO","InChI")]
		#repackData<- repackData[cUniDuplInx, ]
		#if(!is.matrix(repackData))
			#{repackData<-t(as.matrix(repackData))}
		#}
		
	newNames<-c(colnames(repackData),"notHInchi")
	repackData<-cbind(repackData,notHydrInchi)
	colnames(repackData)<-newNames
	
	return(repackData)
	}


#colnames(repackData)
#repackData[,"CAS_NO"]

howManyComp<- function(currMolec)
	{mainForm<- unlist(currMolec)[1]
	appForm<- unlist(strsplit(mainForm, "/", fixed=TRUE))[1]
	appForm<- unlist(strsplit(appForm, ".", fixed=TRUE))
	numComp<- length(appForm)
	return(numComp)
	}

ex_repackInChI<- function(inchiDataLs, diagData, missData, tautom, stereoChem, reduce, standard, extendCAS, bestMatch) 
	{#VER/(chqpbtms/i(hbtms)/f(hqbtms/i(btms)o))/r(chqpbtms/i(hbtms)/f(hqbtms/i(btms)o))
	resLs<- list()
	ID_names<- names(inchiDataLs)
	diagVars<- colnames(diagData) 
	fullColNames<- c("CAS_NO", "InChI", "StereoDiagn")
	if(stereoChem)
		{allVars<- union(fullColNames, diagVars)
		repackData<- matrix("", nrow=0, ncol=length(allVars))
		colnames(repackData)<- allVars
		}
	if(!stereoChem)
		{allVars<- union(fullColNames[-3], diagVars)
		repackData<- matrix("", nrow=0, ncol=length(allVars))
		colnames(repackData)<- allVars
		}
	addedRecord<- rep("", ncol(repackData))
	cumInchi <- c()
	for (h in 1: length(inchiDataLs))#
		{#h<- 1; h<- h+1
		repackData<- rbind(repackData, addedRecord)
		repackData[h,"CAS_NO"] <- ID_names[h]
		if(standard) 
			{repackInchi<- "InChI=1S"
			repackInchi2<- "InChI=1S"}
		if(!standard) 
			{repackInchi<- "InChI=1"
			repackInchi2<- "InChI=1"
			}
		nameOfLayer<- "/"
		currPack<- inchiDataLs[[h]]
		namesOfPack<- names(currPack)
		diagnSterRec<- "0"
		diagnStereo<-  "0"
		for (j in 1: length(currPack))
			{#j<- 1; j<- j+1
			appNam<- unlist(strsplit(namesOfPack[j], "_", fixed=TRUE))
			testSplit<- length(appNam)>1 
			if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }			
			subPack<- currPack[[j]]
			namesOfSubPack<- names(subPack)
			for (k in 1: length(subPack))
				{#k<- 1; k<- k+1
				if(!is.null(namesOfSubPack))
					{appNam<- unlist(strsplit(namesOfSubPack[k], "_", fixed=TRUE))
					testSplit<- length(appNam)>1
					}
				if(is.null(namesOfSubPack)){testSplit<- FALSE}
				if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }
				sub2Pack<- subPack[[k]]
				namesOfSub2Pack<- names(sub2Pack)
				for (z in 1: length(sub2Pack))
					{#z<- 1; z<- z+1
					if(!is.null(namesOfSub2Pack))
						{appNam<- unlist(strsplit(namesOfSub2Pack[z], "_", fixed=TRUE))
						testSplit<- length(appNam)>1
						}
					if(is.null(namesOfSub2Pack)){testSplit<- FALSE}
					if(testSplit) {nameOfLayer<- paste("/", appNam[1], collapse="", sep="") }
					sub3Pack<- sub2Pack[[z]]
					if(sub3Pack!="")
						{repackInchi2<- paste(repackInchi2, nameOfLayer, sub3Pack, collapse="", sep="")
						if(tautom==0)
							{repackInchi<- paste(repackInchi, nameOfLayer, sub3Pack, collapse="", sep="") }
						if(tautom==1)
							{if(nameOfLayer!="/f")
								{repackInchi<- paste(repackInchi, nameOfLayer, sub3Pack, collapse="", sep="") 
								}							
							}
						if(tautom==2)
							{if(nameOfLayer!="/f")
								{appFragms<- unlist(strsplit(sub3Pack, "/", fixed=TRUE))
								appNams<- substr(appFragms, 1,1)
								locOfNotLayer<- which(appNams!="h" & appNams!="")
								sub4Pack<- paste(appFragms[locOfNotLayer], "/", collapse="", sep="")
								repackInchi<- paste(repackInchi, nameOfLayer, sub4Pack, collapse="", sep="")
								}
							}
						if(stereoChem)
							{appFragms<- unlist(strsplit(repackInchi, "/", fixed=TRUE))
							appNams<- substr(appFragms, 1,1)
							
							locOfStereoLayers<- which(appNams=="b" | appNams=="t" | appNams=="m" | appNams=="s" )
							if(length(locOfStereoLayers)>0)
								{for (d in 1: length(locOfStereoLayers))
									{#d<- 1; d<- d+1
									uncerStereo<- appFragms[locOfStereoLayers[d]]
									appUn1<- unlist(strsplit(uncerStereo, "?", fixed=TRUE))
									appUn2<- unlist(strsplit(uncerStereo, "u", fixed=TRUE))
									testUn1<- (length(appUn1)>1 | (nchar(appUn1[1])<nchar(uncerStereo)))
									testUn2<- (length(appUn2)>1 | (nchar(appUn2[1])<nchar(uncerStereo)))
									diagnStereo<- ifelse(!testUn1 & !testUn2, "0", ifelse(testUn1, "?", "u"))
									}
								}
							if(diagnSterRec=="0")
								{diagnSterRec<- diagnStereo}
							if(diagnSterRec=="?" & diagnStereo=="u")
								{diagnSterRec<- "?u"}
							if(diagnSterRec=="u" & diagnStereo=="?")
								{diagnSterRec<- "?u"}
							}
						if(!stereoChem)
							{appFragms<- unlist(strsplit(repackInchi, "/", fixed=TRUE))
							appNams<- substr(appFragms, 1,1)
							locOfNotStereoLayers<- which(appNams!="b" & appNams!="t" & appNams!="m" & appNams!="s" )
							repackInchi<- paste(appFragms[locOfNotStereoLayers], "/", collapse="", sep="")							
							appFragms2<- unlist(strsplit(repackInchi2, "/", fixed=TRUE))
							appNams2<- substr(appFragms2, 1,1)
							locOfNotStereoLayers2<- which(appNams2!="b" & appNams2!="t" & appNams2!="m" & appNams2!="s" & appNams2!="")
							repackInchi2<- paste(appFragms2[locOfNotStereoLayers2],  "/",  collapse="", sep="")
							}								
						}
					}				
				}
			}
		repackData[h,"InChI"] <- repackInchi2
		##repackData[h,diagVars] <- as.matrix(diagData[h,diagVars,with=FALSE])
		repackData[h,diagVars] <- as.matrix(diagData[h,diagVars])
		cumInchi <- c(cumInchi, repackInchi)
		if(stereoChem) 
			{appFragms3<- unlist(strsplit(repackInchi2, "/", fixed=TRUE))
			appNams3<- substr(appFragms3, 1,1)
			locOfNotStereoLayers3<- which(appNams3=="b" | appNams3=="t" | appNams3=="m" | appNams3=="s")
			if(length(locOfNotStereoLayers3)==0 & diagnSterRec==0)
				{diagnSterRec<- -1}
			repackData[h,"StereoDiagn"] <- diagnSterRec
			}
		}
	rownames(repackData)<- NULL
	if(extendCAS)
		{monlyCAS<- ID_names
		unOnlyCAS<- unique(ID_names)
		newId<- c()
		cumInx<- c()
		for (j in 1: length(unOnlyCAS))
			{#j<- 1; j<- j+1
			currCASinx<- which(monlyCAS==unOnlyCAS[j])
			currInchiSet<- cumInchi[currCASinx] 
			currIDSet<- repackData[currCASinx,"CAS_NO"]
			findNotDuplInx<- which(!duplicated(currInchiSet))
			testOneFound<- length(findNotDuplInx)==1
			if(testOneFound | bestMatch)
				{##OKKIO##>>> newId<- c(newId, unOnlyCAS[j])
				newId<- c(newId, currIDSet[findNotDuplInx])
				cumInx<- c(cumInx, currCASinx[1])
				}
			if(!testOneFound & !bestMatch)
				{findNotDuplIDs<- currIDSet[findNotDuplInx]
				newId<- c(newId, findNotDuplIDs)
				cumInx<- c(cumInx, currCASinx[findNotDuplInx])
				}			
			}
		elimInx<- setdiff(seq(1,nrow(repackData)), cumInx)
		elimCAS<- repackData[elimInx, c("CAS_NO","InChI")]
		names(newId)<- NULL
		repackData<- repackData[cumInx, ]
		repackData[,"CAS_NO"]<- newId
		cumInchi<- repackData[,"InChI"]
		}
	if(!reduce)
		{duplInchi<- rep("", length(cumInchi))
		for (s in 1: length(cumInchi))
			{#s<- 1; s<- s+1
			currInch<- cumInchi[s]
			duplTest<-  currInch==cumInchi
			duplInx<-  which(currInch==cumInchi)
			if(length(which(duplInchi=="" & duplTest))>0)
				{duplInchi[which(duplInchi=="" & duplTest)]<- repackData[duplInx[1], "InChI"] }
			}
		repackData[, "InChI"] <- duplInchi
		}
	appDupl<- missData
	if(nrow(appDupl)>0){appDupl[,2]<- rep("inChI_failed", nrow(appDupl))}
	if(nrow(appDupl)==0){appDupl<- matrix("", nrow=0, ncol=2)}
 	if(reduce) 
		{uInChi<- unique(repackData[, "InChI"])
		appDupl2<- matrix(nrow=0, ncol=2)
		colnames(appDupl2)<- c("ID", "ass_ID")
		cumUniqInx<- c()
		for (u in 1: length(uInChi)) #4301
			{#u<-4301 1; u<- 1; u<- +1
			selUinx<- which(repackData[, "InChI"]==uInChi[u])
			sourceSet<- repackData[selUinx, "CAS_NO"]
			appSourceSet<- gsub("PubChem", "", sourceSet, fixed=TRUE)
			pubChemSubInx<- which(nchar(appSourceSet)<nchar(sourceSet))
			firSourInx<- selUinx[1]
			if(length(pubChemSubInx)>0)
				{firSourInx<- selUinx[pubChemSubInx[1]]}
			cumUniqInx<- c(cumUniqInx, firSourInx)
			if(length(selUinx)>1)
				{appenDupl<- unique(repackData[setdiff(selUinx, firSourInx), "CAS_NO"])
				ass_ID<- rep(repackData[firSourInx,"CAS_NO"], length(appenDupl))
				appDuplX<- cbind(appenDupl, ass_ID)
				appDupl2<- rbind(appDupl2, appDuplX)
				}
			}
		appDupl<- rbind(appDupl, appDupl2)
		repackData<- repackData[cumUniqInx,]
		rownames(appDupl)<- NULL
		}

	resLs[[1]]<- repackData[!duplicated(repackData[,"CAS_NO"]),]
	resLs[[2]]<- appDupl[!duplicated(appDupl[,"ID"]),]
	return(resLs)
	}


#ls<- currMolec
#currPart<- "C10H13N2O8.3K/h13H;;;/q-3;3m"
#currPart<- "2C6H11O7.Zn/q2*-1;m"
#currPart<- "C15H34N.Br/h;1h/qm;-1"

abbrevFinder<- function(currPart)
	{abbrevSet<- c("m", "M", "n", "N")#, "i"
	# Different letters are used to refer to different locations of the first instance of the same layer information:
	# m – item in the first section of the Identifier, but not in the isotopic segment
	# M – item in the isotopic part of the first section
	# n – item in the fixed-H section, but not in the isotopic segment
	# N – item in the isotopic part of fixed-H section
	# i – prefix to m, M, n, or N – indicates that sp3-stereo has been inverted.

	infoLs<- list("testSplit"=c(), "foundAbbrev"=c(), "nameOfLayer"=c(), "contentLocation"=c(), "toChange"=c(), "locatOfComp"=c())
	currHeadPart<- paste("/", unlist(strsplit(currPart, "/", fixed=TRUE))[1],  collapse="", sep="")
	currPart<- paste("/", unlist(strsplit(currPart, "/", fixed=TRUE))[-1],  collapse="", sep="")
	testSplit<- FALSE
	maxAbbrev<- FALSE
	cnt<- 0
	foundAbb<- ""
	while(!testSplit & !maxAbbrev)
		{cnt<- cnt+1
		currAbbrev<- abbrevSet[cnt]
		mayBeSplit<- unlist(strsplit(currPart, currAbbrev, fixed=TRUE))
		testSplit<- (length(mayBeSplit)>1 | (length(mayBeSplit)==1 & nchar(mayBeSplit[1])<nchar(currPart)))
		checkSplit<- substr(mayBeSplit[1], nchar(mayBeSplit[1]), nchar(mayBeSplit[1]))!="/"
		if(!checkSplit){testSplit<- FALSE}
		if(testSplit){foundAbb<- currAbbrev}
		if(cnt==length(abbrevSet)){maxAbbrev<- TRUE}
		}
	if(testSplit)	
		{layers<- unlist(strsplit(mayBeSplit[1], "/", fixed=TRUE))
		nameOfLay<- substr(layers[length(layers)], 1, 1)
		infoLs$nameOfLayer<- nameOfLay
		if(length(mayBeSplit)==1){mayBeSplit<- c(mayBeSplit, "")}
		mayBeSplit<- c(paste(currHeadPart, mayBeSplit[1], collapse="", sep=""),  mayBeSplit[2])
		infoLs$toChange<- mayBeSplit
		appJustLayer<- unlist(strsplit(mayBeSplit[1], nameOfLay, fixed=TRUE))
		if(length(appJustLayer)==1 & nchar(appJustLayer[1])<nchar(mayBeSplit[1])) {appJustLayer<- c(appJustLayer, "")}
		justLayer<- appJustLayer[2]
		appJust<- unlist(strsplit(justLayer, ";", fixed=TRUE))
		if(length(appJust)==0) {locOfComp<- 1}
		if(length(appJust)>0)
			{if(length(appJust)==1 & nchar(appJust[1])<nchar(justLayer)){justLayer<- c(justLayer, "")}
			appMultipl<-  unlist(strsplit(appJust, "*", fixed=TRUE))
			fct<- 1
			if(length(appMultipl)>1)
				{multipl<- suppressWarnings(as.numeric(appMultipl[1]))
				if(is.numeric(multipl)){fct<- multipl}
				}
			locOfComp<-  length(justLayer)+ fct-1
			}
		infoLs$locatOfComp<- locOfComp
		}
	infoLs$testSplit<- testSplit
	infoLs$foundAbbrev<- foundAbb
	if(foundAbb!="")
		{if(foundAbb=="m")
			{conLoc<- c("disconnected", "main", "nonIsotopic")}
		if(foundAbb=="M")
			{conLoc<- c("disconnected", "main", "i_sotopic")}
		if(foundAbb=="n")
			{conLoc<- c("disconnected", "f_ixed", "nonIsotopic")}
		if(foundAbb=="N")
			{conLoc<- c("disconnected", "f_ixed", "i_sotopic")}
		#if(foundAbb=="i")
		#	{conLoc<- ""}
		infoLs$contentLocation<- conLoc
		}
	return(infoLs)
	}

#inChIStr_a<- inChIStr
#inChIStr<- inChIStr_a
abbrevSubst<- function(inChIStr)
	{for (i in 1: length(inChIStr))
		{#i<- 1; i<- i+1
		if(length(inChIStr[[i]])>0)
			{if(is.list(inChIStr[[i]]))
				{for (j in 1: length(inChIStr[[i]]))
					{#j<- 1; j<- j+1
					if(length(inChIStr[[i]][[j]])>0)
						{if(is.list(inChIStr[[i]][[j]]))
							{for (k in 1: length(inChIStr[[i]][[j]]))
								{#k<- 1; k<- k+1
								if(length(inChIStr[[i]][[j]][[k]])>0)
									{if(is.list(inChIStr[[i]][[j]][[k]]))
										{for (d in 1: length(inChIStr[[i]][[j]][[k]]))
											{#d<- 1; d<- d+1
											#inChIStr[[i]][[j]][[k]][[d]]<- list()
											currPart<- inChIStr[[i]][[j]][[k]][[d]]
											foundInfo<- abbrevFinder(currPart)
											testSplit<- foundInfo[[1]]	
											if(testSplit)
												{currAbbrev<- foundInfo[[2]]
												nameOfLayer<- foundInfo[[3]]												
												contentLocation<- foundInfo[[4]]
												mayBeSplit<- foundInfo[[5]]
												locOfComp<- foundInfo[[6]]
												if(contentLocation[1]=="disconnected"){fir<- 1}
												if(contentLocation[1]!="disconnected"){fir<- 2}
												if(contentLocation[2]=="main"){sec<- 1}
												if(contentLocation[2]!="main"){sec<- 2}
												if(contentLocation[3]=="nonIsotopic"){ter<- 1}
												if(contentLocation[3]!="nonIsotopic"){ter<- 2}	
												prevLayers<- unlist(strsplit(inChIStr[[fir]][[sec]][[ter]], "/", fixed=TRUE))
												for (g in 1: length(prevLayers))
													{#g<-4 1; g<- g+1
													findPrevLayer<-  unlist(strsplit(prevLayers[g], nameOfLayer, fixed=TRUE))
													if(nchar(findPrevLayer[1])<nchar(prevLayers[g]))
														{contentOfLayer<- findPrevLayer[2]
														contentOfLayer<- unlist(strsplit(contentOfLayer, ";", fixed=TRUE))
														contentOfLayer<- contentOfLayer[locOfComp]
														newString<- paste(mayBeSplit[1], contentOfLayer, mayBeSplit[2], sep="", collapse="")
														}
													}
												inChIStr[[i]][[j]][[k]][[d]]<- newString
												}
											}
										}
									if(!is.list(inChIStr[[i]][[j]][[k]]))
										{#inChIStr[[i]][[j]][[k]]<- list()
										currPart<- inChIStr[[i]][[j]][[k]]
										foundInfo<- abbrevFinder(currPart)
										testSplit<- foundInfo[[1]]
										newString<- currPart	
										#ciao<- se							
										if(testSplit)
											{currAbbrev<- foundInfo[[2]]
											nameOfLayer<- foundInfo[[3]]											
											contentLocation<- foundInfo[[4]]
											mayBeSplit<- foundInfo[[5]]
											locOfComp<- foundInfo[[6]]
											if(contentLocation[1]=="disconnected"){fir<- 1}
											if(contentLocation[1]!="disconnected"){fir<- 2}
											if(contentLocation[2]=="main"){sec<- 1}
											if(contentLocation[2]!="main"){sec<- 2}
											if(contentLocation[3]=="nonIsotopic"){ter<- 1}
											if(contentLocation[3]!="nonIsotopic"){ter<- 2}	
											prevLayers<- unlist(strsplit(inChIStr[[fir]][[sec]][[ter]], "/", fixed=TRUE))
											prevLayers<- prevLayers[which(prevLayers!="" & !is.na(prevLayers))]
											for (g in 1: length(prevLayers))
												{#g<- 1; g<- g+1
												findPrevLayer<-  unlist(strsplit(prevLayers[g], nameOfLayer, fixed=TRUE))
												if(nchar(findPrevLayer[1])<nchar(prevLayers[g]))
													{#if(j==2) {alt<- qu}
													contentOfLayer<- findPrevLayer[2]
													contentOfLayer<- unlist(strsplit(contentOfLayer, ";", fixed=TRUE))
													contentOfLayer<- contentOfLayer[locOfComp]
													newString<- paste(mayBeSplit[1], contentOfLayer, mayBeSplit[2], sep="", collapse="")
													}
												}
											inChIStr[[i]][[j]][[k]]<- newString
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	return(inChIStr)
	}

emptyListCmp<- function(ls)
	{for (i in 1: length(ls))
		{#i<- 1; i<- i+1		
		if(length(ls[[i]])>0)
			{if(is.list(ls[[i]]))
				{for (j in 1: length(ls[[i]]))
					{#j<- 1; j<- j+1
					if(length(ls[[i]][[j]])>0)
						{if(is.list(ls[[i]][[j]]))
							{for (k in 1: length(ls[[i]][[j]]))
								{#k<- 1; k<- k+1
								if(length(ls[[i]][[j]][[k]])>0)
									{if(is.list(ls[[i]][[j]][[k]]))
										{for (d in 1: length(ls[[i]][[j]][[k]]))
											{#d<- 1; d<- d+1
											ls[[i]][[j]][[k]][[d]]<- list()
											}
										}
									if(!is.list(ls[[i]][[j]][[k]]))
										{ls[[i]][[j]][[k]]<- list()}
									}
								}
							}
						}
					}
				}
			}
		}
	return(ls)
	}

#numCmp<- g
#ls<- currLs[[g]]
#contLs<- currMolec

fillListCmp<- function(numCmp, ls, contLs)
	{newLs<- ls
	for (i in 1: length(ls))
		{#i<- 1; i<- i+1		
		if(length(ls[[i]])>0)
			{if(is.list(ls[[i]]))
				{for (j in 1: length(ls[[i]]))
					{#j<- 1; j<- j+1
					if(length(ls[[i]][[j]])>0)
						{if(is.list(ls[[i]][[j]]))
							{for (k in 1: length(ls[[i]][[j]]))
								{#k<- 1; k<- k+1
								if(length(ls[[i]][[j]][[k]])>0)
									{if(is.list(ls[[i]][[j]][[k]]))
										{for (d in 1: length(ls[[i]][[j]][[k]]))
											{#d<- 1; d<- d+1
											appCont<- unlist(contLs[[i]][[j]][[k]][[d]])
											if(appCont!="")
												{newLs[[i]][[j]][[k]][[d]]<- distribComp(appCont, numComp)}
											if(appCont=="")
												{newLs[[i]][[j]][[k]][[d]]<- NULL}
											}
										}
									if(!is.list(ls[[i]][[j]][[k]]))
										{appCont<- unlist(contLs[[i]][[j]][[k]])
										if(appCont!="")
											{newLs[[i]][[j]][[k]]<- distribComp(appCont, numCmp)}
										if(appCont=="")
											{newLs[[i]][[j]][[k]]<- NULL}
										}
									}
								if(length(ls[[i]][[j]][[k]])==0)
									{appCont<-  unlist(contLs[[i]][[j]][[k]])
									if(appCont!="")
										{newLs[[i]][[j]][[k]]<- distribComp(appCont, numCmp)}
									if(appCont=="")
										{newLs[[i]][[j]][[k]]<- NULL}
									}
								}
							}
						if(length(ls[[i]][[j]])==0)
							{appCont<- unlist(contLs[[i]][[j]])
							if(appCont!="")
								{newLs[[i]][[j]]<- distribComp(appCont, numCmp)}
							if(appCont=="")
								{newLs[[i]][[j]]<- NULL}
							}
						}
					}
				}
			}
		if(length(ls[[i]])==0)
			{appCont<-  unlist(contLs[[i]])
			if(appCont!="")
				{newLs[[i]]<- distribComp(appCont, numCmp)}
			if(appCont=="")
				{newLs[[i]]<- NULL}
			}
		}		
	return(newLs)
	}

#fullStruct<- inchiDataLs
#collect_CAS[i]
singleCompInChI<- function(fullStruct)
	{collect_CAS<- names(fullStruct)
	for (i in 1: length(fullStruct))
		{#i<- 1; i<- i+1
		currMolec<- fullStruct[[i]]
		numComp<- howManyComp(currMolec)
		currLs<- list()
		for (g in 1: numComp)
			{#g<-1; g<- g+1
			emptLs<- emptyListCmp(currMolec)
			currLs[[g]]<- emptLs
			}
		names(currLs)<- paste("comp_", seq(1, numComp), sep="" )
		currSingLs<- currLs
		for (g in 1: numComp)
			{#g<-1; g<- g+1
			currSingLs[[g]]<- fillListCmp(g, currLs[[g]], currMolec)			
			}
		fullStruct[[i]]<- currSingLs		
		}
	return(fullStruct)
	}

#appCont<- "C19H39O7S.Na/q-1"
#numComp<- numCmp
x_distribComp<- function(appCont, numComp)
	{dcmpFormul <- unlist(strsplit(appCont, ".", fixed=TRUE))
	cmpSize<- length(dcmpFormul)
	dcmParts <- unlist(strsplit(appCont, "/", fixed=TRUE))
	relevParts<- dcmParts[-1]
	relevPartsSize<- length(relevParts)
	x_dcmpFormul<- dcmParts[1]
	f_dcmpFormul<- unlist(strsplit(x_dcmpFormul, ".", fixed=TRUE))
	dcmpVect<- rep("", cmpSize)
	for (q in 1: length(dcmpFormul))
		{#q<- 1; q<- q+1
		currPiece<- f_dcmpFormul[q]		
		synt<- c()
		for (h in 1: relevPartsSize)
			{#h<- 1; h<- h+1
			currStrata<- relevParts[h]
			cmpOfStrata<- unlist(strsplit(currStrata, ";", fixed=TRUE))
			foundCmp<- FALSE
			for (z in 1: length(cmpOfStrata))
				{#z<- 1; z<- z+1
				if(z==q)
					{synt<- paste(synt, currPiece, "/", cmpOfStrata[z], sep="", collapse="") 
					if(!foundCmp){foundCmp<- TRUE}
					}
				}
			if(!foundCmp) {synt<- f_dcmpFormul[q]}
			}		
		dcmpVect[q]<- synt
		}
	return(dcmpVect[numComp])
	}

#stratum<- currStrata
nameOfStratum<- function(stratum)
	{firLett<-  substr(stratum,1,1)
	testUpperCase<- firLett==toupper(firLett)
	if(testUpperCase)
		{nameOfStra<- ""}
	if(!testUpperCase)
		{if(firLett=="s")
			{secLett<-  substr(stratum,2,2)
			if(secLett=="1" | secLett=="2" | secLett=="3")
				{nameOfStra<- paste("/", firLett, secLett, sep="", collapse="")}
			}
		if(firLett!="s")
			{nameOfStra<- paste("/", firLett, sep="", collapse="")
			}
		}
	return(nameOfStra)
	}
 
#numComp<- numCmp
distribComp<- function(appCont, numComp)
	{dcpStrata<- unlist(strsplit(appCont, "/", fixed=TRUE))
	if(length(dcpStrata)==1 & nchar(dcpStrata[1])<nchar(appCont))
		{dcpStrata<- c(dcpStrata, "")}
	numStrata<- length(dcpStrata)
	collInfo<- c()
	preSize<- 0
	for (s in 1: numStrata)
		{#s<- 1; s<- s+1
		currStrata<- dcpStrata[s]
		currName<- nameOfStratum(currStrata)
		if(currName=="")
			{compOfStrata<- unlist(strsplit(currStrata, ".", fixed=TRUE))
			dotSep<- TRUE
			preSize<- length(compOfStrata)
			}
		if(currName!="")
			{compOfStrata<- unlist(strsplit(currStrata, ".", fixed=TRUE))
			dotSep<- TRUE
			if(nchar(compOfStrata[1])<nchar(currStrata))
				{compOfStrata<- c(compOfStrata[1], paste(rep("", preSize-1), sep=""))}
			if(length(compOfStrata)==1)
				{compOfStrata<- unlist(strsplit(currStrata, ";", fixed=TRUE))
				dotSep<- FALSE}
			if((currName=="s1" | currName=="s2" | currName=="s3") & length(compOfStrata)==1)
				compOfStrata<- rep(compOfStrata, preSize)
			}
		cmpSize<-  length(compOfStrata)
		if(cmpSize>0)
			{for (q in 1: cmpSize)
				{#q<- 1; q<- q+1
				currQ<- compOfStrata[q]
				if(currName!="")
					{if(currName=="/s1" | currName=="/s2" | currName=="/s3")
						{if(nchar(currQ)==2) {currQ<-""}
						else {currQ<- substr(currQ, 3, nchar(currQ))}
						}
					else {if(nchar(currQ)==1) {currQ<-""}
						else {currQ<- substr(currQ, 2, nchar(currQ))}
						}
					}
				if(q==numComp)
					{appString<- paste(currName, currQ, collapse="", sep="")
					collInfo<- paste(collInfo, appString, collapse="", sep="")
					}
				}
			}
		if(cmpSize==0)
			{collInfo<- paste(collInfo, "", collapse="", sep="") }
		}
	return(collInfo)
	} 
