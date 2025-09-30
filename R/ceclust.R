#' Clustering par Entropie Composite sur données mixtes
#'
#' @description
#' Apprend un modèle de mélange (gaussien pour les variables numériques et
#' multinomial pour les variables catégorielles) et renvoie l’assignation de
#' cluster ainsi que les paramètres estimés.
#'
#' @param Z `data.frame` ou `matrix` d'observations (colonnes = variables).
#'   Les colonnes numériques sont traitées via une densité gaussienne, les
#'   colonnes factor via des probabilités discrètes.
#' @param lambda `numeric` (>0). Paramètre d’entropie composite.
#' @param C `numeric` (>0). Constante de borne pour stabiliser les densités.
#' @param r0 `integer` ou `NULL`. Nombre de classes initial si fourni.
#' @param Nshots `integer`. Nombre de redémarrages aléatoires (multi-start).
#' @param Nloop `integer`. Max itérations internes par redémarrage.
#' @param familyType `character`. Famille cible (`"gaussAndDiscreteVector"`, etc.).
#' @param sizeMaxOutlier `integer`. Seuil (taille min) d’un groupe considéré
#'   outlier à regrouper (si > 0).
#' @param autoRegroupOutliers `logical`. Regroupement auto des petits groupes.
#' @param displayRemainingTime `logical`. Affiche le temps estimé.
#' @param focus `integer` ou `character` ou `NULL`. Si une colonne factor est
#'   fournie, le clustering est appris séparément par niveau puis agrégé.
#'
#' @return Une `list` avec au minimum :
#' \itemize{
#'   \item `phi` : vecteur d'appartenance (one-hot) de longueur `n*r`.
#'   \item `params` : liste des paramètres estimés (nu, m, Sigma, etc.).
#'   \item `Hphi` : valeur de l’objectif (entropie composite).
#'   \item `clusters` : vecteur des clusters prédits (longueur `n`).
#'   \item `REO` : nombre final de classes retenues.
#' }
#'
#' @seealso [CECclassifNewData()], [CECpredict()]
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' # Deux clusters gaussiens sur 2 variables
#' x1 <- rnorm(n, 0, 1); x2 <- rnorm(n, 0, 1)
#' x3 <- rnorm(n, 3, 1); x4 <- rnorm(n, 3, 1)
#' X  <- rbind(cbind(x1, x2), cbind(x3, x4))
#' # Une variable catégorielle corrélée au cluster
#' lab <- factor(rep(c("A","B"), each = n))
#' Z   <- data.frame(X1 = X[,1], X2 = X[,2], Cat = lab)
#'
#' fit <- CECclassif(Z, Nshots = 10, Nloop = 200)
#' table(fit$clusters, lab)
#'
#' @export

CECclassif				<- function(Z,lambda=1,C=1,r0=NULL,Nshots = 100,Nloop=1000,familyType="gaussAndDiscreteVector" ,sizeMaxOutlier = 0,autoRegroupOutliers=FALSE,displayRemainingTime = FALSE,focus=NULL)
{
	displayPlotEntropy = FALSE
	if(!is.null(focus))
	{
		if(is.factor(Z[,focus]))
		{
			focus_fac 	<- unique(Z[,focus])
			results   	<- list()
			NU 			<- c()
			R 			<- c()
			H 			<- c()
			for(f in 1:length(focus_fac))
			{
				fac_f 			<- focus_fac[f]
				pos_f 			<- which(Z[,focus]==fac_f)
				results[[f]] 	<- CECclassif(Z[pos_f,],lambda=lambda,C=C,r0=r0,Nloop=Nloop,Nshots = Nshots,familyType=familyType,sizeMaxOutlier =sizeMaxOutlier,displayRemainingTime = displayRemainingTime,focus=NULL)
				NU[f]			<- length(pos_f)
				R[f] 			<- length(results[[f]]$params$states)
				H[f] 			<- results[[f]]$Hphi	
			}
			
			NU <- NU/sum(NU)
			
			results_tot 				<- list()
			results_tot$params 			<- list()
			
			
			r_tot 		<- sum(R)
			states_tot  <- 1:r_tot				
			results_tot$params$states <- 1:r_tot 
			
			results_tot$params$familyType <-familyType
			
			
			
		
			results_tot$params$factors 	<- results[[1]]$params$factors 
			results_tot$params$lambda 	<- lambda
			results_tot$params$C 		<- C
			
			if(length(results[[1]]$params$discreteProbList)>0)
				results_tot$params$discreteProbList <- list()
			
			if(length(results[[1]]$params$m)>0)
			{
				results_tot$params$m <- matrix(NA,r_tot,dim(results[[f]]$params$m)[2])
				colnames(results_tot$params$m) <- colnames(results[[1]]$params$m)
			}
				
			if(length(results[[1]]$params$Sigma)>0)
				results_tot$params$Sigma <- list()
				
			if(length(results[[1]]$params$colFactor)>0)
				results_tot$params$colFactor <- results[[1]]$params$colFactor
				
			if(length(results[[1]]$params$colNum)>0)
				results_tot$params$colNum <- results[[1]]$params$colNum
	 
	 

			
			nu_tot   	<- rep(0,r_tot)
			
			phiM_tot 	<-  matrix(0,dim(Z)[1],r_tot)
			
			rind    <-  1
			
			for(f in 1:length(focus_fac))
			{
				fac_f 			<- focus_fac[f]
			
				pos_f 			<- which(Z[,focus]==fac_f)
			
				phi_f			<- results[[f]]$phi
				phiM_f 			<- vectToMat(phi_f,R[f])
				
				phiM_tot[pos_f,rind:(rind+R[f]-1)] <-  phiM_f
				
				
				
				nu_tot[rind:(rind+R[f]-1)] = NU[f]*results[[f]]$params$nu
				
				if(length(results[[f]]$params$discreteProbList)>0)
				{
					for(k in 1:R[f])
					{
						results_tot$params$discreteProbList[[rind+k-1]] <- results[[f]]$params$discreteProbList[[k]]
					}
				}
				
				if(length(results[[f]]$params$m)>0)
				{
					results_tot$params$m[rind:(rind+R[f]-1),] <- results[[f]]$params$m
				}
				
				if(length(results[[f]]$params$Sigma)>0)
				{
					for(k in 1:R[f])
					{
						results_tot$params$Sigma[[rind+k-1]] <- results[[f]]$params$Sigma[[k]]
					}
				}
				
				
				rind <- rind +R[f]
			}
			
			phi_tot 				<- matToVect(phiM_tot)
			results_tot$phi 		<- phi_tot
			results_tot$params$phi 	<- phi_tot
			results_tot$params$nu	<- nu_tot
			
			
			Hphi_tot <- evalCompositeEntropy(phi=phi_tot ,Z=Z,lambda=lambda,C=C,familyType=familyType)
			
			results_tot$Hphi <- Hphi_tot
			
			return(results_tot)
		}				
	}
	
	n <- length(Z)
	if(is.matrix(Z)|is.data.frame(Z) )
	{
		n <- dim(Z)[1] 
	}
		
	
	PHI 		<- c() 
	H 			<- Inf
	bestClassif <- NULL
	
	T <- Sys.time()
	for(k in 1:Nshots)
	{
		if(k < 2)
		{
			phi0 <- NULL 
		}else{
			U <- runif(1)
			if(U<0.25)
			{
				phi0 <- c(PHI*runif(length(PHI)) ,runif(n) )
			}else if(U>=0.25 & U<0.5){
				if(is.null(r0) | runif(1)<0.5)
				{
					rr <- 1 
				}else{
					rr <- sample(r0,1)
				}
				phi0 <- runif(rr*n)
			}else{
				phi0 <- NULL 
			}
			
			if(k %in%((Nshots- floor(Nshots)/10 ):Nshots) & H<Inf & n >100)
			{
			 
				r 			<- length(bestClassif$params$states)
				nu 			<- bestClassif$params$nu
				phi 		<- bestClassif$phi 
			
				phiM 		<- vectToMat(phi,r)
			
				toKeep      <- which(nu*n>=10)
				phiM 		<-  phiM[,toKeep,drop=FALSE]
				
				tryCatch({
				  phiM 		<- 	phiM + matrix(0.1*runif(prod(dim(phiM))),dim(phiM)[1],dim(phiM)[2])
				}, error = function(e) {
				  message("error phi M: ", phiM)
				})

				
				sumColsPhi 	<- apply(phiM,1,sum) 
				phiM		<- sweep(phiM, 1, sumColsPhi, "/")
				phi0        <- matToVect(phiM)
			}
			
		}
		
		resultList 	<- CECclassifOneShot(Z=Z,lambda=lambda,C=C,r0=r0,Nloop=Nloop,phi0=phi0,familyType=familyType,displayPlotEntropy=displayPlotEntropy)
		
		if(sizeMaxOutlier<n | autoRegroupOutliers)
		{
			# autoRegroupOutliers écrase sizeMaxOutlier
			
			continueOutlier = TRUE 
			while(continueOutlier){
			
				sizeGroups 			<- resultList$params$nu*n
				if(autoRegroupOutliers)
				{
					groupRank 				<- rank(sizeGroups,ties.method="min")
					functionBoundReached 	<- resultList$params$functionBoundReached
					
					candidateToRegroup		<- which(groupRank<=floor(length(groupRank)/2))
					 
					
					candidateToRegroup		<- candidateToRegroup[which(candidateToRegroup%in%functionBoundReached )]
					outliersToRegroup  		<-  c()
					if(length(candidateToRegroup)>0  &  length(sizeGroups)>1)
					{
						outliersToRegroup  	<- candidateToRegroup[which.min(groupRank[candidateToRegroup])]
					}
					 
					
				}else{
					outliersToRegroup  	<- which(sizeGroups<sizeMaxOutlier)
				}
				
			
				if(length(outliersToRegroup)>0 )
				{
					r 							<- length(sizeGroups)
					phi 						<- resultList$phi
					phiM 						<- vectToMat(phi,r)
					if(length(outliersToRegroup)<r)
					{
						phiM0  						<- phiM[,-outliersToRegroup,drop=FALSE]
					}else{
						phiM0  						<- phiM[,1,drop=FALSE]	
					}
						
					
					dataToReassign 				<- which(apply(phiM[,outliersToRegroup,drop=FALSE],1,sum)>0)
					phiM0[dataToReassign,] 		<- runif(length(dataToReassign)*dim(phiM0)[2])
					sumColsphiM0Reassign  		<- apply(phiM0[dataToReassign,,drop=FALSE],1,sum) 
					
					phiM0[dataToReassign,] 		<- sweep(phiM0[dataToReassign,,drop=FALSE] , 1, sumColsphiM0Reassign, "/")
					phi0						<- matToVect(phiM0)
					resultList 					<- CECclassifOneShot(Z=Z,lambda=lambda,C=C,r0=r0,Nloop=Nloop,phi0=phi0,familyType=familyType,displayPlotEntropy=displayPlotEntropy)
			
				}else{
				
					continueOutlier = FALSE
				}
			}
			
		}
		
		
		if(resultList$Hphi <H )
		{
			H 			<- resultList$Hphi 
			bestClassif <- resultList
			if(displayRemainingTime)
				print(paste("newBest shot :",k, "H = ", round(H,2)))
		}
		
		#PHI			<- c(PHI,resultList$phi)
		PHI				<- resultList$phi
		if(displayRemainingTime)
			remainingTime(T,k,Nshots)
		
	}
	
	REO 				<- length(bestClassif$params$states)
	bestClassif$REO 	<- REO
	phiM 				<- vectToMat(bestClassif$phi,REO)
	clusters 			<- apply(phiM,1,which.max)
	
	bestClassif$clusters 	<- clusters
	return(bestClassif)
	
}


#' Affectation de nouveaux individus à des clusters existants
#'
#' @description
#' À partir d’un objet `params` appris par [CECclassif()], assigne à chaque
#' ligne de `Zpred` le cluster le plus probable.
#'
#' @param Zpred `data.frame`/`matrix` de nouvelles observations
#'   (même schéma de variables que `Z` lors de l’apprentissage, sauf
#'   éventuellement certaines colonnes à prédire).
#' @param params `list`. Paramètres du modèle retournés par [CECclassif()].
#' @param idColToPred `integer` ou `character`. Indices/noms des colonnes
#'   considérées comme à prédire (donc ignorées pour l’affectation).
#'
#' @return `integer` vector des indices de cluster (longueur `nrow(Zpred)`).
#'
#' @seealso [CECclassif()], [CECpredict()]
#'
#' @examples
#' # suite de l'exemple CECclassif()
#' Znew <- Z[1:5, ]
#' Znew$Cat <- factor(NA, levels = levels(Z$Cat))  # ex: colonne catégorielle manquante
#' pred_clusters <- CECclassifNewData(Znew, fit$params, idColToPred = "Cat")
#' pred_clusters
#'
#' @export

CECclassifNewData 		<- function(Zpred,params,idColToPred)
{
	if(params$familyType	== "gaussAndDiscreteVector")
	{ 
		r 			<-  length(params$states)
		colFactor	<- 	params$colFactor 
		colNum		<- 	params$colNum
		
		lambda 		<- params$lambda
		
		colFactorMissing		<- which(colFactor%in%idColToPred)
		colNumMissing 			<- which(colNum%in%idColToPred)
	
		paramsPred 				<- params
		
		paramsPred$m  			<- params$m[,-colNumMissing,drop=FALSE]
		if(length(colNumMissing)==0)
			paramsPred$m  			<- params$m
		
			
		colFactorPred 			<- colFactor[-colFactorMissing]
		if(length(colFactorMissing)==0)
			colFactorPred 			<- colFactor
		
		
		colNumPred 				<- colNum[-colNumMissing]
		if(length(colNumMissing)==0)
			colNumPred 			<- colNum
			
		regroupColNum 	<- rank(c(colFactorPred,colNumPred ))
		
		if(length(colFactorPred)>0)
		{
			colFactorPred 	<- regroupColNum[1:length(colFactorPred)]
			
		}else{
			colFactorPred 	<- c()
		}
		
		if(length(regroupColNum)>length(colFactorPred))
		{
			colNumPred 		<- regroupColNum[(1+length(colFactorPred)):length(regroupColNum)]
		}else{
			colNumPred 		<- c()
		}
		
		
		
		
		paramsPred$colFactor	<- colFactorPred
		paramsPred$colNum		<- colNumPred
		
		
		
		for(x in 1:r)
		{	
			paramsPred$discreteProbList[[x]] <- params$discreteProbList[[x]][,-colFactorMissing,drop=FALSE]
			paramsPred$Sigma[[x]]			 <- params$Sigma[[x]][,-colNumMissing,drop=FALSE]
			paramsPred$Sigma[[x]]			 <- paramsPred$Sigma[[x]][-colNumMissing,,drop=FALSE]
			if(length(colFactorMissing)==0)
				paramsPred$discreteProbList[[x]] <- params$discreteProbList[[x]]
			 
			if(length(colNumMissing)==0)
			{
				paramsPred$Sigma[[x]]			 <- params$Sigma[[x]] 
				paramsPred$Sigma[[x]]			 <- paramsPred$Sigma[[x]] 
			}
		}
	
		phiPred <- optPhi(Z=Zpred,param = paramsPred,lambda=lambda)
		phiPredM<- vectToMat(phiPred,r)
		
		XpredHat 	<- apply(phiPredM,1,which.max)
		
		
		
		
		
		return(XpredHat)
	}

}

#' Prédiction de colonnes manquantes conditionnellement au cluster
#'
#' @description
#' En utilisant les paramètres `params` appris par [CECclassif()], prédit les
#' valeurs des colonnes `idColToPred` pour `Zpred`. Les colonnes factor sont
#' prédites par l’argmax des probabilités discrètes; les colonnes numériques
#' via la moyenne conditionnelle gaussienne.
#'
#' @param Zpred `data.frame`/`matrix` de nouvelles observations.
#' @param params `list`. Paramètres du modèle retournés par [CECclassif()].
#' @param idColToPred `integer` ou `character`. Colonnes à prédire.
#'
#' @return `data.frame` de taille `nrow(Zpred) × length(idColToPred)` avec les
#'   valeurs prédites (facteurs et/ou numériques).
#'
#' @seealso [CECclassif()], [CECclassifNewData()]
#'
#' @examples
#' # suite de l'exemple CECclassif()
#' Zmiss <- Z[1:5, ]
#' Zmiss$X2 <- NA_real_
#' preds <- CECpredict(Zmiss, fit$params, idColToPred = "X2")
#' preds
#'
#' @export

CECpredict				<- function(Zpred,params,idColToPred)
{
	if(params$familyType	== "gaussAndDiscreteVector")
	{ 
		r 			<-  length(params$states)
		colFactor	<- 	params$colFactor 
		colNum		<- 	params$colNum
		
		lambda 		<- params$lambda
		
		colFactorMissing		<- which(colFactor%in%idColToPred)
		colNumMissing 			<- which(colNum%in%idColToPred)
	
		paramsPred 				<- params
		
		paramsPred$m  			<- params$m[,-colNumMissing,drop=FALSE]
		if(length(colNumMissing)==0)
			paramsPred$m  			<- params$m
		
			
		colFactorPred 			<- colFactor[-colFactorMissing]
		if(length(colFactorMissing)==0)
			colFactorPred 			<- colFactor
		
		
		colNumPred 				<- colNum[-colNumMissing]
		if(length(colNumMissing)==0)
			colNumPred 			<- colNum
			
		regroupCols		<- rank(c(colFactorPred,colNumPred ))
	 
		if(length(colFactorPred)>0)
			colFactorPred 	<- regroupCols[1:length(colFactorPred)]
		
		 
		colNumPred 		<- regroupCols[(1+length(colFactorPred)):length(regroupCols)]
		
		paramsPred$colFactor	<- colFactorPred
		paramsPred$colNum		<- colNumPred
		
		
		
		for(x in 1:r)
		{	
			paramsPred$discreteProbList[[x]] <- params$discreteProbList[[x]][,-colFactorMissing,drop=FALSE]
			paramsPred$Sigma[[x]]			 <- params$Sigma[[x]][,-colNumMissing,drop=FALSE]
			paramsPred$Sigma[[x]]			 <- paramsPred$Sigma[[x]][-colNumMissing,,drop=FALSE]
			if(length(colFactorMissing)==0)
				paramsPred$discreteProbList[[x]] <- params$discreteProbList[[x]]
			 
			if(length(colNumMissing)==0)
			{
				paramsPred$Sigma[[x]]			 <- params$Sigma[[x]] 
				paramsPred$Sigma[[x]]			 <- paramsPred$Sigma[[x]] 
			}
		}
		
		phiPred <- optPhi(Z=Zpred,param = paramsPred,lambda=lambda)
		phiPredM<- vectToMat(phiPred,r)
		
		XpredHat 	<- apply(phiPredM,1,which.max)
		ZpredVals	<-data.frame(matrix(NA,dim(Zpred)[1],length(idColToPred)))
		
		if(length(colFactorMissing)>0)
		{
			factors <- params$factors
			predFact = as.data.frame(matrix(NA,dim(Zpred)[1],length(colFactorMissing)))
			for(j in 1:length(colFactorMissing))
			{
				predFact[,j] = factor(predFact[,j] ,levels=factors)
			}
			
			for(x in 1:r)
			{
				fact_x <- factors[apply(params$discreteProbList[[x]][,colFactorMissing,drop=FALSE],2,which.max)]
				for(j in 1:length(colFactorMissing))
				{
					predFact[XpredHat==x,j] <- fact_x[j]
				}						
			}
			
			
			
			ZpredVals[,which(idColToPred%in%colFactor)] <- predFact
			
		}
		
		
		if(length(colNumMissing)>0)
		{
			if(length(colNum)>length(colNumMissing))
			{
				# si j'ai un m et un Sigma pour Z1, Z2, Z3
				# je veux prédire Z3 sachant Z1 Z2 
				# f(z1,z2,z3) = 1/sqrt(2*pi)^3  . 1/sqrt(det(Sigma)) exp( - 1/2  ((z1,z2,23) - m)^t Sigma^(-1) ((z1,z2,23) - m) ) 
				# f(z3|z1,z2 )  =  f(z1,z2,z3)/f(z1,z2) 
				# 				= sqrt(2*pi) sqrt(det(Sigma12))/sqrt(det(Sigma))  exp( - 1/2  ((z1,z2,z3) - m)^t Sigma^(-1) ((z1,z2,23) - m)  + 1/2 ((z1,z2) - m12)^t Sigma12^(-1) ((z1,z2) - m12) )
				# Maximiser ça en z3 revient à minimiser ((z1,z2,z3) - m)^t Sigma^(-1) ((z1,z2,z3) - m) en z3  -> Polynôme de degrès 2 
				# P(z3) = sum_isum_j Sigma^(-1)_{i,j} (z_i - m_i) Sigma^(-1)_{i,j}(z_j - m_j)
				# P'(z3) =  2(z_3 - m_3) Sigma^(-1)_{3,3} + 2 sum_i (z_i - m_i) Sigma^(-1)_{i,3}
				# P'(z3) = 0 ssi z_3 = m_3 - 1/( Sigma^(-1)_{3,3} ) sum_i (z_i - m_i) Sigma^(-1)_{i,3}
				# z_3 = m_3  - 1/( Sigma^(-1)_{3,3} )  Sigma^(-1)_{pred,3} (zpred- mpred)
				m 			<- params$m
				
				predNum 	<- as.data.frame(matrix(NA,dim(Zpred)[1],length(colNumMissing)))
				idNumToPred <- which(idColToPred%in%colNum)
				for(x in 1:r)
				{
					Sigmax 		<- params$Sigma[[x]]
					mx  		<- m[x,]
					Sigmainv 	<- solve(Sigmax)
					posx 		<- which(XpredHat==x)
					if(length(posx)>0)
					{
						Zpredx 		<- as.matrix( Zpred[posx, colNumPred ,drop=FALSE])
							
							
						for (j in 1:length(colNumMissing))
						{
							col_j 				<- colNumMissing[j]
							Sigmainv_jj 		<- Sigmainv[col_j,col_j] 
							m_jj				<- mx[col_j]
							predNum[posx,j] 	<- m_jj - Sigmainv_jj^(-1) * t(t(Zpredx) -mx[-colNumMissing] ) %*%Sigmainv[-col_j,col_j,drop=FALSE]
						}
					
					}
					
				}
			
			}else{
				m 			<- params$m
				predNum 	<- as.data.frame(matrix(NA,dim(Zpred)[1],length(colNumMissing)))
				{
					 
					mx  		<- m[x,]
					 
					posx 		<- which(XpredHat==x)
				 
					for (j in 1:length(colNumMissing))
					{
						col_j 		<- colNumMissing[j]
						
						m_jj		<- mx[col_j]
						predNum[posx,j] 		<- m_jj  
					}
				
				}
			
			}
			ZpredVals[,which(idColToPred%in%colNum)] <- predNum
		
		}
		
		
		
		return(ZpredVals)
		
	}	
}


# calcul temps écoulé et à venir 
#' @keywords internal
remainingTime <- function(T,indice,nbIndices)
{
	deltaT = difftime(Sys.time(),T,units="mins")
	remaining = deltaT*(nbIndices-indice)/indice
	print(paste0("elapsed = ",round(deltaT)," mins "," remaining = ",round(remaining)," mins "))
	return(list(elapsed= deltaT, remaining = remaining))
}

#' @keywords internal
vectToMat 		<- function(phi,r)
{
	n = length(phi)/r
	return(matrix(phi,n,r,byrow=FALSE))
}

#' @keywords internal
matToVect 		<- function(phiM)
{
	return(as.vector(phiM))
}

#' @keywords internal
vectMatCoords 	<- function(n,r)
{
	return(cbind(matToVect(matrix(1:n,n,r,byrow=FALSE)),matToVect(matrix(1:r,n,r,byrow=TRUE))))						
}

#' @keywords internal
em_mixture_mixed <- function(df, r = 2, max_iter = 100, tol = 1e-6){
  #set.seed(seed)
  n <- nrow(df)
  
  # Séparation des variables numériques et catégorielles
  is_fact <- sapply(df, is.factor)
  Xnum <- df[, !is_fact, drop = FALSE]
  Xcat <- df[, is_fact, drop = FALSE]
  
  # Transformation des catégories en entiers (1:k)
  Xcat_enc <- lapply(Xcat, function(col) as.integer(col))
  levels_list <- lapply(Xcat, nlevels)
  Xcat_mat <- as.data.frame(Xcat_enc)
  
  # Initialisation aléatoire
  phi <- matrix(runif(n * r), nrow = n)
  phi <- phi / rowSums(phi)
  pi_k <- colMeans(phi)
  
  # Moyennes et variances pour les numériques
  means <- matrix(0, r, ncol(Xnum))
  vars <- matrix(1, r, ncol(Xnum))  # éviter variance nulle
  for (k in 1:r){
	means[k, ] <- colSums(phi[, k] * Xnum) / sum(phi[, k])
	vars[k, ] <- apply(Xnum, 2, function(x) sum(phi[, k] * (x - means[k, ])^2)) / sum(phi[, k])
  }
  
  # Probas multinomiales pour les catégorielles
  probs_cat <- lapply(1:ncol(Xcat_mat), function(j){
	matrix(1 / levels_list[[j]], nrow = levels_list[[j]], ncol = r)
  })
  
  loglik_prev <- -Inf
  
  for (iter in 1:max_iter){
	# E-step
	log_phi <- matrix(0, n, r)
	
	for (k in 1:r){
	  # Numérique : densité normale
	  if (ncol(Xnum) > 0){
		logdens_num <- rowSums(dnorm(as.matrix(Xnum), 
									 mean = matrix(means[k, ], n, ncol(Xnum), byrow = TRUE),
									 sd = sqrt(matrix(vars[k, ], n, ncol(Xnum), byrow = TRUE)),
									 log = TRUE))
	  } else {
		logdens_num <- 0
	  }
	  
	  # Catégoriel : produit des probabilités
	  logdens_cat <- rep(0, n)
	  if (ncol(Xcat_mat) > 0){
		for (j in 1:ncol(Xcat_mat)){
		  cat_vals <- Xcat_mat[[j]]
		  probs <- probs_cat[[j]][cat_vals, k]
		  logdens_cat <- logdens_cat + log(probs)
		}
	  }
	  
	  log_phi[, k] <- log(pi_k[k]) + logdens_num + logdens_cat
	}
	
	# Normalisation (log-sum-exp)
	max_log_phi <- apply(log_phi, 1, max)
	log_phi <- log_phi - max_log_phi
	phi <- exp(log_phi)
	phi <- phi / rowSums(phi)
	
	# M-step
	pi_k <- colMeans(phi)
	
	for (k in 1:r){
	  if (ncol(Xnum) > 0){
		means[k, ] <- colSums(phi[, k] * Xnum) / sum(phi[, k])
		vars[k, ] <- apply(Xnum, 2, function(x) sum(phi[, k] * (x - means[k, ])^2)) / sum(phi[, k])
	  }
	}
	
	for (j in 1:ncol(Xcat_mat)){
	  L <- levels_list[[j]]
	  probs_cat[[j]] <- matrix(0, L, r)
	  for (l in 1:L){
		for (k in 1:r){
		  probs_cat[[j]][l, k] <- sum(phi[Xcat_mat[[j]] == l, k])
		}
	  }
	  # Normaliser
	  probs_cat[[j]] <- probs_cat[[j]] / matrix(colSums(probs_cat[[j]]), L, r, byrow = TRUE)
	}
	
	# Log-vraisemblance
	loglik <- sum(log(rowSums(exp(log_phi + max_log_phi))))
	if (abs(loglik - loglik_prev) < tol){
	  break
	}
	loglik_prev <- loglik
  }
  
  phiM2 <- phi
  return(phiM2)
}
						   
#' @keywords internal
phiInitR			<- function(n,r)
{
	phi  		<- runif(n*r)
	phiM 		<- vectToMat(phi,r)
	sumByRow 	<- apply(phiM,1,sum)
	phiM2 		<- (1/sumByRow)*phiM  # Divise chaque colonne par le même vecteur 
	 
	phi 		<- matToVect(phiM2)
	return(phi)
}
						   
#' @keywords internal
phiInit			<- function(n,r,Z = NULL)
{
	if(is.null(Z))
	{
		phi 		<- phiInitR(n,r)
	}else{
		phiM2 <- em_mixture_mixed(df = Z, r = r)
		phi 		<- matToVect(phiM2)
	}
 
	

	return(phi)
}
						   
#' @keywords internal
phiToNu 		<- function(phi,lambda,C,r,isPhiAlreadyMat = FALSE)
{
	if(!isPhiAlreadyMat)
		phi <- vectToMat(phi=phi,r=r)
		
	return(apply(phi,2,sum)/(dim(phi)[1]))
}	
						   
#' @keywords internal
CECdens	 		<- function(Z,param,lambda=1,applyLog=FALSE)
{
	familyType <- param$familyType
	eval(parse(text = paste0("S <- dens_",familyType,"(Z=Z,param=param,lambda=lambda,applyLog=applyLog)")))
	 
	return(S)
	
}

						   
#' @keywords internal
optPhi			<- function(Z,param,lambda=1)
{
	n 				<- length(Z)
	if(is.matrix(Z)|is.data.frame(Z) )
		n <- dim(Z)[1]
		
	nu 				<- param$nu 
	r 				<- length(nu)
	logG 			<- CECdens(Z=Z,param=param,lambda=lambda,applyLog=TRUE)
	logG$density 	<- rep(log(nu),each=n)  +logG$density
	densityM 		<- vectToMat(logG$density,r)
	loc				<- apply(densityM,1,which.max)
	phi				<- rep(0,n*r)
	phi[(1:n)+ n*(loc-1) ] <- 1
	
	return(phi)			
}


# Fonctions propres à classe fonction 
{
	# gaussUniv
	{
		#' @keywords internal
		phiToMean 		<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
				
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			 
			n <- dim(phi)[1] 
			
			return(apply(Z*phi,2,sum)/(n*nuPhi))
		}
		
		#' @keywords internal
		phiToVar 		<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL,mPhi=NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			if(is.null(mPhi))
				mPhi 	<- phiToMean(phi=phi,Z=Z,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE,nuPhi =nuPhi)
			
			msquarePhi	<- phiToMean(phi=phi,Z=Z^2,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE,nuPhi =nuPhi)
			 
			return(msquarePhi - mPhi^2)
		}
		
		#' @keywords internal		
		dens_gaussUniv 				<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			states <- param$states
			m <- param$m
			s <- param$s
			
			r <- length(m)
			n <- length(Z) 
			
			S <- data.frame(state = rep(states,each = n) ,obs =  rep(1:n, r),density = rep(1/C,n))        
			for(i in 1:r)
			{
				state_i <- states[i]
				m_i		<- m[i]
				s_i 	<- s[i]
				if(applyLog)
				{
					S$density[1:n + (i-1)*n]<- dnorm(Z,mean=m_i,sd=s_i,log=TRUE) -log(lambda)
				}else{
					S$density[1:n + (i-1)*n]<- (dnorm(Z,mean=m_i,sd=s_i))^(1/lambda)
				}
				
			}
		 
			return(S)
			
		}
		 
		
		#' @keywords internal
		optParam_gaussUniv			<- function(Z,phi,lambda=1,C=1)
		{
			nr 		<- length(phi)
			n 		<- length(Z)
			r 		<- nr/n
			
			nu 		<- phiToNu(phi,lambda,C,r,isPhiAlreadyMat = FALSE)
			
			statesToKeep 	<- which(nu>0) 
			if(length(statesToKeep)<r)
			{
				
				indToKeep <- which(rep(1:r,each=n) %in%statesToKeep)
				phi <- phi[indToKeep]
				r 	<- length(statesToKeep)
				nu 	<- nu[statesToKeep]
			}
				
			states 	<- 1:r
			m 		<- phiToMean(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi = nu)
			s		<- sqrt(phiToVar(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi = nu,mPhi=m))
			
			# NON C est la borne sur le max densité, indépendemment de lambda : 1/(sqrt(2*pi)*s) <=C ssi s >= 1/(sqrt(2*pi)*C)
			
			functionBoundReached <- which(s<1/(sqrt(2*pi)*C))
			s[functionBoundReached] =  1/(sqrt(2*pi)*C) 
			
			params	<- list(states=states,nu=nu,m=m,s=s,lambda=lambda,C=C,familyType ="gaussUniv",functionBoundReached=functionBoundReached)
			
			return(params)
		}

		#' @keywords internal
		evalCompositeEntropy_gaussUniv	<-function(phi ,Z,lambda,C,includeClassEntropy = TRUE)
		{
			n 		<- length(Z)
			nr 		<- length(phi)
			r 		<- nr/n
			
			s2Min 	<- 1/(sqrt(2*pi)*C^lambda) 
			
			
			nuPhi 	<- phiToNu(phi,lambda,C,r,isPhiAlreadyMat = FALSE)
			varPhi 	<- phiToVar(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =nuPhi,mPhi=NULL)
			
			toKeep 	<- which(nuPhi>0)
			nuPhi	<- nuPhi[toKeep]
			varPhi	<- varPhi[toKeep]
			varSupS2min	<-  varPhi
			varSupS2min[varSupS2min<s2Min] = s2Min
			
			H <- 0 
			if(includeClassEntropy)
				H <- -sum(nuPhi*log(nuPhi)) 
			
			H <- H + (1/(2*lambda))* 	(	log(2*pi) +
															sum(nuPhi*log( varSupS2min  )) +
															sum(nuPhi*varPhi/varSupS2min) 
														) 
			return( H )
		}

	}

	# gaussVector
	{
		#' @keywords internal
		weighted_cov <- function(res, w) {
		  # Normaliser les poids
		  w <- w / sum(w)
		  
		  # Moyenne pondérée par colonne
		  mean_w <- colSums(w * res)
		  
		  # Centrer les données
		  res_centered <- sweep(res, 2, mean_w, "-")
		  
		  # Calcul de la matrice de covariance pondérée
		  cov_w <- t(res_centered) %*% (res_centered * w)
		  
		  return(cov_w)
		}

		#' @keywords internal
		adjust_covariance 				<- function(res, sigmaZmin,w=NULL,sizeMinSigmaComp=3)
		{
		 
		  # res matrice de résidius 
		  nRes <- dim(res)[1]
		  # estimation empirique
		  
		  if(is.null(w))
		  {
			w =rep(1/nRes,nRes)
		  
		  }
		  
		  if(sum(w)>sizeMinSigmaComp)
		  {
			  
			  Sigma_emp <- weighted_cov(res,w)# on enlève la correction
			  Sigma_emp[is.na(Sigma_emp)] <- 0 
			  # variances empiriques
			  variances <- diag(Sigma_emp)
			  
			  
			  # correction des variances
			  variances_corrected <- pmax(variances, sigmaZmin^2)
			  diag(Sigma_emp) <- variances_corrected
			  
			  
			  detMin <- prod(sigmaZmin^2)
			  
			  detSigmaEmp <- det(Sigma_emp)
			  boundReached <- FALSE
			  if(detSigmaEmp<detMin)
			  {
				boundReached <- TRUE
				cor_emp <- cov2cor(Sigma_emp)
				D_sqrt <- sqrt(diag(diag(Sigma_emp)))
			
				eig_R <- eigen(cor_emp, symmetric = TRUE)
				lambda_vals0 <- eig_R$values
				lambda_vals  <- lambda_vals0
				vectors 	<- eig_R$vectors
				
				
				lambda_vals[lambda_vals<0] = abs(lambda_vals[lambda_vals<0])
				
				det_R 	<- prod(lambda_vals)
				notNull <- which(lambda_vals>0)
				if(det_R == 0 )
				{
					if(length(notNull)>0)
					{
						lambda_vals[lambda_vals==0] = min(lambda_vals[notNull])/2
					}else{
						lambda_vals[lambda_vals==0] = 1 
					}					
				}
				
				det_R 	<- prod(lambda_vals)
				det_Rmin <- detMin/ prod(diag(Sigma_emp))
				
				if(det_R<det_Rmin)
				{
					lambda_vals <- lambda_vals*(det_Rmin/det_R)^(1/dim(Sigma_emp)[1])
					
					
					cor_emp_corrected <- (vectors)%*%diag(lambda_vals,nrow=dim(Sigma_emp)[1],ncol=dim(Sigma_emp)[1])%*%t(vectors)
				
					Sigma_corrected <- D_sqrt %*% cor_emp_corrected %*% D_sqrt
				
				}else{
				
					Sigma_corrected <- Sigma_emp
				}
				
				
 
				 
			  }else{
				Sigma_corrected = Sigma_emp
			  }
			  
			#Contrôle SigmaCorrected
			{
				
				# controle NA
				contrNA <- length(which(is.na(Sigma_corrected)))>0
				if(contrNA)
				{
					boundReached		<- TRUE
					variances 			<- diag(Sigma_emp)
					varNa 				<- is.na(variances)
					variances[varNa] 	<- sigmaZmin[varNa]^2
			  
					# correction des variances
					variances_corrected <- pmax(variances, sigmaZmin^2)
					Sigma_corrected <- diag(variances,nrow=length(variances),ncol=length(variances))
				}else{
					Sigma_corrected = (Sigma_corrected +t(Sigma_corrected))/2
					
					if(det(Sigma_corrected)<=0)
					{
						boundReached		<- TRUE
						variances 			<- diag(Sigma_emp)
						varNa 				<- which(is.na(variances))
						variances[varNa] 	<- sigmaZmin[varNa]^2
				  
						# correction des variances
						variances_corrected <- pmax(variances, sigmaZmin^2)
						Sigma_corrected <- diag(variances,nrow=length(variances),ncol=length(variances))
					}
				}
			}
			   
		  }else{
			boundReached		<- TRUE
			Sigma_corrected 	<-  diag(sigmaZmin,ncol=length(sigmaZmin),nrow=length(sigmaZmin))
		  }
		  
		 
		  
		  return(list(Sigma_corrected=Sigma_corrected,boundReached=boundReached))
		}

		#' @keywords internal
		phiToMeanVec 					<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
				
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			 
			n <- dim(phi)[1] 
			
			return( t(phi)%*%Z /(n*nuPhi))
		}

		#' @keywords internal
		phiToCovMat 					<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL,mPhi=NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			if(is.null(mPhi))
				mPhi 	<- phiToMeanVec(phi=phi,Z=Z,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE,nuPhi =nuPhi)
			 
			r <- length( nuPhi)
			
			Sigma <- list()
			for(i in 1:r)
			{
				phi_i 			<- phi[,i]/sum(phi[,i])
				Z_centered_i	<- sweep(Z, 2, mPhi[i,], "-")
				Sigma[[i]] 		<- t(Z_centered_i) %*% (Z_centered_i * phi_i)
			}
			return(Sigma)
		}

		#' @keywords internal
		dens_gaussVector 				<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			if(!is.matrix(Z))
				Z <- as.matrix(Z)
				
			l 		<- dim(Z)[2] 
			n 		<- dim(Z)[1]  
			
		 
			states 	<- param$states 
			
			r 		<- length(states)
			m 		<- param$m    		# matrice de taille r x l 
			Sigma	<- param$Sigma	    # list de taille r contenant pour tout i = 1,...,r, une matrice de corrélation de taille lxl
			C 		<- param$C 
			S 		<- data.frame(state = rep(states,each = n) ,obs =  rep(1:n, r),density = rep(1/C,n))        
			for(i in 1:r)
			{
				state_i 	<- states[i]
				m_i			<- m[i,]
				Sigma_i 	<- Sigma[[i]]
				
				if(applyLog)
				{
					S$density[1:n + (i-1)*n]<-  dmvnorm(Z,mean=m_i, sigma = Sigma_i ,log=TRUE) -log(lambda)
				}else{
					S$density[1:n + (i-1)*n]<- ( dmvnorm(Z,mean=m_i, sigma = Sigma_i ))^(1/lambda)
				}
				
			}
		 
			return(S)
			
		}

	 
		#' @keywords internal
		optParam_gaussVector 			<- function(Z,phi,lambda=1,C=1)
		{
			if(!is.matrix(Z))
				Z <- as.matrix(Z)
			
			n 		<- dim(Z)[1]
			l 		<- dim(Z)[2]
			
			nr 		<- length(phi)
		
			r 		<- nr/n
				
			
			nu 		<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = FALSE)
			
			statesToKeep 	<- which(nu>0) 
			if(length(statesToKeep)<r)
			{
				
				indToKeep <- which(rep(1:r,each=n) %in%statesToKeep)
				phi <- phi[indToKeep]
				r 	<- length(statesToKeep)
				nu 	<- nu[statesToKeep]
			}
			
		
			states 	<- 1:r
			
			
			
			m 		<- phiToMeanVec(phi=phi,Z=Z,lambda=lambda,C=C,r=r,isPhiAlreadyMat = FALSE,nuPhi = nu)
			
			
		 
			sigmaZmin <- rep(1/(sqrt(2*pi)*C^(1/l)) ,l)
			
			functionBoundReached <-  c()
			
			phiM <- vectToMat(phi=phi,r=r)
			Sigma <-list()
			for(i in 1:r)
			{	
				w 			<-  phiM[,i]
				sigList     <- adjust_covariance(res=Z, sigmaZmin=sigmaZmin,w=w,sizeMinSigmaComp=3)
				Sigma[[i]] 	<- sigList$Sigma_corrected
				if(sigList$boundReached)
				{
					functionBoundReached <- c(functionBoundReached,i)
				}
					
			}
			 
			
			params	<- list(states=states,nu=nu,m=m,Sigma=Sigma,lambda=lambda,C=C,familyType ="gaussVector",phi=phi,functionBoundReached=functionBoundReached)
			
			return(params)
		}

		#' @keywords internal
		evalCompositeEntropy_gaussVector 	<- function(phi ,Z,lambda,C,includeClassEntropy = TRUE)
		{	
			if(!is.matrix(Z))
				Z <- as.matrix(Z)
			
			n 		<- dim(Z)[1]
			
			
			params <- optParam_gaussVector(Z= Z,phi=phi,lambda=lambda,C=C)
			  
			nuPhi 	<- params$nu
			
			phiNew 	<- params$phi
			nr 		<- length(phiNew)
			r 		<- nr/n
			
			toKeep 	<- which(nuPhi>0) 
			 
			phiM <- vectToMat(phi=phiNew,r=r)
		
		
		
			logg 		<- CECdens(Z,param=params ,lambda=lambda,applyLog=TRUE)
			 
			H 		<- 0
			if(includeClassEntropy)
				H 		<- -sum(nuPhi[toKeep]*log(nuPhi[toKeep]))
			
			for(i in toKeep)
			{
				H  <- H - (1/lambda)*nuPhi[i]* sum( phiM[,i] *logg$density[logg$state==i])/sum(phiM[,i])
			}
			
			return(H)
		}

	}
	
	# discreteVector
	{
		#' @keywords internal
		dens_discreteVector 				<- function(Z,param,lambda=1,applyLog=FALSE)
		{
			
			factors <- param$factors
			
			if(!is.data.frame(Z))
			{
				Z <- as.data.frame(Z)
			}
			
			l 		<- dim(Z)[2] 
			for(j in 1:l)
			{
				Z[,j]= factor(Z[,j],levels= factors)
			}
			
			l 		<- dim(Z)[2] 
			n 		<- dim(Z)[1]  
			
		 
			states 	<- param$states 
			
			r 		<- length(states)
			
			discreteProbList <- param$discreteProbList
			
			
			 
			S 		<- data.frame(state = rep(states,each = n) ,obs =  rep(1:n, r),density = rep(NA,n))        
			for(i in 1:r)
			{
				state_i 	 	<- states[i]
				discreteProb_i 	<- discreteProbList[[i]] 
				if(applyLog)
				{
					S$density[1:n + (i-1)*n] <- 0 
					for(j in 1:l)
					{
						S$density[1:n + (i-1)*n] <- S$density[1:n + (i-1)*n] + log(discreteProb_i[Z[,j],j] )
						
						
					}
				}else{
					S$density[1:n + (i-1)*n] <- 1
					for(j in 1:l)
					{
						S$density[1:n + (i-1)*n] <- S$density[1:n + (i-1)*n] * discreteProb_i[Z[,j],j] 
						
						
					}
				}
				 
				
				 
				
			}
		 
			return(S)
			
		}

		#' @keywords internal
		optParam_discreteVector 			<- function(Z,phi,lambda=1,C=1)
		{
			 
			if(!is.data.frame(Z))
			{
				Z 		<- as.data.frame(Z)
			}
			
			{
				l 		<- dim(Z)[2] 
				factors = c()
				for(j in 1:l)
				{
					Z[,j] 	<- as.factor(Z[,j])
					factors <- c(factors,levels(Z[,j]))
				}
				
				factors <- sort(unique(factors))
				for(j in 1:l)
				{
					Z[,j] 	<-  factor(Z[,j],levels=factors)
				}						
			}
			# Ici il faut que toutes les coordonnées de Z soient de type factor et de mêmes niveaux
			
			n 		<- dim(Z)[1]
			l 		<- dim(Z)[2]
			
			nr 		<- length(phi)
		
			r 		<- nr/n
				
			
			nu 		<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = FALSE)
			
			statesToKeep 	<- which(nu>0) 
			if(length(statesToKeep)<r)
			{
				
				indToKeep <- which(rep(1:r,each=n) %in%statesToKeep)
				phi <- phi[indToKeep]
				r 	<- length(statesToKeep)
				nu 	<- nu[statesToKeep]
			}
			
			C_by_coord = rep(1/0.0005,l)
			if(TRUE)
			{
				uniqueFactors  = list()
				for(j in 1:l)
				{
					uniqueFactors[[j]] <- unique(Z[,j])
					C_by_coord[j]  = max(c(10*length(uniqueFactors[[j]]),C))
				}
				
			}
			
			
			states 	<- 1:r
			phiM <- vectToMat(phi,r)
			 
			
			factors = levels(Z[,1])
			discreteProbList <- list()
			
			functionBoundReached <-  c()
			for(i in 1:r)
			{
				discreteProbList[[i]] <- matrix(NA,length(factors),l)
				for(j in 1:l)
				{
					tbl_ij 					<- tapply(phiM[,i], Z[,j], sum)
					tbl_ij 					<- tbl_ij/sum(tbl_ij,na.rm=TRUE)
					
					toLift					<- which(!is.na(tbl_ij) & tbl_ij<1/C_by_coord[j])
					toShrink    			<- which(!is.na(tbl_ij) & tbl_ij>=1/C_by_coord[j])
					
					if(length(toLift)>0)
						functionBoundReached <- c(functionBoundReached,i)
						
					tbl_ij[toLift]  		<- 1/C_by_coord[j] 
					tbl_ij[toShrink] 		<- tbl_ij[toShrink]/(sum(tbl_ij[toShrink]))*(1-sum(tbl_ij[toLift]) )
					
					tbl_ij[is.na(tbl_ij)] 	<- 0 
					discreteProbList[[i]][,j] <- tbl_ij
				}
			}
			
			functionBoundReached <- unique(functionBoundReached)
			 
			
			params	<- list(states=states,nu=nu,factors=factors,discreteProbList=discreteProbList,lambda=lambda,C=C,familyType ="discreteVector",phi=phi)
			
			return(params)
		}

		#' @keywords internal
		evalCompositeEntropy_discreteVector 	<- function(phi ,Z,lambda,C,includeClassEntropy = TRUE)
		{	
			if(!is.data.frame(Z))
			{
				Z 		<- as.data.frame(Z)
				
				l 		<- dim(Z)[2] 
				factors = c()
				for(j in 1:l)
				{
					Z[,j] 	<- as.factor(Z[,j])
					factors <- c(factors,levels(Z[,j]))
					
				}
				
				factors <- sort(unique(factors))
				for(j in 1:l)
				{
					Z[,j] 	<-  factor(Z[,j],levels=factors)
				}						
			}
			
			n 		<- dim(Z)[1]
			
			
			params <- optParam_discreteVector(Z= Z,phi=phi,lambda=lambda,C=C)
			  
			nuPhi 	<- params$nu
			
			phiNew 	<- params$phi
			nr 		<- length(phiNew)
			r 		<- nr/n
			
			toKeep 	<- which(nuPhi>0) 
			 
			phiM <- vectToMat(phi=phiNew,r=r)
		
		
		
			logg 		<- CECdens(Z,param=params ,lambda=lambda,applyLog=TRUE)
			
			H 		<- 0
			if(includeClassEntropy)
				H 		<- -sum(nuPhi[toKeep]*log(nuPhi[toKeep]))
			
			for(i in toKeep)
			{
				H  <- H - (1/lambda)*nuPhi[i]* sum( phiM[,i] *logg$density[logg$state==i])/sum(phiM[,i])
			}
			
			return(H)
		}

	
	}

	# gaussAndDiscreteVector
	{
		#' @keywords internal
		dens_gaussAndDiscreteVector 		<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			 
			if(!is.data.frame(Z))
				Z <- as.data.frame(Z)
			 
			
			l 		<- dim(Z)[2] 
			n 		<- dim(Z)[1]  
			
			states 	<- param$states 
			
			r 		<- length(states)
			
			S 		<- data.frame(state = rep(states,each = n) ,obs =  rep(1:n, r),density = rep(NA,n))   
			if(applyLog)
			{
				S$density = 0
				
			}else{
				S$density = 1
			}
			
			
			Zclass <- sapply(Z, class)
			
			if(length(which(Zclass=="numeric"))>0)
			{
				SNum 		<- dens_gaussVector(Z=Z[,which(Zclass=="numeric"),drop=FALSE],param=param,lambda=lambda,applyLog=applyLog)
				S$density 	<- SNum$density  
			}
				
			
			if(length(which(Zclass=="factor"))>0)
			{
				SFac <- dens_discreteVector(Z=Z[,which(Zclass=="factor"),drop=FALSE],param=param,lambda=lambda,applyLog=applyLog)
				
				if(applyLog)
				{
					S$density = S$density  + SFac$density
					
				}else{
					S$density = S$density  * SFac$density
				}
			}
			 
			 
			return(S)
			
		}

		#' @keywords internal
		optParam_gaussAndDiscreteVector		<- function(Z,phi,lambda=1,C=1)
		{
			  
			if(!is.data.frame(Z))
				Z <- as.data.frame(Z)
			 
		
			Zclass <- sapply(Z, class)
			params	<- list()
			
			colFactor 	<- which(Zclass=="factor")
			colNum  	<- which(Zclass=="numeric")
			functionBoundReached <- c()
			
			if(length(colFactor)>0)
			{
				paramsFact  <- optParam_discreteVector(Z= Z[,which(Zclass=="factor"),drop=FALSE],phi=phi,lambda=lambda,C=C)
				
				params$states 	<- paramsFact$states
				params$nu		<- paramsFact$nu
				params$factors	<- paramsFact$factors
				params$discreteProbList	<- paramsFact$discreteProbList
				params$lambda	<- paramsFact$lambda
				params$C	 	<- paramsFact$C
				params$phi	 	<- paramsFact$phi
				
				functionBoundReached <- paramsFact$functionBoundReached
				
			}
			
			if(length(colNum)>0)
			{
				paramsNum  		<- optParam_gaussVector(Z= Z[,which(Zclass=="numeric"),drop=FALSE],phi=phi,lambda=lambda,C=C)
				
				params$states 	<- paramsNum$states
				params$nu		<- paramsNum$nu
				params$m		<- paramsNum$m
				params$Sigma	<- paramsNum$Sigma
				params$lambda	<- paramsNum$lambda
				params$C	 	<- paramsNum$C
				params$phi	 	<- paramsNum$phi
				
				functionBoundReached <-c( functionBoundReached,paramsNum$functionBoundReached)
				
			}
			functionBoundReached <- unique(functionBoundReached)
			
			params$functionBoundReached <- functionBoundReached
			params$familyType	<- "gaussAndDiscreteVector"
			params$colFactor	<- colFactor
			params$colNum		<- colNum
			
			return(params)
		}

		#' @keywords internal
		evalCompositeEntropy_gaussAndDiscreteVector 	<- function(phi ,Z,lambda,C,includeClassEntropy = TRUE)
		{	
			if(!is.data.frame(Z))
				Z <- as.data.frame(Z)
	
			n 		<- dim(Z)[1]
			
			
			params <- optParam_gaussAndDiscreteVector(Z= Z,phi=phi,lambda=lambda,C=C)
			  
			nuPhi 	<- params$nu
			 
			toKeep 	<- which(nuPhi>0) 
			 
			 
			H 		<- 0
			if(includeClassEntropy)
				H 		<- -sum(nuPhi[toKeep]*log(nuPhi[toKeep]))
			
			Zclass <- sapply(Z, class)
			
			if(length(which(Zclass=="factor"))>0)
			{
				H <- H + evalCompositeEntropy_discreteVector(phi=phi,Z=Z[,which(Zclass=="factor"),drop=FALSE],lambda=lambda,C=C,includeClassEntropy=FALSE)
			}
			
			if(length(which(Zclass=="numeric"))>0)
			{
				H <- H + evalCompositeEntropy_gaussVector(phi=phi,Z=Z[,which(Zclass=="numeric"),drop=FALSE],lambda=lambda,C=C,includeClassEntropy=FALSE)
			}
			
			return(H)
		}

		
	}
}
						   
#' @keywords internal
optParam 		<- function(Z,phi,lambda=1,C=1,familyType ="gaussAndDiscreteVector")
{		
	eval(parse(text = paste0("params <- optParam_",familyType,"(Z=Z,phi=phi,lambda=lambda,C=C)")))
	return(params)	
}

#' @keywords internal
evalCompositeEntropy <- function(phi ,Z,lambda,C,familyType="gaussAndDiscreteVector")
{
	eval(parse(text = paste0("H <- evalCompositeEntropy_",familyType,"(phi =phi,Z=Z,lambda=lambda,C=C)")))
	return(H)
}

#' @keywords internal
CECclassifOneShot 		<- function(Z,lambda=1,C=1,r0=NULL,Nloop=1000,phi0 = NULL,familyType="gaussAndDiscreteVector",displayPlotEntropy = FALSE)
{
	n <- length(Z)
	if(is.matrix(Z)|is.data.frame(Z) )
		n <- dim(Z)[1] 
	
	if(is.null(phi0))
	{
		if(is.null(r0))
		{
			r <- floor(log(n) + 1)
		}else{
			r <- r0			
		}
		 
		 
		phi		<- phiInit(n,r)
	}else{
		nr		<- length(phi0)
		r0		<- nr/n 
		r  		<- r0 
		phi 	<- phi0 
	}
	

	phiFingerPrint <- sum(phi*cos(2*pi*(0:(n*r-1))/n))
	
	H 		<- c()
	
	for(i in 1:Nloop)
	{
		params 	<- optParam(Z=Z,phi=phi,lambda=lambda,C=C,familyType=familyType)
		if(displayPlotEntropy)
		{
			 H[i]<- evalCompositeEntropy(phi=phi ,Z=Z,lambda=lambda,C=C,familyType=familyType)
			 plot(H,main="entropy",xlab = "nLoop",ylab="H")
		}
		
		phi 	<- optPhi(Z=Z,param=params,lambda=lambda)	
		nr 		<- length(phi) 
		r 		<- nr/n
		 
		if(r==1)
			break()
			
		phiFingerPrint_i <- sum(phi*cos(2*pi*(0:(n*r-1))/n))
		if(phiFingerPrint_i%in%phiFingerPrint)
			break()
		
		phiFingerPrint <-c(phiFingerPrint,phiFingerPrint_i)
		
	}
	
	params 	<- optParam(Z=Z,phi=phi,lambda=lambda,C=C,familyType=familyType)
	phi 	<- optPhi(Z=Z,param=params,lambda=lambda)	
	nr 		<- length(phi) 
	r 		<- nr/n
	
	Hphi 	<- evalCompositeEntropy(phi=phi,Z=Z,lambda=lambda,C=C,familyType=params$familyType)
	return(list(phi=phi,params=params,Hphi=Hphi))
}
 
