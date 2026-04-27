# Internal helpers for the CEClust package.

.CEC_fast_backend <- new.env(parent = emptyenv())
.CEC_fast_backend$available <- NULL
.CEC_fast_backend$force_disable <- FALSE

CECdetect_source_file <- function() {
	frames <- sys.frames()
	ofile <- NULL
	if (length(frames) > 0) {
		for (idx in rev(seq_along(frames))) {
			if (!is.null(frames[[idx]]$ofile)) {
				ofile <- frames[[idx]]$ofile
				break
			}
		}
	}
	if (is.null(ofile)) {
		NA_character_
	} else {
		normalizePath(ofile, winslash = "/", mustWork = FALSE)
	}
}

CECget_package_root <- function() {
	pkg_root_opt <- getOption("CEClust.package_root", default = "")
	if (is.character(pkg_root_opt) && length(pkg_root_opt) >= 1L && nzchar(pkg_root_opt[1])) {
		return(normalizePath(pkg_root_opt[1], winslash = "/", mustWork = FALSE))
	}

	pkg_root <- tryCatch(system.file(package = "CEClust"), error = function(e) "")
	if (nzchar(pkg_root)) {
		return(normalizePath(pkg_root, winslash = "/", mustWork = FALSE))
	}
	source_file <- CECdetect_source_file()
	if (!is.na(source_file)) {
		return(dirname(dirname(source_file)))
	}
	normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

CECcompiled_backend_fun_names <- function() {
	c(
		"cec_cpp_choose_clusters",
		"cec_cpp_phi_from_clusters",
		"cec_cpp_gaussian_stats",
		"cec_cpp_logdens_gaussian",
		"cec_cpp_discrete_counts",
		"cec_cpp_logdens_discrete",
		"cec_cpp_cluster_fingerprint"
	)
}

CECcompiled_backend_available <- function() {
	ns <- environment(CECcompiled_backend_available)
	all(vapply(
		CECcompiled_backend_fun_names(),
		exists,
		logical(1),
		envir = ns,
		mode = "function",
		inherits = FALSE
	))
}

#' Reset cached information about the optional fast backend.
#'
#' @rdname CECconfigure_runtime
#' @param rebuild Logical. If `TRUE`, backend detection is run immediately after
#'   the cache is cleared.
#' @return `CECreset_fast_backend()` invisibly returns `TRUE`.
#' @export
CECreset_fast_backend <- function(rebuild = FALSE) {
	.CEC_fast_backend$available <- NULL
	if (isTRUE(rebuild)) {
		CECensure_fast_backend()
	}
	invisible(TRUE)
}

#' Enable or disable the optional fast backend.
#'
#' @rdname CECconfigure_runtime
#' @param enabled Logical. If `TRUE`, CEClust is allowed to use the compiled
#'   backend when it is available.
#' @param rebuild Logical. If `TRUE`, refresh backend detection immediately
#'   after toggling the backend state.
#' @return `CECset_fast_backend()` invisibly returns the named list produced by
#'   [CECperformanceInfo()].
#' @export
CECset_fast_backend <- function(enabled = TRUE, rebuild = enabled) {
	.CEC_fast_backend$force_disable <- !isTRUE(enabled)
	.CEC_fast_backend$available <- NULL
	if (isTRUE(enabled) && isTRUE(rebuild)) {
		CECensure_fast_backend()
	}
	invisible(CECperformanceInfo())
}

CECfast_backend_available <- function() {
	isTRUE(CECensure_fast_backend())
}

CECensure_fast_backend <- function() {
	if (isTRUE(.CEC_fast_backend$force_disable)) {
		.CEC_fast_backend$available <- FALSE
		return(FALSE)
	}

	if (isTRUE(.CEC_fast_backend$available)) {
		return(TRUE)
	}

	.CEC_fast_backend$available <- isTRUE(CECcompiled_backend_available())
	.CEC_fast_backend$available
}

CECget_fast_fun <- function(name) {
	if (!CECensure_fast_backend()) {
		stop("Fast backend is not available.")
	}
	get(name, envir = environment(CECcompiled_backend_available), inherits = FALSE)
}

#' Inspect the runtime capabilities detected by CEClust.
#'
#' @rdname CECconfigure_runtime
#' @return `CECperformanceInfo()` returns a named list containing
#'   `fast_backend_available`, `fast_backend_forced_disabled`,
#'   `make_available`, `BLAS`, `LAPACK`, and `detected_cores`.
#' @export
CECperformanceInfo <- function() {
	ext_versions <- extSoftVersion()
	get_ext_version <- function(name) {
		if (name %in% names(ext_versions)) {
			value <- unname(ext_versions[[name]])
			if (length(value) == 0 || is.na(value) || !nzchar(value)) {
				return(NA_character_)
			}
			return(value)
		}
		NA_character_
	}

	list(
		fast_backend_available = CECfast_backend_available(),
		fast_backend_forced_disabled = isTRUE(.CEC_fast_backend$force_disable),
		make_available = nzchar(Sys.which("make")),
		BLAS = get_ext_version("BLAS"),
		LAPACK = get_ext_version("LAPACK"),
		detected_cores = parallel::detectCores(logical = FALSE)
	)
}

CECis_fast_family <- function(familyType) {
	familyType %in% c("gaussVector", "discreteVector", "gaussAndDiscreteVector")
}

CECprepare_discrete_block <- function(Z) {
	if (is.null(Z) || ncol(Z) == 0) {
		return(NULL)
	}

	if (!is.data.frame(Z)) {
		Z <- as.data.frame(Z)
	}

	l <- ncol(Z)
	factors <- c()
	for (j in seq_len(l)) {
		Z[[j]] <- as.factor(Z[[j]])
		factors <- c(factors, levels(Z[[j]]))
	}
	factors <- sort(unique(factors))
	for (j in seq_len(l)) {
		Z[[j]] <- factor(Z[[j]], levels = factors)
	}

	codes <- do.call(cbind, lapply(Z, as.integer))
	if (!is.matrix(codes)) {
		codes <- matrix(codes, ncol = l)
	}
	storage.mode(codes) <- "integer"

	n_unique_by_coord <- vapply(
		seq_len(l),
		function(j) length(unique(codes[codes[, j] > 0L, j])),
		integer(1)
	)
	level_present_by_coord <- lapply(
		seq_len(l),
		function(j) tabulate(codes[, j], nbins = length(factors)) > 0L
	)

	list(
		df = Z,
		codes = codes,
		factors = factors,
		n_levels = length(factors),
		n_unique_by_coord = n_unique_by_coord,
		level_present_by_coord = level_present_by_coord
	)
}

CECprepare_backend_data <- function(Z, familyType = "gaussAndDiscreteVector") {
	n <- CECget_n_obs(Z)

	if (familyType == "gaussVector") {
		X <- if (is.matrix(Z)) Z else as.matrix(Z)
		return(list(
			familyType = familyType,
			n = n,
			optimized = TRUE,
			raw = Z,
			X_num = X
		))
	}

	if (familyType == "discreteVector") {
		discrete <- CECprepare_discrete_block(Z)
		return(list(
			familyType = familyType,
			n = n,
			optimized = TRUE,
			raw = Z,
			discrete = discrete
		))
	}

	if (familyType == "gaussAndDiscreteVector") {
		if (!is.data.frame(Z)) {
			Z <- as.data.frame(Z)
		}

		Zclass <- vapply(Z, function(x) class(x)[1], character(1))
		colFactor <- which(Zclass == "factor")
		colNum <- which(Zclass == "numeric")

		discrete <- if (length(colFactor) > 0) CECprepare_discrete_block(Z[, colFactor, drop = FALSE]) else NULL
		X_num <- if (length(colNum) > 0) as.matrix(Z[, colNum, drop = FALSE]) else NULL

		return(list(
			familyType = familyType,
			n = n,
			optimized = TRUE,
			raw = Z,
			colFactor = colFactor,
			colNum = colNum,
			discrete = discrete,
			X_num = X_num
		))
	}

	list(
		familyType = familyType,
		n = n,
		optimized = FALSE,
		raw = Z
	)
}

CECclusters_to_phi <- function(clusters, r) {
	clusters <- as.integer(clusters)
	if (CECfast_backend_available()) {
		return(CECget_fast_fun("cec_cpp_phi_from_clusters")(clusters, as.integer(r)))
	}

	n <- length(clusters)
	phi <- numeric(n * r)
	idx <- which(clusters >= 1L & clusters <= r)
	phi[idx + n * (clusters[idx] - 1L)] <- 1
	phi
}

CECclusters_from_phi_if_hard <- function(phi, n = NULL, r = NULL, tol = sqrt(.Machine$double.eps)) {
	if (is.null(phi)) {
		return(NULL)
	}

	if (is.null(n) || is.null(r)) {
		if (is.null(r)) {
			stop("r must be provided when converting phi to hard clusters.")
		}
		n <- length(phi) / r
	}

	phiM <- vectToMat(phi, r)
	row_sums <- rowSums(phiM)
	if (any(abs(row_sums - 1) > tol)) {
		return(NULL)
	}

	row_max <- max.col(phiM, ties.method = "first")
	test_mat <- matrix(0, n, r)
	test_mat[cbind(seq_len(n), row_max)] <- 1
	if (max(abs(phiM - test_mat)) > tol) {
		return(NULL)
	}

	as.integer(row_max)
}

CECcompress_clusters <- function(clusters) {
	if (is.null(clusters)) {
		return(NULL)
	}
	as.integer(match(clusters, sort(unique(clusters))))
}

CECcluster_fingerprint <- function(clusters) {
	clusters <- as.integer(clusters)
	if (CECfast_backend_available()) {
		return(CECget_fast_fun("cec_cpp_cluster_fingerprint")(clusters))
	}
	n <- length(clusters)
	if (n == 0) {
		return(0)
	}
	sum(clusters * cos(2 * pi * (seq_len(n) - 1L) / n))
}

CECentropy_from_assignment <- function(nu, assigned_logdens, lambda) {
	toKeep <- which(nu > 0)
	H_class <- 0
	if (length(toKeep) > 0) {
		H_class <- -sum(nu[toKeep] * log(nu[toKeep]))
	}
	H_cond <- -(1 / lambda) * mean(assigned_logdens)
	list(
		H_total = H_class + H_cond,
		H_class = H_class,
		H_cond = H_cond
	)
}

#' Fit a composite entropy clustering model for one value of lambda.
#'
#' `CECclassif()` is the main fixed-lambda fitting routine of the package. It
#' launches several random initialisations, keeps the best solution according to
#' the composite entropy criterion, and augments the returned object with a hard
#' partition and the realised number of occupied clusters (`REO`).
#'
#' @param Z Input data. Use a numeric vector or one-column numeric object for
#'   one-dimensional Gaussian clustering, a numeric matrix/data frame for
#'   Gaussian multivariate clustering, a factor data frame for discrete data, or
#'   a mixed data frame for `"gaussAndDiscreteVector"`.
#' @param lambda Regularisation parameter of the composite entropy criterion.
#' @param C Positive upper bound for fitted component densities; independent
#'   of `lambda`.
#' @param r0 Optional upper bound for the initial number of clusters.
#' @param Nshots Number of repeated random initialisations.
#' @param Nloop Maximum number of optimisation iterations per shot.
#' @param familyType Model family. Supported values are `"gaussUniv"`,
#'   `"gaussVector"`, `"discreteVector"`, and `"gaussAndDiscreteVector"`.
#' @param sizeMaxOutlier Maximum outlier size used by the optional regrouping
#'   logic.
#' @param autoRegroupOutliers Logical. If `TRUE`, small outlier groups may be
#'   merged automatically during fitting.
#' @param displayRemainingTime Logical. If `TRUE`, progress information is
#'   printed across shots.
#' @param focus Optional column index used by the focused workflow for factor
#'   conditioning.
#' @param backend_data Optional preprocessed backend object returned internally
#'   by CEClust to reuse numeric/factor encodings across repeated calls.
#'
#' @return A list describing the best solution found over all shots. The main
#'   components are:
#' - `phi`: partition weights stored in vectorised form;
#' - `params`: fitted component parameters and tuning settings;
#' - `Hphi`: composite entropy criterion at the selected solution;
#' - `REO`: realised number of occupied clusters;
#' - `clusters`: hard cluster assignment for each observation.
#'
#' @examples
#' CECconfigure_runtime("base")
#' Z <- simulate_multidim_benchmark_data(n = 80, p_num = 3, p_fac = 2, seed = 1)
#'
#' fit <- CECclassif(
#'   Z = Z,
#'   lambda = 0.8,
#'   C = 10,
#'   r0 = 6,
#'   Nshots = 3,
#'   Nloop = 15,
#'   familyType = "gaussAndDiscreteVector",
#'   backend_data = list(optimized = FALSE, raw = Z)
#' )
#'
#' fit$REO
#' table(fit$clusters)
#' @export
CECclassif				<- function(Z,lambda=1,C=1,r0=NULL,Nshots = 100,Nloop=1000,familyType="gaussAndDiscreteVector" ,sizeMaxOutlier = 0,autoRegroupOutliers=FALSE,displayRemainingTime = FALSE,focus=NULL,backend_data=NULL)
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

	if (is.null(backend_data)) {
		backend_data <- CECprepare_backend_data(Z, familyType = familyType)
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

				
				sumColsPhi 	<- rowSums(phiM) 
				phiM		<- sweep(phiM, 1, sumColsPhi, "/")
				phi0        <- matToVect(phiM)
			}
			
		}
		
		resultList 	<- CECrun_one_shot_with_regroup(
			Z = Z,
			lambda = lambda,
			C = C,
			r0 = r0,
			Nloop = Nloop,
			phi0 = phi0,
			familyType = familyType,
			displayPlotEntropy = displayPlotEntropy,
			sizeMaxOutlier = sizeMaxOutlier,
			autoRegroupOutliers = autoRegroupOutliers,
			backend_data = backend_data
		)
		
		
		
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
	
	return(CECfinalize_classif_result(bestClassif))
	
}

CECget_n_obs <- function(Z) {
	if (is.matrix(Z) || is.data.frame(Z)) {
		return(dim(Z)[1])
	}
	length(Z)
}

CECfinalize_classif_result <- function(resultList) {
	if (is.null(resultList) || is.null(resultList$phi) || is.null(resultList$params$states)) {
		return(resultList)
	}

	REO <- length(resultList$params$states)
	resultList$REO <- REO
	phiM <- vectToMat(resultList$phi, REO)
	resultList$clusters <- max.col(phiM, ties.method = "first")
	resultList
}

CECrun_one_shot_with_regroup <- function(
	Z,
	lambda = 1,
	C = 1,
	r0 = NULL,
	Nloop = 1000,
	phi0 = NULL,
	familyType = "gaussAndDiscreteVector",
	displayPlotEntropy = FALSE,
	sizeMaxOutlier = 0,
	autoRegroupOutliers = FALSE,
	backend_data = NULL
) {
	n <- CECget_n_obs(Z)

	if (is.null(backend_data)) {
		backend_data <- CECprepare_backend_data(Z, familyType = familyType)
	}

	resultList <- CECclassifOneShot(
		Z = Z,
		lambda = lambda,
		C = C,
		r0 = r0,
		Nloop = Nloop,
		phi0 = phi0,
		familyType = familyType,
		displayPlotEntropy = displayPlotEntropy,
		backend_data = backend_data
	)

	if (sizeMaxOutlier < n || autoRegroupOutliers) {
		continueOutlier <- TRUE
		while (continueOutlier) {
			sizeGroups <- resultList$params$nu * n

			if (autoRegroupOutliers) {
				groupRank <- rank(sizeGroups, ties.method = "min")
				functionBoundReached <- resultList$params$functionBoundReached

				candidateToRegroup <- which(groupRank <= floor(length(groupRank) / 2))
				candidateToRegroup <- candidateToRegroup[which(candidateToRegroup %in% functionBoundReached)]

				outliersToRegroup <- c()
				if (length(candidateToRegroup) > 0 && length(sizeGroups) > 1) {
					outliersToRegroup <- candidateToRegroup[which.min(groupRank[candidateToRegroup])]
				}
			} else {
				outliersToRegroup <- which(sizeGroups < sizeMaxOutlier)
			}

			if (length(outliersToRegroup) > 0) {
				r <- length(sizeGroups)
				phi <- resultList$phi
				phiM <- vectToMat(phi, r)

				if (length(outliersToRegroup) < r) {
					phiM0 <- phiM[, -outliersToRegroup, drop = FALSE]
				} else {
					phiM0 <- phiM[, 1, drop = FALSE]
				}

				dataToReassign <- which(rowSums(phiM[, outliersToRegroup, drop = FALSE]) > 0)
				phiM0[dataToReassign, ] <- runif(length(dataToReassign) * dim(phiM0)[2])
				sumColsphiM0Reassign <- rowSums(phiM0[dataToReassign, , drop = FALSE])
				phiM0[dataToReassign, ] <- sweep(
					phiM0[dataToReassign, , drop = FALSE],
					1,
					sumColsphiM0Reassign,
					"/"
				)

				phi0_next <- matToVect(phiM0)
				resultList <- CECclassifOneShot(
					Z = Z,
					lambda = lambda,
					C = C,
					r0 = r0,
					Nloop = Nloop,
					phi0 = phi0_next,
					familyType = familyType,
					displayPlotEntropy = displayPlotEntropy,
					backend_data = backend_data
				)
			} else {
				continueOutlier <- FALSE
			}
		}
	}

	resultList
}

CECcompose_fit_origin_label <- function(init_origin = NA_character_, path_direction = NA_character_) {
	init_origin <- as.character(init_origin)[1]
	path_direction <- as.character(path_direction)[1]

	if (is.na(init_origin) || !nzchar(init_origin)) {
		return(NA_character_)
	}

	if (identical(init_origin, "fresh")) {
		return("fresh")
	}

	if (identical(init_origin, "warm") && !is.na(path_direction) && nzchar(path_direction)) {
		return(paste0("warm_", path_direction))
	}

	init_origin
}

CECget_fit_origin_label <- function(fit) {
	if (is.null(fit)) {
		return(NA_character_)
	}

	if (!is.null(fit$init_origin_label)) {
		return(as.character(fit$init_origin_label)[1])
	}

	CECcompose_fit_origin_label(
		init_origin = if (!is.null(fit$init_origin)) fit$init_origin else NA_character_,
		path_direction = if (!is.null(fit$path_direction)) fit$path_direction else NA_character_
	)
}

CECget_fit_path_direction <- function(fit) {
	if (is.null(fit) || is.null(fit$path_direction)) {
		return(NA_character_)
	}

	label <- CECget_fit_origin_label(fit)
	if (!is.na(label) && identical(label, "fresh")) {
		return(NA_character_)
	}

	as.character(fit$path_direction)[1]
}

CECtag_fit_origin <- function(
	fit,
	init_origin = NA_character_,
	path_direction = NA_character_,
	candidate_group = NA_character_,
	candidate_index = NA_integer_
) {
	if (is.null(fit)) {
		return(NULL)
	}

	fit$init_origin <- as.character(init_origin)[1]
	fit$path_direction <- as.character(path_direction)[1]
	fit$init_origin_label <- CECcompose_fit_origin_label(init_origin, path_direction)
	fit$candidate_group <- as.character(candidate_group)[1]
	fit$candidate_index <- as.integer(candidate_index)[1]
	fit
}

CECcount_fit_origins <- function(fits) {
	labels <- vapply(fits, CECget_fit_origin_label, character(1))

	c(
		fresh = sum(labels == "fresh", na.rm = TRUE),
		warm_forward = sum(labels == "warm_forward", na.rm = TRUE),
		warm_backward = sum(labels == "warm_backward", na.rm = TRUE)
	)
}

CECbest_fit_origin_info <- function(fits, Hvals) {
	valid <- which(!is.na(Hvals) & vapply(fits, function(x) !is.null(x), logical(1)))

	if (length(valid) == 0) {
		return(list(
			init_origin = NA_character_,
			path_direction = NA_character_,
			index = NA_integer_
		))
	}

	best_idx <- valid[which.min(Hvals[valid])]

	list(
		init_origin = CECget_fit_origin_label(fits[[best_idx]]),
		path_direction = CECget_fit_path_direction(fits[[best_idx]]),
		index = best_idx
	)
}

CECperturb_best_warm_phi <- function(bestClassif, n) {
	r <- length(bestClassif$params$states)
	nu <- bestClassif$params$nu
	phi <- bestClassif$phi
	phiM <- vectToMat(phi, r)

	toKeep <- which(nu * n >= 10)
	if (length(toKeep) == 0) {
		toKeep <- which.max(nu)
	}

	phiM <- phiM[, toKeep, drop = FALSE]
	phiM <- phiM + matrix(0.1 * runif(prod(dim(phiM))), dim(phiM)[1], dim(phiM)[2])

	sumColsPhi <- rowSums(phiM)
	sumColsPhi[sumColsPhi <= 0] <- 1
	phiM <- sweep(phiM, 1, sumColsPhi, "/")
	matToVect(phiM)
}

CECgenerate_warm_shot_phi0 <- function(phi_seed, PHI, bestClassif, shot_index, Nshots, n) {
	if (shot_index < 2) {
		return(phi_seed)
	}

	if (!is.null(bestClassif) &&
		shot_index %in% seq.int(max(1L, Nshots - floor(Nshots / 10)), Nshots) &&
		is.finite(bestClassif$Hphi) &&
		n > 100) {
		return(CECperturb_best_warm_phi(bestClassif, n))
	}

	if (is.null(PHI) || length(PHI) == 0) {
		PHI <- phi_seed
	}

	U <- runif(1)
	if (U < 0.5) {
		return(c(PHI * runif(length(PHI)), runif(n)))
	}

	if (!is.null(bestClassif) && !is.null(bestClassif$phi)) {
		return(bestClassif$phi)
	}

	phi_seed
}

CECfit_warm_multi <- function(
	Z,
	lambda,
	phi0,
	C = 1,
	r0 = NULL,
	Nshots_warm = 5,
	Nloop = 300,
	familyType = "gaussAndDiscreteVector",
	sizeMaxOutlier = 0,
	autoRegroupOutliers = FALSE,
	path_direction = NA_character_,
	silent = TRUE,
	backend_data = NULL
) {
	if (is.null(phi0)) {
		return(NULL)
	}

	if (Nshots_warm < 1) {
		stop("Nshots_warm must be >= 1 when phi0 is provided.")
	}

	n <- CECget_n_obs(Z)
	if (is.null(backend_data)) {
		backend_data <- CECprepare_backend_data(Z, familyType = familyType)
	}
	phi_seed <- phi0
	PHI <- phi_seed
	H <- Inf
	bestClassif <- NULL

	for (k in seq_len(Nshots_warm)) {
		phi0_k <- CECgenerate_warm_shot_phi0(
			phi_seed = phi_seed,
			PHI = PHI,
			bestClassif = bestClassif,
			shot_index = k,
			Nshots = Nshots_warm,
			n = n
		)

		resultList <- tryCatch(
			CECrun_one_shot_with_regroup(
				Z = Z,
				lambda = lambda,
				C = C,
				r0 = r0,
				Nloop = Nloop,
				phi0 = phi0_k,
				familyType = familyType,
				displayPlotEntropy = FALSE,
				sizeMaxOutlier = sizeMaxOutlier,
				autoRegroupOutliers = autoRegroupOutliers,
				backend_data = backend_data
			),
			error = function(e) {
				if (!silent) {
					message(
						"Warm shot failed at lambda = ", lambda,
						" (shot ", k, "/", Nshots_warm, "): ",
						e$message
					)
				}
				NULL
			}
		)

		if (is.null(resultList)) {
			next
		}

		if (is.finite(resultList$Hphi) && resultList$Hphi < H) {
			H <- resultList$Hphi
			bestClassif <- resultList
		}

		PHI <- resultList$phi
	}

	bestClassif <- CECfinalize_classif_result(bestClassif)
	CECtag_fit_origin(
		bestClassif,
		init_origin = "warm",
		path_direction = path_direction,
		candidate_group = "warm",
		candidate_index = NA_integer_
	)
}

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
		
		XpredHat 	<- max.col(phiPredM, ties.method = "first")
		
		
		
		
		
		return(XpredHat)
	}

}
 
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
		
		if(length(regroupCols)>length(colFactorPred))
		{
			colNumPred 		<- regroupCols[(1+length(colFactorPred)):length(regroupCols)]
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
		
		XpredHat 	<- max.col(phiPredM, ties.method = "first")
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
				fact_x <- factors[max.col(t(params$discreteProbList[[x]][,colFactorMissing,drop=FALSE]), ties.method = "first")]
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
				# If m and Sigma are known for Z1, Z2, Z3,
				# predict Z3 conditionally on Z1 and Z2.
				# f(z1,z2,z3) = 1/sqrt(2*pi)^3  . 1/sqrt(det(Sigma)) exp( - 1/2  ((z1,z2,23) - m)^t Sigma^(-1) ((z1,z2,23) - m) ) 
				# f(z3|z1,z2 )  =  f(z1,z2,z3)/f(z1,z2) 
				# 				= sqrt(2*pi) sqrt(det(Sigma12))/sqrt(det(Sigma))  exp( - 1/2  ((z1,z2,z3) - m)^t Sigma^(-1) ((z1,z2,23) - m)  + 1/2 ((z1,z2) - m12)^t Sigma12^(-1) ((z1,z2) - m12) )
				# Maximizing this in z3 is equivalent to minimizing the quadratic form in z3.
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
				for(x in 1:r)
				{
					mx  		<- m[x,]

					posx 		<- which(XpredHat==x)

					if(length(posx)>0)
					{
						for (j in 1:length(colNumMissing))
						{
							col_j 		<- colNumMissing[j]

							m_jj		<- mx[col_j]
							predNum[posx,j] 		<- m_jj
						}
					}

				}
			
			}
			ZpredVals[,which(idColToPred%in%colNum)] <- predNum
		
		}
		
		
		
		return(ZpredVals)
		
	}	
}


# Display elapsed and remaining time.
remainingTime <- function(T,indice,nbIndices)
{
	deltaT = difftime(Sys.time(),T,units="mins")
	remaining = deltaT*(nbIndices-indice)/indice
	print(paste0("elapsed = ",round(deltaT)," mins "," remaining = ",round(remaining)," mins "))
	return(list(elapsed= deltaT, remaining = remaining))
}

vectToMat 		<- function(phi,r)
{
	n = length(phi)/r
	return(matrix(phi,n,r,byrow=FALSE))
}

matToVect 		<- function(phiM)
{
	return(as.vector(phiM))
}

vectMatCoords 	<- function(n,r)
{
	return(cbind(matToVect(matrix(1:n,n,r,byrow=FALSE)),matToVect(matrix(1:r,n,r,byrow=TRUE))))						
}

em_mixture_mixed <- function(df, r = 2, max_iter = 100, tol = 1e-6){
  #set.seed(seed)
  n <- nrow(df)
  
  # Split numeric and categorical variables.
  is_fact <- sapply(df, is.factor)
  Xnum <- df[, !is_fact, drop = FALSE]
  Xcat <- df[, is_fact, drop = FALSE]
  
  # Encode categories as integers in 1:k.
  Xcat_enc <- lapply(Xcat, function(col) as.integer(col))
  levels_list <- lapply(Xcat, nlevels)
  Xcat_mat <- as.data.frame(Xcat_enc)
  
  # Random initialization.
  phi <- matrix(runif(n * r), nrow = n)
  phi <- phi / rowSums(phi)
  pi_k <- colMeans(phi)
  
  # Means and variances for numeric variables.
  means <- matrix(0, r, ncol(Xnum))
  vars <- matrix(1, r, ncol(Xnum))  # Avoid zero variance.
  for (k in 1:r){
	means[k, ] <- colSums(phi[, k] * Xnum) / sum(phi[, k])
	vars[k, ] <- apply(Xnum, 2, function(x) sum(phi[, k] * (x - means[k, ])^2)) / sum(phi[, k])
  }
  
  # Multinomial probabilities for categorical variables.
  probs_cat <- lapply(1:ncol(Xcat_mat), function(j){
	matrix(1 / levels_list[[j]], nrow = levels_list[[j]], ncol = r)
  })
  
  loglik_prev <- -Inf
  
  for (iter in 1:max_iter){
	# E-step
	log_phi <- matrix(0, n, r)
	
	for (k in 1:r){
	  # Numeric block: Gaussian density.
	  if (ncol(Xnum) > 0){
		logdens_num <- rowSums(dnorm(as.matrix(Xnum), 
									 mean = matrix(means[k, ], n, ncol(Xnum), byrow = TRUE),
									 sd = sqrt(matrix(vars[k, ], n, ncol(Xnum), byrow = TRUE)),
									 log = TRUE))
	  } else {
		logdens_num <- 0
	  }
	  
	  # Categorical block: product of category probabilities.
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

phiInitR			<- function(n,r)
{
	phi  		<- runif(n*r)
	phiM 		<- vectToMat(phi,r)
	sumByRow 	<- rowSums(phiM)
	phiM2 		<- (1/sumByRow)*phiM  # Divise chaque colonne par le mÃƒÆ’Ã‚Âªme vecteur 
	 
	phi 		<- matToVect(phiM2)
	return(phi)
}

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

phiToNu 		<- function(phi,lambda,C,r,isPhiAlreadyMat = FALSE)
{
	if(!isPhiAlreadyMat)
		phi <- vectToMat(phi=phi,r=r)
		
	return(colSums(phi) / nrow(phi))
}	

getFamilyFunction <- function(prefix, familyType)
{
	get(paste0(prefix, familyType), mode = "function")
}

CECdens	 		<- function(Z,param,lambda=1,applyLog=FALSE)
{
	familyType <- param$familyType
	S <- getFamilyFunction("dens_", familyType)(Z = Z, param = param, lambda = lambda, applyLog = applyLog)
	 
	return(S)
	
}

optPhi			<- function(Z,param,lambda=1)
{
	n 				<- length(Z)
	if(is.matrix(Z)|is.data.frame(Z) )
		n <- dim(Z)[1]
		
	nu 				<- param$nu 
	r 				<- length(nu)
	logG 			<- CECdens(Z=Z,param=param,lambda=lambda,applyLog=TRUE)
	logG$density 	<- rep(log(nu),each=n)  + logG$density / lambda
	densityM 		<- vectToMat(logG$density,r)
	loc				<- max.col(densityM, ties.method = "first")
	phi				<- rep(0,n*r)
	phi[(1:n)+ n*(loc-1) ] <- 1
	
	return(phi)			
}


# Family-specific functions.
{
	# gaussUniv
	{
				
		phiToMean 		<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
				
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			 
			n <- dim(phi)[1] 
			
			return(apply(Z*phi,2,sum)/(n*nuPhi))
		}
	 
		phiToVar 		<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL,mPhi=NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			if(is.null(mPhi))
				mPhi 	<- phiToMean(phi=phi,Z=Z,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE,nuPhi =nuPhi)

			n <- dim(phi)[1]
			varPhi <- rep(0, length(nuPhi))
			pos <- which(nuPhi > 0)

			if (length(pos) > 0) {
				phi_pos <- phi[, pos, drop = FALSE]
				centered <- sweep(
					matrix(Z, nrow = n, ncol = length(pos)),
					2,
					mPhi[pos],
					"-"
				)
				varPhi[pos] <- colSums(phi_pos * centered^2) / (n * nuPhi[pos])
			}

			return(varPhi)
		}
					
		dens_gaussUniv 				<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			states <- param$states
			m <- param$m
			s <- param$s
			C <- param$C
			
			r <- length(m)
			n <- length(Z) 
			
			density <- rep(1 / C, n * r)
			for(i in 1:r)
			{
				state_i <- states[i]
				m_i		<- m[i]
				s_i 	<- s[i]
				idx_i 	<- 1:n + (i-1)*n
				if(applyLog)
				{
					density[idx_i] <- dnorm(Z,mean=m_i,sd=s_i,log=TRUE)
				}else{
					density[idx_i] <- dnorm(Z,mean=m_i,sd=s_i)
				}
				
			}
		 
			return(list(density = density))
			
		}
		
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
			
			# C bounds the maximum density independently of lambda:
			# 1 / (sqrt(2*pi) * s) <= C if and only if s >= 1 / (sqrt(2*pi) * C).
			functionBoundReached <- which(s<1/(sqrt(2*pi)*C))
			s[functionBoundReached] =  1/(sqrt(2*pi)*C) 
			
			params	<- list(states=states,nu=nu,m=m,s=s,lambda=lambda,C=C,familyType ="gaussUniv",phi=phi,functionBoundReached=functionBoundReached)
			
			return(params)
		}

		 
		evalCompositeEntropy_gaussUniv	<-function(phi ,Z,lambda,C,includeClassEntropy = TRUE)
		{
			n 		<- length(Z)
			nr 		<- length(phi)
			r 		<- nr/n
			
			# C bounds the Gaussian density independently of lambda.
			# Since varPhi stores variances, the lower bound is s_min^2.
			s2Min 	<- (1/(sqrt(2*pi)*C))^2
			
			
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
		weighted_cov <- function(res, w) {
		  # Normaliser les poids
		  w <- w / sum(w)
		  
		  # Weighted mean by column.
		  mean_w <- colSums(w * res)
		  
		  # Center the data.
		  res_centered <- sweep(res, 2, mean_w, "-")
		  
		  # Compute the weighted covariance matrix.
		  cov_w <- t(res_centered) %*% (res_centered * w)
		  
		  return(cov_w)
		}
										
		adjust_covariance 				<- function(res, sigmaZmin,w=NULL,sizeMinSigmaComp=3)
		{
		 
		  # Residual matrix.
		  nRes <- dim(res)[1]
		  # estimation empirique
		  
		  if(is.null(w))
		  {
			w =rep(1/nRes,nRes)
		  
		  }
		  
		  if(sum(w)>sizeMinSigmaComp)
		  {
			  
			  Sigma_emp <- weighted_cov(res,w)# on enlÃƒÆ’Ã‚Â¨ve la correction
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
			  
			# Check Sigma_corrected.
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

		CECadjust_covariance_from_matrix <- function(Sigma_emp, sigmaZmin, count, sizeMinSigmaComp = 3) {
			if (count > sizeMinSigmaComp) {
				Sigma_emp[is.na(Sigma_emp)] <- 0
				variances <- diag(Sigma_emp)
				variances_corrected <- pmax(variances, sigmaZmin^2)
				diag(Sigma_emp) <- variances_corrected

				detMin <- prod(sigmaZmin^2)
				detSigmaEmp <- det(Sigma_emp)
				boundReached <- FALSE

				if (detSigmaEmp < detMin) {
					boundReached <- TRUE
					cor_emp <- cov2cor(Sigma_emp)
					D_sqrt <- sqrt(diag(diag(Sigma_emp)))

					eig_R <- eigen(cor_emp, symmetric = TRUE)
					lambda_vals0 <- eig_R$values
					lambda_vals <- lambda_vals0
					vectors <- eig_R$vectors

					lambda_vals[lambda_vals < 0] <- abs(lambda_vals[lambda_vals < 0])

					det_R <- prod(lambda_vals)
					notNull <- which(lambda_vals > 0)
					if (det_R == 0) {
						if (length(notNull) > 0) {
							lambda_vals[lambda_vals == 0] <- min(lambda_vals[notNull]) / 2
						} else {
							lambda_vals[lambda_vals == 0] <- 1
						}
					}

					det_R <- prod(lambda_vals)
					det_Rmin <- detMin / prod(diag(Sigma_emp))

					if (det_R < det_Rmin) {
						lambda_vals <- lambda_vals * (det_Rmin / det_R)^(1 / dim(Sigma_emp)[1])
						cor_emp_corrected <- vectors %*%
							diag(lambda_vals, nrow = dim(Sigma_emp)[1], ncol = dim(Sigma_emp)[1]) %*%
							t(vectors)
						Sigma_corrected <- D_sqrt %*% cor_emp_corrected %*% D_sqrt
					} else {
						Sigma_corrected <- Sigma_emp
					}
				} else {
					Sigma_corrected <- Sigma_emp
				}

				contrNA <- any(is.na(Sigma_corrected))
				if (contrNA) {
					boundReached <- TRUE
					variances <- diag(Sigma_emp)
					varNa <- is.na(variances)
					variances[varNa] <- sigmaZmin[varNa]^2
					variances_corrected <- pmax(variances, sigmaZmin^2)
					Sigma_corrected <- diag(variances, nrow = length(variances), ncol = length(variances))
				} else {
					Sigma_corrected <- (Sigma_corrected + t(Sigma_corrected)) / 2
					if (det(Sigma_corrected) <= 0) {
						boundReached <- TRUE
						variances <- diag(Sigma_emp)
						varNa <- which(is.na(variances))
						variances[varNa] <- sigmaZmin[varNa]^2
						variances_corrected <- pmax(variances, sigmaZmin^2)
						Sigma_corrected <- diag(variances, nrow = length(variances), ncol = length(variances))
					}
				}
			} else {
				boundReached <- TRUE
				Sigma_corrected <- diag(sigmaZmin, ncol = length(sigmaZmin), nrow = length(sigmaZmin))
			}

			list(Sigma_corrected = Sigma_corrected, boundReached = boundReached)
		}

		CECassigned_logdens_from_matrix <- function(logdens_mat, clusters) {
			logdens_mat[cbind(seq_along(clusters), clusters)]
		}

		CECoptParam_fast_from_clusters <- function(backend_data, clusters, lambda = 1, C = 1) {
			clusters <- CECcompress_clusters(clusters)
			n <- length(clusters)
			r <- length(unique(clusters))
			familyType <- backend_data$familyType
			use_cpp <- CECfast_backend_available()

			if (familyType == "gaussVector") {
				if (use_cpp) {
					stats <- CECget_fast_fun("cec_cpp_gaussian_stats")(backend_data$X_num, clusters, r)
					counts <- as.numeric(stats$counts)
					m <- stats$means
					covs <- stats$covs
				} else {
					counts <- as.numeric(tabulate(clusters, nbins = r))
					m <- rowsum(backend_data$X_num, group = clusters, reorder = FALSE) / counts
					covs <- vector("list", r)
					for (i in seq_len(r)) {
						Xi <- backend_data$X_num[clusters == i, , drop = FALSE]
						centered <- sweep(Xi, 2, m[i, ], "-")
						covs[[i]] <- crossprod(centered) / nrow(Xi)
					}
				}
				nu <- counts / n
				if (!is.null(colnames(backend_data$X_num))) {
					colnames(m) <- colnames(backend_data$X_num)
				}

				l <- ncol(backend_data$X_num)
				sigmaZmin <- rep(1 / (sqrt(2 * pi) * C^(1 / l)), l)
				Sigma <- vector("list", r)
				functionBoundReached <- integer()
				for (i in seq_len(r)) {
					sigList <- CECadjust_covariance_from_matrix(
						Sigma_emp = covs[[i]],
						sigmaZmin = sigmaZmin,
						count = counts[i],
						sizeMinSigmaComp = 3
					)
					Sigma[[i]] <- sigList$Sigma_corrected
					if (isTRUE(sigList$boundReached)) {
						functionBoundReached <- c(functionBoundReached, i)
					}
				}

				return(list(
					states = seq_len(r),
					nu = nu,
					m = m,
					Sigma = Sigma,
					lambda = lambda,
					C = C,
					familyType = "gaussVector",
					functionBoundReached = unique(functionBoundReached)
				))
			}

			if (familyType == "discreteVector") {
				discrete <- backend_data$discrete
				if (use_cpp) {
					counts_obj <- CECget_fast_fun("cec_cpp_discrete_counts")(
						discrete$codes,
						as.integer(clusters),
						r,
						discrete$n_levels
					)
					cluster_sizes <- as.numeric(counts_obj$cluster_sizes)
					counts_list <- counts_obj$counts
				} else {
					cluster_sizes <- as.numeric(tabulate(clusters, nbins = r))
					counts_list <- vector("list", r)
					for (i in seq_len(r)) {
						idx_i <- which(clusters == i)
						counts_i <- matrix(0, nrow = discrete$n_levels, ncol = ncol(discrete$codes))
						for (j in seq_len(ncol(discrete$codes))) {
							counts_i[, j] <- tabulate(discrete$codes[idx_i, j], nbins = discrete$n_levels)
						}
						counts_list[[i]] <- counts_i
					}
				}
				nu <- cluster_sizes / n
				C_by_coord <- pmax(10 * discrete$n_unique_by_coord, C)
				discreteProbList <- vector("list", r)
				functionBoundReached <- integer()

				for (i in seq_len(r)) {
					counts_i <- counts_list[[i]]
					discreteProbList[[i]] <- matrix(NA_real_, discrete$n_levels, ncol(discrete$codes))
					for (j in seq_len(ncol(discrete$codes))) {
						valid_levels_j <- if (!is.null(discrete$level_present_by_coord)) {
							discrete$level_present_by_coord[[j]]
						} else {
							tabulate(discrete$codes[, j], nbins = discrete$n_levels) > 0L
						}
						tbl_ij <- counts_i[, j]
						tbl_ij[!valid_levels_j] <- NA_real_
						tbl_ij <- tbl_ij / sum(tbl_ij, na.rm = TRUE)

						toLift <- which(!is.na(tbl_ij) & tbl_ij < 1 / C_by_coord[j])
						toShrink <- which(!is.na(tbl_ij) & tbl_ij >= 1 / C_by_coord[j])

						if (length(toLift) > 0) {
							functionBoundReached <- c(functionBoundReached, i)
						}

						tbl_ij[toLift] <- 1 / C_by_coord[j]
						if (length(toShrink) > 0) {
							tbl_ij[toShrink] <- tbl_ij[toShrink] / sum(tbl_ij[toShrink]) * (1 - sum(tbl_ij[toLift]))
						}
						tbl_ij[is.na(tbl_ij)] <- 0
						discreteProbList[[i]][, j] <- tbl_ij
					}
				}

				return(list(
					states = seq_len(r),
					nu = nu,
					factors = discrete$factors,
					discreteProbList = discreteProbList,
					lambda = lambda,
					C = C,
					familyType = "discreteVector",
					functionBoundReached = unique(functionBoundReached)
				))
			}

			if (familyType == "gaussAndDiscreteVector") {
				params <- list(
					states = seq_len(r),
					nu = tabulate(clusters, nbins = r) / n,
					lambda = lambda,
					C = C,
					familyType = "gaussAndDiscreteVector",
					colFactor = backend_data$colFactor,
					colNum = backend_data$colNum
				)
				functionBoundReached <- integer()

				if (!is.null(backend_data$discrete)) {
					paramsFact <- CECoptParam_fast_from_clusters(
						list(
							familyType = "discreteVector",
							n = n,
							discrete = backend_data$discrete
						),
						clusters = clusters,
						lambda = lambda,
						C = C
					)
					params$factors <- paramsFact$factors
					params$discreteProbList <- paramsFact$discreteProbList
					functionBoundReached <- c(functionBoundReached, paramsFact$functionBoundReached)
				}

				if (!is.null(backend_data$X_num)) {
					paramsNum <- CECoptParam_fast_from_clusters(
						list(
							familyType = "gaussVector",
							n = n,
							X_num = backend_data$X_num
						),
						clusters = clusters,
						lambda = lambda,
						C = C
					)
					params$m <- paramsNum$m
					params$Sigma <- paramsNum$Sigma
					functionBoundReached <- c(functionBoundReached, paramsNum$functionBoundReached)
				}

				params$functionBoundReached <- unique(functionBoundReached)
				return(params)
			}

			stop("Fast parameter update is not implemented for familyType = ", familyType)
		}

		CECcompute_logdens_matrix_fast <- function(backend_data, params, lambda = params$lambda) {
			familyType <- params$familyType
			use_cpp <- CECfast_backend_available()

			if (familyType == "gaussVector") {
				if (use_cpp) {
					return(CECget_fast_fun("cec_cpp_logdens_gaussian")(
						backend_data$X_num,
						params$m,
						params$Sigma,
						lambda
					))
				}
				r <- length(params$states)
				n <- nrow(backend_data$X_num)
				out <- matrix(0, n, r)
				for (i in seq_len(r)) {
					out[, i] <- mvtnorm::dmvnorm(
						backend_data$X_num,
						mean = params$m[i, ],
						sigma = params$Sigma[[i]],
						log = TRUE
					)
				}
				return(out)
			}

			if (familyType == "discreteVector") {
				if (use_cpp) {
					return(CECget_fast_fun("cec_cpp_logdens_discrete")(
						backend_data$discrete$codes,
						params$discreteProbList
					))
				}
				codes <- backend_data$discrete$codes
				n <- nrow(codes)
				r <- length(params$states)
				l <- ncol(codes)
				out <- matrix(0, n, r)
				for (i in seq_len(r)) {
					for (j in seq_len(l)) {
						out[, i] <- out[, i] + log(params$discreteProbList[[i]][codes[, j], j])
					}
				}
				return(out)
			}

			if (familyType == "gaussAndDiscreteVector") {
				n <- backend_data$n
				r <- length(params$nu)
				logdens <- matrix(0, n, r)

				if (!is.null(backend_data$X_num)) {
					logdens <- logdens + CECcompute_logdens_matrix_fast(
						list(
							familyType = "gaussVector",
							n = n,
							X_num = backend_data$X_num
						),
						list(
							states = params$states,
							m = params$m,
							Sigma = params$Sigma,
							lambda = lambda,
							familyType = "gaussVector"
						),
						lambda = lambda
					)
				}

				if (!is.null(backend_data$discrete)) {
					logdens <- logdens + CECcompute_logdens_matrix_fast(
						list(
							familyType = "discreteVector",
							n = n,
							discrete = backend_data$discrete
						),
						list(
							states = params$states,
							discreteProbList = params$discreteProbList,
							familyType = "discreteVector"
						),
						lambda = lambda
					)
				}

				return(logdens)
			}

			stop("Fast density computation is not implemented for familyType = ", familyType)
		}

		CECoptPhi_backend <- function(Z, param, lambda = 1, backend_data = NULL) {
			if (!is.null(backend_data) &&
				isTRUE(backend_data$optimized) &&
				CECis_fast_family(param$familyType)) {
				logdens_mat <- CECcompute_logdens_matrix_fast(backend_data, param, lambda = lambda)
				score_logdens_mat <- logdens_mat / lambda
				if (CECfast_backend_available()) {
					choice <- CECget_fast_fun("cec_cpp_choose_clusters")(score_logdens_mat, param$nu)
					clusters <- as.integer(choice$clusters)
					assigned_logdens <- CECassigned_logdens_from_matrix(logdens_mat, clusters)
				} else {
					score_mat <- sweep(score_logdens_mat, 2, log(param$nu), "+")
					clusters <- max.col(score_mat, ties.method = "first")
					assigned_logdens <- CECassigned_logdens_from_matrix(logdens_mat, clusters)
				}
				return(list(
					phi = CECclusters_to_phi(clusters, length(param$nu)),
					clusters = clusters,
					logdens_mat = logdens_mat,
					assigned_logdens = assigned_logdens
				))
			}

			phi <- optPhi(Z = Z, param = param, lambda = lambda)
			r <- length(param$nu)
			clusters <- max.col(vectToMat(phi, r), ties.method = "first")
			logg <- CECdens(Z = Z, param = param, lambda = lambda, applyLog = TRUE)
			logdens_mat <- vectToMat(logg$density, r)
			list(
				phi = phi,
				clusters = clusters,
				logdens_mat = logdens_mat,
				assigned_logdens = CECassigned_logdens_from_matrix(logdens_mat, clusters)
			)
		}
	
		phiToMeanVec 					<- function(phi,Z,lambda,C,r,isPhiAlreadyMat = FALSE,nuPhi =NULL)
		{
			if(!isPhiAlreadyMat)
				phi <- vectToMat(phi=phi,r=r)
				
			if(is.null(nuPhi))
				nuPhi 	<- phiToNu(phi=phi,lambda=lambda,C=C,r=r,isPhiAlreadyMat = TRUE)
			 
			n <- dim(phi)[1] 
			
			return( t(phi)%*%Z /(n*nuPhi))
		}
	 
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
						
		dens_gaussVector 				<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			if(!is.matrix(Z))
				Z <- as.matrix(Z)
				
			l 		<- dim(Z)[2] 
			n 		<- dim(Z)[1]  
			
		 
			states 	<- param$states 
			
			r 		<- length(states)
			m 		<- param$m    		# matrice de taille r x l 
			Sigma	<- param$Sigma	    # list de taille r contenant pour tout i = 1,...,r, une matrice de corrÃƒÆ’Ã‚Â©lation de taille lxl
			C 		<- param$C 
			density <- rep(1 / C, n * r)
			for(i in 1:r)
			{
				state_i 	<- states[i]
				m_i			<- m[i,]
				Sigma_i 	<- Sigma[[i]]
				idx_i 		<- 1:n + (i-1)*n
				
				if(applyLog)
				{
					density[idx_i] <- mvtnorm::dmvnorm(Z,mean=m_i, sigma = Sigma_i ,log=TRUE)
				}else{
					density[idx_i] <- mvtnorm::dmvnorm(Z,mean=m_i, sigma = Sigma_i)
				}
				
			}
		 
			return(list(density = density))
			
		}
		
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
				idx_i <- 1:n + (i - 1) * n
				H  <- H - (1/lambda)*nuPhi[i]* sum( phiM[,i] *logg$density[idx_i])/sum(phiM[,i])
			}
			
			return(H)
		}

	}
	
	# discreteVector
	{
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
			
			
			 
			density <- rep(NA_real_, n * r)
			for(i in 1:r)
			{
				state_i 	 	<- states[i]
				discreteProb_i 	<- discreteProbList[[i]] 
				idx_i 			<- 1:n + (i-1)*n
				if(applyLog)
				{
					density[idx_i] <- 0 
					for(j in 1:l)
					{
						density[idx_i] <- density[idx_i] + log(discreteProb_i[Z[,j],j] )
						
						
					}
				}else{
					density[idx_i] <- 1
					for(j in 1:l)
					{
						density[idx_i] <- density[idx_i] * discreteProb_i[Z[,j],j] 
						
						
					}
				}
				 
				
				 
				
			}
		 
			return(list(density = density))
			
		}
		
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
			# All coordinates of Z must be factors with aligned levels.
			
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
			 
			
			params	<- list(states=states,nu=nu,factors=factors,discreteProbList=discreteProbList,lambda=lambda,C=C,familyType ="discreteVector",phi=phi,functionBoundReached=functionBoundReached)
			
			return(params)
		}
		
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
				idx_i <- 1:n + (i - 1) * n
				H  <- H - (1/lambda)*nuPhi[i]* sum( phiM[,i] *logg$density[idx_i])/sum(phiM[,i])
			}
			
			return(H)
		}

	
	}

	# gaussAndDiscreteVector
	{
		dens_gaussAndDiscreteVector 		<-  function(Z,param,lambda=1,applyLog=FALSE)
		{
			 
			if(!is.data.frame(Z))
				Z <- as.data.frame(Z)
			 
			
			l 		<- dim(Z)[2] 
			n 		<- dim(Z)[1]  
			
			states 	<- param$states 
			
			r 		<- length(states)
			
			density <- if(applyLog) rep(0, n * r) else rep(1, n * r)
			
			
			Zclass <- vapply(Z, function(x) class(x)[1], character(1))
			
			if(length(which(Zclass=="numeric"))>0)
			{
				SNum 		<- dens_gaussVector(Z=Z[,which(Zclass=="numeric"),drop=FALSE],param=param,lambda=lambda,applyLog=applyLog)
				density 	<- SNum$density  
			}
				
			
			if(length(which(Zclass=="factor"))>0)
			{
				SFac <- dens_discreteVector(Z=Z[,which(Zclass=="factor"),drop=FALSE],param=param,lambda=lambda,applyLog=applyLog)
				
				if(applyLog)
				{
					density = density  + SFac$density
					
				}else{
					density = density  * SFac$density
				}
			}
			 
			 
			return(list(density = density))
			
		}
	
		optParam_gaussAndDiscreteVector		<- function(Z,phi,lambda=1,C=1)
		{
			  
			if(!is.data.frame(Z))
				Z <- as.data.frame(Z)
			 
		
			Zclass <- vapply(Z, function(x) class(x)[1], character(1))
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
			
			Zclass <- vapply(Z, function(x) class(x)[1], character(1))
			
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

optParam 		<- function(Z,phi,lambda=1,C=1,familyType ="gaussAndDiscreteVector")
{		
	params <- getFamilyFunction("optParam_", familyType)(Z = Z, phi = phi, lambda = lambda, C = C)
	return(params)	
}

evalCompositeEntropy <- function(phi ,Z,lambda,C,familyType="gaussAndDiscreteVector")
{
	H <- getFamilyFunction("evalCompositeEntropy_", familyType)(phi = phi, Z = Z, lambda = lambda, C = C)
	return(H)
}
	
CECclassifOneShot_legacy <- function(
	Z,
	lambda = 1,
	C = 1,
	r0 = NULL,
	Nloop = 1000,
	phi0 = NULL,
	familyType = "gaussAndDiscreteVector",
	displayPlotEntropy = FALSE
) {
	n <- length(Z)
	if (is.matrix(Z) | is.data.frame(Z))
		n <- dim(Z)[1]

	if (is.null(phi0)) {
		if (is.null(r0)) {
			r <- floor(log(n) + 1)
		} else {
			r <- r0
		}
		phi <- phiInit(n, r)
	} else {
		nr <- length(phi0)
		r0 <- nr / n
		r <- r0
		phi <- phi0
	}

	phiFingerPrint <- sum(phi * cos(2 * pi * (0:(n * r - 1)) / n))
	H <- c()

	for (i in 1:Nloop) {
		params <- optParam(Z = Z, phi = phi, lambda = lambda, C = C, familyType = familyType)
		if (displayPlotEntropy) {
			H[i] <- evalCompositeEntropy(phi = phi, Z = Z, lambda = lambda, C = C, familyType = familyType)
			plot(H, main = "entropy", xlab = "nLoop", ylab = "H")
		}

		phi <- optPhi(Z = Z, param = params, lambda = lambda)
		nr <- length(phi)
		r <- nr / n

		if (r == 1)
			break()

		phiFingerPrint_i <- sum(phi * cos(2 * pi * (0:(n * r - 1)) / n))
		if (phiFingerPrint_i %in% phiFingerPrint)
			break()

		phiFingerPrint <- c(phiFingerPrint, phiFingerPrint_i)
	}

	params <- optParam(Z = Z, phi = phi, lambda = lambda, C = C, familyType = familyType)
	phi <- optPhi(Z = Z, param = params, lambda = lambda)
	Hphi <- evalCompositeEntropy(phi = phi, Z = Z, lambda = lambda, C = C, familyType = params$familyType)
	list(phi = phi, params = params, Hphi = Hphi)
}

CECclassifOneShot 		<- function(Z,lambda=1,C=1,r0=NULL,Nloop=1000,phi0 = NULL,familyType="gaussAndDiscreteVector",displayPlotEntropy = FALSE,backend_data=NULL)
{
	n <- CECget_n_obs(Z)

	if (is.null(backend_data)) {
		backend_data <- CECprepare_backend_data(Z, familyType = familyType)
	}

	use_fast <- isTRUE(backend_data$optimized) &&
		CECis_fast_family(familyType)

	if (!use_fast) {
		return(
			CECclassifOneShot_legacy(
				Z = Z,
				lambda = lambda,
				C = C,
				r0 = r0,
				Nloop = Nloop,
				phi0 = phi0,
				familyType = familyType,
				displayPlotEntropy = displayPlotEntropy
			)
		)
	}

	if (is.null(phi0)) {
		if (is.null(r0)) {
			r <- floor(log(n) + 1)
		} else {
			r <- r0
		}
		phi_current <- phiInit(n, r)
	} else {
		phi_current <- phi0
		r <- length(phi0) / n
	}

	clusters_current <- CECclusters_from_phi_if_hard(phi_current, n = n, r = r)
	fingerprints <- numeric()
	if (!is.null(clusters_current)) {
		clusters_current <- CECcompress_clusters(clusters_current)
		fingerprints <- CECcluster_fingerprint(clusters_current)
	}

	H_plot <- c()
	for (i in seq_len(Nloop)) {
		if (is.null(clusters_current)) {
			params <- optParam(Z = Z, phi = phi_current, lambda = lambda, C = C, familyType = familyType)
			if (displayPlotEntropy) {
				H_plot[i] <- evalCompositeEntropy(phi = phi_current, Z = Z, lambda = lambda, C = C, familyType = familyType)
				plot(H_plot, main = "entropy", xlab = "nLoop", ylab = "H")
			}
		} else {
			params <- CECoptParam_fast_from_clusters(
				backend_data = backend_data,
				clusters = clusters_current,
				lambda = lambda,
				C = C
			)
			if (displayPlotEntropy) {
				logdens_current <- CECcompute_logdens_matrix_fast(backend_data, params, lambda = lambda)
				assigned_current <- CECassigned_logdens_from_matrix(logdens_current, clusters_current)
				H_plot[i] <- CECentropy_from_assignment(params$nu, assigned_current, lambda)$H_total
				plot(H_plot, main = "entropy", xlab = "nLoop", ylab = "H")
			}
		}

		opt_step <- CECoptPhi_backend(Z = Z, param = params, lambda = lambda, backend_data = backend_data)
		clusters_next <- CECcompress_clusters(opt_step$clusters)
		phi_current <- opt_step$phi

		if (length(unique(clusters_next)) == 1L) {
			clusters_current <- clusters_next
			break
		}

		fp_next <- CECcluster_fingerprint(clusters_next)
		clusters_current <- clusters_next
		if (fp_next %in% fingerprints) {
			break
		}
		fingerprints <- c(fingerprints, fp_next)
	}

	if (is.null(clusters_current)) {
		return(
			CECclassifOneShot_legacy(
				Z = Z,
				lambda = lambda,
				C = C,
				r0 = r0,
				Nloop = Nloop,
				phi0 = phi_current,
				familyType = familyType,
				displayPlotEntropy = displayPlotEntropy
			)
		)
	}

	clusters_current <- CECcompress_clusters(clusters_current)
	params_final <- CECoptParam_fast_from_clusters(
		backend_data = backend_data,
		clusters = clusters_current,
		lambda = lambda,
		C = C
	)
	logdens_final <- CECcompute_logdens_matrix_fast(backend_data, params_final, lambda = lambda)
	assigned_final <- CECassigned_logdens_from_matrix(logdens_final, clusters_current)
	H_final <- CECentropy_from_assignment(params_final$nu, assigned_final, lambda)
	phi_final <- CECclusters_to_phi(clusters_current, length(params_final$states))
	params_final$phi <- phi_final

	list(
		phi = phi_final,
		params = params_final,
		Hphi = H_final$H_total
	)
}


# ============================================================
# ============================================================
# Lambda diagnostics for CEClust
# Main workflow kept in this file: linked lambda paths
# ============================================================

# ------------------------------------------------------------
# Distance between partitions (match by default, or ARI).
# ------------------------------------------------------------
partition_distance <- function(x, y, method = c("match", "ARI")) {
  
  method <- match.arg(method)
  
  x <- as.integer(as.factor(x))
  y <- as.integer(as.factor(y))
  
  if (length(x) != length(y)) stop("x and y must have the same length.")
  n <- length(x)
  if (n <= 1) return(0)
  
  # ----------------------------------------------------------
  # Method 1: optimal matching distance (default).
  # ----------------------------------------------------------
  if (method == "match") {
    
    r <- length(unique(x))
    s <- length(unique(y))
    # Enforce r <= s.
    
    # on impose r <= s
    if (r > s) return(partition_distance(y, x, method = "match"))
    
    tab <- table(x, y)
    tab <- as.matrix(tab)
    
    # Pad to a square matrix when needed.
    if (r < s) {
      tab <- rbind(tab, matrix(0, nrow = s - r, ncol = s))
    }
    
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Package 'clue' required for match distance. Install it.")
    }
    
    assignment <- clue::solve_LSAP(tab, maximum = TRUE)
    matched <- sum(tab[cbind(seq_len(nrow(tab)), assignment)])
    
    return(1 - matched / n)
  }
  
  # ----------------------------------------------------------
  # Method 2: adjusted Rand index.
  # ----------------------------------------------------------
  if (method == "ARI") {
    
    tab <- table(x, y)
    a <- rowSums(tab)
    b <- colSums(tab)
    
    choose2 <- function(z) {
      z <- as.numeric(z)
      z * (z - 1) / 2
    }
    
    sum_nij <- sum(choose2(tab))
    sum_ai  <- sum(choose2(a))
    sum_bj  <- sum(choose2(b))
    total   <- choose2(n)
    
    if (total == 0) return(0)
    
    expected <- (sum_ai * sum_bj) / total
    max_ind  <- 0.5 * (sum_ai + sum_bj)
    denom    <- max_ind - expected
    
    if (denom == 0) return(0)
    
    ari <- (sum_nij - expected) / denom
    return(1 - ari)
  }
}


# ------------------------------------------------------------
# Partition similarity = 1 - distance.
# ------------------------------------------------------------
partition_similarity <- function(x, y, method = c("match", "ARI")) {
  method <- match.arg(method)
  1 - partition_distance(x, y, method = method)
}


# ------------------------------------------------------------
# Bootstrap rows
# ------------------------------------------------------------
bootstrap_rows <- function(Z, replace = TRUE) {
  n <- if (is.matrix(Z) || is.data.frame(Z)) nrow(Z) else length(Z)
  sample.int(n, size = n, replace = replace)
}

CECrun_fixed_lambda_path_task <- function(
  rep_idx,
  Z,
  lambda_grid,
  C,
  r0,
  Nshots_fresh,
  Nshots_warm,
  Nloop,
  familyType,
  sizeMaxOutlier,
  autoRegroupOutliers,
  focus,
  silent,
  backend_data = NULL
) {
  dir_rep <- if (rep_idx %% 2 == 1) "forward" else "backward"
  if (is.null(backend_data)) {
    backend_data <- CECprepare_backend_data(Z, familyType = familyType)
  }
  path_rep <- CECfollowLambdaPath(
    Z = Z,
    lambda_grid = lambda_grid,
    direction = dir_rep,
    C = C,
    r0 = r0,
    Nshots_fresh = Nshots_fresh,
    Nshots_warm = Nshots_warm,
    Nloop = Nloop,
    familyType = familyType,
    sizeMaxOutlier = sizeMaxOutlier,
    autoRegroupOutliers = autoRegroupOutliers,
    focus = focus,
    silent = silent,
    backend_data = backend_data
  )

  per_lambda <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    fit_i <- path_rep$fits[[i]]
    if (!is.null(fit_i)) {
      dec_i <- CECdecompose_fit_on_data(fit_i, Z = Z)
      proj_i <- list(
        clusters = fit_i$clusters,
        K = sum(fit_i$params$nu > 0),
        H = dec_i$H_total,
        H_class = dec_i$H_class,
        H_cond = dec_i$H_cond
      )
      per_lambda[[i]] <- list(
        fit = fit_i,
        proj = proj_i,
        part = proj_i$clusters,
        H = proj_i$H,
        Hraw = fit_i$Hphi,
        REO = fit_i$REO,
        K = proj_i$K
      )
    } else {
      per_lambda[[i]] <- list(
        fit = NULL,
        proj = NULL,
        part = NULL,
        H = NA_real_,
        Hraw = NA_real_,
        REO = NA_integer_,
        K = NA_integer_
      )
    }
  }

  list(rep = rep_idx, direction = dir_rep, per_lambda = per_lambda)
}

CECrun_bootstrap_lambda_path_task <- function(
  b_idx,
  Z,
  lambda_grid,
  C,
  r0,
  Nshots_fresh,
  Nshots_warm,
  Nloop,
  familyType,
  sizeMaxOutlier,
  autoRegroupOutliers,
  focus,
  replace,
  silent
) {
  dir_b <- if (b_idx %% 2 == 1) "forward" else "backward"
  n <- if (is.matrix(Z) || is.data.frame(Z)) nrow(Z) else length(Z)
  ib <- sample.int(n, size = n, replace = replace)
  Zb <- if (is.matrix(Z) || is.data.frame(Z)) Z[ib, , drop = FALSE] else Z[ib]

  path_b <- CECfollowLambdaPath(
    Z = Zb,
    lambda_grid = lambda_grid,
    direction = dir_b,
    C = C,
    r0 = r0,
    Nshots_fresh = Nshots_fresh,
    Nshots_warm = Nshots_warm,
    Nloop = Nloop,
    familyType = familyType,
    sizeMaxOutlier = sizeMaxOutlier,
    autoRegroupOutliers = autoRegroupOutliers,
    focus = focus,
    silent = silent,
    backend_data = CECprepare_backend_data(Zb, familyType = familyType)
  )

  per_lambda <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    fit_i <- path_b$fits[[i]]
    if (!is.null(fit_i)) {
      proj_i <- CECproject_fit_on_reference(
        fit = fit_i,
        Zref = Z,
        lambda = lambda_grid[i]
      )
      per_lambda[[i]] <- list(
        fit = fit_i,
        proj = proj_i,
        part = proj_i$clusters,
        H = proj_i$H,
        Hraw = fit_i$Hphi,
        REO = fit_i$REO,
        K = proj_i$K
      )
    } else {
      per_lambda[[i]] <- list(
        fit = NULL,
        proj = NULL,
        part = NULL,
        H = NA_real_,
        Hraw = NA_real_,
        REO = NA_integer_,
        K = NA_integer_
      )
    }
  }

  list(b = b_idx, direction = dir_b, boot_indices = ib, per_lambda = per_lambda)
}

CECmake_diagnostics_cluster <- function(
  n_cores,
  seed = NULL,
  context = NULL,
  backend_enabled = !isTRUE(.CEC_fast_backend$force_disable)
) {
  cl <- parallel::makePSOCKcluster(n_cores)
  package_root <- CECget_package_root()
  lib_paths <- .libPaths()
  parallel::clusterCall(
    cl,
    function(pkg_root, lib_paths_i, ctx, backend_flag) {
      package_name <- "CEClust"
      .libPaths(unique(c(lib_paths_i, .libPaths())))
      options(CEClust.package_root = pkg_root)
      load_ok <- tryCatch(
        {
          loadNamespace(package_name)
          TRUE
        },
        error = function(e) FALSE
      )

      if (!isTRUE(load_ok)) {
        if (requireNamespace("pkgload", quietly = TRUE) && file.exists(file.path(pkg_root, "DESCRIPTION"))) {
          pkgload::load_all(
            path = pkg_root,
            export_all = TRUE,
            helpers = FALSE,
            quiet = TRUE,
            reset = FALSE
          )
        } else {
          stop(
            "Unable to load the CEClust package in a worker. ",
            "Install the package or install 'pkgload' for development mode."
          )
        }
      }

      ns <- asNamespace(package_name)
      fast_backend <- get(".CEC_fast_backend", envir = ns, inherits = FALSE)
      fast_backend$force_disable <- !isTRUE(backend_flag)
      fast_backend$available <- NULL
      assign(".CEC_diagnostics_context", ctx, envir = .GlobalEnv)
      NULL
    },
    package_root,
    lib_paths,
    context,
    backend_enabled
  )
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }
  cl
}

CECsimple_hash <- function(x) {
  raw <- serialize(x, NULL, version = 2)
  ints <- as.integer(raw)
  weights <- (seq_along(ints) %% 251L) + 1L
  hash_val <- sum(ints * weights) %% 2147483647
  sprintf("%010d", as.integer(hash_val))
}

CECquick_data_signature <- function(Z) {
  if (is.data.frame(Z) || is.matrix(Z)) {
    nr <- nrow(Z)
    nc <- ncol(Z)
    row_idx <- if (nr > 0) unique(c(seq_len(min(2L, nr)), seq.int(max(1L, nr - 1L), nr))) else integer(0)
    col_idx <- if (nc > 0) unique(c(seq_len(min(3L, nc)), seq.int(max(1L, nc - 2L), nc))) else integer(0)
    sample_vals <- if (length(row_idx) > 0L && length(col_idx) > 0L) {
      as.character(unlist(as.data.frame(Z[row_idx, col_idx, drop = FALSE]), use.names = FALSE))
    } else {
      character(0)
    }
    classes <- if (is.data.frame(Z)) {
      vapply(Z, function(x) class(x)[1], character(1))
    } else {
      rep(class(Z)[1], nc)
    }
    return(list(
      type = "rectangular",
      dim = c(nr, nc),
      names = colnames(Z),
      classes = classes,
      sample = sample_vals
    ))
  }

  n <- length(Z)
  idx <- if (n > 0) unique(c(seq_len(min(5L, n)), seq.int(max(1L, n - 4L), n))) else integer(0)
  list(
    type = "vector",
    length = n,
    class = class(Z)[1],
    sample = if (length(idx) > 0L) as.character(Z[idx]) else character(0)
  )
}

CECformat_elapsed_time <- function(start_time, end_time = Sys.time()) {
  elapsed <- max(0, as.numeric(difftime(end_time, start_time, units = "secs")))
  hours <- floor(elapsed / 3600)
  minutes <- floor((elapsed %% 3600) / 60)
  seconds <- floor(elapsed %% 60)
  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

CECtemp_checkpoint_dir <- function(prefix = "cecdiagnose_lambda_grid_linked", signature = NULL) {
  root <- Sys.getenv("TEMP", unset = "")
  if (!nzchar(root)) {
    root <- tempdir()
  }
  suffix <- if (is.null(signature) || !nzchar(signature)) {
    paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  } else {
    signature
  }
  path <- file.path(root, "CEClust_checkpoints", paste0(prefix, "_", suffix))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

CECwrite_rds_atomic <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid(),
    "_",
    as.integer(stats::runif(1, min = 1, max = 1e9))
  )
  saveRDS(object, file = tmp_path, compress = compress)
  if (file.exists(path)) {
    file.remove(path)
  }
  if (!file.rename(tmp_path, path)) {
    ok <- file.copy(tmp_path, path, overwrite = TRUE)
    unlink(tmp_path, force = TRUE)
    if (!ok) {
      stop("Could not atomically write checkpoint file: ", path)
    }
  }
  invisible(path)
}

CECdiagnostics_task_key <- function(task) {
  paste(task$type, task$idx, sep = ":")
}

CECdiagnostics_task_label <- function(task) {
  paste(task$type, task$idx, task$direction)
}

CECdiagnostics_build_tasks <- function(
  k0,
  B,
  checkpoint_dir = NULL,
  checkpoint_compress = FALSE,
  seed = NULL
) {
  fixed_dir <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, "fixed") else NULL
  boot_dir <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, "boot") else NULL

  make_task <- function(type, idx, offset = 0L) {
    direction <- if (idx %% 2 == 1L) "forward" else "backward"
    result_file <- NULL
    if (!is.null(checkpoint_dir)) {
      subdir <- if (identical(type, "fixed")) fixed_dir else boot_dir
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      result_file <- file.path(
        subdir,
        sprintf("%s_%03d_%s.rds", type, idx, direction)
      )
    }

    list(
      type = type,
      idx = idx,
      direction = direction,
      result_file = result_file,
      checkpoint_compress = checkpoint_compress,
      seed = if (is.null(seed)) NULL else as.integer(seed + offset + idx - 1L)
    )
  }

  fixed_tasks <- if (k0 > 0) {
    lapply(seq_len(k0), function(idx) make_task("fixed", idx, offset = 0L))
  } else {
    list()
  }
  boot_tasks <- if (B > 0) {
    lapply(seq_len(B), function(idx) make_task("boot", idx, offset = as.integer(k0)))
  } else {
    list()
  }

  list(
    fixed = fixed_tasks,
    boot = boot_tasks,
    all = c(fixed_tasks, boot_tasks)
  )
}

CECdiagnostics_make_task_meta <- function(
  Z,
  lambda_grid,
  k0,
  B,
  C,
  r0,
  Nshots_fresh,
  Nshots_warm,
  Nloop,
  familyType,
  sizeMaxOutlier,
  autoRegroupOutliers,
  focus,
  replace,
  seed
) {
  payload <- list(
    data_signature = CECquick_data_signature(Z),
    lambda_grid = lambda_grid,
    k0 = k0,
    B = B,
    C = C,
    r0 = r0,
    Nshots_fresh = Nshots_fresh,
    Nshots_warm = Nshots_warm,
    Nloop = Nloop,
    familyType = familyType,
    sizeMaxOutlier = sizeMaxOutlier,
    autoRegroupOutliers = autoRegroupOutliers,
    focus = focus,
    replace = replace,
    seed = seed
  )

  list(
    version = 1L,
    created_at = as.character(Sys.time()),
    payload = payload,
    signature = CECsimple_hash(payload)
  )
}

CECdiagnostics_prepare_checkpoint <- function(
  checkpoint_dir,
  task_meta,
  auto_checkpoint = FALSE,
  resume = TRUE
) {
  if (isFALSE(checkpoint_dir) || identical(checkpoint_dir, "none")) {
    return(list(enabled = FALSE, dir = NULL, auto = FALSE))
  }

  if (isTRUE(checkpoint_dir) || identical(checkpoint_dir, "temp")) {
    checkpoint_dir <- CECtemp_checkpoint_dir(signature = task_meta$signature)
    auto_checkpoint <- TRUE
  } else if (is.null(checkpoint_dir)) {
    if (!isTRUE(auto_checkpoint)) {
      return(list(enabled = FALSE, dir = NULL, auto = FALSE))
    }
    checkpoint_dir <- CECtemp_checkpoint_dir(signature = task_meta$signature)
  } else {
    checkpoint_dir <- normalizePath(path.expand(as.character(checkpoint_dir)[1]), winslash = "/", mustWork = FALSE)
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(checkpoint_dir, "fixed"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(checkpoint_dir, "boot"), recursive = TRUE, showWarnings = FALSE)

  meta_path <- file.path(checkpoint_dir, "meta.rds")
  if (file.exists(meta_path)) {
    saved_meta <- readRDS(meta_path)
    saved_sig <- saved_meta$signature
    current_sig <- task_meta$signature
    if (!identical(saved_sig, current_sig)) {
      stop(
        "Existing checkpoint directory does not match the current run: ",
        checkpoint_dir,
        ". Choose another checkpoint_dir or remove the old checkpoints."
      )
    }
  } else {
    CECwrite_rds_atomic(task_meta, meta_path, compress = FALSE)
  }

  list(
    enabled = TRUE,
    dir = checkpoint_dir,
    auto = auto_checkpoint,
    meta_path = meta_path,
    resume = resume
  )
}

CECdiagnostics_completed_task_keys <- function(tasks) {
  keys <- character(0)
  for (task in tasks) {
    if (!is.null(task$result_file) && file.exists(task$result_file)) {
      keys <- c(keys, CECdiagnostics_task_key(task))
    }
  }
  keys
}

CECdiagnostics_partition_lengths_in_result <- function(res) {
  if (is.null(res) || is.null(res$per_lambda)) {
    return(integer(0))
  }

  lens <- vapply(
    res$per_lambda,
    function(slot_i) {
      if (is.null(slot_i) || is.null(slot_i$part)) {
        return(NA_integer_)
      }
      length(slot_i$part)
    },
    integer(1)
  )

  unique(stats::na.omit(lens))
}

CECdiagnostics_validate_result_lengths <- function(results, expected_n, label = "checkpoint") {
  expected_n <- as.integer(expected_n)
  bad <- character(0)

  for (i in seq_along(results)) {
    res <- results[[i]]
    if (is.null(res)) {
      next
    }

    lens <- CECdiagnostics_partition_lengths_in_result(res)
    if (length(lens) > 0L && any(lens != expected_n)) {
      result_name <- names(results)[i]
      if (is.null(result_name) || is.na(result_name) || !nzchar(result_name)) {
        result_name <- paste0(label, " result #", i)
      }
      bad <- c(bad, paste0(result_name, " has partition length(s) ", paste(lens, collapse = ", ")))
    }
  }

  if (length(bad) > 0L) {
    stop(
      "Existing diagnostic checkpoint results are incompatible with the current data length (n = ",
      expected_n,
      "): ",
      paste(head(bad, 5L), collapse = "; "),
      ". Choose another checkpoint_dir or remove the old checkpoints."
    )
  }

  invisible(TRUE)
}

CECdiagnostics_new_progress_state <- function(k0, B, show_progress = TRUE) {
  state <- new.env(parent = emptyenv())
  state$show_progress <- isTRUE(show_progress)
  state$start_time <- Sys.time()
  state$total_fixed <- as.integer(k0)
  state$total_boot <- as.integer(B)
  state$done_fixed <- 0L
  state$done_boot <- 0L
  state
}

CECdiagnostics_print_progress <- function(progress_state, prefix = "CECdiagnose progress", task = NULL) {
  if (!isTRUE(progress_state$show_progress)) {
    return(invisible(NULL))
  }

  done_total <- progress_state$done_fixed + progress_state$done_boot
  total <- progress_state$total_fixed + progress_state$total_boot
  task_txt <- if (is.null(task)) "" else paste0(" | last=", CECdiagnostics_task_label(task))
  cat(
    prefix,
    ": fixed ", progress_state$done_fixed, "/", progress_state$total_fixed,
    " | boot ", progress_state$done_boot, "/", progress_state$total_boot,
    " | total ", done_total, "/", total,
    task_txt,
    " | elapsed ", CECformat_elapsed_time(progress_state$start_time),
    "\n",
    sep = ""
  )
  invisible(NULL)
}

CECdiagnostics_set_progress_from_keys <- function(progress_state, completed_keys, tasks) {
  fixed_done <- 0L
  boot_done <- 0L
  if (length(completed_keys) > 0L) {
    completed_keys <- unique(completed_keys)
    for (task in tasks) {
      if (CECdiagnostics_task_key(task) %in% completed_keys) {
        if (identical(task$type, "fixed")) {
          fixed_done <- fixed_done + 1L
        } else {
          boot_done <- boot_done + 1L
        }
      }
    }
  }
  progress_state$done_fixed <- fixed_done
  progress_state$done_boot <- boot_done
  invisible(progress_state)
}

CECdiagnostics_mark_task_completed <- function(progress_state, task) {
  if (identical(task$type, "fixed")) {
    progress_state$done_fixed <- progress_state$done_fixed + 1L
  } else {
    progress_state$done_boot <- progress_state$done_boot + 1L
  }
  CECdiagnostics_print_progress(progress_state, task = task)
  invisible(progress_state)
}

CECrun_diagnostic_task_worker <- function(task, context = NULL) {
  if (is.null(context)) {
    context <- get(".CEC_diagnostics_context", envir = .GlobalEnv, inherits = TRUE)
  }

  if (!is.null(task$result_file) && file.exists(task$result_file)) {
    return(list(
      ok = TRUE,
      task = task,
      stored = TRUE,
      from_checkpoint = TRUE
    ))
  }

  if (!is.null(task$seed)) {
    set.seed(task$seed)
  }

  out <- tryCatch(
    {
      if (identical(task$type, "fixed")) {
        CECrun_fixed_lambda_path_task(
          rep_idx = task$idx,
          Z = context$Z,
          lambda_grid = context$lambda_grid,
          C = context$C,
          r0 = context$r0,
          Nshots_fresh = context$Nshots_fresh,
          Nshots_warm = context$Nshots_warm,
          Nloop = context$Nloop,
          familyType = context$familyType,
          sizeMaxOutlier = context$sizeMaxOutlier,
          autoRegroupOutliers = context$autoRegroupOutliers,
          focus = context$focus,
          silent = context$silent,
          backend_data = context$fixed_backend_data
        )
      } else {
        CECrun_bootstrap_lambda_path_task(
          b_idx = task$idx,
          Z = context$Z,
          lambda_grid = context$lambda_grid,
          C = context$C,
          r0 = context$r0,
          Nshots_fresh = context$Nshots_fresh,
          Nshots_warm = context$Nshots_warm,
          Nloop = context$Nloop,
          familyType = context$familyType,
          sizeMaxOutlier = context$sizeMaxOutlier,
          autoRegroupOutliers = context$autoRegroupOutliers,
          focus = context$focus,
          replace = context$replace,
          silent = context$silent
        )
      }
    },
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "CECdiagnose_task_error"
      )
    }
  )

  if (inherits(out, "CECdiagnose_task_error")) {
    return(list(
      ok = FALSE,
      task = task,
      error = out$message
    ))
  }

  if (!is.null(task$result_file)) {
    CECwrite_rds_atomic(out, task$result_file, compress = isTRUE(task$checkpoint_compress))
    return(list(
      ok = TRUE,
      task = task,
      stored = TRUE,
      from_checkpoint = FALSE
    ))
  }

  list(
    ok = TRUE,
    task = task,
    stored = FALSE,
    result = out
  )
}

CECrun_diagnostic_task_local <- function(task, context) {
  CECrun_diagnostic_task_worker(task = task, context = context)
}

CECrun_diagnostic_tasks_sequential <- function(
  tasks,
  context,
  progress_state,
  completed_keys = character(0)
) {
  results <- list()
  completed_keys <- unique(completed_keys)

  for (task in tasks) {
    key <- CECdiagnostics_task_key(task)
    if (key %in% completed_keys) {
      next
    }

    res <- CECrun_diagnostic_task_local(task, context = context)
    if (!isTRUE(res$ok)) {
      stop("Task failed (", CECdiagnostics_task_label(task), "): ", res$error)
    }
    results[[key]] <- res
    completed_keys <- c(completed_keys, key)
    CECdiagnostics_mark_task_completed(progress_state, task)
  }

  list(results = results, completed_keys = unique(completed_keys))
}

CECrun_diagnostic_tasks_parallel_batch <- function(
  tasks,
  context,
  n_cores,
  progress_state
) {
  if (length(tasks) == 0L) {
    return(list(results = list(), completed_keys = character(0), failed = character(0)))
  }

  cl <- CECmake_diagnostics_cluster(
    n_cores = min(n_cores, length(tasks)),
    seed = NULL,
    context = context,
    backend_enabled = !isTRUE(.CEC_fast_backend$force_disable)
  )
  on.exit(parallel::stopCluster(cl), add = TRUE)

  results <- vector("list", length(tasks))
  names(results) <- vapply(tasks, CECdiagnostics_task_key, character(1))
  completed_keys <- character(0)
  send_call <- get("sendCall", envir = asNamespace("parallel"), inherits = FALSE)
  recv_one_result <- get("recvOneResult", envir = asNamespace("parallel"), inherits = FALSE)

  next_idx <- 1L
  active <- 0L
  for (worker_idx in seq_len(min(length(cl), length(tasks)))) {
    send_call(
      cl[[worker_idx]],
      function(task) {
        get("CECrun_diagnostic_task_worker", envir = asNamespace("CEClust"), inherits = FALSE)(task)
      },
      list(tasks[[next_idx]]),
      tag = next_idx
    )
    next_idx <- next_idx + 1L
    active <- active + 1L
  }

  while (active > 0L) {
    recv <- recv_one_result(cl)
    task_idx <- as.integer(recv$tag)
    task <- tasks[[task_idx]]
    key <- CECdiagnostics_task_key(task)
    results[[key]] <- recv$value

    if (isTRUE(recv$value$ok)) {
      completed_keys <- c(completed_keys, key)
      CECdiagnostics_mark_task_completed(progress_state, task)
    } else {
      stop("Task failed (", CECdiagnostics_task_label(task), "): ", recv$value$error)
    }

    active <- active - 1L
    if (next_idx <= length(tasks)) {
      send_call(
        recv$node,
        function(task) {
          get("CECrun_diagnostic_task_worker", envir = asNamespace("CEClust"), inherits = FALSE)(task)
        },
        list(tasks[[next_idx]]),
        tag = next_idx
      )
      next_idx <- next_idx + 1L
      active <- active + 1L
    }
  }

  list(results = results, completed_keys = unique(completed_keys), failed = character(0))
}

# ------------------------------------------------------------
# Safe wrapper around CECclassif.
# ------------------------------------------------------------
CECfit_safe <- function(
  Z,
  lambda,
  C = 1,
  r0 = NULL,
  Nshots = 100,
  Nloop = 1000,
  familyType = "gaussAndDiscreteVector",
  sizeMaxOutlier = 0,
  autoRegroupOutliers = FALSE,
  displayRemainingTime = FALSE,
  focus = NULL,
  silent = TRUE,
  backend_data = NULL
) {
  if (is.null(backend_data)) {
    backend_data <- CECprepare_backend_data(Z, familyType = familyType)
  }
  out <- tryCatch(
    {
      CECclassif(
        Z = Z,
        lambda = lambda,
        C = C,
        r0 = r0,
        Nshots = Nshots,
        Nloop = Nloop,
        familyType = familyType,
        sizeMaxOutlier = sizeMaxOutlier,
        autoRegroupOutliers = autoRegroupOutliers,
        displayRemainingTime = displayRemainingTime,
        focus = focus,
        backend_data = backend_data
      )
    },
    error = function(e) {
      if (!silent) message("CECclassif failed at lambda = ", lambda, ": ", e$message)
      NULL
    }
  )
  out
}

# ------------------------------------------------------------
# Warm start helper used when following a lambda path
#
# Reuses a previous partition matrix as initialization.
#
# ------------------------------------------------------------
CECfit_warm <- function(
  Z, lambda, phi0 = NULL,
  C = 1, r0 = NULL,
  Nshots = 3, Nloop = 30,
  familyType = "gaussUniv"
) {
  if (is.null(phi0)) {
    return(
      CECclassif(
        Z = Z, lambda = lambda, C = C, r0 = r0,
        Nshots = Nshots, Nloop = Nloop,
        familyType = familyType
      )
    )
  }

  one <- CECclassifOneShot(
    Z = Z,
    lambda = lambda,
    C = C,
    r0 = r0,
    Nloop = Nloop,
    phi0 = phi0,
    familyType = familyType,
    displayPlotEntropy = FALSE
  )

  REO <- length(one$params$states)
  phiM <- vectToMat(one$phi, REO)
  one$clusters <- max.col(phiM, ties.method = "first")
  one$REO <- REO
  one
}

CECproject_fit_on_reference <- function(fit, Zref, lambda = NULL) {
  if (is.null(fit)) return(NULL)

  params <- fit$params
  if (is.null(lambda)) lambda <- params$lambda

  phi_ref <- optPhi(Z = Zref, param = params, lambda = lambda)
  r_ref <- length(params$states)
  phiM_ref <- vectToMat(phi_ref, r_ref)
  clusters_ref <- max.col(phiM_ref, ties.method = "first")
  nu_ref <- colMeans(phiM_ref)
  K_ref <- sum(nu_ref > 0)

  dec_ref <- evalCompositeEntropyDecomposed(
    phi = phi_ref,
    Z = Zref,
    lambda = lambda,
    C = params$C,
    familyType = params$familyType
  )

  list(
    phi = phi_ref,
    phiM = phiM_ref,
    clusters = clusters_ref,
    nu = nu_ref,
    K = K_ref,
    H = dec_ref$H_total,
    H_class = dec_ref$H_class,
    H_cond = dec_ref$H_cond
  )
}


# ------------------------------------------------------------
# Pairwise ARI stability summary for a list of partitions
# ------------------------------------------------------------
# ------------------------------------------------------------
# Pairwise partition stability summary
# similarity = 1 - partition_distance
# ------------------------------------------------------------
partition_stability_summary <- function(partitions, partition_metric = c("match", "ARI")) {
  partition_metric <- match.arg(partition_metric)
  
  ok <- vapply(partitions, function(x) !is.null(x), logical(1))
  parts <- partitions[ok]
  m <- length(parts)

  if (m == 0) {
    return(list(
      n_valid = 0L,
      sim_matrix = matrix(NA_real_, 0, 0),
      mean_sim = NA_real_,
      sd_sim = NA_real_,
      median_sim = NA_real_,
      min_sim = NA_real_,
      max_sim = NA_real_,
      medoid_index = NA_integer_,
      medoid_mean_sim = NA_real_,
      metric = partition_metric
    ))
  }

  if (m == 1) {
    return(list(
      n_valid = 1L,
      sim_matrix = matrix(1, 1, 1),
      mean_sim = 1,
      sd_sim = 0,
      median_sim = 1,
      min_sim = 1,
      max_sim = 1,
      medoid_index = 1L,
      medoid_mean_sim = 1,
      metric = partition_metric
    ))
  }

  S <- matrix(NA_real_, m, m)
  diag(S) <- 1

  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      S[i, j] <- partition_similarity(parts[[i]], parts[[j]], method = partition_metric)
      S[j, i] <- S[i, j]
    }
  }

  upper_vals <- S[upper.tri(S)]
  row_mean_sim <- rowMeans(S, na.rm = TRUE)
  medoid_index <- which.max(row_mean_sim)

  list(
    n_valid = m,
    sim_matrix = S,
    mean_sim = mean(upper_vals, na.rm = TRUE),
    sd_sim = sd(upper_vals, na.rm = TRUE),
    median_sim = median(upper_vals, na.rm = TRUE),
    min_sim = min(upper_vals, na.rm = TRUE),
    max_sim = max(upper_vals, na.rm = TRUE),
    medoid_index = medoid_index,
    medoid_mean_sim = row_mean_sim[medoid_index],
    metric = partition_metric
  )
}


# ------------------------------------------------------------
# Stability summary relative to a reference partition
# The reference is typically the partition with smallest H
# among the projected partitions on Z.
# ------------------------------------------------------------
partition_reference_stability_summary <- function(
  partitions,
  H,
  partition_metric = c("match", "ARI"),
  ref_index = NULL
) {
  partition_metric <- match.arg(partition_metric)
  
  ok_part <- vapply(partitions, function(x) !is.null(x), logical(1))
  ok_H <- !is.na(H)
  ok <- ok_part & ok_H
  
  parts <- partitions[ok]
  H_ok <- H[ok]
  
  m <- length(parts)
  
  if (m == 0) {
    return(list(
      n_valid = 0L,
      ref_index_valid = NA_integer_,
      ref_index_original = NA_integer_,
      ref_H = NA_real_,
      similarities = numeric(0),
      mean_sim = NA_real_,
      sd_sim = NA_real_,
      median_sim = NA_real_,
      min_sim = NA_real_,
      max_sim = NA_real_,
      metric = partition_metric
    ))
  }
  
  if (is.null(ref_index)) {
    ref_index_valid <- which.min(H_ok)
  } else {
    idx_map <- which(ok)
    pos <- match(ref_index, idx_map)
    if (is.na(pos)) stop("ref_index does not correspond to a valid partition/H pair.")
    ref_index_valid <- pos
  }
  
  ref_partition <- parts[[ref_index_valid]]
  sims <- vapply(
    parts,
    function(p) partition_similarity(ref_partition, p, method = partition_metric),
    numeric(1)
  )
  
  idx_map <- which(ok)
  
  list(
    n_valid = m,
    ref_index_valid = ref_index_valid,
    ref_index_original = idx_map[ref_index_valid],
    ref_H = H_ok[ref_index_valid],
    similarities = sims,
    mean_sim = mean(sims),
    sd_sim = if (length(sims) > 1) sd(sims) else 0,
    median_sim = median(sims),
    min_sim = min(sims),
    max_sim = max(sims),
    metric = partition_metric
  )
}


# ------------------------------------------------------------
# Numeric summary helper
# ------------------------------------------------------------
numeric_summary <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(list(
      n = 0L,
      mean = NA_real_,
      sd = NA_real_,
      min = NA_real_,
      q25 = NA_real_,
      median = NA_real_,
      q75 = NA_real_,
      max = NA_real_
    ))
  }

  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)

  list(
    n = length(x),
    mean = mean(x),
    sd = if (length(x) > 1) sd(x) else 0,
    min = min(x),
    q25 = qs[1],
    median = qs[2],
    q75 = qs[3],
    max = max(x)
  )
}

.safe_ratio <- function(num, den) {
  if (is.null(num) || is.null(den) || is.na(num) || is.na(den) || den == 0) {
    return(NA_real_)
  }
  num / den
}

# ------------------------------------------------------------
# Summary row helper for one lambda result
# ------------------------------------------------------------
CECdiagnose_lambda_summary_row <- function(obj) {
  sH0 <- numeric_summary(obj$H0)
  sHB <- numeric_summary(obj$HB)
  sH0raw <- numeric_summary(obj$H0_raw)
  sHBraw <- numeric_summary(obj$HB_raw)
  origin0 <- CECcount_fit_origins(obj$fits0)
  originB <- CECcount_fit_origins(obj$fitsB)
  best0_origin <- CECbest_fit_origin_info(obj$fits0, obj$H0)
  bestB_origin <- CECbest_fit_origin_info(obj$fitsB, obj$HB)

  sREO0 <- numeric_summary(obj$REO0)
  sREOB <- numeric_summary(obj$REOB)
  sKproj0 <- numeric_summary(obj$Kproj0)
  sKprojB <- numeric_summary(obj$KprojB)

  out <- data.frame(
    lambda = obj$lambda,
    partition_metric = obj$partition_metric,
    stability_approach = obj$stability_approach,

    stab0_mean = obj$stability0$mean_sim,
    stab0_sd = obj$stability0$sd_sim,
    stab0_median = obj$stability0$median_sim,
    stab0_min = obj$stability0$min_sim,
    stab0_max = obj$stability0$max_sim,
    n_valid0 = obj$stability0$n_valid,

    stabB_mean = obj$stabilityB$mean_sim,
    stabB_sd = obj$stabilityB$sd_sim,
    stabB_median = obj$stabilityB$median_sim,
    stabB_min = obj$stabilityB$min_sim,
    stabB_max = obj$stabilityB$max_sim,
    n_validB = obj$stabilityB$n_valid,

    stab_ratio_mean = obj$stab_ratio_mean,
    stab_ratio_medoid = obj$stab_ratio_medoid,

    meanH_0 = sH0$mean,
    sdH_0 = sH0$sd,
    minH_0 = sH0$min,
    q25H_0 = sH0$q25,
    medH_0 = sH0$median,
    q75H_0 = sH0$q75,
    maxH_0 = sH0$max,

    meanH_B = sHB$mean,
    sdH_B = sHB$sd,
    minH_B = sHB$min,
    q25H_B = sHB$q25,
    medH_B = sHB$median,
    q75H_B = sHB$q75,
    maxH_B = sHB$max,

    meanH0_raw = sH0raw$mean,
    sdH0_raw = sH0raw$sd,
    minH0_raw = sH0raw$min,

    meanHB_raw = sHBraw$mean,
    sdHB_raw = sHBraw$sd,
    minHB_raw = sHBraw$min,

    meanREO_0 = sREO0$mean,
    sdREO_0 = sREO0$sd,
    minREO_0 = sREO0$min,
    maxREO_0 = sREO0$max,

    meanREO_B = sREOB$mean,
    sdREO_B = sREOB$sd,
    minREO_B = sREOB$min,
    maxREO_B = sREOB$max,

    meanKproj_0 = sKproj0$mean,
    sdKproj_0 = sKproj0$sd,
    minKproj_0 = sKproj0$min,
    maxKproj_0 = sKproj0$max,

    meanKproj_B = sKprojB$mean,
    sdKproj_B = sKprojB$sd,
    minKproj_B = sKprojB$min,
    maxKproj_B = sKprojB$max,

    n0_fresh = origin0[["fresh"]],
    n0_warm_forward = origin0[["warm_forward"]],
    n0_warm_backward = origin0[["warm_backward"]],
    nB_fresh = originB[["fresh"]],
    nB_warm_forward = originB[["warm_forward"]],
    nB_warm_backward = originB[["warm_backward"]],

    best0_init_origin = best0_origin$init_origin,
    best0_path_direction = best0_origin$path_direction,
    best0_index = best0_origin$index,
    bestB_init_origin = bestB_origin$init_origin,
    bestB_path_direction = bestB_origin$path_direction,
    bestB_index = bestB_origin$index
  )

  if (obj$stability_approach == "pairwise") {
    out$stab0_medoid <- obj$stability0$medoid_mean_sim
    out$stabB_medoid <- obj$stabilityB$medoid_mean_sim
  } else {
    out$ref0_index <- obj$stability0$ref_index_original
    out$ref0_H <- obj$stability0$ref_H
    out$refB_index <- obj$stabilityB$ref_index_original
    out$refB_H <- obj$stabilityB$ref_H
  }

  out
}


# ------------------------------------------------------------
# Diagnostics over a lambda grid
# ------------------------------------------------------------
CECmovingAverage <- function(x, k = 3) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0) return(x)
  if (k <= 1) return(x)

  out <- rep(NA_real_, n)
  h <- floor(k / 2)

  for (i in seq_len(n)) {
    lo <- max(1, i - h)
    hi <- min(n, i + h)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }

  out
}


# ------------------------------------------------------------
# Enrich lambda-grid summary with smoothed columns
# ------------------------------------------------------------
#' Add smoothed summary curves to a linked-lambda diagnostic object.
#'
#' @rdname CECdiagnose_lambda_grid_linked
#' @param lambda_diag Object returned by [CECdiagnose_lambda_grid_linked()].
#' @param k Width of the moving-average smoother applied to the selected summary
#'   columns.
#' @param ratio_col,stab0_col,stabB_col,H_col Column names in
#'   `lambda_diag$summary` used for smoothing.
#' @return `CECaddSmoothedDiagnostics()` returns the input `lambda_diag` object
#'   with additional `*_smooth` columns in `summary` and a `smoothing`
#'   component describing the chosen settings.
#' @export
CECaddSmoothedDiagnostics <- function(
  lambda_diag,
  k = 3,
  ratio_col = "stab_ratio_mean",
  stab0_col = "stab0_mean",
  stabB_col = "stabB_mean",
  H_col = "meanH_B"
) {
  df <- lambda_diag$summary

  if (!(ratio_col %in% names(df))) stop("ratio_col not found in summary.")
  if (!(stab0_col %in% names(df))) stop("stab0_col not found in summary.")
  if (!(stabB_col %in% names(df))) stop("stabB_col not found in summary.")
  if (!(H_col %in% names(df))) stop("H_col not found in summary.")

  df[[paste0(ratio_col, "_smooth")]] <- CECmovingAverage(df[[ratio_col]], k = k)
  df[[paste0(stab0_col, "_smooth")]] <- CECmovingAverage(df[[stab0_col]], k = k)
  df[[paste0(stabB_col, "_smooth")]] <- CECmovingAverage(df[[stabB_col]], k = k)
  df[[paste0(H_col, "_smooth")]]     <- CECmovingAverage(df[[H_col]], k = k)

  lambda_diag$summary <- df
  lambda_diag$smoothing <- list(
    k = k,
    ratio_col = ratio_col,
    stab0_col = stab0_col,
    stabB_col = stabB_col,
    H_col = H_col
  )
  lambda_diag
}


# ------------------------------------------------------------
# Identify stable lambda regions from diagnostics
#
# rule = "ratio"       : ratio >= threshold
# rule = "stabB"       : stab_B >= threshold
# rule = "both"        : both conditions simultaneously
#
# min_consecutive = require the condition to hold on c consecutive lambdas
# for a lambda to be declared stable
# ------------------------------------------------------------
CECstableRunMask <- function(cond, min_consecutive = 1) {
  cond <- as.logical(cond)
  cond[is.na(cond)] <- FALSE

  n <- length(cond)
  min_consecutive <- max(1L, as.integer(min_consecutive[1]))

  empty_runs <- data.frame(
    start = integer(0),
    end = integer(0),
    length = integer(0)
  )

  if (n == 0) {
    return(list(mask = logical(0), runs = empty_runs))
  }

  if (min_consecutive > n) {
    return(list(mask = rep(FALSE, n), runs = empty_runs))
  }

  rr <- rle(cond)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  keep <- which(rr$values & rr$lengths >= min_consecutive)

  mask <- rep(FALSE, n)
  if (length(keep) > 0) {
    for (k in keep) {
      mask[starts[k]:ends[k]] <- TRUE
    }
  }

  runs <- empty_runs
  if (length(keep) > 0) {
    runs <- data.frame(
      start = starts[keep],
      end = ends[keep],
      length = rr$lengths[keep]
    )
  }

  list(mask = mask, runs = runs)
}

#' Identify stable lambda regions from a diagnostic summary.
#'
#' @rdname CECdiagnose_lambda_grid_linked
#' @param rule Stability rule. `"ratio"` uses the bootstrap-to-fixed stability
#'   ratio, `"stabB"` uses bootstrap stability alone, and `"both"` requires the
#'   two conditions simultaneously.
#' @param ratio_threshold,stabB_threshold Thresholds used by the selected
#'   stability rule.
#' @param use_smoothed Logical. If `TRUE`, the smoothed columns generated by
#'   [CECaddSmoothedDiagnostics()] are used when available.
#' @param min_consecutive Minimum run length required for a lambda value to be
#'   declared part of a stable region.
#' @return `CECidentifyStableLambdas()` returns a list with the first stable
#'   lambda (`lambda_min`), the complete stable mask, detected stable runs, and
#'   a compact summary data frame.
#' @export
CECidentifyStableLambdas <- function(
  lambda_diag,
  rule = c("ratio", "stabB", "both"),
  ratio_threshold = 0.8,
  stabB_threshold = 0.8,
  use_smoothed = TRUE,
  min_consecutive = 1,
  ratio_col = "stab_ratio_mean",
  stabB_col = "stabB_mean"
) {
  rule <- match.arg(rule)

  df <- lambda_diag$summary
  ratio_name <- if (use_smoothed) paste0(ratio_col, "_smooth") else ratio_col
  stabB_name <- if (use_smoothed) paste0(stabB_col, "_smooth") else stabB_col

  if (!(ratio_name %in% names(df))) stop("ratio column not found. Maybe run CECaddSmoothedDiagnostics first.")
  if (!(stabB_name %in% names(df))) stop("stabB column not found. Maybe run CECaddSmoothedDiagnostics first.")

  cond_ratio <- df[[ratio_name]] >= ratio_threshold
  cond_stabB <- df[[stabB_name]] >= stabB_threshold

  cond <- switch(
    rule,
    ratio = cond_ratio,
    stabB = cond_stabB,
    both = cond_ratio & cond_stabB
  )
  cond[is.na(cond)] <- FALSE

  stable_info <- CECstableRunMask(cond, min_consecutive = min_consecutive)
  stable_idx <- which(stable_info$mask)
  lambda_min <- if (length(stable_idx) == 0) NA_real_ else df$lambda[min(stable_idx)]

  runs <- stable_info$runs
  if (nrow(runs) > 0) {
    runs$lambda_start <- df$lambda[runs$start]
    runs$lambda_end <- df$lambda[runs$end]
  } else {
    runs$lambda_start <- numeric(0)
    runs$lambda_end <- numeric(0)
  }

  list(
    lambda_min = lambda_min,
    stable_lambda = df$lambda[stable_idx],
    stable_index = stable_idx,
    stable = stable_info$mask,
    unstable = !stable_info$mask,
    raw_condition = cond,
    runs = runs,
    rule = rule,
    use_smoothed = use_smoothed,
    min_consecutive = min_consecutive,
    ratio_threshold = ratio_threshold,
    stabB_threshold = stabB_threshold,
    summary = data.frame(
      lambda = df$lambda,
      ratio = df[[ratio_name]],
      stabB = df[[stabB_name]],
      condition = cond,
      stable = stable_info$mask
    )
  )
}

# ------------------------------------------------------------
# Backward-compatible wrapper:
# returns the first stable lambda, but the condition now reflects
# the lambdas declared stable on the whole grid.
# ------------------------------------------------------------
CECidentifyLambdaMin <- function(
  lambda_diag,
  rule = c("ratio", "stabB", "both"),
  ratio_threshold = 0.8,
  stabB_threshold = 0.8,
  use_smoothed = TRUE,
  min_consecutive = 1,
  ratio_col = "stab_ratio_mean",
  stabB_col = "stabB_mean"
) {
  out <- CECidentifyStableLambdas(
    lambda_diag = lambda_diag,
    rule = rule,
    ratio_threshold = ratio_threshold,
    stabB_threshold = stabB_threshold,
    use_smoothed = use_smoothed,
    min_consecutive = min_consecutive,
    ratio_col = ratio_col,
    stabB_col = stabB_col
  )

  out$condition <- out$stable
  out$summary$selected <- out$summary$stable
  out
}


# ------------------------------------------------------------
# Extract best partition for one lambda
#
# source = "fixed"     -> among repeated runs on Z
# source = "bootstrap" -> among bootstrap fits projected on Z
# source = "all"       -> among both
#
# criterion = "projected_H" (recommended)
# criterion = "raw_H"
# ------------------------------------------------------------
CECbestPartitionOneLambda <- function(
  lambda_detail,
  source = c("fixed", "bootstrap", "all"),
  criterion = c("projected_H", "raw_H")
) {
  source <- match.arg(source)
  criterion <- match.arg(criterion)

  candidates <- list()
  rows <- list()

  add_block <- function(parts, fits, Hproj, Hraw, tag) {
    m <- length(parts)
    if (m == 0) return(NULL)

    for (i in seq_len(m)) {
      if (!is.null(parts[[i]]) && !is.na(Hproj[i])) {
        candidates[[length(candidates) + 1L]] <<- list(
          partition = parts[[i]],
          fit = fits[[i]],
          H_projected = Hproj[i],
          H_raw = Hraw[i],
          source = tag,
          source_index = i,
          init_origin = CECget_fit_origin_label(fits[[i]]),
          path_direction = CECget_fit_path_direction(fits[[i]])
        )
        rows[[length(rows) + 1L]] <<- data.frame(
          source = tag,
          source_index = i,
          init_origin = CECget_fit_origin_label(fits[[i]]),
          path_direction = CECget_fit_path_direction(fits[[i]]),
          H_projected = Hproj[i],
          H_raw = Hraw[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (source %in% c("fixed", "all")) {
    add_block(
      parts = lambda_detail$partitions0,
      fits = lambda_detail$fits0,
      Hproj = lambda_detail$H0,
      Hraw = lambda_detail$H0_raw,
      tag = "fixed"
    )
  }

  if (source %in% c("bootstrap", "all")) {
    add_block(
      parts = lambda_detail$partitionsB,
      fits = lambda_detail$fitsB,
      Hproj = lambda_detail$HB,
      Hraw = lambda_detail$HB_raw,
      tag = "bootstrap"
    )
  }

  if (length(candidates) == 0) return(NULL)

  tab <- do.call(rbind, rows)
  score_values <- if (criterion == "projected_H") tab$H_projected else tab$H_raw
  idx <- .best_finite_index(score_values)
  if (length(idx) == 0L) {
    fallback_values <- if (criterion == "projected_H") tab$H_raw else tab$H_projected
    idx <- .best_finite_index(fallback_values)
  }
  if (length(idx) == 0L) return(NULL)
  best <- candidates[[idx]]

  REO_best <- if (!is.null(best$fit)) best$fit$REO else length(unique(best$partition))

  list(
    lambda = lambda_detail$lambda,
    source = best$source,
    source_index = best$source_index,
    init_origin = best$init_origin,
    path_direction = best$path_direction,
    criterion = criterion,
    H_projected = best$H_projected,
    H_raw = best$H_raw,
    partition = best$partition,
    fit = best$fit,
    REO = REO_best,
    candidates_table = tab[order(tab$H_projected), , drop = FALSE]
  )
}


# ------------------------------------------------------------
# Extract best partition for each lambda of a grid diagnosis
# ------------------------------------------------------------
CECrelabelToReference <- function(ref, target) {
  ref <- as.integer(as.factor(ref))
  target <- as.integer(as.factor(target))

  if (length(ref) != length(target)) stop("ref and target must have same length.")

  r <- length(unique(ref))
  s <- length(unique(target))

  tab <- table(target, ref)
  tab <- as.matrix(tab)

  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' required for relabeling. Install it.")
  }

  # rows = target labels, cols = ref labels
  nr <- nrow(tab)
  nc <- ncol(tab)
  m <- max(nr, nc)

  M <- matrix(0, nrow = m, ncol = m)
  M[1:nr, 1:nc] <- tab

  assign <- clue::solve_LSAP(M, maximum = TRUE)

  map <- rep(NA_integer_, nr)
  for (i in seq_len(nr)) {
    if (assign[i] <= nc) map[i] <- assign[i]
  }

  out <- target
  uniq_t <- sort(unique(target))

  used_ref <- unique(na.omit(map))
  next_new <- if (length(used_ref) == 0) 1L else max(used_ref) + 1L

  for (j in seq_along(uniq_t)) {
    lab <- uniq_t[j]
    pos <- which(uniq_t == lab)
    if (!is.na(map[pos])) {
      out[target == lab] <- map[pos]
    } else {
      out[target == lab] <- next_new
      next_new <- next_new + 1L
    }
  }

  out
}


# ------------------------------------------------------------
# Detect partition changes along lambda
#
# A change is flagged if:
# - REO changes, OR
# - similarity drops below sim_threshold.
# By default the comparison is local between consecutive lambdas. With
# method = "regime_reference", each lambda is compared to the first partition
# of the current regime; a new regime starts whenever a change is detected.
# ------------------------------------------------------------
#' Detect change points along a best-partition trajectory.
#'
#' @rdname CECextractBestPartitions
#' @param best_partitions_obj Object returned by [CECextractBestPartitions()].
#' @param partition_metric Similarity measure used to compare consecutive
#'   partitions.
#' @param sim_threshold Threshold below which a similarity-based change is
#'   detected.
#' @param criterion Rule used to declare a change. In
#'   `CECextractBestPartitions()` this argument selects the criterion used to
#'   choose the representative partition at each lambda. In
#'   `CECdetectPartitionChanges()` it selects the rule combining REO changes and
#'   partition similarity.
#' @param include_first Logical. If `TRUE`, the first lambda value is always
#'   reported as the start of a new regime.
#' @param method Change-detection method. `"consecutive"` compares each
#'   partition to the previous lambda. `"regime_reference"` compares each
#'   partition to the reference partition at the beginning of the current
#'   regime and updates the reference when a change is detected.
#' @param lambda_subset Optional subset of lambda values on which changes are
#'   detected. This is useful for detecting regimes only on stable lambdas.
#' @return `CECdetectPartitionChanges()` returns a list with a summary data
#'   frame, the lambda values flagged as change points, and the settings used to
#'   detect them.
#' @export
CECdetectPartitionChanges <- function(
  best_partitions_obj,
  partition_metric = c("match", "ARI"),
  sim_threshold = 0.80,
  criterion = c("REO_or_threshold", "REO_only", "threshold_only", "REO_and_threshold"),
  include_first = TRUE,
  method = c("consecutive", "regime_reference"),
  lambda_subset = NULL
) {
  partition_metric <- match.arg(partition_metric)
  criterion <- match.arg(criterion)
  method <- match.arg(method)

  best_list <- best_partitions_obj$best
  ok <- vapply(best_list, function(x) !is.null(x), logical(1))
  best_list <- best_list[ok]

  if (length(best_list) == 0) {
    return(list(
      summary = data.frame(),
      change_lambdas = numeric(0),
      change_indices = integer(0),
      criterion = criterion,
      partition_metric = partition_metric,
      sim_threshold = sim_threshold,
      include_first = include_first,
      method = method,
      lambda_subset = lambda_subset
    ))
  }

  lambda <- vapply(best_list, function(x) x$lambda, numeric(1))
  if (!is.null(lambda_subset)) {
    keep <- lambda %in% lambda_subset
    best_list <- best_list[keep]
    lambda <- lambda[keep]
  }

  if (length(best_list) == 0) {
    return(list(
      summary = data.frame(),
      change_lambdas = numeric(0),
      change_indices = integer(0),
      criterion = criterion,
      partition_metric = partition_metric,
      sim_threshold = sim_threshold,
      include_first = include_first,
      method = method,
      lambda_subset = lambda_subset
    ))
  }

  REO <- vapply(best_list, function(x) x$REO, numeric(1))
  parts <- lapply(best_list, function(x) x$partition)
  part_lengths <- vapply(parts, length, integer(1))
  expected_length <- as.integer(names(sort(table(part_lengths), decreasing = TRUE))[1])
  ok_length <- part_lengths == expected_length

  if (any(!ok_length)) {
    warning(
      "Dropping ", sum(!ok_length),
      " partition(s) with length inconsistent with the main trajectory in ",
      "CECdetectPartitionChanges().",
      call. = FALSE
    )
    best_list <- best_list[ok_length]
    lambda <- lambda[ok_length]
    REO <- REO[ok_length]
    parts <- parts[ok_length]
  }

  n <- length(best_list)
  if (n == 0) {
    return(list(
      summary = data.frame(),
      change_lambdas = numeric(0),
      change_indices = integer(0),
      criterion = criterion,
      partition_metric = partition_metric,
      sim_threshold = sim_threshold,
      include_first = include_first,
      method = method,
      lambda_subset = lambda_subset
    ))
  }

  sim_prev <- rep(NA_real_, n)
  sim_ref <- rep(NA_real_, n)
  ref_index <- rep(NA_integer_, n)
  REO_change <- rep(FALSE, n)
  threshold_change <- rep(FALSE, n)
  is_change <- rep(FALSE, n)
  current_ref <- 1L
  ref_index[1] <- current_ref

  for (i in 2:n) {
    sim_prev[i] <- partition_similarity(
      parts[[i - 1]],
      parts[[i]],
      method = partition_metric
    )

    if (method == "consecutive") {
      ref_index[i] <- i - 1L
      sim_ref[i] <- sim_prev[i]
      REO_change[i] <- (REO[i] != REO[i - 1])
      threshold_change[i] <- (!is.na(sim_prev[i]) && sim_prev[i] < sim_threshold)
    } else {
      ref_index[i] <- current_ref
      sim_ref[i] <- partition_similarity(
        parts[[current_ref]],
        parts[[i]],
        method = partition_metric
      )
      REO_change[i] <- (REO[i] != REO[current_ref])
      threshold_change[i] <- (!is.na(sim_ref[i]) && sim_ref[i] < sim_threshold)
    }

    is_change[i] <- switch(
      criterion,
      REO_only = REO_change[i],
      threshold_only = threshold_change[i],
      REO_or_threshold = REO_change[i] || threshold_change[i],
      REO_and_threshold = REO_change[i] && threshold_change[i]
    )

    if (method == "regime_reference" && is_change[i]) {
      current_ref <- i
    }
  }

  if (include_first) {
    is_change[1] <- TRUE
  }

  summary_df <- data.frame(
    index = seq_len(n),
    lambda = lambda,
    REO = REO,
    sim_prev = sim_prev,
    ref_index = ref_index,
    ref_lambda = lambda[ref_index],
    sim_ref = sim_ref,
    REO_change = REO_change,
    threshold_change = threshold_change,
    is_change = is_change
  )

  list(
    summary = summary_df,
    change_lambdas = lambda[is_change],
    change_indices = which(is_change),
    criterion = criterion,
    partition_metric = partition_metric,
    sim_threshold = sim_threshold,
    include_first = include_first,
    method = method,
    lambda_subset = lambda_subset
  )
}

# ------------------------------------------------------------
# Helpers for shading stable / unstable lambda regions on plots
# ------------------------------------------------------------
CECcomputeLambdaIntervals <- function(lambda) {
  lambda <- as.numeric(lambda)
  n <- length(lambda)

  if (n == 0) {
    return(data.frame(lambda = numeric(0), lower = numeric(0), upper = numeric(0)))
  }

  if (n == 1) {
    delta <- max(0.5, abs(lambda[1]) * 0.05)
    return(data.frame(
      lambda = lambda,
      lower = lambda - delta,
      upper = lambda + delta
    ))
  }

  mids <- (lambda[-1] + lambda[-n]) / 2
  first_step <- mids[1] - lambda[1]
  last_step <- lambda[n] - mids[length(mids)]

  data.frame(
    lambda = lambda,
    lower = c(lambda[1] - first_step, mids),
    upper = c(mids, lambda[n] + last_step)
  )
}

CECextractStableMask <- function(lambda, lambda_stable_obj) {
  if (is.null(lambda_stable_obj)) return(NULL)

  if (!is.null(lambda_stable_obj$summary) && "lambda" %in% names(lambda_stable_obj$summary)) {
    summary_df <- lambda_stable_obj$summary
    stable_name <- NULL

    if ("stable" %in% names(summary_df)) {
      stable_name <- "stable"
    } else if ("selected" %in% names(summary_df)) {
      stable_name <- "selected"
    }

    if (!is.null(stable_name)) {
      idx <- match(lambda, summary_df$lambda)
      out <- rep(NA, length(lambda))
      ok <- !is.na(idx)
      out[ok] <- as.logical(summary_df[[stable_name]][idx[ok]])
      return(out)
    }
  }

  if (!is.null(lambda_stable_obj$stable) && length(lambda_stable_obj$stable) == length(lambda)) {
    return(as.logical(lambda_stable_obj$stable))
  }

  if (!is.null(lambda_stable_obj$condition) && length(lambda_stable_obj$condition) == length(lambda)) {
    return(as.logical(lambda_stable_obj$condition))
  }

  NULL
}

CECshadeLambdaBands <- function(
  lambda,
  stable_mask,
  axis = c("x", "y"),
  col = grDevices::adjustcolor("grey70", alpha.f = 0.35)
) {
  axis <- match.arg(axis)

  if (is.null(stable_mask) || length(stable_mask) != length(lambda)) {
    return(invisible(NULL))
  }

  stable_mask <- as.logical(stable_mask)
  stable_mask[is.na(stable_mask)] <- FALSE
  unstable_idx <- which(!stable_mask)

  if (length(unstable_idx) == 0) {
    return(invisible(NULL))
  }

  bands <- CECcomputeLambdaIntervals(lambda)
  usr <- par("usr")

  for (i in unstable_idx) {
    if (axis == "x") {
      rect(bands$lower[i], usr[3], bands$upper[i], usr[4], col = col, border = NA)
    } else {
      rect(usr[1], bands$lower[i], usr[2], bands$upper[i], col = col, border = NA)
    }
  }

  invisible(bands[unstable_idx, , drop = FALSE])
}

# ------------------------------------------------------------
# Plot smoothed diagnostics + stable / unstable lambda regions
# ------------------------------------------------------------
#' Plot diagnostic summaries and interactive views for CEClust workflows.
#'
#' The plotting helpers provide complementary views of a CEClust analysis:
#' static summaries of the lambda grid, one-dimensional inspection of extracted
#' partitions, PCA-based views of partition trajectories, and an interactive
#' Shiny explorer for mixed-data results.
#'
#' @rdname plotCECdiagnoseInteractive
#' @param Z Original data set used for plotting. For
#'   `plotCECBestPartitions1D()`, this should be one-dimensional and numeric.
#' @param lambda_diag Object returned by [CECdiagnose_lambda_grid_linked()].
#' @param ratio_col,stab0_col,stabB_col,H_col Column names used in the lambda
#'   summary plot.
#' @param use_smoothed Logical. If `TRUE`, smoothed summary columns are used when
#'   available.
#' @param lambda_min_obj Deprecated compatibility object containing
#'   `lambda_min`.
#' @param lambda_stable_obj Optional output of [CECidentifyStableLambdas()] used
#'   to shade stable and unstable regions.
#' @param unstable_col Colour used to shade unstable lambda bands.
#' @return `plotCECdiagnoseLambdaGrid()` is called for its side effects and
#'   returns `NULL` invisibly.
#' @export
plotCECdiagnoseLambdaGrid <- function(
  lambda_diag,
  ratio_col = "stab_ratio_mean",
  stab0_col = "stab0_mean",
  stabB_col = "stabB_mean",
  H_col = "meanH_B",
  use_smoothed = TRUE,
  lambda_min_obj = NULL,
  lambda_stable_obj = NULL,
  unstable_col = grDevices::adjustcolor("grey70", alpha.f = 0.35)
) {
  df <- lambda_diag$summary

  ratio_name <- if (use_smoothed) paste0(ratio_col, "_smooth") else ratio_col
  stab0_name <- if (use_smoothed) paste0(stab0_col, "_smooth") else stab0_col
  stabB_name <- if (use_smoothed) paste0(stabB_col, "_smooth") else stabB_col
  H_name     <- if (use_smoothed) paste0(H_col, "_smooth") else H_col

  if (is.null(lambda_stable_obj)) {
    lambda_stable_obj <- lambda_min_obj
  }

  stable_mask <- CECextractStableMask(df$lambda, lambda_stable_obj)

  lambda_ref <- NA_real_
  if (!is.null(lambda_stable_obj) && !is.null(lambda_stable_obj$lambda_min)) {
    lambda_ref <- lambda_stable_obj$lambda_min
  } else if (!is.null(lambda_min_obj) && !is.null(lambda_min_obj$lambda_min)) {
    lambda_ref <- lambda_min_obj$lambda_min
  }

  ratio_threshold <- if (!is.null(lambda_stable_obj) && !is.null(lambda_stable_obj$ratio_threshold)) {
    lambda_stable_obj$ratio_threshold
  } else {
    0.8
  }

  stabB_threshold <- if (!is.null(lambda_stable_obj) && !is.null(lambda_stable_obj$stabB_threshold)) {
    lambda_stable_obj$stabB_threshold
  } else {
    0.8
  }

  plot_metric <- function(y, ylab, main, ylim = NULL, h = NULL) {
    yy <- as.numeric(y)
    finite_y <- yy[is.finite(yy)]

    if (is.null(ylim)) {
      ref_y <- finite_y
      if (!is.null(h)) ref_y <- c(ref_y, h)
      if (length(ref_y) == 0) ref_y <- c(0, 1)
      ylim <- range(ref_y)
      if (diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
    }

    plot(
      df$lambda, yy,
      type = "n",
      ylim = ylim,
      xlab = expression(lambda),
      ylab = ylab,
      main = main
    )

    CECshadeLambdaBands(df$lambda, stable_mask, axis = "x", col = unstable_col)
    lines(df$lambda, yy, type = "b", pch = 19)

    if (!is.null(h)) {
      for (hh in h) abline(h = hh, lty = 2)
    }

    if (!is.na(lambda_ref)) {
      abline(v = lambda_ref, lty = 2)
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 1), mgp = c(2.2, 0.7, 0))

  plot_metric(df[[stab0_name]], ylab = "stability", main = expression(stab[0](lambda)), ylim = c(0, 1))
  plot_metric(df[[stabB_name]], ylab = "stability", main = expression(stab[B](lambda)), ylim = c(0, 1), h = stabB_threshold)
  plot_metric(df[[ratio_name]], ylab = "ratio", main = expression(stab[B](lambda) / stab[0](lambda)), h = ratio_threshold)
  plot_metric(df[[H_name]], ylab = "mean projected H", main = "meanH_B(lambda)")

  invisible(df)
}

# ------------------------------------------------------------
# 1D visualization of best partition evolution with lambda
#
# x-axis: sorted Z
# y-axis: lambda
# color: cluster label after relabeling for continuity
# ------------------------------------------------------------
plotCECPartitionEvolution1D <- function(
  Z,
  best_partitions_obj,
  lambda_subset = NULL,
  reorder_Z = TRUE,
  relabel_continuously = TRUE,
  main = "Partition evolution with lambda",
  col = NULL,
  use_rainbow = FALSE,
  piecewise_constant = TRUE,
  lambda_stable_obj = NULL,
  unstable_col = grDevices::adjustcolor("grey70", alpha.f = 0.35)
) {
  if (!is.numeric(Z)) stop("This plot is only for 1D numeric Z.")

  best_list <- best_partitions_obj$best
  ok <- vapply(best_list, function(x) !is.null(x), logical(1))
  best_list <- best_list[ok]

  if (length(best_list) == 0) stop("No valid best partitions.")

  lambda_all <- vapply(best_list, function(x) x$lambda, numeric(1))
  part_all   <- lapply(best_list, function(x) x$partition)
  stable_mask <- CECextractStableMask(lambda_all, lambda_stable_obj)

  ord <- if (reorder_Z) order(Z) else seq_along(Z)
  z_ord <- Z[ord]

  # Continuous relabeling over the full lambda grid.
  aligned_all <- vector("list", length(part_all))
  for (i in seq_along(part_all)) {
    p <- part_all[[i]][ord]
    if (i > 1 && relabel_continuously) {
      p <- CECrelabelToReference(aligned_all[[i - 1]], p)
    }
    aligned_all[[i]] <- p
  }

  Kmax <- max(vapply(aligned_all, function(p) length(unique(p)), integer(1)))
  n <- length(Z)
  m_all <- length(lambda_all)

  if (is.null(col)) {
    if (use_rainbow) {
      col <- rainbow(Kmax)
    } else {
      base_cols <- c(
        "#D73027", "#4575B4", "#1A9850", "#984EA3", "#FF7F00",
        "#A65628", "#F781BF", "#4D4D4D", "#66C2A5", "#E6AB02"
      )
      if (Kmax <= length(base_cols)) {
        col <- base_cols[1:Kmax]
      } else {
        col <- grDevices::colorRampPalette(base_cols)(Kmax)
      }
    }
  }

  # ----------------------------------------------------------
  # Cas 1 : pas de sous-ensemble -> comportement standard
  # ----------------------------------------------------------
  if (is.null(lambda_subset)) {
    M <- matrix(NA_integer_, nrow = m_all, ncol = n)
    for (i in seq_len(m_all)) {
      M[i, ] <- aligned_all[[i]]
    }

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    image(
      x = z_ord,
      y = lambda_all,
      z = t(M),
      col = col,
      xlab = "sorted Z",
      ylab = expression(lambda),
      main = main,
      axes = FALSE
    )
    CECshadeLambdaBands(lambda_all, stable_mask, axis = "y", col = unstable_col)
    axis(1)
    axis(2, at = lambda_all, labels = lambda_all)
    box()

    return(invisible(list(
      matrix = M,
      lambda = lambda_all,
      z_sorted = z_ord,
      partitions_sorted = aligned_all,
      stable = stable_mask,
      col = col
    )))
  }

  # ----------------------------------------------------------
  # Cas 2 : sous-ensemble de lambdas
  # -> propagation par morceaux constants
  # ----------------------------------------------------------
  keep_idx <- which(lambda_all %in% lambda_subset)
  if (length(keep_idx) == 0) {
    stop("No lambda from lambda_subset was found in best_partitions_obj.")
  }

  keep_idx <- sort(unique(keep_idx))

  # Mfull[i, ] is the partition displayed on row i.
  Mfull <- matrix(NA_integer_, nrow = m_all, ncol = n)

  if (piecewise_constant) {
    for (h in seq_along(keep_idx)) {
      i_start <- if (h == 1L) 1L else keep_idx[h]
      i_end <- if (h < length(keep_idx)) keep_idx[h + 1] - 1L else m_all
      for (i in i_start:i_end) {
        Mfull[i, ] <- aligned_all[[keep_idx[h]]]
      }
    }
  } else {
    for (i in keep_idx) {
      Mfull[i, ] <- aligned_all[[i]]
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  image(
    x = z_ord,
    y = lambda_all,
    z = t(Mfull),
    col = col,
    xlab = "sorted Z",
    ylab = expression(lambda),
    main = main,
    axes = FALSE
  )
  CECshadeLambdaBands(lambda_all, stable_mask, axis = "y", col = unstable_col)
  axis(1)
  axis(2, at = lambda_all, labels = lambda_all)
  box()

  invisible(list(
    matrix = Mfull,
    lambda = lambda_all,
    lambda_subset = lambda_all[keep_idx],
    z_sorted = z_ord,
    partitions_sorted = aligned_all,
    displayed_partitions = aligned_all[keep_idx],
    stable = stable_mask,
    col = col
  ))
}

# ------------------------------------------------------------
# Composite entropy criterion decomposition.
# H = H(nu) + sum_x nu_x H(G_x | g_{theta_x})
# ------------------------------------------------------------
evalCompositeEntropyDecomposed <- function(phi, Z, lambda, C, familyType = "gaussAndDiscreteVector") {
  params <- optParam(Z = Z, phi = phi, lambda = lambda, C = C, familyType = familyType)
  
  nu <- params$nu
  toKeep <- which(nu > 0)
  nu <- nu[toKeep]
  
  # terme d'entropie de classes
  H_class <- -sum(nu * log(nu))
  
  # terme conditionnel
  logg <- CECdens(Z = Z, param = params, lambda = lambda, applyLog = TRUE)
  
  n <- if (is.matrix(Z) || is.data.frame(Z)) nrow(Z) else length(Z)
  phi_use <- if (!is.null(params$phi)) params$phi else phi
  r <- length(params$states)
  if (length(phi_use) != n * r) {
    stop("Inconsistent phi length in evalCompositeEntropyDecomposed.")
  }
  phiM <- vectToMat(phi_use, r)
  
  H_cond <- 0
  for (i in toKeep) {
    idx_i <- ((i - 1) * n + 1):(i * n)
    phi_i <- phiM[, i]
    active_i <- which(phi_i > 0)

    if (length(active_i) > 0L) {
      phi_active <- phi_i[active_i]
      dens_active <- logg$density[idx_i][active_i]
      sphi_i <- sum(phi_active)
      H_cond <- H_cond - (1 / lambda) * nu[which(toKeep == i)] *
        sum(phi_active * dens_active) / sphi_i
    }
  }
  
  list(
    H_total = H_class + H_cond,
    H_class = H_class,
    H_cond = H_cond,
    nu = params$nu,
    params = params
  )
}


# ------------------------------------------------------------
# Hard partition helpers
# ------------------------------------------------------------
CECpartitionToPhi <- function(partition) {
  partition <- as.integer(as.factor(partition))
  n <- length(partition)
  r <- length(unique(partition))

  phiM <- matrix(0, nrow = n, ncol = r)
  phiM[cbind(seq_len(n), partition)] <- 1
  matToVect(phiM)
}

CECevaluatePartitionOnLambda <- function(
  partition,
  Z,
  lambda,
  C = 1,
  familyType = "gaussAndDiscreteVector"
) {
  phi <- CECpartitionToPhi(partition)

  dec <- evalCompositeEntropyDecomposed(
    phi = phi,
    Z = Z,
    lambda = lambda,
    C = C,
    familyType = familyType
  )

  list(
    H_total = dec$H_total,
    H_class = dec$H_class,
    H_cond = dec$H_cond,
    phi = phi,
    params = dec$params
  )
}

# ------------------------------------------------------------
# Check coherence of best partitions along the lambda grid
#
# For each k:
# H_{lambda_{k-1}}(P_k) >= H_{lambda_{k-1}}(P_{k-1})
# H_{lambda_{k+1}}(P_k) >= H_{lambda_{k+1}}(P_{k+1})
# ------------------------------------------------------------
#' Check whether extracted best partitions remain locally coherent.
#'
#' @rdname CECextractBestPartitions
#' @param Z Original data set used to evaluate neighbouring partitions.
#' @param tol Numerical tolerance used when comparing criterion values.
#' @return `CECcheckBestPartitionCoherence()` returns a list with a row per
#'   lambda, a global coherence flag, and the indices/lambda values of
#'   incoherent positions.
#' @export
CECcheckBestPartitionCoherence <- function(
  best_partitions_obj,
  Z,
  tol = 1e-10
) {
  best_list <- best_partitions_obj$best
  ok <- vapply(best_list, function(x) !is.null(x), logical(1))
  best_list <- best_list[ok]

  if (length(best_list) <= 1) {
    summary_df <- data.frame(
      index = seq_along(best_list),
      lambda = vapply(best_list, function(x) x$lambda, numeric(1)),
      prev_lambda = NA_real_,
      H_prev_of_Pk = NA_real_,
      H_prev_of_Pkm1 = NA_real_,
      prev_gap = NA_real_,
      prev_ok = NA,
      next_lambda = NA_real_,
      H_next_of_Pk = NA_real_,
      H_next_of_Pkp1 = NA_real_,
      next_gap = NA_real_,
      next_ok = NA,
      coherent = NA
    )

    return(list(
      summary = summary_df,
      all_coherent = TRUE,
      n_incoherent = 0L,
      incoherent_indices = integer(0),
      incoherent_lambdas = numeric(0),
      tol = tol
    ))
  }

  ord <- order(vapply(best_list, function(x) x$lambda, numeric(1)))
  best_list <- best_list[ord]

  n <- length(best_list)
  lambda <- vapply(best_list, function(x) x$lambda, numeric(1))
  partitions <- lapply(best_list, function(x) x$partition)

  get_fit_param <- function(obj, name, default = NULL) {
    if (!is.null(obj$fit) && !is.null(obj$fit$params) && !is.null(obj$fit$params[[name]])) {
      return(obj$fit$params[[name]])
    }
    default
  }

  familyType <- vapply(
    best_list,
    function(x) {
      val <- get_fit_param(x, "familyType", default = NA_character_)
      if (length(val) == 0) NA_character_ else as.character(val[1])
    },
    character(1)
  )

  Cvec <- vapply(
    best_list,
    function(x) {
      val <- get_fit_param(x, "C", default = NA_real_)
      if (length(val) == 0) NA_real_ else as.numeric(val[1])
    },
    numeric(1)
  )

  familyType_default <- familyType[which(!is.na(familyType))[1]]
  C_default <- Cvec[which(!is.na(Cvec))[1]]

  if (length(familyType_default) == 0 || is.na(familyType_default)) {
    stop("Unable to infer familyType from best_partitions_obj.")
  }
  if (length(C_default) == 0 || is.na(C_default)) {
    stop("Unable to infer C from best_partitions_obj.")
  }

  prev_lambda <- rep(NA_real_, n)
  H_prev_of_Pk <- rep(NA_real_, n)
  H_prev_of_Pkm1 <- rep(NA_real_, n)
  prev_gap <- rep(NA_real_, n)
  prev_ok <- rep(NA, n)

  next_lambda <- rep(NA_real_, n)
  H_next_of_Pk <- rep(NA_real_, n)
  H_next_of_Pkp1 <- rep(NA_real_, n)
  next_gap <- rep(NA_real_, n)
  next_ok <- rep(NA, n)

  for (k in seq_len(n)) {
    if (k > 1) {
      lam_prev <- lambda[k - 1]
      C_prev <- if (is.na(Cvec[k - 1])) C_default else Cvec[k - 1]
      family_prev <- if (is.na(familyType[k - 1])) familyType_default else familyType[k - 1]

      prev_lambda[k] <- lam_prev
      H_prev_of_Pk[k] <- CECevaluatePartitionOnLambda(
        partition = partitions[[k]],
        Z = Z,
        lambda = lam_prev,
        C = C_prev,
        familyType = family_prev
      )$H_total
      H_prev_of_Pkm1[k] <- CECevaluatePartitionOnLambda(
        partition = partitions[[k - 1]],
        Z = Z,
        lambda = lam_prev,
        C = C_prev,
        familyType = family_prev
      )$H_total
      prev_gap[k] <- H_prev_of_Pk[k] - H_prev_of_Pkm1[k]
      prev_ok[k] <- isTRUE(prev_gap[k] >= -tol)
    }

    if (k < n) {
      lam_next <- lambda[k + 1]
      C_next <- if (is.na(Cvec[k + 1])) C_default else Cvec[k + 1]
      family_next <- if (is.na(familyType[k + 1])) familyType_default else familyType[k + 1]

      next_lambda[k] <- lam_next
      H_next_of_Pk[k] <- CECevaluatePartitionOnLambda(
        partition = partitions[[k]],
        Z = Z,
        lambda = lam_next,
        C = C_next,
        familyType = family_next
      )$H_total
      H_next_of_Pkp1[k] <- CECevaluatePartitionOnLambda(
        partition = partitions[[k + 1]],
        Z = Z,
        lambda = lam_next,
        C = C_next,
        familyType = family_next
      )$H_total
      next_gap[k] <- H_next_of_Pk[k] - H_next_of_Pkp1[k]
      next_ok[k] <- isTRUE(next_gap[k] >= -tol)
    }
  }

  coherent <- rep(TRUE, n)
  for (k in seq_len(n)) {
    checks_k <- c(prev_ok[k], next_ok[k])
    checks_k <- checks_k[!is.na(checks_k)]
    coherent[k] <- if (length(checks_k) == 0) TRUE else all(checks_k)
  }

  summary_df <- data.frame(
    index = seq_len(n),
    lambda = lambda,
    prev_lambda = prev_lambda,
    H_prev_of_Pk = H_prev_of_Pk,
    H_prev_of_Pkm1 = H_prev_of_Pkm1,
    prev_gap = prev_gap,
    prev_ok = prev_ok,
    next_lambda = next_lambda,
    H_next_of_Pk = H_next_of_Pk,
    H_next_of_Pkp1 = H_next_of_Pkp1,
    next_gap = next_gap,
    next_ok = next_ok,
    coherent = coherent
  )

  incoherent_idx <- which(!summary_df$coherent)

  list(
    summary = summary_df,
    all_coherent = length(incoherent_idx) == 0,
    n_incoherent = length(incoherent_idx),
    incoherent_indices = incoherent_idx,
    incoherent_lambdas = summary_df$lambda[incoherent_idx],
    tol = tol
  )
}

CECbestPartitionSummaryRow <- function(obj) {
  if (is.null(obj)) {
    return(data.frame(
      lambda = NA_real_,
      source = NA_character_,
      source_index = NA_integer_,
      init_origin = NA_character_,
      path_direction = NA_character_,
      criterion = NA_character_,
      H_projected = NA_real_,
      H_class_projected = NA_real_,
      H_cond_projected = NA_real_,
      H_raw = NA_real_,
      H_class_raw = NA_real_,
      H_cond_raw = NA_real_,
      REO = NA_integer_,
      repaired = NA,
      propagated_from_lambda = NA_real_,
      repair_origin_lambda = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    lambda = .scalar_or_na(obj$lambda),
    source = as.character(.scalar_or_na(obj$source, default = NA_character_)),
    source_index = .scalar_or_na(obj$source_index, default = NA_integer_),
    init_origin = as.character(.scalar_or_na(obj$init_origin, default = NA_character_)),
    path_direction = as.character(.scalar_or_na(obj$path_direction, default = NA_character_)),
    criterion = as.character(.scalar_or_na(obj$criterion, default = NA_character_)),
    H_projected = .scalar_or_na(obj$H_projected),
    H_class_projected = .scalar_or_na(obj$H_class_projected),
    H_cond_projected = .scalar_or_na(obj$H_cond_projected),
    H_raw = .scalar_or_na(obj$H_raw),
    H_class_raw = .scalar_or_na(obj$H_class_raw),
    H_cond_raw = .scalar_or_na(obj$H_cond_raw),
    REO = .scalar_or_na(obj$REO, default = NA_integer_),
    repaired = if (is.null(obj$repaired)) FALSE else isTRUE(obj$repaired),
    propagated_from_lambda = .scalar_or_na(obj$propagated_from_lambda),
    repair_origin_lambda = .scalar_or_na(obj$repair_origin_lambda),
    stringsAsFactors = FALSE
  )
}

CECbuildBestPartitionObjectAtLambda <- function(
  partition,
  lambda,
  Z,
  C = 1,
  familyType = "gaussAndDiscreteVector",
  template_obj = NULL,
  eval_obj = NULL,
  criterion = NULL
) {
  partition <- as.integer(as.factor(partition))

  if (is.null(eval_obj)) {
    eval_obj <- CECevaluatePartitionOnLambda(
      partition = partition,
      Z = Z,
      lambda = lambda,
      C = C,
      familyType = familyType
    )
  }

  r <- length(unique(partition))
  phiM <- vectToMat(eval_obj$phi, r)

  fit_obj <- list(
    phi = eval_obj$phi,
    params = eval_obj$params,
    Hphi = eval_obj$H_total,
    REO = r,
    clusters = partition
  )

  projection_obj <- list(
    phi = eval_obj$phi,
    phiM = phiM,
    clusters = partition,
    nu = eval_obj$params$nu,
    K = sum(eval_obj$params$nu > 0),
    H = eval_obj$H_total,
    H_class = eval_obj$H_class,
    H_cond = eval_obj$H_cond
  )

  source_val <- if (!is.null(template_obj) && !is.null(template_obj$source)) {
    template_obj$source
  } else {
    NA_character_
  }

  source_index_val <- if (!is.null(template_obj) && !is.null(template_obj$source_index)) {
    template_obj$source_index
  } else {
    NA_integer_
  }

  init_origin_val <- if (!is.null(template_obj) && !is.null(template_obj$init_origin)) {
    template_obj$init_origin
  } else if (!is.null(template_obj) && !is.null(template_obj$fit)) {
    CECget_fit_origin_label(template_obj$fit)
  } else {
    NA_character_
  }

  path_direction_val <- if (!is.null(template_obj) && !is.null(template_obj$path_direction)) {
    template_obj$path_direction
  } else if (!is.null(template_obj) && !is.null(template_obj$fit)) {
    CECget_fit_path_direction(template_obj$fit)
  } else {
    NA_character_
  }

  fit_obj <- CECtag_fit_origin(
    fit_obj,
    init_origin = if (!is.na(init_origin_val) && startsWith(init_origin_val, "warm")) "warm" else init_origin_val,
    path_direction = path_direction_val
  )

  criterion_val <- if (!is.null(criterion)) {
    criterion
  } else if (!is.null(template_obj) && !is.null(template_obj$criterion)) {
    template_obj$criterion
  } else {
    "projected_H"
  }

  repair_origin_lambda <- if (!is.null(template_obj) && !is.null(template_obj$repair_origin_lambda)) {
    template_obj$repair_origin_lambda
  } else if (!is.null(template_obj) && !is.null(template_obj$lambda)) {
    template_obj$lambda
  } else {
    NA_real_
  }

  repair_origin_source <- if (!is.null(template_obj) && !is.null(template_obj$repair_origin_source)) {
    template_obj$repair_origin_source
  } else if (!is.null(template_obj) && !is.null(template_obj$source)) {
    template_obj$source
  } else {
    NA_character_
  }

  repair_origin_index <- if (!is.null(template_obj) && !is.null(template_obj$repair_origin_index)) {
    template_obj$repair_origin_index
  } else if (!is.null(template_obj) && !is.null(template_obj$source_index)) {
    template_obj$source_index
  } else {
    NA_integer_
  }

  list(
    lambda = lambda,
    partition = partition,
    fit = fit_obj,
    projection = projection_obj,
    source = source_val,
    source_index = source_index_val,
    init_origin = init_origin_val,
    path_direction = path_direction_val,
    criterion = criterion_val,
    H_projected = eval_obj$H_total,
    H_class_projected = eval_obj$H_class,
    H_cond_projected = eval_obj$H_cond,
    H_raw = eval_obj$H_total,
    H_class_raw = eval_obj$H_class,
    H_cond_raw = eval_obj$H_cond,
    REO = r,
    repaired = TRUE,
    propagated_from_lambda = if (!is.null(template_obj) && !is.null(template_obj$lambda)) template_obj$lambda else NA_real_,
    repair_origin_lambda = repair_origin_lambda,
    repair_origin_source = repair_origin_source,
    repair_origin_index = repair_origin_index
  )
}

#' Repair an incoherent best-partition trajectory.
#'
#' @rdname CECextractBestPartitions
#' @param max_iter Maximum number of forward/backward propagation rounds used by
#'   the repair heuristic.
#' @return `CECrepairBestPartitionTrajectory()` returns an updated
#'   `best_partitions_obj` enriched with `repair_log`, coherence summaries
#'   before/after repair, and convergence information.
#' @export
CECrepairBestPartitionTrajectory <- function(
  best_partitions_obj,
  Z,
  tol = 1e-10,
  max_iter = 100
) {
  best_full <- best_partitions_obj$best
  ok <- vapply(best_full, function(x) !is.null(x), logical(1))
  valid_positions <- which(ok)

  if (length(valid_positions) == 0) {
    out <- best_partitions_obj
    out$repair_log <- data.frame()
    out$coherence_before <- CECcheckBestPartitionCoherence(best_partitions_obj, Z, tol = tol)
    out$coherence_after <- out$coherence_before
    out$n_repairs <- 0L
    out$converged <- TRUE
    out$tol <- tol
    out$max_iter <- max_iter
    return(out)
  }

  best_valid <- best_full[valid_positions]
  lambda_valid <- vapply(best_valid, function(x) x$lambda, numeric(1))
  ord_valid <- order(lambda_valid)
  best_sorted <- best_valid[ord_valid]
  lambda_sorted <- lambda_valid[ord_valid]
  n <- length(best_sorted)

  if (n <= 1) {
    out <- best_partitions_obj
    if (!is.null(out$summary) && n == 1) {
      row_new <- CECbestPartitionSummaryRow(best_sorted[[1]])
      row_new$repaired <- FALSE
      row_new$propagated_from_lambda <- NA_real_
      row_new$repair_origin_lambda <- row_new$lambda
      for (nm in names(row_new)) {
        if (!(nm %in% names(out$summary))) out$summary[[nm]] <- NA
      }
      out$summary[valid_positions[ord_valid], names(row_new)] <- row_new[1, names(row_new), drop = FALSE]
    }
    out$repair_log <- data.frame()
    out$coherence_before <- CECcheckBestPartitionCoherence(best_partitions_obj, Z, tol = tol)
    out$coherence_after <- out$coherence_before
    out$n_repairs <- 0L
    out$converged <- TRUE
    out$tol <- tol
    out$max_iter <- max_iter
    return(out)
  }

  get_fit_param <- function(obj, name, default = NULL) {
    if (!is.null(obj$fit) && !is.null(obj$fit$params) && !is.null(obj$fit$params[[name]])) {
      return(obj$fit$params[[name]])
    }
    default
  }

  familyType_vec <- vapply(
    best_sorted,
    function(x) {
      val <- get_fit_param(x, "familyType", default = NA_character_)
      if (length(val) == 0) NA_character_ else as.character(val[1])
    },
    character(1)
  )

  C_vec <- vapply(
    best_sorted,
    function(x) {
      val <- get_fit_param(x, "C", default = NA_real_)
      if (length(val) == 0) NA_real_ else as.numeric(val[1])
    },
    numeric(1)
  )

  familyType_default <- familyType_vec[which(!is.na(familyType_vec))[1]]
  C_default <- C_vec[which(!is.na(C_vec))[1]]

  if (length(familyType_default) == 0 || is.na(familyType_default)) {
    stop("Unable to infer familyType from best_partitions_obj.")
  }
  if (length(C_default) == 0 || is.na(C_default)) {
    stop("Unable to infer C from best_partitions_obj.")
  }

  get_slot_family <- function(idx) {
    if (is.na(familyType_vec[idx])) familyType_default else familyType_vec[idx]
  }

  get_slot_C <- function(idx) {
    if (is.na(C_vec[idx])) C_default else C_vec[idx]
  }

  coherence_before <- CECcheckBestPartitionCoherence(best_partitions_obj, Z, tol = tol)

  cache_list <- vector("list", n)

  evaluate_partition_at_slot <- function(partition, idx) {
    partition <- as.integer(as.factor(partition))
    bucket <- cache_list[[idx]]

    if (!is.null(bucket)) {
      for (j in seq_along(bucket$partitions)) {
        if (identical(bucket$partitions[[j]], partition)) {
          return(bucket$results[[j]])
        }
      }
    }

    out <- CECevaluatePartitionOnLambda(
      partition = partition,
      Z = Z,
      lambda = lambda_sorted[idx],
      C = get_slot_C(idx),
      familyType = get_slot_family(idx)
    )

    if (is.null(bucket)) {
      bucket <- list(partitions = list(), results = list())
    }
    bucket$partitions[[length(bucket$partitions) + 1L]] <- partition
    bucket$results[[length(bucket$results) + 1L]] <- out
    cache_list[[idx]] <<- bucket

    out
  }

  current_slot_H <- function(obj, idx) {
    h <- .scalar_or_na(obj$H_projected)
    if (!is.na(h)) {
      return(h)
    }
    evaluate_partition_at_slot(obj$partition, idx)$H_total
  }

  repair_log <- list()
  converged <- TRUE

  for (iter in seq_len(max_iter)) {
    changed_iter <- FALSE

    for (k in seq_len(n - 1L)) {
      eval_right <- evaluate_partition_at_slot(best_sorted[[k]]$partition, k + 1L)
      H_new <- eval_right$H_total
      H_old <- current_slot_H(best_sorted[[k + 1L]], k + 1L)

      if (isTRUE(H_new < H_old - tol)) {
        old_obj <- best_sorted[[k + 1L]]
        best_sorted[[k + 1L]] <- CECbuildBestPartitionObjectAtLambda(
          partition = best_sorted[[k]]$partition,
          lambda = lambda_sorted[k + 1L],
          Z = Z,
          C = get_slot_C(k + 1L),
          familyType = get_slot_family(k + 1L),
          template_obj = best_sorted[[k]],
          eval_obj = eval_right,
          criterion = if (!is.null(best_partitions_obj$criterion)) best_partitions_obj$criterion else NULL
        )

        repair_log[[length(repair_log) + 1L]] <- data.frame(
          iteration = iter,
          direction = "right",
          from_index = k,
          to_index = k + 1L,
          from_lambda = lambda_sorted[k],
          to_lambda = lambda_sorted[k + 1L],
          old_H = H_old,
          new_H = H_new,
          improvement = H_old - H_new,
          replaced_lambda = .scalar_or_na(old_obj$lambda),
          propagated_from_lambda = .scalar_or_na(best_sorted[[k]]$lambda)
        )

        changed_iter <- TRUE
      }
    }

    for (k in n:2) {
      eval_left <- evaluate_partition_at_slot(best_sorted[[k]]$partition, k - 1L)
      H_new <- eval_left$H_total
      H_old <- current_slot_H(best_sorted[[k - 1L]], k - 1L)

      if (isTRUE(H_new < H_old - tol)) {
        old_obj <- best_sorted[[k - 1L]]
        best_sorted[[k - 1L]] <- CECbuildBestPartitionObjectAtLambda(
          partition = best_sorted[[k]]$partition,
          lambda = lambda_sorted[k - 1L],
          Z = Z,
          C = get_slot_C(k - 1L),
          familyType = get_slot_family(k - 1L),
          template_obj = best_sorted[[k]],
          eval_obj = eval_left,
          criterion = if (!is.null(best_partitions_obj$criterion)) best_partitions_obj$criterion else NULL
        )

        repair_log[[length(repair_log) + 1L]] <- data.frame(
          iteration = iter,
          direction = "left",
          from_index = k,
          to_index = k - 1L,
          from_lambda = lambda_sorted[k],
          to_lambda = lambda_sorted[k - 1L],
          old_H = H_old,
          new_H = H_new,
          improvement = H_old - H_new,
          replaced_lambda = .scalar_or_na(old_obj$lambda),
          propagated_from_lambda = .scalar_or_na(best_sorted[[k]]$lambda)
        )

        changed_iter <- TRUE
      }
    }

    if (!changed_iter) {
      break
    }

    if (iter == max_iter) {
      converged <- FALSE
    }
  }

  repaired_valid_original_order <- vector("list", length(best_valid))
  repaired_valid_original_order[ord_valid] <- best_sorted

  best_out <- best_full
  best_out[valid_positions] <- repaired_valid_original_order

  summary_out <- if (!is.null(best_partitions_obj$summary)) {
    best_partitions_obj$summary
  } else {
    data.frame(lambda = vapply(best_out, function(x) if (is.null(x)) NA_real_ else x$lambda, numeric(1)))
  }

  for (pos in valid_positions) {
    row_new <- CECbestPartitionSummaryRow(best_out[[pos]])
    for (nm in names(row_new)) {
      if (!(nm %in% names(summary_out))) summary_out[[nm]] <- NA
    }
    summary_out[pos, names(row_new)] <- row_new[1, names(row_new), drop = FALSE]
  }

  coherence_after <- CECcheckBestPartitionCoherence(
    list(best = best_out, summary = summary_out),
    Z = Z,
    tol = tol
  )

  out <- best_partitions_obj
  out$best <- best_out
  out$summary <- summary_out
  out$repair_log <- if (length(repair_log) > 0) do.call(rbind, repair_log) else data.frame()
  out$coherence_before <- coherence_before
  out$coherence_after <- coherence_after
  out$n_repairs <- if (length(repair_log) > 0) nrow(out$repair_log) else 0L
  out$converged <- converged
  out$tol <- tol
  out$max_iter <- max_iter
  out
}


# ------------------------------------------------------------
# 1D plot of best partitions individually on histogram
# useful for inspecting change points
# ------------------------------------------------------------
#' Plot extracted one-dimensional partitions on a histogram.
#'
#' @rdname plotCECdiagnoseInteractive
#' @param best_partitions_obj Object returned by [CECextractBestPartitions()].
#' @param lambda_subset Optional subset of lambda values to display.
#' @param max_plots Maximum number of panels shown in the histogram overview.
#' @param digits_lambda,digits_H Number of digits used in panel annotations.
#' @return `plotCECBestPartitions1D()` is called for its side effects and returns
#'   `NULL` invisibly.
#' @export
plotCECBestPartitions1D <- function(
  Z,
  best_partitions_obj,
  lambda_subset = NULL,
  max_plots = 6,
  digits_lambda = 3,
  digits_H = 4
) {
  if (!is.numeric(Z)) stop("This plot is only for 1D numeric Z.")

  best_list <- best_partitions_obj$best
  ok <- vapply(best_list, function(x) !is.null(x), logical(1))
  best_list <- best_list[ok]

  if (!is.null(lambda_subset)) {
    keep <- vapply(best_list, function(x) x$lambda %in% lambda_subset, logical(1))
    best_list <- best_list[keep]
  }

  if (length(best_list) == 0) stop("No partition to plot.")

  if (length(best_list) > max_plots) {
    idx <- round(seq(1, length(best_list), length.out = max_plots))
    best_list <- best_list[idx]
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  nr <- ceiling(length(best_list) / 2)
  par(mfrow = c(nr, 2), mar = c(4, 4, 5, 1))

  for (obj in best_list) {
    cl <- obj$partition
    K <- length(unique(cl))
    cols <- seq_len(K)

    Htot  <- .scalar_or_na(obj$H_projected)
    Hnu   <- .scalar_or_na(obj$H_class_projected)
    Hcond <- .scalar_or_na(obj$H_cond_projected)
    REO   <- .scalar_or_na(obj$REO, default = length(unique(cl)))
    lam   <- .scalar_or_na(obj$lambda)

    main_txt <- paste0(
      "lambda = ", round(lam, digits_lambda),
      " | REO = ", REO,
      "\nH = ", round(Htot, digits_H),
      " ; H(nu) = ", round(Hnu, digits_H),
      " ; sum nu_x H(G_x|g_theta) = ", round(Hcond, digits_H)
    )

    hist(
      Z,
      breaks = 40,
      freq = FALSE,
      col = "grey85",
      border = "white",
      main = main_txt,
      xlab = "Z"
    )

    y0 <- par("usr")[3]
    points(
      Z,
      rep(y0, length(Z)),
      col = cols[as.integer(as.factor(cl))],
      pch = 20,
      cex = 0.6
    )

    if (!is.null(obj$fit) && !is.null(obj$fit$params$m)) {
      abline(
        v = obj$fit$params$m,
        lty = 2,
        col = cols[seq_along(obj$fit$params$m)]
      )
    }
  }
}

# ------------------------------------------------------------
# Petit helper : force un scalaire, sinon NA
# ------------------------------------------------------------
.scalar_or_na <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(default)
  x[1]
}

.best_finite_index <- function(values) {
  valid <- which(is.finite(values))
  if (length(valid) == 0L) return(integer(0))
  valid[which.min(values[valid])]
}

# ------------------------------------------------------------
# Decompose one fit on a given data set.
# ------------------------------------------------------------
CECdecompose_fit_on_data <- function(fit, Z) {
  if (is.null(fit)) return(NULL)
  
  dec <- evalCompositeEntropyDecomposed(
    phi = fit$phi,
    Z = Z,
    lambda = fit$params$lambda,
    C = fit$params$C,
    familyType = fit$params$familyType
  )
  
  list(
    H_total = .scalar_or_na(dec$H_total),
    H_class = .scalar_or_na(dec$H_class),
    H_cond  = .scalar_or_na(dec$H_cond)
  )
}

#' Extract one representative partition for each lambda in a diagnostic object.
#'
#' `CECextractBestPartitions()` condenses the output of
#' [CECdiagnose_lambda_grid_linked()] into one representative partition per
#' lambda value. Downstream helpers can then assess coherence, repair
#' trajectories, and identify change points along the selected lambda path.
#'
#' @rdname CECextractBestPartitions
#' @param lambda_diag Object returned by [CECdiagnose_lambda_grid_linked()].
#' @param source Candidate pool used when selecting the representative partition:
#'   `"fixed"` uses runs on the original data, `"bootstrap"` uses projected
#'   bootstrap fits, and `"all"` combines both.
#' @param Z Original data set. It is required for coherence checks and repair
#'   operations, and recommended when extracting partitions so that raw criterion
#'   decompositions can be recovered when available.
#' @return `CECextractBestPartitions()` returns a list with:
#' - `best`: a list of per-lambda representative partitions;
#' - `summary`: a compact data frame with criterion values, REO, and origin
#'   labels;
#' - `source` and `criterion`: the selection settings used.
#'
#' @examples
#' CECconfigure_runtime("base")
#' Z <- simulate_multidim_benchmark_data(n = 80, p_num = 3, p_fac = 2, seed = 2)
#' lambda_diag <- CECdiagnose_lambda_grid_linked(
#'   Z = Z,
#'   lambda_grid = seq(0.2, 0.8, by = 0.2),
#'   k0 = 2,
#'   B = 2,
#'   C = 10,
#'   r0 = 6,
#'   Nshots_fresh = 1,
#'   Nshots_warm = 1,
#'   Nloop = 12,
#'   familyType = "gaussAndDiscreteVector",
#'   seed = 2,
#'   silent = TRUE,
#'   verbose = FALSE,
#'   n_cores = 1L,
#'   checkpoint_dir = FALSE,
#'   auto_checkpoint = FALSE,
#'   show_progress = FALSE
#' )
#'
#' best_parts <- CECextractBestPartitions(
#'   lambda_diag,
#'   source = "fixed",
#'   criterion = "projected_H",
#'   Z = Z
#' )
#' coherence <- CECcheckBestPartitionCoherence(best_parts, Z)
#' changes <- CECdetectPartitionChanges(best_parts)
#'
#' head(best_parts$summary)
#' coherence$all_coherent
#' changes$change_lambdas
#' @export
CECextractBestPartitions <- function(
  lambda_diag,
  source = c("all", "fixed", "bootstrap"),
  criterion = c("projected_H", "raw_H"),
  Z = NULL
) {
  source <- match.arg(source)
  criterion <- match.arg(criterion)
  
  details <- lambda_diag$details
  out_best <- vector("list", length(details))
  out_rows <- vector("list", length(details))
  
  for (i in seq_along(details)) {
    obj <- details[[i]]
    lam <- obj$lambda
    
    candidates <- list()
    
    # -------------------------
    # Candidats issus des runs sur Z
    # -------------------------
    if (source %in% c("all", "fixed")) {
      for (j in seq_along(obj$fits0)) {
        fit_j  <- obj$fits0[[j]]
        proj_j <- obj$projections0[[j]]
        
        if (!is.null(fit_j) && !is.null(proj_j)) {
          raw_dec <- tryCatch(
            CECdecompose_fit_on_data(fit_j, Z = Z),
            error = function(e) NULL
          )
          
          candidates[[length(candidates) + 1]] <- list(
            lambda = lam,
            partition = proj_j$clusters,
            fit = fit_j,
            projection = proj_j,
            source = "fixed",
            source_index = j,
            init_origin = CECget_fit_origin_label(fit_j),
            path_direction = CECget_fit_path_direction(fit_j),
            H_projected = .scalar_or_na(proj_j$H),
            H_class_projected = .scalar_or_na(proj_j$H_class),
            H_cond_projected = .scalar_or_na(proj_j$H_cond),
            H_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_total else NULL),
            H_class_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_class else NULL),
            H_cond_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_cond else NULL),
            REO = .scalar_or_na(fit_j$REO, default = NA_integer_)
          )
        }
      }
    }
    
    # -------------------------
    # Candidats issus des bootstraps
    # -------------------------
    if (source %in% c("all", "bootstrap")) {
      for (j in seq_along(obj$fitsB)) {
        fit_j  <- obj$fitsB[[j]]
        proj_j <- obj$projectionsB[[j]]
        
        if (!is.null(fit_j) && !is.null(proj_j)) {
          raw_dec <- NULL
          
          if (!is.null(Z) && !is.null(obj$boot_indices[[j]])) {
            ib <- obj$boot_indices[[j]]
            Zb <- if (is.matrix(Z) || is.data.frame(Z)) Z[ib, , drop = FALSE] else Z[ib]
            
            raw_dec <- tryCatch(
              CECdecompose_fit_on_data(fit_j, Z = Zb),
              error = function(e) NULL
            )
          }
          
          candidates[[length(candidates) + 1]] <- list(
            lambda = lam,
            partition = proj_j$clusters,
            fit = fit_j,
            projection = proj_j,
            source = "bootstrap",
            source_index = j,
            init_origin = CECget_fit_origin_label(fit_j),
            path_direction = CECget_fit_path_direction(fit_j),
            H_projected = .scalar_or_na(proj_j$H),
            H_class_projected = .scalar_or_na(proj_j$H_class),
            H_cond_projected = .scalar_or_na(proj_j$H_cond),
            H_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_total else NULL),
            H_class_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_class else NULL),
            H_cond_raw = .scalar_or_na(if (!is.null(raw_dec)) raw_dec$H_cond else NULL),
            REO = .scalar_or_na(fit_j$REO, default = NA_integer_)
          )
        }
      }
    }
    
    # aucun candidat valide
    if (length(candidates) == 0) {
      out_best[[i]] <- NULL
      out_rows[[i]] <- data.frame(
        lambda = lam,
        source = NA_character_,
        source_index = NA_integer_,
        init_origin = NA_character_,
        path_direction = NA_character_,
        criterion = criterion,
        H_projected = NA_real_,
        H_class_projected = NA_real_,
        H_cond_projected = NA_real_,
        H_raw = NA_real_,
        H_class_raw = NA_real_,
        H_cond_raw = NA_real_,
        REO = NA_integer_
      )
      next
    }
    
    crit_values <- vapply(
      candidates,
      function(x) {
        if (criterion == "projected_H") .scalar_or_na(x$H_projected) else .scalar_or_na(x$H_raw)
      },
      numeric(1)
    )

    best_id <- .best_finite_index(crit_values)
    if (length(best_id) == 0L) {
      fallback_values <- vapply(
        candidates,
        function(x) {
          if (criterion == "projected_H") .scalar_or_na(x$H_raw) else .scalar_or_na(x$H_projected)
        },
        numeric(1)
      )
      best_id <- .best_finite_index(fallback_values)
    }
    if (length(best_id) == 0L) {
      out_best[[i]] <- NULL
      out_rows[[i]] <- data.frame(
        lambda = lam,
        source = NA_character_,
        source_index = NA_integer_,
        init_origin = NA_character_,
        path_direction = NA_character_,
        criterion = criterion,
        H_projected = NA_real_,
        H_class_projected = NA_real_,
        H_cond_projected = NA_real_,
        H_raw = NA_real_,
        H_class_raw = NA_real_,
        H_cond_raw = NA_real_,
        REO = NA_integer_
      )
      next
    }
    best <- candidates[[best_id]]
    out_best[[i]] <- best
    
    out_rows[[i]] <- data.frame(
      lambda = .scalar_or_na(best$lambda),
      source = as.character(.scalar_or_na(best$source, default = NA_character_)),
      source_index = .scalar_or_na(best$source_index, default = NA_integer_),
      init_origin = as.character(.scalar_or_na(best$init_origin, default = NA_character_)),
      path_direction = as.character(.scalar_or_na(best$path_direction, default = NA_character_)),
      criterion = criterion,
      H_projected = .scalar_or_na(best$H_projected),
      H_class_projected = .scalar_or_na(best$H_class_projected),
      H_cond_projected = .scalar_or_na(best$H_cond_projected),
      H_raw = .scalar_or_na(best$H_raw),
      H_class_raw = .scalar_or_na(best$H_class_raw),
      H_cond_raw = .scalar_or_na(best$H_cond_raw),
      REO = .scalar_or_na(best$REO, default = NA_integer_)
    )
  }
  
  summary_df <- do.call(rbind, out_rows)
  rownames(summary_df) <- NULL
  
  list(
    best = out_best,
    summary = summary_df,
    source = source,
    criterion = criterion
  )
}





CECfit_linked_one_lambda <- function(
  Z,
  lambda,
  phi_prev = NULL,
  direction = c("forward", "backward"),
  C = 1,
  r0 = NULL,
  Nloop = 300,
  Nshots_warm = 5,
  Nshots_fresh = 5,
  familyType = "gaussAndDiscreteVector",
  sizeMaxOutlier = 0,
  autoRegroupOutliers = FALSE,
  focus = NULL,
  silent = TRUE,
  backend_data = NULL
) {
  direction <- match.arg(direction)
  if (is.null(backend_data)) {
    backend_data <- CECprepare_backend_data(Z, familyType = familyType)
  }
  candidates <- list()

  # candidat warm-start multi-shot
  if (!is.null(phi_prev) && Nshots_warm > 0) {
    fit_warm <- tryCatch(
      CECfit_warm_multi(
        Z = Z,
        lambda = lambda,
        phi0 = phi_prev,
        C = C,
        r0 = r0,
        Nshots_warm = Nshots_warm,
        Nloop = Nloop,
        familyType = familyType,
        sizeMaxOutlier = sizeMaxOutlier,
        autoRegroupOutliers = autoRegroupOutliers,
        path_direction = direction,
        silent = silent,
        backend_data = backend_data
      ),
      error = function(e) NULL
    )
    candidates[[length(candidates) + 1]] <- fit_warm
  }

  # Small independent candidate to avoid tracking a poor branch.
  if (Nshots_fresh > 0) {
    fit_fresh <- CECfit_safe(
      Z = Z,
      lambda = lambda,
      C = C,
      r0 = r0,
      Nshots = Nshots_fresh,
      Nloop = Nloop,
      familyType = familyType,
      sizeMaxOutlier = sizeMaxOutlier,
      autoRegroupOutliers = autoRegroupOutliers,
      focus = focus,
      silent = silent,
      backend_data = backend_data
    )
    fit_fresh <- CECtag_fit_origin(
      fit_fresh,
      init_origin = "fresh",
      path_direction = direction,
      candidate_group = "fresh",
      candidate_index = 1L
    )
    candidates[[length(candidates) + 1]] <- fit_fresh
  }

  ok <- !vapply(candidates, is.null, logical(1))
  if (!any(ok)) return(NULL)

  candidates <- candidates[ok]
  Hvals <- vapply(candidates, function(x) x$Hphi, numeric(1))
  best_idx <- .best_finite_index(Hvals)
  if (length(best_idx) == 0L) return(NULL)
  candidates[[best_idx]]
}

CECfollowLambdaPath <- function(
  Z,
  lambda_grid,
  direction = c("forward", "backward"),
  C = 1,
  r0 = NULL,
  Nshots_fresh = 5,
  Nshots_warm = 5,
  Nshots_init = NULL,
  Nloop = 300,
  familyType = "gaussAndDiscreteVector",
  sizeMaxOutlier = 0,
  autoRegroupOutliers = FALSE,
  focus = NULL,
  silent = TRUE,
  backend_data = NULL
) {
  direction <- match.arg(direction)
  if (is.null(backend_data)) {
    backend_data <- CECprepare_backend_data(Z, familyType = familyType)
  }

  if (!is.null(Nshots_init) && missing(Nshots_fresh)) {
    Nshots_fresh <- Nshots_init
  }
  if (!is.null(Nshots_init) && missing(Nshots_warm)) {
    Nshots_warm <- Nshots_init
  }

  if (Nshots_fresh < 1) {
    stop("Nshots_fresh must be >= 1 to initialize a lambda path.")
  }

  if (Nshots_warm < 0) {
    stop("Nshots_warm must be >= 0.")
  }

  lambda_grid <- sort(unique(lambda_grid))
  L <- length(lambda_grid)

  ord <- if (direction == "forward") seq_len(L) else rev(seq_len(L))

  fits <- vector("list", L)
  phi_prev <- NULL

  for (jj in seq_along(ord)) {
    i <- ord[jj]
    lam <- lambda_grid[i]

    if (is.null(phi_prev)) {
      fit_i <- CECfit_safe(
        Z = Z,
        lambda = lam,
        C = C,
        r0 = r0,
        Nshots = Nshots_fresh,
        Nloop = Nloop,
        familyType = familyType,
        sizeMaxOutlier = sizeMaxOutlier,
        autoRegroupOutliers = autoRegroupOutliers,
        focus = focus,
        silent = silent,
        backend_data = backend_data
      )
      fit_i <- CECtag_fit_origin(
        fit_i,
        init_origin = "fresh",
        path_direction = direction,
        candidate_group = "fresh",
        candidate_index = 1L
      )
    } else {
      fit_i <- CECfit_linked_one_lambda(
        Z = Z,
        lambda = lam,
        phi_prev = phi_prev,
        direction = direction,
        C = C,
        r0 = r0,
        Nloop = Nloop,
        Nshots_warm = Nshots_warm,
        Nshots_fresh = Nshots_fresh,
        familyType = familyType,
        sizeMaxOutlier = sizeMaxOutlier,
        autoRegroupOutliers = autoRegroupOutliers,
        focus = focus,
        silent = silent,
        backend_data = backend_data
      )
    }

    fits[i] <- list(fit_i)
    if (!is.null(fit_i)) {
      phi_prev <- fit_i$phi
    }
  }

  list(
    lambda_grid = lambda_grid,
    direction = direction,
    fits = fits
  )
}

# ------------------------------------------------------------
# Main entry point for lambda diagnostics based on linked paths
# ------------------------------------------------------------
#' Diagnose a grid of lambda values with linked forward and backward paths.
#'
#' `CECdiagnose_lambda_grid_linked()` is the main workflow for studying how the
#' fitted partition evolves across a sequence of regularisation values. For each
#' `lambda`, the function combines repeated runs on the original data with
#' bootstrap trajectories projected back on the reference data set. The result
#' can be checkpointed and resumed, which is helpful for long multi-core runs.
#'
#' @param Z Input data used for all fits and projections.
#' @param lambda_grid Numeric vector of candidate `lambda` values. Duplicates are
#'   removed and the grid is sorted increasingly.
#' @param k0 Number of linked paths fitted on the original data.
#' @param B Number of linked bootstrap paths.
#' @param C Positive upper bound for fitted component densities; independent
#'   of `lambda`.
#' @param r0 Optional upper bound for the initial number of clusters.
#' @param Nshots_fresh Number of fresh random starts used at the first lambda of
#'   each path.
#' @param Nshots_warm Number of warm starts used for neighbouring lambda values
#'   along a path.
#' @param Nshots_init Deprecated compatibility alias that can still initialise
#'   `Nshots_fresh` and `Nshots_warm`.
#' @param Nloop Maximum number of optimisation iterations per fit.
#' @param familyType Model family passed to [CECclassif()].
#' @param sizeMaxOutlier,autoRegroupOutliers,focus Additional controls forwarded
#'   to [CECclassif()].
#' @param replace Logical. If `TRUE`, bootstrap samples are drawn with
#'   replacement.
#' @param seed Optional random seed used for reproducibility.
#' @param silent,verbose Logical controls for fit-level and workflow-level
#'   messaging.
#' @param n_cores Number of worker processes used for independent path tasks.
#' @param checkpoint_dir Directory used to store per-task `.rds` checkpoints.
#'   Use `TRUE` to request an automatically generated temporary folder and
#'   `FALSE` to disable checkpointing explicitly.
#' @param auto_checkpoint Logical. If `NULL`, checkpointing is enabled
#'   automatically when `n_cores > 1`.
#' @param resume Logical. If `TRUE`, completed tasks already present in
#'   `checkpoint_dir` are reused.
#' @param batch_size Number of tasks launched together in each parallel batch.
#' @param show_progress Logical. If `TRUE`, progress lines are printed while
#'   tasks complete.
#' @param rerun_failed_serial Logical. If `TRUE`, unfinished tasks are retried
#'   sequentially after a parallel batch failure.
#' @param checkpoint_compress Logical passed to `saveRDS()` for checkpoint
#'   files.
#' @param partition_metric Similarity measure used for partition comparisons.
#'   `"match"` relies on an optimal label matching; `"ARI"` uses the adjusted
#'   Rand index.
#' @param stability_approach Strategy used to summarise bootstrap stability.
#'
#' @return `CECdiagnose_lambda_grid_linked()` returns a list with:
#' - `summary`: one row per lambda with stability and criterion summaries;
#' - `details`: the full per-lambda objects used to build downstream summaries;
#' - `lambda_grid`: the analysed lambda grid;
#' - `partition_metric` and `stability_approach`: the comparison settings used;
#' - `linked`: a compact description of the linked fixed and bootstrap paths;
#' - `execution`: runtime metadata, including checkpoint information.
#'
#' @examples
#' Z <- simulate_multidim_benchmark_data(n = 80, p_num = 3, p_fac = 2, seed = 1)
#' runtime <- CECconfigure_runtime("base")
#'
#' lambda_diag <- CECdiagnose_lambda_grid_linked(
#'   Z = Z,
#'   lambda_grid = seq(0.2, 0.8, by = 0.2),
#'   k0 = 2,
#'   B = 2,
#'   C = 10,
#'   r0 = 6,
#'   Nshots_fresh = 1,
#'   Nshots_warm = 1,
#'   Nloop = 12,
#'   familyType = "gaussAndDiscreteVector",
#'   seed = 1,
#'   silent = TRUE,
#'   verbose = FALSE,
#'   n_cores = runtime$n_cores,
#'   checkpoint_dir = FALSE,
#'   auto_checkpoint = FALSE,
#'   show_progress = FALSE
#' )
#'
#' lambda_diag <- CECaddSmoothedDiagnostics(lambda_diag, k = 3)
#' stable <- CECidentifyStableLambdas(lambda_diag, rule = "ratio")
#'
#' head(lambda_diag$summary)
#' stable$lambda_min
#' @export
CECdiagnose_lambda_grid_linked <- function(
  Z,
  lambda_grid = seq(0.1, 2, by = 0.1),
  k0 = 10,
  B = 10,
  C = 1,
  r0 = NULL,
  Nshots_fresh = 5,
  Nshots_warm = 5,
  Nshots_init = NULL,
  Nloop = 300,
  familyType = "gaussAndDiscreteVector",
  sizeMaxOutlier = 0,
  autoRegroupOutliers = FALSE,
  focus = NULL,
  replace = TRUE,
  seed = NULL,
  silent = TRUE,
  verbose = TRUE,
  n_cores = 1L,
  checkpoint_dir = NULL,
  auto_checkpoint = NULL,
  resume = TRUE,
  batch_size = NULL,
  show_progress = NULL,
  rerun_failed_serial = TRUE,
  checkpoint_compress = FALSE,
  partition_metric = c("match", "ARI"),
  stability_approach = c("pairwise", "reference_bestH")
) {
  partition_metric <- match.arg(partition_metric)
  stability_approach <- match.arg(stability_approach)
  n_cores <- max(1L, as.integer(n_cores[1]))
  if (is.null(auto_checkpoint)) {
    auto_checkpoint <- n_cores > 1L
  }
  if (is.null(show_progress)) {
    show_progress <- isTRUE(verbose)
  }
  if (is.null(batch_size)) {
    batch_size <- min(n_cores, max(1L, k0 + B))
  }
  batch_size <- max(1L, as.integer(batch_size[1]))

  if (!is.null(Nshots_init) && missing(Nshots_fresh)) {
    Nshots_fresh <- Nshots_init
  }
  if (!is.null(Nshots_init) && missing(Nshots_warm)) {
    Nshots_warm <- Nshots_init
  }

  if (!is.null(seed)) set.seed(seed)

  lambda_grid <- sort(unique(lambda_grid))
  L <- length(lambda_grid)
  fixed_backend_data <- CECprepare_backend_data(Z, familyType = familyType)

  task_meta <- CECdiagnostics_make_task_meta(
    Z = Z,
    lambda_grid = lambda_grid,
    k0 = k0,
    B = B,
    C = C,
    r0 = r0,
    Nshots_fresh = Nshots_fresh,
    Nshots_warm = Nshots_warm,
    Nloop = Nloop,
    familyType = familyType,
    sizeMaxOutlier = sizeMaxOutlier,
    autoRegroupOutliers = autoRegroupOutliers,
    focus = focus,
    replace = replace,
    seed = seed
  )

  checkpoint_info <- CECdiagnostics_prepare_checkpoint(
    checkpoint_dir = checkpoint_dir,
    task_meta = task_meta,
    auto_checkpoint = auto_checkpoint,
    resume = resume
  )

  tasks <- CECdiagnostics_build_tasks(
    k0 = k0,
    B = B,
    checkpoint_dir = checkpoint_info$dir,
    checkpoint_compress = checkpoint_compress,
    seed = seed
  )
  all_tasks <- tasks$all
  all_task_keys <- vapply(all_tasks, CECdiagnostics_task_key, character(1))

  existing_completed_keys <- if (checkpoint_info$enabled) {
    CECdiagnostics_completed_task_keys(all_tasks)
  } else {
    character(0)
  }

  if (checkpoint_info$enabled && !isTRUE(resume) && length(existing_completed_keys) > 0L) {
    stop(
      "checkpoint_dir already contains completed tasks. Use resume = TRUE or choose another checkpoint_dir."
    )
  }

  completed_keys <- if (checkpoint_info$enabled && isTRUE(resume)) unique(existing_completed_keys) else character(0)
  pending_tasks <- Filter(
    function(task) !(CECdiagnostics_task_key(task) %in% completed_keys),
    all_tasks
  )

  progress_state <- CECdiagnostics_new_progress_state(k0 = k0, B = B, show_progress = show_progress)
  CECdiagnostics_set_progress_from_keys(progress_state, completed_keys, all_tasks)

  if (isTRUE(show_progress)) {
    if (checkpoint_info$enabled) {
      cat("CECdiagnose checkpoints: ", checkpoint_info$dir, "\n", sep = "")
    }
    if (length(completed_keys) > 0L) {
      CECdiagnostics_print_progress(progress_state, prefix = "CECdiagnose resume")
    } else {
      CECdiagnostics_print_progress(progress_state, prefix = "CECdiagnose start")
    }
  }

  task_context <- list(
    Z = Z,
    lambda_grid = lambda_grid,
    C = C,
    r0 = r0,
    Nshots_fresh = Nshots_fresh,
    Nshots_warm = Nshots_warm,
    Nloop = Nloop,
    familyType = familyType,
    sizeMaxOutlier = sizeMaxOutlier,
    autoRegroupOutliers = autoRegroupOutliers,
    focus = focus,
    replace = replace,
    silent = silent,
    fixed_backend_data = fixed_backend_data
  )

  runtime_results <- list()

  if (length(pending_tasks) > 0L) {
    if (n_cores > 1L) {
      task_batches <- split(
        pending_tasks,
        ceiling(seq_along(pending_tasks) / batch_size)
      )

      for (batch_idx in seq_along(task_batches)) {
        batch_tasks <- task_batches[[batch_idx]]
        if (isTRUE(show_progress)) {
          cat(
            "CECdiagnose batch ", batch_idx, "/", length(task_batches),
            " (", length(batch_tasks), " task(s))\n",
            sep = ""
          )
        }

        batch_res <- tryCatch(
          CECrun_diagnostic_tasks_parallel_batch(
            tasks = batch_tasks,
            context = task_context,
            n_cores = n_cores,
            progress_state = progress_state
          ),
          error = function(e) e
        )

        if (inherits(batch_res, "error")) {
          if (isTRUE(show_progress)) {
            cat("CECdiagnose parallel batch failed: ", conditionMessage(batch_res), "\n", sep = "")
          }

          if (checkpoint_info$enabled) {
            recovered_keys <- setdiff(
              CECdiagnostics_completed_task_keys(batch_tasks),
              completed_keys
            )
            if (length(recovered_keys) > 0L) {
              completed_keys <- unique(c(completed_keys, recovered_keys))
              CECdiagnostics_set_progress_from_keys(progress_state, completed_keys, all_tasks)
              CECdiagnostics_print_progress(progress_state, prefix = "CECdiagnose recovered")
            }
          }

          remaining_batch_tasks <- Filter(
            function(task) !(CECdiagnostics_task_key(task) %in% completed_keys),
            batch_tasks
          )

          if (length(remaining_batch_tasks) > 0L) {
            if (!isTRUE(rerun_failed_serial)) {
              stop(batch_res)
            }
            if (isTRUE(show_progress)) {
              cat(
                "CECdiagnose fallback to sequential for ",
                length(remaining_batch_tasks),
                " remaining task(s)\n",
                sep = ""
              )
            }
            serial_res <- CECrun_diagnostic_tasks_sequential(
              tasks = remaining_batch_tasks,
              context = task_context,
              progress_state = progress_state,
              completed_keys = completed_keys
            )
            completed_keys <- unique(c(completed_keys, serial_res$completed_keys))
            if (!checkpoint_info$enabled) {
              runtime_results[names(serial_res$results)] <- serial_res$results
            }
          }
        } else {
          completed_keys <- unique(c(completed_keys, batch_res$completed_keys))
          if (!checkpoint_info$enabled) {
            runtime_results[names(batch_res$results)] <- batch_res$results
          }
        }
      }
    } else {
      serial_res <- CECrun_diagnostic_tasks_sequential(
        tasks = pending_tasks,
        context = task_context,
        progress_state = progress_state,
        completed_keys = completed_keys
      )
      completed_keys <- unique(c(completed_keys, serial_res$completed_keys))
      runtime_results[names(serial_res$results)] <- serial_res$results
    }
  }

  if (checkpoint_info$enabled) {
    completed_keys <- unique(c(completed_keys, CECdiagnostics_completed_task_keys(all_tasks)))
  }

  missing_keys <- setdiff(all_task_keys, completed_keys)
  if (length(missing_keys) > 0L) {
    stop(
      "Some diagnostic tasks did not complete: ",
      paste(head(missing_keys, 10L), collapse = ", ")
    )
  }

  # stockage par lambda
  fits0_by_lambda <- vector("list", L)
  proj0_by_lambda <- vector("list", L)
  parts0_by_lambda <- vector("list", L)
  H0_by_lambda <- vector("list", L)
  H0raw_by_lambda <- vector("list", L)
  REO0_by_lambda <- vector("list", L)
  K0_by_lambda <- vector("list", L)

  fitsB_by_lambda <- vector("list", L)
  projB_by_lambda <- vector("list", L)
  partsB_by_lambda <- vector("list", L)
  HB_by_lambda <- vector("list", L)
  HBraw_by_lambda <- vector("list", L)
  REOB_by_lambda <- vector("list", L)
  KB_by_lambda <- vector("list", L)
  boot_indices_by_lambda <- vector("list", L)

  for (i in seq_len(L)) {
    fits0_by_lambda[[i]] <- vector("list", k0)
    proj0_by_lambda[[i]] <- vector("list", k0)
    parts0_by_lambda[[i]] <- vector("list", k0)
    H0_by_lambda[[i]] <- rep(NA_real_, k0)
    H0raw_by_lambda[[i]] <- rep(NA_real_, k0)
    REO0_by_lambda[[i]] <- rep(NA_integer_, k0)
    K0_by_lambda[[i]] <- rep(NA_integer_, k0)

    fitsB_by_lambda[[i]] <- vector("list", B)
    projB_by_lambda[[i]] <- vector("list", B)
    partsB_by_lambda[[i]] <- vector("list", B)
    HB_by_lambda[[i]] <- rep(NA_real_, B)
    HBraw_by_lambda[[i]] <- rep(NA_real_, B)
    REOB_by_lambda[[i]] <- rep(NA_integer_, B)
    KB_by_lambda[[i]] <- rep(NA_integer_, B)
    boot_indices_by_lambda[[i]] <- vector("list", B)
  }

  if (checkpoint_info$enabled) {
    fixed_results <- lapply(tasks$fixed, function(task) {
      if (is.null(task$result_file) || !file.exists(task$result_file)) {
        return(NULL)
      }
      readRDS(task$result_file)
    })
    names(fixed_results) <- vapply(
      tasks$fixed,
      function(task) if (is.null(task$result_file)) "" else task$result_file,
      character(1)
    )
    boot_results <- lapply(tasks$boot, function(task) {
      if (is.null(task$result_file) || !file.exists(task$result_file)) {
        return(NULL)
      }
      readRDS(task$result_file)
    })
    names(boot_results) <- vapply(
      tasks$boot,
      function(task) if (is.null(task$result_file)) "" else task$result_file,
      character(1)
    )
    expected_partition_length <- CECget_n_obs(Z)
    CECdiagnostics_validate_result_lengths(fixed_results, expected_partition_length, "fixed")
    CECdiagnostics_validate_result_lengths(boot_results, expected_partition_length, "bootstrap")
  } else {
    fixed_results <- lapply(tasks$fixed, function(task) {
      key <- CECdiagnostics_task_key(task)
      res <- runtime_results[[key]]
      if (is.null(res) || !isTRUE(res$ok)) {
        return(NULL)
      }
      res$result
    })
    boot_results <- lapply(tasks$boot, function(task) {
      key <- CECdiagnostics_task_key(task)
      res <- runtime_results[[key]]
      if (is.null(res) || !isTRUE(res$ok)) {
        return(NULL)
      }
      res$result
    })
  }

  for (res in fixed_results) {
    rep <- res$rep
    for (i in seq_along(lambda_grid)) {
      slot_i <- res$per_lambda[[i]]
      fits0_by_lambda[[i]][rep] <- list(slot_i$fit)
      proj0_by_lambda[[i]][rep] <- list(slot_i$proj)
      parts0_by_lambda[[i]][rep] <- list(slot_i$part)
      H0_by_lambda[[i]][rep] <- slot_i$H
      H0raw_by_lambda[[i]][rep] <- slot_i$Hraw
      REO0_by_lambda[[i]][rep] <- slot_i$REO
      K0_by_lambda[[i]][rep] <- slot_i$K
    }
  }

  for (res in boot_results) {
    b <- res$b
    for (i in seq_along(lambda_grid)) {
      slot_i <- res$per_lambda[[i]]
      fitsB_by_lambda[[i]][b] <- list(slot_i$fit)
      projB_by_lambda[[i]][b] <- list(slot_i$proj)
      partsB_by_lambda[[i]][b] <- list(slot_i$part)
      HB_by_lambda[[i]][b] <- slot_i$H
      HBraw_by_lambda[[i]][b] <- slot_i$Hraw
      REOB_by_lambda[[i]][b] <- slot_i$REO
      KB_by_lambda[[i]][b] <- slot_i$K
      boot_indices_by_lambda[[i]][b] <- list(res$boot_indices)
    }
  }

  details <- vector("list", L)
  rows <- vector("list", L)

  for (i in seq_along(lambda_grid)) {
    if (stability_approach == "pairwise") {
      stab0 <- partition_stability_summary(
        parts0_by_lambda[[i]],
        partition_metric = partition_metric
      )
      stabB <- partition_stability_summary(
        partsB_by_lambda[[i]],
        partition_metric = partition_metric
      )
      stab_ratio_mean <- .safe_ratio(stabB$mean_sim, stab0$mean_sim)
      stab_ratio_medoid <- .safe_ratio(stabB$medoid_mean_sim, stab0$medoid_mean_sim)
    } else {
      stab0 <- partition_reference_stability_summary(
        partitions = parts0_by_lambda[[i]],
        H = H0_by_lambda[[i]],
        partition_metric = partition_metric
      )
      stabB <- partition_reference_stability_summary(
        partitions = partsB_by_lambda[[i]],
        H = HB_by_lambda[[i]],
        partition_metric = partition_metric
      )
      stab_ratio_mean <- .safe_ratio(stabB$mean_sim, stab0$mean_sim)
      stab_ratio_medoid <- NA_real_
    }

    obj_i <- list(
      lambda = lambda_grid[i],
      partition_metric = partition_metric,
      stability_approach = stability_approach,

      fits0 = fits0_by_lambda[[i]],
      projections0 = proj0_by_lambda[[i]],
      partitions0 = parts0_by_lambda[[i]],
      stability0 = stab0,
      H0 = H0_by_lambda[[i]],
      H0_raw = H0raw_by_lambda[[i]],
      REO0 = REO0_by_lambda[[i]],
      Kproj0 = K0_by_lambda[[i]],

      fitsB = fitsB_by_lambda[[i]],
      projectionsB = projB_by_lambda[[i]],
      partitionsB = partsB_by_lambda[[i]],
      stabilityB = stabB,
      HB = HB_by_lambda[[i]],
      HB_raw = HBraw_by_lambda[[i]],
      REOB = REOB_by_lambda[[i]],
      KprojB = KB_by_lambda[[i]],
      boot_indices = boot_indices_by_lambda[[i]],

      stab_ratio_mean = stab_ratio_mean,
      stab_ratio_medoid = stab_ratio_medoid
    )

    details[[i]] <- obj_i
    rows[[i]] <- CECdiagnose_lambda_summary_row(obj_i)
  }

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL

  list(
    summary = summary_df,
    details = details,
    lambda_grid = lambda_grid,
    partition_metric = partition_metric,
    stability_approach = stability_approach,
    linked = TRUE,
    execution = list(
      n_cores = n_cores,
      batch_size = batch_size,
      checkpoint_enabled = checkpoint_info$enabled,
      checkpoint_dir = checkpoint_info$dir,
      auto_checkpoint = checkpoint_info$auto,
      resume = isTRUE(resume)
    )
  )
}

#' Build a PCA-ready data frame describing partition trajectories.
#'
#' @rdname plotCECdiagnoseInteractive
#' @param best_parts Object returned by [CECextractBestPartitions()].
#' @param scale_data Logical. If `TRUE`, numeric variables are scaled before PCA.
#' @return `build_partition_path_df()` returns a long data frame with PCA
#'   coordinates, lambda values, cluster labels, and optional species labels for
#'   faceted plotting across lambdas.
#' @export
build_partition_path_df <- function(Z, best_parts, scale_data = TRUE, lambda_subset = NULL) {
  Z_df <- as.data.frame(Z)
  num_cols <- vapply(Z_df, is.numeric, logical(1))
  if (sum(num_cols) < 2L) {
    stop("Z must contain at least two numeric columns for PCA path visualisation.")
  }

  Z_num <- as.matrix(Z_df[, num_cols, drop = FALSE])
  n <- nrow(Z_num)

  pca <- stats::prcomp(Z_num, scale. = scale_data)
  coords <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  colnames(coords) <- c("PC1", "PC2")
  coords$id <- seq_len(n)

  if ("Species" %in% names(Z_df)) {
    coords$species <- Z_df$Species
  }

  best_list <- best_parts$best
  ok <- vapply(best_list, function(x) {
    !is.null(x) && length(x$partition) == n
  }, logical(1))

  if (!any(ok)) {
    stop("No valid partition has length equal to nrow(Z).")
  }

  best_list <- best_list[ok]

  if (!is.null(lambda_subset)) {
    keep <- vapply(best_list, function(x) x$lambda %in% lambda_subset, logical(1))
    best_list <- best_list[keep]
  }

  if (length(best_list) == 0) {
    stop("No valid partition left after filtering by lambda_subset.")
  }

  out <- vector("list", length(best_list))

  for (i in seq_along(best_list)) {
    obj <- best_list[[i]]
    df_i <- coords
    df_i$lambda <- obj$lambda
    df_i$cluster <- factor(obj$partition)
    df_i$REO <- obj$REO
    out[[i]] <- df_i
  }

  do.call(rbind, out)
}

#' Plot PCA coordinates of extracted partitions across lambdas.
#'
#' @rdname plotCECdiagnoseInteractive
#' @param path_df Data frame returned by [build_partition_path_df()].
#' @return `plot_partition_path_pca()` returns a `ggplot2` object.
#' @export
plot_partition_path_pca <- function(path_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for PCA path plots. Install it.")
  }
  PC1 <- PC2 <- cluster <- species <- NULL

  if ("species" %in% names(path_df)) {
    ggplot2::ggplot(path_df, ggplot2::aes(PC1, PC2, color = cluster, shape = species)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::facet_wrap(~ lambda) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "Partition path along lambda",
        x = "PC1",
        y = "PC2",
        color = "Cluster",
        shape = "Species"
      )
  } else {
    ggplot2::ggplot(path_df, ggplot2::aes(PC1, PC2, color = cluster)) +
      ggplot2::geom_point(size = 2, alpha = 0.85) +
      ggplot2::facet_wrap(~ lambda) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "Partition path along lambda",
        x = "PC1",
        y = "PC2",
        color = "Cluster"
      )
  }
}

CECresolve_best_parts_for_plot <- function(
  obj,
  Z = NULL,
  source = c("fixed", "bootstrap", "all"),
  criterion = c("projected_H", "raw_H")
) {
  source <- match.arg(source)
  criterion <- match.arg(criterion)

  if (is.list(obj) && !is.null(obj$best) && !is.null(obj$summary)) {
    return(obj)
  }

  if (is.list(obj) && !is.null(obj$details) && !is.null(obj$summary)) {
    if (is.null(Z)) {
      stop("Z must be provided when plotting directly from a lambda_diag object.")
    }
    return(
      CECextractBestPartitions(
        lambda_diag = obj,
        source = source,
        criterion = criterion,
        Z = Z
      )
    )
  }

  stop("obj must be either a best_partitions object or a lambda_diag object.")
}

CECinteractive_lambda_key <- function(lambda) {
  sprintf("%.12g", as.numeric(lambda))
}

CECinteractive_lambda_label <- function(lambda) {
  paste0("lambda = ", formatC(as.numeric(lambda), format = "fg", digits = 4))
}

CECinteractive_valid_best_list <- function(best_parts, n) {
  best_list <- best_parts$best
  ok <- vapply(best_list, function(x) {
    !is.null(x) && !is.null(x$partition) && length(x$partition) == n
  }, logical(1))

  if (!any(ok)) {
    stop("No valid partition has length equal to nrow(Z).")
  }

  best_list <- best_list[ok]
  names(best_list) <- vapply(best_list, function(x) CECinteractive_lambda_key(x$lambda), character(1))
  best_list
}

CECinteractive_selected_best_list <- function(best_list, lambda_keys = NULL) {
  if (is.null(lambda_keys)) {
    return(best_list)
  }
  keep <- names(best_list) %in% lambda_keys
  best_list[keep]
}

CECinteractive_extract_mean_for_var <- function(fit, var_name) {
  if (is.null(fit) || is.null(fit$params$m)) {
    return(numeric(0))
  }

  m <- fit$params$m
  if (is.null(dim(m))) {
    if (identical(var_name, names(fit$params$m)[1])) {
      return(as.numeric(m))
    }
    return(as.numeric(m))
  }

  if (!is.null(colnames(m)) && var_name %in% colnames(m)) {
    return(as.numeric(m[, var_name]))
  }

  if (ncol(m) == 1L) {
    return(as.numeric(m[, 1]))
  }

  numeric(0)
}

CECinteractive_sort_labels <- function(x) {
  x <- unique(as.character(x))
  if (length(x) == 0L) {
    return(character(0))
  }

  x_num <- suppressWarnings(as.numeric(x))
  if (all(!is.na(x_num))) {
    return(x[order(x_num)])
  }

  sort(x)
}

CECinteractive_cluster_levels <- function(best_list) {
  CECinteractive_sort_labels(unlist(lapply(best_list, function(x) x$partition), use.names = FALSE))
}

CECinteractive_cluster_palette <- function(cluster_levels) {
  cluster_levels <- as.character(cluster_levels)
  if (length(cluster_levels) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }

  hues <- seq(15, 375, length.out = length(cluster_levels) + 1L)
  cols <- grDevices::hcl(
    h = hues[seq_len(length(cluster_levels))],
    l = 65,
    c = 100
  )
  stats::setNames(cols, cluster_levels)
}

CECinteractive_shape_values <- function(shape_levels) {
  shape_levels <- as.character(shape_levels)
  if (length(shape_levels) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }

  base_shapes <- c(0:20, 22:25)
  stats::setNames(rep(base_shapes, length.out = length(shape_levels)), shape_levels)
}

CECinteractive_add_manual_scales <- function(p, cluster_colors = NULL, shape_values = NULL) {
  if (!is.null(cluster_colors) && length(cluster_colors) > 0L) {
    p <- p + ggplot2::scale_color_manual(values = cluster_colors)
  }

  if (!is.null(shape_values) && length(shape_values) > 0L) {
    p <- p + ggplot2::scale_shape_manual(values = shape_values)
  }

  p
}

CECinteractive_filter_cluster_data <- function(data, cluster_filter = "__all__") {
  if (is.null(data) || identical(cluster_filter, "__all__")) {
    return(data)
  }

  data[data$cluster == cluster_filter, , drop = FALSE]
}

CECinteractive_single_lambda_clusters <- function(best_list) {
  if (length(best_list) != 1L) {
    return(character(0))
  }

  CECinteractive_sort_labels(best_list[[1]]$partition)
}

CECinteractive_default_lambda_keys <- function(best_list) {
  if (length(best_list) == 0L) {
    return(character(0))
  }

  lambda_values <- suppressWarnings(as.numeric(vapply(best_list, function(x) x$lambda, numeric(1))))
  if (anyNA(lambda_values)) {
    return(names(best_list)[1])
  }

  names(best_list)[order(lambda_values, seq_along(lambda_values))][1]
}

CECinteractive_add_cluster_weight_to_title <- function(base_title, data, cluster_filter = "__all__") {
  if (
    is.null(data) ||
    identical(cluster_filter, "__all__") ||
    !"cluster" %in% names(data) ||
    nrow(data) == 0L
  ) {
    return(base_title)
  }

  cluster_chr <- as.character(data$cluster)
  cluster_filter_chr <- as.character(cluster_filter)
  cluster_n <- sum(cluster_chr == cluster_filter_chr, na.rm = TRUE)
  total_n <- length(cluster_chr)

  if (total_n < 1L || cluster_n < 1L) {
    return(base_title)
  }

  weight_pct <- 100 * cluster_n / total_n
  paste0(
    base_title,
    " - cluster ",
    cluster_filter_chr,
    " (",
    formatC(weight_pct, format = "f", digits = 1),
    "%)"
  )
}

CECinteractive_build_cluster_bar_data <- function(
  data,
  qual_var = NULL,
  cluster_filter = "__all__",
  shape_levels = NULL
) {
  if (
    is.null(data) ||
    is.null(qual_var) ||
    !"shape_value" %in% names(data)
  ) {
    return(NULL)
  }

  data_cluster <- CECinteractive_filter_cluster_data(data, cluster_filter = cluster_filter)
  if (nrow(data_cluster) == 0L) {
    return(NULL)
  }

  if (is.null(shape_levels)) {
    shape_levels <- levels(factor(data_cluster$shape_value))
  }

  counts <- table(factor(as.character(data_cluster$shape_value), levels = shape_levels))
  total_n <- sum(counts)

  data.frame(
    shape_value = factor(shape_levels, levels = shape_levels),
    n = as.numeric(counts),
    pct = if (total_n > 0) 100 * as.numeric(counts) / total_n else 0,
    label = if (total_n > 0) sprintf("%.1f%%", 100 * as.numeric(counts) / total_n) else "0.0%",
    stringsAsFactors = FALSE
  )
}

CECinteractive_plot_cluster_bar <- function(
  bar_data,
  qual_var,
  cluster_filter = "__all__",
  fill_color = "#4D4D4D"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for interactive plotting. Install it.")
  }
  shape_value <- n <- label <- NULL

  title_text <- if (identical(cluster_filter, "__all__")) {
    paste("Distribution of", qual_var, "across all clusters")
  } else {
    paste("Distribution of", qual_var, "in cluster", cluster_filter)
  }

  ggplot2::ggplot(
    bar_data,
    ggplot2::aes(x = shape_value, y = n)
  ) +
    ggplot2::geom_col(fill = fill_color, alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.35, size = 3.5) +
    ggplot2::coord_cartesian(ylim = c(0, max(1, bar_data$n, na.rm = TRUE) * 1.15)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    ) +
    ggplot2::labs(
      title = title_text,
      x = qual_var,
      y = "Count"
    )
}

CECinteractive_draw_plot_stack <- function(top_plot, bottom_plot = NULL, top_ratio = 0.72) {
  if (is.null(bottom_plot)) {
    print(top_plot)
    return(invisible(NULL))
  }

  top_ratio <- max(0.1, min(0.9, as.numeric(top_ratio)[1]))
  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 2,
        ncol = 1,
        heights = grid::unit(c(top_ratio, 1 - top_ratio), "null")
      )
    )
  )
  print(top_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(bottom_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  invisible(NULL)
}

CECinteractive_build_hist_plot_data <- function(
  Z,
  best_list,
  xvar,
  qual_var = NULL,
  bins = 30L
) {
  Z_df <- as.data.frame(Z)
  x <- Z_df[[xvar]]
  n <- length(x)
  cluster_levels <- CECinteractive_cluster_levels(best_list)
  cluster_colors <- CECinteractive_cluster_palette(cluster_levels)
  shape_levels <- NULL
  shape_values <- NULL

  if (!is.null(qual_var)) {
    shape_levels <- if (is.factor(Z_df[[qual_var]])) {
      levels(Z_df[[qual_var]])
    } else {
      levels(factor(Z_df[[qual_var]]))
    }
    shape_values <- CECinteractive_shape_values(shape_levels)
  }

  hist_template <- hist(x, breaks = bins, plot = FALSE)
  max_count <- max(1, hist_template$counts, na.rm = TRUE)
  base_y <- -0.04 * max_count

  point_data <- vector("list", length(best_list))
  hist_data <- vector("list", length(best_list))
  mean_data <- vector("list", length(best_list))

  for (i in seq_along(best_list)) {
    obj <- best_list[[i]]
    lambda_label <- CECinteractive_lambda_label(obj$lambda)
    clusters <- factor(obj$partition)
    point_y <- base_y + (((seq_len(n) - 1L) %% 7L) - 3L) / 7L * abs(base_y) * 0.25

    df_i <- data.frame(
      value = x,
      y = point_y,
      lambda = lambda_label,
      cluster = factor(clusters, levels = cluster_levels),
      stringsAsFactors = FALSE
    )

    if (!is.null(qual_var)) {
      df_i$shape_value <- factor(Z_df[[qual_var]], levels = shape_levels)
    }
    point_data[[i]] <- df_i

    hist_data[[i]] <- data.frame(
      value = x,
      lambda = lambda_label,
      stringsAsFactors = FALSE
    )

    means_i <- CECinteractive_extract_mean_for_var(obj$fit, xvar)
    if (length(means_i) > 0L) {
      mean_data[[i]] <- data.frame(
        value = means_i,
        lambda = lambda_label,
        cluster = factor(seq_along(means_i), levels = cluster_levels),
        stringsAsFactors = FALSE
      )
    } else {
      mean_data[[i]] <- NULL
    }
  }

  list(
    hist = do.call(rbind, hist_data),
    points = do.call(rbind, point_data),
    means = if (any(vapply(mean_data, Negate(is.null), logical(1)))) do.call(rbind, mean_data) else NULL,
    xlim = range(x, na.rm = TRUE),
    ylim = c(base_y * 1.6, max_count * 1.05),
    cluster_levels = cluster_levels,
    cluster_colors = cluster_colors,
    shape_levels = shape_levels,
    shape_values = shape_values
  )
}

CECinteractive_plot_histogram <- function(
  plot_data,
  xvar,
  qual_var = NULL,
  bins = 30L,
  cluster_filter = "__all__"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for interactive plotting. Install it.")
  }
  value <- y <- cluster <- shape_value <- NULL

  points_plot <- CECinteractive_filter_cluster_data(
    plot_data$points,
    cluster_filter = cluster_filter
  )
  means_plot <- CECinteractive_filter_cluster_data(
    plot_data$means,
    cluster_filter = cluster_filter
  )

  if (is.null(qual_var)) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = plot_data$hist,
        ggplot2::aes(x = value),
        bins = bins,
        fill = "grey85",
        color = "white"
      ) +
      ggplot2::geom_point(
        data = points_plot,
        ggplot2::aes(x = value, y = y, color = cluster),
        size = 1.5,
        alpha = 0.8
      )
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = plot_data$hist,
        ggplot2::aes(x = value),
        bins = bins,
        fill = "grey85",
        color = "white"
      ) +
      ggplot2::geom_point(
        data = points_plot,
        ggplot2::aes(x = value, y = y, color = cluster, shape = shape_value),
        size = 1.7,
        alpha = 0.85
      )
  }

  if (!is.null(means_plot) && nrow(means_plot) > 0L) {
    p <- p + ggplot2::geom_vline(
      data = means_plot,
      ggplot2::aes(xintercept = value, color = cluster),
      linetype = 2,
      linewidth = 0.7,
      alpha = 0.7,
      show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::facet_wrap(~ lambda, scales = "free_y") +
    ggplot2::coord_cartesian(xlim = plot_data$xlim, ylim = plot_data$ylim) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = CECinteractive_add_cluster_weight_to_title(
        base_title = paste("Partition path on", xvar),
        data = plot_data$points,
        cluster_filter = cluster_filter
      ),
      x = xvar,
      y = "Count",
      color = "Cluster",
      shape = if (is.null(qual_var)) NULL else qual_var
    )

  CECinteractive_add_manual_scales(
    p,
    cluster_colors = plot_data$cluster_colors,
    shape_values = plot_data$shape_values
  )
}

CECinteractive_build_2d_plot_data <- function(
  Z,
  best_list,
  xvar = NULL,
  yvar = NULL,
  qual_var = NULL,
  mode = c("pair", "pca"),
  scale_pca = TRUE
) {
  mode <- match.arg(mode)
  Z_df <- as.data.frame(Z)
  num_cols <- names(Z_df)[vapply(Z_df, is.numeric, logical(1))]
  cluster_levels <- CECinteractive_cluster_levels(best_list)
  cluster_colors <- CECinteractive_cluster_palette(cluster_levels)
  shape_levels <- NULL
  shape_values <- NULL

  if (!is.null(qual_var)) {
    shape_levels <- if (is.factor(Z_df[[qual_var]])) {
      levels(Z_df[[qual_var]])
    } else {
      levels(factor(Z_df[[qual_var]]))
    }
    shape_values <- CECinteractive_shape_values(shape_levels)
  }

  if (mode == "pair") {
    if (is.null(xvar) || is.null(yvar)) {
      stop("xvar and yvar must be provided in pair mode.")
    }
    coords <- data.frame(
      x = Z_df[[xvar]],
      y = Z_df[[yvar]],
      stringsAsFactors = FALSE
    )
    xlab <- xvar
    ylab <- yvar
    title <- paste("Partition path in the plane", paste0("(", xvar, ", ", yvar, ")"))
  } else {
    if (length(num_cols) < 2L) {
      stop("PCA display requires at least two numeric variables.")
    }
    pca <- stats::prcomp(as.matrix(Z_df[, num_cols, drop = FALSE]), scale. = scale_pca)
    coords <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    names(coords) <- c("x", "y")
    xlab <- "PC1"
    ylab <- "PC2"
    title <- "Partition path on the first two PCA axes"
  }

  out <- vector("list", length(best_list))
  for (i in seq_along(best_list)) {
    obj <- best_list[[i]]
    df_i <- coords
    df_i$lambda <- CECinteractive_lambda_label(obj$lambda)
    df_i$cluster <- factor(obj$partition, levels = cluster_levels)
    if (!is.null(qual_var)) {
      df_i$shape_value <- factor(Z_df[[qual_var]], levels = shape_levels)
    }
    out[[i]] <- df_i
  }

  list(
    data = do.call(rbind, out),
    xlab = xlab,
    ylab = ylab,
    title = title,
    xlim = range(coords$x, na.rm = TRUE),
    ylim = range(coords$y, na.rm = TRUE),
    cluster_levels = cluster_levels,
    cluster_colors = cluster_colors,
    shape_levels = shape_levels,
    shape_values = shape_values
  )
}

CECinteractive_plot_2d <- function(
  plot_data,
  qual_var = NULL,
  cluster_filter = "__all__"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for interactive plotting. Install it.")
  }
  x <- y <- cluster <- shape_value <- NULL

  data_plot <- CECinteractive_filter_cluster_data(
    plot_data$data,
    cluster_filter = cluster_filter
  )

  if (is.null(qual_var)) {
    p <- ggplot2::ggplot(
      data_plot,
      ggplot2::aes(x = x, y = y, color = cluster)
    )
  } else {
    p <- ggplot2::ggplot(
      data_plot,
      ggplot2::aes(x = x, y = y, color = cluster, shape = shape_value)
    )
  }

  p <- p +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~ lambda) +
    ggplot2::coord_cartesian(xlim = plot_data$xlim, ylim = plot_data$ylim) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = CECinteractive_add_cluster_weight_to_title(
        base_title = plot_data$title,
        data = plot_data$data,
        cluster_filter = cluster_filter
      ),
      x = plot_data$xlab,
      y = plot_data$ylab,
      color = "Cluster",
      shape = if (is.null(qual_var)) NULL else qual_var
    )

  CECinteractive_add_manual_scales(
    p,
    cluster_colors = plot_data$cluster_colors,
    shape_values = plot_data$shape_values
  )
}

#' Launch an interactive explorer for CEClust partitions.
#'
#' @rdname plotCECdiagnoseInteractive
#' @param obj Either a `lambda_diag` object returned by
#'   [CECdiagnose_lambda_grid_linked()] or a `best_parts` object returned by
#'   [CECextractBestPartitions()].
#' @param Z Original data set used to build the diagnostic object.
#' @param source Candidate pool used when `obj` is a `lambda_diag` object and
#'   best partitions still need to be extracted.
#' @param criterion Criterion used when extracting best partitions from a
#'   `lambda_diag` object.
#' @param launch.browser Logical passed to `shiny::runApp()`.
#' @param scale_pca Logical. If `TRUE`, numeric variables are scaled before PCA
#'   is computed inside the interactive application.
#' @details `plotCECdiagnoseInteractive()` requires the optional packages
#'   \pkg{shiny} and \pkg{ggplot2}. For one-dimensional numeric data, the app
#'   switches automatically to a histogram-based display.
#' @return `plotCECdiagnoseInteractive()` launches a local Shiny application and
#'   returns the value produced by `shiny::runApp()`.
#'
#' @examples
#' Z <- simulate_multidim_benchmark_data(n = 80, p_num = 3, p_fac = 2, seed = 3)
#' lambda_diag <- CECdiagnose_lambda_grid_linked(
#'   Z = Z,
#'   lambda_grid = seq(0.2, 0.8, by = 0.2),
#'   k0 = 2,
#'   B = 2,
#'   C = 10,
#'   r0 = 6,
#'   Nshots_fresh = 1,
#'   Nshots_warm = 1,
#'   Nloop = 12,
#'   familyType = "gaussAndDiscreteVector",
#'   seed = 3,
#'   silent = TRUE,
#'   verbose = FALSE,
#'   n_cores = 1L,
#'   checkpoint_dir = FALSE,
#'   auto_checkpoint = FALSE,
#'   show_progress = FALSE
#' )
#' best_parts <- CECextractBestPartitions(lambda_diag, source = "fixed", Z = Z)
#'
#' plotCECdiagnoseLambdaGrid(lambda_diag, use_smoothed = FALSE)
#'
#' path_df <- build_partition_path_df(iris[, -5], CECextractBestPartitions(
#'   CECdiagnose_lambda_grid_linked(
#'     iris[, -5],
#'     lambda_grid = seq(0.2, 0.4, by = 0.2),
#'     k0 = 2,
#'     B = 1,
#'     C = 10,
#'     r0 = 4,
#'     Nshots_fresh = 1,
#'     Nshots_warm = 1,
#'     Nloop = 10,
#'     familyType = "gaussVector",
#'     seed = 1,
#'     silent = TRUE,
#'     verbose = FALSE,
#'     n_cores = 1L,
#'     checkpoint_dir = FALSE,
#'     auto_checkpoint = FALSE,
#'     show_progress = FALSE
#'   ),
#'   source = "fixed",
#'   Z = iris[, -5]
#' ))
#' plot_partition_path_pca(path_df)
#'
#' @examplesIf interactive() && requireNamespace("shiny", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
#' plotCECdiagnoseInteractive(best_parts, Z, launch.browser = FALSE)
#' @export
plotCECdiagnoseInteractive <- function(
  obj,
  Z,
  source = c("fixed", "bootstrap", "all"),
  criterion = c("projected_H", "raw_H"),
  launch.browser = TRUE,
  scale_pca = TRUE
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' required for interactive plotting. Install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for interactive plotting. Install it.")
  }

  source <- match.arg(source)
  criterion <- match.arg(criterion)

  best_parts <- CECresolve_best_parts_for_plot(
    obj = obj,
    Z = Z,
    source = source,
    criterion = criterion
  )

  Z_df <- as.data.frame(Z)
  n <- nrow(Z_df)
  best_list <- CECinteractive_valid_best_list(best_parts, n = n)

  num_cols <- names(Z_df)[vapply(Z_df, is.numeric, logical(1))]
  fac_cols <- names(Z_df)[vapply(Z_df, is.factor, logical(1))]

  if (length(num_cols) < 1L) {
    stop("Interactive plotting currently requires at least one numeric variable in Z.")
  }

  lambda_keys <- names(best_list)
  lambda_labels <- vapply(best_list, function(x) CECinteractive_lambda_label(x$lambda), character(1))
  names(lambda_labels) <- lambda_keys
  default_lambda_keys <- CECinteractive_default_lambda_keys(best_list)

  has_hist_mode <- length(num_cols) == 1L
  default_mode <- if (has_hist_mode) "hist" else "pair"
  default_x <- num_cols[1]
  default_y <- if (length(num_cols) >= 2L) num_cols[2] else num_cols[1]
  qual_choices <- c("(none)" = "__none__")
  if (length(fac_cols) > 0L) {
    qual_choices <- c(qual_choices, stats::setNames(fac_cols, fac_cols))
  }
  default_qual <- if (length(fac_cols) == 1L) fac_cols[1] else "__none__"

  app <- shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::titlePanel("CEClust interactive partition explorer"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::helpText("Color represents the cluster. A factor can be shown with point shapes."),
          shiny::selectizeInput(
            inputId = "lambda_keys",
            label = "Lambdas to display",
            choices = stats::setNames(lambda_keys, lambda_labels),
            selected = default_lambda_keys,
            multiple = TRUE,
            options = list(plugins = list("remove_button"))
          ),
          shiny::uiOutput("cluster_filter_ui"),
          if (has_hist_mode) {
            shiny::helpText("One numeric variable detected: histogram display is used.")
          } else {
            shiny::radioButtons(
              inputId = "display_mode",
              label = "Numeric display",
              choices = c(
                "Variable pair" = "pair",
                "PCA axes (PC1, PC2)" = "pca"
              ),
              selected = default_mode
            )
          },
          if (!has_hist_mode) {
            shiny::conditionalPanel(
              condition = "input.display_mode == 'pair'",
              shiny::selectInput("xvar", "X variable", choices = num_cols, selected = default_x),
              shiny::selectInput("yvar", "Y variable", choices = num_cols, selected = default_y)
            )
          },
          if (has_hist_mode) {
            shiny::sliderInput("bins", "Histogram bins", min = 10, max = 80, value = 30, step = 1)
          },
          shiny::selectInput(
            inputId = "qual_var",
            label = "Categorical variable shown by shape",
            choices = qual_choices,
            selected = default_qual
          ),
          shiny::uiOutput("detail_plot_ratio_ui")
        ),
        shiny::mainPanel(
          shiny::plotOutput("cec_partition_plot", height = "820px"),
          shiny::verbatimTextOutput("cec_partition_info")
        )
      )
    ),
    server = function(input, output, session) {
      if (!has_hist_mode) {
        shiny::observeEvent(input$xvar, {
          y_choices <- num_cols[num_cols != input$xvar]
          if (length(y_choices) == 0L) {
            y_choices <- input$xvar
          }
          selected_y <- if (!is.null(input$yvar) && input$yvar %in% y_choices) input$yvar else y_choices[1]
          shiny::updateSelectInput(session, "yvar", choices = y_choices, selected = selected_y)
        }, ignoreNULL = FALSE)
      }

      qual_var_reactive <- shiny::reactive({
        if (is.null(input$qual_var) || identical(input$qual_var, "__none__")) {
          return(NULL)
        }
        input$qual_var
      })

      selected_best_reactive <- shiny::reactive({
        shiny::validate(
          shiny::need(length(input$lambda_keys) > 0L, "Select at least one lambda.")
        )
        best_sel <- CECinteractive_selected_best_list(best_list, lambda_keys = input$lambda_keys)
        shiny::validate(
          shiny::need(length(best_sel) > 0L, "No valid partition available for the selected lambdas.")
        )
        best_sel
      })

      single_lambda_clusters_reactive <- shiny::reactive({
        best_sel <- selected_best_reactive()
        clusters <- CECinteractive_single_lambda_clusters(best_sel)
        if (length(clusters) == 0L) {
          return(NULL)
        }
        clusters
      })

      output$cluster_filter_ui <- shiny::renderUI({
        clusters <- single_lambda_clusters_reactive()
        if (is.null(clusters)) {
          return(NULL)
        }

        shiny::selectInput(
          inputId = "cluster_filter",
          label = "Cluster subset",
          choices = c("All clusters" = "__all__", stats::setNames(clusters, paste("Cluster", clusters))),
          selected = "__all__"
        )
      })

      selected_cluster_reactive <- shiny::reactive({
        clusters <- single_lambda_clusters_reactive()
        if (is.null(clusters)) {
          return("__all__")
        }

        value <- input$cluster_filter
        valid <- c("__all__", clusters)
        if (is.null(value) || !value %in% valid) {
          return("__all__")
        }

        value
      })

      output$detail_plot_ratio_ui <- shiny::renderUI({
        qual_var <- qual_var_reactive()
        best_sel <- selected_best_reactive()

        if (length(best_sel) != 1L || is.null(qual_var)) {
          return(NULL)
        }

        shiny::sliderInput(
          inputId = "detail_plot_ratio",
          label = "Upper plot height share",
          min = 0.35,
          max = 0.9,
          value = 0.72,
          step = 0.05
        )
      })

      output$cec_partition_plot <- shiny::renderPlot({
        best_sel <- selected_best_reactive()
        qual_var <- qual_var_reactive()
        cluster_filter <- if (length(best_sel) == 1L) selected_cluster_reactive() else "__all__"
        top_ratio <- if (is.null(input$detail_plot_ratio)) 0.72 else input$detail_plot_ratio
        bottom_plot <- NULL

        if (has_hist_mode) {
          plot_data <- CECinteractive_build_hist_plot_data(
            Z = Z_df,
            best_list = best_sel,
            xvar = num_cols[1],
            qual_var = qual_var,
            bins = input$bins
          )

          top_plot <- CECinteractive_plot_histogram(
            plot_data = plot_data,
            xvar = num_cols[1],
            qual_var = qual_var,
            bins = input$bins,
            cluster_filter = cluster_filter
          )

          if (!is.null(qual_var) && length(best_sel) == 1L) {
            bar_data <- CECinteractive_build_cluster_bar_data(
              data = plot_data$points,
              qual_var = qual_var,
              cluster_filter = cluster_filter,
              shape_levels = plot_data$shape_levels
            )
            if (!is.null(bar_data)) {
              bottom_plot <- CECinteractive_plot_cluster_bar(
                bar_data = bar_data,
                qual_var = qual_var,
                cluster_filter = cluster_filter,
                fill_color = if (identical(cluster_filter, "__all__")) {
                  "#9E9E9E"
                } else {
                  unname(plot_data$cluster_colors[[as.character(cluster_filter)]])
                }
              )
            }
          }
        } else {
          mode <- if (is.null(input$display_mode)) default_mode else input$display_mode
          plot_data <- CECinteractive_build_2d_plot_data(
            Z = Z_df,
            best_list = best_sel,
            xvar = if (identical(mode, "pair")) input$xvar else NULL,
            yvar = if (identical(mode, "pair")) input$yvar else NULL,
            qual_var = qual_var,
            mode = mode,
            scale_pca = scale_pca
          )

          top_plot <- CECinteractive_plot_2d(
            plot_data = plot_data,
            qual_var = qual_var,
            cluster_filter = cluster_filter
          )

          if (!is.null(qual_var) && length(best_sel) == 1L) {
            bar_data <- CECinteractive_build_cluster_bar_data(
              data = plot_data$data,
              qual_var = qual_var,
              cluster_filter = cluster_filter,
              shape_levels = plot_data$shape_levels
            )
            if (!is.null(bar_data)) {
              bottom_plot <- CECinteractive_plot_cluster_bar(
                bar_data = bar_data,
                qual_var = qual_var,
                cluster_filter = cluster_filter,
                fill_color = if (identical(cluster_filter, "__all__")) {
                  "#9E9E9E"
                } else {
                  unname(plot_data$cluster_colors[[as.character(cluster_filter)]])
                }
              )
            }
          }
        }

        CECinteractive_draw_plot_stack(
          top_plot = top_plot,
          bottom_plot = bottom_plot,
          top_ratio = top_ratio
        )
      })

      output$cec_partition_info <- shiny::renderPrint({
        qual_var <- qual_var_reactive()
        best_sel <- selected_best_reactive()
        cluster_filter <- if (length(best_sel) == 1L) selected_cluster_reactive() else "__all__"
        cat("Lambdas displayed:", length(best_sel), "\n")
        cat("Numeric variables:", paste(num_cols, collapse = ", "), "\n")
        cat(
          "Cluster subset:",
          if (identical(cluster_filter, "__all__")) "all clusters" else paste("cluster", cluster_filter),
          "\n"
        )
        cat(
          "Categorical variable displayed:",
          if (is.null(qual_var)) "(none)" else qual_var,
          "\n"
        )
      })
    }
  )

  esc_hint <- "Press Esc in the R console to return to the prompt and stop the Shiny app.\n"
  launch_browser_internal <- launch.browser
  if (isTRUE(launch.browser)) {
    launch_browser_internal <- function(url) {
      utils::browseURL(url)
      cat(esc_hint)
    }
  } else if (is.function(launch.browser)) {
    launch_browser_internal <- function(url) {
      launch.browser(url)
      cat(esc_hint)
    }
  }

  invisible(
    shiny::runApp(
      app,
      launch.browser = launch_browser_internal
    )
  )
}



# ------------------------------------------------------------
# Distance de partition par matching optimal
# ------------------------------------------------------------
partition_distance_match <- function(z_true, z_pred) {
  z_true <- as.integer(as.factor(z_true))
  z_pred <- as.integer(as.factor(z_pred))
  
  n <- length(z_true)
  stopifnot(length(z_pred) == n)
  
  tab <- table(z_true, z_pred)
  
  # solve_LSAP minimizes a cost, so convert gains to costs first.
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required. Install it with install.packages('clue').")
  }
  
  nr <- nrow(tab)
  nc <- ncol(tab)
  m  <- max(nr, nc)
  
  # Square gain matrix completed with zeros.
  gain <- matrix(0, nrow = m, ncol = m)
  gain[1:nr, 1:nc] <- tab
  
  # Cost = max(gain) - gain.
  cost <- max(gain) - gain
  
  perm <- clue::solve_LSAP(cost)
  matched <- sum(gain[cbind(seq_len(m), perm)])
  
  proximity <- matched / n
  distance  <- 1 - proximity
  
  list(
    distance = distance,
    proximity = proximity,
    matched = matched,
    contingency = tab
  )
}

compute_path_distance <- function(best_parts, partitions_other) {
  
  res <- data.frame(
    lambda = best_parts$summary$lambda,
    k_hat = NA_integer_,
    distance = NA_real_
  )
  
  for (i in seq_along(best_parts$best)) {
    
    obj <- best_parts$best[[i]]
    if (is.null(obj)) next
    
    z_hat <- obj$partition
    k <- length(unique(z_hat))
    
    res$k_hat[i] <- k
    
    key <- paste0("k", k)
    
    if (key %in% names(partitions_other)) {
      z_other <- partitions_other[[key]]
      res$distance[i] <- partition_distance_match(z_hat, z_other)$distance
    }
  }
  
  res
}

