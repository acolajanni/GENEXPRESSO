#' ebayesGEO
#'
#'
#' @export eBayesGEO



#  EMPIRICAL BAYES FUNCTIONS

eBayesGEO <- function(fit,proportion=0.01,stdev.coef.lim=c(0.1,4),trend=FALSE,robust=FALSE,winsor.tail.p=c(0.05,0.1),proportion_vector=0.01)
  #fit = fit2; proportion_vector=prop;
  #proportion=0.01;stdev.coef.lim=c(0.1,4);
  #trend=FALSE;robust=FALSE;winsor.tail.p=c(0.05,0.1);

  #	Empirical Bayes statistics to select differentially expressed genes
#	Object orientated version
#	Gordon Smyth
#	4 August 2003.  Last modified 20 November 2012.
{
	# reordering the proportion vector and filling NA's with
	# for unfound genes in the p-vector make there proportion the average
	if (class(proportion_vector) == "data.frame" || class(proportion_vector) == "matrix"){
		genes <- names(fit$sigma)
		p <- vector()
		if (dim(proportion_vector)[2] == 1){
			#p_avg <- median(proportion_vector[,1])
			p_avg <- 0.06932688
			print(p_avg)

			p = sapply(1:length(genes), function(x) {ifelse( genes[x] %in% rownames(proportion_vector), proportion_vector[genes[x],], p_avg )} )

			#for (i in 1:length(genes)){

			  #print(c("gene ID in expression matrix:", genes[i]))

			  #print(c(genes[i] %in% rownames(proportion_vector)))


				#if(genes[i] %in% rownames(proportion_vector)){

				  #print(c("Id in probability:", proportion_vector[genes[i],]))

					#p[i]<-proportion_vector[genes[i],]

				#} else{
					#p[i]<-p_avg
				#}
			#}

		} else if (dim(proportion_vector)[2] == 2){
			warning(paste("A matrix or dataframe with 2 columns the first column will be used as gene names.",
				"To avoid this in the future provide the proportions in a matrix, dataframe, or vector where the gene names are the row names.",
				"This might take sometime"))
			# needs work
			p_avg <- median(proportion_vector[,2])
			plist <- vector()
			for (i in 1:dim(p4)[1]){
				plist[as.character(p4[i,1])]<-p4[i,2]
			}
			for (i in 1:length(genes)){
				if(genes[i] %in% names(plist)){
					p[i]<-pplist[genes[i],]
				} else{
					p[i]<-p_avg
				}
			}
		} else{
			exit(paste("A matrix or dataframe with undeterminable columns was provided.",
				"To avoid this in the future provide the proportions in a matrix, dataframe, or vector where the gene names are the row names."))
		}
	} else if (class(proportion_vector) == "numeric" && length(proportion_vector) != 1){
		genes <- names(fit$sigma)
		p <- vector()
		p_avg <- median(proportion_vector)
		for (i in 1:length(genes)){
			if(genes[i] %in% names(proportion_vector)){
				p[i]<-proportion_vector[genes[i],]
			} else{
				p[i]<-p_avg
			}
		}
	} else{
		warning(paste("A matrix, dataframe, or named was not provided.",
				"To avoid this in the future provide the proportions in a matrix, dataframe, or vector where the gene names are the row names/names."))
		p <- proportion
	}

	if(trend) if(is.null(fit$Amean)) stop("Need Amean component in fit to estimate trend")
	eb <- ebayesGEO(fit=fit,proportion=proportion,stdev.coef.lim=stdev.coef.lim,trend=trend,robust=robust,winsor.tail.p=winsor.tail.p, proportion_vector=p)
	fit$df.prior <- eb$df.prior
	fit$s2.prior <- eb$s2.prior
	fit$var.prior <- eb$var.prior
	fit$proportion <- proportion
	fit$s2.post <- eb$s2.post
	fit$t <- eb$t
	fit$df.total <- eb$df.total
	fit$p.value <- eb$p.value
	fit$lods <- eb$lods
	if(!is.null(fit$design) && is.fullrank(fit$design)) {
		F.stat <- limma:::classifyTestsF(fit,fstat.only=TRUE)
		fit$F <- as.vector(F.stat)
		df1 <- attr(F.stat,"df1")
		df2 <- attr(F.stat,"df2")
		if(df2[1] > 1e6){ # Work around bug in R 2.1
			fit$F.p.value <- pchisq(df1*fit$F,df1,lower.tail=FALSE)
		} else{
			fit$F.p.value <- pf(fit$F,df1,df2,lower.tail=FALSE)
		}
	}
	fit
}

ebayesGEO <- function(fit,proportion=0.01,stdev.coef.lim=c(0.1,4),trend=FALSE,robust=FALSE,winsor.tail.p=c(0.05,0.1), proportion_vector=0.01)
#	Empirical Bayes statistics to select differentially expressed genes
#	Gordon Smyth
#	8 Sept 2002.  Last revised 1 May 2013.
{
  #print("using Kevinses")
  coefficients <- fit$coefficients
	stdev.unscaled <- fit$stdev.unscaled
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
	if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
	if(all(!is.finite(sigma))) stop("No finite residual standard deviations")
	if(trend) {
		covariate <- fit$Amean
		if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
	} else {
		covariate <- NULL
	}

#	Moderated t-statistic
	out <- limma:::squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust, winsor.tail.p=winsor.tail.p)
	out$s2.prior <- out$var.prior
	out$s2.post <- out$var.post
	out$var.prior <- out$var.post <- NULL
	out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
	df.total <- df.residual + out$df.prior
	df.pooled <- sum(df.residual,na.rm=TRUE)
	df.total <- pmin(df.total,df.pooled)
	out$df.total <- df.total
	out$p.value <- 2*pt(-abs(out$t),df=df.total)

#	B-statistic
	#out$s2.prior <- prop$probability
	#out$var.prior <- prop$probability
	var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
	out$var.prior <- limma:::tmixture.matrix(out$t,stdev.unscaled,df.total,proportion,var.prior.lim)

	if(any(is.na(out$var.prior))) {
		out$var.prior[ is.na(out$var.prior) ] <- 1/out$s2.prior
		warning("Estimation of var.prior failed - set to default value")
	}
	r <- rep(1,NROW(out$t)) %o% out$var.prior

    #r <- r[1,]

	r <- (stdev.unscaled^2+r) / stdev.unscaled^2
	t2 <- out$t^2
	Infdf <- out$df.prior > 10^6
	if(any(Infdf)) {
		kernel <- t2*(1-1/r)/2
		if(any(!Infdf)) {
			t2.f <- t2[!Infdf]
			r.f <- r[!Infdf]
			df.total.f <- df.total[!Infdf]
			kernel[!Infdf] <- (1+df.total.f)/2*log((t2.f+df.total.f) / (t2.f/r.f+df.total.f))
		}
	} else{
		kernel <- (1+df.total)/2*log((t2+df.total) / (t2/r+df.total))
	}
	out$lods <- log(proportion_vector/(1-proportion_vector))-log(r)/2+kernel

	#out$oldB <- log(proportion/(1-proportion))-log(r)/2+kernel

    #out$logodds <- log(proportion_vector/(1-proportion_vector))

	out
}

# tmixture.matrix <- function(tstat,stdev.unscaled,df,proportion,v0.lim=NULL) {
# #	Estimate the prior variance of the coefficients for DE genes
# #	Gordon Smyth
# #	18 Nov 2002. Last modified 12 Dec 2003.

# 	tstat <- as.matrix(tstat)
# 	stdev.unscaled <- as.matrix(stdev.unscaled)
# 	if(any(dim(tstat) != dim(stdev.unscaled))) stop("Dims of tstat and stdev.unscaled don't match")
# 	if(!is.null(v0.lim)) if(length(v0.lim) != 2) stop("v0.lim must have length 2")
# 	ncoef <- ncol(tstat)
# 	v0 <- rep(0,ncoef)
# 	for (j in 1:ncoef) v0[j] <- tmixture.vector(tstat[,j],stdev.unscaled[,j],df,proportion,v0.lim)
# 	v0
# }

# tmixture.vector <- function(tstat,stdev.unscaled,df,proportion,v0.lim=NULL) {
# #	Estimate scale factor in mixture of two t-distributions
# #	tstat is assumed to follow sqrt(1+v0/v1)*t(df) with probability proportion and t(df) otherwise
# #	v1 is stdev.unscaled^2 and v0 is to be estimated
# #	Gordon Smyth
# #	18 Nov 2002.  Last modified 13 Dec 2003.

# 	if(any(is.na(tstat))) {
# 		o <- !is.na(tstat)
# 		tstat <- tstat[o]
# 		stdev.unscaled <- stdev.unscaled[o]
# 		df <- df[o]
# 	}
# 	ngenes <- length(tstat)
# 	ntarget <- ceiling(proportion/2*ngenes)
# 	if(ntarget < 1) return(NA)

# #	If ntarget is v small, ensure p at least matches selected proportion
# #	This ensures ptarget < 1
# 	p <- max(ntarget/ngenes,proportion)

# 	tstat <- abs(tstat)
# 	if(ngenes>1)
# 		ttarget <- quantile(tstat,(ngenes-ntarget)/(ngenes-1))
# 	else
# 		ttarget <- tstat
# 	top <- (tstat >= ttarget)
# 	tstat <- tstat[top]
# 	v1 <- stdev.unscaled[top]^2
# 	df <- df[top]
# 	r <- ntarget-rank(tstat)+1
# 	p0 <- pt(-tstat,df=df)
# 	ptarget <- ( (r-0.5)/2/ngenes - (1-p)*p0 ) / p
# 	pos <- ptarget > p0
# 	v0 <- rep(0,ntarget)
# 	if(any(pos)) {
# 		qtarget <- qt(ptarget[pos],df=df[pos])
# 		v0[pos] <- v1[pos]*((tstat[pos]/qtarget)^2-1)
# 	}
# 	if(!is.null(v0.lim)) v0 <- pmin(pmax(v0,v0.lim[1]),v0.lim[2])
# 	mean(v0)
# }
