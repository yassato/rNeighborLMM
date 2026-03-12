#' Mixed models for testing self and neighbor effects with their covariance
#'
#' A function to provide coefficients and p-values of self and neighbor effects for each marker.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param g_nei An output of \code{nei_coval()} object, namely an individual x marker matrix including neighbor genotypic identity.
#' @param pheno A numeric vector including phenotypes for individuals
#' @param X0 An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param asym If TRUE, beta12 is also tested. Default is FALSE, which test beta2 without beta12 following the original Neighbor GWAS.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @param ... Arguments to be passed to \code{EM3.cov()}.
#' @return A data.frame including coefficients and p-values of self and neighbor effects, without the chromosome numbers and marker position.
#' \itemize{
#'  \item beta1 coefficient for self effects.
#'  \item beta12 coefficient for asymmetric (one-sided) neighbor effects. Available only when asym = TRUE.
#'  \item beta2 coefficient for symmetric (interactive) neighbor effects.
#'  \item chisq1 likelihood ratio between a null and standard GWAS model.
#'  \item chisq12 likelihood ratio between models with or without asymmetric neighbor effects. Available only when asym = TRUE.
#'  \item chisq2 likelihood ratio between models with or without symmetric neighbor effects.
#'  \item p1 p-value for self effects by a likelihood ratio test between a null and standard GWAS model.
#'  \item p12 p-value for neighbor effects by a likelihood ratio test between models with or without asymmetric neighbor effects. Available only when asym = TRUE.
#'  \item p2 p-value for neighbor effects by a likelihood ratio test models with or without symmetric neighbor effects.
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp}) and Kosuke Hamazaki
#' @import Matrix gaston parallel RAINBOWR
#' @export
neiGWAS.cov = function(geno, g_nei, pheno, X0=matrix(1,nrow=length(pheno)), asym=FALSE, n_core=1L, ...){
  if((min(geno)==-1)&(max(geno)==1)){
    if((min(g_nei)==-1)&(max(g_nei)==1)){
      message("genotype data, ok")
    }
  } else {
    stop("genotype values must range from -1 to +1")
  }

  message("LMM is being solved...")
  args <- list(...)
  args$geno <- geno
  args$g_nei <- g_nei
  args$pheno <- pheno
  args$X0 <- X0
  EM3res <- do.call(neiEM3.cov, args)
  message("LMM, done.")

  message("Calculating weighted kernels...")
  q <- ncol(geno)

  selfCross <- tcrossprod(geno)
  K1 <- (q+selfCross)/(2*(q-1))

  neiCross <- tcrossprod(g_nei)
  K2 <- neiCross/(q-1)

  K12 <- tcrossprod(geno,g_nei)/(q-1)
  K21 <- t(K12)

  NAs <- (is.na(pheno)==FALSE)
  y <- pheno[NAs]
  geno <- geno[NAs,]
  g_nei <- g_nei[NAs,]
  X0 <- as.matrix(X0[NAs,])
  selfK <- K1[NAs,NAs]
  neiK <- K2[NAs,NAs]
  K12 <- K12[NAs,NAs]
  K21 <- K21[NAs,NAs]

  message("GWAS is running...")
  K_prime <- EM3res$weights[1]*selfK +
    EM3res$weights[2]*neiK +
    EM3res$rhosMat[1,2]*sqrt(EM3res$weights[1]*EM3res$weights[2])*(K12 + K21)

  eigenK <- eigen(K_prime)
  lmm0 <- gaston::lmm.diago(Y=y,X=X0,eigenK=eigenK,verbose=FALSE)
  LL0 <- gaston::lmm.diago.profile.likelihood(tau=lmm0$tau,s2=lmm0$sigma2,h2=lmm0$tau/(lmm0$tau+lmm0$sigma2),Y=y,X=X0,eigenK=eigenK)

  if(asym==FALSE){
    test_i = function(i){
      X1 <- cbind(X0,geno[,i])
      lmm1 <- gaston::lmm.diago(Y=y,X=X1,eigenK=eigenK,verbose=FALSE)
      LL1 <- gaston::lmm.diago.profile.likelihood(tau=lmm1$tau,s2=lmm1$sigma2,h2=lmm1$tau/(lmm1$tau+lmm1$sigma2),Y=y,X=X1,eigenK=eigenK)

      X2 <- cbind(X1,g_nei[,i])
      lmm2 <- gaston::lmm.diago(Y=y,X=X2,eigenK=eigenK,verbose=FALSE)
      LL2 <- gaston::lmm.diago.profile.likelihood(tau=lmm2$tau,s2=lmm2$sigma2,h2=lmm2$tau/(lmm2$tau+lmm2$sigma2),Y=y,X=X2,eigenK=eigenK)
      b2 <- lmm2$BLUP_beta[(length(lmm2$BLUP_beta)-1):length(lmm2$BLUP_beta)]

      res_i <- c(b2,2*(LL1[1,1]-LL0[1,1]),2*(LL2[1,1]-LL1[1,1]))
      names(res_i) <- c("beta1","beta2","chisq1","chisq2")

      return(res_i)
    }

    res <- parallel::mcmapply(test_i,1:ncol(geno),mc.cores=n_core)
    res <- t(res)
    res <- as.data.frame(res)

    p1 <- stats::pchisq(res$chisq1,1,lower.tail=FALSE)
    p2 <- stats::pchisq(res$chisq2,1,lower.tail=FALSE)
    res <- data.frame(res,p1,p2)
  } else if(asym==TRUE){
    g_12 <- geno*g_nei

    test_i = function(i){
      X1 <- cbind(X0,geno[,i])
      lmm1 <- lmm.diago(Y=y,X=X1,eigenK=eigenK,verbose=FALSE)
      LL1 <- lmm.diago.profile.likelihood(tau=lmm1$tau,s2=lmm1$sigma2,h2=lmm1$tau/(lmm1$tau+lmm1$sigma2),Y=pheno,X=X1,eigenK=eigenK)

      X12 <- cbind(X1,g_12[,i])
      lmm12 <- lmm.diago(Y=y,X=X12,eigenK=eigenK,verbose=FALSE)
      LL12 <- lmm.diago.profile.likelihood(tau=lmm12$tau,s2=lmm12$sigma2,h2=lmm12$tau/(lmm12$tau+lmm12$sigma2),Y=pheno,X=X12,eigenK=eigenK)

      X2 <- cbind(X12,g_nei[,i])
      lmm2 <- lmm.diago(Y=y,X=X2,eigenK=eigenK,verbose=FALSE)
      LL2 <- lmm.diago.profile.likelihood(tau=lmm2$tau,s2=lmm2$sigma2,h2=lmm2$tau/(lmm2$tau+lmm2$sigma2),Y=pheno,X=X2,eigenK=eigenK)
      b2 <- lmm2$BLUP_beta[(length(lmm2$BLUP_beta)-2):length(lmm2$BLUP_beta)]

      res_i <- c(b2,2*(LL1[1,1]-LL0[1,1]),2*(LL12[1,1]-LL1[1,1]),2*(LL2[1,1]-LL12[1,1]))
      names(res_i) <- c("beta1","beta12","beta2","chisq1","chisq12","chisq2")

      return(res_i)
    }

    res <- parallel::mcmapply(test_i,1:ncol(geno),mc.cores=n_core)
    res <- t(res)
    res <- as.data.frame(res)

    p1 <- stats::pchisq(res$chisq1,1,lower.tail=FALSE)
    p12 <- stats::pchisq(res$chisq12,1,lower.tail=FALSE)
    p2 <- stats::pchisq(res$chisq2,1,lower.tail=FALSE)
    res <- data.frame(res,p1,p12,p2)
  }
  message("GWAS, done.")
  return(res)
}
