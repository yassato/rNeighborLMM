#' Calculating neighbor genotypic identity
#'
#' A function to calculate neighbor genotypic identity, with a given reference scale and a degree of distance decay.
#' @param geno An individual x marker matrix. Bialleles (i.e., A or a) must be converted into -1 or 1 digit.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial points along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial (Euclidean) distance between a focal individual and neighbors to define neighbor effects.
#' @param dist.func An option to set a distance decay function given as a function of scale. Default is NULL, meaning no distance decay.
#' @param grouping A positive integer vector assigning each individual to a group. This argument can be useful when a "smap" contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param n_core No. of cores for a multi-core computation. This does not work for Windows OS. Default is a single-core computation.
#' @return A numeric matrix for neighbor covariates, with no. of individuals x markers.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' The argument \code{dist.func} allows a user-defined function, but default setting is recommended unless spatial distance decay of neighbor effects needs to be modeled.
#' If \code{dist.func} is not NULL, this argument must be pre-defined and given as a function of x. Then output variables are weighted on the basis of Euclidean distance from a focal individual up to \code{scale}.
#' @import parallel
#' @examples
#' set.seed(1)
#' g <- matrix(sample(c(-1,1),100*1000,replace = TRUE),100,1000)
#' gmap <- cbind(c(rep(1,nrow(g)/2),rep(2,nrow(g)/2)),c(1:ncol(g)))
#' x <- rep(c(1:10),10)
#' y <- rep(c(1:10),each=10)
#' smap <- cbind(x,y)
#' grouping <- c(rep(1,nrow(g)/2), rep(2,nrow(g)/2))
#'
#' g_nei <- nei.coval(geno=g,smap=smap,scale=sqrt(2),grouping=grouping)
#'
#' # Distance decay with f(d) = 1/d^2 as suggested by Muir (2005) Genetics 170:1247–1259.
#' dist.func = function(x) return((1/x^2))
#' g_nei <- nei.coval(geno=g,smap=smap,scale=sqrt(2),dist.func=dist.func,grouping=grouping)
#' @export
nei.coval = function(geno, smap, scale, dist.func=NULL, grouping=rep(1,nrow(smap)), n_core=1L) {
  if((min(geno)==-1)&(max(geno)==1)){
    message("calculation starts...")
  } else {
    stop("genotype values must range from -1 to +1")
  }

  p <- nrow(smap)

  g.d2 <- lapply(unique(grouping), function(gi, s2) {
    id <- which(grouping == gi)
    d2 <- outer(smap[id,1], smap[id,1], "-")^2 + outer(smap[id,2], smap[id,2], "-")^2
    d2[d2>s2] <- 0
    if (is.null(dist.func)==TRUE){
      wd <- NULL
    } else if(is.function(dist.func)==TRUE){
      wd <- dist.func(x=sqrt(d2))
    } else {
      stop("dist.func must be a function of x")
    }
    list(id=id,d2=d2,wd=wd)
  }, scale^2)

  n_div <- ceiling(ncol(geno)/3000/n_core)*n_core
  div.i <- div.seq(ncol(geno), max(n_core, n_div))

  coval_i = function(i1, i2) {
    res <- matrix(0, p, i2-i1+1L)
    for (g.d2_i in g.d2) {
      for (i in seq_along(g.d2_i$id)) {
        id <- g.d2_i$id[i]
        j_id <- g.d2_i$id[g.d2_i$d2[,i] != 0]
        if (length(j_id) > 0) {
          if (is.null(dist.func)==TRUE) {
            res_i12 = geno[j_id,i1:i2,drop=FALSE]
          } else {
            w_i <- g.d2_i$wd[g.d2_i$d2[,i] != 0, i]
            res_i12 <- w_i * geno[j_id,i1:i2,drop=FALSE]
          }
          if(length(j_id)==1) {
            res[id,] <- res_i12
          } else {
            res[id,] <- colSums(res_i12, na.rm = TRUE)/colSums(!is.na(res_i12))
          }
        }
      }
    }
    res
  }
  g_nei <- geno*do.call(cbind,
                        parallel::mcmapply(coval_i, div.i$i1, div.i$i2,
                                           mc.cores=getOption("mc.cores", n_core),
                                           USE.NAMES = FALSE, SIMPLIFY = FALSE))
  g_nei[is.na(g_nei)] <- 0

  colnames(g_nei) <- colnames(geno)
  rownames(g_nei) <- rownames(geno)
  message("done.")

  return(g_nei)
}

div.seq <- function(n, n.div) {
  if (n.div >= n) return(list(i1 = 1L:n, i2 = 1L:n))
  c.div <- rep(n%/%n.div, n.div)
  if ((r <- n%%n.div) > 0) {
    c.div[1:r] <- c.div[1:r] + 1L
  }
  i2 <- cumsum(c.div)
  list(i1 = c(1L, i2[-length(i2)]+1L), i2 = i2)
}
