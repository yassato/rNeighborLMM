#' Calculating x-axis for Manhattan plot
#'
#' A utility function to calculate relative genomic positions to depict the x-axis of Manhattan plot
#' @param chr Chromosomal positions. Non-numeric allowed.
#' @param pos Genomic positions in base pairs within chromosomes.
#' @return A list object including relative genomic positions and mid-points of each chromosome.
#'  \itemize{
#'  \item coord relative x-axis positions of Manhattan plot coefficient for self effects.
#'  \item tic mid-point of each chromosome.
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
coord = function(chr, pos) {
  if(length(pos)!=length(chr)) stop("chr and pos length differ")
  chr <- as.factor(chr)
  coord <- 0
  M <- 0
  tic <- numeric(nlevels(chr))
  for (i in 1:nlevels(chr)) {
            w <- (chr == levels(chr)[i])
            pos.c <- pos[w]
            coord[w] <- M + pos.c
            mx <- max(pos.c)
            tic[i] <- M + mx/2
            M <- M + mx
            }
  coord <- coord/M
  tic <- tic/M
  return(list(coord=coord,tic=tic))
}
