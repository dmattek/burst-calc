
#' Calculate burst size function
#' 
#' Calculates burst size function \code{Phi(tau)} according to Eq.4 of
#' Dobrzynski & Bruggeman (2009) Elongation dynamics shape bursty transcription and translation. PNAS, 106(8)
#' Formula: \code{Phi (tau) = total no. of intervals / number of intervals longer than tau}
#' 
#' The input should contain a single-column data with intervals between events such as 
#' firings of a neuron, or production of a protein. The function searches for the number of events
#' longer than a threshold. The procedure is repeated for a series of threshold lengths.
#' The output data frame contains two columns with the number of events for each of the thresholds.
#'
#' @param in.v Vector with intervals
#' @param in.nthr Total number of thresholds
#' @param in.sc Scaling factor for threshold interval. Provide as 1/x to obtain the maximum threshold size = x.
#'
#' @return Data frame with threshold lengths and the value of burst size function.
#' @export
#'
#' @examples
#' df = read.csv('example-data-1/output1.prel', header = F)
#' phi=calcBurstSzFn(v$V1, 100, 1/2)
#' plot(phi$t, phi$b, type = 'l')
#' 
calcBurstSzFn = function (in.v, in.nthr = 100, in.sc = 1/2) {
  if(!is.vector(in.v))
    stop('Input in.v must be a vector.')
  
  loc.vlen = length(in.v)
  loc.maxthr = 1. / (in.sc * in.nthr)
  
  # allocate data
  loc.out = data.frame(t = rep(0, in.nthr),
                       b = rep(0, in.nthr))
  
  # Search for the number of events longer than a threshold
  # Loop over a range of thresholds
  for(ii in 1:in.nthr) {
    loc.sublen = length(in.v[in.v > ii * loc.maxthr])
    loc.out$t[ii] = ii / (in.sc * in.nthr)
    loc.out$b[ii] = loc.vlen / loc.sublen
  }
  
  return(loc.out)
}
