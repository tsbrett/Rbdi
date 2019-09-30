#' Likelihood estimation using the Birth-Death-Immigration (BDI) process
#'
#' \code{bdi_likelihood} uses the analytical solution for the
#' BDI process to calculate the log-likelihood of a time series. The BDI model
#' used assumes that the immigration rate and death rate are fixed, and the
#' birth-rate can vary approximately linearly with time (see below).
#'
#' The log-likelihood is calculated for each permutation of four input
#' parameters. These parameters are: i) gamm, the death rate; ii) eta, the
#' immigration rate; iii) R00, the initial value of the reproductive number
#' \eqn{R_0}; iv) dR0, the rate of increase in \eqn{R_0}. The function assumes
#' that the input data x are an equally-spaced time series with time step delta_t.
#'
#'
#' The birth-rate is specified via the reproductive  number \eqn{R_0} (the ratio
#' of birth rate to death rate) as is standard in epidemiological and
#' ecological applications. It is modelled as a piecewise constant function,
#' \deqn{
#' R_0(t) = R00 +  dR0 \sum_{i=1}^{N} \delta \chi_{\delta i \le t},
#' }
#' where \eqn{N} is the number of data points \chi_{\delta i \le t} is an
#' indicator function. For \eqn{dR0 \ll \delta} this approximates a linear
#' variation in \eqn{R_0(t)}. Setting dR0 = 0 recovers the BDI model with
#' constant parameters.
#'
#' As the dynamics may be non-stationary, there is no general expression for the
#' unconditional probability of observing the first data point,
#' \eqn{P(x_0; \theta)}, that do not require additional assumptions on the
#' dynamics. The log-likelihood returned is therefore the conditional log-likelihood,
#' \deqn{\ell(\{x_i\}_{i=0}^N; \theta) = \sum_{i=1}^{N}\log P(x_i|x_{i-1};\theta),}
#' where \eqn{\theta} are the model parameters.
#' @references \url{https://en.wikipedia.org/wiki/Birth-death_process}
#'
#' @param x An equally spaced time series (numeric vector) of the number of
#' individuals present.
#' @param delta_t Time step of the data
#' @param eta Numerical vector listing the set of immigration rates
#' @param gamm Numerical vector listing the set of death (recovery) rates
#' @param R00 Numerical vector listing the set of initial reproductive number
#'  values
#' @param dR0 Numerical vector listing the set of rates of change in
#' reproductive number
#' @return A data frame that contains the log-likelihood for each unique
#' parameter combinations.
#' @export
bdi_likelihood <- function(z, delta_t = 1, eta = 1, gamm = 1,
                                        dR0 = 0.0, R00 = 0.9){

  #library(foreach)
  #library(Rcpp)

  #library("doParallel")
  #doParallel::registerDoParallel(cores=3)
  #RNGkind("L'Ecuyer-CMRG")
  #set.seed(3234)
  #parallel::mc.reset.stream()

  #This needs fixing and might not be the best way of including C functions
  #Rcpp::sourceCpp("./src/bdi.cpp")

  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  levs <- list()
  levs$eta <- eta
  levs$gamm <-  gamm
  levs$dR0 <- dR0
  levs$R00 <- R00
  parameters <- do.call(expand.grid, levs)

  LogLikelihood <- foreach::foreach(i=seq(1, nrow(parameters)),
                           .options.multicore=list(set.seed=TRUE)) %do%
    do.call(C_bdi_ll_linear_R0, append(list(z,delta_t = delta_t), parameters[i,]))

  tm <- parameters
  tm$LogLike <- unlist(LogLikelihood)
  tm$delta <- (1- tm$R00)/(tm$dR0*365)
  return(tm)
}


#' MLE
#'
#' @export
bdi_likelihood_ratio_test <- function(likelihood_data){

  # Calculate number of free parameters
  get_fp <- function(vlist){
    fp=0
    for(v in vlist){
      i = length(unique(ll_test[, v]))
      if(i > 1){
        fp = fp + 1
      }
    }
    return(fp)
   }

  fp_s = get_fp(c("gamm","eta","R00"))
  fp_e = fp_s + get_fp(c("dR0"))



  # could use alternative optimisation methods to find mle
  ll <- na.omit(likelihood_data)
  # Just non-negative (dR0 >= 0)
  ll_e <- ll[ll$dR0>=0, ]
  ll_s <- ll[ll$dR0==0, ]
  if (dim(ll_s)[1] ==0){
    stop("Log-likelihood not calculated for dR0=0")
  }
  ll_max_e <- max(ll$LogLike)
  ll_max_s <- max(ll_s$LogLike)

  mle_e <- ll[ll$LogLike == ll_max_e, ]
  mle_s <- ll_s[ll_s$LogLike == ll_max_s, ]

  AIC_e = 2*fp_e - 2*ll_max_e
  AIC_s = 2*fp_s - 2*ll_max_s
  is_emerging = AIC_e < AIC_s
  cox_delta <- 2*(ll_max_e - ll_max_s)
  sum_stats <- list("ML_e" = ll_max_e, "ML_s" = ll_max_s, "fp_e" = fp_e,
                    "fp_s" = fp_s, "AIC_e" = AIC_e, "AIC_s" = AIC_s,
                    "cox_threshold" = 2*(fp_e-fp_s),
                    "cox_delta" = cox_delta, "emerging" = is_emerging)
  return(sum_stats)
  #Tomorrow: return summary of stats, mle_e and mle_s as two data frames
  #Tomorrow2: add moving window
}



#' Moving window
#'
#' @export
bdi_lrt_moving_window <- function(z, delta_t = 1, eta = 1, gamm = 1,
                                  dR0 = 0.0, R00 = 0.9, window=10){


  get_window_stats <- function(i){
    zt = z[max(i-window,0):(i)]
    bdi_l = bdi_likelihood(zt, delta_t = delta_t, eta = eta, gamm = gamm,
                           dR0 = dR0, R00 = R00)
    bdi_lrt = bdi_likelihood_ratio_test(bdi_l)
    bdi_lrt$i = i
    return(data.frame(bdi_lrt))
  }

  N = length(z)
  df = get_window_stats(1)
  for(i in 2:N){
    df = rbind(df, get_window_stats(i))
  }
  row.names(df) <- df$i

  return(df)
}
