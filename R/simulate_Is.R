#' Simulate Incidence Time Series of an Epidemic based on two-refion model
#'
#' \code{simulate_Is} simulate incidence time series of an epidemic, given
#' the time series of \eqn{D_{Im}}, \eqn{R_{Im}} and \eqn{R_{Lc}}, and the
#' serial interval distribution.
#'
#' @param R_im A vector of non-negative integers containing the time series of
#' the mean Instantaneous Reproduction Number for imported infectors.
#'
#' @param R_lc A vector of non-negative integers containing the time series of
#' the mean Instantaneous Reproduction Number for local infectors.
#'
#' @param D_im A vector of non-negative integers containing the time series of
#' the imported incidences.
#'
#' @param mean_im The time series of the imported incidences will be drawn from
#' a Poisson distribution with mean \code{mean_im} without given \code{D_im}
#'
#' @param si_distr A vector of probabilities giving the discrete distribution of
#' the serial interval, starting with si_distr[1] (probability that the serial
#' interval is zero), which should be zero. If not given, \code{si_distr} will
#' be generated using \code{\link{discr_si}} with parameters \code{mean_si} and
#' \code{std_si}.
#'
#' @param mean_si The mean serial interval.
#'
#' @param std_si The stadard deviation of the serial interval.
#' 
#' @param I_past The historical incidence data, if not NULL, must be a dataframe 
#' with 3 columns called 'D_im', 'I_im' and 'I_lc'.
#'
#' @return A list of vectors
#'
#' @seealso \code{\link{estimate_Rs}} \code{\link{discr_si}}
#'
#' @author Jinshan Wu \email{jinshanw@bnu.edu.cn}
#'
#' @export
#'
simulate_Is <- function (R_im,
                         R_lc,
                         D_im = NULL,
                         mean_im = NULL,
                         si_distr = NULL,
                         mean_si = NULL,
                         std_si = NULL,
                         I_past = NULL,
                         window = 6) {
  if (!is.vector(R_im)) {
    stop("R_im must be a vector.")
  }
  if (!is.vector(R_lc)) {
    stop("R_lc must be a vector.")
  }
  if (length(R_im) != length(R_lc)) {
    stop("R_im and R_lc must have the same length.")
  }
  if (!is.null(D_im)) {
    if (!is.vector(D_im) || length(D_im) != length(R_im)) {
      stop("D_im must be either NULL or a vector with the same length as R_im.")
    }
  } else if (!is.null(mean_im)) {
    D_im <- rpois(length(R_im), mean_im)
  } else {
    stop("D_im and mean_im cannot both be null.")
  }
  
  T <- length(D_im)
  T_all <- T
  if (!is.null(I_past)) {
    T_all <- T + nrow(I_past)
  }
  
  if (!is.null(si_distr)) {
    if (!is.vector(si_distr)) {
      stop("si_distr must be either NULL or a vector.")
    }
  } else if (!is.null(mean_si) && !is.null(std_si)) {
    si_distr <- discr_si(seq(0, T_all - 1), mean_si, std_si)
  } else {
    stop("si_distr and (mean_si, std_si) cannot both be null.")
  }
  if (length(si_distr) < T_all + 1) {
    si_distr[seq(length(si_distr) + 1, T_all + 1)] <- 0
  }
  
  if (!is.null(I_past)) {
    T_past <- nrow(I_past)
    lambda_im <- vector()
    lambda_im[1] <- NA
    lambda_lc <- vector()
    lambda_lc[1] <- NA
    
    for (t in seq(2, T_past)) {
      lambda_im[t] <-
        sum(si_distr[seq_len(t)] * I_past$D_im[seq(t, 1)], na.rm = TRUE)
      lambda_lc[t] <-
        sum(si_distr[seq_len(t)] * (I_past$I_im + I_past$I_lc)[seq(t, 1)], na.rm = TRUE)
    }
    
    I_im <- rep(0, T)
    I_lc <- rep(0, T)
    for (t in seq(1, T)) {
      t1 <- t + T_past
      I <-
        rbind(I_past, data.frame(
          D_im = D_im,
          I_im = I_im,
          I_lc = I_lc
        ))
      lambda_im[t1] <-
        sum(si_distr[seq_len(t1)] * I$D_im[seq(t1, 1)], na.rm = TRUE)
      lambda_lc[t1] <-
        sum(si_distr[seq_len(t1)] * (I$I_im + I$I_lc)[seq(t1, 1)], na.rm = TRUE)
      
      if (window <= 0) {
        I_im[t] <- rpois(1, R_im[t] * lambda_im[t1])
        I_lc[t] <- rpois(1, R_lc[t] * lambda_lc[t1])
      } else {
        t0 <- max(1, t1 - window)
        mean_I_im <-
          R_im[t] * sum(lambda_im[seq(t0, t1)], na.rm = TRUE) - sum(I$I_im[seq(t0, t1 - 1)], na.rm = TRUE)
        if (mean_I_im < 0) {
          mean_I_im <- 0
        }
        mean_I_lc <-
          R_lc[t] * sum(lambda_lc[seq(t0, t1)], na.rm = TRUE) - sum(I$I_lc[seq(t0, t1 - 1)], na.rm = TRUE)
        if (mean_I_lc < 0) {
          mean_I_lc <- 0
        }
        I_im[t] <- rpois(1, mean_I_im)
        I_lc[t] <- rpois(1, mean_I_lc)
      }
    }
  } else {
    lambda_im <- vector()
    lambda_im[1] <- NA
    lambda_lc <- vector()
    lambda_lc[1] <- NA
    
    I_im <- rep(0, T)
    I_lc <- rep(0, T)
    for (t in seq(2, T)) {
      lambda_im[t] <-
        sum(si_distr[seq_len(t)] * D_im[seq(t, 1)], na.rm = TRUE)
      lambda_lc[t] <-
        sum(si_distr[seq_len(t)] * (I_im + I_lc)[seq(t, 1)],
            na.rm = TRUE)
      
      if (window <= 0) {
        I_im[t] <- rpois(1, R_im[t] * lambda_im[t])
        I_lc[t] <- rpois(1, R_lc[t] * lambda_lc[t])
      } else {
        t0 <- max(1, t - window)
        mean_I_im <-
          R_im[t] * sum(lambda_im[seq(t0, t)], na.rm = TRUE) - sum(I$I_im[seq(t0, t - 1)], na.rm = TRUE)
        if (mean_I_im < 0) {
          mean_I_im <- 0
        }
        mean_I_lc <-
          R_lc[t] * sum(lambda_lc[seq(t0, t)], na.rm = TRUE) - sum(I$I_lc[seq(t0, t - 1)], na.rm = TRUE)
        if (mean_I_lc < 0) {
          mean_I_lc <- 0
        }
        I_im[t] <- rpois(1, mean_I_im)
        I_lc[t] <- rpois(1, mean_I_lc)
      }
    }
  }
  
  result = list(
    R_im = R_im,
    R_lc = R_lc,
    D_im = D_im,
    I_im = I_im,
    I_lc = I_lc,
    I = D_im + I_im + I_lc,
    si_distr = si_distr,
    lambda_im = lambda_im,
    lambda_lc = lambda_lc
  )
  return(result)
}
