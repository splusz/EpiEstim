#' Simulate Incidence Time Series of an Epidemic
#'
#' \code{simulate_Is} simulate incidence time series of an epidemic, given
#' the time series of \code{R} and \code{I_imported}, and the
#' serial interval distribution.
#'
#' @param R A vector of non-negative integers containing the time series of
#' the mean Instantaneous Reproduction Number.
#'
#' @param I_imported A vector of non-negative integers containing the time series of
#' the imported incidences.
#'
#' @param mean_im The time series of the imported incidences will be drawn from
#' a Poisson distribution with mean \code{mean_im} without given \code{I_imported}
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
#' @param I_past The historical incidence data, if not NULL, must be a dataframe with 2 
#' columns called 'local' and 'imported'.
#'
#' @return A list of vectors
#'
#' @seealso \code{\link{estimate_R}} \code{\link{discr_si}}
#'
#' @author Jinshan Wu \email{jinshanw@bnu.edu.cn}
#'
#' @export
#'
simulate_I <- function (R,
                        I_imported = NULL,
                        mean_im = NULL,
                        si_distr = NULL,
                        mean_si = NULL,
                        std_si = NULL,
                        I_past = NULL) {
  if (!is.vector(R)) {
    stop("R must be a vector.")
  }
  if (!is.null(I_imported)) {
    if (!is.vector(I_imported) || length(I_imported) != length(R)) {
      stop("I_imported must be either NULL or a vector with the same length as R.")
    }
  } else if (!is.null(mean_im)) {
    I_imported <- rpois(length(R), mean_im)
  } else {
    stop("I_imported and mean_im cannot both be null.")
  }
  
  T <- length(I_imported)
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
    lambda <- vector()
    lambda[1] <- NA
    for (t in seq(2, T_past)) {
      lambda[t] <- sum(si_distr[seq_len(t)] * rowSums(I_past[seq(t, 1), c("local", "imported")]),
                      na.rm = TRUE)
    }

    I_local <- rep(0, T)
    for (t in seq(1, T)) {
      t1 <- t + T_past
      I <- rbind(I_past, data.frame(imported = I_imported, local = I_local))
      lambda[t1] <- sum(si_distr[seq_len(t1)] * rowSums(I[seq(t1, 1), c("local", "imported")]),
                      na.rm = TRUE)
      I_local[t] = rpois(1, R[t] * lambda[t1])
    }
  } else {
    lambda <- vector()
    lambda[1] <- NA

    I_local <- rep(0, T)
    for (t in seq(2, T)) {
      lambda[t] <- sum(si_distr[seq_len(t)] * (I_imported + I_local)[seq(t, 1)],
                      na.rm = TRUE)
      I_local[t] = rpois(1, R[t] * lambda[t])
    }
  }

  result = list(
    R = R,
    I = I_imported + I_local,
    I_imported = I_imported,
    I_local = I_local,
    si_distr = si_distr,
    lambda = lambda
  )
  return(result)
}
