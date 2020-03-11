#' Estimated Instantaneous Reproduction Numbers by distinguishing
#' imported, infected by imported, and infected by local cases
#'
#' \code{estimate_Rs} estimates the reproduction numbers of an epidemic, given
#' the incidence time series and the serial interval distribution. By default,
#' the incidences are divided into three groups: imported, infected by imported,
#' and infected by local cases. Therefore the function needs three time series
#' and estimates two Instantaneous Reproduction Numbers. In a more general
#' framework, the incidences are divided into five groups: imported from a
#' designated region, imported from other regions, infected by imported from the
#' designated region, infected by imported from other regoins, and infected by
#' local cases. Therefore the function needs five time series and estimates
#' three Instantaneous Reproduction Numbers.
#'
#' @param incid One of the following
#' \itemize{
#'
#' \item{A dataframe of non-negative integers with three columns, so that
#' \code{incid$D_im} contains the incidence of imported cases, \code{incid$I_im}
#' contains the incidence infected by the imported cases, and \code{incid$I_lc}
#' contains the incidence infected by the local cases. If the dataframe
#' contains a column \code{incid$dates}, this is used for plotting.
#' \code{incid$dates} must contains only dates in a row.}
#'
#' \item{A dataframe of non-negative integers with five columns, so that
#' \code{incid$D_im} contains the incidence imported from a designated region,
#' \code{incid$O_im} contains the incidence imported from a other regions,
#' \code{incid$I_imd} contains the incidence infected by the imported from the
#' designated region, \code{incid$I_imo} contains the incidence infected by the
#' imported from other regions, and \code{incid$I_lc} contains the incidence
#' infected by the local cases. If the dataframe contains a column
#' \code{incid$dates}, this is used for plotting. \code{incid$dates} must
#' contains only dates in a row.}
#'
#' \item{An object of class \code{\link{incidence}}}
#'
#' }
#'
#' @param method One of "non_parametric_si", "parametric_si", "uncertain_si",
#'   "si_from_data" or "si_from_sample" (see details).
#'
#' @param si_sample For method "si_from_sample" ; a matrix where each column
#'   gives one distribution of the serial interval to be explored (see details).
#'
#' @param si_data For method "si_from_data" ; the data on dates of symptoms of
#'   pairs of infector/infected individuals to be used to estimate the serial
#'   interval distribution (see details).
#'
#' @param config An object of class \code{estimate_R_config}, as returned by
#' function \code{make_config}.
#'
#' @return {
#' an object of class \code{estimate_Rs}, with components:
#' \itemize{
#'
#' \item{model}{: The model applied to the incidence data, one of "2-region",
#' "3-region"}
#'
#' \item{R_im}{: an object of class \code{estimate_R} (2-region model)}
#'
#' \item{R_imd}{: an object of class \code{estimate_R} (3-region model)}
#'
#' \item{R_imo}{: an object of class \code{estimate_R} (3-region model)}
#'
#' \item{R_lc}{: an object of class \code{estimate_R} (both models)}
#'
#' \item{method}{: the method used to estimate R, one of "non_parametric_si",
#' "parametric_si", "uncertain_si", "si_from_data" or "si_from_sample"}
#'
#' \item{si_distr}{: a vector or dataframe (depending on the method) containing
#'  the discrete serial interval distribution(s) used for estimation}
#'
#' \item{SI.Moments}{: a vector or dataframe (depending on the method)
#' containing the mean and std of the discrete serial interval distribution(s)
#' used for estimation}
#'
#' \item{I}{: the time series of total incidence (sum of all the cases)}
#'
#' \item{D_im}{: the time series of incidence imported from a designated region
#' (in both models)}
#'
#' \item{O_im}{: the time series of incidence imported from a other regions
#' (3-region model)}
#'
#' \item{I_im}{: the time series of incidence infected by the imported cases
#' (2-region model)}
#'
#' \item{I_imd}{: the time series of incidence infected by the imported from the
#' designated region (3-region model)}
#'
#' \item{I_imo}{: the time series of incidence infected by the imported from the
#' other regions (2-region model)}
#'
#' \item{I_lc}{: the time series of incidence infected by the local cases
#' (in both model)}
#'
#' \item{dates}{: a vector of dates corresponding to the incidence time series}
#' }
#' }
#'
#' @details
#' This function extends the original \code{estimate_Rs} by introducing a
#' multi-region model on the incidence data. By default (two-region model), the
#' incidences are divided into three groups: imported, infected by imported,
#' and infected by local cases. Therefore the function needs three time series
#' and estimates two Instantaneous Reproduction Numbers. In a more general
#' framework (three-region model), the incidences are divided into five groups:
#' imported from a designated region, imported from other regions, infected by
#' imported from the designated region, infected by imported from other regoins,
#' and infected by local cases. Therefore the function needs five time series
#' and estimates three Instantaneous Reproduction Numbers.
#'
#' @seealso \code{\link{estimate_R}} \code{\link{discr_si}} \code{\link{make_config}}
#'
#' @author Jinshan Wu \email{jinshanw@bnu.edu.cn}
#'
#' @references {
#' Furthre improved EpiEstim by distinguishing imported,infected by imported,
#' infected by local cases and its application to COVID-19 in China.
#' Cori, A. et al. A new framework and software to estimate time-varying
#' reproduction numbers during epidemics (AJE 2013).
#' Wallinga, J. and P. Teunis. Different epidemic curves for severe acute
#' respiratory syndrome reveal similar impacts of control measures (AJE 2004).
#' Reich, N.G. et al. Estimating incubation period distributions with coarse
#' data (Statis. Med. 2009)
#' }
#'
#' @importFrom incidence incidence
#' @export
#' @examples
#' ##
# SI_MEAN = 8.4
# SI_STD = 3.8
# T <- 50
# R_im <- rep(6, T)
# R_lc <- rep(2, T)
# sim_data <- simulate_Is(R_im, R_lc,
#                                  mean_im = 10,
#                                  mean_si = SI_MEAN, std_si = SI_STD)
# incid <- data.frame(D_im = sim_data$D_im,
#                     I_im = sim_data$I_im,
#                     I_lc = sim_data$I_lc)
# result <- estimate_Rs(incid,
#                       method="parametric_si",
#                       config = make_config(list(
#                         mean_si = SI_MEAN,
#                         std_si = SI_STD)))
# plot(result)
#'
estimate_Rs <- function(incid,
                        method = c(
                          "non_parametric_si",
                          "parametric_si",
                          "uncertain_si",
                          "si_from_data",
                          "si_from_sample"
                        ),
                        si_data = NULL,
                        si_sample = NULL,
                        config = make_config(incid = incid, method = method)) {
  ret <- process_Is(incid)
  incid <- ret$incid
  model <- ret$model
  
  if (model == "2-region") {
    incid_im <- data.frame(
      dates = incid$dates,
      local = incid$I_im,
      imported = incid$D_im
    )
    config$group <- "imported"
    result_R_im <- estimate_R(
      incid = incid_im,
      method = method,
      si_data = si_data,
      si_sample = si_sample,
      config = config
    )
    
    incid_lc <- data.frame(
      dates = incid$dates,
      local = incid$I_lc,
      imported = incid$I_im
    )
    config$group <- "all"
    result_R_lc <- estimate_R(
      incid = incid_lc,
      method = method,
      si_data = si_data,
      si_sample = si_sample,
      config = config
    )
    
    results <- list(R_im = result_R_im)
    results$R_lc = result_R_lc
    results$model = model
    results$method <- result_R_lc$method
    results$si_distr <- result_R_lc$si_distr
    results$SI.Moments <- result_R_lc$SI.Moments
    results$dates <- result_R_lc$dates
    
    results$I <- rowSums(incid[, c("D_im", "I_im", "I_lc")])
    results$D_im <- incid$D_im
    results$I_im <- incid$I_im
    results$I_lc <- incid$I_lc
  } else if (model == "3-region") {
    incid_imd <- data.frame(
      dates = incid$dates,
      local = incid$I_imd,
      imported = incid$D_im
    )
    config$group <- "imported"
    result_R_imd <- estimate_R(
      incid = incid_imd,
      method = method,
      si_data = si_data,
      si_sample = si_sample,
      config = config
    )
    
    incid_imo <- data.frame(
      dates = incid$dates,
      local = incid$I_imo,
      imported = incid$O_im
    )
    config$group <- "imported"
    result_R_imo <- estimate_R(
      incid = incid_imo,
      method = method,
      si_data = si_data,
      si_sample = si_sample,
      config = config
    )
    
    incid_lc <- data.frame(
      dates = incid$dates,
      local = incid$I_lc,
      imported = incid$I_imd + incid$I_imo
    )
    config$group <- "all"
    result_R_lc <- estimate_R(
      incid = incid_lc,
      method = method,
      si_data = si_data,
      si_sample = si_sample,
      config = config
    )
    
    results <- list(R_imd = result_R_imd)
    results$R_imo = result_R_imd
    results$R_lc = result_R_lc
    results$model = model
    results$method <- method
    results$dates <- result_R_imd$dates
    
    results$I <-
      rowSums(incid[, c("D_im", "O_im", "I_imd", "I_imo", "I_lc")])
    results$D_im <- incid$D_im
    results$O_im <- incid$O_im
    results$I_imd <- incid$I_imd
    results$I_imo <- incid$I_imo
    results$I_lc <- incid$I_lc
  } else {
    stop("model must be either '2-region' or '3-region'.")
  }
  
  class(results) <- "estimate_Rs"
  return(results)
}


# check and process incidence data
process_Is <- function(incid) {
  if (inherits(incid, "incidence")) {
    I_inc   <- incid
    incid   <- as.data.frame(I_inc)
    incid$I <- rowSums(incidence::get_counts(I_inc))
  }
  
  if (!is.data.frame(incid)) {
    stop("incid must be a dataframe.")
  }
  
  model = ""
  if (all(c("D_im", "I_im", "I_lc") %in% names(incid))) {
    model = "2-region"
  } else if (all(c("D_im", "O_im", "I_imd", "I_imo", "I_lc") %in% names(incid))) {
    model = "3-region"
  } else {
    stop(
      "incid must be a dataframe with 3 columns called
         'D_im', 'I_im' and 'I_lc', or with 5 columns called
         'D_im', 'O_im', 'I_imd', 'I_imo', and 'I_lc'."
    )
  }
  
  incid[which(is.na(incid))] <- 0
  date_col <- names(incid) == "dates"
  if (any(date_col)) {
    if (any(incid[, !date_col] < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  } else {
    if (any(incid < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  }
  
  if (!is.null(incid$dates)) {
    incid$dates <- check_dates(incid)
  } else {
    incid$dates <- as.numeric(seq_len(nrow(incid)))
  }
  
  return(list(incid = incid, model = model))
}
