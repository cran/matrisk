#' Historical distributions
#'
#' @param qt_trgt Numeric vector, dim k, of k quantiles for different qt-estimations
#' @param v_dep Numeric vector of the dependent variable
#' @param v_expl Numeric vector of the (k) explanatory variable(s)
#' @param type_function String argument : "gaussian" for normal distribution or "skew-t" for t-student distribution
#' @param starting_values Numeric vector with initial values for optimization
#' @param step Numeric argument for accuracy graphics abscissa
#' @param x_min Numeric optional argument (default value = -15)
#' @param x_max Numeric optional argument (default value = 10)
#'
#' @return
#' A list with:
#' \item{distrib_histo}{Numeric matrix with historical values of x, y and t}
#' \item{param_histo}{Numeric matrix containing the parameters of the distribution for each period}
#'
#' @importFrom stats dnorm
#' @importFrom sn dst
#' @import plot3D
#' @export
#' @description
#' This function is based on f_distrib function (Adrian et al., 2019; Adrian et al., 2022) and is used to get historical estimation of empirical distributions and associated parameters. Results allow to realize a 3D graphical representation.
#'
#' @references Adrian, Tobias, Nina Boyarchenko, and Domenico Giannone. "Vulnerable growth." American Economic Review 109.4 (2019): 1263-89.
#' @references Adrian, Tobias, et al. "The term structure of growth-at-risk." American Economic Journal: Macroeconomics 14.3 (2022): 283-323.
#'
#' @examples
#' \donttest{# Import data
#' data("data_euro")
#'
#' # Data process
#' PIB_euro_forward_4 = data_euro["GDP"][c(5:length(data_euro["GDP"][,1])),]
#' FCI_euro_lag_4 = data_euro["FCI"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#' CISS_euro_lag_4 = data_euro["CISS"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#'
#' results_histo <- f_distrib_histo(qt_trgt=c(0.10,0.25,0.75,0.90), v_dep=PIB_euro_forward_4,
#' v_expl=cbind(FCI_euro_lag_4,CISS_euro_lag_4),
#' type_function="skew-t",
#' starting_values=c(0, 1, -0.5, 1.3),
#' step=5, x_min=-10, x_max=5)
#'
#' library(plot3D) # load
#' scatter3D(results_histo$distrib_histo[,3],
#' results_histo$distrib_histo[,1],
#' results_histo$distrib_histo[,2],
#' pch = 10,  theta = 70, phi = 10,
#' main = "Distribution of GDP Growth over time - Euro Area",
#' xlab = "Date",
#' ylab ="Pib",
#' zlab="", cex = 0.3)}

f_distrib_histo <- function(qt_trgt, v_dep, v_expl, type_function, starting_values, step, x_min, x_max){

  # error management
  if (missing(x_min)){
    x_min <- -15
  }else{
  }
  if (missing(x_max)){
    x_max <- 10
  }else{
  }
  # initialization abscissa for distrib
  x <- as.vector(seq(from=x_min, to =x_max, by=step))
  lg_x <- length(x)
  if (is.vector(v_dep)==FALSE){
    stop("'v_dep' has to be a vector")
  }else{
    nb_T <- length(v_dep)
  }

  # initialization matrix results for historical distributions
  res_histo <- matrix(data=0,ncol=3, nrow=lg_x*nb_T)
  if(type_function=="gaussian"){
    param_histo <- matrix(data=0,ncol=3, nrow=nb_T)
  }else if(type_function=="skew-t"){
    param_histo <- matrix(data=0,ncol=5, nrow=nb_T)
  }else{
  }
  # estimation distrib for each T
  for (ct_period in 1:nb_T){
    res_qt_reg <- f_compile_quantile(qt_trgt, v_dep, v_expl, t_trgt = ct_period)

    results_distrib <- f_distrib(type_function, compile_qt=res_qt_reg, starting_values)

    if(type_function=="gaussian"){
      y <- dnorm(x, mean=results_distrib$mean, sd=results_distrib$sd)
    }else if(type_function=="skew-t"){
      y <- dst(x, xi=results_distrib$xi, omega=results_distrib$omega,
               alpha=results_distrib$alpha, nu=results_distrib$nu, dp=NULL, log=FALSE)
    }else{
    }

    res_histo[(1+(ct_period-1)*lg_x):(ct_period*lg_x),1] <- x
    res_histo[(1+(ct_period-1)*lg_x):(ct_period*lg_x),2] <- y
    res_histo[(1+(ct_period-1)*lg_x):(ct_period*lg_x),3] <- ct_period

    if(type_function=="gaussian"){
      param_histo[ct_period,1] <- ct_period
      param_histo[ct_period,2] <- results_distrib$mean
      param_histo[ct_period,3] <- results_distrib$sd
    }else if(type_function=="skew-t"){
      param_histo[ct_period,1] <- ct_period
      param_histo[ct_period,2] <- results_distrib$xi
      param_histo[ct_period,3] <- results_distrib$omega
      param_histo[ct_period,4] <- results_distrib$alpha
      param_histo[ct_period,5] <- results_distrib$nu
    }else{
    }
  }
  results <- list("distrib_histo"=res_histo, "param_histo"=param_histo)
  return(results)
}
