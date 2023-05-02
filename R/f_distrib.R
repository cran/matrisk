#' Distribution
#'
#' @param type_function String argument : "gaussian" for normal distribution or "skew-t" for t-student distribution
#' @param compile_qt Numeric matrix containing different quantiles and associated values
#' @param starting_values Numeric vector with initial values for optimization
#'
#' @return a data.frame with the parameters of the distribution
#' @importFrom dfoptim nmkb
#' @importFrom stats qnorm
#' @importFrom sn qst
#' @description
#' This function is used to estimate the parameters of the distribution (mean and standard deviation for Gaussian, xi, omega, alpha, and nu for skew-t) based on the quantile regression results (Koenker and Basset, 1978). See Adrian et al. (2019) and Adrian et al. (2022) for more details on the estimation steps.
#'
#' @references Adrian, Tobias, Nina Boyarchenko, and Domenico Giannone. "Vulnerable growth." American Economic Review 109.4 (2019): 1263-89.
#' @references Adrian, Tobias, et al. "The term structure of growth-at-risk. " American Economic Journal: Macroeconomics 14.3 (2022): 283-323.
#' @references Koenker, Roger, and Gilbert Bassett Jr. "Regression quantiles." Econometrica: journal of the Econometric Society (1978): 33-50.
#'
#' @export
#'
#' @examples
#' # Import data
#' data("data_euro")
#'
#' # Data process
#' PIB_euro_forward_4 = data_euro["GDP"][c(5:length(data_euro["GDP"][,1])),]
#' FCI_euro_lag_4 = data_euro["FCI"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#' CISS_euro_lag_4 = data_euro["CISS"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#'
#' # for a gaussian
#' quantile_target <- as.vector(c(0.25,0.75))
#' results_quantile_reg <- f_compile_quantile(qt_trgt=quantile_target,
#' v_dep=PIB_euro_forward_4,
#' v_expl=cbind(FCI_euro_lag_4, CISS_euro_lag_4),
#' t_trgt = 30)
#'
#' results_g <- f_distrib(type_function="gaussian",
#' compile_qt=results_quantile_reg,
#' starting_values=c(0, 1))
#'
#' # for a skew-t
#' quantile_target <- as.vector(c(0.10,0.25,0.75,0.90))
#' results_quantile_reg <- f_compile_quantile(qt_trgt=quantile_target,
#' v_dep=PIB_euro_forward_4,
#' v_expl=cbind(FCI_euro_lag_4, CISS_euro_lag_4),
#' t_trgt = 30)
#'
#' results_s <- f_distrib(type_function="skew-t",
#' compile_qt=results_quantile_reg,
#' starting_values=c(0, 1, -0.5, 1.3))
#'
f_distrib <- function(type_function, compile_qt, starting_values){

  # for a gaussian function
  if(type_function=="gaussian"){
    # error management
    if (length(starting_values)!=2 && is.matrix(compile_qt)==FALSE){
      stop("for a gaussian function, 'starting_values' has to be of dimension 2 and 'compile_qt' has to be a matrix")
    }else if(length(starting_values)!=2){
      if(nrow(compile_qt)<2){
        stop("for a gaussian function, 'starting_values' has to be of dimension 2 and 'compile_qt' has to be a matrix with a minimum of 2 rows")
      }else{
        stop("for a gaussian function, 'starting_values' has to be of dimension 2")
      }
    }else if(is.matrix(compile_qt)==FALSE){
      stop("'compile_qt' has to be a matrix")
    }else if(nrow(compile_qt)<2){
      error_results <- TRUE
      stop("'compile_qt' has to be a matrix with a minimum of 2 rows")
    }else{
      # objective function
      f_objective <- function(X, par){
        # initialization
        sum <- 0

        # Loop on each elements of X
        for (compteur in 1:nrow(X)){
          sum <- sum + (qnorm(X[compteur,1], mean=par[1], sd=par[2]) - X[compteur,2])^2
        }
        return(sum)
      }
      # optimization
      param <-nmkb(par=starting_values, fn=f_objective,
                   lower=c(-Inf,0),
                   upper=c(+Inf, +Inf), X=compile_qt)
    }

    results <- data.frame("mean"=param$par[1], "sd"=param$par[2])
    return(results)

  # for a skew-t function
  }else if(type_function=="skew-t"){
    # error management
    if (length(starting_values)!=4 && is.matrix(compile_qt)==FALSE){
      stop("for a skew-t function, 'starting_values' has to be of dimension 4 and 'compile_qt' has to be a matrix")
    }else if(length(starting_values)!=4){
      if(nrow(compile_qt)<4){
        stop("for a skew-t function, 'starting_values' has to be of dimension 4 and 'compile_qt' has to be a matrix with a minimum of 4 rows")
      }else{
        stop("for a skew-t function, 'starting_values' has to be of dimension 4")
      }
    }else if(is.matrix(compile_qt)==FALSE){
      stop("'compile_qt' has to be a matrix")
    }else if(nrow(compile_qt)<4){
      stop("'compile_qt' has to be a matrix with a minimum of 4 rows")
    }else{
      # objective function
      f_objective <- function(X, par){
        # initialization
        sum <- 0
        # Loop on each elements of X
        for (compteur in 1:nrow(X)){
          sum <- sum + (qst(X[compteur,1], xi=par[1], omega=par[2], alpha=par[3], nu=par[4], tol=1e-08, method=0) - X[compteur,2])^2
        }
        return(sum)
      }
      # optimization
      param <-nmkb(par=starting_values, fn=f_objective,
               lower=c(-Inf,10e-6, -1, 10e-6),
               upper=c(+Inf, +Inf, 1, +Inf), X=compile_qt)
      results <- data.frame("xi"=param$par[1], "omega"=param$par[2], "alpha"=param$par[3], "nu"=param$par[4])
      return(results)
    }
  }else{
  }

}
