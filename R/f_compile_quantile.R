#' Estimation of quantiles
#'
#' @param qt_trgt Numeric vector, dim k, of k quantiles for different qt-estimations
#' @param v_dep Numeric vector of the dependent variable
#' @param v_expl Numeric vector of the (k) explanatory variable(s)
#' @param t_trgt Numeric time target (optional)
#'
#' @return Numeric matrix with the predicted values based on each quantile regression, at time fixed in input
#' @importFrom stats predict
#' @importFrom quantreg rq
#' @export
#' @description
#' Predicted values based on each quantile regression (Koenker and Basset, 1978), at time=t_trgt, for each quantile in qt_trgt.
#'
#' @references Koenker, Roger, and Gilbert Bassett Jr. "Regression quantiles." Econometrica: journal of the Econometric Society (1978): 33-50.
#'
#' @examples
#' # Import data
#' data("data_euro")
#'
#' #' # Data process
#' PIB_euro_forward_4 = data_euro["GDP"][c(5:length(data_euro["GDP"][,1])),]
#' FCI_euro_lag_4 = data_euro["FCI"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#' CISS_euro_lag_4 = data_euro["CISS"][c(1:(length(data_euro["GDP"][,1]) - 4)),]
#'
#' quantile_target <- as.vector(c(0.10,0.25,0.75,0.90))
#' results_quantile_reg <- f_compile_quantile(qt_trgt=quantile_target,
#' v_dep=PIB_euro_forward_4,
#' v_expl=cbind(FCI_euro_lag_4, CISS_euro_lag_4),
#' t_trgt = 30)
#'
f_compile_quantile <- function(qt_trgt, v_dep, v_expl, t_trgt){

  # number of quantile regressions (for k quantile regressions)
  nb_qt <- length(qt_trgt)

  # initialization of matrix results
  results_qt <- matrix(data=0, ncol=2, nrow=nb_qt)

  # loop on each quantile regression
  for (ct_qt in 1:nb_qt){

    reg_qt <- rq(v_dep ~ cbind(v_expl), tau=qt_trgt[ct_qt]) # quantile regression
    pred_qt <- predict(reg_qt, newdata=as.data.frame(v_expl)) # prediction

    # store the value that corresponds to t_trgt, time target
    if(missing(t_trgt)){
      results_qt[ct_qt,2] <- pred_qt[length(pred_qt)]
    }else{
      results_qt[ct_qt,2] <- pred_qt[t_trgt]
    }

  }
  results_qt[,1] <- qt_trgt
  return(results_qt)
}
