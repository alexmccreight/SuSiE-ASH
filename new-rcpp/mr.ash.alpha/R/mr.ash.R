#' @title Multiple Regression with Adaptive Shrinkage
#'
#' @description Model fitting algorithms for Multiple Regression with
#'   Adaptive Shrinkage ("Mr.ASH"). Mr.ASH is a variational empirical
#'   Bayes (VEB) method for multiple linear regression. The fitting
#'  algorithms (locally) maximize the approximate marginal likelihood
#'   (the "evidence lower bound", or ELBO) via coordinate-wise updates.
#'
#' @details Mr.ASH is a statistical inference method for the following
#' multiple linear regression model: \deqn{y | X, \beta, \sigma^2 ~
#' N(X \beta, \sigma I_n),} in which the regression coefficients
#' \eqn{\beta} admit a mixture-of-normals prior, \deqn{\beta | \pi,
#' \sigma ~ g = \sum_{k=1}^K N(0, \sigma^2 \sigma_k^2).} Each mixture
#' component in the prior, \eqn{g}, is a normal density centered at
#' zero, with variance \eqn{\sigma^2 \sigma_k^2}.
#'
#' The fitting algorithm, it run for a large enough number of
#' iterations, will find an approximate posterior for the regression
#' coefficients, denoted by \eqn{q(\beta)}, residual variance
#' parameter \eqn{sigma^2}, and prior mixture weights \eqn{\pi_1,
#' \ldots, \pi_K} maximizing the evidence lower bound, \deqn{F(q, \pi,
#' \sigma^2) = E_q \log p(y | X, \beta, \sigma^2) - \sum_{j=1}^p
#' D_{KL}(q_j || g),} where \eqn{D_{KL}(q_j || g)} denotes the
#' Kullback-Leibler (KL) divergence, a measure of the "distance"
#' between (approximate) posterior \eqn{q_j(\beta_j)} and prior
#' \eqn{g(\beta_j)}. The fitting algorithm iteratively updates the
#' approximate posteriors \eqn{q_1, \ldots, q_p}, separately for each
#' \eqn{j = 1, \ldots, p} (in an order determined by
#' \code{update.order}), then separately updates the mixture weights
#' and \eqn{\pi} and residual variance \eqn{\sigma^2}. This
#' coordinate-wise update scheme iterates until the convergence
#' criterion is met, or until the algorithm hits an upper bound on
#' the number of iterations (specified by \code{max.iter}).
#'
#' See \sQuote{References} for more details about the model and
#' algorithm.
#'
#' @param X The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. The matrix
#'   cannot be sparse.
#'
#' @param y The observed continuously-valued responses, a vector of
#'   length p.
#'
#' @param Z The covariate matrix, of dimension (n,k), where k is the
#'   number of covariates. This feature is not yet implemented;
#'   \code{Z} must be set to \code{NULL}.
#'
#' @param sa2 The vector of prior mixture component variances. The
#'   variances should be in increasing order, starting at zero; that is,
#'   \code{sort(sa2)} should be the same as \code{sa2}. When \code{sa2}
#'   is \code{NULL}, the default setting is used, \code{sa2[k] =
#'   (2^(0.05*(k-1)) - 1)^2}, for \code{k = 1:20}. For this default
#'   setting, \code{sa2[1] = 0}, and \code{sa2[20]} is roughly 1.
#'
#' @param method_q The algorithm used to update the variational
#'   approximation to the posterior distribution of the regression
#'   coefficients, \code{method = "sigma_dep_q"}, \code{method =
#'   "sigma_indep_q"} and \code{"sigma_scaled_beta"}, take different
#'   approaches to updating the residual variance \eqn{sigma^2}.
#'
#' @param method_g \code{method = "caisa"}, an abbreviation of
#'   "Cooridinate Ascent Iterative Shinkage Algorithm", fits the model
#'   by approximate EM; it iteratively updates the variational
#'   approximation to the posterior distribution of the regression
#'   coefficients (the approximate E-step) and the model parameters
#'   (mixture weights and residual covariance) in an approximate
#'   M-step. Settings \code{method = "block"} and
#'   \code{method = "accelerate"} are considered experimental.
#'
#' @param max.iter The maximum number of outer loop iterations allowed.
#'
#' @param min.iter The minimum number of outer loop iterations allowed.
#'
#' @param beta.init The initial estimate of the (approximate)
#'   posterior mean regression coefficients. This should be \code{NULL},
#'   or a vector of length p. When \code{beta.init} is \code{NULL}, the
#'   posterior mean coefficients are all initially set to zero.
#'
#' @param update.pi If \code{update.pi = TRUE}, the mixture
#'   proportions in the mixture-of-normals prior are estimated from the
#'   data. In the manuscript, \code{update.pi = TRUE}.
#'
#' @param pi The initial estimate of the mixture proportions
#'   \eqn{\pi_1, \ldots, \pi_K}. If \code{pi} is \code{NULL}, the
#'   mixture weights are initialized to \code{rep(1/K,K)}}, where
#'   \code{K = length(sa2).
#'
#' @param update.sigma2 If \code{update.sigma2 = TRUE}, the residual
#'   variance \eqn{sigma^2} is estimated from the data.  In the manuscript,
#'   \code{update.sigma = TRUE}.
#'
#' @param sigma2 The initial estimate of the residual variance,
#'   \eqn{\sigma^2}. If \code{sigma2 = NULL}, the residual variance is
#'   initialized to the empirical variance of the residuals based on the
#'   initial estimates of the regression coefficients, \code{beta.init},
#'   after removing linear effects of the intercept and any covariances.
#'
#' @param update.order The order in which the co-ordinate ascent
#'   updates for estimating the posterior mean coefficients are
#'   performed. \code{update.order} can be \code{NULL}, \code{"random"},
#'   or any permutation of \eqn{(1,\ldots,p)}, where \code{p} is the number
#'   of columns in the input matrix \code{X}. When \code{update.order}
#'   is \code{NULL}, the co-ordinate ascent updates are performed in
#'   order in which they appear in \code{X}; this is equivalent to
#'   setting \code{update.order = 1:p}. When \code{update.order =
#'   "random"}, the co-ordinate ascent updates are performed in a
#'   randomly generated order, and this random ordering is different at
#'   each outer-loop iteration.
#'
#' @param standardize The logical flag for standardization of the
#'   columns of X variable, prior to the model fitting. The coefficients
#'   are always returned on the original scale.
#'
#' @param intercept When \code{intercept = TRUE}, an intercept is
#'   included in the regression model.
#'
#' @param tol Additional settings controlling behaviour of the model
#'   fitting algorithm. \code{tol$convtol} controls the termination
#'   criterion for the model fitting. When \code{update.pi = TRUE}, the
#'   outer-loop updates stop when the largest change in the mixture
#'   weights is less than \code{convtol*K}; when \code{update.pi =
#'   FALSE}, the outer-loop updates stop when the largest change in the
#'   estimates of the posterior mean coefficients is less than
#'   \code{convtol*K}. \code{tol$epstol} is a small, positive number
#'   added to the likelihoods to avoid logarithms of zero.
#'
#' @return A list object with the following elements:
#'
#' \item{intercept}{The estimated intercept.}
#'
#' \item{beta}{A vector containing posterior mean estimates of the
#'   regression coefficients for all predictors.}
#'
#' \item{sigma2}{The estimated residual variance.}
#'
#' \item{pi}{A vector of containing the estimated mixture
#'   proportions.}
#'
#' \item{iter}{The number of outer-loop iterations that were
#'   performed.}
#'
#' \item{update.order}{The ordering used for performing the
#'   coordinate-wise updates. For \code{update.order = "random"}, the
#'   orderings for outer-loop iterations are provided in a vector of
#'   length \code{p*max.iter}, where \code{p} is the number of predictors.}
#'
#' \item{varobj}{A vector of length \code{numiter}, containing the
#'   value of the variational objective (equal to the negative "evidence
#'   lower bound") attained at each (outer-loop) model fitting
#'   iteration. Note that the objective does not account for the
#'   intercept term, even when \code{intercept = TRUE}; therefore, this
#'   value shoudl be interpreted as being an approximation to the
#'   marginal likelihood \emph{conditional} on the estimate of the
#'   intercept.}
#'
#' \item{data}{The preprocessed data (X, Z, y) provided as input to the model
#'   fitting algorithm. \code{data$w} is equal to
#'   \code{diag(crossprod(X))}, in which \code{X} is the preprocessed
#'   data matrix. Additionally, \code{data$sa2} gives the prior variances
#'   used.}
#'
#' @seealso \code{\link{get.full.posterior}}, \code{\link{predict.mr.ash}}
#'
#' @references
#'
#' Y. Kim (2020), Bayesian shrinkage methods for high dimensional
#' regression. Ph.D. thesis, University of Chicago.
#'
#' @useDynLib mr.ash.alpha
#'
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' @importFrom stats var
#'
#' @examples
#' ### generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#'
#' ### fit Mr.ASH
#' fit.mr.ash  = mr.ash(X,y, method_q = "sigma_indep_q")
#' fit.mr.ash  = mr.ash(X,y, method_q = "sigma_scaled_beta")
#' fit.mr.ash  = mr.ash(X,y, method_q = "sigma_dep_q")
#'
#' ### prediction routine
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ynew        = Xnew %*% beta + rnorm(n)
#' ypred       = predict(fit.mr.ash, Xnew)
#'
#' ### test error
#' rmse        = norm(ynew - ypred, '2') / sqrt(n)
#'
#' ### coefficients
#' betahat     = predict(fit.mr.ash, type = "coefficients")
#' # this equals c(fit.mr.ash$intercept, fit.mr.ash$beta)
#'
#' @export
#'
mr.ash                      = function(X, y, Z = NULL, sa2 = NULL,
                                       method_q = c("sigma_dep_q","sigma_indep_q"),
                                       method_g = c("caisa","accelerate","block"),
                                       max.iter = 1000, min.iter = 1,
                                       beta.init = NULL,
                                       update.pi = TRUE, pi = NULL,
                                       update.sigma2 = TRUE, sigma2 = NULL,
                                       update.order = NULL,
                                       standardize = FALSE, intercept = TRUE,
                                       tol = set_default_tolerance(),
                                       diagXtOmegaX =  diagXtOmegaX, #CHANGED!
                                       XtOmega = XtOmega, #CHANGED!
                                       V = V, # CHANGED!
                                       tausq = tausq,
                                       sum_Dsq = sum_Dsq,
                                       Dsq = Dsq,
                                       VtXt = VtXt
                                       ){

  # get sizes
  n                 = nrow(X)
  p                 = ncol(X)

  # check necessary conditions
  if (!is.null(sa2)) {
    if (any(sa2 < 0)) {
      stop ("all the mixture component variances must be non-negative.")
    }
    if (sa2[1] != 0) {
      stop ("the first mixture component variance sa2[1] must be 0.")
    }
  }

  # check Z
  if (!is.null(Z)) {
    stop("covariates are not currently fully implemented; Z should be set to NULL")
  }

  # match method
  method_q          = match.arg(method_q)
  method_g          = match.arg(method_g)

  # set default tolerances unless specified
  tol0              = set_default_tolerance()
  tol               = modifyList(tol0,tol,keep.null = TRUE)

  # remove covariates
  data              = remove_covariate(X, y, Z, standardize, intercept)

  # initialize beta
  if ( is.null(beta.init) ){
    data$beta       = as.vector(double(p))
  } else {
    if (standardize) {
      data$beta     = as.vector(beta.init) * attr(data$X,"scaled:scale")
    } else {
      data$beta     = as.vector(beta.init)
    }
  }
  data$beta[1]      = data$beta[1] + 0   # to make sure beta.init is not modified

  # initialize r
  r                 = data$y - data$X %*% data$beta

  # sigma2
  if (is.null(sigma2))
    sigma2 = c(var.n(r))

  # precompute x_j^T x_j
  # w = colSums(data$X^2) #CHANGED!
  w = 1/diagXtOmegaX
  data$w  = w

  # set sa2 if missing
  if ( is.null(sa2) ) {
    sa2             = (2^((0:19) / 20) - 1)^2
    sa2             = sa2 / median(data$w) * n
  }
  K                 = length(sa2)
  data$sa2          = sa2

  # initialize other parameters
  if ( is.null(pi) ) {
    if ( is.null(beta.init) ){

      Phi           = matrix(1,p,K)/K
      pi            = rep(1,K)/K

    } else {

      S             = outer(1/w, sa2, '+') * sigma2
      Phi           = -data$beta^2/S/2 - log(S)/2
      Phi           = exp(Phi - apply(Phi,1,max))
      Phi           = Phi / rowSums(Phi)
      pi            = colMeans(Phi)

    }
  } else
    Phi             = matrix(rep(pi, each = p), nrow = p)
  pi[1]            <- pi[1] + 0

  # verbose = TRUE (TO DO LIST)
  verbose           = TRUE

  # run algorithm

  if ( is.null(update.order) ) {
    o               = rep(0:(p-1), max.iter)
  } else if (is.numeric(update.order)) {
    o               = rep(update.order - 1, max.iter)
  } else if (update.order == "random") {
    o               = random_order(p, max.iter)
  }

  out               = caisa_rcpp (data$X, data$y, w, sa2, pi, data$beta, r, sigma2, o,
                                  max.iter, min.iter, tol$convtol, tol$epstol,
                                  method_q, update.pi, update.sigma2, verbose,
                                  XtOmega,
                                  tausq,
                                  sum_Dsq,
                                  V,
                                  Dsq,
                                  VtXt
                                  ) # CHANGED! added: XtOmega (pxn matrix)

  if (method_q == "sigma_scaled_beta") {
    out$beta        = out$beta * sqrt(out$sigma2)
  }

  ## polish return object
  out$intercept     = c(data$ZtZiZy - data$ZtZiZX %*% out$beta)
  data["beta"]      = NULL
  out$data          = data
  out$update.order  = o

  ## rescale beta is needed
  if (standardize)
    out$beta        = out$beta / attr(data$X, "scaled:scale")
  class(out)       <- c("mr.ash", "list")

  ## warn if necessary
  if (update.pi & out$pi[K] > 1/K) {
    warning(sprintf(paste("The mixture proportion associated with the",
                          "largest prior variance is greater than %0.2e;",
                          "this indicates that the model fit could be",
                          "improved by using a larger setting of the",
                          "prior variance. Consider increasing the range",
                          "of the variances \"sa2\"."),1/K))
  }

  return(out)
}
