#' UPoSI confidence intervals for penalized G-estimator
#' @description This function provides UPoSI confidence intervals (considering random design) for the
#' blip coefficients selected by penalized G-estimation in the context of a structural nested mean
#' model (SNMM) for repeated outcomes.
#' @param data A data frame containing the variables in longitudinal format. In the data, the
#' outcome should be continuous and the treatment/exposure should be binary. The data must have 
#' an additional column containing the propensity scores.
#' @param wc.str A character string specifying the working correlation structure. The
#' following are currently allowed: "independence", "exchangeable", "ar1", and "unstructured".
#' @param id.var The column name in data that corresponds to the variable id (unique identifier).
#' @param response.var The column name in data that corresponds to the response variable.
#' @param treat.var The column name in data that corresponds to the treatment/exposure variable.
#' @param tf.model A single formula object specifying the covariates of a (linear)
#' treatment-free model.
#' @param theta.tilde The vector of regression estimates for the selected model. i.e., the
#' estimates for the selected blip coefficients and the treatment-free model parameters.
#' @param alpha.hat The estimated correlation parameter(s) alpha(s) if the provided structure is either
#' "exchangeable", "ar1, or "unstructured". For unstructured, the elements of alpha.hat correspond
#' to the upper triangular portion of working correlation matrix having dimension equal to
#' the largest cluster size.
#' @param sigma2.hat The estimated variance parameter sigma^2.
#' @param test.size The significance level. For a 95\% confidence interval test.size should be 0.05.
#' @param nrep Number of bootstrap replicates.
#' @param M A logical vector of TRUE/FALSE values indicating which variables are selected for
#' the blip model.
#' @param r A list of length 'nrep', where each element is a numeric vector of length 'n'
#' containing random draws from an Exponential(1) distribution, with 'n' equal to the number
#' of subjects in the data.
#' @param continuous.covs  A logical vector of TRUE/FALSE values identifying the potential
#' effect modifiers that are continuous.
#'
#' @return A matrix showing confidence interval estimates for the selected blip coefficients.
#' @export
#'
#' @importFrom stats model.matrix as.formula binomial glm predict.glm toeplitz
#'
#' @examples
#' library(mvtnorm)
#' expit <- function(x) exp(x)/(1+exp(x))
#'
#' ## data.gen is a function that generates a longitudinal data set for a specific correlation
#' ## structure. Available structures in this function are: independence, exchangeable and AR1.
#'
#' ## Arguments(data.gen):
#' #     n = Number of subjects
#' #     ni = A vector containing number of time points for each subject
#' #     sigma2.e = Error variance
#' #     alpha = Correlation parameter
#' #     corstr = The correlation structure among the repeated outcomes
#' #     autocorr.coef = The autocorrelation coefficient for inducing correlation among the
#' #     continuous confounders and the noise covariates
#'
#' data.gen <- function(n, ni, sigma2.e, alpha, corstr, autocorr.coef){
#'  ncovs  <-  2+4+45 # 2+4=6 confounders and 45 noise covariates
#'   beta <- c(0, 1, -1.1, 1.2, 0.75, -0.9, 1.2) # treatment model parameters
#'   delta <- c(1, 1, 1.2, 1.2, -0.9, 0.8, -1,
#'              rep(1, 20), rep(0, ncovs-6-20),
#'              -0.8, 1, 1.2, -1.5) # treatment-free model parameters
#'   psi <- c(1, 1, -1, -0.9, 0.8, 1, 0, rep(0, 20), rep(0, ncovs-6-20)) # blip parameters
#'
#'   # generating two continuous baseline covariates
#'   l1 <- rnorm(n, 0, 1)
#'   l2 <- rnorm(n, 0, 1)
#'
#'   # V is the covariance matrix of the time-varying confounders (l3,..,l6) and
#'   # noise covariates (x1,...)
#'   V <- toeplitz(autocorr.coef^(0:(ncovs-2-1)))
#'
#'   lx <- a <- y <- vector(mode="list", length=n)
#'   lx.mat <- NULL
#'   for(i in 1:n){
#'     a[[i]] <- y[[i]] <- rep(NA, ni[i])
#'     lx[[i]] <- matrix(NA, ni[i], ncovs)
#'     lx[[i]][,1] <- rep(l1[i], ni[i])
#'     lx[[i]][,2] <- rep(l2[i], ni[i])
#'
#'     corr.mat <- switch(corstr,
#'                        "exchangeable" = toeplitz(c(1, rep(alpha, ni[i]-1))),
#'                        "ar1" = toeplitz(alpha^(0:(ni[i]-1))),
#'                        "independence" = diag(ni[i])
#'     )
#'     cov.mat <- diag(sqrt(sigma2.e), ni[i]) %*% corr.mat %*% diag(sqrt(sigma2.e), ni[i])
#'     e <- rmvnorm(1, sigma = cov.mat)
#'
#'     # j=1
#'     mu.lx <- c(NA, NA, # for l1 and l2
#'                rep(0, 4), rep(0, ncovs-6)) #rep(0.3*lx[[i]][1, 1]+0.3*lx[[i]][1, 2],4)
#'     lx[[i]][1,3:ncovs] <- rmvnorm(1, mean=mu.lx[3:ncovs], sigma = V)
#'     a[[i]][1] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][1,1:6])*beta)))  # no correlation
#'     tf.mean <- sum(c(1, lx[[i]][1,1:ncovs],
#'                      lx[[i]][1,1]*lx[[i]][1,5], lx[[i]][1,3]*lx[[i]][1,4], sin(lx[[i]][1,3]-lx[[i]][1,4]),
#'                      cos(2*lx[[i]][1,5])) * delta)
#'     blip <- sum(c(a[[i]][1], a[[i]][1]*c(lx[[i]][1,1:ncovs])) * psi)
#'     y[[i]][1] <- (tf.mean + blip) + e[1]
#'
#'
#'     # j=2:ni
#'     for(j in 2:ni[i]){
#'       mu.lx <- c(NA, NA, # for l1 and l2
#'                  0.3*lx[[i]][j-1, 3:6] + 0.3*a[[i]][j-1], 0.5*lx[[i]][j-1,7:ncovs])
#'       lx[[i]][j,3:ncovs] <- rmvnorm(1, mean=mu.lx[3:ncovs], sigma = V)
#'       a[[i]][j] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][j,1:6])*beta)))  # no correlation
#'       tf.mean <- sum(c(1, lx[[i]][j,1:ncovs],
#'                        lx[[i]][j,1]*lx[[i]][j,5], lx[[i]][j,3]*lx[[i]][j,4], sin(lx[[i]][j,3]-lx[[i]][j,4]),
#'                        cos(2*lx[[i]][j,5])) * delta)
#'       blip <- sum(c(a[[i]][j], a[[i]][j]*c(lx[[i]][j,1:ncovs])) * psi)
#'       y[[i]][j] <- (tf.mean + blip) + e[j]
#'     }
#'     lx.mat <- rbind(lx.mat, lx[[i]])
#'   }
#'
#'   colnames(lx.mat) <- c(paste("l", 1:6, sep=""), paste("x", 1:(ncovs-6), sep=""))
#'   data <- data.frame(id=rep(1:n, times=ni), a=unlist(a), lx.mat, y=round(unlist(y),3))
#'   return(data)
#' }
#'
#' data.s <- data.gen(n = 500, ni = rep(6, 500), sigma2.e = 1, alpha = 0.8,
#'                    corstr = "exchangeable", autocorr.coef = 0.25)
#'
#' ncovs  <-  2+4+45
#' #treatment-free model is misspecified
#' tf.model <- as.formula(paste("~",paste(c(paste("l",1:6,sep=""), paste("x",c(1:9,11:(ncovs-6)),sep="")),
#'                                        collapse = "+"), collapse=""))
#' treat.model <- ~l1+l2+l3+l4+l5+l6
#'
#' lam_max <- 1
#' lam_min <- 0.01*lam_max
#' lambda.seq <- sort(seq(lam_min,lam_max,(lam_max-lam_min)/99), decreasing=TRUE)
#'
#' # library(devtools) # if already installed, otherwise need to install it first
#' # install_github("ajmeryjaman/penalizedG") # run this command if this package is not already installed
#' library(penalizedG)
#'
#' # Run penalized G-estimation for a sequence of tuning parameters (lambda)
#' out <- penalizedG(data = data.s, wc.str = "exchangeable", id.var="id", response.var="y",
#'                    treat.var="a", tf.model=tf.model, treat.model = treat.model,
#'                    lambda.seq = lambda.seq, maxitr = 100, penalty = "SCAD")
#' out$selected.EMs
#'
#' data <- out$data # have an additional column containing the propensity scores
#' n <- length(unique(data$id))
#' p <- length(model.matrix(tf.model, data=data)[1,])
#' effect.modification <- out$estimate[2:p]
#' delta <- out$estimate[(p+1):(2*p)]
#' theta.tilde <- c(out$estimate[1], effect.modification[abs(effect.modification) > 0.001], delta)
#' nrep.boot <- 500
#' selected.EM <- c(TRUE, abs(effect.modification) > 0.001)
#' set.seed(100)
#' random.numb.exp <- lapply(1:nrep.boot, function(B) rexp(n,1))
#'
#' # UPoSI confidence interval (random design)
#' CI.boot.psi <- uposi_confint(data = data, wc.str = "exchangeable", id.var = "id", response.var = "y",
#'                              treat.var = "a", tf.model = tf.model, theta.tilde = theta.tilde,
#'                              alpha.hat = out$alpha.hat, sigma2.hat = out$sigma2.hat,
#'                              test.size = 0.05, nrep = nrep.boot, M = selected.EM, r = random.numb.exp,
#'                              continuous.covs = rep(TRUE, length(all.vars(tf.model))))
#' CI.boot.psi
#'
#' # Naive CI based on sandwich variance
#' M <- selected.EM
#' CI.naive <- cbind(theta.tilde[1:sum(M)] - qnorm(0.975)*sqrt(diag(out$asymp.var.psi[M,M])),
#'                  theta.tilde[1:sum(M)] + qnorm(0.975)*sqrt(diag(out$asymp.var.psi[M,M])))
#' CI.naive

uposi_confint <- function(data, wc.str, id.var, response.var, treat.var, tf.model, theta.tilde,
                          alpha.hat, sigma2.hat, test.size, nrep, M, r, continuous.covs){
  ## theta.tilde := (only the selected psi's, delta)
  ## continuous.covs : vector with T/F indicating continuous.covs covs in tf.model
  names(data)[names(data)==id.var] <- "id"
  names(data)[names(data)==treat.var] <- "a"
  names(data)[names(data)==response.var] <- "y"

  # n, K >> scalar
  # ni, l, a, y, l.mat.split, e >> list of length n

  n <- length(split(data, data$id))
  p <- length(model.matrix(tf.model, data)[1,])
  dat <- split(data, data$id)
  ni <- lapply(1:n, function(i) length(dat[[i]]$id))
  a <- lapply(1:n, function(i) as.matrix(dat[[i]]$a))
  E.a <- lapply(1:n, function(i) as.matrix(dat[[i]]$E.a))
  y <- lapply(1:n, function(i) as.matrix(dat[[i]]$y))
  l.mat.full <- data.frame(id=data$id,model.matrix(tf.model,data)) # cov+treat history
  l.mat.split.full <- split(l.mat.full, l.mat.full$id)
  l.full <- lapply(1:n, function(i) as.matrix(l.mat.split.full[[i]][,-1]))
  l.mat <- data.frame(id=data$id,model.matrix(tf.model,data)[,M]) # cov+treat history
  l.mat.split <- split(l.mat, l.mat$id)
  l <- lapply(1:n, function(i) as.matrix(l.mat.split[[i]][,-1])) # according to selected model
  K <- dim(l.mat.full[,-1])[2] # dimension of psi/delta including intercept
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))

  sum3.n <- Reduce("+", lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]]))) #(c(a[[i]]-E.a[[i]])*l[[i]])))
  sum4.n <- Reduce("+", lapply(1:n, function(i) t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%l.full[[i]]))
  sum6.n <- Reduce("+", lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%l.full[[i]]))
  sum5.n <- Reduce("+", lapply(1:n, function(i) t(l.full[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]]))) ##t(sum4.n)
  W_M <- rbind(cbind(sum3.n,sum4.n),cbind(sum5.n,sum6.n))
  for(k in all.vars(tf.model)[continuous.covs]){
    data[,k] <- (data[,k] - mean(data[,k]))/sd(data[,k])
  }

  ### bivariate quantiles should be calculated with standardized continuous covariates
  l.mat.full <- data.frame(id=data$id, model.matrix(tf.model,data))
  l.mat.split.full <- split(l.mat.full, l.mat.full$id)
  l.full <- lapply(1:n, function(i) as.matrix(l.mat.split.full[[i]][,-1]))
  QB <- bivar_quantile(n, p, a, E.a, y, l.full, K, V.inv, theta.tilde, test.size, nrep, M, r)

  e_l <- diag(rep(1, dim(W_M)[1]))
  mltp <- apply(t(e_l)%*%solve(W_M), 1, function(x) sum(abs(x)))
  half.length <- mltp*(QB$C_G + QB$C_W*sum(abs(theta.tilde)))
  CI.boot.blip <- matrix(cbind(theta.tilde-half.length,
                               theta.tilde+half.length)[-tail(1:length(half.length), p),],sum(M),2)

  colnames(CI.boot.blip) <- c("lower", "upper")
  selected.EMs <- colnames(model.matrix(tf.model, data))[M][-1]
  row.names(CI.boot.blip) <- c("a", paste("a*", selected.EMs, sep=""))
  return(CI.boot.blip)
}
