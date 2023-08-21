#### Functions for WMM and ARMM ####
library(ddalpha)
library(forecast)
library(LaplacesDemon)
library(matrixcalc)
library(mclust)

#### Fill in missing value of a time series with NA ####
fill_missing <- function(t, y){ 
  
  # t: time point index e.g. 1 to n
  # y: time series, repeated observation
  
  df1 <- data.frame("t"=t,
                    "y"=y)
  df2 <- data.frame("t"=c(1:max(df1$t)),
                    "temp"=rep(NA, max(df1$t)) )
  df3 <- merge(x = df1, y = df2, by="t", all.y=T )[c("t","y")]
  return(df3)
}

#### convert a time series into a data matrix ####
data_matrix <- function(y, lag){
  
  # y: time series, repeated observation
  # lag: number of past time points used in model e.g k-1
  
  lag <- lag+1
  df = NULL
  for (i in 1:lag) {
    temp <-c( rep(NA,i-1), y[1:(length(y)-i+1)] )
    df <- cbind( df, temp )
  }
  return(df)
}

#### turn autocorrelation vector into autocorrelation matrix ####
cov_to_matrix <- function(covar){
  
  # covar: autocovariance values, can be obtained using acf()
  
  mat <- matrix(covar, nrow = 1)
  dim <- length(covar)
  for (j in 2:dim){
    mat <- rbind( mat, 
                  matrix( c(rep(0,j-1), covar[1:(dim-j+1)] ), nrow = 1) 
    )
  }
  mat <- mat + t(mat) - diag(diag(mat))
  return(mat)
}

#### calculates AR coefficient phi using Yule-Walker ####
AR_coeff <- function(matrix){
  
  # matrix: autocorrelation/autocovariance/sigma matrix
  
  d <- dim( matrix )[1]
  return( as.numeric( solve( matrix[2:d,2:d] ) %*% matrix[2:d,1] ) )
}

#### calculate mean from a list of matrices ####
matrix_sum <- function(matrix_list){
  out <- 0*matrix_list[[1]]
  for(i in 1:length(matrix_list) ){
    out = out + matrix_list[[i]]
  }
  return( out )
}

#### EM functions for all algorithms ####
#### E step ####

# E-step for z, returns a list of length G, with vector z of length I estimate for each group ####
E_step_z <- function(matrix_list, pi_list, sigma_list,
                     n_vector, z_list, lambda_list = NULL,
                     method = "EM1"){
  
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # pi_list: list of length G, estimated pi values for each group
  # sigma_list: list of length G, estimated sigma matrices for each group
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  # z_list: list of length G, with estimated vector z of length I for each group
  # lambda_list: lambda value for each group
  # method: WMM model, EM1, EM2, EM3
  
  #
  I <- length(matrix_list)
  G <- length(pi_list)
  #
  if(method == "EM1"){
    new_list <- list()
    for(i in 1:I){
      new_list[[i]] <- list( "matrix" = matrix_list[[i]], "df"=n_vector[i])
    }
    #
    denom <- rep(0,I)
    # calculate the denominator
    for( g in 1:G ){
      denom <- denom + pi_list[[g]] * 
        sapply(new_list, function(x) LaplacesDemon::dwishart(x$matrix, x$df, sigma_list[[g]]) )
    }
    # calculate the numerator and estimate z for each group 1,...,G
    for( g in 1:G ){
      numer <- pi_list[[g]] * 
        sapply(new_list, function(x) LaplacesDemon::dwishart(x$matrix, x$df, sigma_list[[g]]) )
      out <- numer/denom
      out[denom==0] <- 1/G
      z_list[[g]] <- out
    }    
  }
  if(method == "EM2"){
    new_list <- list()
    for(i in 1:I){
      new_list[[i]] <- list( "matrix" = matrix_list[[i]], "df"=n_vector[i])
    }
    #
    denom <- rep(0,I)
    for( g in 1:G ){
      denom <- denom +
        pi_list[[g]] * 
        sapply(new_list, function(x) LaplacesDemon::dwishart(x$matrix, x$df*lambda_list[[g]], sigma_list[[g]]) )
    }
    for( g in 1:G ){
      numer <- pi_list[[g]] * 
        sapply(new_list, function(x) LaplacesDemon::dwishart(x$matrix, x$df*lambda_list[[g]], sigma_list[[g]]) )
      out <- numer/denom
      out[denom==0] <- 1/G
      z_list[[g]] <- out
    }
  }
  #
  return(z_list)
}

#### M step ####
# M-step for pi, return a list of length G, with scalar pi estimates for each group ####
M_step_pi <- function(z_list, method = NULL){
  
  # z_list: list of length G, with estimated vector z of length I for each group
  # method: WMM model, EM1, EM2, EM3
  
  G <- length(z_list)
  pi_list <- list()
  #
  for(g in 1:G){
    pi_list[[g]] <- mean(z_list[[g]])
  }
  return(pi_list)
}

# M-step for sigma, return a list of length G, with matrix sigma estimate for each group ####
M_step_sigma <- function(matrix_list, z_list, 
                         n_vector = NULL, lambda_list = NULL,
                         method = "EM1"){
  
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # z_list: list of length G, with estimated vector z of length I for each group
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  # lambda_list: lambda value for each group
  # method: WMM model, EM1, EM2, EM3
  
  G <- length(z_list)
  I <- length(matrix_list)
  sigma_list <- list()
  #
  for(g in 1:G){
    numer <- 0
    denom <- 0
    if(method == "EM1"){
      for(i in 1:I){
        numer <- numer + matrix_list[[i]] * z_list[[g]][i]
        denom <- denom + n_vector[i] * z_list[[g]][i]
      }
    }
    #
    if(method == "EM2"){
      for(i in 1:I){
        numer <- numer + matrix_list[[i]] * z_list[[g]][i]
        denom <- denom + (n_vector[i] * lambda_list[[g]]) * z_list[[g]][i]
      }
    }
    #
    sigma_list[[g]] <- numer/denom
  }
  return(sigma_list)
}

# Sum of digamma function with lambda used in the M-step of lambda, EM2 ####
sum_digamma_lambda <- function(dim, z, n_vector, lambda){
  
  # dim: dimension of wishart 
  # z: z estimate for group g
  # n_vector: degree of freedom vector for subjects
  # lambda: scalar value for group g
  
  C <- 0
  for( j in 1:dim){
    C <- C + sum( z * n_vector * digamma(0.5*(n_vector*lambda-j+1) ) )
  }
  return( C )
}

# M-step for lambda_g, return a list of length G, with scalar lambda estimate for each group ####
M_step_lambda <- function(matrix_list, z_list, sigma_list, n_vector, lambda_list, upper=2){
  
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # z_list: list of length G, with estimated vector z of length I for each group
  # sigma_list: list of length G, estimated sigma matrices for each group
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  # lambda_list: lambda value for each group
  
  I <- length(matrix_list)
  G <- length(z_list)
  dim <- dim(sigma_list[[1]])[1]
  lower = (dim) / min(n_vector)
  # lambda_list
  output_list <- list()
  #
  for(g in 1:G){
    A <- 0
    for(i in 1:I){
      A <- A + z_list[[g]][i]*n_vector[i]*log( det( matrix_list[[i]] %*% solve(sigma_list[[g]])/2 ) )
    }
    output <- optim( lambda_list[[g]], # use square loss
                     function(x) (A-sum_digamma_lambda(dim, z_list[[g]], n_vector, x))**2, 
                     lower = lower, upper = upper, method = "L-BFGS-B")
    output_list[[g]] <- as.numeric(output[1])
  }
  return( output_list )
}

#### EM algorithm for WMM ####
WMM <- function(matrix_list, n_vector = NULL, G = 2,
                method = "EM1", upper=2, tol = 1e-5,
                burn_in = 100, max_iter = 1000,
                pi_list_init = NULL, sigma_list_init = NULL, 
                z_list_init = NULL, lambda_list_init = NULL){
  
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  # G: total number of groups
  # method: method, EM1, EM2, EM3
  # upper: upper bound for n_g or lambda_g
  # tol: log likelihood gain convergence threshold
  # burn_in: minimum number of iterations
  # max_iter: maximum number of iterations
  
  # pi_list_init: initial values
  # sigma_list_init: initial values
  # z_list_init: initial values
  # lambda_list_init: initial values
  
  # initialize pi values
  if( is.null(pi_list_init) ){
    pi_list <- as.list( LaplacesDemon::rdirichlet( 1,rep(100,G) ) )
  } else {  
    pi_list <- pi_list_init
  }
  
  # initialize z vector
  if( is.null(z_list_init) ){
    
    indices <- split(1:length(matrix_list),
                     f = sample(1:G, size = length(matrix_list), replace = T) )
    z_list <- list()
    
    for(g in 1:G){
      
      ind <- array( indices[ as.character(g) ] )[[1]]
      z_list[[g]] <- 1*( 1:length(matrix_list) %in% ind )
      
    }
  } else {
    z_list <- z_list_init
  }
  
  # initialize sigma matrices
  if( is.null(sigma_list_init) ){
    
    sigma_list <- list()
    
    for(g in 1:G){
      
      ind <- z_list[[g]]==1
      
      sigma_list[[g]] <- matrix_sum( matrix_list[ ind ] )/sum(n_vector[ind])
    }
  } else {
    sigma_list <- sigma_list_init
  }
  
  # initialize lambda_g values
  if( is.null(lambda_list_init) ){
    lambda_list <- as.list( rep( 1, G ) )
  } else {
    lambda_list <- lambda_list_init
  }
  
  #
  sigma_old <- sigma_list_init
  converged <- F
  for(iter in 1:max_iter){
    
    # E step #
    z_list <- E_step_z(matrix_list=matrix_list, 
                       pi_list=pi_list, sigma_list=sigma_list,
                       n_vector = n_vector, z_list = z_list, 
                       lambda_list = lambda_list,
                       method = method)

    # M step #
    pi_list <- M_step_pi(z_list, method = NULL)
    #
    sigma_list <-  M_step_sigma(matrix_list=matrix_list, z_list=z_list, 
                                n_vector = n_vector, lambda_list = lambda_list,
                                method = method)
    #
    if(method == "EM2"){
      lambda_list <- M_step_lambda(matrix_list=matrix_list, z_list=z_list, 
                                   sigma_list=sigma_list, n_vector=n_vector, 
                                   lambda_list=lambda_list, upper=upper)
    }
    # Stopping Criteria
    sigma_delta <- 0
    for(g in 1:G){
      sigma_delta <- sigma_delta + 
        sum( abs( as.vector( sigma_old[[g]] - sigma_list[[g]] ) ) )
    }
    sigma_old <- sigma_list
    #
    #if(iter %% 10 == 0){
    if(F){
      print("Iteration number:")
      print(iter)
      print("lambda list")
      print(lambda_list)
      print("sigma list")
      print(sigma_list)
    }
    if(iter > burn_in & sigma_delta <= tol){
      converged <- T
      break
    }
  }
  output <- list(
    "method" = method,
    "z_list" = z_list,
    "pi_list" = pi_list,
    "sigma_list" = sigma_list,
    "lambda_list" = lambda_list,
    
    "converged" = converged,
    "iter" = iter
  )
  return(output)
}

#### calculates maximum a posteriori of z ####
z_max <- function(z_list){
  
  # z_list: list of length G, with estimated vector z of length I for each group
  
  I <- length(z_list[[1]])
  G <- length(z_list)
  z_map <- diag( rep(1,G) )
  #
  z_df <- NULL
  for (g in 1:G){
    # convert z_list to a I x G df
    z_df <- cbind( z_df, z_list[[g]] )
  }
  for(i in 1:I){
    # obtain argmax group index, from I x G df
    z_df[i,] <- z_map[which.max(z_df[i,]),]
  }
  return(z_df)
}

#### calculates estimate phi and asymptotic dist. ####
summary_coeff <- function(z_list,matrix_list,sigma_list,n_vector){
  
  # z_list: list of length G, with estimated vector z of length I for each group
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # sigma_list: list of length G, estimated sigma matrices for each group
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  
  I <- length(z_list[[1]])
  G <- length(z_list)
  K <- dim(sigma_list[[1]])[1]
  #
  out_list <- list()
  for(g in 1:G){
    #
    phi_g <- solve( sigma_list[[g]][2:K,2:K] ) %*% sigma_list[[g]][2:K]
    #
    bread <- sigma_list[[1]][2:K,2:K]*0
    center <- sigma_list[[1]][2:K,2:K]*0
    for(i in 1:I){
      bread <- bread + z_list[[g]][i]*matrix_list[[i]][2:K,2:K]
      #
      v_i <- matrix_list[[i]][1,1]/n_vector[i] * 
        (1 - t(sigma_list[[g]][2:K]) %*% 
           solve( sigma_list[[g]][2:K,2:K] ) %*% 
           sigma_list[[g]][2:K] / (sigma_list[[g]][1,1]) )
      if(v_i < 0 ){
        print(i)
        print(g)
      }
      center <- center + as.numeric(v_i) * z_list[[g]][i]**2 * matrix_list[[i]][2:K,2:K]
    }
    out_list[[g]] <- list("coeff"=phi_g,
                          "var" = solve(bread) %*% center %*% solve(bread) )
  }
  return(out_list)
}

#### calculates log likelihood of ARMM model using sigma_list ####
logLik_ARMM <- function(z_list, matrix_list, sigma_list, n_vector){
  
  # z_list: list of length G, with estimated vector z of length I for each group
  # matrix_list: list of length I, estimated autocorrelation matrices for each subject
  # sigma_list: list of length G, estimated sigma matrices for each group
  # n_vector: degree of freedom values for each subject i as a vector, i.e. n_i
  
  # initial log likelihood
  ll <- 0
  I <- length(z_list[[1]])
  G <- length(z_list)
  K <- dim(sigma_list[[1]])[1]
  z_map <- diag( rep(1,G) )
  #
  z_df <- NULL
  for (g in 1:G){
    # convert z_list to a I x G df
    z_df <- cbind( z_df, z_list[[g]] )
  }
  for(i in 1:I){
    # obtain argmax group index, from I x G df
    z_df[i,] <- z_map[which.max(z_df[i,]),]
    for (g in 1:G){
      # estimate alpha_i and tau_i for each subject
      if(z_df[i,g]==1){
        v_i <- matrix_list[[i]][1,1]/n_vector[i] * (1 - t(sigma_list[[g]][2:K]) %*% 
          solve( sigma_list[[g]][2:K,2:K] ) %*% 
          sigma_list[[g]][2:K] / (sigma_list[[g]][1,1]) )
        if(v_i < 0 ){
          print(i)
          print(g)
        }
        # add to the log likelihood for each subject
        ll <- ll + n_vector[i] * log(v_i)
      }
    }
  }
  return(ll)
}

#### TSclust functions ####

diss.ACF <- function(x, y, lag.max = 50){
  rhox <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
  rhoy <- acf(y, lag.max = lag.max, plot = FALSE)$acf[-1]
  return(sqrt(t(rhox - rhoy) %*% diag(length(rhox)) %*% (rhox - rhoy)))
}

diss.PACF <- function(x, y, lag.max = 50){
  rhox <- as.vector(pacf(x, lag.max = lag.max, plot = FALSE)$acf)
  rhoy <- as.vector(pacf(y, lag.max = lag.max, plot = FALSE)$acf)
  return(sqrt(t(rhox - rhoy) %*% diag(length(rhox)) %*% (rhox - rhoy)))
}

diss.AR.PIC <- function(x.mat, y.mat){
  return( as.numeric(dist(rbind(AR_coeff(x.mat), AR_coeff(y.mat)))) )
}
#
diss.ACF2 <- function(rhox, rhoy){
  return(sqrt(t(rhox - rhoy) %*% diag(length(rhox)) %*% (rhox - rhoy)))
}

diss.PACF2 <- function(rhox, rhoy){
  return(sqrt(t(rhox - rhoy) %*% diag(length(rhox)) %*% (rhox - rhoy)))
}
diss.AR.PIC2 <- function(coef1, coef2){
  return(sqrt(t(coef1 - coef2) %*% diag(length(coef1)) %*% (coef1 - coef2)))
}

#### Simulation functions ####
sim <- function(ar=c(0,0), ma=c(0,0), sigma_i, n_i, I=50){
  y <- NULL
  for(i in 1:I){
    y_i <- arima.sim( model=list(ar=ar,ma=ma, mean=0, sd=sigma_i), n = n_i, n.start = 100 )
    y <- cbind(y,y_i)
  }
  colnames(y) <- NULL
  return(y)
}

sim_fit <- function(z_true,
                    y, # time series data as columns
                    return_data=T){
  #
  matrix_list <- list()
  coeff_df <- NULL
  acf_df <- NULL
  pacf_df <- NULL
  acf_dist <- diag(rep(0,dim(y)[2]))
  pacf_dist <- diag(rep(0,dim(y)[2]))
  pic_dist <- diag(rep(0,dim(y)[2]))
  n_vector <- c()
  #
  for(i in 1:dim(y)[2]){
    n_vector <- c(n_vector,sum(!is.na(y[,i])))
    matrix_list[[i]] <- sum(!is.na(y[,i]))*
      cov_to_matrix(acf(y[,i], lag.max = 2, plot = F, 
                        type = "covariance",
                        na.action = na.pass)$acf)
    #
    coeff <- arima(y[,i],order = c(2,0,0),method = "ML",include.mean = F)$coef[1:2]
    coeff_df <- rbind(coeff_df, coeff)
    #
    acf_i <- acf(y[,i],lag.max = 2,plot = F,na.action = na.pass)$acf[2:3]
    acf_df <- rbind(acf_df, acf_i)
    #
    pacf_i <- pacf(y[,i],lag.max = 2,plot = F,na.action = na.pass)$acf
    pacf_df <- rbind(pacf_df,pacf_i)
  }
  #
  for(r in 1:dim(y)[2] ){
    for( c in 1:dim(y)[2] ){
      pic_dist[r,c] <- diss.AR.PIC2(coef1 = coeff_df[r,],
                                    coef2 = coeff_df[c,])
      acf_dist[r,c] <- diss.ACF2(acf_df[r,],acf_df[c,])
      pacf_dist[r,c] <- diss.PACF2(pacf_df[r,],pacf_df[c,])
    }
  }
  #
  acf_dist <- as.dist(acf_dist)
  pacf_dist <- as.dist(pacf_dist)
  pic_dist <- as.dist(pic_dist)
  #
  out_acf <- cutree(hclust(acf_dist), k = 2)
  out_pacf <- cutree(hclust(pacf_dist), k = 2)
  out_pic <- cutree(hclust(pic_dist), k = 2)
  GMM <- Mclust(data = coeff_df, G=2)
  
  out1 <- WMM(matrix_list=matrix_list, n_vector = n_vector, G = 2,
              method = "EM1",
              burn_in = 10, max_iter = 1000,
              pi_list_init = NULL, sigma_list_init = NULL, 
              z_list_init = NULL, 
              lambda_list_init = NULL)
  out2 <- WMM(matrix_list=matrix_list, n_vector = n_vector, G = 2,
              method = "EM2",
              burn_in = 10, max_iter = 200,
              pi_list_init = NULL, 
              sigma_list_init = out1$sigma_list, 
              z_list_init = NULL, 
              lambda_list_init = list(1.0,1.0) )
  #####
  pGMM <- mean( (GMM$z[,1] > 0.5) == z_true )
  pGMM <- max(pGMM,1-pGMM)
  
  pACF <- mean( (out_acf==1) == z_true )
  pACF <- max(pACF,1-pACF)
  
  pPACF <- mean( (out_pacf==1) == z_true )
  pPACF <- max(pPACF,1-pPACF)
  
  pPIC <- mean( (out_pic==1) == z_true )
  pPIC <- max(pPIC,1-pPIC)
  
  pEM1 <- mean( (out1$z_list[[1]] > 0.5) == z_true )
  pEM1 <- max(pEM1,1-pEM1)
  
  pEM2 <- mean( (out2$z_list[[1]] > 0.5) == z_true )
  pEM2 <- max(pEM2,1-pEM2)
  
  output <- list(
    "pACF" = pACF,
    "pPACF" = pPACF,
    "pPIC" = pPIC,
    "pGMM" = pGMM,
    "pEM1" = pEM1,
    "pEM2" = pEM2
  )
  if(return_data){
    output$data <- y
  }
  return(output)
}
