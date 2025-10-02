#################Defining some kernels and kernel estimators##############

#' Gamma kernel 'of y' of bandwidth h estimated at x
#' (parameter of gamma distribution depend on y)
#' @param y array of jumping times
#' @param x value at which it is estimated
#' @param h bandwidth
#' @returns value of y,h-dependent gamma kernel evaluated at point x
gamma_kernel  = function(y, x, h) {
  indices = (y >= 2*h)
  y_new = y/h*indices + (1/4*(y/h)^2 + 1)*(!indices)
  h*dgamma(x, y_new , scale = h)
} #gamma kernel for non-parametric estimation

gamma_ker = function(y,x,h)(gamma_kernel(y,x,h)/h)


#' Gaussian kernel estimator of hazard rate of jumping times in X with bandwidth
#' h at point t
#' @param t point at which to evaluate kernel estimation
#' @param X array of jumping times
#' @param h bandwidth
#' @returns kernel estimator with bandwidth h of hazard rate of jumping times X at t 
ker_estg <- function(t,X,h){
  SF = length(X):(1+length(X)-length(sort(X))) #empirical survival function
  X2 = (t-sort(X))/h
  expX = exp(-X2^2/2)/(sqrt(2*pi))
  return(1/h*sum(1/SF*expX))
} 


#' Auxiliary function for lognormal kernel estimator of hazard rate 
#' of jumping time x with bandwidth h at point t
#' @param x jumping time
#' @param t point at which to evaluate kernel estimation
#' @param h bandwidth
#' @returns lognormal kernel with bandwidth h and jumping time x at t 
log_ker = function(x,t,h) (1/x*(dnorm(log(x), log(t), h)))

#' Lognormal kernel density estimator of jumping times in X with bandwidth h at point t
#' @param t point at which to evaluate kernel estimation
#' @param h bandwidth
#' @param X array of jumping times
#' @returns lognormal kernel estimation of hazard rate of the jumping times 
#' in X at time t
log_ker_dens= function(t,h,X) (
  sum(sapply(sort(X), function(u) (log_ker(t,u,h)))/(length(sort(X)))))


#' Estimation of the distribution function, usefull for the ratio estimator
#' of the hazard rate
#' @param t point at which to evaluate kernel estimation
#' @param h bandwidth
#' @param X array of jumping times
#' @returns estimation of the distribution function of jumping times in X at t
est_surv = function(t,h,X)(
  sum(sapply(X,function(u)( integrate(function(y)(log_ker(y,u,h)), 0,t)$value) ))/length(X))



#' Gamma kernel estimator of hazard rate of jumping times in X with fixed bandwidth
#' h at point t
#' @param t point at which to evaluate kernel estimation
#' @param X array of jumping times
#' @param h bandwidth
#' @returns gamma kernel estimator with fixed bandwidth h of hazard rate of jumping times X at t 
ker_est_gamma_c <- function(t,X,h){
  return(sum(gamma_kernel(t, sort(X),h)/(length(X):(1+length(X)-length(sort(X)))))/h)
} # gamma kernel constant bandwidth


#' Gaussian kernel estimator of hazard rate of jumping times in X with k-nearest
#' neighbor bandwidth
#' @param t point at which to evaluate kernel estimation
#' @param X array of jumping times
#' @param k number of neighbors for bandwidth
#' @returns gaussian kernel estimator with k-nearest neighbor bandwidth of hazard rate of jumping times X at t 
ker_est_neigh <- function(t,X,k){
  X = X[!is.na(X)]
  SF = length(X):1 #empirical survival function
  dis = sort((X-t)^2)[k] #distance of the kth nearest neighbor
  h = sqrt(dis)
  X = (t-sort(X))/h
  expX = exp(-X^2/2)/(sqrt(2*pi))
  return(1/h*sum(1/SF*expX))
}
#########Simulation function#############
#'Returns the death time for a population of size init_size in a one phase model
#'with an exponential hazard rate a+bexp(c*t)
#'@param a the first parameter in the hazard rate
#'@param b the second parameter
#'@param c the third parameter
#'@param init_size the population size
#'@param n the time of simulation divided by 10
#'@returns an array of death times 
sim_pop_simple <- function(a,b,c,init_size,n=200){
  pop_init <- population(data.frame( birth = rep(0, init_size), 
                                     death = NA)) #initial pop
  death = mk_event_individual(type = 'death', 
                              intensity_code ="{result = a+b*exp(c*I.age(t));}"
  )
  params =  list('a'=a,'b'=b,'c'=c) 
  birth_death <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = params)
  sim_out <- popsim(model = birth_death,
                    initial_population = pop_init,
                    events_bounds = c(
                      'death' = 1),
                    parameters = params,
                    time = (0:n)*10)
  return(sim_out$population[[n]]$death)
}


#'Returns the death time for a population of size init_size in a one phase model
#'with a hazard rate defined by three different parameters a,b,c whose expression
#'is given in text of the form '{result = expression;}' 
#'@param a the first parameter in the hazard rate
#'@param b the second parameter
#'@param c the third parameter
#'@param text a string of the expression of the hazard rate in C++ of the form 
#''{result = expression;}' where the variable is I.age(t)
#'@param init_size the population size
#'@param n the time of simulation divided by 10
#'@returns an array of n death times/events with hazard defined by text
sim_pop_text <- function(params, text,init_size,n=200){
  pop_init <- population(data.frame( birth = rep(0, init_size), 
                                     death = NA))
  death = mk_event_individual(type = 'death', 
                              intensity_code =text
  )
  birth_death <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = params)
  sim_out <- popsim(model = birth_death,
                    initial_population = pop_init,
                    events_bounds = c(
                      'death' = 1),
                    parameters = params,
                    time = n*10)
  return(sim_out$population$death)
}
##################Cross validation and max likelihood ################
#Both of the following functions are partly taken from the ones i the package kernhaz
#https://github.com/cran/kernhaz/blob/master/R/crossval_haz.R


#' s-block cross-validation of an estimator of hazard rate of jumping times 
#' @param Y array of jumping times
#' @param h bandwidth of kernel estimator OR number of neighbors for nearest neighbor bandwidth
#' @param ker function taking as input a point at which to evaluate, an array of 
#' jumping times and a parameter (bandwidth or number of neighbors)
#' @param del number of neighbors for s-block cross validation (s = 2*del + 1)
#' @returns cross val of kernel estimator ker withcparameter h of hazard rate of jumping times Y 
crossval <- function(Y,h,ker,del){
  n<-length(Y)
  nh<-length(h)
  a = 0.1
  Y = sort(Y)
  n <- length(Y)
  haz_sqr<-function(t,Y,h){
    ker(t,Y,h)^2
  }
  CV<-numeric(nh)
  for (i in 1:nh) { 
    S<-0
    for (j in 1:n) {
      deb = max(1,j-del)
      fin = min(j+del,n)
      S<-S+n*ker(Y[j],Y[-(deb:fin)],h[i])/((n-j+1)) 
    }
    b=max(Y)+h[i]
    trap = trapezoid(haz_sqr,a,b,1000,Y,h[i])
    CV[i]<-trap-2/n*S
  }
  return(CV)
} #crossvalidation for bandwidth h

#' trapezoid function for quick integral computation
#' @param f function of three variables, a time, an array of jumping times and a parameter
#' (bandwidth or number of neighbors) typically f is a kernel estimator
#' @param a lower border of interval for t
#' @param b upper border of interval for t
#' @param mx number of intervals in direction t
#' @param Y array of jumping times
#' @param h bandwidth of kernel estimator OR number of neighbors for nearest neighbor bandwidth
#' @returns aproximation of integral of function f with parameters Y and h
#'  on interval [a,b] by trapezoid method
trapezoid<-function(f,a,b,mx,Y,h){
  hx<-(b-a)/mx
  x<-seq(from=a,to=b,by=hx)
  nx<-length(x)
  H<-sapply(x, function(t)(f(t,Y,h)))
  M<-rep(1,nx)
  M[c(1,nx)]<-1/2
  aprox_integral<-H%*%M*hx 
  return(aprox_integral)
}

#####################Minimax bandwidth estimation########


#'Function V_0 for local minimax bandwidth choice
#'@param t time of estimation 
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns the value of V_0 at t for the bandwidth b (estimation of the variance 
#'of the estimator) 
V_0 = function(t,b,B,max,m, kappa_0 = NULL, lambda = NULL){
  if(is.null(kappa_0)){kappa_0 = 0.03}
  if(is.null(lambda)){lambda = 4}
  return(kappa_0*max*log(m)/(m*sqrt(b))) 
}

#'Function A_0 for local minimax bandwidth choice
#'@param t time of estimation 
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param Ttest array of jumping times for hazard estimation
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns the value of A_0 at t for the bandwidth b (estimation of the bias
#'of the estimator) 
A_0 = function(t,b,B,max,m,Ttest,kappa_0 = NULL, lambda = NULL){
  A = c()
  for (i in 1:length(B)){
    inter = (ker_est_gamma_c(t,Ttest,B[i]) - ker_est_gamma_c(t,Ttest,max(b,B[i])))^2 - V_0(t,B[i],B, max,m,kappa_0 , lambda )
    A[i] = inter}
  return(max(A,0))   
}

#'Function V for global minimax bandwidth choice
#'@param Grid array on which to perform estimation
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns the value of V on Grid for the bandwidth b (estimation of the variance 
#'of the estimator) 
V = function(Grid,b,B,max,m, kappa_2 = NULL, lambda = NULL, epsilon = NULL) {
  if(is.null(kappa_2)){kappa_2 = 20}
  if(is.null(lambda)){lambda = 4}
  if(is.null(epsilon)){epsilon = 0.5}
  return((1+epsilon)^2*kappa_2*max/(m*sqrt(b)))
}

#'Function A for global minimax bandwidth choice
#'@param Grid array on which to perform estimation
#'@param b bandwidth 
#'@param B bandwidth set
#'@param max maximum of the hazard rate on the estimation interval
#'@param m sample size
#'@param Ttest array of jumping times for hazard estimation
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns the value of A on Grid for the bandwidth b (estimation of the bias
#'of the estimator) 
A = function(Grid,b,B,max,m,Ttest,kappa_2 = NULL, lambda = NULL, epsilon = NULL){
  A = c()
  for (i in 1:length(B)){
    integrand = function(t)((ker_est_gamma_c(t,Ttest,B[i]) - ker_est_gamma_c(t,Ttest,max(b,B[i])))^2)
    integrand = Vectorize(integrand)
    inter = integrate(integrand,0, max(Grid),subdivisions=2000)$value - V(Grid,B[i],B, max,m,kappa_2, lambda, epsilon )
    A[i] = inter}
  return(max(A,0))   
}


#'Simulates a sample for a given hazard rate and returns  the adaptive bandwidth 
#'estimator of the hazard rate and the chosen bandwidth for each point
#'@param Times the data sample
#'@param Grid the grid for the estimation 
#'@param seed a number to use as seed for randomness
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns a list with three arrays : K, the estimated hazard rate on Grid 
#'@p B the chosen bandwidth for each point of Grid and T, the simulated sample
minimax_pointwise <- function(Times,Grid,B,kappa_0 = NULL,lambda = NULL){

  m=length(Times)
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated
  #upper bound of hazard rate
  Kopt = c()
  Bopt=c()
  for (j in 1:length(Grid)){ 
    t = Grid[j]
    C = sapply(B, function(b)(A_0(t,b,B,max,m,Times,kappa_0, lambda) + V_0(t,b,B,max,m,kappa_0, lambda)))
    bopt = B[which.min(C)]
    Bopt[j] = bopt #local minimax bandwidth
    Kopt[j] = ker_est_gamma_c(t,Times,bopt) #estimator 
  }
  return(list('K' = Kopt, 'B'= Bopt,'T' =Times))
}


#'Simulates a sample for a given hazard rate and returns  the adaptive bandwidth 
#'estimator of the hazard rate and the chosen bandwidth for each point
#'@param Times a vector containing the observations
#'@param Grid the grid for the estimation 
#'@param B the bandwidth set to choose from
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns a list with three arrays : K, the estimated hazard rate on Grid 
#', B, the chosen bandwidth for each point of Grid and T, the simulated sample
minimax_global <- function(Times,Grid,B,kappa_2=NULL, lambda=NULL,epsilon=NULL){
  m = length(Times)
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated 
  #upper bound of hazard rate
  C = sapply(B, function(b)(A(Grid,b,B,max,m,Times,kappa_2, lambda,epsilon) + V(Grid,b,B,max,m,kappa_2, lambda,epsilon)))
  bopt = B[which.min(C)] #global minimax bandwidth
  return(list('B'= bopt,'T' =Times))
}




#'Returns a set of bandwidth depending for adaptive bandwidth choice
#'for a given sample size
#'@param m a sample size
#'@returns a set of bandwidths
Bandwidth_set <- function(m){
  B = c()
  B[1] = 400*(log(m)/m)^2
  i = 1
  while (i <= 20*log(m)){
    B[(i-1)/4+1] = log(m)^2*i/m
    i = i+4
  }
  B = B[!is.na(B)]
  B = B[sqrt(B) <= min(1,6/log(m))]
  B = B[sqrt(B) >= B[1]]
  return(B)
}



#'Returns a set of bandwidth depending for adaptive bandwidth choice
#'for a given sample size for the global bandwidth choice procedure
#'@param m a sample size
#'@returns a set of bandwidths
Bandwidth_set_global <- function(m){
  B = c()
  B[1] = 1/m^(2/3)
  i = 1
  while (i <= (10*sqrt(m))){
    B[(i-1)/10+1] = i/m^(2/3)
    i = i+10
  }
  B = B[!is.na(B)]
  B = B[sqrt(B) <= min(1,6/log(m))]
  B = B[sqrt(B) >= sqrt(B[1])]
  return(B)
}


#'Gives the adaptive bandwidth estimator of the hazard rate for a data sample on  
#'a Grid and the chosen bandwidth for each point of the Grid
#'@param Times the data sample 
#'@param Grid the grid for the estimation 
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@returns a list with two arrays : K, the estimated hazard rate on Grid 
#', B, the chosen bandwidth for each point of Grid 
minimax_pointwise_data <- function(Times, Grid,kappa_0 = NULL,lambda = NULL){
  m = length(Times)
  B = Bandwidth_set(m)
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated
  #upper bound of hazard rate
  Kopt = c()
  Bopt=c()
  for (j in 1:length(Grid)){ 
    t = Grid[j]
    C = sapply(B, function(b)(A_0(t,b,B,max,m,Times,kappa_0, lambda) + V_0(t,b,B,max,m,kappa_0, lambda)))
    bopt = B[which.min(C)]
    Bopt[j] = bopt #local minimax bandwidth
    Kopt[j] = ker_est_gamma_c(t,Times,bopt) #local estimator value
  }
  
  return(list('K' = Kopt, 'B'= Bopt))
}

#'Gives the adaptive bandwidth estimator of the hazard rate for a data sample on  
#'a Grid and the chosen global bandwidth 
#'@param Times the data sample 
#'@param Grid the grid for the estimation 
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@returns a list with two arrays : K, the estimated hazard rate on Grid 
#', B, the chosen bandwidth 
minimax_global_data <- function(Times, Grid,kappa_2=NULL, lambda=NULL,epsilon=NULL){
  m = length(Times)
  B = Bandwidth_set(m)
  max = max(sapply(Grid, function(t)(ker_est_gamma_c(t,Times,B[1])))) #estimated
  #upper bound of hazard rate
  Bopt=c()
  C = sapply(B, function(b)(A(Grid,b,B,max,m,Times,kappa_2, lambda,epsilon) + V(Grid,b,B,max,m,kappa_2, lambda,epsilon)))
  bopt = B[which.min(C)] #global minimax bandwidth
  Kopt = sapply(Grid, function(t)(ker_est_gamma_c(t,Times,bopt))) #global minmax bandwidth gamma estimator
  return(list('K' = Kopt, 'B'= bopt))
}

#'Computes an empirical approximation of the MISE on a grid and at 0 (or the first point of the Grid)
#'for the local minimax bandwidth choice and nearest neighbour bandwidth with the 
#'gamma kernel estimator
#'@param m sample size
#'@param reps number of repetitions to compute empirical MISE
#'@param param named list of parameters for the hazard rate 
#'@param hazard string of the form '{result= expression;}'
#'@param th_haz theoretical hazard evaluated on Grid
#'@param Grid a grid on which to conduct estimations
#'@param kappa_0 (optional) value of kappa_0
#'@param lambda (optional) value of lambda
#'@param nneigh number of neighbours for nearest neighbour bandwidth 
#'@returns a list with 2 elements mx and nn with the SE value on the Grid
#' for the reps repetitions for each method
MISE_approx_local <- function(m, reps ,param, hazard,th_haz,Grid,nneigh,kappa_0 = NULL,lambda = NULL){
  B1 = Bandwidth_set(m)
  m2 = length(Grid)
  MISE_mx = matrix(0L, nrow = reps, ncol = m2)
  MISE_nn = matrix(0L, nrow = reps, ncol = m2) 
  pop_init <- population(data.frame( birth = rep(0, m), 
                                     death = NA))
death <- mk_event_individual(type = 'death', 
                              intensity_code = hazard)

model <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = param)
  for (i in 1:reps){
    sim_out <- popsim(model = model, initial_population = pop_init,
                    events_bounds = c('death' = 1),
                    parameters = param,
                    time = 2000)

   data<- sim_out$population$death #simulated data
    RES = minimax_pointwise(data, Grid, B1,kappa_0, lambda)
    Kopt1 = RES$K #adaptive hazard estimator
    Bopt1 = RES$B #chosen bandwidth
    est_band1 = sapply(Grid, function(x)( ker_est_neigh(x,data,nneigh))) #nearest neighbour bandwidth gaussian kernel estimator  
    MISE_mx[i,] = (Kopt1-th_haz)^2
    MISE_nn[i,]  =(est_band1-th_haz)^2
  }
  return(list('mx' = MISE_mx,'nn'=MISE_nn))
}

#'Computes an empirical approximation of the MISE on a grid and at 0 (or the first point of the Grid)
#'for the global minimax bandwidth choice and CV bandwidth with the 
#'gamma kernel estimator
#'@param m sample size
#'@param reps number of repetitions to compute empirical MISE
#'@param param named list of parameters for the hazard rate 
#'@param hazard string of the form '{result= expression;}'
#'@param th_haz theoretical hazard evaluated on Grid
#'@param Grid a grid on which to conduct estimations
#'@param kappa_2 (optional) value of kappa_2
#'@param lambda (optional) value of lambda
#'@param epsilon (optional) value of epsilon
#'@param BW_CV set of bandwidths to choose from for cross-validation
#'@returns a list with 2 elements mx_g and cv with the SE value on the Grid
#'for the reps repetitions for each method
MISE_approx_global <- function(m, reps ,param, hazard,th_haz,Grid,BW_CV,kappa_2=NULL, lambda=NULL,epsilon=NULL){
  Bset = Bandwidth_set_global(m)
  m2 = length(Grid)
  MISE_mx_g = matrix(0L, nrow = reps, ncol = m2)
  MISE_cv = matrix(0L, nrow = reps, ncol = m2)
pop_init <- population(data.frame( birth = rep(0, m), 
                                     death = NA))
death <- mk_event_individual(type = 'death', 
                              intensity_code = hazard)

model <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = param)
  for (i in 1:reps){ 
   sim_out <- popsim(model = model, initial_population = pop_init, events_bounds = c('death' = 1),
                    parameters = param,
                    time = 2000)
    data<- sim_out$population$death #simulated data
    RES = minimax_global(data ,  Grid, Bset,kappa_2, lambda,epsilon)
    Bopt_g= RES$B #chosen bandwidth
    CV =crossval(sort(data),BW_CV,ker_est_gamma_c,0)
    Bopt_cv = BW_CV[which.min(CV)]
    res_adapt = sapply(Grid, function(t)(ker_est_gamma_c(t,data,Bopt_g)))
    res_CV = sapply(Grid, function(t)(ker_est_gamma_c(t,data,Bopt_cv)))
    MISE_mx_g[i,] =  (res_adapt-th_haz)^2
    MISE_cv[i,]  = (res_CV-th_haz)^2
  }
  return(list('mx_g' = MISE_mx_g,'cv'=MISE_cv))
}


#'Computes an empirical approximation of the MISE on a grid and at 0 (or the first point of the Grid)
#'for the CV bandwidth and nearest neighbour bandwidth with the gaussian kernel estimator
#'@param m sample size
#'@param reps number of repetitions to compute empirical MISE
#'@param param named list of parameters for the hazard rate 
#'@param hazard string of the form '{result= expression;}'
#'@param th_haz theoretical hazard evaluated on Grid
#'@param Grid a grid on which to conduct estimations
#'@param nneigh number of neighbours for nearest neighbour bandwidth 
#'@param BW_CV set of bandwidths to choose from for cross-validation
#'@returns a list with 2 elements  cv and nn with the SE value on the Grid
#' for the reps repetitions for each method
MISE_approx_gaussian <- function(m, reps ,param, hazard,th_haz,Grid,nneigh,BW_CV){
  with(as.list(param),{
    m2 = length(Grid)
    MISE_cv = matrix(0L, nrow = reps, ncol = m2)
    MISE_nn = matrix(0L, nrow = reps, ncol = m2)
    pop_init <- population(data.frame( birth = rep(0, m), 
                                     death = NA))
    pop_init <- population(data.frame( birth = rep(0, m), 
                                     death = NA))
    death <- mk_event_individual(type = 'death', 
                              intensity_code = hazard)

    model <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = param)
    for (i in 1:reps){
       sim_out <- popsim(model = model, initial_population = pop_init, events_bounds = c('death' = 1), parameters = param, time = 2000)
      Ttest = sim_out$population$death 
      CV =crossval(sort(Ttest),BW_CV,ker_estg,0)
      Bopt_cv = BW_CV[which.min(CV)]
      est_band1 = sapply(Grid, function(x)(ker_est_neigh(x,Ttest,nneigh))) #nearest neighbour bandwidth gaussian kernel estimator  
      res_CV = sapply(Grid, function(t)(ker_estg(t,Ttest,Bopt_cv ))) #gaussian kernel estimator         
      MISE_nn[i,]  =(est_band1-th_haz)^2
      MISE_cv[i,]  = (res_CV-th_haz)^2
    }
    return(list('cv' = MISE_cv,'nn'=MISE_nn))
  })
}
