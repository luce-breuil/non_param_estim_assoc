library("ggplot2")
library('IBMPopSim')
library("viridis")
library('logKDE')
library('TeachingDemos')

source('Kernel_estimation_article.R')#R file with base functions

BW_CV = (1:150)/3 #set of bandwidths for cross validation

a1 = 7*10^(-3)
b1 = 3*10^(-2)
c1 = 7*10^(-2)
param1 = list('a1'=a1, 'b1' = b1, 'c1' = c1) #parameters as list
Grid = c(0.1,1:600)
th_haz = a1+b1*exp(-c1*Grid) #theoretical hazard
reps = 50 #number of repetitions for empirical MISE computation
hazard1 = '{result = a1+b1*exp(-c1*pow(I.age(t),1));}' #theoretical hazard rate expression 
m = 500
n=200
b3 = 0.1 #bandwidth for log estimator

pop_init <- population(data.frame( birth = rep(0, m), 
                                     death = NA))
death <- mk_event_individual(type = 'death', 
                              intensity_code = hazard1)
model_death <- mk_model(characteristics = get_characteristics(pop_init),
                          events = list(death),
                          parameters = param1)


    m2 = length(Grid)
    MISE_log = matrix(0L, nrow = reps, ncol = m2)
    for (i in 1:reps){
        print(i)
       sim_out <- popsim(model = model_death,
                    initial_population = pop_init,
                    events_bounds = c(
                      'death' = 1),
                    parameters = param1,
                    time = n*10)
      Ttest = sim_out$population$death

        est_log = sapply(Grid, function(x)(log_ker_dens(x,b3,Ttest))) #lognormal kernel estimator of density
        Surv_log =1- sapply(Grid, function(x)(est_surv(x,b3,Ttest))) # empirical survival function for lognormal kernel ratio estimation
        MISE_log[i,]  = (est_log/Surv_log-th_haz)^2
}

    # Construct the file name
    file_name <- paste0("MISE_log/MISE_log_m_", m, "_b3_", b3,"rep",reps,".csv")
    
    # Save the matrix
    write.csv(MISE_log, file = file_name, row.names = FALSE)
