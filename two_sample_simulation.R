# Description -------------------------------------------------------------

### Two Sample (or called "Split Sample") Simulation
### Sample Size is controlled by "num.replicate"
### Effect Sizes are controllend by "coef.multi.u", "coef.multi.y"


# Setwd & Read Data -------------------------------------------------------
library(ivpack)
load("sample_gene_expression_and_snp.Rdata")

# Correct SE For TS2SLS ---------------------------------------------------

Correct_SE_TS2SLS <- function(n1,n2,lm1,lm2)
{
  #lm1 is the first stage, n1 is its sample size
  #lm2 is the second stage, n2 is its sample size
  const = 1 + n2/n1 * (lm2$coefficients[2])^2 * (summary(lm1)$sigma)^2 / 
    (summary(lm2)$sigma)^2
  return(const * vcov(lm2))
}


# Correct 2SRI var --------------------------------------------------------
# Function "calculate_A" is for calculating matrix "A"
# which is used for correcting standard error of GLM 2SRI

calculate_A <- function(Z,y,SNP,alpha,beta)
{
  
  #=====================
  #SNP is n by k
  #alpha is of length 1+k
  #beta is of length 3
  #y is of length n
  #=====================
  
  n = nrow(SNP)
  k = ncol(SNP)
  A_final = matrix(0,nrow = 3,ncol = 1+k)
  for(i in 1:n)
  {
    u_hat = y[i] - sum(c(1,SNP[i,])*alpha)
    v2 = c(1,y[i],u_hat)
    logit_p = sum(beta*v2)
    p = exp(logit_p)/(1+exp(logit_p))
    A_final = A_final - beta[3]* (Z[i]-p)^2 * (t(t(v2)) %*% t(c(1,SNP[i,])))
  }
  return(A_final)
}





# Calculating Simulation Coefs --------------------------------------------
# This section is aimed at generating coefficients for simulation,
# as described in the paper

# Linear Regression
colnames(cleaning_result)[4:33] = paste("SNP",1:30,sep="")
lm_formula = paste("y ~ ", paste("SNP",1:30,sep = "",collapse =" + "))
lm_formula = as.formula(lm_formula)
lm1 = lm(lm_formula,data = cleaning_result)

# Predict y and u
y_hat = predict(lm1)
u_hat = cleaning_result$y - y_hat
sd.uhat = sd(u_hat)

# Logistic Regression Only u_hat
glm1 = glm(cleaning_result$Status ~ u_hat, family = binomial(link = "logit"))

# Logistic Regression with y and u_hat
glm2 = glm(cleaning_result$Status ~ cleaning_result$y + u_hat,
           family = binomial(link = "logit"))



# SIMULATION ----------------------------------------------------------------
# Function for simulation

simulation_split <- function(coef.multi.u,coef.multi.y,coef.multi.const,
                             num.replicate,sd.uhat,glm1,glm2,lm1,
                             cleaning_result,null_case)
{
  set.seed(1)
  simulation_ind = sample(1: (num.replicate*712) )
  simulation_stage1_ind = simulation_ind[1 : (num.replicate*712/2)]
  simulation_stage2_ind = simulation_ind[(num.replicate*712/2 + 1) : 
                                           (num.replicate*712)]
  
  # SPLIT -------------------------------------------------------------------
  result = NULL
  set.seed(1)
  for(num_simulation in 1:1000)
  {
    # Generate u and y
    u = rnorm(num.replicate*712, mean = 0 , sd = sd.uhat)
    SNP = NULL
    for(i in 1:num.replicate)
    {
      SNP = rbind(SNP, cleaning_result[,4:33])
    }
    y = as.matrix(cbind(1,SNP[,c(5,15,25)])) %*% 
      t(t(lm1$coefficients[c(1,6,16,26)] )) + u
    
    # Generate p
    if(null_case)
    {
      logit.p = coef.multi.const * glm1$coefficients[1] + 
        coef.multi.u * glm1$coefficients[2] * u
    } else {
      logit.p = coef.multi.const * glm2$coefficients[1] + 
        coef.multi.y * glm2$coefficients[2] * y + 
        coef.multi.u * glm2$coefficients[3] * u
    }

    p = exp(logit.p) / (1 + exp(logit.p))
    
    # Generate Z
    Z = rbinom(num.replicate*712,size = 1, prob = p)
    
    # Merge Data
    simulation_data = cbind(Z,y,SNP)
    simulation_data_stage1 = simulation_data[simulation_stage1_ind,]
    simulation_data_stage2 = simulation_data[simulation_stage2_ind,]
    
    # First Stage
    SNP_stage1_ind = order(abs(cor(y[simulation_stage1_ind],SNP[simulation_stage1_ind,])), 
                           decreasing = T)[1:10]
    SNP_stage1 = SNP[simulation_stage1_ind,SNP_stage1_ind]
    SNP_stage1 = as.matrix(SNP_stage1)
    stage1_lm_formula = paste("y ~ ", paste("SNP",SNP_stage1_ind,sep = "",collapse =" + "))
    stage1_lm_formula = as.formula(stage1_lm_formula)
    lm_stage1 = lm(stage1_lm_formula,data = simulation_data_stage1)
    
    # Predict
    y_hat = predict(lm_stage1,newdata = simulation_data_stage2)
    u_hat = y[simulation_stage2_ind] - y_hat
    Z = Z[simulation_stage2_ind]
    y = y[simulation_stage2_ind]
    u = u[simulation_stage2_ind]
    
    # Second Stage
    glm_2sps = glm(Z ~ y_hat, family = binomial(link = "logit"))
    glm_2sri = glm(Z ~ y + u_hat, family = binomial(link = "logit"))
    glm_2sri_2 = glm(Z ~ y_hat + u_hat, family = binomial(link = "logit"))
    
    lm_2sps = lm(Z ~ y_hat)
    lm_2sri = lm(Z ~ y + u_hat)
    lm_2sri_2 = lm(Z ~ y_hat + u_hat)
    
    # Corrected lm_2sps
    #iv.model = ivreg( Z ~ y | SNP_stage1)
    #invisible(capture.output( corrected_se <- (robust.se(iv.model))[2,2]))
   
    
    # Correct se for 2SRI
    A = calculate_A(Z,y,SNP_stage1,lm_stage1$coefficients,glm_2sri$coefficients)
    v_alpha = vcov(lm_stage1)
    v_beta = vcov(glm_2sri)
    var_corrected = v_beta %*% A %*% v_alpha %*% t(A) %*% v_beta + v_beta
    se_corrected_2sri = sqrt(var_corrected[2,2])
    
    # Correct se for 2sps lm
    correct_2sps_se = Correct_SE_TS2SLS(num.replicate*712/2,num.replicate*712/2,lm_stage1,lm_2sps)
    
    # Oracle
    glm_oracle = glm(Z ~ y + u, family = binomial(link = "logit"))
    oracle = summary(glm_oracle)$coefficients[2,1:2]
    
    # Keep result
    result = rbind(result, 
                   c(summary(glm_2sps)$coefficients[2,1:2],
                     summary(glm_2sri)$coefficients[2,1:2],
                     se_corrected_2sri,
                     summary(glm_2sri_2)$coefficients[2,1:2],
                     summary(lm_2sps)$coefficients[2,1:2],
                     summary(lm_2sri)$coefficients[2,1:2],
                     summary(lm_2sri_2)$coefficients[2,1:2],
                     #iv.model$coefficients[2],
                     sqrt(correct_2sps_se[2,2]),
                     oracle
                   )
    )
  }
  return(result)
}


# DO SIMULATION -----------------------------------------------------------
# Here you can change "coef.multi.u", "coef.multi.y" 

coef.multi.const = 1

for(coef.multi.u in c(1))
{
  for(coef.multi.y in c(0))
  {
    for(num.replicate in c(1,3,5,7,9))
    {
      null_case = (coef.multi.y == 0)
      
      result1 = simulation_split(coef.multi.u,coef.multi.y,coef.multi.const,
                                 num.replicate,sd.uhat,glm1,glm2,lm1,
                                 cleaning_result,null_case)
      name_result = paste("u",coef.multi.u,"y",coef.multi.y,
                          "rep",num.replicate,"_split.txt",sep="")
      
      write.table(result1, name_result,row.names = F,
                  quote = F)
      
    }
  }
}

