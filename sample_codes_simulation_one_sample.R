# Description -------------------------------------------------------------

### One Sample (or called "Whole Sample") Simulation
### Sample Size is controlled by "num.replicate"
### Effect Sizes are controllend by "coef.multi.u", "coef.multi.y"


# Setwd & Read Data -------------------------------------------------------
library(ivpack)

load("sample_gene_expression_and_snp.Rdata")


# correct 2sri var --------------------------------------------------------

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



# SIMULATION WHOLE --------------------------------------------------------
simulation_whole <- function(coef.multi.u,coef.multi.y,coef.multi.const,
                             num.replicate,sd.uhat,glm1,glm2,lm1,
                             cleaning_result,null_case)
{
  # WHOLE -------------------------------------------------------------------
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
    
    # First Stage
    SNP_stage1_ind = order(abs(cor(y,SNP)), decreasing = T)[1:10]
    SNP_stage1 = SNP[,SNP_stage1_ind]
    SNP_stage1 = as.matrix(SNP_stage1)
    lm_stage1 = lm(y ~ SNP_stage1)
    
    # Predict
    y_hat = predict(lm_stage1)
    u_hat = y - y_hat
    
    # Second Stage
    glm_2sps = glm(Z ~ y_hat, family = binomial(link = "logit"))
    glm_2sri = glm(Z ~ y + u_hat, family = binomial(link = "logit"))
    glm_2sri_2 = glm(Z ~ y_hat + u_hat, family = binomial(link = "logit"))
    
    lm_2sps = lm(Z ~ y_hat)
    lm_2sri = lm(Z ~ y + u_hat)
    lm_2sri_2 = lm(Z ~ y_hat + u_hat)
    
    # Corrected lm_2sps
    iv.model = ivreg( Z ~ y | SNP_stage1)
    invisible(capture.output( corrected_se <- (robust.se(iv.model))[2,2]))
   
    
    # Correct se for 2SRI
    A = calculate_A(Z,y,SNP_stage1,lm_stage1$coefficients,glm_2sri$coefficients)
    v_alpha = vcov(lm_stage1)
    v_beta = vcov(glm_2sri)
    var_corrected = v_beta %*% A %*% v_alpha %*% t(A) %*% v_beta + v_beta
    se_corrected_2sri = sqrt(var_corrected[2,2])
    
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
                     iv.model$coefficients[2],
                     corrected_se,
                     oracle
                   )
    )
  }
}


# DO SIMULATION -----------------------------------------------------------
coef.multi.const = 1

for(coef.multi.u in c(1))
{
  for(coef.multi.y in c(0))
  {
    for(num.replicate in c(1,3,5,7,9))
    {
      null_case = (coef.multi.y == 0)
      
      result1 = simulation_whole(coef.multi.u,coef.multi.y,coef.multi.const,
                                 num.replicate,sd.uhat,glm1,glm2,lm1,
                                 cleaning_result,null_case)
      name_result = paste("u",coef.multi.u,"y",coef.multi.y,
                          "rep",num.replicate,"_whole.txt",sep="")
      
      write.table(result1, name_result,row.names = F,
                  quote = F)
      
    }
  }
}

