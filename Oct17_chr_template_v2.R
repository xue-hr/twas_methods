library(plink2R)
library(ivpack)
# Correct SE For TS2SLS ---------------------------------------------------

Correct_SE_TS2SLS <- function(n1,n2,lm1,lm2)
{
  #lm1 is the first stage, n1 is its sample size
  #lm2 is the second stage, n2 is its sample size
  const = 1 + n2/n1 * (lm2$coefficients[2])^2 * (summary(lm1)$sigma)^2 / 
    (summary(lm2)$sigma)^2
  return(const * vcov(lm2))
}

#==========correct 2sri var
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
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### read data
### 1.gene_exp
gene_exp_raw = read.csv("Sep_data/ADNI_Gene_Expression_Profile.csv", header = FALSE, sep = ",", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")
gene_exp = gene_exp_raw[10:49395,3:747]
colnames(gene_exp) = c("Symbol",as.character(t(gene_exp_raw[3,4:747])))
row.names(gene_exp) = NULL
for(i in 1:745)
{
  gene_exp[,i] = as.character(gene_exp[,i])
}
for(i in 2:745)
{
  gene_exp[,i] = as.numeric(gene_exp[,i])
}
### 2.glist data
glist = read.table("Sep_data/glist-hg19.txt")
glist$V1 = as.character(glist$V1)
glist$V4 = as.character(glist$V4)
dup_glist = glist$V4[which(duplicated(glist$V4))]
dup_ind_glist = which(is.element(glist$V4,dup_glist))
glist = glist[-dup_ind_glist,]
### 3.case control data
case_control_raw = read.csv("Sep_data/keyADNItables.csv", header = FALSE, sep = ",", quote = "\"",
                            dec = ".", fill = TRUE, comment.char = "")
case_control = case_control_raw[,c(2,8)]
case_control = case_control[-1,]
case_control[,1] = as.character(case_control[,1])
case_control[,2] = as.character(case_control[,2])
colnames(case_control) = c("ID","Status")
case_control = case_control[!duplicated(case_control$ID),]

###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### For Chromosome INDEX_CHR
glist_sub = glist[glist$V1==INDEX_CHR,]
length(glist_sub$V4) #1266, no duplication
gene_list = intersect(gene_exp$Symbol,glist_sub$V4)
final_result_2SPS = NULL
final_result_2SRI = NULL
final_result_2SPS_sep = NULL
final_result_2SRI_sep = NULL
step1_result = NULL
old_dir = getwd()
setwd(paste(old_dir,"/PLINK",sep=""))
for( gene_ind_chr in 1:length(gene_list))
{
  ###### First Step
  gene_name = gene_list[gene_ind_chr]
  gene_exp_sub = gene_exp[gene_exp$Symbol == gene_name,]
  gene_exp_sub = sapply(gene_exp_sub[,2:745],median)
  ID = names(gene_exp_sub)
  y = gene_exp_sub
  names(y) = NULL
  gene_exp_sub = cbind(ID,y)
  gene_exp_sub = as.data.frame(gene_exp_sub)
  gene_exp_sub$ID = as.character(gene_exp_sub$ID)
  gene_exp_sub$y = as.numeric(as.character(gene_exp_sub$y))
  
  ###### Second Step
  window_size = 100 #in kb
  gene_ind = which(glist$V4 == gene_name)
  start = floor(glist[gene_ind,2]/1000) - window_size
  end = ceiling(glist[gene_ind,3]/1000) + window_size
  ###
  plink_command = paste("./plink --bfile chromosomes/GATK_chrINDEX_CHR --chr INDEX_CHR --from-kb ",start," --to-kb ",end," --make-bed --out ",gene_name)
  system(plink_command)
  ###
  genotype_windowsize = read_plink(gene_name)
  genotype_windowsize.bed = genotype_windowsize$bed
  genotype_windowsize.fam = genotype_windowsize$fam
  genotype_windowsize.bim = genotype_windowsize$bim
  ID = row.names(genotype_windowsize.bed)
  genotype_windowsize.bed = cbind(ID,genotype_windowsize.bed)
  genotype_windowsize.bed = as.data.frame(genotype_windowsize.bed)
  for(i in 2:ncol(genotype_windowsize.bed))
  {
    genotype_windowsize.bed[,i] = as.numeric(as.character(genotype_windowsize.bed[,i]))
  }
  genotype_windowsize.bed[,1] = as.character(genotype_windowsize.bed[,1])
  row.names(genotype_windowsize.bed) = NULL
  genotype_windowsize.bed$ID = substr(genotype_windowsize.bed$ID,3,20)
  
  #######Now after 2 steps, we finish data cleaning, merge data by ID
  gene_data = merge(gene_exp_sub,genotype_windowsize.bed, by = "ID")
  gene_data = merge(case_control,gene_data, by = "ID")
  gene_data$Status = as.numeric(gene_data$Status != "CN")
  
  ######Remove SNPs
  gene_data_profile = gene_data[,1:3]
  gene_data_snp = gene_data[,4:ncol(gene_data)]
  ###NA
  ind_na = which(is.na(colSums(gene_data_snp)))
  gene_data_snp = gene_data_snp[,-ind_na]
  ###MAF
  MAF = 0.05
  ind_MAF = which(colSums(gene_data_snp) > (2*nrow(gene_data)*MAF))
  gene_data_snp = gene_data_snp[,ind_MAF]
  ###
  gene_data_snp = as.data.frame(gene_data_snp)
  if(ncol(gene_data_snp)<2)
  {
    next()
  }
  ###Delete No Variations
  ind_novariation1 = which(colSums(gene_data_snp==1) == nrow(gene_data_snp))
  ind_novariation2 = which(colSums(gene_data_snp==2) == nrow(gene_data_snp))
  ind_novariation = c(ind_novariation1,ind_novariation2)
  if(length(ind_novariation)>0)
  {
    gene_data_snp = gene_data_snp[,-ind_novariation]
  }
  ###
  gene_data_snp = as.data.frame(gene_data_snp)
  if(ncol(gene_data_snp)<2)
  {
    next()
  }
  ###correlation
  cor_cutoff = 0.9
  cor_matrix = (abs(cor(gene_data_snp)) < cor_cutoff)^2
  cor_matrix = cor_matrix + diag(1,ncol(cor_matrix))
  i = 1
  while(i < nrow(cor_matrix) )
  {
    ind = which(cor_matrix[i,] == 1)
    cor_matrix = cor_matrix[ind,ind]
    i = i + 1
  }
  cor_name = colnames(cor_matrix)
  ind_cor = which(is.element(colnames(gene_data_snp),cor_name))
  gene_data_snp = gene_data_snp[,ind_cor]
  ###
  gene_data_snp = as.data.frame(gene_data_snp)
  if(ncol(gene_data_snp)<2)
  {
    next()
  }
  ###top 30 with maximum cor with y
  y_snp_cor = abs(cor(gene_data$y,gene_data_snp))
  y_snp_cor = order(y_snp_cor,decreasing = T)
  gene_data_snp = gene_data_snp[,y_snp_cor]
  gene_data_snp = gene_data_snp[,1:min(ncol(gene_data_snp),30)]
  ###
  if( (ncol(gene_data_snp)>0) & (ncol(gene_data_snp)<nrow(gene_data_profile)) )
  {
      
      lm_pre = lm(gene_data_profile$y ~ as.matrix(gene_data_snp))
      na_ind = which(is.na(lm_pre$coefficients))
      if(length(na_ind)>0)
      {
          gene_data_snp = gene_data_snp[,-(na_ind-1)]
      }
      
    ###run lm
    lm1 = lm(gene_data_profile$y ~ as.matrix(gene_data_snp))
    summary_lm1 = summary(lm1)
    step1_result = rbind(step1_result,
                         c(summary_lm1$r.squared,
                           summary_lm1$fstatistic,
                           summary_lm1$df
                           )
                         )
    ybar = predict(lm1)
    yresidual = gene_data$y - ybar
    ###run 2SPS
    glm_2SPS = glm(cbind(gene_data$Status,1-gene_data$Status)~ybar,family = binomial(link = "logit"))
    
    ###run 2SPS lm
    lm_2SPS = lm(gene_data$Status ~ ybar)
    
    ###run 2SRI
    glm_2SRI = glm(cbind(gene_data$Status,1-gene_data$Status)~gene_data$y + yresidual,family = binomial(link = "logit"))
    
    ###run 2SRI lm
    lm_2SRI = lm(gene_data$Status ~ gene_data$y + yresidual)
    
    #==========
    #adjust var of 2sri
    A = calculate_A(gene_data$Status,gene_data$y,as.matrix(gene_data_snp),
                    lm1$coefficients,glm_2SRI$coefficients)
    v_alpha = vcov(lm1)
    v_beta = vcov(glm_2SRI)
    var_corrected = v_beta %*% A %*% v_alpha %*% t(A) %*% v_beta + v_beta
    var_corrected_2sri = sqrt(var_corrected[2,2])
    #==========
    #adjust var for 2sps
    # Corrected lm_2sps
    iv.model = ivreg( gene_data$Status ~ gene_data$y | as.matrix(gene_data_snp))
    invisible(capture.output( corrected_se <- (robust.se(iv.model))[2,2]))
    #==========
    final_result_2SPS = rbind(final_result_2SPS,
                              c(gene_name,summary(glm_2SPS)$coefficients[2,1:2],
                                ncol(gene_data_snp),
                                summary(lm1)$r.squared,
                                summary(lm1)$adj.r.squared,
                                summary(lm_2SPS)$coefficients[2,1:2],
                                corrected_se
                                )
                              )
    
    final_result_2SRI = rbind(final_result_2SRI,
                              c(summary(glm_2SRI)$coefficients[2,1:2],
                                var_corrected_2sri,
                                summary(lm_2SRI)$coefficients[2,1:2]
                                )
                              )
  } else {
    final_result_2SPS = rbind(final_result_2SPS,
                              c(gene_name,0,0,ncol(gene_data_snp),
                                summary(lm1)$r.squared,
                                summary(lm1)$adj.r.squared,0,0,0
                                )
                              )
    final_result_2SRI = rbind(final_result_2SRI,
                              c(0,0,0,0,0))
  }
  
  ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
  ######split
  gene_data_final_snp = as.data.frame(gene_data_snp)
  gene_data_final_profile = as.data.frame(gene_data_profile)
  gene_data_final = cbind(gene_data_final_profile,gene_data_final_snp)
  gene_data_final = gene_data_final[,-c(1,2)]
  ### need to change name!
  colnames(gene_data_final) = c("y",paste("SNP",1:(ncol(gene_data_final)-1),sep = ""))
  ###Seperate
  set.seed(1)
  sample_size = nrow(gene_data_final)
  stage12_sample = sample(1:sample_size)
  stage1_ind = stage12_sample[1:(sample_size/2)]
  stage2_ind = stage12_sample[(sample_size/2+1):sample_size]
  gene_data_final_stage1 = gene_data_final[stage1_ind,]
  gene_data_final_stage2 = gene_data_final[stage2_ind,]
  lm_pre = lm(y ~ ., data = gene_data_final_stage1)
  na_ind = which(is.na(lm_pre$coefficients))
  if(length(na_ind)>0)
  {
    gene_data_final_stage1 = gene_data_final_stage1[,-na_ind]
    gene_data_final_stage2 = gene_data_final_stage2[,-na_ind]
    gene_data_snp = gene_data_snp[,-(na_ind-1)]
  }

  ###run lm
  if( (ncol(gene_data_snp)>0) & (ncol(gene_data_snp)<(nrow(gene_data_profile)/2)) )
  {
      lm1 = lm(y ~ ., data = gene_data_final_stage1)
      ybar = predict(lm1,newdata = gene_data_final_stage2)
      yresidual = gene_data_final_stage2$y - ybar
      ###run 2SPS
      glm_2SPS = glm(cbind(gene_data$Status,1-gene_data$Status)[stage2_ind,]~ybar,family = binomial(link = "logit"))
      
      ###run 2SPS lm
      lm_2SPS = lm(gene_data$Status[stage2_ind]~ybar)
      
      ###run 2SRI
      glm_2SRI = glm(cbind(gene_data$Status,1-gene_data$Status)[stage2_ind,]~gene_data_final_stage2$y + yresidual,family = binomial(link = "logit"))
      
      ###run lm 2SRI
      lm_2SRI = lm(gene_data$Status[stage2_ind]~gene_data_final_stage2$y + yresidual)
      
      
      #==========
      #adjust var of 2sri
      A = calculate_A(gene_data$Status[stage2_ind],gene_data_final_stage2$y,
                      as.matrix(gene_data_snp[stage2_ind,]),
                      lm1$coefficients,glm_2SRI$coefficients)
      v_alpha = vcov(lm1)
      v_beta = vcov(glm_2SRI)
      var_corrected = v_beta %*% A %*% v_alpha %*% t(A) %*% v_beta + v_beta
      var_corrected_2sri = sqrt(var_corrected[2,2])
      #==========
      # Correct se for 2sps lm
      correct_2sps_se = Correct_SE_TS2SLS(712/2,712/2,lm1,lm_2SPS)
      
      #==========
      final_result_2SPS_sep = rbind(final_result_2SPS_sep,
                                    c(gene_name,summary(glm_2SPS)$coefficients[2,1:2],
                                      ncol(gene_data_snp),
                                      summary(lm1)$r.squared,
                                      summary(lm1)$adj.r.squared,
                                      summary(lm_2SPS)$coefficients[2,1:2],
                                      sqrt(correct_2sps_se[2,2])
                                    )
      )
      
      final_result_2SRI_sep = rbind(final_result_2SRI_sep,
                                    c(summary(glm_2SRI)$coefficients[2,1:2],
                                      var_corrected_2sri,
                                      summary(lm_2SRI)$coefficients[2,1:2]
                                    )
      )      
      
  } else{
    final_result_2SPS_sep = rbind(final_result_2SPS_sep,
                              c(gene_name,0,0,ncol(gene_data_snp),
                                summary(lm1)$r.squared,
                                summary(lm1)$adj.r.squared,0,0,0
                              )
    )
    final_result_2SRI_sep = rbind(final_result_2SRI_sep,
                              c(0,0,0,0,0))
  }
  remove_command = paste("rm ",gene_name,".*",sep = "")
  system(remove_command)
}

colnames(final_result_2SPS) = c("Gene_Name","glm2sps_est","glm2sps_se",
                                "num_gene","r2","adjr2","lm2sps_est",
                                "lm2sps_se","corr")

colnames(final_result_2SPS_sep) = c("Gene_Name_sep","glm2sps_est_sep","glm2sps_se_sep",
                                    "num_gene_sep","r2_sep","adjr2_sep","lm2sps_est_sep",
                                    "lm2sps_se_sep","corr")

colnames(final_result_2SRI) = c("glm2sri_est","glm2sri_se","glm2sri_adjse",
                                "lm2sri_est","lm2sri_se")

colnames(final_result_2SRI_sep) = c("glm2sri_est_sep","glm2sri_se_sep","glm2sri_adjse_sep",
                                "lm2sri_est_sep","lm2sri_se_sep")


final_result = cbind(final_result_2SPS,
                     final_result_2SPS_sep,
                     final_result_2SRI,
                     final_result_2SRI_sep)
write.table(final_result, "Oct_17_chrINDEX_CHR_v2.txt",
            quote = F,
            row.names = F)

write.table(step1_result, "Oct_17_chrINDEX_CHR_v2_step1.txt",
            quote = F,
            row.names = F)
###36,37,68,217,221

