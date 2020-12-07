#############################################################################################
library(tidyverse)

## Hyperparameters
#n                            # Population size
#sens                         # Sensitivity of PCR
#spec                         # Specificity of PCR
#prob_expo                    # Probability of exposure in the general population
#inf                          # Probability of disease prevalence after exposure
#allele_freq_inf              # Frequency of the allele affecting disease infection after exposure
#inf_effect                   # Log odds ratio after exposure with an additional allele 
#allele_freq_risk             # Frequency of the risk allele
#risk_effect                  # Log odds ratio of hospitalization after infection with an additional allele 
#baseline_risk                # Baseline risk of hospitalization with no risk alleles
#unrelated_symptoms           # Probability of showing symptoms but not exposed to COVID. Assumed random. 
#prob_tested_mild             # Probability of those infeced with COVID and with mild symptoms get tested

######## Generate the population
generate_pop <- function(n, sens, spec, inf, prob_expo, allele_freq_inf, allele_freq_risk, inf_effect, risk_effect, 
                         baseline_risk, unrelated_symptoms, prob_tested_mild){
  df <- data.frame(snp_risk=rbinom(n=n, size=2, prob=allele_freq_risk),
                   snp_inf=rbinom(n=n, size=2, prob=allele_freq_inf))
  
  # Define Exposure. Dominant Model
  df <- df %>%
    mutate(exposure=rbinom(n=n, size=1, prob=prob_expo))
  
  
  # Define disease status. Dominant Model
  df <- df %>%
    mutate(status = case_when(
      exposure==0 ~ as.integer(0),
      snp_inf==0 & exposure==1 ~ rbinom(n(), size=1, prob=inf),
      snp_inf==1 & exposure==1 ~ rbinom(n(), size=1, prob=plogis(qlogis(inf)+1*inf_effect)),
      snp_inf==2 & exposure==1 ~ rbinom(n(), size=1, prob=plogis(qlogis(inf)+2*inf_effect))
    )
    )
  
  
  # Define prob of hospitalization and actual hospitalization
  df <- df %>%
    mutate(severity = case_when(
      status==0 ~ 0, 
      status==1 ~ plogis(qlogis(baseline_risk)+snp_risk*risk_effect)
    )
    ) %>%
    mutate(hospitalization = rbinom(n=n, size=1, prob=severity))
  
  # Define those who has symptoms unrelated to COVID and got tested negative.
  # Assume random in the population
  df <- df %>%
    mutate(get_tested = case_when(
      status==0 ~ rbinom(n(), size=1, prob=unrelated_symptoms),
      hospitalization==1 ~ as.integer(1),
      status==1 & hospitalization==0 ~ rbinom(n(), size=1, prob=prob_tested_mild)
    )
    )
  
  ## Those who are hospitalized should have test result equal to one
  df <- df %>%
    mutate(test_result = case_when(
      status == 0 ~ rbinom(n(), size=1, prob=1-spec),
      status == 1 & hospitalization == 0 ~ rbinom(n(), size=1, prob=sens),
      hospitalization == 1 ~ as.integer(1)
    )
    )
  
  
  df
}

##########################################################################
## Cases: Tested Positive, Controls: Tested Negative.
## dim: Number of iterations
simulation <- function(dim=500, n_case=1000, n, sens, spec, inf, prob_expo, allele_freq_inf, allele_freq_risk, inf_effect, risk_effect, 
                       baseline_risk, unrelated_symptoms, prob_tested_mild){
  pval_1 <- array(dim=dim)
  
  for(i in 1:dim){
    df <- generate_pop(n, sens, spec, inf, prob_expo, allele_freq_inf, allele_freq_risk, inf_effect, risk_effect, 
                       baseline_risk, unrelated_symptoms, prob_tested_mild)
    sd1 <- df %>%
      filter(get_tested==1) %>%
      mutate(case_control=ifelse(test_result==1, 1, 0))
    
    sd1 <- rbind(sample_n(sd1[sd1$case_control==0,], n_case, replace=FALSE), sample_n(sd1[sd1$case_control==1,], n_case, replace=FALSE))
    pval_1[i] <- summary(glm(case_control ~ snp_inf, data=sd1, family=binomial(link="logit")))$coef[2,4]
    
  }
  
  df <- data.frame(value=c(sum(pval_1<5e-8)/dim))
  
  df
}


