library(dplyr) # for data-cleaning
library(metafor) # for multivariate meta-analytic working models 
library(clubSandwich)

#-------------------------------------------------------------------------------
# Prepare data
#-------------------------------------------------------------------------------

TSL15 <- 
  readRDS("Brief alcohol interventions/Tanner-Smith-Lipsey-2015-subset.rds") %>%
  # exclude observations missing control variables
  filter(!is.na(percoll), !is.na(attrition_all), !is.na(permale))

levels(TSL15$dv_cat) <- c("freq","heavy","quantity","peak","BAC","combined")
levels(TSL15$Ctype) <- c("Strawman","Attention","TAU","nontreatment")

TSL15_cent <- 
  TSL15 %>%
  group_by(studyid) %>%
  mutate(
    n_grps = length(unique(Ctype)),
    n_obs = n_tx_ob + n_ct_ob,
    postwks_c = pmin(postwks, 26) - 12,
    postwks_long = as.numeric(postwks > 26),
    percoll_c = percoll - 1,
    permale_c = permale - 0.5,
    attrition_c = attrition_all - median(attrition_all),
    study_dv = paste(studyid, dv_cat, sep = "-"),
    study_ctype = paste(studyid, Ctype, sep = "-")
  ) %>%
  # drop non-treatment conditions from studies with multiple control groups
  # so that Ctype is a between-study covariate
  filter(n_grps == 1 | Ctype != "nontreatment")


#-------------------------------------------------------------------------------
# dvcat model
#-------------------------------------------------------------------------------

covariates <- model.matrix(~ 0 + dv_cat + postwks_c + postwks_long + percoll_c + permale_c + attrition_c, data = TSL15_cent)
studyID <- as.factor(TSL15_cent$studyid)
subgroup <- as.factor(TSL15_cent$dv_cat)
sample_sizes <- tapply(TSL15_cent$n_obs, studyID, mean)

# Fit CMVE model to use as reference parameters

V_mat <- impute_covariance_matrix(TSL15_cent$V, 
                                  cluster = TSL15_cent$studyid, 
                                  r = 0.6,
                                  smooth_vi = TRUE)

dv_multivariate <- rma.mv(es ~ 0 + dv_cat + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                          V = V_mat, 
                          random = ~ dv_cat | studyid, struct = "HCS",
                          data = TSL15_cent, sparse = TRUE)

beta <- dv_multivariate$beta[,1]
tau_covariate <- TSL15_cent$dv_cat
tau_corr <- dv_multivariate$rho
tau_scalar <- sqrt(dv_multivariate$tau2 / dv_multivariate$tau2[3])
tau_scalar <- sqrt(seq(2.0, 0.5, length.out = length(tau_scalar)))

save(covariates, studyID, subgroup, sample_sizes, beta, tau_covariate, tau_corr, tau_scalar,
     file = "simulations/dv-model-data.Rdata")

#-------------------------------------------------------------------------------
# Ctype model
#-------------------------------------------------------------------------------

covariates <- model.matrix(~ 0 + Ctype + postwks_c + postwks_long + percoll_c + permale_c + attrition_c, data = TSL15_cent)
studyID <- as.factor(TSL15_cent$studyid)
subgroup <- as.factor(TSL15_cent$Ctype)
sample_sizes <- tapply(TSL15_cent$n_obs, studyID, mean)

# Fit SCE model to use as reference parameters

V_mat <- impute_covariance_matrix(TSL15_cent$V, 
                                  cluster = TSL15_cent$studyid, 
                                  r = 0.6,
                                  smooth_vi = TRUE)

dv_Ctype <- rma.mv(es ~ 0 + Ctype + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                          V = V_mat, 
                          random = ~ Ctype | studyid, struct = "DIAG",
                          data = TSL15_cent, sparse = TRUE)

beta <- dv_Ctype$beta[,1]
tau_covariate <- TSL15_cent$Ctype
tau_corr <- dv_Ctype$rho
tau_scalar <- sqrt(dv_Ctype$tau2 / dv_Ctype$tau2[4])
tau_scalar <- sqrt(seq(0.5, 2.0, length.out = length(tau_scalar)))

save(covariates, studyID, subgroup, sample_sizes, beta, tau_covariate, tau_corr, tau_scalar,
     file = "simulations/Ctype-model-data.Rdata")
