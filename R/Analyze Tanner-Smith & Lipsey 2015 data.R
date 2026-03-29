library(tidyverse) # for data-cleaning
library(kableExtra) # for making formatted table
library(robumeta) # for original implementation of RVE
library(metafor) # for multivariate meta-analytic working models 
library(clubSandwich) # for RVE with metafor; requires version 0.5.1 or higher

TSL15 <- 
  readRDS("Brief alcohol interventions/Tanner-Smith-Lipsey-2015-subset.rds") %>%
  # exclude observations missing control variables
  filter(!is.na(percoll), !is.na(attrition_all), !is.na(permale))

# Number of studies, total number of effect sizes, and distribution of effect sizes

TSL15 %>%
  group_by(studyid) %>%
  summarise(es = n()) %>%
  summarise(
    studies = n(), 
    effects = sum(es),
    es_min = min(es),
    es_max = max(es),
    es_median = median(es),
    es_Q1 = quantile(es, .25),
    es_q3 = quantile(es, .75)
  )

# Center control variables

TSL15_cent <- 
  TSL15 %>%
  mutate(
    postwks_c = pmin(postwks, 26) - 12,
    postwks_long = as.numeric(postwks > 26),
    percoll_c = percoll - 1,
    permale_c = permale - 0.5,
    attrition_c = attrition_all - median(attrition_all),
    study_dv = paste(studyid, dv_cat, sep = "-"),
    study_ctype = paste(studyid, Ctype, sep = "-")
  )

# constant sampling correlation assumption
rho <- 0.6


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#
# Moderator analysis by type of alcohol consumption measure
#
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

# count studies and effects
N_dv_cat <- 
  TSL15_cent %>%
  group_by(dv_cat, studyid) %>%
  summarise(effects = n()) %>%
  summarise(
    studies = n(),
    effects = sum(effects)
  )

#------------------------------------------------------------------------------------------
# CE: Original RVE correlated effects working model

dv_robu <- robu(es ~ 0 + dv_cat + postwks_c + postwks_long + percoll_c + permale_c + attrition_c, 
                    var.eff.size = V, studynum = studyid, 
                    data = TSL15_cent, 
                    rho = rho)
dv_robu

# Robust F-test
Wald_dv_robu <- Wald_test(dv_robu, 
                          constraints = constrain_equal(1:6), 
                          vcov = "CR2")
Wald_dv_robu

#------------------------------------------------------------------------------------------
# CHE: multilevel random effects model with constant sampling correlation working model

# constant sampling correlation working model
V_mat <- impute_covariance_matrix(TSL15_cent$V, 
                                  cluster = TSL15_cent$studyid, 
                                  r = rho, 
                                  smooth_vi = TRUE)

# fit random effects working model in metafor
dv_multilevel <- rma.mv(es ~ 0 + dv_cat + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                        V = V_mat, 
                        random = ~ 1 | studyid / esid,
                        data = TSL15_cent, sparse = TRUE)

dv_multilevel # Note that this reports model-based (not robust) standard errors

# RVE standard errors
CI_dv_multilevel <- conf_int(dv_multilevel, vcov = "CR2")
CI_dv_multilevel

# Robust F-test
Wald_dv_multilevel <- Wald_test(dv_multilevel,
                                constraints = constrain_equal(1:6), 
                                vcov = "CR2")
Wald_dv_multilevel


#------------------------------------------------------------------------------------------
# SCE: subgroup random effects model with constant sampling correlation working model

# constant sampling correlation working model within dv_cat subgroups
V_dv_subgroup <- impute_covariance_matrix(TSL15_cent$V, 
                                          cluster = TSL15_cent$studyid, 
                                          r = rho,
                                          smooth_vi = TRUE,
                                          subgroup = TSL15_cent$dv_cat)

# fit random effects working model in metafor
dv_subgroup <- rma.mv(es ~ 0 + dv_cat + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                      V = V_dv_subgroup, 
                      random = list(~ dv_cat | studyid), struct = "DIAG",
                      data = TSL15_cent, sparse = TRUE)

dv_subgroup # Note that this reports model-based (not robust) standard errors

# RVE standard errors
CI_dv_subgroup <- conf_int(dv_subgroup, vcov = "CR2")
CI_dv_subgroup

# Robust F-test
Wald_dv_subgroup <- Wald_test(dv_subgroup, 
                              constraints = constrain_equal(1:6), 
                              vcov = "CR2")
Wald_dv_subgroup


#------------------------------------------------------------------------------------------
# Assemble Table of DV results

CI_dv_robu <- conf_int(dv_robu, vcov = "CR2")
vcomp_dv_robu <- data.frame(term = "tau", beta = sqrt(as.numeric(dv_robu$mod_info$tau.sq)))

vcomp_dv_multilevel <- data.frame(
  term = c("tau","omega"),
  beta = sqrt(dv_multilevel$sigma2)
)

vcomp_dv_subgroup <- data.frame(
  term = rownames(CI_dv_subgroup)[1:6],
  tau = sqrt(dv_subgroup$tau2)
)

dv_Wald_tests <- 
  list(
    robumeta = Wald_dv_robu,
    multilevel = Wald_dv_multilevel,
    subgroup = Wald_dv_subgroup
  ) %>%
  map(as_tibble, rownames = "term") %>%
  bind_rows(.id = "Model") %>%
  select(Model, term, Est = p_val) %>%
  mutate(term = "Wald test p-value")


dv_results <- 
  list(
    robumeta = bind_rows(as_tibble(CI_dv_robu, rownames = "term"), vcomp_dv_robu),
    multilevel = bind_rows(as_tibble(CI_dv_multilevel, rownames = "term"), vcomp_dv_multilevel),
    subgroup = left_join(as_tibble(CI_dv_subgroup, rownames = "term"), vcomp_dv_subgroup)
  )%>%
  bind_rows(.id = "Model") %>%
  select(Model, term, Est = beta, SE, tau) %>%
  filter(str_detect(term, "^dv_cat") | term %in% c("tau","omega")) %>%
  mutate(
    term = str_remove(term, "^dv_cat"),
    term = str_replace_all(term, "\\.", " ")
  ) %>%
  bind_rows(dv_Wald_tests)

dv_table <- 
  dv_results %>%
  mutate(
    Est_SE = if_else(is.na(SE), 
                     formatC(Est, digits = 3, format = "f"), 
                     paste0(formatC(Est, digits = 3, format = "f"), "\n[", formatC(SE, digits = 3, format = "f"), "]"))
  ) %>%
  pivot_wider(id_cols = term, names_from = Model, values_from = c(Est, SE, Est_SE, tau)) %>%
  left_join(N_dv_cat, by = c("term" = "dv_cat"))

options(knitr.kable.NA = " ")

dv_table %>%
  select(term, studies, effects, starts_with("Est_SE"), starts_with("tau")) %>%
  select(term, studies, effects, ends_with("robumeta"), ends_with("multilevel"), ends_with("subgroup")) %>%
  select(-tau_robumeta, -tau_multilevel) %>%
  kable(
    digits = 3, 
    escape = FALSE,
    col.names = c("Coef","Studies","Effect sizes","Est. [SE]","Est. [SE]","Est. [SE]","tau")
  ) %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, 
                     "Correlated effects" = 1, 
                     "Correlated hierarchical effects" = 1, 
                     "Sub-group correlated effects" = 2)) %>%
  save_kable(file = "Table-DV.html")


# SEs versus correlated effects working model

dv_table %>%
  select(term, starts_with("SE")) %>%
  mutate(
    multilevel_robu = SE_multilevel / SE_robumeta - 1,
    subgroup_robu = SE_subgroup / SE_robumeta - 1
  )

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#
# Moderator analysis by type of control group
#
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

# count studies and effects
N_Ctype <- 
  TSL15_cent %>%
  group_by(Ctype, studyid) %>%
  summarise(effects = n()) %>%
  summarise(
    studies = n(),
    effects = sum(effects)
  )

#------------------------------------------------------------------------------------------
# CE: Original RVE correlated effects working model

Ctype_robu <- robu(es ~ 0 + Ctype + postwks_c + postwks_long + percoll_c + permale_c + attrition_c, 
                var.eff.size = V, studynum = studyid, 
                data = TSL15_cent, 
                rho = rho)
Ctype_robu

# Robust F-test
Wald_Ctype_robu <- Wald_test(Ctype_robu, 
                             constraints = constrain_equal(1:4), 
                             vcov = "CR2")
Wald_Ctype_robu

#------------------------------------------------------------------------------------------
# CHE: multilevel random effects model with constant sampling correlation working model

# constant sampling correlation working model
V_mat <- impute_covariance_matrix(TSL15_cent$V, 
                                  cluster = TSL15_cent$studyid, 
                                  r = rho,
                                  smooth_vi = TRUE)

# fit random effects working model in metafor
Ctype_multilevel <- rma.mv(es ~ 0 + Ctype + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                            V = V_mat, 
                            random = ~ 1 | studyid / esid,
                            data = TSL15_cent, sparse = TRUE)

Ctype_multilevel # Note that this reports model-based (not robust) standard errors

# RVE standard errors
CI_Ctype_multilevel <- conf_int(Ctype_multilevel, vcov = "CR2")
CI_Ctype_multilevel

# Robust F-test
Wald_Ctype_multilevel <- Wald_test(Ctype_multilevel, 
                                   constraints = constrain_equal(1:4), 
                                   vcov = "CR2")
Wald_Ctype_multilevel


#------------------------------------------------------------------------------------------
# SCE: subgroup random effects model with constant sampling correlation working model

# constant sampling correlation working model within dv_cat subgroups
V_Ctype_subgroup <- impute_covariance_matrix(TSL15_cent$V, 
                                             cluster = TSL15_cent$studyid, 
                                             r = rho,
                                             smooth_vi = TRUE,
                                             subgroup = TSL15_cent$Ctype)

# fit random effects working model in metafor
Ctype_subgroup <- rma.mv(es ~ 0 + Ctype + postwks_c + postwks_long + percoll_c + permale_c + attrition_c,
                          V = V_Ctype_subgroup, 
                          random = list(~ Ctype | studyid), struct = "DIAG",
                          data = TSL15_cent, sparse = TRUE)

Ctype_subgroup # Note that this reports model-based (not robust) standard errors

# RVE standard errors
CI_Ctype_subgroup <- conf_int(Ctype_subgroup, vcov = "CR2")
CI_Ctype_subgroup

# Robust F-test
Wald_Ctype_subgroup <- Wald_test(Ctype_subgroup, 
                                 constraints = constrain_equal(1:4), 
                                 vcov = "CR2")
Wald_Ctype_subgroup


#------------------------------------------------------------------------------------------
# Assemble Table of control group type results

CI_Ctype_robu <- conf_int(Ctype_robu, vcov = "CR2")
vcomp_Ctype_robu <- data.frame(term = "tau", beta = sqrt(as.numeric(Ctype_robu$mod_info$tau.sq)))

vcomp_Ctype_multilevel <- data.frame(
  term = c("tau","omega"),
  beta = sqrt(Ctype_multilevel$sigma2)
)

vcomp_Ctype_subgroup <- data.frame(
  term = rownames(CI_Ctype_subgroup)[1:4],
  tau = sqrt(Ctype_subgroup$tau2)
)

Ctype_Wald_tests <- 
  list(
    robumeta = Wald_Ctype_robu,
    multilevel = Wald_Ctype_multilevel,
    subgroup = Wald_Ctype_subgroup
  ) %>%
  map(as_tibble, rownames = "term") %>%
  bind_rows(.id = "Model") %>%
  select(Model, term, Est = p_val) %>%
  mutate(term = "Wald test p-value")


Ctype_results <- 
  list(
    robumeta = bind_rows(as_tibble(CI_Ctype_robu, rownames = "term"), vcomp_Ctype_robu),
    multilevel = bind_rows(as_tibble(CI_Ctype_multilevel, rownames = "term"), vcomp_Ctype_multilevel),
    subgroup = left_join(as_tibble(CI_Ctype_subgroup, rownames = "term"), vcomp_Ctype_subgroup)
  )%>%
  bind_rows(.id = "Model") %>%
  select(Model, term, Est = beta, SE, tau) %>%
  filter(str_detect(term, "^Ctype") | term %in% c("tau","omega")) %>%
  mutate(
    term = str_remove(term, "^Ctype"),
    term = recode(term, "Attention.sham" = "Attention/sham"),
    term = str_replace_all(term, "\\.", " ")
  ) %>%
  bind_rows(Ctype_Wald_tests)

Ctype_table <- 
  Ctype_results %>%
  mutate(
    Est_SE = if_else(is.na(SE), 
                     formatC(Est, digits = 3, format = "f"), 
                     paste0(formatC(Est, digits = 3, format = "f"), "\n[", formatC(SE, digits = 3, format = "f"), "]"))
  ) %>%
  pivot_wider(id_cols = term, names_from = Model, values_from = c(Est, SE, Est_SE, tau)) %>%
  left_join(N_Ctype, by = c("term" = "Ctype"))

options(knitr.kable.NA = " ")

Ctype_table %>%
  select(term, studies, effects, starts_with("Est_SE"), starts_with("tau")) %>%
  select(term, studies, effects, ends_with("robumeta"), ends_with("multilevel"),
         ends_with("subgroup")) %>%
  select(-tau_robumeta, -tau_multilevel) %>%
  kable(
    digits = 3, 
    escape = FALSE,
    col.names = c("Coef","Studies","Effect sizes","Est. [SE]","Est. [SE]","Est. [SE]","tau")
  ) %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, "Correlated effects" = 1, "Correlated hierarchical effects" = 1, "Sub-group correlated effects" = 2)) %>%
  save_kable(file = "Table-Ctype.html")

# SEs versus correlated effects working model

Ctype_table %>%
  select(term, starts_with("SE")) %>%
  mutate(
    multilevel_robu = SE_multilevel / SE_robumeta - 1,
    subgroup_robu = SE_subgroup / SE_robumeta - 1
  )
