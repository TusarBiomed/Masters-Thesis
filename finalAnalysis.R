
load(file = "FinalData.RData")

library(dplyr)
library(ggplot2)

## Select variables for analysis; 
selectedVar <- c("in68genes", "syncount", "syncount.group", "chr9_139395259_G_A", "chr13_28609651_C_G")

## Variables on complete remission of patients;
## Statistical approach: Logistic regression

glm_fits <- lapply(selectedVar, function(e){
  tmp <- data %>% select(e, casesOfremission, ELN_2017) 
    tmp <- tmp %>% mutate_at(c(e, "casesOfremission"), function(e){ as.numeric(as.character(e)) })
  # }
  fit <- glm(casesOfremission ~ get(e) + ELN_2017, data = tmp, family = "binomial")
  fit_coe <- summary(fit)$coefficients
  fit_coe2 <- exp(cbind(OR = coef(fit), confint(fit)))
  res <- data.frame(Variables = e, OR = fit_coe2["get(e)", "OR"], lower = fit_coe2["get(e)", "2.5 %"], upper = fit_coe2["get(e)", "97.5 %"], p = fit_coe["get(e)", "Pr(>|z|)"], pheno = "CR")
  return(res)
}
  )

res <- bind_rows(glm_fits)
res_cr <- res

## Plot the significance level of each variables on 'Complete remission' of patients.
ggplot(data = res_cr, aes(x = OR, y = Variables)) +
  geom_point(aes(color = ifelse(p < 0.05, TRUE, FALSE))) +
  geom_errorbar(aes(xmin = lower, xmax = upper, color = ifelse(p < 0.05, TRUE, FALSE))) +
  geom_vline(xintercept = 1, lty = 6) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "black"), name = 'p-value < 0.05') +
  scale_x_log10() +
  ggtitle("95% CI of Variables to predict the complete remission of patients") + xlab("Confidence intervals") + ylab("Variables adjusting with ELN 2017")



##############
## Survival analysis

dataOS <- data %>% filter(!is.na(dead))

library(survival)

## Survival analysis by the Cox-proportional hazard regression;

Coxph_res <- lapply(selectedVar, function(e){
  tmp <- dataOS %>% select(e, overallSurvival, dead, ELN_2017) 
  tmp <- tmp %>% mutate_at(c(e), function(e){ as.numeric(as.character(e)) })
  fit <- coxph(Surv(overallSurvival, dead) ~ get(e) + strata(ELN_2017), data = tmp)
  fit_coef <- summary(fit)$coefficients
  fit_conf <- summary(fit)$conf.int
  res <- data.frame(Variables = e, HR = fit_conf['get(e)', 'exp(coef)'], lower = fit_conf['get(e)', 'lower .95'], upper = fit_conf['get(e)', 'upper .95'], p = fit_coef['get(e)', 'Pr(>|z|)'], pheno = 'OS')
  return(res)
})
res_os <- bind_rows(Coxph_res)


## Plot the significance level of each variables on survival analysis.
ggplot(data = res_os, aes(x = HR, y = Variables)) +
  geom_point(aes(color = ifelse(p < 0.05, TRUE, FALSE))) +
  geom_errorbar(aes(xmin = lower, xmax = upper, color = ifelse(p < 0.05, TRUE, FALSE))) + geom_vline(xintercept = 1, lty = 6) +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "black"), name = 'p-value < 0.05') +
  scale_x_log10() +
  ggtitle("95% CI of Variables in CPH regression") + xlab("Confidence interval") + ylab("Variables adjusting with stratified ELN 2017")


library(survminer)

# Unified 'major' and 'medium' allele frequency.

dataOS <- dataOS %>% mutate_at(2:28, function(e){
  ifelse(e == 2, 1, ifelse(e == 1, 1, 0))
})


# The Kaplan-Meier survival curve of the variables.
for(e in selectedVar) {
  if (e == "syncount" | "syncount.group") {
    next
  }
  tmp <- dataOS %>% select(e, overallSurvival, dead, ELN_2017) 
  tmp <- tmp %>% mutate_at(c(e), function(e){ as.numeric(as.character(e)) })
  fit <- survfit(Surv(overallSurvival, dead) ~ get(e), data = tmp)
  p <- ggsurvplot(fit1, data = dataOS,
                  surv.median.line = "hv", 
                  legend.title = "Variable",
                  pval = TRUE, conf.int = TRUE,  risk.table = TRUE,
                  tables.height = 0.2, tables.theme = theme_cleantable(), 
                  palette = c("blue", "red"),  ggtheme =theme_bw())
  ggsave(
    filename =paste0('OS_', e, '.pdf'),
    plot = print(p, newpage = FALSE),
    device = 'pdf',
    path = 'data/output',
    width = 8,
    height = 9
  )
}


## The Kaplan-Meier survival curve of each mutations-
reqsnps <- grep("^chr", names(dataOS), value = T)

for(e in reqsnps) {
  tmp <- dataOS %>% select(e, overallSurvival, dead, ELN_2017) 
  tmp <- tmp %>% mutate_at(c(e), function(e){ as.numeric(as.character(e)) })
  fit <- survfit(Surv(overallSurvival, dead) ~ get(e), data = tmp)
  p <- ggsurvplot(fit1, data = dataOS,
                  surv.median.line = "hv", 
                  legend.title = "Variable",
                  # Add p-value and tervals
                  pval = FALSE, conf.int = TRUE,  risk.table = TRUE,
                  tables.height = 0.2, tables.theme = theme_cleantable(), 
                  palette = c("blue", "red"),  ggtheme =theme_bw())
  ggsave(
    filename =paste0('OS_', e, '.pdf'),
    plot = print(p, newpage = FALSE),
    device = 'pdf',
    path = 'data/output',
    width = 8,
    height = 9
  )
}
