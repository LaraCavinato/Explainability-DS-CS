library("readxl")
library(ggplot2)
library(tidyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(cowplot)
library(ggpubr)
library(ggpmisc)
library("readxl")

library("survival")
library("survminer")

library(randomForest)
library(sirus)
library(caret)

setwd("~/Desktop/POLIMI/PhD/progetti/S2GC-master")

# ------------------ CLINCAL DATASET

hrh <- read_excel("DATA/humanitas_linfoma/summary_clinical_hum.xls")
int <- read_excel("DATA/int_linfoma/summary_clinical_int.xls")

# ------------------ SURVIVAL SUMMARY STATISTICS
km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ STADIO, data=hrh)
ggsurvplot(km_trt_fit)

km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ SINTOMI_B, data=hrh)
ggsurvplot(km_trt_fit)

km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ EXTRANODAL, data=hrh)
ggsurvplot(km_trt_fit)

km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ BONE, data=hrh)
ggsurvplot(km_trt_fit)

km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ RADIOTERAPIA, data=hrh)
ggsurvplot(km_trt_fit)

km_trt_fit <- survfit(Surv(REF_TIME, REF) ~ Pt_ID, data=hrh)
ggsurvplot(km_trt_fit)

# -------
fit_strata_hrh <- survfit( Surv(REF_TIME, REF) ~ BONE + RADIOTERAPIA + SINTOMI_B,
                 data = hrh )

ggsurv <- ggsurvplot(fit_strata_hrh, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())

ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(BONE ~ RADIOTERAPIA)

# -------
fit_strata_int <- survfit( Surv(REF_TIME, REF) ~ BONE + RADIOTERAPIA + SINTOMI_B,
                       data = int )

ggsurv <- ggsurvplot(fit_strata_int, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())

ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(BONE ~ RADIOTERAPIA)

# ------------------ EXPLAINABILITY

# RANDOM FOREST WITH SIRUS
X_train <- data.frame(hum_ds[,1:45])
y_train <- hum_ds[,48]
X_test <- data.frame(int_ds[,1:45])
model.cv <- sirus.cv(X_train,
                     y_train,
                     type = "auto",
                     nfold = 10,
                     ncv = 10,
                     num.rule.max = 25,
                     q = 10,
                     discrete.limit = 10,
                     num.trees.step = 1000,
                     alpha = 0.05,
                     mtry = NULL,
                     max.depth = 5,
                     num.trees = NULL,
                     num.threads = NULL,
                     replace = TRUE,
                     sample.fraction = NULL,
                     verbose = TRUE,
                     seed = NULL
)

plot.error <- sirus.plot.cv(model.cv)$error
plot(plot.error)

model <- sirus.fit(X_train,
                   y_train,
                   type = "classif",
                   num.rule = 10,
                   p0 = model.cv$p0.pred, # NULL
                   num.rule.max = 25,
                   q = 10,
                   discrete.limit = 10,
                   num.trees.step = 1000,
                   alpha = 0.05,
                   mtry = NULL,
                   max.depth = 2,
                   num.trees = NULL,
                   num.threads = NULL,
                   replace = TRUE,
                   sample.fraction = 0.632, # 1
                   verbose = TRUE,
                   seed = NULL
)

sirus.print(model)

sp <- ggplot(X_train, aes(x=dispersion_all, y=RADIOTERAPIA, color=factor(y_train))) + 
  geom_point(size=3)
sp + scale_color_manual(values=c("#FEC50F", "#56B4E9"))

sp <- ggplot(X_train, aes(x=GLZLM_ZLNU, y=GLRLM_RLNU, color=factor(y_train))) + 
  geom_point(size=3)
sp + scale_color_manual(values=c("#FEC50F", "#56B4E9"))

sp <- ggplot(X_train, aes(x=dispersion_all, y=Volume_max, color=factor(y_train))) + 
  geom_point(size=3)
sp + scale_color_manual(values=c("#FEC50F", "#56B4E9"))

label <- sirus.predict(model, X_train)
prob.train <- sum(y_train == 1) / length(y_train)
label <- factor(ifelse(label > prob.train, "1", "0"))
train.confusion <- confusionMatrix(label, factor(y_train), positive = NULL)

y_pred <- sirus.predict(model, X_test)

