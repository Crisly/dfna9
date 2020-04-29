## -----------------------------------------------------------------------------
library(ggplot2)
library(drc)
library(sjPlot)
library("readxl")
library("nlme")
library(lme4)
library(knitr)
library(forcats)
library(tidyr)
library(dplyr)


## -----------------------------------------------------------------------------
data_raw <- read_excel("../data/raw_data/database_20-04-2020.xlsx")
data_raw$group = factor(data_raw$Domain)
#leave out data with only n=1 dataset per domain/certain unpublished data.
data_subset <- subset(data_raw, Smits=='no' & Domainrec!=1 & group!="Ivd1")  
data <- subset(data_subset, select=c('pid', 'group', 'Leeftijd', 'PTA54ADS'))
data <- droplevels(data)
#data_raw %>% group_by(group) %>% summarize(count=n())
#head(data)
#check for NaNs in PTA and Leeftijd, should be 0.
nrow(data[is.na(data$PTA54ADS) | is.na(data$Leeftijd), ])
#number (n) of counts (i.e. audiograms) per subject (pid)
summarytable <- data %>% count(group, pid)


## -----------------------------------------------------------------------------
num_meas_per_id <- aggregate(PTA54ADS ~ pid , data, function(x) length(unique(x)))
t1 <- table(num_meas_per_id$PTA54ADS)


## -----------------------------------------------------------------------------
kable(t1, caption = "Table 1. The number of subjects that each have n audiograms", col.names = c("# audiograms","# subjects"))


## -----------------------------------------------------------------------------
ggplot(data=summarytable, aes(x=n,fill=group)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~group) +
  xlab("# of measurements") +
  ylab("Frequency (# patients)") 
  dev.print(pdf, '../results/histogram_number_meas_pid.pdf')


## -----------------------------------------------------------------------------
ggplot(data, aes(x = Leeftijd, y = PTA54ADS, group = pid, color = group)) + 
  #xlab("Age (years)") +
  #ylab("PTA (dB HL)") +
  geom_point(aes(colour = factor(group))) +
  geom_line(data=data, size=1, alpha = .3) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light() 
  dev.print(pdf, '../results/pta_age_pid_groups.pdf')


## -----------------------------------------------------------------------------
lccl = subset(data, group=="LCCL")
ggplot(lccl, aes(x=Leeftijd, y=PTA54ADS), group=pid, color=group) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  geom_point(aes(colour = factor(group))) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  theme_light()
  dev.print(pdf, '../results/pta_age_pid_lccl.pdf')


## -----------------------------------------------------------------------------
lin_fit <- nls(PTA54ADS ~ a*Leeftijd + b, lccl)
summary(lin_fit)

## -----------------------------------------------------------------------------
nls_fit <- nls(PTA54ADS ~ a*Leeftijd^b, lccl, start = list(a = 0.05, b = 1.5))
summary(nls_fit)

## -----------------------------------------------------------------------------
startvec <- c(Asym = 120, xmid = 50, scal = 15)
nls_logis<- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal),
                 data=lccl,
                 start = startvec)
summary(nls_logis)


## -----------------------------------------------------------------------------
anova(nls_fit,nls_logis)


## -----------------------------------------------------------------------------
newdat = expand.grid(Leeftijd = seq(0, 100, by = 1))
newdat$lin <-predict(lin_fit, newdata = newdat)
newdat$pta_logistic <- predict(nls_logis,newdata = newdat)
newdat$pta_power <- predict(nls_fit,newdata = newdat)
#newdat
ggplot(lccl, aes(x=Leeftijd, y=PTA54ADS)) +
  geom_point(aes(colour = factor(group))) +
  geom_line(data = newdat, aes(y = lin), size = 1, col='blue') +
  geom_line(data = newdat, aes(y = pta_logistic), size = 1) +
  geom_line(data = newdat, aes(y = pta_power), size = 1, col='red') +
 
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light() 
  dev.print(pdf, '../results/pta_age_lccl_fits.pdf')



## -----------------------------------------------------------------------------
fit0 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal), data=data)
summary(fit0)
coef(fit0)


## -----------------------------------------------------------------------------
# https://stats.stackexchange.com/questions/27273/how-do-i-fit-a-nonlinear-mixed-effects-model-for-repeated-measures-data-using-nl
# https://stats.stackexchange.com/questions/316801/how-to-compare-logistic-regression-curves
fit1 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid[group], scal), 
            data=data,
            start=list(Asym=rep(120,1), xmid=rep(50,3), scal=rep(15,1)))
summary(fit1)


## -----------------------------------------------------------------------------
fit2 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid[group], scal[group]), 
            data=data,
            start=list(Asym=rep(120,1), xmid=rep(50,3), scal=rep(15,3)))
summary(fit2)

## -----------------------------------------------------------------------------
plot(profile(fit2))


## -----------------------------------------------------------------------------
fit3 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym[group], xmid[group], scal[group]), 
            data=data,
            start=list(Asym=rep(120,3), xmid=rep(50,3), scal=rep(15,3)))
summary(fit3)


## -----------------------------------------------------------------------------
anova(fit0,fit1,fit2,fit3)


## -----------------------------------------------------------------------------
newdat = expand.grid(Leeftijd = seq(0, 100, by = 1), group = c("LCCL","vWFA1","vWFA2"))
newdat$fit <- predict(fit2, newdata = newdat)

## -----------------------------------------------------------------------------

ggplot(data, aes(x=Leeftijd, y=PTA54ADS, group = pid, color = group)) +
  geom_point(aes(colour = factor(group)),alpha = .4) +
  geom_line(data=data, size=1, alpha = .2) +
  geom_line(data = newdat, aes(y = fit, group = group, colour = factor(group)), size = 1) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light() 
  dev.print(pdf, '../results/pta_age_pid_groups_fits.pdf')


## -----------------------------------------------------------------------------
#https://stackoverflow.com/questions/14439770/filter-rows-in-dataframe-by-number-of-rows-per-level-of-a-factor
pidlengths <- ave(as.numeric(data$pid), 
                     data$pid, FUN = length)
#df2 <- lccl[pidlengths > 5, ]
df2 <- data[pidlengths > 2, ]
with(df2, table(group))

## -----------------------------------------------------------------------------
models <- nlsList(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal) | pid, data = df2)


## -----------------------------------------------------------------------------
data_id <- subset(df2, pid=="147")
data_id
ggplot(data=data_id,  aes(x=Leeftijd, y=PTA54ADS)) +
  geom_point() +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light() 


## -----------------------------------------------------------------------------
df2$Pred <-predict(models)
df2_na <- na.omit(df2)
with(df2_na, table(group, pid))


## -----------------------------------------------------------------------------
le <- unique(df2_na$pid)
#le
newdat = expand.grid(Leeftijd = seq(0, 100, by = 1), pid=le)
newdat$prednlm <-predict(models, newdata=newdat)


## ----fig.height = 7, fig.width = 8--------------------------------------------
#https://stackoverflow.com/questions/37122994/plotting-a-list-of-non-linear-regressions-with-ggplot
#https://aosmith.rbind.io/2018/11/16/plot-fitted-lines/
ip <- ggplot(data=df2_na,  aes(x=Leeftijd, y=PTA54ADS, colour = pid)) +
  geom_point() +
  geom_line(data=newdat,aes(y=prednlm)) +
  facet_wrap(~pid) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,25)) +
  scale_y_reverse(breaks=seq(0,120,40), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light()
ip
dev.print(pdf, '../results/pta_age_pid_lccl_ind_fits.pdf')


## -----------------------------------------------------------------------------
ggplot(data=df2_na,  aes(x=Leeftijd, y=PTA54ADS, colour = pid)) +
  geom_point() +
  geom_line(data=newdat,aes(y=prednlm, group = pid)) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_x_continuous(breaks=seq(0,100,20)) +
  scale_y_reverse(breaks=seq(-10,130,20), limits=c(130,-10)) +
  #scale_y_reverse(limits=c(130,-10)) +
  xlab("Age (years)") +
  ylab("PTA (dB HL)") +
  theme_light() 
  dev.print(pdf, '../results/pta_age_pid_lccl_ind_fits_overlay.pdf')


## -----------------------------------------------------------------------------
nm1 <- nlmer(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal) ~ Asym + xmid + scal | pid, df2_na, start = c(Asym = 100, xmid = 60, scal = 15), corr = FALSE)
summary(nm1)
coef(nm1)


## -----------------------------------------------------------------------------
require(lattice) 
qqmath(ranef(nm1, condVar=TRUE))



## -----------------------------------------------------------------------------
data_all <- subset(data_subset, select=c('pid', 'group', 'Leeftijd', '250.AD','500.AD', '1000.AD','2000.AD','4000.AD','8000.AD', '250.AS','500.AS', '1000.AS','2000.AS','4000.AS','8000.AS'))
head(data_all)


## -----------------------------------------------------------------------------
tidier <- data_all %>%
  gather(f, dB, -pid, -group, -Leeftijd)
data_all_l <- tidier %>%
  separate(f, into = c("frequency", "ear"), sep = "\\.")
#head(data_all_l)
data_all_l$frequency = factor(data_all_l$frequency)
data_all_l$ear= factor(data_all_l$ear)
data_all_l <- na.omit(data_all_l)
#head(data_all_l)


## -----------------------------------------------------------------------------
p <- data_all_l%>%
  mutate(frequency = fct_relevel(frequency,"250", "500", "1000", "2000", "4000", "8000")) %>%
  ggplot(aes(x = frequency, y = dB)) +
    #geom_bar(stat="identity") +
    #geom_histogram() +
    geom_violin() + 
    facet_wrap(~ group, ncol = 3) +
    #geom_point()
    xlab("Frequency (Hz)") +
    ylab("Hearing level (dB)") +
    scale_y_reverse(limits=c(130,-10)) 
    #theme_classic()
p


## -----------------------------------------------------------------------------
random_intercept <-lme(dB ~ frequency , 
            random = ~1|pid,   #p. 896
            method = "ML", 
            na.action = na.exclude, 
            control = list(opt="optim"),
            correlation = corAR1(),  #see p.897; timepoints are not equally spaced;use corCAR1 
            data = data_all_l)
summary(random_intercept)
anova(random_intercept)


## -----------------------------------------------------------------------------
timeRI <- update(random_intercept, .~. + Leeftijd)
summary(timeRI)


## -----------------------------------------------------------------------------
sessionInfo()


## ----code = readLines(knitr::purl("DFNA9_notebook.Rmd", documentation = 1)), echo = T, eval = F----
## NA

