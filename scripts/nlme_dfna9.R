library(ggplot2)
library(drc)
library(dplyr)
library(sjPlot)
library("readxl")
library("nlme")
library(lme4)

data <- read_excel("/Users/cplanting/Nextcloud/ProjectsScience/Hearing_Genes/DFNA9_SV/databasevoorpythonRid.xlsx")
datazjeroen <- subset(data, Smits == 'no',Domainrec!=1)
datazjeroenzIvd1 <-subset(datazjeroen, Domainrec !=1) #leave out data with only n=1 dataset per domain
nuttigedata<-subset(datazjeroenzIvd1,select=c('Leeftijd','Domainrec','PTA54ADS','Hz8000AD','Hz8000AS','Hz1000AD','Hz1000AS','Familynrrec','pid','Domain'))
nuttigedata$group=factor(nuttigedata$Domain)
#check for nans in PTA and Leeftijd
nrow(nuttigedata[is.na(nuttigedata$PTA54ADS) | is.na(nuttigedata$Leeftijd), ])

#todo: make fit on binned data

# non-linear fit of data with logistic function, self starting. See e.g.
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/SSlogis
# model_ss<-nls(PTA54ADS~SSlogis(Leeftijd,phi1,phi2,phi3),data=nuttigedata)
# summary(model_ss)

# https://stats.stackexchange.com/questions/316801/how-to-compare-logistic-regression-curves
# to do: scale PTA54ADS to range [0,1]
# glm_model<- glm(PTA54ADS ~ Leeftijd + Domain, family = binomial, data = nuttigedata)
# summary(glm_model)


# function fits with  DRC
df_lccl <- nuttigedata %>% filter(Domain == "LCCL")
df_wfa1 <- nuttigedata %>% filter(Domain == "vWFA1")
df_wfa2 <- nuttigedata %>% filter(Domain == "vWFA2")

# https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# I've tried the LL.2 and G.2 functions, see e.g. https://cran.r-project.org/web/packages/drc/drc.pdf
# for now on 1000 Hz to compare with Janssens de Varebeke

lccl <- drm(Hz1000AD ~ Leeftijd, data = df_lccl, fct = LL.2(upper=130))
summary(lccl)
lccl.fits <- expand.grid(l=seq(1,100, length=100))
pm_lccl <- predict(lccl, newdata=lccl.fits, interval="confidence") 
lccl.fits$p <- pm_lccl[,1]
lccl.fits$pmin <- pm_lccl[,2]
lccl.fits$pmax <- pm_lccl[,3]

wfa2 <- drm(Hz1000AD ~ Leeftijd, data = df_wfa2, fct = LL.2(upper=130))
summary(wfa2)
wfa2.fits <- expand.grid(l=seq(1,100, length=100))
pm_wfa2 <- predict(wfa2, newdata=wfa2.fits, interval="confidence") 
wfa2.fits$p <- pm_wfa2[,1]
wfa2.fits$pmin <- pm_wfa2[,2]
wfa2.fits$pmax <- pm_wfa2[,3]

ggplot(nuttigedata, aes(x = Leeftijd, y = Hz1000AD)) + 
  geom_point(aes(colour = factor(Domain))) +
  geom_ribbon(data=lccl.fits, aes(x=l, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill='red') +
  geom_line(data=lccl.fits, aes(x=l, y=p),color='red') +
  geom_ribbon(data=wfa2.fits, aes(x=l, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill='blue') +
  geom_line(data=wfa2.fits, aes(x=l, y=p),color='blue') +
  ylim(0, 125) +
  xlim(0,80) +
  #geom_line(data=df2, size=1, alpha = .5) +
  theme_classic()


# fit SSlogis to data to get some init coefs
# https://stats.stackexchange.com/questions/27273/how-do-i-fit-a-nonlinear-mixed-effects-model-for-repeated-measures-data-using-nl
# https://stats.stackexchange.com/questions/316801/how-to-compare-logistic-regression-curves

fit1 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal), data=nuttigedata)
summary(fit1)
coef(fit1)

# compare parameters across groups; first with xmid and slope/scale different across domains (fit2); later with all params fitted
fit2 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid[group], scal[group]), 
            data=nuttigedata,
            start=list(Asym=rep(120,1), xmid=rep(50,3), scal=rep(15,3)))
summary(fit2)
fit3 <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym[group], xmid[group], scal[group]), 
            data=nuttigedata,
            start=list(Asym=rep(120,3), xmid=rep(50,3), scal=rep(15,3)))
summary(fit3)

anova(fit1,fit2,fit3)
#apparently, giving each group a different value for Asym does not increase the Res.Sum Sq.

#good start for plots/ideas: https://asancpt.github.io/nlme/chapter-6.html
plot(fit2, group ~ resid(.), abline = 0)


#https://datascienceplus.com/second-step-with-non-linear-regression-adding-predictors/
z <- nls(PTA54ADS ~ a*Leeftijd^b+c, start = list(a=0.1, b=2, c=0), 
         data=nuttigedata)
summary(z)

plot(nuttigedata$Leeftijd, nuttigedata$PTA54ADS)
lines(predict(z))
xfit <- expand.grid(l=seq(1,100, length=100))
yfit <- predict(z,newdata=xfit)
line(xfit,yfit)

z.fits <- expand.grid(l=seq(1,100, length=100))
fit <- predict(z)
plot()
lines(predict(z))
z.fits <- expand.grid(l=seq(1,100, length=100))
z.fits$yfit <- predict(z,newdata = z_fits)

ggplot(nuttigedata, aes(x = Leeftijd, y = Hz1000AD)) + 
  geom_point(aes(colour = factor(Domain))) +
  geom_line(predict(z))

pm_all_data <- predict(nuttigedata, newdata=z.fits, interval="confidence") 
z.fits$p <- pm_wfa2[,1]
z.fits$pmin <- pm_wfa2[,2]
z.fits$pmax <- pm_wfa2[,3]
# https://rdrr.io/cran/lme4/man/nlmer.html
startvec <- c(Asym = 120, xmid = 50, scal = 15)



nl_null <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal),
                 data=nuttigedata,
                 start = startvec)
summary(nl_null)

domain.list <- nlsList(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal) | pid, 
                       data = df_lccl,
                       na.action = na.exclude,
                       start = startvec,
                       control = list(maxiter = 1000))
summary(domain.list)


dom.list <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal), 
                    data = nuttigedata,
                    na.action = na.exclude,
                    start = startvec,
                    fixed = list(Asym ~ group, xmid ~ group))
summary(nm1)
# Likelihood ratio test
anova(fit,fit2, test = "Chisq")


#label pids with >1 datapoints
#https://stackoverflow.com/questions/14439770/filter-rows-in-dataframe-by-number-of-rows-per-level-of-a-factor
cts <- count(nuttigedata, vars = "pid")
keep <- as.character(subset(cts, freq > 5)$pid)
id = keep[8]
keep2 <- nuttigedata$pid %in% id
df2 <- nuttigedata[keep2,]


ggplot(nuttigedata, aes(x = Leeftijd, y = PTA54ADS, group = pid, color = Domain)) + 
  geom_point(aes(colour = factor(Domain))) +
  geom_line(data=df2, size=1, alpha = .5) +
  theme_classic()


single <- nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal), 
                       data = df2)
summary(single)

# simple linear model: PTA is a function of the affected domain; here are two levels in this mixed model; 
# 1: timepoints for each patient; 2: genetic domain, with domain the fixed effect and patient the random effect allowing
# the intercept to vary across patient (~1|pid).
random_intercept <-lme(PTA54ADS ~ Domain , 
            random = ~1|pid,   #p. 896
            method = "ML", 
            na.action = na.exclude, 
            control = list(opt="optim"),
            correlation = corAR1(),  #see p.897; timepoints are not equally spaced;use corCAR1 
            data = nuttigedata)
summary(random_intercept)


#add Leeftijd as fixed effect; PTA ~ Domain + Leeftijd` (see. e.g. p.897)
timeRI <- update(random_intercept, .~. + Leeftijd)
summary(timeRI)

#perhaps even an fixed effect of family
#timeRI_fam <- update(timeRI, .~. + Familynrrec)
# this results in problems with fitting the model; leaving it out for now
#check whether adding the fixed effect of time and family improves the fit 
anova(random_intercept,timeRI)

# -> it does; adding both lowers the AIC/BIC
# Now, update the model such that it includes a random slope; the intercepts and the effects of 
# time (Leeftijd) vary across people (pid); this allows for a random effect
timeRS <- update(timeRI, random=~Leeftijd|pid)
summary(timeRS)

#finally, one can expect an interaction of domain x time (leeftijd)
timeRS_interact <- update(timeRS, .~. + Domain * Leeftijd)
summary(timeRS_interact)

#full_model <- update(timeRS_interact, .~. + Familynrrec) does not converge
anova(random_intercept,timeRI,timeRS,timeRS_interact) #compare all models
intervals(timeRS_interact)
#random effect of Leeftijd is 0.71 (0.57, 0.87); indicative of a PTA that varies significantly across people

#plot parameter estimates for the finam model
plot_model(timeRS_interact)





# might not work beyond this point---------
newdat <- expand.grid(Domain=unique(nuttigedata$Domain),
                      Leeftijd=c(min(nuttigedata$Leeftijd),
                            max(nuttigedata$Leeftijd)))

ggplot(nuttigedata, aes(x=Leeftijd, y=PTA54ADS, factor=Domain)) +
  geom_point(aes(colour = factor(Domain))) +
  geom_line(data=newdat,aes(y=predict(timeRS_interact,level=0, newdata=newdat)))


mydf <- ggpredict(timeRS_interact, terms = c("Leeftijd", "Domain"))


fit <- lm(PTA54ADS ~ Domain + Leeftijd + Domain*Leeftijd, data = nuttigedata)
mydf = ggpredict(fit, terms = c("Domain","Leeftijd"))

ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)

plot_model(dat, type = "eff", terms = "Leeftijd")
plot(dat)
random_intercept_time2 <-lme(PTA54ADS ~ Domain*I(Leeftijd^2) , 
                       random = ~1|pid,   #p. 896
                       method = "ML", 
                       na.action = na.exclude, 
                       control = list(opt="optim"),
                       correlation = corAR1(),  #see p.897; timepoints are not equally spaced;use corCAR1 
                       data = nuttigedata)
summary(random_intercept_time2)
random_intercept_time2_rs <- update(random_intercept_time2, random=~Leeftijd|pid)

#plot

ggplot(nuttigedata, aes(x=Leeftijd, y=PTA54ADS, group=Domain)) +
  geom_point()+
  geom_line(color="grey")


#newdat <- expand.grid(Domain=unique(nuttigedata$Domain),
#                      Leeftijd=c(min(nuttigedata$Leeftijd),
#                            max(nuttigedata$Leeftijd)))

newdat <- expand.grid(Leeftijd=c(10,20,30,40,50,60,70,80), Domain=c("LCCL","vWFA1","vWFA2"), Familynrrec=unique(nuttigedata$Familynrrec))
newdat$pred <- predict(model, newdat, level=0)

ggplot(df2, aes(x=Leeftijd, y=PTA54ADS)) +
      geom_point(aes(colour = factor(Domain))) +
      geom_line(aes(group = Domain)) +
      geom_line(aes(y=predict(timeRS_interact), group=Domain, size="Subjects"))
      geom_line(newdat,aes(x=leeftijd,y=pred,colour=Domain))
      geom_line(newdata,aes(x=leeftijd,y=predict(model), group=Domain))
  
print(p)


#----
nm1 <- nlmer(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal) ~ Asym | Domain , nuttigedata, start = c(Asym = 1,xmid = 60, scal = 120))
fm1 = nls(PTA54ADS ~ SSlogis(Leeftijd, Asym, xmid, scal), data=nuttigedata)
fm2 = nls(PTA54ADS ~ SSlogis(Leeftijd, Asym[Domain], xmid[Domain], scal[Domain]), data=nuttigedata)

  
#non-linear mixed models
library(saemix)
saemix.data <- saemixData(name.data       = nuttigedata,
                          name.group      = "Domain",
                          name.predictors = "Leeftijd",
                          name.response   = "PTA54ADS")

model_ll4 <- function(psi,Domain,x) { 
  t   <-x[,1] 
  b   <-psi[Domain,1]
  c   <-psi[Domain,2]
  d   <-psi[Domain,3]
  e   <-psi[Domain,4]
  fpred <- c + (d-c)/(1+exp(b * log(t)-log(e)))
  return(fpred)
}

ll4 <- function(psi,x) { 
  t   <-x 
  b   <-psi[1]
  c   <-psi[2]
  d   <-psi[3]
  e   <-psi[4]
  fpred <- c + (d-c)/(1+exp(b * (log(t)-log(e))))  
  return(fpred)}

x = seq(1, 100, by=1)
psi = c(b=10,c=-90,d=130,e=60)

saemix.model <- saemixModel(model = model_ll4, 
                           psi0  = c(b=10,c=0,d=130,e=60))
saemix.options <- list(map=TRUE, fim=TRUE, ll.is=FALSE, displayProgress=FALSE, seed=632545)
saemix.fit1    <- saemix(saemix.model, saemix.data, saemix.options)

saemix.fit <- saemix.predict(saemix.fit1)
saemix.plot.fits(saemix.fit1)