library(nlme)
library(bbmle)
library(lme4)
library(multcomp)
library(ggplot2)
library(mgcv)
library(mgcViz)
library(lubridate)
library(readxl)
library(plotrix)
library(egg)

##turtle data
HABs_CPUE_01_18_Net_Km <- read_excel("Dropbox/HABs CPUE 01-18_Net_Km.xlsx")
attach(HABs_CPUE_01_18_Net_Km)
StudyDay<-HABs_CPUE_01_18_Net_Km$`Study Day`
HABs_CPUE_01_18_Net_Km$StudyDay<-StudyDay-min(StudyDay)+1

##season for visualization
HABs_CPUE_01_18_Net_Km$Season<-HABs_CPUE_01_18_Net_Km$Month
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="12"]<-"Winter"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="1"]<-"Winter"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="2"]<-"Winter"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="3"]<-"Spring"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="4"]<-"Spring"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="5"]<-"Spring"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="6"]<-"Summer"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="7"]<-"Summer"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="8"]<-"Summer"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="9"]<-"Fall"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="10"]<-"Fall"
HABs_CPUE_01_18_Net_Km$Season[HABs_CPUE_01_18_Net_Km$Season=="11"]<-"Fall"
HABs_CPUE_01_18_Net_Km$Season<-factor(HABs_CPUE_01_18_Net_Km$Season, order = TRUE, levels =c("Winter", "Spring","Summer","Fall"))



###testing variance vs mean relationship for turtle captures##
par(mfrow=c(1,1))

HABs_CPUE_01_18_Net_Km$cuts<-cut(HABs_CPUE_01_18_Net_Km$Green_Captures, c(-1,3,7,11,15,19,23,27,31,35,39,43))
HABs_CPUE_01_18_Net_Km$cuts<-as.factor(HABs_CPUE_01_18_Net_Km$cuts)

meanStats = function(x) c(mean = mean(x))
SEStats = function(x) c(se = std.error(x))

meancuts<-tapply(HABs_CPUE_01_18_Net_Km$Green_Captures, HABs_CPUE_01_18_Net_Km$cuts, meanStats)
SEcuts<-tapply(HABs_CPUE_01_18_Net_Km$Green_Captures, HABs_CPUE_01_18_Net_Km$cuts, SEStats)
plot(meancuts, SEcuts, xlim=c(0,40), ylim=c(0,1))
##variance increases as mean increases, but not exponentially
##In addition, and more subjectively, we want to more heavily weight 
##high capture days as they drive trend
##and account for most of the captures
##Therefore, quasipoisson is best distribution 
##to use for abundance analyses 
##(see Ver Hoef and Boveng 2007)

##Turtle gam model code##
##k's here are just for demonstration purposes
##not the final model

gam(Green_Captures~s(StudyDay, bs="tp", k=9)+s(Month, bs="cc", k=12)+
      ti(StudyDay, Month, bs=c("tp","cc"), k=c(9,12))+s(Net_Km_Hours),
    select=TRUE, family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)

##iteratively increasing k to find optimal k for Study Day

GAM3s<-gam(Green_Captures~s(StudyDay, bs="tp", k=3)+s(Month, bs="cc", k=12)+
                 ti(StudyDay, Month, bs=c("tp","cc"), k=c(3,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM4s<-gam(Green_Captures~s(StudyDay, bs="tp", k=4)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(4,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM5s<-gam(Green_Captures~s(StudyDay, bs="tp", k=5)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(5,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM6s<-gam(Green_Captures~s(StudyDay, bs="tp", k=6)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(6,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM7s<-gam(Green_Captures~s(StudyDay, bs="tp", k=7)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(7,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM8s<-gam(Green_Captures~s(StudyDay, bs="tp", k=8)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(8,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM9s<-gam(Green_Captures~s(StudyDay, bs="tp", k=9)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(9,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM10s<-gam(Green_Captures~s(StudyDay, bs="tp", k=10)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(10,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM11s<-gam(Green_Captures~s(StudyDay, bs="tp", k=11)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(11,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM12s<-gam(Green_Captures~s(StudyDay, bs="tp", k=12)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(12,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM13s<-gam(Green_Captures~s(StudyDay, bs="tp", k=13)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(13,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
GAM14s<-gam(Green_Captures~s(StudyDay, bs="tp", k=14)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(14,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)
##comparison of RMSE of iterated modelss
thr<-sqrt(mean(residuals(GAM3s)^2))
fo<-sqrt(mean(residuals(GAM4s)^2))
fi<-sqrt(mean(residuals(GAM5s)^2))
si<-sqrt(mean(residuals(GAM6s)^2))
seven<-sqrt(mean(residuals(GAM7s)^2))
eig<-sqrt(mean(residuals(GAM8s)^2))
nin<-sqrt(mean(residuals(GAM9s)^2))
te<-sqrt(mean(residuals(GAM10s)^2))
el<-sqrt(mean(residuals(GAM11s)^2))
tw<-sqrt(mean(residuals(GAM12s)^2))
thi<-sqrt(mean(residuals(GAM13s)^2))
plot(c(3,4,5,6,7,8,9,10,11,12,13),c(thr,fo,fi,si,seven,eig,nin,te,el,tw,thi),
     xlab="k", ylab="RMSE")

##once you get to k=7 and higher, minimal changes in RMSE 
##but there are some small differences so need another metric

##Use -REML score from summary/fitting
three<-summary(GAM3s)
four<-summary(GAM4s)
five<-summary(GAM5s)
six<-summary(GAM6s)
sev<-summary(GAM7s)
eight<-summary(GAM8s)
nine<-summary(GAM9s)
ten<-summary(GAM10s)
eleven<-summary(GAM11s)
twelve<-summary(GAM12s)
thirteen<-summary(GAM13s)

##extract -REML scores, looking for minimum
REMLS<-c(three$sp.criterion, four$sp.criterion, five$sp.criterion, six$sp.criterion,
       sev$sp.criterion, eight$sp.criterion, nine$sp.criterion, ten$sp.criterion, eleven$sp.criterion,
       twelve$sp.criterion, thirteen$sp.criterion)
plot(c(3,4,5,6,7,8,9,10,11,12,13), REMLS,
     xlab="k", ylab="-REML")

##Again small differences but minimum -REML is with k=7
##And has the benefit of not overfitting with a small dataset 
##(for a GAM anyway)
##Takeaway: GAM7s is the model to go with

summary(GAM7s)
11GAM7s<-gam(Green_Captures~s(StudyDay, bs="tp", k=7)+s(Month, bs="cc", k=12)+
             ti(StudyDay, Month, bs=c("tp","cc"), k=c(7,12))+s(Net_Km_Hours), select=TRUE,family=quasipoisson(), method="REML", data=HABs_CPUE_01_18_Net_Km)

##individual factor plots
labels=c("2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")

turtviz<-getViz(GAM7s)
turtsd<-plot(sm(turtviz, 1))+ l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  theme_classic()+
  ylim(c(-1,1))+
  ylab("Green Turtle Captures (GAM scale)")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  annotate("text",x=6000, y=0.9, label="a")
  
turtsd

mon<-plot(sm(turtviz, 2))+
  l_fitLine(colour = "red")+
  ylim(c(-1,1))+
  ylab("Green Turtle Captures (GAM scale)")+
  xlab("Month")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+theme_classic()+
  annotate("text",x=12, y=0.9, label="b")

intmsd<-plot(sm(turtviz, 3))+
  l_fitRaster()+
  labs(title = NULL)+
  ylab("Month")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme_classic()+
  annotate("text",x=6800, y=12, label="c")

effort<-plot(sm(turtviz, 4))+ylim(c(-1,1))+l_fitLine(colour = "red")+
  ylab("Green Turtle Captures (GAM scale)")+
  xlab("Net Km Hours")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  theme_classic()+
  annotate("text",x=3, y=0.9, label="d")

gridPrint(turtsd, mon, intmsd, effort)

##plots with raw data
GAM7pred<-predict(GAM7s, se.fit=TRUE)
GAM7fit<-GAM7pred$fit
GAM7se<-GAM7pred$se.fit
crit.t<-qt(0.975, df=df.residual(GAM7s))
Turtfit<-exp(GAM7fit)
turtse<-exp(GAM7se)
Turtupper<-Turtfit+crit.t*turtse
Turtlower<-Turtfit-crit.t*turtse

Turtplot<-ggplot(data=HABs_CPUE_01_18_Net_Km, aes(x=StudyDay, y=Green_Captures))+
  geom_point()+
  facet_grid(~Season)+
  geom_smooth(aes(x=StudyDay, y=Turtfit), se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=Turtlower), se=FALSE, linetype="dashed")+
  geom_smooth(aes(x=StudyDay, y=Turtupper), se=FALSE, linetype="dashed")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  ylab("Green Turtle Captures")
Turtplot
##monthly plot for Appendix 2
Turtplotmonth<-ggplot(data=HABs_CPUE_01_18_Net_Km, aes(x=StudyDay, y=Green_Captures))+
  geom_point()+
  facet_wrap(~Month)+
  geom_smooth(aes(x=StudyDay, y=Turtfit), se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=Turtlower), se=FALSE, linetype="dashed")+
  geom_smooth(aes(x=StudyDay, y=Turtupper), se=FALSE, linetype="dashed")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  ylab("Green Turtle Captures")
Turtplotmonth<-tag_facet(Turtplotmonth, tag_pool = letters)
Turtplotmonth
##mean model-predicted abundances
##beginning and end of study period
Turtvalues<-fitted(GAM7s)
mean2001<-mean(Turtvalues[1:22])
mean2001
std.error(Turtvalues[1:22])

mean2018<-mean(Turtvalues[402:419])
mean2018
std.error(Turtvalues[402:419])

###Turtle growth rates
Cm_turt_data_2001_2018 <- read_excel("Dropbox/Cm BC turt data 2001-2018.xlsx")
GRFP2001<-subset(PapsGRdata, GRStudyDay>0)
GRFP2001$TID<-as.factor(GRFP2001$TID)

##increasing k to find optimum value
##first set holds k for SCL constant, increases for Study Day
GRBAM33i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=3)+ti(MeanSCL,GRStudyDay,k=c(3,3))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM34i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=4)+ti(MeanSCL,GRStudyDay,k=c(3,4))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM35i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=5)+ti(MeanSCL,GRStudyDay,k=c(3,5))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM36i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=6)+ti(MeanSCL,GRStudyDay,k=c(3,6))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM37i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=7)+ti(MeanSCL,GRStudyDay,k=c(3,7))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM38i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=8)+ti(MeanSCL,GRStudyDay,k=c(3,8))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM39i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=9)+ti(MeanSCL,GRStudyDay,k=c(3,9))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)

##Second set holds k for Study Day constant, increases k for SCL
GRBAM43i<-bam(Growth~s(MeanSCL, k=4)+s(GRStudyDay, k=3)+ti(MeanSCL,GRStudyDay,k=c(4,3))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM53i<-bam(Growth~s(MeanSCL, k=5)+s(GRStudyDay, k=3)+ti(MeanSCL,GRStudyDay,k=c(5,3))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM63i<-bam(Growth~s(MeanSCL, k=6)+s(GRStudyDay, k=3)+ti(MeanSCL,GRStudyDay,k=c(6,3))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)
GRBAM73i<-bam(Growth~s(MeanSCL, k=7)+s(GRStudyDay, k=3)+ti(MeanSCL,GRStudyDay,k=c(7,3))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)

###RMSE evaluation
thrthri<-sqrt(mean(residuals(GRBAM33i)^2))
thrfoi<-sqrt(mean(residuals(GRBAM34i)^2))
thrfii<-sqrt(mean(residuals(GRBAM35i)^2))
thrsixi<-sqrt(mean(residuals(GRBAM36i)^2))
thrsevi<-sqrt(mean(residuals(GRBAM37i)^2))
threigi<-sqrt(mean(residuals(GRBAM38i)^2))
thrnini<-sqrt(mean(residuals(GRBAM39i)^2))

fothri<-sqrt(mean(residuals(GRBAM43i)^2))
fithri<-sqrt(mean(residuals(GRBAM53i)^2))
sixthri<-sqrt(mean(residuals(GRBAM63i)^2))
sevthri<-sqrt(mean(residuals(GRBAM73i)^2))



plot(c(thrthri,thrfoi,thrfii,thrsixi,thrsevi,threigi,thrnini,fothri,fithri,sixthri,sevthri),ylab="RMSE")

##increasing k for Mean SCL has little effect on RMSE (Index >=8 is flat)
##increasing k for Study Day has large effect on RMSE (Index <8 declines)
##at k=7-9 for Study Day, relatively few changes but not no changes in RMSE
##But at k=10, too many coefficients for data
##need to be careful

##evaluating changes in Study k only for REML

sum33i<-summary(GRBAM33i)
sum34i<-summary(GRBAM34i)
sum35i<-summary(GRBAM35i)
sum36i<-summary(GRBAM36i)
sum37i<-summary(GRBAM37i)
sum38i<-summary(GRBAM38i)
sum39i<-summary(GRBAM39i)
GRreml<-c(sum33i$sp.criterion,sum34i$sp.criterion,sum35i$sp.criterion,
          sum36i$sp.criterion,sum37i$sp.criterion,sum38i$sp.criterion,
          sum39i$sp.criterion)
plot(GRreml)

###REML minimizes at k=7 for Study Day.
##Considering that degrees of freedom max out at k=9, and that the model results
##Are, qualitatively, almost exactly the same for k=7-9, we go with k=7
summary(GRBAM37i)
GRBAM37i<-bam(Growth~s(MeanSCL, k=3)+s(GRStudyDay, k=7)+ti(MeanSCL,GRStudyDay,k=c(3,7))+s(TID, bs="re")+s(Interval), method="REML",family=gaussian, data=GRFP2001)

##Growth rate visualization
##component smooths for Appendix 2
plot(GRBAM37i, scheme=1)
GRviz<-getViz(GRBAM37i)

GRstdy<-plot(sm(GRviz, 2))+l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  ylab("Turtle Growth Rate (GAM scale)")+
  theme_classic()+
  ylim(-4,4)+
  scale_x_continuous(name="Mean Year", breaks=(seq(0,6560,365)), labels=labels1)+
  annotate("text",x=6000, y=4, label="a")

SCL<-plot(sm(GRviz, 1))+l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  ylab("Turtle Growth Rate (GAM scale)")+
  theme_classic()+
  ylim(-4,4)+
  annotate("text",x=70, y=4, label="b")

intGR<-plot(sm(GRviz, 3))+
  l_fitRaster()+
  scale_y_continuous(name="Mean Year", breaks=(seq(0,6560,365)), labels=labels1)+
  labs(title = NULL)+
  theme_classic()+
  annotate("text",x=73, y=6700, label="c")

turtid<-plot(sm(GRviz, 4))+
  l_points() + labs(title = NULL)+theme_classic()+
  annotate("text",x=3.1, y=3, label="d")

gridPrint(GRstdy, SCL, intGR, turtid)

##3-d for SI?
vis.gam(GRBAM37i, theta=45, phi=30,color="topo", ticktype="detailed",zlab="cm/year")
##same fig but rotated 180 to show low SCL/High study day pattern
vis.gam(GRBAM37i, theta=225, phi=30,color="topo", ticktype="detailed",zlab="cm/year")

##2D fig of just study day
GR<-predict(GRBAM37i, se.fit=TRUE)
GRfit<-GR$fit
GRse.fit<-GR$se.fit
GRqt<-qt(0.975, df=df.residual(GRBAM37i))
GRminCI<-GRfit-GRqt*GRse.fit
GRmaxCI<-GRfit+GRqt*GRse.fit

ggplot(data=GRFP2001, aes(x=GRStudyDay, y=Growth))+
  geom_point()+
  geom_smooth(aes(x=GRStudyDay, y=GRfit), se=FALSE)+
  geom_smooth(aes(x=GRStudyDay, y=GRminCI),linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=GRStudyDay, y=GRmaxCI), linetype="dashed", se=FALSE)+
  geom_rug(sides="r", position="identity", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  scale_x_continuous(name="Mean Year", breaks=(seq(0,6560,365)), labels=labels)+
  theme(axis.text.x=element_text(angle=60, size=18, vjust=0.5, color="black"),
        axis.title.x= element_text(size=18))+
  theme(strip.text.x = element_blank())+
  theme(axis.title.y= element_text(size=18),
        axis.text.y=element_text(size=16, color="black"))+
  theme_classic()+
  ylab("Growth Rate (cm/year)")

##Note that CIs (and therefore r-squared) are wide because of the high number of random effects 
##(each turtle is its own, n=164) for this amount of data (n=202).
##with many turtles only having one recapture
##Without Turtle ID random effect, variation explained is more like 0.18


##SEAGRASS
##Seagrass data
Sebastian_area_transects_1994_2018 <- read_excel("Dropbox/Sebastian area transects_1994-2018.xlsx")
SAT<-Sebastian_area_transects_1994_2018

##taking raw data and taking mean of the data at each site on each sampled day
SiteDailyAverage<-aggregate.data.frame(SAT$TOTVIS, by=list(SAT$`Site#`, SAT$SampleDate), FUN=mean)
SiteDailyAverage$Month<-month(SiteDailyAverage$Group.2)
SiteDailyAverage$StudySite<-as.character(SiteDailyAverage$Group.1)
SiteDailyAverage$StudySite<-as.factor(SiteDailyAverage$StudySite)
SiteDailyAverage$Date<-SiteDailyAverage$Group.2
SiteDailyAverage$StudyDay<-as.numeric(SiteDailyAverage$Date-min(SiteDailyAverage$Date))/(24*60*60)

##adding a constant to TOTVIS, the total seagrass cover, so that it can be transformed
SiteDailyAverage$TOTVIS<-SiteDailyAverage$x/100+0.001
##logit transformation
SiteDailyAverage$logitTOTVIS<-log(SiteDailyAverage$TOTVIS/(1-SiteDailyAverage$TOTVIS))


###creating a season factor for visualization
SiteDailyAverage$Season<-SiteDailyAverage$Month
SiteDailyAverage$Season[SiteDailyAverage$Season=="12"]<-"Winter"
SiteDailyAverage$Season[SiteDailyAverage$Season=="1"]<-"Winter"
SiteDailyAverage$Season[SiteDailyAverage$Season=="2"]<-"Winter"
SiteDailyAverage$Season[SiteDailyAverage$Season=="3"]<-"Spring"
SiteDailyAverage$Season[SiteDailyAverage$Season=="4"]<-"Spring"
SiteDailyAverage$Season[SiteDailyAverage$Season=="5"]<-"Spring"
SiteDailyAverage$Season[SiteDailyAverage$Season=="6"]<-"Summer"
SiteDailyAverage$Season[SiteDailyAverage$Season=="7"]<-"Summer"
SiteDailyAverage$Season[SiteDailyAverage$Season=="8"]<-"Summer"
SiteDailyAverage$Season[SiteDailyAverage$Season=="9"]<-"Fall"
SiteDailyAverage$Season[SiteDailyAverage$Season=="10"]<-"Fall"
SiteDailyAverage$Season[SiteDailyAverage$Season=="11"]<-"Fall"
SiteDailyAverage$Season<-factor(SiteDailyAverage$Season, order = TRUE, levels =c("Winter", "Spring","Summer","Fall"))

###subsetting seagrass data to after 2001 only
##Then resetting the study Day to match turtles
SDA2001_2018<-subset(SiteDailyAverage, StudyDay>=2375)
SDA2001_2018$StudyDay<-SDA2001_2018$StudyDay-min(SDA2001_2018$StudyDay)+31

##seagrass model formulation
SGgam3<-gam(logitTOTVIS~s(StudyDay, k=3)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(3,12))+s(StudySite, bs="re"),
               family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam4<-gam(logitTOTVIS~s(StudyDay, k=4)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(4,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam5<-gam(logitTOTVIS~s(StudyDay, k=5)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(5,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam6<-gam(logitTOTVIS~s(StudyDay, k=6)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(6,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam7<-gam(logitTOTVIS~s(StudyDay, k=7)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(7,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam8<-gam(logitTOTVIS~s(StudyDay, k=8)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(8,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam9<-gam(logitTOTVIS~s(StudyDay, k=9)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(9,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam10<-gam(logitTOTVIS~s(StudyDay, k=10)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(10,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam11<-gam(logitTOTVIS~s(StudyDay, k=11)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(11,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam12<-gam(logitTOTVIS~s(StudyDay, k=12)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(12,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam13<-gam(logitTOTVIS~s(StudyDay, k=13)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(13,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
SGgam14<-gam(logitTOTVIS~s(StudyDay, k=14)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(14,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=SDA2001_2018)
##RMSE eval
sg3<-sqrt(mean(residuals(SGgam3)^2))
sg4<-sqrt(mean(residuals(SGgam4)^2))
sg5<-sqrt(mean(residuals(SGgam5)^2))
sg6<-sqrt(mean(residuals(SGgam6)^2))
sg7<-sqrt(mean(residuals(SGgam7)^2))
sg8<-sqrt(mean(residuals(SGgam8)^2))
sg9<-sqrt(mean(residuals(SGgam9)^2))
sg10<-sqrt(mean(residuals(SGgam10)^2))
sg11<-sqrt(mean(residuals(SGgam11)^2))
sg12<-sqrt(mean(residuals(SGgam12)^2))
sg13<-sqrt(mean(residuals(SGgam13)^2))
sg14<-sqrt(mean(residuals(SGgam14)^2))
plot(c(sg3,sg4,sg5,sg6,sg7,sg8,sg9,sg10,sg11,sg12,sg13,sg14))
##minimal changes past k=9

##REML evaluation
sgsum3<-summary(SGgam3)
sgsum4<-summary(SGgam4)
sgsum5<-summary(SGgam5)
sgsum6<-summary(SGgam6)
sgsum7<-summary(SGgam7)
sgsum8<-summary(SGgam8)
sgsum9<-summary(SGgam9)
sgsum10<-summary(SGgam10)
sgsum11<-summary(SGgam11)
sgsum12<-summary(SGgam12)
sgsum13<-summary(SGgam13)
sgsum14<-summary(SGgam14)
sgreml<-c(sgsum3$sp.criterion,sgsum4$sp.criterion,sgsum5$sp.criterion,sgsum6$sp.criterion,
          sgsum7$sp.criterion,sgsum8$sp.criterion,sgsum9$sp.criterion,sgsum10$sp.criterion,
          sgsum11$sp.criterion,sgsum12$sp.criterion,sgsum13$sp.criterion, sgsum14$sp.criterion)
plot(c(3,4,5,6,7,8,9,10,11,12,13,14), sgreml,
     xlab="k", ylab="-REML")

##small but noticeable differences after k=9
##after k=12 it levels off
##use k=12 to avoid overfitting

summary(SGgam12)

##plotting seagrass
##individual effects
plot(SGgam12, scheme=1)
SGviz<-getViz(SGgam12)
SGsd<-plot(sm(SGviz, 1))+ l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  ylab("Seagrass Cover (GAM scale)")+
  theme_classic()+
  ylim(c(-4,4))+
  scale_x_continuous(name="Mean Year", breaks=(seq(0,6560,365)), labels=labels1)+
  annotate("text",x=6000, y=4, label="a")

SGmon<-plot(sm(SGviz, 2))+
  ylab("Seagrass Cover (GAM scale)")+
  l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  theme_classic()+
  annotate("text",x=12, y=1, label="b")


SGint<-plot(sm(SGviz, 3))+
  l_fitRaster()+
  labs(title = NULL)+
  ylab("Month")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme_classic()+
  annotate("text",x=6700, y=12, label="c")

SGsite<-plot(sm(SGviz, 4))+
  l_points()+
  labs(title = NULL)+
  annotate("text",x=1.5, y=1.3, label="d")



gridPrint(SGsd, SGmon, SGint, SGsite)

##long-term plot
SGpred<-predict(SGgam12, se.fit=TRUE)
SGfit<-SGpred$fit
SGse.fit<-SGpred$se.fit
sgqt<-qt(0.975, df=df.residual(SGgam12))
BTfit<-exp(SGfit)/(1+exp(SGfit))
BTminCI<-exp(SGfit-sgqt*SGse.fit)/(1+exp(SGfit-sgqt*SGse.fit))
BTmaxCI<-exp(SGfit+sgqt*SGse.fit)/(1+exp(SGfit+sgqt*SGse.fit))

SGplot<-ggplot(data=SDA2001_2018, aes(x=StudyDay, y=TOTVIS))+
  geom_point(size=0.6)+
  facet_grid(~Season)+
  geom_smooth(aes(x=StudyDay, y=BTfit), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTminCI), method="loess", linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTmaxCI), method="loess", linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels)+
  ylab("Seagrass Percent Cover")

SGplot
##Monthly seagrass plot for Appendix 2
SGplotmonth<-ggplot(data=SDA2001_2018, aes(x=StudyDay, y=TOTVIS))+
  geom_point(size=0.6)+
  facet_wrap(~Month)+
  geom_smooth(aes(x=StudyDay, y=BTfit), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTminCI), method="loess", linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTmaxCI), method="loess", linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  ylab("Seagrass Percent Cover")
SGplotmonth<-tag_facet(SGplotmonth)
SGplotmonth
##averages at beginning and end
SGvalues<-exp(fitted(SGgam12))/(1+exp(fitted(SGgam12)))
mean2001<-mean(SGvalues[1:22])
mean2001
std.error(SGvalues[1:22])


mean2018<-mean(SGvalues[375:288])
mean2018
std.error(SGvalues[375:288])


##Algaedata
Algaedata<-subset(SAT, ALGCOV>=0)
AlgaeDailyAverage<-aggregate.data.frame(Algaedata$ALGCOV, by=list(Algaedata$`Site#`, Algaedata$SampleDate), FUN=mean)
AlgaeDailyAverage$Month<-month(AlgaeDailyAverage$Group.2)
AlgaeDailyAverage$StudySite<-as.character(AlgaeDailyAverage$Group.1)
AlgaeDailyAverage$Date<-AlgaeDailyAverage$Group.2
AlgaeDailyAverage$StudyDay<-as.numeric(AlgaeDailyAverage$Date-min(AlgaeDailyAverage$Date))/(24*60*60)
AlgaeDailyAverage$ALGCOV<-AlgaeDailyAverage$x

##creating Season for visualization
AlgaeDailyAverage$Season<-AlgaeDailyAverage$Month
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="12"]<-"Winter"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="1"]<-"Winter"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="2"]<-"Winter"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="3"]<-"Spring"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="4"]<-"Spring"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="5"]<-"Spring"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="6"]<-"Summer"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="7"]<-"Summer"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="8"]<-"Summer"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="9"]<-"Fall"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="10"]<-"Fall"
AlgaeDailyAverage$Season[AlgaeDailyAverage$Season=="11"]<-"Fall"
AlgaeDailyAverage$Season<-factor(AlgaeDailyAverage$Season, order = TRUE, levels =c("Winter", "Spring","Summer","Fall"))
AlgaeDailyAverage$StudySite<-as.factor(AlgaeDailyAverage$StudySite)

##subsetting algae data to after 2001
Algae2001_2018<-subset(AlgaeDailyAverage, StudyDay>=2375)
Algae2001_2018$StudyDay<-Algae2001_2018$StudyDay-min(Algae2001_2018$StudyDay)+31

##Transformation of algal cover data for logit
Algae2001_2018$ALGCOV<-Algae2001_2018$ALGCOV/100+0.001
Algae2001_2018$logitALGCOV<-log(Algae2001_2018$ALGCOV/(1-Algae2001_2018$ALGCOV))

##algae model
AlggamTI<-gam(logitALGCOV~s(StudyDay, k=9)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(9,12))+s(StudySite, bs="re"),
              family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
##increasing k
Algam3<-gam(logitALGCOV~s(StudyDay, k=3)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(3,12))+s(StudySite, bs="re"),
              family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam4<-gam(logitALGCOV~s(StudyDay, k=4)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(4,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam5<-gam(logitALGCOV~s(StudyDay, k=5)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(5,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam6<-gam(logitALGCOV~s(StudyDay, k=6)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(6,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam7<-gam(logitALGCOV~s(StudyDay, k=7)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(7,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam8<-gam(logitALGCOV~s(StudyDay, k=8)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(8,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam9<-gam(logitALGCOV~s(StudyDay, k=9)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(9,12))+s(StudySite, bs="re"),
           family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam10<-gam(logitALGCOV~s(StudyDay, k=10)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(10,12))+s(StudySite, bs="re"),
            family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam11<-gam(logitALGCOV~s(StudyDay, k=11)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(11,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam12<-gam(logitALGCOV~s(StudyDay, k=12)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(12,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam13<-gam(logitALGCOV~s(StudyDay, k=13)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(13,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam14<-gam(logitALGCOV~s(StudyDay, k=14)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(14,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam15<-gam(logitALGCOV~s(StudyDay, k=15)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(15,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)
Algam16<-gam(logitALGCOV~s(StudyDay, k=16)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(16,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)

AIC(Algam3,Algam4,Algam5,Algam6,Algam7,Algam8,Algam9,Algam10,Algam11,Algam12,Algam13,Algam14,Algam15,Algam16)

##RMSE
Alg3<-sqrt(mean(residuals(Algam3)^2))
Alg4<-sqrt(mean(residuals(Algam4)^2))
Alg5<-sqrt(mean(residuals(Algam5)^2))
Alg6<-sqrt(mean(residuals(Algam6)^2))
Alg7<-sqrt(mean(residuals(Algam7)^2))
Alg8<-sqrt(mean(residuals(Algam8)^2))
Alg9<-sqrt(mean(residuals(Algam9)^2))
Alg10<-sqrt(mean(residuals(Algam10)^2))
Alg11<-sqrt(mean(residuals(Algam11)^2))
Alg12<-sqrt(mean(residuals(Algam12)^2))
Alg13<-sqrt(mean(residuals(Algam13)^2))
Alg14<-sqrt(mean(residuals(Algam14)^2))
Alg15<-sqrt(mean(residuals(Algam15)^2))
Alg16<-sqrt(mean(residuals(Algam16)^2))
Alg17<-sqrt(mean(residuals(Algam17)^2))
Alg18<-sqrt(mean(residuals(Algam18)^2))

plot(c(Alg3,Alg4,Alg5,Alg6,Alg7,Alg8,Alg9,Alg10,Alg11,Alg12,Alg13,Alg14,Alg15,Alg16,Alg17,Alg18))

##REML evaluation for algae
Algsum3<-summary(Algam3)
Algsum4<-summary(Algam4)
Algsum5<-summary(Algam5)
Algsum6<-summary(Algam6)
Algsum7<-summary(Algam7)
Algsum8<-summary(Algam8)
Algsum9<-summary(Algam9)
Algsum10<-summary(Algam10)
Algsum11<-summary(Algam11)
Algsum12<-summary(Algam12)
Algsum13<-summary(Algam13)
Algsum14<-summary(Algam14)
Algsum15<-summary(Algam15)
Algsum16<-summary(Algam16)

Algreml<-c(Algsum3$sp.criterion,Algsum4$sp.criterion,Algsum5$sp.criterion,Algsum6$sp.criterion,
          Algsum7$sp.criterion,Algsum8$sp.criterion,Algsum9$sp.criterion,Algsum10$sp.criterion,
          Algsum11$sp.criterion,Algsum12$sp.criterion,Algsum13$sp.criterion, Algsum14$sp.criterion,
          Algsum15$sp.criterion,Algsum16$sp.criterion)

plot(c(3,4,5,6,7,8,9,10,11,12,13,14,15,16), Algreml,
     xlab="k", ylab="-REML")

##Going with Algam16. Has similar but slightly RMSE than higher ks
##but is a local minimum for REML score. 
##Again, wary of overfitting so
##going with lower k
summary(Algam16)
Algam16<-gam(logitALGCOV~s(StudyDay, k=16)+s(Month, bs="cc", k=12)+ti(StudyDay, Month, bs=c("tp", "cc"), k=c(16,12))+s(StudySite, bs="re"),
             family=gaussian(link="identity"), method="REML", data=Algae2001_2018)

###Algae visualization
##component smooths
plot(Algam16, scheme=1)

Algviz<-getViz(Algam16)
Algsd<-plot(sm(Algviz, 1))+ l_fitLine(colour = "red")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  theme_classic()+
  ylim(c(-4,4))+
  ylab("Algae Cover (GAM scale)")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  annotate("text",x=6000, y=4, label="a")

Algmon<-plot(sm(Algviz, 2))+ l_fitLine(colour = "red")+
  ylab("Algae Cover (GAM scale)")+
  l_ciLine(mul = 2, colour = "blue", linetype = 2)+
  theme_classic()+
  annotate("text",x=12, y=0.7, label="b")

Algint<-plot(sm(Algviz, 3))+
  l_fitRaster()+
  labs(title = NULL)+
  ylab("Month")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme_classic()+
  annotate("text",x=6700, y=12, label="c")

Algsite<-plot(sm(Algviz, 4))+l_points()+labs(title = NULL)+
  annotate("text",x=1.5, y=0.35, label="d")

gridPrint(Algsd, Algmon, Algint, Algsite)
##Long-term
Algpred<-predict(Algam16, se=TRUE)
Algfit<-Algpred$fit
Algse.fit<-Algpred$se.fit
BTalg<-exp(Algfit)/(1+exp(Algfit))
AlgminCI<-exp(Algfit-1.96*Algse.fit)/(1+exp(Algfit-1.96*Algse.fit))
AlgmaxCI<-exp(Algfit+1.96*Algse.fit)/(1+exp(Algfit+1.96*Algse.fit))

Algplot<-ggplot(data=Algae2001_2018, aes(x=StudyDay, y=ALGCOV))+
  geom_point(size=0.6)+
  geom_smooth(aes(x=StudyDay, y=BTalg), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgminCI), linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgmaxCI), linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  ylab("Algae Percent Cover")

Algplot

##Monthly Algae plot for Appendix 2
Algplotmonth<-ggplot(data=Algae2001_2018, aes(x=StudyDay, y=ALGCOV))+
  facet_wrap(~Month)+
  geom_point(size=0.6)+
  geom_smooth(aes(x=StudyDay, y=BTalg), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgminCI), linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgmaxCI), linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  ylab("Algae Percent Cover")
Algplotmonth<-tag_facet(Algplotmonth)
Algplotmonth
###Algae averages and SEs
View(Algae2001_2018)
Algvalues<-exp(fitted(Algam16))/(1+exp(fitted(Algam16)))
mean2001<-mean(Algvalues[1:22])
mean2001

mean2008<-mean(Algvalues[105:118])
mean2008
std.error(Algvalues[105:118])

mean2012<-mean(Algvalues[204:234])
mean2012
std.error(Algvalues[204:234])

##Combined plots of turt, algae, seagrass
labels1=c("2001", "", "2003", "", "2005", "", "2007", "", "2009", "", "2011", "", "2013", "", "2015", "", "2017", "")
tag_facet <- function(p, open = "", close = "", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}

Turtplotall<-ggplot(data=HABs_CPUE_01_18_Net_Km, aes(x=StudyDay, y=Green_Captures))+
  geom_point(size=0.6)+
  facet_grid(~Season)+
  geom_smooth(aes(x=StudyDay, y=Turtfit), se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=Turtlower), se=FALSE, linetype="dashed")+
  geom_smooth(aes(x=StudyDay, y=Turtupper), se=FALSE, linetype="dashed")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y= element_text(size=18),
        axis.text.y=element_text(size=18, color="black"))+
  theme(strip.text.x= element_text(size=18))+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  ylab("Green Turtle Captures")
  
Turtplotall<-tag_facet(Turtplotall, tag_pool=c("a","b","c","d"))

SGplotall<-ggplot(data=SDA2001_2018, aes(x=StudyDay, y=TOTVIS))+
  geom_point(size=0.6)+
  facet_grid(~Season)+
  geom_smooth(aes(x=StudyDay, y=BTfit), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTminCI), method="loess", linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=BTmaxCI), method="loess", linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  theme(axis.text.x=element_text(angle=60, vjust=0.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(strip.text.x = element_blank())+
  theme(axis.title.y= element_text(size=18),
        axis.text.y=element_text(size=16, color="black"))+
  ylab("Seagrass Percent Cover")
SGplotall<-tag_facet(SGplotall, tag_pool=c("e","f","g","h"))
Algplotall<-ggplot(data=Algae2001_2018, aes(x=StudyDay, y=ALGCOV))+
  geom_point(size=0.6)+
  facet_grid(~Season)+
  geom_smooth(aes(x=StudyDay, y=BTalg), method="loess", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgminCI), linetype="dashed", se=FALSE)+
  geom_smooth(aes(x=StudyDay, y=AlgmaxCI), linetype="dashed", se=FALSE)+
  theme_classic()+
  geom_rug(sides="r", position="jitter", alpha=0.5)+
  geom_rug(sides="b", position="identity", alpha=0.5)+
  scale_x_continuous(name="Year", breaks=(seq(0,6560,365)), labels=labels1)+
  theme(axis.text.x=element_text(angle=60, size=18, vjust=0.5, color="black"),
        axis.title.x= element_text(size=18))+
  theme(strip.text.x = element_blank())+
  theme(axis.title.y= element_text(size=18),
        axis.text.y=element_text(size=16, color="black"))+
  ylab("Algae Percent Cover")
Algplotall<-tag_facet(Algplotall, tag_pool=c("i","j","k","l"))

Allplot<-gridPrint(Turtplotall, SGplotall, Algplotall)
##Combined plot of study day effect



