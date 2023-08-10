#####Scripts for "Spatiotemporal organization of cryptic North American Culex species along an urbanization gradient" Arsenault-Benoit and Fritz 2023


#Load packages
library(readxl)
library(vegan)
require(ggplot2)
require(grid)
require(tidyverse)
library(dplyr)
library(ggrepel)
library(permute)
library(lattice)
library(vcd)
library(lme4)
library(car)
library(HardyWeinberg)
library(stats)
library(MASS)
library(eulerr)
library(comecol)

#Read in data *Change path based on data location
catch <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/MD_3yr_totalcatch.xlsx", 
                    
                    col_types = c("text", "text", "date", 
                                  "text", "numeric", "text", "text", 
                                  "text", "numeric", "numeric", "numeric"))

summary_matrix <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/summary.matrix.xlsx")
summatrix21 <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/summatrix21.xlsx")
md_3 <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/MD_3yr_finals.xlsx")
env_phen <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/env_phen.xlsx")
dat_phen<-read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/dat_phen.xlsx")
pq_HWE <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/pq HWE.xlsx")


#############Trap Catch Analysis#######
#Summary analysis of all unfed cryptic Culex collected (before molecular ID)
#head(catch)
catch$WOY<-ordered(catch$WOY)
catch$Year<-as.factor(catch$Year)
catch$log_coll<-log(catch$Unfed_PerTrap+1)
catch$Date<-as.Date(catch$Date)
catch$Site<-as.factor(catch$Site)

#June-October 2019-2021
catch_reg<-catch%>%filter(Year!=2021)
catch_reg$WOY<-ordered(catch_reg$WOY)
catch_reg$Year<-as.factor(catch_reg$Year)
catch_reg$log_coll<-log(catch_reg$Unfed_PerTrap+1)
catch_reg$Date<-as.Date(catch_reg$Date)
catch_reg$Site<-as.factor(catch_reg$Site)
catch_reg$Season<-factor(catch_reg$Season, levels=c("Early","Mid","Late"))


sum.matrix<-summary_matrix
sum.matrix$Season<-factor(sum.matrix$Season, ordered=TRUE, levels=c("Early", "Mid", "Late"))
sum.matrix$Class<-factor(sum.matrix$Class)

#April-October, 2021
catch21<-catch%>%filter(Year==2021)
#View(catch21)
catch21$Season<-factor(catch21$Season, levels=c("Very Early","Early","Mid","Late"))

#generate mean =/- se to summarize catch. SE was calculated by generating sd, then divided by sqrt(N) for each group. 
table(catch_reg$Year, catch_reg$Prim_Class)
table(catch_reg$Season, catch_reg$Prim_Class)

table(catch21$Prim_Class)
table(catch21$Season)
table(catch21$Season, catch21$Prim_Class)

aggregate(Unfed_PerTrap ~ Prim_Class + Season, data = catch21, FUN = mean)
aggregate(Unfed_PerTrap~Prim_Class,data = catch_reg, FUN = mean )
aggregate(Unfed_PerTrap~Season,data = catch_reg, FUN = mean )
aggregate(Unfed_PerTrap~Year,data = catch_reg, FUN = mean )
aggregate(Unfed_PerTrap~Year,data = catch_reg, FUN = sd )
aggregate(Unfed_PerTrap~Prim_Class*Season,data = catch_reg, FUN = sd )


aggregate(Unfed_PerTrap~Prim_Class,data = catch21, FUN = mean )
aggregate(Unfed_PerTrap~Season,data = catch21, FUN = mean )
aggregate(Unfed_PerTrap~Season,data = catch21, FUN = sd )
aggregate(Unfed_PerTrap~Prim_Class,data = catch21, FUN = sd )
aggregate(Unfed_PerTrap~Season*Prim_Class,data = catch21, FUN = mean )
aggregate(Unfed_PerTrap~Season*Prim_Class,data = catch21, FUN = sd )

#Modeling for total catch analysis. Begin with June-October, 3 years
ctsy<-glmer.nb(Unfed_PerTrap~Season*Prim_Class+(1|Year),  data=catch_reg)
summary(ctsy)

scy<-glmer.nb(Unfed_PerTrap~Season+Prim_Class+(1|Year),  data=catch_reg)
summary(scy)

cts<-glmer.nb(Unfed_PerTrap~Season+(1|Year),  data=catch_reg)
summary(cts)

ctc<-glmer.nb(Unfed_PerTrap~Prim_Class+(1|Year),  data=catch_reg)
summary(ctc)



##Full model is ctsy (interaction)
##Then scy (addititve)
##Then cts (just season)
##Then ctc (justclass)
anova(scy,ctsy) #0.7539
anova(cts, ctsy) #0.452
anova(ctsy, ctc) #0.099
anova(cts,ctc) #0.0



##2021 only

ctsy21<-glm.nb(Unfed_PerTrap~Season*Prim_Class,  data=catch21)
summary(ctsy21)

add21<-glm.nb(Unfed_PerTrap~Season+Prim_Class,  data=catch21)
summary(add21)

ctc21<-glm.nb(Unfed_PerTrap~Prim_Class,  data=catch21)
summary(ctc21)

cts21<-glm.nb(Unfed_PerTrap~Season,  data=catch21)
summary(cts21)

#model comparisons
anova(ctsy21, add21)#***

#temporal autocorrelation assessment
plot(acf(resid(ctsy), type="covariance", lag = 6))


#figures for total trap catch summary

library(viridis)

unfed.box<-ggplot(catch_reg,aes(x=Prim_Class, y=Unfed_PerTrap, fill=Season)) +
  geom_boxplot(outlier.shape=NA, width=0.8) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(x="Site Class", y="Culex Collected per Trap")+
  theme_classic() +scale_fill_brewer(palette = "Accent")
unfed.box

summatrix21$Season<-factor(summatrix21$Season, ordered=TRUE, levels=c("VeryEarly","Early", "Mid", "Late"))
summatrix21$Class<-factor(summatrix21$Class)


unfed.box21<-ggplot(catch21,aes(x=Prim_Class, y=Unfed_PerTrap, fill=Season)) +
  geom_boxplot(outlier.size=0.4, width=0.8) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.5, alpha=0.5) +
  labs(x="Site Class", y="Culex Collected per Trap") +
  theme_classic() +scale_fill_brewer(palette = "Accent")
unfed.box21


###HWE Analysis####

X <-as.matrix(pq_HWE[-c(1:3)])
colnames(X) <- c("PP","PQ","QQ")

Results.wc <- HWChisqMat(X)
Output.wc <- cbind(Results.wc$chisqvec,Results.wc$pvalvec)
colnames(Output.wc)<-c("X2", "pval")
HWE.output<-cbind(pq_HWE[,1:3], Output.wc)

HWE.output #find any significant rows- start by eliminating rows with no quinqs, where X2=NA and pval=1
HWE.poss<-which(HWE.output$pval!=1)
HWE.poss
HWE.output.rdx<-HWE.output %>% filter(pval!=1)
View(HWE.output.rdx)
which(HWE.output.rdx$pval>0.05)


####Community Analyses####
#Clean up datasheet of molecular IDs
md_3 $Prim_Class<-factor(md_3 $Prim_Class)
md_3 $Sec_Class<-factor(md_3 $Sec_Class)
md_3 $Final_ID<-factor(md_3 $Final_ID)
md_3 $Coll_Event<-factor(md_3 $Coll_Event, ordered=TRUE, levels=c(1,2,3,4,5,6,7,8,9,10,11))
md_3 $Season<-factor(md_3 $Season, ordered=TRUE, levels=c("VeryEarly","Early", "Mid", "Late"))
md_3 $Year<-factor(md_3 $Year)
md_3 $Location<-factor(md_3 $Location)
summary(md_3)

#Start building env. matrix; make site x species matrix(dat_p)

uni_sp = unique(md_3$Final_ID)
uni_ev = unique(md_3$TotID)

dat_p<-table(md_3$TotID_Yr,md_3$Final_ID)
#dat_p
head(dat_p)

cols_to_keep = c('Prim_Class', 'Sec_Class')
Tot_env = aggregate(md_3[ , cols_to_keep], by = list(md_3$TotID), function(x) x[1])
row.names(Tot_env) = Tot_env[ , 1]
Tot_env = Tot_env[ , -1]
head(Tot_env)

env_phen$Prim_Class<-factor(env_phen$Prim_Class)
env_phen$Sec_Class<-factor(env_phen$Sec_Class)
env_phen$Tot_ID_Yr<-as.character(env_phen$Tot_ID_Yr)
env_phen$Coll_Event<-factor(env_phen$Coll_Event, ordered=TRUE, levels=c(1,2,3,4,5,6,7,8,9,10,11)) 
env_phen$Season<-factor(env_phen$Season, ordered=TRUE, levels=c("VeryEarly", "Early", "Mid", "Late"))
env_phen$Location<-factor(env_phen$Location)
dim(env_phen)

uni_env<-unique(env_phen$Tot_ID_Yr)
env_phen[order(env_phen$Tot_ID_Yr),]
all.equal(rownames(dat_p), env_phen$Tot_ID_Yr) ### Make sure order of events is same

#Set up as long format for unconstrained analyses

table(md_3$TotID,md_3$Final_ID)
#Create TotID_m sheet from above table
TotID_m <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/TotID_m.xlsx")
Tot_m<-as.data.frame(TotID_m)
Tot_m$pip<-as.numeric(Tot_m$pip)
Tot_m$quinq<-as.numeric(Tot_m$quinq)
Tot_m$rest<-as.numeric(Tot_m$rest)
Tot_dat <- Tot_m[,-1]
rownames(Tot_dat) <- Tot_m[,1]


#####Contingency Analysis##### Do counts of species vary by site class and season?
###All Years###
dat_phen$Season=ordered(dat_phen$Season, levels=c("Early", "Mid", "Late"))


##Count by class 
dat_type<-table(dat_phen$Prim_Class,dat_phen$Final_ID)
head(dat_type)
chisq.test(dat_type)

##Count by Season
dat_seas<-table(dat_phen$Season,dat_phen$Final_ID)
head(dat_seas)
chisq.test(dat_seas)

Table = table(dat_phen$Season, dat_phen$Prim_Class, dat_phen$Final_ID)
ftable(Table)                     # Display a flattened table
mantelhaen.test(Table)

##2021 only
md21<-md_3%>%filter(Year=="2021")
md21$Season=ordered(md21$Season, levels=c("VeryEarly", "Early", "Mid", "Late"))
dat_type<-table(md21$Prim_Class,md21$Final_ID)
head(dat_type)

assoc(dat_type, shade=TRUE,legend=TRUE)

chisq.test(dat_type)

##Count by Season
dat_seas<-table(md21$Season,md21$Final_ID)
head(dat_seas)
chisq.test(dat_seas)


Table21 = table(md21$Season, md21$Prim_Class, md21$Final_ID)
ftable(Table21)                     # Display a flattened table
mantelhaen.test(Table21)



###NMDS##### Unconstrained ordination 
#####NMDS- all years but exclude very early season

dat_reg<-read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/dat_reg.xlsx")
dat_reg$pip<-as.numeric(dat_reg$pip)
dat_reg$quinq<-as.numeric(dat_reg$quinq)
dat_reg$rest<-as.numeric(dat_reg$rest)

dat_reg<-as.data.frame(dat_reg)
Tot_dat_reg <- dat_reg[,-1]
rownames(Tot_dat_reg) <- dat_reg[,1]
Tot_dat_reg<-Tot_dat_reg%>%filter_all(any_vars(. != 0))

env_3<-env_phen[,-1] #remove Tot_ID_Yr
env_3<-env_3[,-10]#removeYear
env_3<-env_3[,-2] #remove CollNum
env_3<-env_3[,-8]
env_3<-env_3[,-6]
head(env_3)
env_3<-env_3%>% distinct(TotID, .keep_all = TRUE)
env_reg<-env_3%>%filter(Season!="VeryEarly")
env_reg<-env_reg%>% distinct(TotID, .keep_all = TRUE)

env_reg[order(env_reg$TotID),]
dat_reg[order(dat_reg$TotID),]

library(vegan)
nmds_reg<-metaMDS(Tot_dat_reg,"bray", k=3) #change values of K to generate pval in Table S3
plot(nmds_reg)
nmds_reg$stress
en_reg = envfit(nmds_reg, env_reg, permutations = 999, na.rm = TRUE)
en_reg
plot(en_reg)

data.scores_reg = as.data.frame(scores(nmds_reg)$sites)
data.scores_reg$Prim_Class<-env_reg$Prim_Class
data.scores_reg$Season<-env_reg$Season
data.scores_reg$combo<-paste(data.scores_reg$Prim_Class, "_", data.scores_reg$Season)
env_reg$combo<-paste(env_reg$Prim_Class, "_", env_reg$Season)

ggcombo = ggplot(data = data.scores_reg, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores_reg, aes(colour = combo), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c( "lightskyblue",  "darkslateblue", "steelblue3","plum3",  "purple4","darkorchid3", "orange1",  "brown1","darkorange3")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Class x Season")
ggcombo + theme(legend.position = "bottom")

ggcombo



library(devtools)
install_github("fawda123/ggord")
library("ggord")

p.s <- ggord(nmds_reg, data.scores_reg$season, poly = FALSE, polylntyp = data.scores_reg$season, grp_title="Season", xlims=c(-1.5, 2), ylims=c(-1.5,1.5), show.legend=FALSE)
p.s + scale_linetype_manual(values = c('dashed', 'dashed', 'dashed'))

p.c<- ggord(nmds_reg, data.scores_reg$Prim_Class, poly = FALSE, polylntyp = data.scores_reg$Prim_Class,grp_title="Site Class", xlims=c(-1.5, 2), ylims=c(-1.5,1.5))
p.c + scale_linetype_manual(values = c('dashed', 'dashed', 'dashed'))

ordiareatest(nmds_reg, groups=env_reg$combo, area = c("hull", "ellipse"), kind = "sd", permutations = 999, parallel = getOption("mc.cores"))


###NMDS 2021
table(md21$TotID,md21$Final_ID)
dat_21<-read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/dat_21.xlsx")
dat_21$pip<-as.numeric(dat_21$pip)
dat_21$quinq<-as.numeric(dat_21$quinq)
dat_21$rest<-as.numeric(dat_21$rest)

dat_21<-as.data.frame(dat_21)

dat_21$pip<-as.numeric(dat_21$pip)
dat_21$quinq<-as.numeric(dat_21$quinq)
dat_21$rest<-as.numeric(dat_21$rest)


Tot_dat21 <- dat_21[,-1]
rownames(Tot_dat21) <- dat_21[,1]
Tot_dat21<-Tot_dat21%>%filter_all(any_vars(. != 0))

env_21 <- read_excel("~/Desktop/Research/MD_Culex_manuscript/Ch2-to share/MD_Culex_Submit/Data and Scripts/env_21.xlsx")
head(env_21)

env_21<-env_21[,-10]#removeYear
env_21<-env_21[,-2] #remove CollNum
#env_3<-env_3[,-8]
#env_21<-env_21[,-6]


nmds21<-metaMDS(Tot_dat21,"bray")
plot(nmds21)
nmds21$stress
en21 = envfit(nmds21, env_21, permutations = 999, na.rm = TRUE)
en21
plot(en21)

##plot with ggplot
data.scores.21 = as.data.frame(scores(nmds21)$sites)
data.scores.21$ID=rownames(data.scores.21)
data.scores.21$Season=env_21$Season[match(data.scores.21$ID, env_21$TotID)]
data.scores.21$Prim_Class=env_21$Prim_Class[match(data.scores.21$ID, env_21$TotID)]
env_21$Combo<-paste(env_21$Prim_Class, "_", env_21$Season)
data.scores.21$Combo<-paste(data.scores.21$Prim_Class, "_", data.scores.21$Season)
head(data.scores.21)

#View(data.scores)

en_coord_cont = as.data.frame(scores(en21, "vectors")) * ordiArrowMul(en21)
en_coord_cat = as.data.frame(scores(en21, "factors")) * ordiArrowMul(en21)
en_coord_cat<-en_coord_cat[-c(1:148),] #remove sites and sites x events

ggcombo21 = ggplot(data = data.scores.21, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores.21, aes(colour = Combo), size = 3, alpha = 0.5) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Combo")

ggcombo21



###CCA#### Constrained ordination
md_3$Season<-ordered(md_3$Season, levels=c("VeryEarly", "Early", "Mid", "Late"))
md_3$Combo<-paste(md_3$Location,"-",md_3$Season)
md.reg<-md_3%>%filter(Season!="VeryEarly")
dat_seas<-table(md.reg$Combo, md.reg$Final_ID)
dat_seas
target2<-as.vector(rownames(dat_seas))

env_phen$Combo<-paste(env_phen$Location,"-", env_phen$Season)
env_seas<-env_phen[ !duplicated(env_phen$Combo), ]
env_seas<-env_seas%>%filter(Season!="VeryEarly")
env_seas<-env_seas %>% arrange(factor(Combo, levels = target2))
#View(env_seas)
#str(env_seas)

#Remove columns too specific for the site by environment matrix
env_seas<-env_seas[,-25]
env_seas<-env_seas[,-17]
env_seas<-env_seas[,-14]
env_seas<-env_seas[,-11]
env_seas<-env_seas[,-10]
env_seas<-env_seas[,-9]
env_seas<-env_seas[,-8]
env_seas<-env_seas[,-3]
env_seas<-env_seas[,-2]
env_seas<-env_seas[,-1]

rownames(env_seas)<-paste(env_seas$Location,"-",env_seas$Season)
all.equal(rownames(dat_seas), rownames(env_seas))

#regularseason


cca_min=cca(dat_seas~1, env_seas)
cca_full=cca(dat_seas~., env_seas)

ordistep(cca_min, scope=formula(cca_full))

cca_tot6=cca(dat_seas~ Prim_Class+ DCC+Y, env_seas)
cca_tot6
RsquareAdj(cca_tot6)
anova(cca_tot6)
anova.cca(cca_tot6, by="term")
anova.cca(cca_tot6, by="margin")
vif.cca(cca_tot6)

with(env_seas, levels(Prim_Class))
scl<-3
colvec<-c("blue", "purple", "orange")

with(env_seas, levels(Sec_Class))
scl<-3
colvec2<-c("darkgoldenrod", "azure4", "darkgreen", "brown3")

plot(cca_tot6, type='n', scaling=scl, xlim=c(-8,8), ylim=c(-3, 3))
ordiellipse(cca_tot6, groups = env_seas$Season, draw = "polygon",lty = 0, col = c( "palegreen","slategray3","lightsalmon1"),label=FALSE)
with(env_seas, points(cca_tot6, display = "sites", col = colvec[Prim_Class],
                      scaling = scl, pch = 21, cex=0.5, bg = colvec[Prim_Class]))
text(cca_tot6, display='cn', col='black')
orditorp(cca_tot6, display='sp', cex=1, scaling=3, col='blue')
with(env_seas, legend(3, -4, cex=0.8,legend = levels(Prim_Class), y.intersp = 0.2, x.intersp=0.2, bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))


#plot(cca_tot6, type='n', scaling=3,  xlim=c(-8,8), ylim=c(-5,5))
#ordiellipse(cca_tot6, groups = env_seas$Season, draw = "polygon",lty = 0, col = c( "palegreen","slategray3","lightsalmon1"),label=FALSE)
#orditorp(cca_tot6, display='sp', cex=1, scaling=3, col='blue')
#orditorp(cca_tot6,display ="sites",  cex=0.5, scaling=3, col="red") ## Orditop includes text for some sites to avoid overcrowding, others are plotted as data points
#text(cca_tot6, display='cn', col='black')


###Multiscale ordination for spatial autocorrelation
xy.full<-cbind(env_seas$X,env_seas$Y)
mso.full1 <- mso(cca_tot6, xy.full , grain=0.1, perm=999)
msoplot(mso.full1, legend=FALSE)
mso.full1$vario


cca_tot7=cca(dat_seas~ Sec_Class+PercentImp_Mean+DCC+Y, env_seas)
cca_tot7
RsquareAdj(cca_tot7)
anova.cca(cca_tot7)
anova.cca(cca_tot7, by="term")
anova.cca(cca_tot7, by="margin")
vif.cca(cca_tot7)

p.cca7<-plot(cca_tot7, type='n', scaling=3, ylim=c(-5,5), xlim=c(-10,10))
orditorp(cca_tot7, display='sp', cex=1, scaling=1, col='blue')
#orditorp(cca_tot7,display ="sites",cex=0.5, scaling=2, col='black') ## Orditop includes text for some sites to avoid overcrowding, others are plotted as data points
text(cca_tot7, display='cn', col='black')
ordiellipse(cca_tot7, groups = env_seas$Season, draw = "polygon",lty = 0, col = c("palegreen","slategray3","lightsalmon1"),label=FALSE)

#Revise plots with site scores
plot(cca_tot7, type='n', scaling=scl, xlim=c(-8,8), ylim=c(-3, 3))
ordiellipse(cca_tot7, groups = env_seas$Season, draw = "polygon",lty = 0, col = c( "palegreen","slategray3","lightsalmon1"),label=FALSE)
with(env_seas, points(cca_tot6, display = "sites", col = colvec2[Sec_Class],
                      scaling = scl, pch = 21, cex=0.5, bg = colvec2[Sec_Class]))
text(cca_tot7, display='cn', col='black')
orditorp(cca_tot7, display='sp', cex=1, scaling=3, col='blue')
with(env_seas, legend(3, 1, cex=0.8,legend = levels(Sec_Class), y.intersp = 0.2, x.intersp=0.2, bty = "n",
                      col = colvec2, pch = 21, pt.bg = colvec2))

#Variance partitioning

pip3<-as.vector(dat_seas[,1])
quinq3<-as.vector(dat_seas[,2])
rest3<-as.vector(dat_seas[,3])

dat.seas<-as.matrix(cbind(pip3, quinq3, rest3))

vp<-varpart(dat.seas,~Season,~X+Y, ~DCC+NDVI+ Pop+PercentImp_Mean+WaterTbl+Elev, data=env_seas, chisquare = TRUE)
vp
plot(vp, Xnames=c("Season", "Space", "Environment"))

plot_vp_euler3(
  vp,
  names = c("Season", "Space", "Environment"),
  col = c("brown3", "skyblue3", "orange")
)

mso.full <- mso(cca_tot7,xy.full , grain=0.1, perm=999)
msoplot(mso.full)
mso.full$vario

###2021 Only###
md_3$Combo<-paste(md_3$Location,"-",md_3$Season)
md.21<-md_3%>%filter(Year==2021)
View(md.21)
dat_seas.21<-table(md.21$Combo, md.21$Final_ID)
dat_seas.21


target<-as.vector(rownames(dat_seas.21))
env_seas.21<-env_phen[ !duplicated(env_phen$Combo), ]
env_seas.21<-env_seas.21 %>% arrange(factor(Combo, levels = target))

#remove columns that are too specific for this step and make rownames; rownames need to be the same
env_seas.21<-env_seas.21[,-25]
env_seas.21<-env_seas.21[,-11]
env_seas.21<-env_seas.21[,-10]
env_seas.21<-env_seas.21[,-9]
env_seas.21<-env_seas.21[,-8]
env_seas.21<-env_seas.21[,-3]
env_seas.21<-env_seas.21[,-2]
env_seas.21<-env_seas.21[,-1]
env_seas.21<-env_seas.21[-60,] #remove NatHist Late-no data in 2021

rownames(env_seas.21)<-paste(env_seas.21$Location,"-",env_seas.21$Season)
all.equal(rownames(dat_seas.21), rownames(env_seas.21))

cca_min=cca(dat_seas.21~1, env_seas.21)
cca_full=cca(dat_seas.21~., env_seas.21)

ordistep(cca_min, scope=formula(cca_full))

cca4_tot21=cca(dat_seas.21~Prim_Class+Y+DCC, env_seas.21)
cca4_tot21
RsquareAdj(cca4_tot21)
anova(cca4_tot21)
anova(cca4_tot21, by="terms")
anova(cca4_tot21, by="margin")


#MSO for autocorrelation
xy.21<-cbind(env_seas.21$X,env_seas.21$Y)
mso1.21 <- mso(cca4_tot21, xy.21, grain = 0.1, perm=999)
msoplot(mso1.21,ylim=c(0,1), xlim=c(0,0.5))
mso1.21$vario

#Old version of plot
#plot(cca4_tot21, type='n', scaling=3, ylim=c(-4,4))
#orditorp(cca4_tot21, display='sp', cex=1, scaling=3, col='blue')
#orditorp(cca_tot7,display ="sites",cex=0.5, scaling=2, col='black') ## Orditop includes text for some sites to avoid overcrowding, others are plotted as data points
#text(cca4_tot21, display='cn', col='black')
#ordiellipse(cca4_tot21, groups = env_seas.21$Season, draw = "polygon",lty = 0, col = c("thistle", "slategray3", "lightsalmon1", "palegreen"),label=TRUE)

#Updated figure code
plot(cca4_tot21, type='n', scaling=scl, xlim=c(-8,8), ylim=c(-3, 3))
ordiellipse(cca4_tot21, groups = env_seas.21$Season, draw = "polygon",lty = 0, col = c( "thistle","slategray3","lightsalmon1","palegreen"),label=FALSE)
with(env_seas, points(cca4_tot21, display = "sites", col = colvec[Prim_Class],
                      scaling = scl, pch = 21, cex=0.5, bg = colvec[Prim_Class]))
text(cca4_tot21, display='cn', col='black')
orditorp(cca4_tot21, display='sp', cex=1, scaling=3, col='blue')
with(env_seas.21, legend(3, -1.5, cex=0.8,legend = levels(Prim_Class), y.intersp = 0.2, x.intersp=0.2, bty = "n",
                         col = colvec, pch = 21, pt.bg = colvec))


cca5_tot21=cca(dat_seas.21~PercentImp_Mean+Y+DCC+Sec_Class, env_seas.21)
cca5_tot21
RsquareAdj(cca5_tot21)
anova(cca5_tot21)
anova(cca5_tot21, by="terms")
anova(cca5_tot21, by="margin")
vif.cca(cca5_tot21) 


#Revise plots with site scores
plot(cca5_tot21, type='n', scaling=scl, xlim=c(-8,8), ylim=c(-3, 3))
ordiellipse(cca5_tot21, groups = env_seas.21$Season, draw = "polygon",lty = 0, col = c( "thistle","slategray3","lightsalmon1","palegreen"),label=FALSE)
with(env_seas, points(cca5_tot21, display = "sites", col = colvec2[Sec_Class],
                      scaling = scl, pch = 21, cex=0.5, bg = colvec2[Sec_Class]))
text(cca5_tot21, display='cn', col='black')
orditorp(cca5_tot21, display='sp', cex=1, scaling=3, col='blue')
with(env_seas.21, legend(3, -1.5, cex=0.8,legend = levels(Sec_Class), y.intersp = 0.2, x.intersp=0.2, bty = "n",
                         col = colvec2, pch = 21, pt.bg = colvec2))

#Need data in a matrix for variance partitioning. Extract each column as a vector and bind
pip21<-as.vector(dat_seas.21[,1])
quinq21<-as.vector(dat_seas.21[,2])
rest21<-as.vector(dat_seas.21[,3])

dat.seas.21<-as.matrix(cbind(pip21, quinq21, rest21))

vp21<-varpart(dat.seas.21,~Season,~X+Y, ~DCC+NDVI+ Pop+PercentImp_Mean+WaterTbl+Elev, data=env_seas.21, chisquare = TRUE)
vp21
plot(vp21, Xnames=c("Season", "Space", "Environment"))

plot_vp_euler3(
  vp21,
  names = c("Season", "Space", "Environment"),
  col = c("brown3", "skyblue3", "orange")
)

#Function to summarize data and generate Mean, SE, CI
###create function to summarize
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Correlation plots for env. variables
pd<-summarySE(env_phen,measurevar="DCC", groupvars="Prim_Class")
pd
pi<-summarySE(env_phen,measurevar= "PercentImp_Mean", groupvars="Prim_Class")
pi
px<-summarySE(env_phen,measurevar= "X", groupvars="Prim_Class")
px
py<-summarySE(env_phen,measurevar= "Y",groupvars="Prim_Class")
py
pn<-summarySE(env_phen,measurevar="NDVI",groupvars="Prim_Class")
pn

P.D<-ggplot(pd, aes(x=Prim_Class,y=DCC))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=DCC-ci,ymax=DCC+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  theme_bw()+scale_x_discrete(labels=c("Rural" = "Rur", "Suburban" = "Sub",
                                       "Urban" = "Urb"))+
  labs(x="Site Class", y = "Distance to City Center (km)")+
  theme(legend.position = "none", text = element_text(size = 12))
P.D

P.I<-ggplot(pi, aes(x=Prim_Class,y=PercentImp_Mean))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=PercentImp_Mean-ci,ymax=PercentImp_Mean+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  theme_bw()+scale_x_discrete(labels=c("Rural" = "Rur", "Suburban" = "Sub",
                                       "Urban" = "Urb"))+
  labs(x="Site Class", y = "Mean Percent Impervious Surface (%)")+
  theme(legend.position = "none", text = element_text(size = 12))
P.I

P.X<-ggplot(px, aes(x=Prim_Class,y=X))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=X-ci,ymax=X+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  theme_bw()+scale_x_discrete(labels=c("Rural" = "Rur", "Suburban" = "Sub",
                                       "Urban" = "Urb"))+
  labs(x="Site Class", y = "Longitude")+
  theme(legend.position = "none", text = element_text(size = 12))
P.X

P.Y<-ggplot(py, aes(x=Prim_Class,y=Y))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=Y-ci,ymax=Y+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  theme_bw()+scale_x_discrete(labels=c("Rural" = "Rur", "Suburban" = "Sub",
                                       "Urban" = "Urb"))+
  labs(x="Site Class", y = "Latitude")+
  theme(legend.position = "none", text = element_text(size = 12))
P.Y

P.N<-ggplot(pn, aes(x=Prim_Class,y=NDVI))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=NDVI-ci,ymax=NDVI+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+scale_x_discrete(labels=c("Rural" = "Rur", "Suburban" = "Sub",
                                       "Urban" = "Urb"))+
  labs(x="Site Class", y = "NDVI")+
  theme(legend.position = "none", text = element_text(size = 12))
P.N

sd<-summarySE(env_phen,measurevar="DCC", groupvars="Sec_Class")
sd
si<-summarySE(env_phen,measurevar= "PercentImp_Mean", groupvars="Sec_Class")
si
sx<-summarySE(env_phen,measurevar= "X", groupvars="Sec_Class")
sx
sy<-summarySE(env_phen,measurevar= "Y",groupvars="Sec_Class")
sy
sn<-summarySE(env_phen,measurevar="NDVI",groupvars="Sec_Class")
sn

S.D<-ggplot(sd, aes(x=Sec_Class,y=DCC))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=DCC-ci,ymax=DCC+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  scale_x_discrete(labels=c("Agricultural" = "Agr", "Commercial" = "Com",
                            "Natural" = "Nat", "Residential"="Res"))+
  labs(x="Secondary Class", y = "Distance to City Center (km)")+
  theme(legend.position = "none", text = element_text(size = 12))
S.D

S.I<-ggplot(si, aes(x=Sec_Class,y=PercentImp_Mean))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=PercentImp_Mean-ci,ymax=PercentImp_Mean+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  scale_x_discrete(labels=c("Agricultural" = "Agr", "Commercial" = "Com",
                            "Natural" = "Nat", "Residential"="Res"))+
  labs(x="Secondary Class", y = "Mean Percent Impervious Surface (%)")+
  theme(legend.position = "none", text = element_text(size = 12))
S.I

S.X<-ggplot(sx, aes(x=Sec_Class,y=X))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=X-ci,ymax=X+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  scale_x_discrete(labels=c("Agricultural" = "Agr", "Commercial" = "Com",
                            "Natural" = "Nat", "Residential"="Res"))+
  labs(x="Secondary Class", y = "Longitude")+
  theme(legend.position = "none", text = element_text(size = 12))
S.X

S.Y<-ggplot(sy, aes(x=Sec_Class,y=Y))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=Y-ci,ymax=Y+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  scale_x_discrete(labels=c("Agricultural" = "Agr", "Commercial" = "Com",
                            "Natural" = "Nat", "Residential"="Res"))+
  labs(x="Secondary Class", y = "Latitude")+
  theme(legend.position = "none", text = element_text(size = 12))
S.Y

S.N<-ggplot(sn, aes(x=Sec_Class,y=NDVI))+
  geom_point(size=1, alpha=0.7)   +
  geom_errorbar(aes(ymin=NDVI-ci,ymax=NDVI+ci), width=0.2) +
  theme(text = element_text(size = 15))+      
  theme_bw()+
  scale_x_discrete(labels=c("Agricultural" = "Agr", "Commercial" = "Com",
                            "Natural" = "Nat", "Residential"="Res"))+
  labs(x="Secondary Class", y = "NDVI")+
  theme(legend.position = "none", text = element_text(size = 12))
S.N

library(gridExtra)
grid.arrange(P.D, P.I, P.N, P.X, P.Y, S.D, S.I, S.N, S.X, S.Y, ncol=5, nrow=2)

cxd<-ggplot(env_phen, aes(x=X, y=DCC))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Longitude", y = "DCC")

cyd<-ggplot(env_phen, aes(x=Y, y=DCC))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Latitude", y = "DCC")

cxi<-ggplot(env_phen, aes(x=X, y=PercentImp_Mean))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Longitude", y = "MPIS")

cyi<-ggplot(env_phen, aes(x=Y, y=PercentImp_Mean))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Latitude", y = "MPIS")

cxn<-ggplot(env_phen, aes(x=X, y=NDVI))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Longitude", y = "NDVI")

cyn<-ggplot(env_phen, aes(x=Y, y=NDVI))+
  geom_point()+
  theme(text = element_text(size = 15))+      
  theme_bw()+
  labs(x="Latitude", y = "NDVI")

grid.arrange(cxd, cxi, cxn, cyd, cyi, cyd, ncol=3, nrow=2)


