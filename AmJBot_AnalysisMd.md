# PCA 
Environmental PCA for the populations with all 18 Bioclim variables at 2.5 minutes resolution.

Variables with the highest contribution to PC1:

* bio3 Isothermality (bio2/bio7)*100 = (Mean of monthly (max temp - min temp))/temp annual range
* bio18 precip of warmest quarter
* bio15 precipitation seasonality (coeff of var)
* bio11 mean temp of coldest quarter
* bio9 Mean Temp warmest month

```{r pca, message=F, warning=F,echo=F}
# packages
library(raster)

# read dataframe
pcadf <- read.csv(file = "PNW_EU_AllPops_meta.csv", sep = ",", header = T)
pcadf2 <- read.csv(file = "2022GH_pops.csv", sep = ",", header = T)

# extract coord and get data from worldclim (lon,lat)
coord <- pcadf[,4:3]

r <- raster::getData("worldclim",var = "bio", res = 2.5)
values <- extract(r,coord)

pcabiodf <- cbind(pcadf,values)

write.csv(pcabiodf,file = "pcabiodf.csv",quote = F)

pcabiodf<-read.csv("pcabiodf.csv",header = T,sep = ",")
pcabiodf <- pcabiodf[,2:27]

################# PCA of climate variables
# biodf is dataframe of pops, fullPop, lat/lon,region, country, state and bio1-19
# [,7:25] are bioclim variables

res.pca <- PCA(pcabiodf[8:26])
eig.val <- get_eigenvalue(res.pca)
head(eig.val)

fviz_eig(res.pca, addlabels = T)

# contribution of variables to pc1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 11)
# bio3 Isothermality (bio2/bio7)*100 = (Mean of monthly (max temp - min temp))/temp annual range
# bio18 precip of warmest quarter
# bio15 precipitation seasonality (coeff of var)
# bio11 mean temp of coldest quarter
# bio9 Mean Temp warmest month

####### use all samples PCs to make figure and highlight GH samples
allpca<-res.pca[["ind"]][["coord"]]
allpca <- as.data.frame(allpca)
alldf <- bind_cols(pcabiodf,allpca)
# write.csv(alldf, file = "allpcadf.csv",quote = F)

pcab<-ggplot(data = alldf, aes(x = Dim.1, y=Dim.2,color = state,shape = state)) +
  geom_point(size = 3.5)+
  scale_shape_manual(values=c(0, 1, 3, 4, 5, 25, 7, 8))+
  xlab("PC 1 (55.12%)")+
  ylab("PC 2 (18.84%)")+
  geom_point(data=pcadf2,
             aes(x=Dim.1,y=Dim.2),
             color='red',
             shape = 17,
             size=3)+
  theme_classic()+
  theme(legend.position= "top")

pcab
pcac<-fviz_pca_var(res.pca)+ 
  theme_classic()+
  xlab("PC 1 (55.12%)")+
  ylab("PC 2 (18.84%)")+
  ggtitle(label = "")

pcac

```

# Map of samples in EU and PNW

```{r map, message=F,warning=F,echo=F}


dfEU <- pcadf2 %>% filter(region == "EU")
dfPNW <- pcadf2 %>% filter(region == "US")

world <- ne_countries(scale = "medium", returnclass = "sf")


pnwmap<-ggplot(data = world) +
  geom_sf()+ 
  coord_sf(xlim = c(-126,-118), ylim = c(42,49.3))+
  geom_point(data=dfPNW, 
             aes(x=lon,y=lat), 
             color='blue',
             shape = 19,
             size=4,
             alpha = 0.5)+
  ylab(expression("Latitude ("*degree*")" )) +  
  xlab(expression("Longitude ("*degree*")" )) +
  theme_classic()
# pnwmap

eumap<-ggplot(data = world) +
  geom_sf()+ 
  coord_sf(xlim = c(-8,15), ylim = c(43,62))+
  geom_point(data=dfEU, 
             aes(x=lon,y=lat), 
             color='blue',
             shape = 19,
             size=4,
             alpha = 0.5)+
  ylab(expression("Latitude ("*degree*")" )) +  
  xlab(expression("Longitude ("*degree*")" )) +
  theme_classic()
# eumap

ggdraw() +
  draw_plot(eumap, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(pnwmap, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
```

# Plot of the moisture measurements over the course of the experiment

```{r moisture, message=F,warning=F,echo=F}
moist <- read.csv("moistMeasure2.csv", header = T)

moist2 <- 
  moist%>%
    group_by(trt, day)%>%
  summarize(mean=mean(sample),)

ggplot(moist2, aes(x=day, y=mean, color=trt)) +
  geom_point()+geom_line()+
  theme_classic()+
  scale_color_manual(values = c("#649cfc", "#e67069"))+
  labs(x="Day", y="Mean VCW across 10 plants", color = "Treatment")
```

# Variables

```{r, dataprep, warning=F,message=F, echo=F}
# read in data
df <- read.csv(file = "GH_data_master.csv", header = T)
moist <- read.csv("moistMeasure2.csv", header = T)
meta <- read.csv("2022GH_pops.csv", header = T)
popSumm <- read.csv("popSumGH2.csv", header = T)

# manipulating df
# closedCap + openCap = totalCap
df <- df %>% mutate(totalCap = closedCap + openCap)

# estSeed = 
df <- df %>% mutate(estSeed = (numSeed/numSeedCap)*totalCap)

# closedCap + openCap = totalCap
df <- df %>% mutate(totalCap = closedCap + openCap)

# estSeed = 
df <- df %>% mutate(estSeed = (numSeed/numSeedCap)*totalCap)

# remove nigrescens samples
df2 <- df %>%filter(nigrescens == 0)

# remove unbolted ind
df_bolt <- df2 %>%filter(bolted == 1)

# remove nigrescens samples
df2 <- df %>%filter(nigrescens == 0)

###time from bolting to flowering
df2$boltToFlower <- df2$flowerExpDays - df2$boltExpDays
```

## Transforming variables

#### First looking at the variables themselves, checking mean and sample size of the two treatments/regions and transforming variables to improve normality.

* Log transformed: dry weight, stem width, Tleaf.

* Square root transformed: stem height, total capitula, totalCap, gsw

* Didn't need/want to transform: boltToFlower, estSeed

* Didn't work when I transformed (tried both log and sqrt): RWC, boltExpDays, flowerExpDays, PhiPS2

### Dry weight with all plants

```{r dryweight, warning=F, message=F,echo=F}
#dryWeight
temp <- df2 %>% 
  filter_at(vars(dryWeight), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(dryWeight), n = n())

# hist(temp$dryWeight)

shapiro.test(temp$dryWeight)

## Transform

df2$dryWeightLog <- log10(df2$dryWeight)

temp <- df2 %>% 
  filter_at(vars(dryWeightLog), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(dryWeightLog), n = n())

hist(temp$dryWeightLog)

shapiro.test(temp$dryWeightLog)
```

### Dry weight with only bolted plants

```{r,bolted,message=F,warning=F,echo=F}
#dryWeight of only bolted
tempb <- df_bolt %>% 
  filter_at(vars(dryWeight), all_vars(!is.na(.)))
tempb %>%
  group_by(trt) %>%
  summarise(mean = mean(dryWeight), n = n())

# hist(temp$dryWeight)

shapiro.test(temp$dryWeight)

## Transform

df_bolt$dryWeightLog <- log(df_bolt$dryWeight)

temp <- df_bolt %>% 
  filter_at(vars(dryWeightLog), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(dryWeightLog), n = n())

hist(temp$dryWeightLog)

shapiro.test(temp$dryWeightLog)
```

### Stem height and width

```{r, stem, warning=F, message=F,echo=F}
#stemHeight
temp <- df2 %>% 
  filter_at(vars(stemHeight), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(stemHeight), n = n())
# hist(temp$stemHeight)

shapiro.test(temp$stemHeight)

## Transform

df2$stemHeightsqrt <- sqrt(df2$stemHeight)

temp <- df2 %>% 
  filter_at(vars(stemHeightsqrt), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(stemHeightsqrt), n = n())

hist(temp$stemHeightsqrt)

shapiro.test(temp$stemHeightsqrt)

#stemWidth
temp <- df2 %>% 
  filter_at(vars(stemWidth), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(stemWidth), n = n())
# hist(temp$stemWidth)

shapiro.test(temp$stemWidth)

## Transform

df2$stemWidthLog <- log10(df2$stemWidth)

temp <- df2 %>% 
  filter_at(vars(stemWidthLog), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(stemWidthLog), n = n())

hist(temp$stemWidthLog)

shapiro.test(temp$stemWidthLog)
```

### RWC

```{r, rwc, message=F, warning=F,echo=F}
#RWC
temp <- df2 %>% 
  filter_at(vars(RWC), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(RWC), n = n())

hist(temp$RWC)

shapiro.test(temp$RWC)

shapiro.test(sqrt(temp$RWC))

## Transform
# 
# df2$ <- sqrt(df2$stemHeight)
# 
# temp <- df2 %>% 
#   filter_at(vars(stemHeightsqrt), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(stemHeightsqrt), n = n())
# 
# hist(temp$stemHeightsqrt)
# 
# shapiro.test(temp$stemHeightsqrt)
```

### bolting and flowering in days since the beginning of the experiment

```{r, boltflower, message=F,warning=F,echo=F}
#boltExpDays
temp <- df2 %>% 
  filter_at(vars(boltExpDays), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(boltExpDays), n = n())
hist(temp$boltExpDays)

shapiro.test(temp$boltExpDays)
shapiro.test(sqrt(temp$boltExpDays))

#flowerExpDays
temp <- df2 %>% 
  filter_at(vars(flowerExpDays), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(flowerExpDays), n = n())
hist(temp$flowerExpDays)

shapiro.test(temp$flowerExpDays)
shapiro.test(sqrt(temp$flowerExpDays))
```


### Time from bolting to flowering in days since start of the experiment

```{r,timebolttoflower,message=F,warning=F,echo=F}
# time from bolting to flowering
#temp <- df2 %>% 
  #filter_at(vars(boltToFlower), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(boltToFlower), n = n())
hist(temp$boltToFlower)

shapiro.test(temp$boltToFlower)
```

### Total capitula and estimated seed

```{r, seedcap, message=F,warning=F,echo=F}
#totalCap
temp <- df2 %>% 
  filter_at(vars(totalCap), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(totalCap), n = n())

hist(temp$totalCap)

shapiro.test(temp$totalCap)
shapiro.test(sqrt(temp$totalCap))

## Transform

df2$totalCapsqrt <- sqrt(df2$totalCap)

temp <- df2 %>% 
  filter_at(vars(totalCapsqrt), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(stemHeightsqrt), n = n())

hist(temp$totalCapsqrt)

shapiro.test(temp$totalCapsqrt)

#estSeed
temp <- df2 %>% 
  filter_at(vars(estSeed), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(estSeed), n = n())

hist(temp$estSeed)
shapiro.test(temp$estSeed)

```

### Stomatal Conductance (gsw)

```{r, gsw, warning=F,message=F,echo=F}
########## licor data
licor1 <- read.csv("LiCordf.csv", header = T, sep = ",")

licor <- left_join(licor1,df,by = "blockID")

# view(licor)

licor <- licor %>% dplyr::select(indID,mat,pop,lat,lon,
                          state,region,pc1,trt,blockNum,
                          blockID,nigrescens,bolted,Date,
                          gsw,Fs,PhiPS2,Tref,Tleaf) %>%
  filter(nigrescens == 0)# remove nigrescens samples
####################

#gsw
temp <- licor %>% 
  filter_at(vars(gsw), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(gsw), n = n())

hist(temp$gsw)

shapiro.test(temp$gsw)

## Transform

licor$gswsqrt <- sqrt(licor$gsw)

temp <- licor %>% 
  filter_at(vars(gswsqrt), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(dryWeightLog), n = n())

hist(temp$gswsqrt)

shapiro.test(temp$gswsqrt)
```

### Fluorescence (phips2)

```{r, phips2, warning=F,message=F,echo=F}
#gsw
temp <- licor %>% 
  filter_at(vars(PhiPS2), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(PhiPS2), n = n())

hist(temp$PhiPS2)

shapiro.test(temp$PhiPS2)

## Transform
# 
# licor$PhiPS2sqrt <- sqrt(max(temp$PhiPS2+1) - temp$PhiPS2)
# 
# temp <- licor %>% 
#   filter_at(vars(PhiPS2sqrt), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(dryWeightLog), n = n())

# hist(temp$PhiPS2sqrt)
# 
# shapiro.test(temp$PhiPS2sqrt)
```

### Temperature of Leaf (Tleaf)

```{r, tleaf, warning=F,message=F,echo=F}
#Tleaf
temp <- licor %>% 
  filter_at(vars(Tleaf), all_vars(!is.na(.)))
temp %>%
  group_by(trt) %>%
  summarise(mean = mean(Tleaf), n = n())

hist(temp$Tleaf)

shapiro.test(temp$Tleaf)

## Transform

licor$TleafLog <- log10(licor$Tleaf)

temp <- licor %>% 
  filter_at(vars(TleafLog), all_vars(!is.na(.)))
# temp %>%
#   group_by(trt) %>%
#   summarise(mean = mean(dryWeightLog), n = n())

hist(temp$TleafLog)

shapiro.test(temp$TleafLog)
```


# Models
```{r,models,warning=F,message=F, echo=F}

```

### Bolted vs unbolted

use glmer, but family = binomial

formula:

##### trait ~ trt + region + trt*region + 1|pop + 1|mat + 1|blockNum

* pc1 is significant when included. When pc1 isn't included, region is significant (p=0.023).

* Not including pc1 because it is confounded with region

```{r,boltedglmer,warning=F,message=F}
# df2$bolted <- as.factor(df2$bolted)

m1 <- glmer(bolted ~ trt + region + trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, family = "binomial")

summary(m1)
tab_model(m1)

a1<-plot_model(m1, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Percent Bolted",title = "Probability of Bolting")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
a1

# m2 <- glmer(bolted ~ pc1 + trt + trt*pc1 + (1|pop) + (1|mat) + (1|blockNum), data=df2, family = "binomial")
# 
# summary(m2)
# tab_model(m2)

# a2<-plot_model(m2, type = "int")+
#   theme_classic()+geom_line(position = position_dodge(width=.1))+
#   labs(x="PC1",y="bolted")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
# a2
```

### Main variables

For the variables: dryWeightLog, stemWidthLog, stemHeightsqrt, totCapsqrt, RWC, boltExpDays, flowerExpDays, and boltToFlower, use formula:

##### trait ~ pc1 + trt + region + trt*region + 1|pop + 1|mat + 1|blockNum

Growth traits:

* Dry Weight: Treatment, region, and interaction are significant. US is higher than EU and Control is higher than Treatment.

* Stem Width: Treatment is significant. Control is higher than Treatment.

* Stem Height: Treatment is significant. Control is higher than Treatment.

Physiology Traits:

* RWC: No significance.

* Stomatal Conductance: Region and region*trt significant.

* Fluorescence: Region and region*trt significant.

* Temperature of Leaf: Region and region*trt significant.

Reproduction Traits:

* Total Capitula: Region is significant. US is higher than EU.

* Bolt since exp start: No significance.

* Flowering since exp start: No significance.

* Bolting to flowering: No significance.


### dryWeightLog

```{r,dryWeightLogglmm,warning=F,message=F}
mod1 = lmer(dryWeightLog~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod1)
tab_model(mod1)

# plot_model(mod1, type = "pred")
# plot_model(mod1, type = "re")

a<-plot_model(mod1, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Dry Weight (g)", title = "Dry Weight")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
a
```

### stemWidthLog

```{r,stemWidthLogglmm,warning=F,message=F}
mod2 = lmer(stemWidthLog~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod2)
tab_model(mod2)


# plot_model(mod2, type = "pred")
# plot_model(mod2, type = "re")

b<-plot_model(mod2, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Stem Width (mm)", title = "Stem Width")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
b
```

### stemHeightsqrt

```{r,stemHeightsqrtglmm,warning=F,message=F}
mod3 = lmer(stemHeightsqrt~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod3)
tab_model(mod3)


# plot_model(mod3, type = "pred")
# plot_model(mod3, type = "re")

c<-plot_model(mod3, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Stem Height (cm)", title = "Stem Height")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
c
```

### totCapsqrt

```{r,totCapsqrtglmm,warning=F,message=F}
mod4 = lmer(totalCapsqrt~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod4)
tab_model(mod4)


# plot_model(mod4, type = "pred")
# plot_model(mod4, type = "re")

d<-plot_model(mod4, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Total Capitula", title = "Total Capitula")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
d
```

### RWC

```{r,RWCglmm,warning=F,message=F}
mod5 = lmer(RWC~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod5)
tab_model(mod5)


# plot_model(mod5, type = "pred")
# plot_model(mod5, type = "re")

e<-plot_model(mod5, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="RWC", title = "RWC")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
e
```

### boltExpDays

```{r,boltExpDaysglmm,warning=F,message=F}
mod6 = lmer(boltExpDays~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod6)
tab_model(mod6)


# plot_model(mod6, type = "pred")
# plot_model(mod6, type = "re")

f<-plot_model(mod6, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Bolting Day", title = "Bolting Day")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
f
```

### flowerExpDays

```{r,flowerExpDaysglmm,warning=F,message=F}
mod7 = lmer(flowerExpDays~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod7)
tab_model(mod7)


# plot_model(mod7, type = "pred")
# plot_model(mod7, type = "re")

g<-plot_model(mod7, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Flowering Day", title = "Flowering Day")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
g
```

### boltToFlower

```{r,boltToFlowerglmm,warning=F,message=F}
mod8 = lmer(boltToFlower~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=df2, REML = F)
summary(mod8)
tab_model(mod8)

# plot_model(mod8, type = "pred")
# plot_model(mod8, type = "re")

h<-plot_model(mod8, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Days from Bolting to Flowering", title = "Days from Bolting to Flowering")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
h
```

### Stomatal Conductance (gsw)

```{r, gswlmm, warning=F,message=F}
mod9 = lmer(gswsqrt~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=licor, REML = F)
summary(mod9)
tab_model(mod9)

# plot_model(mod8, type = "pred")
# plot_model(mod8, type = "re")

i<-plot_model(mod9, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="gsw", title = "Stomatal Conductance (gsw)")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
i

```

### Fluorescence (phips2)

```{r, phips2lmm, warning=F,message=F}
mod10 = lmer(PhiPS2~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=licor, REML = F)
summary(mod10)
tab_model(mod10)

# plot_model(mod8, type = "pred")
# plot_model(mod8, type = "re")

j<-plot_model(mod10, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="PhiPS2", title = "Fluorescence")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
j
```

### Temperature of Leaf (Tleaf)

```{r, tleaflmm, warning=F,message=F}
mod11 = lmer(TleafLog~trt*region + (1|pop) + (1|mat) + (1|blockNum), data=licor, REML = F)
summary(mod11)
tab_model(mod11)

# plot_model(mod8, type = "pred")
# plot_model(mod8, type = "re")

k<-plot_model(mod11, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Tleaf", title = "Leaf Temperature (Tleaf)")
# uncomment below to get adjusted means
# +stat_summary(aes(label=..y..), fun=mean, geom="text", size=2)
  
k
```



## Estimated Seed

Formula:

trait ~ pc1 + trt + region + trt*region + 1|pop + 1|mat + 1|blockNum

use zeroinfl

Note to self: add back in NaNs as 0s

### estSeed

```{r,estSeedzeroinfl,warning=F,message=F}
# zero1 <- zeroinfl(estSeed ~ pc1 + trt + region + trt*region + (1|pop) + (1|mat) + (1|blockNum), data = df2)



seedMod <- glmmTMB(estSeed ~ trt + region + trt*region + (1|pop) + (1|mat) + (1|blockNum), 
                   data = df2,
                   family = nbinom2,
                   ziformula = ~trt*region)
summary(seedMod)
tab_model(seedMod)
plot_model(seedMod, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Estimated Seed", title = "Predicted Values of Estimated Seed")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))
```

### estSeed but with the bolted 0 seed ones set to 0 instead of NaN

```{r, estSeed2,message=F,warning=F}
# change the NaNs in the estimated seed column to 0
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

df2[is.nan(df2)] <- 0

seedMod2 <- glmmTMB(estSeed ~ trt + region + trt*region + (1|pop) + (1|mat) + (1|blockNum), 
                   data = df2,
                   family = nbinom2,
                   ziformula = ~trt*region)
summary(seedMod2)
tab_model(seedMod2)
seedplot<-plot_model(seedMod2, type = "int")+
  theme_classic()+geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="Estimated Seed", title = "Estimated Seed")+scale_color_manual(values = c("#577eb8","#e41a1c"))

seedplot
```

#### figures for GLMMs

```{r figures, message=F,warning=F}
### add asterisks for significance?

# Growth traits
#a # Dry Weight
#b # Stem Width
#c # Stem Height

ggarrange(a, b, c + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2,
          common.legend = T)
# Physiology
#e # RWC
#i # gsw
#j # PhiPS2
#k # Tleaf
ggarrange(e, i, j,k + rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          common.legend = T)

# Reproduction
#d # Total Cap
#f # Bolting Day
#g # Flowering Day
#h # Days from Bolting to Flowering
#i # Estimated Seed
ggarrange(f, g, h, d,seedplot + rremove("x.text"), 
          labels = c("A", "B", "C","D","E"),
          ncol = 2, nrow = 3,
          common.legend = T)

```

## Repeated Measures ANOVA

Use for the SPAD data in spad.csv, need to transform to long format

Use formula with additional fixed effect for week and additional random effect for individual:

##### trait ~ week + trt + region + trt*region week + 1|pop + 1|mat + 1|blockNum + 1|ind

### SPAD Data (absorbance)

```{r SPAD, message=F, warning=F}


spad <- read.csv("spad.csv", header = T, sep = ",")

# remove nigrescens samples
spaddf <- spad %>%filter(nigrescens == 0)
# head(spaddf)
# make into long format

spaddf2<-spaddf %>%
  pivot_longer(
    cols = starts_with("chlor_"),
    names_to = "week",
    names_prefix = "chlor_",
    values_to = "spad",
    values_drop_na = TRUE)

head(spaddf2)
spaddf2$week <- as.integer(spaddf2$week)

# spaddf2 %>%
#   group_by(indID)%>%
#   get_summary_stats(spad, type = "mean_sd")

# ggqqplot(spaddf2, "spad", facet.by = "week")

aov <- lmer(spad ~ week + trt + region + trt * region + week*trt + week*region + week*trt*region + (1 | pop) + (1 | mat) + (1 | blockNum) + (1 | indID), data = spaddf2, REML = F)

summary(aov)
tab_model(aov)



aovplot <- plot_model(aov, type = "pred", terms=c("trt","region"))+
  theme_classic()+
  geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="SPAD Measurements", title = "Absorbance")
aovplot

aovplot2<- plot_model(aov, type = "pred", terms=c("week","trt"))+
  theme_classic()+
  geom_line(position = position_dodge(width=.1))+
  labs(x="Week",y="SPAD Measurements", title = "")
aovplot2

plot_model(aov, type = "pred", terms=c("trt","week"))+
  theme_classic()+
  geom_line(position = position_dodge(width=.1))+
  labs(x="Treatment",y="SPAD Measurements")

plot_model(aov, type = "pred", terms=c("trt"))+
  theme_classic()+
  geom_line(position = position_dodge(width=.1))+
  labs(x="Region",y="SPAD Measurements")

plot_model(aov, type = "pred", terms=c("region","week"))+
  theme_classic()+
  geom_line(position = position_dodge(width=.1))+
  labs(x="Region",y="SPAD Measurements")

# aovdat <- ggpredict(aov, terms = c("trt","week"))
# plot(aovdat, facet=F)
# 
# 
# ggplot(aovdat, aes(x=trt,y=))
# 
# aovplot2 <- ggplot(spaddf2, aes(x=trt, y=spad, fill = region))+
#   geom_violin()+theme_classic()
# aovplot2
```

```{r figures2, message=F,warning=F}
### SPAD figures

ggarrange(aovplot, aovplot2 + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
```

### deltaTrait

* Obtain BLUPs by using trait ~ trt + (trt|pop) + (1|mat) + (blockNum)

* lme's coef() function to extract BLUP for each pop at each treatment level

* Take the proportion of those two values for each pop to get deltaTrait

* deltaTrait ~ pc1 + region pc1*region

#### Dry Weight

```{r blupsDryWeight, message=F,warning=F}
# df2 for the first set of traits
# get info for the pops first
meta <- df2 %>%
  group_by(pop) %>%
  summarize(region = first(region),
            state = first(state),
            pc1 = first(pc1),
            lat = first(lat),
            lon = first(lon),
            n = n())

## dryweight
tempb<-lmer(dryWeightLog ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

dwBLUP<-coef(tempb)$pop

dwBLUP$pop <- rownames(dwBLUP)

names(dwBLUP)[names(dwBLUP) == "(Intercept)"] <- "trtC"
names(dwBLUP)[names(dwBLUP) == "trtT"] <- "trtTcoef"

dwBLUP <- dwBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaDryWeight = trtT/trtC)

# join with pop data
dwBLUP <- right_join(dwBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(dwBLUP$deltaDryWeight)$out

# remove SP5
dwBLUP <- dwBLUP %>%
  filter(pop != "SP5")

# linear model with all
lm1a <- lm(deltaDryWeight ~ pc1*region, data = dwBLUP)
lm1b <- lm(deltaDryWeight ~ 1, data = dwBLUP)
lm1c <- lm(deltaDryWeight ~ region, data = dwBLUP)
lm1d1 <- lm(deltaDryWeight ~ pc1, data = dwBLUP)
lm1d2 <- lm(deltaDryWeight ~ pc1+region, data = dwBLUP)

AIC(lm1a,lm1b,lm1c,lm1d1,lm1d2)

# summary(lm1)

delta1 <- ggplot(dwBLUP, aes(x = pc1, 
                             y = deltaDryWeight,
                             color=region,
                             shape = region,
                             group=1))+theme_classic()+
  labs(title = "Dry Weight", x = "PC1", 
       y="Dry Weight DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",
              # size = 0.5)+
  geom_point()+
  annotate("text",x=4,y=0.97, label="Model D1", size = 4,
           fontface="italic")+
  # annotate("text",x=3.5,y=0.968,label="ANOVA: . p = 0.052",
           # size = 3)+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta1
  
# # EU only
# dwBLUP_EU<- dwBLUP %>% filter(region == "EU")
#   
# lm1EU <- lm(deltaDryWeight ~ pc1, data = dwBLUP_EU)
# summary(lm1EU)
# 
# ggplot(dwBLUP_EU, aes(x = pc1, y = deltaDryWeight, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# dwBLUP_US<- dwBLUP %>% filter(region == "US")
# 
# lm1US <- lm(deltaDryWeight ~ pc1, data = dwBLUP_US)
# summary(lm1US)
# 
# ggplot(dwBLUP_US, aes(x = pc1, y = deltaDryWeight, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaDryWeight ~ region+pc1, data=dwBLUP))
```

#### Dry Weight with unbolted removed

```{r blupsDryWeight2, message=F,warning=F}
# df2 for the first set of traits
# get info for the pops first
# meta2 <- df2 %>%
#   dplyr::filter(bolted == 1)%>%
#   group_by(pop) %>%
#   summarize(region = first(region),
#             state = first(state),
#             pc1 = first(pc1),
#             lat = first(lat),
#             lon = first(lon),
#             n = n())
# 
# 
# df3 <- df2 %>%
#   dplyr::filter(bolted == 1)
# 
# ## dryweight
# tempb<-lmer(dryWeightLog ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df3, REML = F)
# 
# dwBLUP2<-coef(tempb)$pop
# 
# dwBLUP2$pop <- rownames(dwBLUP2)
# 
# names(dwBLUP2)[names(dwBLUP2) == "(Intercept)"] <- "trtC"
# names(dwBLUP2)[names(dwBLUP2) == "trtT"] <- "trtTcoef"
# 
# dwBLUP2 <- dwBLUP2 %>%
#   dplyr::select(pop,trtC,trtTcoef) %>%
#   mutate(trtT = trtC + trtTcoef) %>%
#   mutate(deltaDryWeight = trtT/trtC)
# 
# # join with pop data
# dwBLUP2 <- right_join(dwBLUP2,meta2,by = "pop")
# 
# # check for outliers
# boxplot.stats(dwBLUP2$deltaDryWeight)$out
# 
# # remove SP5
# dwBLUP2 <- dwBLUP2 %>% 
#   filter(pop != "SP5")
# 
# # linear model with all
# lm1bolt <- lm(deltaDryWeight ~ pc1, data = dwBLUP2)
# summary(lm1bolt)
# 
# lm1bolt <- lm(deltaDryWeight ~ pc1+region, data = dwBLUP2)
# summary(lm1bolt)
# 
# delta1bolt <- ggplot(dwBLUP2, aes(x = pc1, 
#                              y = deltaDryWeight,
#                              color=region,
#                              shape = region,
#                              group=1))+theme_classic()+
#   labs(title = "Dry Weight Bolted Only", x = "PC1", 
#        y="Dry Weight Drought Ratio")+
#   geom_point()+
#   scale_color_manual(values = c("#e41a1c","#577eb8"))
# 
# delta1bolt
#   
# # EU only
# dwBLUP_EU2<- dwBLUP2 %>% filter(region == "EU")
#   
# lm1EU2 <- lm(deltaDryWeight ~ pc1, data = dwBLUP_EU2)
# summary(lm1EU2)
# 
# ggplot(dwBLUP_EU2, aes(x = pc1, y = deltaDryWeight, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# dwBLUP_US2<- dwBLUP2 %>% filter(region == "US")
# 
# lm1US2 <- lm(deltaDryWeight ~ pc1, data = dwBLUP_US2)
# summary(lm1US2)
# 
# ggplot(dwBLUP_US2, aes(x = pc1, y = deltaDryWeight, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaDryWeight ~ region, data=dwBLUP2))
```

#### Stem Width

```{r blupStemWidth, message=F,warning=F}

## stemWidth
tempb<-lmer(stemWidthLog ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

swBLUP<-coef(tempb)$pop

swBLUP$pop <- rownames(swBLUP)

names(swBLUP)[names(swBLUP) == "(Intercept)"] <- "trtC"
names(swBLUP)[names(swBLUP) == "trtT"] <- "trtTcoef"

swBLUP <- swBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaStemWidth = trtT/trtC)

# join with pop data
swBLUP <- right_join(swBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(swBLUP$deltaStemWidth)$out

# remove SP5
swBLUP <- swBLUP %>% 
  filter(pop != "SP5")

# linear model with all

lm2a <- lm(deltaStemWidth ~ pc1*region, data = swBLUP)
lm2b <- lm(deltaStemWidth ~ 1, data = swBLUP)
lm2c <- lm(deltaStemWidth ~ region, data = swBLUP)
lm2d1 <- lm(deltaStemWidth ~ pc1, data = swBLUP)
lm2d2 <- lm(deltaStemWidth ~ pc1+region, data = swBLUP)

AIC(lm2a,lm2b,lm2c,lm2d1,lm2d2)

# lm2 <- lm(deltaStemWidth ~ pc1+region, data = swBLUP)
# summary(lm2)

delta2 <- ggplot(swBLUP, aes(x = pc1, y = deltaStemWidth, 
                             color=region, 
                             shape = region,
                             group = 1))+
  theme_classic()+
  labs(title = "Stem Width", x = "PC1", 
       y="Stem Width DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",size = 0.5)+
  geom_point()+
  # annotate("text",x=3.5,y=0.86, label="lm: * p=0.025", size = 3)+
  annotate("text",x=4,y=0.865,label="Model D1",
           size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta2

# EU only
# swBLUP_EU<- swBLUP %>% filter(region == "EU")
#   
# lm2EU <- lm(deltaStemWidth ~ pc1, data = swBLUP_EU)
# summary(lm2EU)
# 
# ggplot(swBLUP_EU, aes(x = pc1, y = deltaStemWidth, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# swBLUP_US<- swBLUP %>% filter(region == "US")
# 
# lm2US <- lm(deltaStemWidth ~ pc1, data = swBLUP_US)
# summary(lm2US)
# 
# ggplot(swBLUP_US, aes(x = pc1, y = deltaStemWidth, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaStemWidth ~ region, data=swBLUP))
```


#### Stem Height

```{r blupStemHeight, message=F,warning=F}

## stemWidth
tempb<-lmer(stemHeightsqrt ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

shBLUP<-coef(tempb)$pop

shBLUP$pop <- rownames(shBLUP)

names(shBLUP)[names(shBLUP) == "(Intercept)"] <- "trtC"
names(shBLUP)[names(shBLUP) == "trtT"] <- "trtTcoef"

shBLUP <- shBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaStemHeight = trtT/trtC)

# join with pop data
shBLUP <- right_join(shBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(shBLUP$deltaStemHeight)$out

# remove SP5 and FR3
shBLUP <- shBLUP %>% 
  filter(pop != "SP5")%>%
  filter(pop != "FR3")

# linear model with all
lm3a <- lm(deltaStemHeight ~ pc1*region, data = shBLUP)
lm3b <- lm(deltaStemHeight ~ 1, data = shBLUP)
lm3c <- lm(deltaStemHeight ~ region, data = shBLUP)
lm3d1 <- lm(deltaStemHeight ~ pc1, data = shBLUP)
lm3d2 <- lm(deltaStemHeight ~ pc1+region, data = shBLUP)

AIC(lm3a,lm3b,lm3c,lm3d1,lm3d2)

# summary(lm3)

delta3 <- ggplot(shBLUP, aes(x = pc1, y = deltaStemHeight,
                             color=region, 
                             shape = region,
                             group = 1))+
  theme_classic()+
  labs(title = "Stem Height", x = "PC1", 
       y="Drought Ratio")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",size = 0.5)+
  geom_point()+
  # annotate("text",x=3.5,y=0.835, label="lm: ** p=0.002", 
           # size = 3)+
  # annotate("text",x=4,y=0.835,label="Model D1",
           # size = 4,
           # fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta3

# # EU only
# shBLUP_EU<- shBLUP %>% filter(region == "EU")
#   
# lm3EU <- lm(deltaStemHeight ~ pc1, data = shBLUP_EU)
# summary(lm3EU)
# 
# ggplot(shBLUP_EU, aes(x = pc1, y = deltaStemHeight, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# shBLUP_US<- shBLUP %>% filter(region == "US")
# 
# lm3US <- lm(deltaStemHeight ~ pc1, data = shBLUP_US)
# summary(lm3US)
# 
# ggplot(shBLUP_US, aes(x = pc1, y = deltaStemHeight, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaStemHeight ~ region, data=shBLUP))
```

#### RWC

```{r blupRWC, message=F,warning=F}

## stemWidth
tempb<-lmer(RWC ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

rwcBLUP<-coef(tempb)$pop

rwcBLUP$pop <- rownames(rwcBLUP)

names(rwcBLUP)[names(rwcBLUP) == "(Intercept)"] <- "trtC"
names(rwcBLUP)[names(rwcBLUP) == "trtT"] <- "trtTcoef"

rwcBLUP <- rwcBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaRWC = trtT/trtC)

# join with pop data
rwcBLUP <- right_join(rwcBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(rwcBLUP$deltaRWC)$out

# remove HO
rwcBLUP <- rwcBLUP %>% 
  filter(pop != "HO")

# linear model with all
lm4a <- lm(deltaRWC ~ pc1*region, data = rwcBLUP)
lm4b <- lm(deltaRWC ~ 1, data = rwcBLUP)
lm4c <- lm(deltaRWC ~ region, data = rwcBLUP)
lm4d1 <- lm(deltaRWC ~ pc1, data = rwcBLUP)
lm4d2 <- lm(deltaRWC ~ pc1+region, data = rwcBLUP)

AIC(lm4a,lm4b,lm4c,lm4d1,lm4d2)

delta4 <- ggplot(rwcBLUP, aes(x = pc1, 
                              y = deltaRWC, 
                              color=region,
                              shape = region))+
  geom_point()+theme_classic()+
  labs(title = "RWC", x = "PC1", y="RWC DR")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))+
  annotate("text",x=3,y=1.06,label="Model B",
           size = 4,
           fontface="italic")

delta4
  
# # EU only
# rwcBLUP_EU<- rwcBLUP %>% filter(region == "EU")
#   
# lm4EU <- lm(deltaRWC ~ pc1, data = rwcBLUP_EU)
# summary(lm4EU)
# 
# ggplot(rwcBLUP_EU, aes(x = pc1, y = deltaRWC, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# rwcBLUP_US<- rwcBLUP %>% filter(region == "US")
# 
# lm4US <- lm(deltaRWC ~ pc1, data = rwcBLUP_US)
# summary(lm4US)
# 
# ggplot(rwcBLUP_US, aes(x = pc1, y = deltaRWC, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaRWC ~ region, data=rwcBLUP))
```

#### Total Capitula

```{r bluptotcap, message=F,warning=F}

## stemWidth
tempb<-lmer(totalCapsqrt ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

tcBLUP<-coef(tempb)$pop

tcBLUP$pop <- rownames(tcBLUP)

names(tcBLUP)[names(tcBLUP) == "(Intercept)"] <- "trtC"
names(tcBLUP)[names(tcBLUP) == "trtT"] <- "trtTcoef"

tcBLUP <- tcBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaTotCap = trtT/trtC)

# join with pop data
tcBLUP <- right_join(tcBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(tcBLUP$deltaTotCap)$out

# linear model with all
lm5a <- lm(deltaTotCap ~ pc1*region, data = tcBLUP)
lm5b <- lm(deltaTotCap ~ 1, data = tcBLUP)
lm5c <- lm(deltaTotCap ~ region, data = tcBLUP)
lm5d1 <- lm(deltaTotCap ~ pc1, data = tcBLUP)
lm5d2 <- lm(deltaTotCap ~ pc1+region, data = tcBLUP)

AIC(lm5a,lm5b,lm5c,lm5d1,lm5d2)

delta5 <- ggplot(tcBLUP, aes(x = pc1, y = deltaTotCap, 
                             color=region, 
                             shape = region,
                             group = 1))+
  theme_classic()+
  labs(title = "Total Capitula", x = "PC1",
       y="Total Capitula DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",
              # size = 0.5)+
  geom_point()+
  # annotate("text",x=3.5,y=0.87, label="lm: * p=0.047", 
           # size = 3)+
  annotate("text",x=4,y=0.865, label="Model C", 
           size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta5

# # EU only
# tcBLUP_EU<- tcBLUP %>% filter(region == "EU")
#   
# lm5EU <- lm(deltaTotCap ~ pc1, data = tcBLUP_EU)
# summary(lm5EU)
# 
# ggplot(tcBLUP_EU, aes(x = pc1, y = deltaTotCap, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# tcBLUP_US<- tcBLUP %>% filter(region == "US")
# 
# lm5US <- lm(deltaTotCap ~ pc1, data = tcBLUP_US)
# summary(lm5US)
# 
# ggplot(tcBLUP_US, aes(x = pc1, y = deltaTotCap, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaTotCap ~ region, data=tcBLUP))
```

#### Bolting Day

```{r blupbd, message=F,warning=F}

## stemWidth
tempb<-lmer(boltExpDays ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

bdBLUP<-coef(tempb)$pop

bdBLUP$pop <- rownames(bdBLUP)

names(bdBLUP)[names(bdBLUP) == "(Intercept)"] <- "trtC"
names(bdBLUP)[names(bdBLUP) == "trtT"] <- "trtTcoef"

bdBLUP <- bdBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaBolt = trtT/trtC)

# join with pop data
bdBLUP <- right_join(bdBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(bdBLUP$deltaBolt)$out

# remove TR
bdBLUP <- bdBLUP %>% 
  # filter(pop != "TR")%>% 
  filter(pop != "HL")

# linear model with all
lm6a <- lm(deltaBolt ~ pc1*region, data = bdBLUP)
lm6b <- lm(deltaBolt ~ 1, data = bdBLUP)
lm6c <- lm(deltaBolt ~ region, data = bdBLUP)
lm6d1 <- lm(deltaBolt ~ pc1, data = bdBLUP)
lm6d2 <- lm(deltaBolt ~ pc1+region, data = bdBLUP)

AIC(lm6a,lm6b,lm6c,lm6d1,lm6d2)

delta6 <- ggplot(bdBLUP, aes(x = pc1, y = deltaBolt, color=region))+
  geom_point()+theme_classic()+
  labs(title = "Bolting Day", x = "PC1",
       y="Bolting Day DR")+
  annotate("text",x=4,y=1.2, label="Model B", 
           size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta6

# # EU only
# bdBLUP_EU<- bdBLUP %>% filter(region == "EU")
#   
# lm6EU <- lm(deltaBolt ~ pc1, data = bdBLUP_EU)
# summary(lm6EU)
# 
# ggplot(bdBLUP_EU, aes(x = pc1, y = deltaBolt, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# bdBLUP_US<- bdBLUP %>% filter(region == "US")
# 
# lm6US <- lm(deltaBolt ~ pc1, data = bdBLUP_US)
# summary(lm6US)
# 
# ggplot(bdBLUP_US, aes(x = pc1, y = deltaBolt, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaBolt ~ region, data= bdBLUP))
```

#### Flower Day

```{r blupfd, message=F,warning=F}

## stemWidth
tempb<-lmer(flowerExpDays ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

fdBLUP<-coef(tempb)$pop

fdBLUP$pop <- rownames(fdBLUP)

names(fdBLUP)[names(fdBLUP) == "(Intercept)"] <- "trtC"
names(fdBLUP)[names(fdBLUP) == "trtT"] <- "trtTcoef"

fdBLUP <- fdBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaFlow = trtT/trtC)

# join with pop data
fdBLUP <- right_join(fdBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(fdBLUP$deltaFlow)$out

# remove SP4 and HL
# fdBLUP <- fdBLUP %>% 
#   filter(pop != "SP4")%>% 
#   filter(pop != "HL")

# linear model with all
lm7a <- lm(deltaFlow ~ pc1*region, data = fdBLUP)
lm7b <- lm(deltaFlow ~ 1, data = fdBLUP)
lm7c <- lm(deltaFlow ~ region, data = fdBLUP)
lm7d1 <- lm(deltaFlow ~ pc1, data = fdBLUP)
lm7d2 <- lm(deltaFlow ~ pc1+region, data = fdBLUP)

AIC(lm7a,lm7b,lm7c,lm7d1,lm7d2)

delta7 <- ggplot(fdBLUP, aes(x = pc1, y = deltaFlow,
                             color=region, shape = region))+
  geom_point()+theme_classic()+
  labs(title = "Flowering Day", x = "PC1", 
       y="Flowering Day DR")+
  annotate("text",x=4,y=1.068, label="Model B", 
           size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta7

# # EU only
# fdBLUP_EU<- fdBLUP %>% filter(region == "EU")
#   
# lm7EU <- lm(deltaFlow ~ pc1, data = fdBLUP_EU)
# summary(lm7EU)
# 
# ggplot(fdBLUP_EU, aes(x = pc1, y = deltaFlow, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# fdBLUP_US<- fdBLUP %>% filter(region == "US")
# 
# lm7US <- lm(deltaFlow ~ pc1, data = fdBLUP_US)
# summary(lm7US)
# 
# ggplot(fdBLUP_US, aes(x = pc1, y = deltaFlow, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaFlow ~ region, data= fdBLUP))
```

#### Days from Bolt to Flower

```{r blupbtf, message=F,warning=F}

## stemWidth
tempb<-lmer(boltToFlower ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = df2, REML = F)

btfBLUP<-coef(tempb)$pop

btfBLUP$pop <- rownames(btfBLUP)

names(btfBLUP)[names(btfBLUP) == "(Intercept)"] <- "trtC"
names(btfBLUP)[names(btfBLUP) == "trtT"] <- "trtTcoef"

btfBLUP <- btfBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltabtf = trtT/trtC)

# join with pop data
btfBLUP <- right_join(btfBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(btfBLUP$deltabtf)$out

# remove SP4 
btfBLUP <- btfBLUP %>% 
  filter(pop != "SP4")

# linear model with all
lm8a <- lm(deltabtf ~ pc1*region, data = btfBLUP)
lm8b <- lm(deltabtf ~ 1, data = btfBLUP)
lm8c <- lm(deltabtf ~ region, data = btfBLUP)
lm8d1 <- lm(deltabtf ~ pc1, data = btfBLUP)
lm8d2 <- lm(deltabtf ~ pc1+region, data = btfBLUP)

AIC(lm8a,lm8b,lm8c,lm8d1,lm8d2)

delta8 <- ggplot(btfBLUP, aes(x = pc1, y = deltabtf,
                              color=region, 
                              shape = region))+
  geom_point()+theme_classic()+
  labs(title = "Days from Bolting to Flowering", 
       x = "PC1", 
       y="Bolting to Flowering DR")+
  annotate("text",x=4,y=0.995,label="Model D1",
           size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta8

# # EU only
# btfBLUP_EU<- btfBLUP %>% filter(region == "EU")
#   
# lm8EU <- lm(deltabtf~ pc1, data = btfBLUP_EU)
# summary(lm8EU)
# 
# ggplot(btfBLUP_EU, aes(x = pc1, y = deltabtf, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# btfBLUP_US<- btfBLUP %>% filter(region == "US")
# 
# lm8US <- lm(deltabtf ~ pc1, data = btfBLUP_US)
# summary(lm8US)
# 
# ggplot(btfBLUP_US, aes(x = pc1, y = deltabtf, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltabtf ~ region, data= btfBLUP))
```

#### GSW

```{r blupgsw, message=F,warning=F}

# get info for the pops first
meta <- licor %>%
  group_by(pop) %>%
  summarize(region = first(region),
            state = first(state),
            pc1 = first(pc1),
            lat = first(lat),
            lon = first(lon),
            n = n())

## stemWidth
tempb<-lmer(gswsqrt ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = licor, REML = F)

gswBLUP<-coef(tempb)$pop

gswBLUP$pop <- rownames(gswBLUP)

names(gswBLUP)[names(gswBLUP) == "(Intercept)"] <- "trtC"
names(gswBLUP)[names(gswBLUP) == "trtT"] <- "trtTcoef"

gswBLUP <- gswBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltagsw = trtT/trtC)

# join with pop data
gswBLUP <- right_join(gswBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(gswBLUP$deltagsw)$out

# # remove SP4 
# fdBLUP <- fdBLUP %>% 
#   filter(pop != "SP4")

# linear model with all
lm9a <- lm(deltagsw ~ pc1*region, data = gswBLUP)
lm9b <- lm(deltagsw ~ 1, data = gswBLUP)
lm9c <- lm(deltagsw ~ region, data = gswBLUP)
lm9d1 <- lm(deltagsw ~ pc1, data = gswBLUP)
lm9d2 <- lm(deltagsw ~ pc1+region, data = gswBLUP)

AIC(lm9a,lm9b,lm9c,lm9d1,lm9d2)

delta9 <- ggplot(gswBLUP, aes(x = pc1, y = deltagsw,
                              color=region, shape =region,
                              group = 1))+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",
              # size = 0.5)+
  geom_point()+theme_classic()+
  labs(title = "Stomatal Conductance", 
       x = "PC1", y="Drought Ratio")+
  # annotate("text",x=4,y=0.89,label="Model C",
  #          size = 4,
  #          fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta9

# # EU only
# gswBLUP_EU<- gswBLUP %>% filter(region == "EU")
#   
# lm9EU <- lm(deltagsw~ pc1, data = gswBLUP_EU)
# summary(lm9EU)
# 
# ggplot(gswBLUP_EU, aes(x = pc1, y = deltagsw, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# gswBLUP_US<- gswBLUP %>% filter(region == "US")
# 
# lm9US <- lm(deltagsw ~ pc1, data = gswBLUP_US)
# summary(lm9US)
# 
# ggplot(gswBLUP_US, aes(x = pc1, y = deltagsw, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltagsw ~ region, data= gswBLUP))
```

#### PhiPS2

```{r blupphi, message=F,warning=F}

## stemWidth
tempb<-lmer(PhiPS2 ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = licor, REML = F)

phiBLUP<-coef(tempb)$pop

phiBLUP$pop <- rownames(phiBLUP)

names(phiBLUP)[names(phiBLUP) == "(Intercept)"] <- "trtC"
names(phiBLUP)[names(phiBLUP) == "trtT"] <- "trtTcoef"

phiBLUP <- phiBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaphi = trtT/trtC)

# join with pop data
phiBLUP <- right_join(phiBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(phiBLUP$deltaphi)$out

# # remove SE
phiBLUP <- phiBLUP %>%
  filter(pop != "SE")

# linear model with all
lm10a <- lm(deltaphi ~ pc1*region, data = phiBLUP)
lm10b <- lm(deltaphi ~ 1, data = phiBLUP)
lm10c <- lm(deltaphi ~ region, data = phiBLUP)
lm10d1 <- lm(deltaphi ~ pc1, data = phiBLUP)
lm10d2 <- lm(deltaphi ~ pc1+region, data = phiBLUP)

AIC(lm10a,lm10b,lm10c,lm10d1,lm10d2)

delta10 <- ggplot(phiBLUP, aes(x = pc1, 
                               y = deltaphi, 
                               color=region, 
                               shape = region))+
  geom_point()+theme_classic()+
  labs(title = "Fluorescence", x = "PC1", 
       y="Drought Ratio")+
  # annotate("text",x=3.8,y=0.94,label="Model B",
  #          size = 4,
           # fontface="italic")+ 
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta10

# # EU only
# phiBLUP_EU<- phiBLUP %>% filter(region == "EU")
#   
# lm10EU <- lm(deltaphi~ pc1, data = phiBLUP_EU)
# summary(lm10EU)
# 
# ggplot(phiBLUP_EU, aes(x = pc1, y = deltaphi, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# phiBLUP_US<- phiBLUP %>% filter(region == "US")
# 
# lm10US <- lm(deltaphi ~ pc1, data = phiBLUP_US)
# summary(lm10US)
# 
# ggplot(phiBLUP_US, aes(x = pc1, y = deltaphi, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaphi ~ region, data= phiBLUP))
```

#### TLeaf

```{r bluptl, message=F,warning=F}

## stemWidth
tempb<-lmer(TleafLog ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data = licor, REML = F)

tlBLUP<-coef(tempb)$pop

tlBLUP$pop <- rownames(tlBLUP)

names(tlBLUP)[names(tlBLUP) == "(Intercept)"] <- "trtC"
names(tlBLUP)[names(tlBLUP) == "trtT"] <- "trtTcoef"

tlBLUP <- tlBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltatl = trtT/trtC)

# join with pop data
tlBLUP <- right_join(tlBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(tlBLUP$deltatl)$out

# linear model with all
lm11a <- lm(deltatl ~ pc1*region, data = tlBLUP)
lm11b <- lm(deltatl ~ 1, data = tlBLUP)
lm11c <- lm(deltatl ~ region, data = tlBLUP)
lm11d1 <- lm(deltatl ~ pc1, data = tlBLUP)
lm11d2 <- lm(deltatl ~ pc1+region, data = tlBLUP)

AIC(lm11a,lm11b,lm11c,lm11d1,lm11d2)

delta11 <- ggplot(tlBLUP, aes(x = pc1, y = deltatl,
                              color=region, 
                              shape = region,
                              group = 1))+
  theme_classic()+
  labs(title = "Leaf Temperature (Tleaf)", 
       x = "PC1", y="Tleaf DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",size = 0.5)+
  geom_point()+
  annotate("text",x=4,y=0.997, label="Model C", size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta11

# # EU only
# tlBLUP_EU<- tlBLUP %>% filter(region == "EU")
#   
# lm11EU <- lm(deltatl~ pc1, data = tlBLUP_EU)
# summary(lm11EU)
# 
# ggplot(tlBLUP_EU, aes(x = pc1, y = deltatl, color=region))+
#   geom_point()+theme_classic()
# 
# # US only
# tlBLUP_US<- tlBLUP %>% filter(region == "US")
# 
# lm11US <- lm(deltatl ~ pc1, data = tlBLUP_US)
# summary(lm11US)
# 
# ggplot(tlBLUP_US, aes(x = pc1, y = deltatl, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltatl ~ region, data= tlBLUP))
```

#### est Seed
```{r esdelta,message=F,warning=F}
tempb<-glmmTMB(estSeed ~ trt + (trt|pop) + (1|mat) + (1|blockNum), 
                   data = df2,
                   family = nbinom2,
                   ziformula = ~trt)

esBLUP<-coef(tempb)$cond$pop

esBLUP$pop <- rownames(esBLUP)

names(esBLUP)[names(esBLUP) == "(Intercept)"] <- "trtC"
names(esBLUP)[names(esBLUP) == "trtT"] <- "trtTcoef"

esBLUP <- esBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaEstSeed = trtT/trtC)

# join with pop data
esBLUP <- right_join(esBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(esBLUP$deltaEstSeed)$out

# remove SP5 and WC
esBLUP <- esBLUP %>% 
  filter(pop != "SP5")%>% 
  filter(pop != "WC")

# linear model with all
lm12a <- lm(deltaEstSeed ~ pc1*region, data = esBLUP)
lm12b <- lm(deltaEstSeed ~ 1, data = esBLUP)
lm12c <- lm(deltaEstSeed ~ region, data = esBLUP)
lm12d1 <- lm(deltaEstSeed ~ pc1, data = esBLUP)
lm12d2 <- lm(deltaEstSeed ~ pc1+region, data = esBLUP)

AIC(lm12a,lm12b,lm12c,lm12d1,lm12d2)

delta12 <- ggplot(esBLUP, aes(x = pc1, y = deltaEstSeed, color=region))+
  geom_point()+theme_classic()+
  labs(y = "Estimated Seed DR",
       x="PC1", title = "Estimated Seed")+
  annotate("text",x=4,y=0.99, label="Model D2", size = 4,
           fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta12

# # EU only
# esBLUP_EU<- esBLUP %>% filter(region == "EU")
#   
# lm12EU <- lm(deltaEstSeed ~ pc1, data = esBLUP_EU)
# summary(lm12EU)
# 
# delta12EU <- ggplot(esBLUP_EU, aes(x = pc1, 
#                                  y = deltaEstSeed,
#                                  color=region,
#                                  shape = region))+
#   geom_point()+theme_classic()+
#   labs(title = "Estimated Seed", x = "PC1", 
#        y="delta Estimated Seed")
# 
# delta12EU
# 
# # US only
# esBLUP_US<- esBLUP %>% filter(region == "US")
# 
# lm1US <- lm(deltaEstSeed ~ pc1, data = esBLUP_US)
# summary(lm1US)
# 
# ggplot(esBLUP_US, aes(x = pc1, y = deltaEstSeed, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaEstSeed ~ region, data=esBLUP))
```

#### est Seed with unbolted set to 0

```{r esdelta2,message=F,warning=F}
df4 <- df2 %>% mutate(estSeed = ifelse(is.na(estSeed), 0, estSeed))

tempb<-glmmTMB(estSeed ~ trt + (trt|pop) + (1|mat) + (1|blockNum), 
                   data = df4,
                   family = nbinom2,
                   ziformula = ~trt)

esBLUP2<-coef(tempb)$cond$pop

esBLUP2$pop <- rownames(esBLUP2)

names(esBLUP2)[names(esBLUP2) == "(Intercept)"] <- "trtC"
names(esBLUP2)[names(esBLUP2) == "trtT"] <- "trtTcoef"

esBLUP2 <- esBLUP2 %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaEstSeed = trtT/trtC)

# join with pop data
esBLUP2 <- right_join(esBLUP2,meta,by = "pop")

# check for outliers
# boxplot.stats(esBLUP2$deltaEstSeed)$out
# 
# # remove SP5 and WC
# esBLUP <- esBLUP %>%
#   filter(pop != "SP5")%>%
#   filter(pop != "WC")%>%
#   filter(pop != "SP2")

# linear model with all
# lm122 <- lm(deltaEstSeed ~ pc1+region, data = esBLUP2)
# summary(lm122)

lm122a <- lm(deltaEstSeed ~ pc1*region, data = esBLUP2)
lm122b <- lm(deltaEstSeed ~ 1, data = esBLUP2)
lm122c <- lm(deltaEstSeed ~ region, data = esBLUP2)
lm122d1 <- lm(deltaEstSeed ~ pc1, data = esBLUP2)
lm122d2 <- lm(deltaEstSeed ~ pc1+region, data = esBLUP2)

AIC(lm122a,lm122b,lm122c,lm122d1,lm122d2)

delta122 <- ggplot(esBLUP2, aes(x = pc1, y = deltaEstSeed, color=region, shape = region))+
  geom_point()+theme_classic()+
  labs(y = "Drought Ratio",
       x="PC1", title = "Estimated Seed")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",size = 0.5)+
  geom_point()+
  # annotate("text",x=4,y=0.55, label="Model A", size = 4,
  #          fontface="italic")+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta122

# # EU only
# esBLUP_EU2<- esBLUP2 %>% filter(region == "EU")
#   
# lm12EU2 <- lm(deltaEstSeed ~ pc1, data = esBLUP_EU2)
# summary(lm12EU2)
# 
# delta12EU2 <- ggplot(esBLUP_EU2, aes(x = pc1, 
#                                  y = deltaEstSeed,
#                                  color=region,
#                                  shape = region))+
#   geom_point()+theme_classic()+
#   labs(title = "Estimated Seed", x = "PC1", 
#        y="delta Estimated Seed")
# 
# delta12EU2
# 
# # US only
# esBLUP_US2<- esBLUP2 %>% filter(region == "US")
# 
# lm1US2 <- lm(deltaEstSeed ~ pc1, data = esBLUP_US2)
# summary(lm1US2)
# 
# ggplot(esBLUP_US2, aes(x = pc1, y = deltaEstSeed, color=region))+
#   geom_point()+theme_classic()
# 
# ###### are means different between regions
# 
# summary(aov(deltaEstSeed ~ region, data=esBLUP2))
```


#### Bolting probability

```{r, boltBLUP, message=F,warning=F}
boltBLUP <- glmer(bolted ~ trt + (trt|pop) + (1|mat) + (1|blockNum), data=df2, family = "binomial")

# log odds logit scale
boltBLUP<-coef(boltBLUP)$pop

boltBLUP$pop <- rownames(boltBLUP)

names(boltBLUP)[names(boltBLUP) == "(Intercept)"] <- "trtC"
names(boltBLUP)[names(boltBLUP) == "trtT"] <- "trtTcoef"

boltBLUP <- boltBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaB = trtT/trtC)

# join with pop data
boltBLUP <- right_join(boltBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(boltBLUP$deltaB)$out

# remove SP5 and WC
boltBLUP <- boltBLUP %>%
  filter(pop != "GER4")

# linear model with all
lm13a <- lm(deltaB ~ pc1*region, data = boltBLUP)
lm13b <- lm(deltaB ~ 1, data = boltBLUP)
lm13c <- lm(deltaB ~ region, data = boltBLUP)
lm13d1 <- lm(deltaB ~ pc1, data = boltBLUP)
lm13d2 <- lm(deltaB ~ pc1+region, data = boltBLUP)

AIC(lm13a,lm13b,lm13c,lm13d1,lm13d2)

delta13 <- ggplot(boltBLUP, aes(x = pc1, 
                             y = deltaB,
                             color=region,
                             shape = region,
                             group=1))+theme_classic()+
  labs(title = "Bolting Probability", x = "PC1", 
       y="Bolt Prob. DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",
              # size = 0.5)+
  geom_point()+
  annotate("text",x=4,y=1.9, label="Model B", size = 4,
           fontface="italic")+
  # annotate("text",x=3.5,y=0.968,label="ANOVA: . p = 0.052",
           # size = 3)+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta13
```

#### SPAD BLUP

```{r spadBLUP, message=F,warning=F}

aovBLUP <- lmer(spad ~ week + trt + week*trt + (trt | pop) + (1 | mat) + (1 | blockNum) + (1 | indID), data = spaddf2, REML = F)

aovBLUP<-coef(aovBLUP)$pop

aovBLUP$pop <- rownames(aovBLUP)

names(aovBLUP)[names(aovBLUP) == "(Intercept)"] <- "trtC"
names(aovBLUP)[names(aovBLUP) == "trtT"] <- "trtTcoef"

aovBLUP <- aovBLUP %>%
  dplyr::select(pop,trtC,trtTcoef) %>%
  mutate(trtT = trtC + trtTcoef) %>%
  mutate(deltaSPAD = trtT/trtC)

# join with pop data
aovBLUP <- right_join(aovBLUP,meta,by = "pop")

# check for outliers
boxplot.stats(aovBLUP$deltaSPAD)$out


# linear model with all
lm14a <- lm(deltaSPAD ~ pc1*region, data = aovBLUP)
lm14b <- lm(deltaSPAD ~ 1, data = aovBLUP)
lm14c <- lm(deltaSPAD ~ region, data = aovBLUP)
lm14d1 <- lm(deltaSPAD ~ pc1, data = aovBLUP)
lm14d2 <- lm(deltaSPAD ~ pc1+region, data = aovBLUP)

AIC(lm14a,lm14b,lm14c,lm14d1,lm14d2)

delta14 <- ggplot(aovBLUP, aes(x = pc1, 
                             y = deltaSPAD,
                             color=region,
                             shape = region,
                             group=1))+theme_classic()+
  labs(title = "Absorbance (SPAD)", x = "PC1", 
       y="SPAD DR")+
  # geom_smooth(formula = y~x,na.rm = T, 
              # method = "lm", color = "black",
              # size = 0.5)+
  geom_point()+
  annotate("text",x=4,y=0.97, label="Model B", size = 4,
           fontface="italic")+
  # annotate("text",x=3.5,y=0.968,label="ANOVA: . p = 0.052",
           # size = 3)+
  scale_color_manual(values = c("#e41a1c","#577eb8"))

delta14
```
#### Delta trait figures

```{r figures3, message=F,warning=F}
### add asterisks for significance?

# Growth traits
#a # Dry Weight
#b # Stem Width
#c # Stem Height

ggarrange(delta1, delta2, delta3 + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2,
          common.legend = T)
# Physiology
#e # RWC
#i # gsw
#j # PhiPS2
#k # Tleaf
ggarrange(delta4, delta9, delta10, delta11 + rremove("x.text"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          common.legend = T)

# Reproduction
#d # Total Cap
#f # Bolting Day
#g # Flowering Day
#h # Days from Bolting to Flowering
#i # Estimated Seed
ggarrange(delta6, delta7, delta8, delta5,delta12,delta122 + rremove("x.text"), 
          labels = c("A", "B", "C","D","E"),
          ncol = 2, nrow = 3,
          common.legend = T)

```