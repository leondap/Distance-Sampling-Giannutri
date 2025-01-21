##########################################################################
#                                                                        #
#  ANALYSIS OF DISTANCE SAMPLIG FOR THE STUDY                            #
#                                                                        #
#  Estimating wild bee population size with validated distance sampling  #
#                                                                        #
#                                                                        #
#                                                                        #
##########################################################################

library(lubridate)
library(Distance)
library(PMCMRplus)
source("prepare.data.R")
library(mgcv)
library(ggplot2)


# Open the areas of regions as obtained after GIS analysis
ar_prep <- read.table("Area_Env_types.txt",sep="\t",h=T)
area <- aggregate(ar_prep$area~ar_prep$categ.DS,FUN="sum")
colnames(area)<-c("habitat","number")


#Open the observations
data<-read.table("Data_transect_sampling.txt", h=T, sep="\t")

#Obtain average wheather per day
wheather <-unique(data[c("Id.trail","Date","Cloud.coverage","Wind")])
wheather_mean<-cbind(aggregate(wheather$Wind~wheather$Date, FUN="mean"),aggregate(wheather$Cloud.coverage~wheather$Date, FUN="mean")[,2])

#Make a matrix with variables used
usevar <- c("Observer", "Id.trail","habitat","Date","Subject","distance.class")
dat <- data[,usevar]

#Create a id code
dat$id<-paste0(dat$Observer,dat$Id.trail,dat$Date)
volte<-unique(dat$id)

#Create a vector for days
day<-as.POSIXlt(as.character(dat$Date),format="%Y%m%d")
dat$day<-yday(day)
study_days<-unique(dat$day)

day2<-as.POSIXlt(as.character(wheather_mean[,1]),format="%Y%m%d")
wheather_mean$day<-yday(day2)
colnames(wheather_mean)<-c("date","wind","clouds","day")

# Prepare the region table
strata<-cbind(unique(dat$habitat),c(1:length(unique(dat$habitat))))
colnames(strata)<-c("habitat","numero")
strata<-merge(strata,area,by="habitat")
colnames(strata)<-c("habitat","numero","area")
region<-as.data.frame(strata[,2:3])
colnames(region)<-c("Region.Label", "Area")

#Create a table to save the models
megatable<-as.data.frame(matrix(NA,1000,20))

#Select the species to be used
bestie<-c("Anthophora", "Apis", "Bombus")
keys<-c("hn", "hr", "unif") 
riga<-1

for  (gior in 1:length(study_days)){
	day<-study_days[gior]
	for(beast in 1:length(bestie)){
		bestia<-bestie[beast]
		for (k in 1:length(keys)){
			if(k==3){
				adjterms<-c(1:3)
			}else{
				adjterms<-c(2:4)
			}
			keyf<-keys[k]
			datipre<-prepare.data(giorni=day,dat=dat,bestia=bestia,strata=strata)
			if(nrow(datipre$datafinal)>1){
				result <- try({
 				 	model <- ds(datipre$datafinal, 
              			region_table = datipre$region, 
              			sample_table = datipre$sample, 
             		 	key = keyf, 
              			adjustment=NULL, 
              			truncation = 5, 
              			cutpoints = c(0, 1, 2, 3, 4, 5), 
              			obs_table = datipre$obs)
				}, silent = TRUE)

				if (!inherits(result, "try-error")) {
  					a <- summary(model)
  					megatable[riga, 1] <- day
 					megatable[riga, 2] <- bestia
  					megatable[riga, 3] <- keyf
  					megatable[riga, 4] <- "none"
  					megatable[riga, 5] <- a$dht$individuals$N[6, 2]
 					megatable[riga, 6] <- a$dht$individuals$N[6, 3]
 				 	megatable[riga, 7] <- AIC(model)[2]
  					megatable[riga, 8] <- a$dht$individuals$D[1, 2]
  					megatable[riga, 9] <- a$dht$individuals$D[1, 3]
  					megatable[riga, 10] <- a$dht$individuals$D[2, 2]
  					megatable[riga, 11] <- a$dht$individuals$D[2, 3]
  					megatable[riga, 12] <- a$dht$individuals$D[3, 2]
  					megatable[riga, 13] <- a$dht$individuals$D[3, 3]
  					megatable[riga, 14] <- a$dht$individuals$D[4, 2]
  					megatable[riga, 15] <- a$dht$individuals$D[4, 3]
  					megatable[riga, 16] <- a$dht$individuals$D[5, 2]
  					megatable[riga, 17] <- a$dht$individuals$D[5, 3]
  
  					riga <- riga + 1
   					model <- NULL
  					a <- NULL

				} else {
  					cat(sprintf("Model for day %s, bestia %s, and keyf %s failed and was skipped.\n", day, bestia, keyf))
				}
				write.table(megatable,"megatable_red.txt")
				print(paste0(day,bestia,keyf,"none"))
	
				models <- list()  
				for (i in 1:3) {
  					order <- adjterms[i]
 					result <- try({
   						ds(datipre$datafinal, 
      					region_table = datipre$region, 
       					sample_table = datipre$sample, 
      					key = keyf, 
       					truncation = 5, 
       					cutpoints = c(0, 1, 2, 3, 4, 5), 
       					obs_table = datipre$obs, 
      					order = order)
  						
					}, silent = TRUE)
  
  					if (!inherits(result, "try-error")) {
   						models[[i]] <- result
  					} else {
    						cat(sprintf("Model %d failed and was skipped.\n", i))
  					}
				}

				models <- models[!sapply(models, is.null)]

				if (length(models) > 0) {
  					aic_values <- sapply(models, function(model) model$ddf$criterion)
  					best_model_index <- which.min(aic_values)
 					a <- summary(models[[best_model_index]])
					megatable[riga,1]<-day
					megatable[riga,2]<-bestia
					megatable[riga,3]<-keyf
					megatable[riga,4]<-"adjust"
					megatable[riga,5]<-a$dht$individuals$N[6,2]
					megatable[riga,6]<-a$dht$individuals$N[6,3]
					megatable[riga,7]<-AIC(models[[best_model_index]])[2]
					megatable[riga,8]<-a$dht$individuals$D[1,2]
					megatable[riga,9]<-a$dht$individuals$D[1,3]
					megatable[riga,10]<-a$dht$individuals$D[2,2]
					megatable[riga,11]<-a$dht$individuals$D[2,3]
					megatable[riga,12]<-a$dht$individuals$D[3,2]
					megatable[riga,13]<-a$dht$individuals$D[3,3]
					megatable[riga,14]<-a$dht$individuals$D[4,2]
					megatable[riga,15]<-a$dht$individuals$D[4,3]
					megatable[riga,16]<-a$dht$individuals$D[5,2]
					megatable[riga,17]<-a$dht$individuals$D[5,3]				

					model<-NULL
					a<-NULL
					riga<-riga+1
					write.table(megatable,"megatable_red.txt")

				} else {
  					cat("No valid models were generated.\n")
				}
				print(paste0(day,bestia,keyf,"adjust"))

				if(!(keyf=="unif")){
					result <- try({
  						model <- ds(datipre$datafinal, 
             				 region_table = datipre$region, 
             				 sample_table = datipre$sample, 
              				key = keyf, 
             				 formula = ~ observers, 
              				truncation = 5, 
              				cutpoints = c(0, 1, 2, 3, 4, 5), 
              				obs_table = datipre$obs)
					}, silent = TRUE)

					if (!inherits(result, "try-error")) {
  						a <- summary(model)
  						megatable[riga, 1] <- day
  						megatable[riga, 2] <- bestia
  						megatable[riga, 3] <- keyf
  						megatable[riga, 4] <- "observer"
  						megatable[riga, 5] <- a$dht$individuals$N[6, 2]
  						megatable[riga, 6] <- a$dht$individuals$N[6, 3]
  						megatable[riga, 7] <- AIC(model)[2]
  						megatable[riga, 8] <- a$dht$individuals$D[1, 2]
 						megatable[riga, 9] <- a$dht$individuals$D[1, 3]
  						megatable[riga, 10] <- a$dht$individuals$D[2, 2]
  						megatable[riga, 11] <- a$dht$individuals$D[2, 3]
  						megatable[riga, 12] <- a$dht$individuals$D[3, 2]
  						megatable[riga, 13] <- a$dht$individuals$D[3, 3]
  						megatable[riga, 14] <- a$dht$individuals$D[4, 2]
  						megatable[riga, 15] <- a$dht$individuals$D[4, 3]
  						megatable[riga, 16] <- a$dht$individuals$D[5, 2]
  						megatable[riga, 17] <- a$dht$individuals$D[5, 3]
  
  						riga <- riga + 1
  						model <- NULL
  						a <- NULL
					} else {
  						cat(sprintf("Model for day %s, bestia %s, and keyf %s failed and was skipped.\n", day, bestia, keyf))
					}
					print(paste0(day, bestia, keyf, "observer"))
					write.table(megatable,"megatable_red.txt")
				}
			}
		}
	}
}


#After cleaning the models to be excluded by hand the table can be saved and re-opened


megatable<-read.table("megatable_red_2.txt",h=T,sep="\t")
colnames(megatable)<-c("day","species","keyf","random","est","SE","AIC", "meadow","mSE","wood","wSE","h_shrub","hsSE","urban","uSE","l_shrub","lsSE")

#Search for models withn 2DAIC and average them
newtable<-as.data.frame(matrix(NA,1000,15))
riga<-1
for  (gior in 1:length(study_days)){
	for(beast in 1:3){
		tabe<-megatable[which(megatable$day==study_days[gior] & megatable$species==bestie[beast]),]
		if(nrow(tabe>0)){	
			delta_aic <- tabe$AIC - min(tabe$AIC)
			akaike_weights <- exp(-0.5 * delta_aic[delta_aic < 2])
			akaike_weights <- akaike_weights / sum(akaike_weights)
			models<-tabe[which(delta_aic < 2),]
				if(gior+beast==2){
					modellini<-models[,1:4]
				}else{
					modellini<-rbind(modellini,models[,1:4])
				}
			avg_estimate <- sum(models$est* akaike_weights)
			variance <- sum(akaike_weights * (models$SE^2 + (models$est- avg_estimate)^2))
			avg_se <- sqrt(variance)
			avg_estimate_m <- sum(models$meadow* akaike_weights)
			variance <- sum(akaike_weights * (models$mSE^2 + (models$meadow- avg_estimate_m)^2))
			avg_se_m <- sqrt(variance)
			avg_estimate_hs <- sum(models$h_shrub* akaike_weights)
			variance <- sum(akaike_weights * (models$hsSE^2 + (models$h_shrub- avg_estimate_hs)^2))
			avg_se_hs <- sqrt(variance)
			avg_estimate_w <- sum(models$wood* akaike_weights)
			variance <- sum(akaike_weights * (models$wSE^2 + (models$wood- avg_estimate_w)^2))
			avg_se_w <- sqrt(variance)
			avg_estimate_u <- sum(models$urban* akaike_weights)
			variance <- sum(akaike_weights * (models$uSE^2 + (models$urban- avg_estimate_u)^2))
			avg_se_u <- sqrt(variance)
			avg_estimate_ls <- sum(models$l_shrub* akaike_weights)
			variance <- sum(akaike_weights * (models$lsSE^2 + (models$l_shrub- avg_estimate_ls)^2))
			avg_se_ls <- sqrt(variance)
			newtable[riga,1]<-study_days[gior]
			newtable[riga,2]<-bestie[beast]
			newtable[riga,3]<-avg_estimate
			newtable[riga,4]<-avg_se
			newtable[riga,5]<-length(which(delta_aic < 2))
			newtable[riga,6]<-avg_estimate_m
			newtable[riga,7]<-avg_se_m
			newtable[riga,8]<-avg_estimate_hs
			newtable[riga,9]<-avg_se_hs
			newtable[riga,10]<-avg_estimate_w
			newtable[riga,11]<-avg_se_w
			newtable[riga,12]<-avg_estimate_u
			newtable[riga,13]<-avg_se_u
			newtable[riga,14]<-avg_estimate_ls
			newtable[riga,15]<-avg_se_ls
			riga<-riga+1
		}
	}
}

colnames(newtable)<-c("day","species","est","se","models", "meadow","mSE","h_shrub","hsSE","wood","wSE","urban","uSE","l_shrub","lsSE")
newtable<-newtable[complete.cases(newtable),]


# Check which models entered the DAIC
modelliant<-modellini[which(modellini$species=="Apis"),]
aggregate(rep(1,nrow(modelliant))~modelliant$keyf, FUN="sum")
aggregate(rep(1,nrow(modelliant))~modelliant$random, FUN="sum")

modelliant<-modellini[which(modellini$species=="Bombus"),]
aggregate(rep(1,nrow(modelliant))~modelliant$keyf, FUN="sum")
aggregate(rep(1,nrow(modelliant))~modelliant$random, FUN="sum")

modelliant<-modellini[which(modellini$species=="Anthophora"),]
aggregate(rep(1,nrow(modelliant))~modelliant$keyf, FUN="sum")
aggregate(rep(1,nrow(modelliant))~modelliant$random, FUN="sum")

experimental_cond<-read.table("experimental_condition.txt", h=T, sep="\t")

newtable2<-merge(newtable, wheather_mean, by="day")
newtable2<-merge(newtable2,experimental_cond, by="day")


#Make the boxplots for environmental type
newtable3<-newtable2[-(which(newtable2$species=="Apis" & (newtable2$Hives=="Half" | newtable2$Hives=="Closed"))),]
boxplot(newtable3$meadow~newtable3$species, ylim=c(0,0.15))
boxplot(newtable3$l_shrub~newtable3$species,ylim=c(0,0.15))
boxplot(newtable3$h_shrub~newtable3$species,ylim=c(0,0.15))
boxplot(newtable3$wood~newtable3$species, ylim=c(0,0.15))
boxplot(newtable3$urban~newtable3$species,ylim=c(0,0.15))


#Prepapre for Friedman Test
tab1<-newtable3[,c(1:5,6:7)]
tab1<-cbind(tab1,rep(colnames(newtable)[6],nrow(newtable3)))
colnames(tab1)[6:8]<-c("Density","SEd","Type")
for(c in 4:7){
	col<-c*2
	tab2<-newtable3[,c(1:5,col,col+1)]
	tab2<-cbind(tab2,rep(colnames(newtable3)[col],nrow(newtable3)))
	colnames(tab2)[6:8]<-c("Density","SEd","Type")
	tab1<-rbind(tab1,tab2)
}


# Select the species and make Friedman and post hoc tests
tabdis<-tab1[which(tab1$species=="Bombus"),]
boxplot(tabdis$Density~tabdis$Type, ylim=c(0, 0.15))


friedman_result <- friedman.test(Density ~ Type| day, data = tabdis)
print(friedman_result)

# Post hoc comparisons with Benjamini-Hochberg
posthoc_result <- frdAllPairsExactTest(
  y = tabdis$Density,       
  groups = tabdis$Type,  
  blocks = tabdis$day,      
  p.adjust.method = "BH"    
)
print(posthoc_result)



# Make the GAMMs

#Bombus
Bombus_data<-newtable2[which(newtable2$species=="Bombus"),]
Bombus_data$Hives_numeric <- match(Bombus_data$Hives, c("Closed", "Half", "Open"))


Bombus_data$Weights <- 1 / (Bombus_data$se^2)
Bombus_data$Weights_std <- Bombus_data$Weights * (nrow(Bombus_data) / sum(Bombus_data$Weights))
gamm_Bombus <- gamm(est ~ s(day, k = 3) + clouds+wind, 
                  correlation = corAR1(form = ~ day),  weights=Weights_std, data = Bombus_data)

summary(gamm_Bombus$gam)



#Make the plot
pred_data <- data.frame(
  day = seq(min(Bombus_data$day), max(Bombus_data$day), length.out = 100),  # Giorni per il fit
  wind = mean(Bombus_data$wind), 
clouds = mean(Bombus_data$clouds),
  Weights_std = 1            
)


pred <- predict(gamm_Bombus$gam, newdata = pred_data, se.fit = TRUE)
pred_data$pred <- pred$fit
pred_data$lower <- pred$fit - 1.96 * pred$se.fit  # Limite inferiore (95% CI)
pred_data$upper <- pred$fit + 1.96 * pred$se.fit  # Limite superiore (95% CI)

bombus_graph <- ggplot(Bombus_data, aes(x = day, y = est)) +
  
  geom_errorbar(aes(ymin = est - se, ymax = est + se), 
                width = 0.2, color = "black", alpha = 0.6) +
  
  geom_point(aes(size = 1), shape = 1, fill = "blue", alpha = 0.8) +
  
  geom_ribbon(data = pred_data, aes(x = day, ymin = lower, ymax = upper), 
              fill = "red", alpha = 0.2, inherit.aes = FALSE) +
  
  geom_line(data = pred_data, aes(x = day, y = pred), 
            color = "red", size = 1, inherit.aes = FALSE) +
   scale_size(range = c(2, 6)) +
  
  scale_y_continuous(breaks = seq(0, max(Bombus_data$est + Bombus_data$se), by = 5000)) +
 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Rimuovi la griglia principale
    panel.grid.minor = element_blank()   # Rimuovi la griglia secondaria
  ) +
  
  labs(
    title = "Bombus trends",
    x = "Day",
    y = "DS Estimate",
     )



plot(bombus_graph)
ggsave("Bombus_trend.svg", width = 12, height = 6, dpi = 300)



#Anthophora
Anthophora_data<-newtable2[which(newtable2$species=="Anthophora"),]
Anthophora_data$Hives_numeric <- match(Anthophora_data$Hives, c("Closed", "Half", "Open"))

Anthophora_data$Weights <- 1 / (Anthophora_data$se^2)
Anthophora_data$Weights_std <- Anthophora_data$Weights * (nrow(Anthophora_data) / sum(Anthophora_data$Weights))
gamm_Anthophora <- gamm(est ~ s(day, k = 3) + clouds+wind, correlation = corAR1(form = ~ day),  weights=Weights_std, data = Anthophora_data)
summary(gamm_Anthophora$gam)


pred_data <- data.frame(
  day = seq(min(Anthophora_data$day), max(Anthophora_data$day), length.out = 100),  # Giorni per il fit
  wind = mean(Anthophora_data$wind), 
clouds = mean(Anthophora_data$clouds),
  weights_std = 1            
)

pred <- predict(gamm_Anthophora$gam, newdata = pred_data, se.fit = TRUE)
pred_data$pred <- pred$fit
pred_data$lower <- pred$fit - 1.96 * pred$se.fit  # Limite inferiore (95% CI)
pred_data$upper <- pred$fit + 1.96 * pred$se.fit  # Limite superiore (95% CI)


Anthophora_graph <- ggplot(Anthophora_data, aes(x = day, y = est)) +
    geom_errorbar(aes(ymin = est - se, ymax = est + se), 
                width = 0.2, color = "black", alpha = 0.6) +
    geom_point(aes(size = 1), shape = 1, fill = "white", alpha = 0.8) +
    geom_ribbon(data = pred_data, aes(x = day, ymin = lower, ymax = upper), 
              fill = "red", alpha = 0.2, inherit.aes = FALSE) +
   geom_line(data = pred_data, aes(x = day, y = pred), 
            color = "red", size = 1, inherit.aes = FALSE) +
   scale_size(range = c(2, 6)) +
    scale_y_continuous(breaks = seq(0, max(Anthophora_data$est + Anthophora_data$se), by = 25000)) +
   theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  ) +
    labs(
    title = "Anthophora trend",
    x = "Day",
    y = "DS Estimate Anthophora",
     )

plot(Anthophora_graph)

ggsave("Anthophora_graph.svg", width = 12, height = 6, dpi = 300)




# Apis mellifera

Apis_data<-newtable2[which(newtable2$species=="Apis"),]


plot(Apis_data$day,Apis_data$est, ylim=c(0000,150000), cex=1.5)
for(days in 1:nrow(Apis_data)){
arrows(Apis_data$day[days],Apis_data$est[days],Apis_data$day[days],Apis_data$est[days]+Apis_data$se[days], length=0)
arrows(Apis_data$day[days],Apis_data$est[days],Apis_data$day[days],Apis_data$est[days]-Apis_data$se[days], length=0)

}
points(Apis_data[which(Apis_data$Hives=="Half"),1],Apis_data[which(Apis_data$Hives=="Half"),3],pch=16, cex=1.5)


Apis_data$Ratio<-Apis_data$est/Apis_data$Estimate
Apis_data$RatioSE<-Apis_data$se/Apis_data$Estimate

Apis_data$Weights <- 1/(Apis_data$RatioSE^2)
Apis_data$Weights_std <- Apis_data$Weights * (nrow(Apis_data) / sum(Apis_data$Weights))
gamm_Apis <- gamm(Ratio ~ s(day, k = 5) + clouds+wind, correlation = corAR1(form = ~ day),  weights=Weights_std, data = Apis_data)
summary(gamm_Apis$gam)


pred_data <- data.frame(
  day = seq(min(Apis_data$day), max(Apis_data$day), length.out = 100),  # Giorni per il fit
  wind = mean(Apis_data$wind), 
clouds = mean(Apis_data$clouds),
  weights_std = 1            
)



pred <- predict(gamm_Apis$gam, newdata = pred_data, se.fit = TRUE)
pred_data$pred <- pred$fit
pred_data$lower <- pred$fit - 1.96 * pred$se.fit  # Limite inferiore (95% CI)
pred_data$upper <- pred$fit + 1.96 * pred$se.fit  # Limite superiore (95% CI)


Apis_graph <- ggplot(Apis_data, aes(x = day, y = Ratio)) +
    geom_errorbar(aes(ymin = Ratio - RatioSE, ymax = Ratio+ RatioSE), 
                width = 0.2, color = "black", alpha = 0.6) +
    geom_point(aes(size = 1), shape = 1, fill = "white", alpha = 0.8) +
    geom_ribbon(data = pred_data, aes(x = day, ymin = lower, ymax = upper), 
              fill = "red", alpha = 0.2, inherit.aes = FALSE) +
   geom_line(data = pred_data, aes(x = day, y = pred), 
            color = "red", size = 1, inherit.aes = FALSE) +
    scale_size(range = c(2, 6)) +
   scale_y_continuous(breaks = seq(0, max(Apis_data$Ratio +Apis_data$RatioSE), by = 0.05)) +
    theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  ) +
  # Etichette e titolo
  labs(
    title = "Apis ratio trend",
    x = "Day",
    y = "DS Estimate/Hives Estimates in Apis",
     )

plot(Apis_graph)

ggsave("Apis_graph.svg", width = 12, height = 6, dpi = 300)



#Concordance_Anthophora
Anthophora_data$day2<-c(Anthophora_data$day[-1],NA)
prendi<-which(Anthophora_data$day2-Anthophora_data$day==1)
ratios<-Anthophora_data$est[prendi]/Anthophora_data$est[prendi+1]
mean(ratios)
sd(ratios)

#Concordance_Bombus
Bombus_data$day2<-c(Bombus_data$day[-1],NA)
prendi<-which(Bombus_data$day2-Bombus_data$day==1)
ratios<-Bombus_data$est[prendi]/Bombus_data$est[prendi+1]
mean(ratios)
sd(ratios)

#Concordance_Apis
Apis_data$day2<-c(Apis_data$day[-1],NA)
prendi<-which(Apis_data$day2-Apis_data$day==1)
ratios<-Apis_data$est[prendi]/Apis_data$est[prendi+1]
mean(ratios)
sd(ratios)
