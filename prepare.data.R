prepare.data<-function(giorni,dat,bestia,strata){
res<-NULL

#provide a number to the transects 
transetti<-cbind(unique(dat$Id.trail),c(1:length(unique(dat$Id.trail))))

region<-as.data.frame(strata[,2:3])
colnames(region)<-c("Region.Label", "Area")

#prepapre the sample table
trans<-as.data.frame(unique(dat[,c("Id.trail","habitat")]))
trans<-merge(trans,strata,by="habitat")
trans$num<-c(1:nrow(trans))
colnames(trans)[2]<-"Sample.Label"
samples1<-trans[,c("Sample.Label","numero")]
colnames(samples1)<-c("Sample.Label", "Region.Label")
samples1$Region.Label<-as.numeric(samples1$Region.Label)
samples<-merge(samples1,trans, by="Sample.Label")[,c("num","Region.Label")]
colnames(samples)<-c("Sample.Label", "Region.Label")
samples$Sample.Label<-as.numeric(samples$Sample.Label)

datb <- dat[which(dat$day%in%giorni),]
fatti <- unique(datb$Id.trail)
datb$object<-c(1:nrow(datb))
datc<-datb[which(datb$Subject==bestia),]


prendi <- trans[which(trans$Sample.Label%in% fatti),5]
qua<-unique(datb[,c("Id.trail","day")])
repliche<-aggregate(rep(1,nrow(qua))~qua$Id.trail,FUN="sum")
repliche[,2]<-repliche[,2]*30
colnames(repliche)<-c("Sample.Label","Effort")
repliche_trans<-merge(repliche,trans, by="Sample.Label")
repliche_trans<-repliche_trans[,-1]
colnames(repliche_trans)[5]<-"Sample.Label"
obs<-cbind(datc$object,datc$habitat,datc$Id.trail)
colnames(obs)<-c("object", "Region.Label", "Sample.Label")
obs1<-merge(obs,trans,by="Sample.Label")
obs<-obs1[,c("object","numero","num")]
colnames(obs)<-c("object", "Region.Label", "Sample.Label")
obs$object<-as.numeric(obs$object)
obs$Region.Label<-as.numeric(obs$Region.Label)

region$Region.Label<-as.numeric(region$Region.Label)
region$Area<-as.numeric(region$Area)



datafinal<-as.data.frame(cbind(datc$object,datc$Observer,rep(1,nrow(datc)),datc$distance.class,rep(1,nrow(datc)),rep(1,nrow(datc))))
colnames(datafinal)<-c("object", "observers", "detected", "distance", "sex", "exposure")
datafinal$object<-as.numeric(datafinal$object)
datafinal$detected<-as.numeric(datafinal$detected)
datafinal$sex<-as.numeric(datafinal$sex)
datafinal$exposure<-as.numeric(datafinal$exposure)
datafinal$distance<-as.numeric(datafinal$distance)
sample1 <- samples[which(samples$Sample.Label%in%prendi),]
sample2<-merge(sample1,repliche_trans,by="Sample.Label")
sample3<-sample2[,1:3]
datafinal2<-merge(datafinal,obs, by="object")
datafinal2$Region.Label<-as.factor(datafinal2$Region.Label)

res$datafinal<-datafinal2
res$region<-region
res$sample<-sample3
res$obs<-obs
return(res)
}


