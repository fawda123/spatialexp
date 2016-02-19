library(maptools)  ### (Bivand and Lewin-Koh 2015)
library(scales)    ### (Wickham 2015)
library(spdep)

MS_0510=readShapeSpatial("ky_licking.shp")
names(MS_0510)  
names(MS_0510)[3]="COMID"
dim(MS_0510)

load(file = 'rdata/NLCD2011_Region05.RData')
NLCD2011=NLCD2011_Region05

###merge variable of interest
MS_0510Decid=merge(MS_0510,NLCD2011[,c(1,13)])
dim(MS_0510Decid)

names(MS_0510Decid)

TotalArea=sum(MS_0510Decid$AreaSqKM)
TotalArea 

MaxArea=max(MS_0510Decid$AreaSqKM)
MaxArea    

MinArea=min(MS_0510Decid$AreaSqKM)
MinArea  

MS_0510.nb=poly2nb(MS_0510,row.names=MS_0510$COMID,queen=TRUE)
print(MS_0510.nb)  

listw_0510=nb2listw(MS_0510.nb,style="W")

moran.plot(MS_0510Decid$PctDecid2011Cat,listw_0510,labels=FALSE,xlab="Percent deciduous forest NLCD 2011",ylab="Spatially lagged percent deciduous forest NLCD 2011")
moran.i=moran.test(MS_0510Decid$PctDecid2011Cat,listw_0510);moran.i

legend("bottomright",legend=c(paste("Moran's I statistic =",round(moran.i[[3]][1],3)),paste("p-value=",moran.i[[2]])),bty="n")
