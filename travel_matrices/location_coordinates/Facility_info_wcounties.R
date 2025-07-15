#Sorting all locations with their counties

#getting csvs from previous coding
FQHC_withcounties<-read.csv("FQHC_Location_information.csv")
PPHC_withcounties<-read.csv("PPHC_Location_information.csv")

#transposing csv from above
FQHC_withcounties<-t(FQHC_withcounties)
PPHC_withcounties<-t(PPHC_withcounties)

#renaming columns and rows
colnames(PPHC_withcounties)<-c('address','census tract','county','coordinates')
rownames(PPHC_withcounties)<-paste('PPHC', 1:nrow(PPHC_withcounties))
            
colnames(FQHC_withcounties)<-c('address','census tract','county','coordinates')
rownames(FQHC_withcounties)<-paste('FQHC', 1:nrow(FQHC_withcounties))

Facilities_withcounties<-rbind(FQHC_withcounties,PPHC_withcounties)
