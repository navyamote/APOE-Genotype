#Author : Navya Mote
#Date   : 5/20/2019
#Title  : WUSTL Challenge - APOE Data 
#---------------------------------------------------------------------------------------------
# Libraries
library(dplyr)
library(tidyr)
#---------------------------------------------------------------------------------------------
# --------------------------------BEGIN OF SECTION 1------------------------------------------
# Part 1
#To read the CSV input data
data <- read.csv('APOE_toy_Data.csv',fileEncoding="UTF-8-BOM")
apoe_type<-read.csv('APOE_Type.csv',fileEncoding="UTF-8-BOM")#Contains the APOE Genotype reference table

# Preprocessing the data 
data_reshape<-reshape(data, idvar = "Barcode", timevar = "SnpAssayName", direction = "wide")

# Building a base table with data required for analysis
data_apoe<-data.frame("DateofApoe" = data_reshape$`DateofApoe.rs429358-C___3084793_20`,
                      "PI_Name" = data_reshape$`PI_Name.rs429358-C___3084793_20`,
                      "Barcode" = data_reshape$Barcode,
                      "UniquePhenoID" = data_reshape$`UniquePhenoID.rs429358-C___3084793_20`,
                      "rs429358-C___3084793_20" = data_reshape$`Call.rs429358-C___3084793_20`,
                      "rs7412-C____904973_10" = data_reshape$`Call.rs7412-C____904973_10`)

# To map the APOE Genotype based on SNP and reference table
data_apoe_type<-merge(data_apoe, apoe_type, by=c('rs429358.C___3084793_20','rs7412.C____904973_10'),
                      all.x=TRUE)

# To map all the records with APOE genotype
data_final<-merge(data_apoe, data_apoe_type, by=c('DateofApoe','PI_Name', 'Barcode', 'UniquePhenoID'))

# To sort the table based on Barcode
data_final<-data_final[order(data_final$Barcode),]

# To build the final table
table_1<-data.frame("DateofApoe"=data_final$DateofApoe, "PI_Name" = data_final$PI_Name,
                    "Barcode" = data_final$Barcode, "UniquePhenoID" = data_final$UniquePhenoID,
                    "APOE_Genotype" = data_final$APOE.Genotype)

write.csv(table_1,"table_1.csv")
# -----------------------------END OF SECTION 1------------------------------------------------
# ---------------------------------------------------------------------------------------------

# --------------------------BEGIN OF SECTION 2-------------------------------------------------
# Part 2
# No data and Unverified for one barcode
df_unid<-group_by(table_1,UniquePhenoID)
df_bar_cnt<-summarise(df_unid, count_barcode = n())
df_tab<-data.frame("PI_Name" = table_1$PI_Name, "UniquePhenoID" = table_1$UniquePhenoID,
                    "Genotype" = table_1$APOE_Genotype)
df_get_cnt<-merge(df_bar_cnt, df_tab, by = c("UniquePhenoID"))
df_get_cnt<-df_get_cnt[order(df_get_cnt$UniquePhenoID),]

df_get_cnt$Status[ is.na(df_get_cnt$Genotype) & df_get_cnt$count_barcode=="1"]<-"NoData"
df_get_cnt$Status[!is.na(df_get_cnt$Genotype) & df_get_cnt$count_barcode=="1"]<-"Unverified"
# To get individuals with status Nodata
# write.csv(df_get_cnt,"NoData.csv")
df_nd_unv1<-filter(df_get_cnt, Status == "NoData")
df_nd_unv2<-filter(df_get_cnt, Status == "Unverified")
df1<-data.frame("UniquePhenoID" = df_nd_unv1$UniquePhenoID, 
                "Status" = df_nd_unv1$Status)
df2<-data.frame("UniquePhenoID" = df_nd_unv2$UniquePhenoID, 
                "Status" = df_nd_unv2$Status)
# write.csv(df1,"df1.csv")
# ------------------------------------------------------------------------------------
# No data for more than one barcode
df_na_typ<-filter(df_get_cnt,is.na(Genotype) & is.na(Status))
df_nd<-group_by(df_na_typ,UniquePhenoID)
df_nd<-summarise(df_nd, count_id = n())
# write.csv(df_nd,"v2.csv")
df_nd$Status[df_nd$count_id > "1"]<-"NoData"
df_nd<-filter(df_nd, Status == "NoData")
df3<-data.frame("UniquePhenoID" = df_nd$UniquePhenoID, 
                "Status" = df_nd$Status)
# -----------------------------------------------
# Individuals that have more than one barcode but only one having the Genotype
df_na<-filter(df_get_cnt,!is.na(Genotype)& is.na(Status))
# write.csv(df_na,"df_na.csv")
df_unver<-group_by(df_na,UniquePhenoID,Genotype)
df_unverct<-summarise(df_unver, count_Genotype = n())
# write.csv(df_unverct,"df_unverct.csv")
df_v<-group_by(df_unverct,UniquePhenoID)
df_v2<-summarise(df_v, count_id = n())
df_x<-merge(df_unverct, df_v2, by=c("UniquePhenoID"))
df_x$Status[df_x$count_Genotype == "1" & df_x$count_id == "1"]<-"Unverified"
df_x$Status[df_x$count_Genotype > "1" & df_x$count_id == "1"]<-"Verified"
df_x$Status[df_x$count_id > "1"]<-"Conflict"

df_unv1<-filter(df_x, Status == "Unverified")
df_unv2<-filter(df_x, Status == "Verified")
df_unv3<-filter(df_x, Status == "Conflict")
# write.csv(df_unv,"df_unv.csv")
d4<-data.frame("UniquePhenoID" = df_unv1$UniquePhenoID, 
               "Status" = df_unv1$Status)
d5<-data.frame("UniquePhenoID" = df_unv2$UniquePhenoID, 
               "Status" = df_unv2$Status)
d6<-data.frame("UniquePhenoID" = df_unv3$UniquePhenoID, 
               "Status" = df_unv3$Status)
# -------------------------------
df_sum<-rbind(df1,df2,df3,d4,d5,d6)
# write.csv(df_sum,"df_sum.csv")
df_all<-merge(df_tab,df_sum,by = c("UniquePhenoID"))
df_all<-df_all[!duplicated(df_all), ]
df_all<-df_all[order(df_all$UniquePhenoID),]
write.csv(df_all,"table_2.csv")
# ------------------------------END OF SECTION 2------------------------------------------------