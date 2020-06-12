df_report = as.data.frame(report)
df_list = as.data.frame(list)

total =  merge(df_report, df_list, by="ID")

write.csv2(total, file="total.csv")

as.data.frame()

df_report = as.data.frame(read.csv("D:/MergeTest/report.csv",sep="\t",header=TRUE))

df_list = as.data.frame(read.csv("D:/MergeTest/list.csv",sep="\t",header=TRUE))
total =  merge(df_report, df_list, by="ID")
write.csv2(total, file="total.csv")

#Merg with Pme annotation

df_report2 = as.data.frame(read.csv("D:/MergeTest/report.csv",sep="\t",header=TRUE))

df_list2 = as.data.frame(read.csv("D:/MergeTest/listpme.csv",sep=";",header=TRUE))
total =  merge(df_report2, df_list2, by="Pme_ID")
write.csv2(total, file="totalPme.csv")


x = c(100,200,300,150,120,110,600,125,120,110,800)