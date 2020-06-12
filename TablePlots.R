#Plot statistics for Assembly table

library(base2grob)
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)

assembly_table = as.data.frame(Quast_stat[,1:length(Quast_stat)])

assembly_table

summary(assembly_table[1,])

boxplot(assembly_table[1,1:12])

hist(assembly_table[1,1:22])

contigs500 = assembly_table[9,]
contigs500

soap_contigs500 = contigs500[,1:22]

p = ggplot(assembly_table_t2)
p



boxplot(contigs500[,1:22], contigs500[,25:50])
boxplot(contigs500[,25:50])
boxplot(c(contigs500[,23:24],contigs500[,51:52]))

plot(contigs500[,1:22])

fivenum(contigs500[,1:22])

as.vector(t(contigs500[1,1:22]))

boxplot(as.vector(t(contigs500[1,1:22])))

tdt <- function(inpdt){
  transposed <- t(inpdt[,-1]);
  colnames(transposed) <- inpdt[[1]];
  transposed <- data.table(transposed, keep.rownames=T);
  setnames(transposed, 1, names(inpdt)[1]);
  return(transposed);
}

tdt(assembly_table)

assembly_table_t = tdt(assembly_table)

assembly_table_t2 = assembly_table_t[,2:18]
rownames(assembly_table_t2) = assembly_table_t[1:52,1]
colnames(assembly_table_t2) = c("contigs (>= 0 bp)",
                               "contigs (>= 1000 bp)",
                               "contigs (>= 5000 bp)",
                               "contigs (>= 10000 bp)",
                               "Total length (>= 0 bp)",
                               "Total length (>= 1000 bp)",
                               "Total length (>= 5000 bp)",
                               "Total length (>= 10000 bp)",
                               "contigs (>= 500bp)",
                               "Largest contig",
                               "Total length",
                               "GC (%)",
                               "N50",
                               "N75",
                               "L50",
                               "L75",
                               "N's per 100 kbp",
                               "Group"
)

assembly_table_t2

contigs500 = assembly_table_t2[,3]
contigs500

plot(contigs500)
boxplot(contigs500)

boxplot(assembly_table_t2[,3])
boxplot(assembly_table_t2[,13])


###Transposed Table with groups

assembly_table = as.data.frame(Quast_transposed_grouped[2:53,2:length(Quast_transposed_grouped)])


assembly_table$Group = as.factor(assembly_table$Group)
#assembly_table = t(assembly_table)
#as.numeric(assembly_table)

head(assembly_table)

rownames(assembly_table) = Quast_transposed_grouped[2:53,1]
colnames(assembly_table) = c(   "contigs (>= 0 bp)",
                                "contigs (>= 1000 bp)",
                                "contigs (>= 5000 bp)",
                                "contigs (>= 10000 bp)",
                                "Total length (>= 0 bp)",
                                "Total length (>= 1000 bp)",
                                "Total length (>= 5000 bp)",
                                "Total length (>= 10000 bp)",
                                "contigs (>= 500bp)",
                                "Largest contig",
                                "Total length",
                                "GC (%)",
                                "N50",
                                "N75",
                                "L50",
                                "L75",
                                "N's per 100 kbp",
                                "Group"
)

assembly_colnames = c(   "contigs (>= 0 bp)",
                           "contigs (>= 1000 bp)",
                           "contigs (>= 5000 bp)",
                           "contigs (>= 10000 bp)",
                           "Total length (>= 0 bp)",
                           "Total length (>= 1000 bp)",
                           "Total length (>= 5000 bp)",
                           "Total length (>= 10000 bp)",
                           "contigs (>= 500bp)",
                           "Largest contig",
                           "Total length",
                           "GC (%)",
                           "N50",
                           "N75",
                           "L50",
                           "L75",
                           "N's per 100 kbp",
                           "Group"
)





##SOAP and Velvet: Contigs >500bp

p <- ggplot(assembly_table[1:50,], aes(x=assembly_table[1:50,]$Group, y=as.numeric(assembly_table[1:50,]$`contigs (>= 500bp)`), color=assembly_table[1:50,]$Group)) + geom_boxplot()

p + coord_flip()

p + scale_x_discrete(limits=c("SOAP", "Trinity", "Velvet-Oases")) + coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2))

contigs500 = assembly_table$`contigs (>= 500bp)`
fivenum(as.numeric(contigs500))

##SOAP and Velvet: Contigs >1000bp

p <- ggplot(assembly_table[1:50,], aes(x=assembly_table[1:50,]$Group, y=as.numeric(assembly_table[1:50,]$`contigs (>= 1000 bp)`), color=assembly_table[1:50,]$Group))+
  geom_boxplot()+
  labs(title="Number of Contigs > 1000bp per Assembler",y="n(contigs)",x="Assembler", color="Assembler")

p = p + scale_x_discrete(limits=c("SOAP", "Trinity", "Velvet-Oases")) + coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2))
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=4)

p

##SOAP and Velvet: Contigs >500bp

p <- ggplot(assembly_table[1:50,], aes(x=assembly_table[1:50,]$Group, y=as.numeric(assembly_table[1:50,]$`contigs (>= 5000 bp)`), color=assembly_table[1:50,]$Group))+
  geom_boxplot()+
  labs(title="Number of Contigs > 500bp per Assembler",y="n(contigs)",x="Assembler", color="Assembler")

p = p + scale_x_discrete(limits=c("SOAP", "Trinity", "Velvet-Oases")) + coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2))
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=4)

p

contigs500 = assembly_table$`contigs (>= 500bp)`
fivenum(as.numeric(contigs500))

##SOAP and Velvet: N50

p <- ggplot(assembly_table[1:50,], aes(x=assembly_table[1:50,]$Group, y=as.numeric(assembly_table[1:50,]$N50), color=assembly_table[1:50,]$Group))+
  geom_boxplot()+
  labs(title="Number of Contigs > 500bp per Assembler",y="n(contigs)",x="Assembler", color="Assembler")

p = p + scale_x_discrete(limits=c("SOAP", "Trinity", "Velvet-Oases")) + coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2))
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=4)

p

##SOAP and Velvet: L50

p <- ggplot(assembly_table[1:50,], aes(x=assembly_table[1:50,]$Group, y=as.numeric(assembly_table[1:50,]$L50), color=assembly_table[1:50,]$Group))+
  geom_boxplot()+
  labs(title="Number of Contigs > 500bp per Assembler",y="n(contigs)",x=" ", color="Assembler")

p = p + scale_x_discrete(limits=c("SOAP", "Velvet-Oases")) + coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2))
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=4)
p <- p +  geom_hline(aes(yintercept=mean(c(as.numeric(assembly_table$L50[23]), as.numeric(assembly_table$L50[24]))), colour="Trinity Mean"))
p

trinity_mean_L50 = mean(c(as.numeric(assembly_table$L50[23]), as.numeric(assembly_table$L50[24])))

# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))



#####Transposed Table with groups and kmer size

assembly_table = as.data.frame(Quast_transposed_grouped2[2:53,2:length(Quast_transposed_grouped2)])
assembly_table = as.data.frame(Quast_transposed_grouped2)
assembly_table[37,20]



#assembly_table = t(assembly_table)
#as.numeric(assembly_table)

head(assembly_table)

rownames(assembly_table) = Quast_transposed_grouped[2:53,1]
colnames(assembly_table) = c(   "contigs (>= 0 bp)",
                                "contigs (>= 1000 bp)",
                                "contigs (>= 5000 bp)",
                                "contigs (>= 10000 bp)",
                                "Total length (>= 0 bp)",
                                "Total length (>= 1000 bp)",
                                "Total length (>= 5000 bp)",
                                "Total length (>= 10000 bp)",
                                "contigs (>= 500bp)",
                                "Largest contig",
                                "Total length",
                                "GC (%)",
                                "N50",
                                "N75",
                                "L50",
                                "L75",
                                "N's per 100 kbp",
                                "Group",
                                "Clustered",
                                "k-mer"
)
assembly_table$Group = as.factor(assembly_table$Group)
assembly_table$Clustered = as.factor(assembly_table$Clustered)

##SOAP and Velvet: Contigs >500bp

bp = ggplot(assembly_table[1:48,], aes(x=assembly_table[1:48,]$Group, y=as.numeric(assembly_table[1:48,]$`contigs (>= 500bp)`), fill=assembly_table[1:48,]$Clustered)) +
  geom_boxplot() + labs(title="Number of Contigs > 500bp",y="n(contigs)",x="Assembler", color="Clustered") + scale_fill_discrete(name="Trinity")+
  geom_jitter(shape=8, position=position_jitter(0.1),aes(colour = assembly_table[1:48,]$Clustered)) +   scale_fill_discrete(name="Clustered")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  geom_hline(aes(yintercept=mean(c(as.numeric(assembly_table$`contigs (>= 500bp)`[49]), as.numeric(assembly_table$`contigs (>= 500bp)`[50]))), colour="Trinity Mean"))
bp = bp + coord_flip()

bp

trinity_mean_L50 = mean(c(as.numeric(assembly_table$L50[23]), as.numeric(assembly_table$L50[24])))

##SOAP and Velvet: Contigs >1000bp

bp = ggplot(assembly_table[1:48,], aes(x=assembly_table[1:48,]$Group, y=as.numeric(assembly_table[1:48,]$`contigs (>= 10000 bp)`), fill=assembly_table[1:48,]$Clustered)) +
  geom_boxplot() + labs(title="Number of Contigs > 10000bp",y="n(contigs)",x="Assembler", color="Clustered") + scale_fill_discrete(name="Trinity")+
  geom_jitter(shape=8, position=position_jitter(0.1),aes(colour = assembly_table[1:48,]$Clustered)) +   scale_fill_discrete(name="Clustered")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  geom_hline(aes(yintercept=mean(c(as.numeric(assembly_table$`contigs (>= 10000 bp)`[49]), as.numeric(assembly_table$`contigs (>= 10000 bp)`[50]))), colour="Trinity Mean"))
bp = bp + coord_flip()

bp

##SOAP and Velvet: Contigs >N50

bp = ggplot(assembly_table[1:48,], aes(x=assembly_table[1:48,]$Group, y=as.numeric(assembly_table[1:48,]$N50), fill=assembly_table[1:48,]$Clustered)) +
  geom_boxplot() + labs(title="N50 length of assemblies",y="N50 length (nt)",x="Assembler", color="Clustered") + scale_fill_discrete(name="Trinity")+
  geom_jitter(shape=8, position=position_jitter(0.1),aes(colour = assembly_table[1:48,]$Clustered)) +   scale_fill_discrete(name="Clustered")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  geom_hline(aes(yintercept=mean(c(as.numeric(assembly_table$N50[49]), as.numeric(assembly_table$N50[50]))), colour="Trinity Mean"))
bp = bp + coord_flip()

bp

##SOAP and Velvet: L50

bp = ggplot(assembly_table[1:48,], aes(x=assembly_table[1:48,]$Group, y=as.numeric(assembly_table[1:48,]$L50), fill=assembly_table[1:48,]$Clustered)) +
  geom_boxplot() + labs(title="L50 count of assemblies",y="n(contigs)",x="Assembler", color="Clustered") + scale_fill_discrete(name="Trinity")+
  geom_jitter(shape=8, position=position_jitter(0.1),aes(colour = assembly_table[1:48,]$Clustered)) +   scale_fill_discrete(name="Clustered")+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  geom_hline(aes(yintercept=mean(c(as.numeric(assembly_table$L50[49]), as.numeric(assembly_table$L50[50]))), colour="Trinity Mean"))
bp = bp + coord_flip()

bp



####Barplots:

dat <- read.table(header=TRUE, text='
     cond group result hline
                  control     A     10     9
                  treatment     A   11.5    12
                  control     B     12     9
                  treatment     B     14    12
                  ')
dat
#>        cond group result hline
#> 1   control     A   10.0     9
#> 2 treatment     A   11.5    12
#> 3   control     B   12.0     9
#> 4 treatment     B   14.0    12

# Define basic bar plot
bar = ggplot(assembly_table[1:22,], aes(x=assembly_table[1:22,]$`k-mer`, y=as.numeric(assembly_table[1:22,]$N50), fill=assembly_table[1:22,]$Clustered)) +
  geom_bar(position=position_dodge(), stat="identity")
bar

bar = ggplot(assembly_table[23:48,], aes(x=assembly_table[23:48,]$`k-mer`, y=as.numeric(assembly_table[23:48,]$N50), fill=assembly_table[23:48,]$Clustered)) +
  geom_bar(position=position_dodge(), stat="identity")
bar

# Add the text 
text(bar, assembly_table[1:22,] , paste("n = ",assembly_table[1:22,]$`k-mer`,sep="") ,cex=1) 

# The error bars get plotted over one another -- there are four but it looks
# like two
bp + geom_errorbar(aes(ymax=hline, ymin=hline), linetype="dashed")



# Data of Velvet N50 L50
name= assembly_table[23:48,]$`k-mer`
average= as.numeric(assembly_table[23:48,]$L50)
number= as.numeric(assembly_table[23:48,]$N50)
data=data.frame(name,average,number)

# Basic Barplot of L50 and N50 vs assembly size
my_bar_Velv=barplot(data$average , border=F , names.arg=data$name , las=3 , col=c(rgb(0.3,0.1,0.4,0.6), rgb(0.3,0.9,0.4,0.6)) , 
    ylim=c(0,12000), xlim=c(0,91), width=3, 
    main="Barplot of L50 count against SOAP k-mer size" )
my_bar_Velv + abline(h=10025, col="darkseagreen", lty = 2, lwd = 3)
my_bar_Velv + abline(h=2144, col="darkseagreen4", lty = 2, lwd = 3)
#Label axis and title



# Add the text 
my_bar_Velv + text(my_bar_Velv, data$average , paste("N50 = ","\n",data$number,sep="") ,cex=0.9) 
#text(my_bar_Velv, data$number , paste("#",sep="") ,cex=0.9) 

my_bar_Velv + text(-1.6, 10025, "L50 \n Trinity25", col = "darkseagreen", cex=0.8) 
my_bar_Velv + text(-1.6, 2144, "N50 \n Trinity25", col = "darkseagreen4", cex=0.8) 

my_bar_Velv + mtext("2144", side = 2, col = "darkseagreen4", cex=1, adj= 0.18) 
#Legende
legend("topright", legend = c("Unclustered","Clustered" ) , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.9,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))


# Data of SOAP N50 L50
name= assembly_table[1:22,]$`k-mer`
average= as.numeric(assembly_table[1:22,]$L50)
number= as.numeric(assembly_table[1:22,]$N50)
data=data.frame(name,average,number)

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_Velv=barplot(data$average , border=F , names.arg=data$name, las=3 , col="white", xaxt='n', yaxt='n',
                    ylim=c(0,12000), xlim=c(0,91), width=3.6, 
                    plot=TRUE
                    
)




abline(h=10025, col="darkseagreen", lty = 2, lwd = 3)
abline(h=2144, col=alpha("darkseagreen4",0.9), lty = 3, lwd = 3)

par(new=TRUE)

#Label axis and title

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_Velv=barplot(data$average , border=F , names.arg=data$name , las=3 , col=c(rgb(0.3,0.1,0.4,0.6), rgb(0.3,0.9,0.4,0.6)) , 
                    ylim=c(0,12000), xlim=c(0,91), width=3.6, 
                    main="Barplot of L50 count plus N50 size \nagainst SOAP k-mer size \n compared to Trinity assembly",
                    xlab="SOAP k-mer size",
                    ylab="L50 count"
)


# Add the text 
my_bar_Velv + text(my_bar_Velv, data$number , paste("N50 = ","\n",data$number,sep="") ,cex=0.9) 
my_bar_Velv + text(my_bar_Velv, data$average , paste("L50 = ","\n",data$average,sep="") ,cex=0.9)
#text(my_bar_Velv, data$number , paste("#",sep="") ,cex=0.9) 

 
my_bar_Velv + text(-1.6, 2144, "N50 \n Trinity25", col = "darkseagreen4", cex=0.8) 
my_bar_Velv + text(-1.6, 10025, "L50 \n Trinity25", col = "darkseagreen", cex=0.8)

my_bar_Velv + mtext("2144nt", side = 2, col = "darkseagreen4", cex=1, adj= 0.18) 
#Legende
legend("topright", legend = c("Unclustered","Clustered" ) , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.9,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.1, 0.05))

########################Velvet#############################
# Data of Velvet N50 L50
name= assembly_table[23:48,]$`k-mer`
average= assembly_table[23:48,]$L50
average = as.numeric(levels(average))[average]
number= assembly_table[23:48,]$N50
number = as.numeric(levels(number))[number]
data=data.frame(name,average,number)

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_Velv=barplot(data$average , border=F , names.arg=data$name, las=3 , col="white", xaxt='n', yaxt='n',
                    ylim=c(0,12000), xlim=c(0,91), width=3, 
                    plot=TRUE
                    
)




abline(h=10025, col=alpha("darkseagreen4",0.7), lty = 2, lwd = 3)
abline(h=2144, col=alpha("darkseagreen4",0.8), lty = 3, lwd = 3)

par(new=TRUE)

#Label axis and title

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_Velv=barplot(data$average , border=F , names.arg=data$name , las=3 , col=c(rgb(0.3,0.1,0.4,0.65), rgb(0.3,0.9,0.4,0.65)) , 
                    ylim=c(0,12000), xlim=c(0,91), width=3, 
                    main="Barplot of L50 count plus N50 size \nagainst Velvet k-mer size \n compared to Trinity assembly",
                    xlab="Velvet k-mer size",
                    ylab="L50 count"
)


# Add the text 
my_bar_Velv + text(my_bar_Velv, data$number , paste("N50 ","\n",data$number,sep="") ,cex=0.8, adj = c(0.55,0.78), font=2) 
my_bar_Velv + text(my_bar_Velv, data$average , paste("L50 ","\n",data$average,sep="") ,cex=0.9, font=2)
#text(my_bar_Velv, data$number , paste("#",sep="") ,cex=0.9) 


my_bar_Velv + text(-1.6, 2144, "N50 \n Trinity25", col = "darkseagreen4", cex=0.8) 
my_bar_Velv + text(-1.6, 10025, "L50 \n Trinity25", col = "darkseagreen4", cex=0.8)

my_bar_Velv + mtext("2144nt", side = 2, col = "darkseagreen4", cex=1, adj= 0.18) 
my_bar_Velv + mtext("10025", side = 2, col = "darkseagreen4", cex=1, adj= 0.91) 
#Legende
legend("topright", legend = c("Unclustered","Clustered" ) , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.9,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.1, 0.05))

###############################################################
###SOAP#############################################################################################
# Data of SOAP N50 L50
Quast_transposed_grouped2 <- read.delim("U:/OneDrive/Master/Tabellen/Quast_transposed_grouped2.csv", header=FALSE)

assembly_table = as.data.frame(Quast_transposed_grouped3[2:53,2:length(Quast_transposed_grouped2)])
#assembly_table = as.data.frame(Quast_transposed_grouped2)
assembly_table[37,20]



#assembly_table = t(assembly_table)
#as.numeric(assembly_table)

head(assembly_table)

rownames(assembly_table) = Quast_transposed_grouped[2:53,1]
colnames(assembly_table) = c(   "contigs (>= 0 bp)",
                                "contigs (>= 1000 bp)",
                                "contigs (>= 5000 bp)",
                                "contigs (>= 10000 bp)",
                                "Total length (>= 0 bp)",
                                "Total length (>= 1000 bp)",
                                "Total length (>= 5000 bp)",
                                "Total length (>= 10000 bp)",
                                "contigs (>= 500bp)",
                                "Largest contig",
                                "Total length",
                                "GC (%)",
                                "N50",
                                "N75",
                                "L50",
                                "L75",
                                "N's per 100 kbp",
                                "Group",
                                "Clustered",
                                "k-mer"
)
assembly_table$Group = as.factor(assembly_table$Group)
assembly_table$Clustered = as.factor(assembly_table$Clustered)

name = assembly_table[1:22,]$`k-mer`
average= assembly_table[1:22,]$L50
average = as.numeric(levels(average))[average]
number = assembly_table[1:22,]$N50
number = as.numeric(levels(number))[number]
data=data.frame(name,average,number)

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_SOAP=barplot(data$average , border=F , names.arg=data$name, las=3 , col="white", xaxt='n', yaxt='n',
                    ylim=c(0,12000), xlim=c(0,91), width=3.5, 
                    plot=TRUE
                    
)




abline(h=10025, col=alpha("darkseagreen4",0.7), lty = 2, lwd = 3)
abline(h=2144, col=alpha("darkseagreen4",0.8), lty = 3, lwd = 3)

par(new=TRUE)

#Label axis and title

# Basic Barplot of L50 and N50 vs assembly size SOAP
my_bar_SOAP=barplot(data$average , border=F , names.arg=data$name , las=3 , col=c(rgb(0.3,0.1,0.4,0.65), rgb(0.3,0.9,0.4,0.65)) , 
                    ylim=c(0,12000), xlim=c(0,91), width=3.5, 
                    main="Barplot of L50 count plus N50 size \nagainst SOAP k-mer size \n compared to Trinity assembly",
                    xlab="SOAP k-mer size",
                    ylab="L50 count"
)


# Add the text 
my_bar_SOAP + text(my_bar_SOAP, data$number , paste("N50 ","\n",data$number,sep="") ,cex=0.8, adj = c(0.55,0.78), font=2) 
my_bar_SOAP + text(my_bar_SOAP, data$average , paste("L50 ","\n",data$average,sep="") ,cex=0.9, font=2)
#text(my_bar_SOAP, data$number , paste("#",sep="") ,cex=0.9) 


my_bar_SOAP + text(-1.6, 2144, "N50 \n Trinity25", col = "darkseagreen4", cex=0.8) 
my_bar_SOAP + text(-1.6, 10025, "L50 \n Trinity25", col = "darkseagreen4", cex=0.8)

my_bar_SOAP + mtext("2144nt", side = 2, col = "darkseagreen4", cex=1, adj= 0.18) 
my_bar_SOAP + mtext("10025", side = 2, col = "darkseagreen4", cex=1, adj= 0.91) 
#Legende
legend("topright", legend = c("Unclustered","Clustered" ) , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.9,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.1, 0.05))






###############################################################
###Stacked Barplot Contig Size SOAP#############################################################################################
# Data of SOAP
name= assembly_table[1:22,]$`k-mer`
cond= assembly_table[1:22,]$Clustered
contigs1000= assembly_table[1:22,]$`contigs (>= 1000 bp)`
contigs1000 = as.numeric(levels(contigs1000))[contigs1000]
contigs500= assembly_table[1:22,]$`contigs (>= 500bp)`
contigs500 = as.numeric(levels(contigs500))[contigs500]
contigs5000= assembly_table[1:22,]$`contigs (>= 5000 bp)`
contigs5000 = as.numeric(levels(contigs5000))[contigs5000]

data=data.frame(name,contigs500,contigs1000,contigs5000,cond)

########GGPLOT and MELT
 
mdata = melt(data, id=c("name","cond"))
colnames(mdata) = c("name","cond","variable","value")

mdata <- with(mdata, mdata[order(name,cond, variable, value, decreasing=TRUE),])
mdata
log_value = log10(mdata$value); log_value[63]=0;log_value[66]=0;log_value


Velv_bargroupedstacked = ggplot() +
  geom_bar(data=mdata, aes(y = value, x = cond, fill = variable), stat="identity", position='stack') +  facet_grid( ~ name) +  
  scale_colour_manual(labels = c("Contigs >= 500bp", "Contigs > 1000bp","Contigs > 5000bp")) + 
  geom_text(aes(y = value, x = cond, label = value),fontface = "bold",data=mdata, size = 3,  position = position_stack(vjust=0.5))+
  theme_bw() + scale_y_log10() + labs(colour="Contigs >=",title="Stacked Barplot of logarithmic contig count \nagainst SOAP k-mer size",x="k-mer size (N = unclustered, Y = Clustered)",y="log10 of count")+ 
  scale_fill_discrete(name = "Contigs", labels=c("Contigs >= 500bp", "Contigs >= 1000bp","Contigs >= 5000bp")) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold")
  )

Velv_bargroupedstacked


###############################################################
###Stacked Barplot Contig Size Velv#############################################################################################
# Data of Velv
name= assembly_table[23:48,]$`k-mer`
cond= assembly_table[23:48,]$Clustered
contigs1000= assembly_table[23:48,]$`contigs (>= 1000 bp)`
contigs1000 = as.numeric(levels(contigs1000))[contigs1000]
contigs500= assembly_table[23:48,]$`contigs (>= 500bp)`
contigs500 = as.numeric(levels(contigs500))[contigs500]
contigs5000= assembly_table[23:48,]$`contigs (>= 5000 bp)`
contigs5000 = as.numeric(levels(contigs5000))[contigs5000]

data=data.frame(name,contigs500,contigs1000,contigs5000,cond)

########GGPLOT and MELT

mdata = melt(data, id=c("name","cond"))
colnames(mdata) = c("name","cond","variable","value")

mdata <- with(mdata, mdata[order(name,cond, variable, value, decreasing=TRUE),])
mdata
log_value = log10(mdata$value); log_value[1]=0;log_value[4]=0;log_value


SOAP_bargroupedstacked = ggplot() +
  geom_bar(data=mdata, aes(y = value, x = cond, fill = variable), stat="identity", position='stack') +  facet_grid( ~ name) +  
  scale_colour_manual(labels = c("Contigs >= 500bp", "Contigs > 1000bp","Contigs > 5000bp")) + 
  geom_text(aes(y = value, x = cond, label = value),fontface = "bold",data=mdata, size = 3,  position = position_stack(vjust=0.5))+
  theme_bw() + scale_y_log10() + labs(colour="Contigs >=",title="Stacked Barplot of logarithmic contig count \nagainst Velvet k-mer size",x="k-mer size (N = unclustered, Y = Clustered)",y="log10 of count")+ 
  scale_fill_discrete(name = "Contigs", labels=c("Contigs >= 500bp", "Contigs >= 1000bp","Contigs >= 5000bp")) + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold")
        )
  
SOAP_bargroupedstacked

grid.arrange(Velv_bargroupedstacked,SOAP_bargroupedstacked)

