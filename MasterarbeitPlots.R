library(gridExtra)
library(grid)
library(plyr)
library(ggplot2)
library(gtable)

tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE),core=list(fg_params=list(hjust=0, x=0.01)))

t = tableGrob(blast_tophit_coverage_tr2aacds[,1:3], rows=NULL, theme=tt)

t <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
t <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))

t2 = tableGrob(Backmapping_results, rows=NULL, theme=tt)

t2 <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
t2 <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))

grid.draw(backmap_table)

##unique blastx hits
hj <- matrix(c(1, 0.5), ncol=2, nrow=nrow(NumberOfUniqueBlastxInClusteredAssemblies), byrow=TRUE)
x <- matrix(c(0.95, 0.5), ncol=2, nrow=nrow(NumberOfUniqueBlastxInClusteredAssemblies), byrow=TRUE)

colnames(NumberOfUniqueBlastxInClusteredAssemblies) = c("Assembly", "Number of unique \n blastx hits")

ttuniq <- ttheme_default(core=list(fg_params=list(hjust = as.vector(hj), 
                                               x = as.vector(x))),
                      colhead=list(fg_params=list(hjust=0.5, x=0.5)))

t_unique <- tableGrob(NumberOfUniqueBlastxInClusteredAssemblies[,1:2], rows=NULL, theme=ttuniq)

t_unique <- gtable_add_grob(t_unique,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                     t = 1, b = nrow(t_unique), l = 1, r = ncol(g))
t_unique <- gtable_add_grob(t_unique,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(t_unique))

grid.draw(t_unique)

##Trinity Stats
g <- (Trinity_Stats[,1:3])


colnames(g) <- c("Statistics based on:","all contigs","only on longes \n isoform per gene")

grid.table(g, rows=NULL)

gg <- tableGrob(g[,1:3], rows=NULL, theme=tt)
gg <- gtable_add_grob(gg,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(gg), l = 1, r = ncol(gg))
gg <- gtable_add_grob(gg,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(gg))
grid.draw(gg)


##Backmapping

g <- (Backmapping_results[,1:5])

colnames(g) <- c("Assembly name","Back mapping (unclustered)"," ", "Back mapping (clustered)", "  ")

grid.table(g, rows=NULL)

gg <- tableGrob(g[,1:5], rows=NULL, theme=tt)
gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                      t = 1, b = nrow(gg), l = 1, r = ncol(gg))
gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = ncol(gg))

gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = 3, b=nrow(gg))

gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, r = 1, b=nrow(gg))

grid.draw(gg)

##Table Assemblies vs Noncode

tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE),core=list(fg_params=list(hjust=1, x=0.98)))

g <- (blastn_1e.20_vsNONCODE_output[,1:2])

colnames(g) <- c("Assembly","Blast hit count \nvs. Noncode")

grid.table(g, rows=NULL)

gg <- tableGrob(g[,1:2], rows=NULL, theme=tt)
gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 2, b = nrow(gg), l = 1, r = ncol(gg))
gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                      t = 1, l = 1, b=1,r = ncol(gg))
gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                      t = 1, l = 1, r = ncol(gg), b=nrow(gg))
grid.draw(gg)



tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE),core=list(fg_params=list(hjust=1, x=0.98)))

g <- (blastn_1e.20_vsNONCODE_output[,1:2])
df = data.frame(g[,1],g[,2])
colnames(df) <- c("Assembly name","Blast hit count \nvs. Noncode")
df

gg <- tableGrob(df[,1:2], rows=NULL, theme=tt)

gg <- gtable_add_grob(gg,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                      t = 1, l = 1, r = ncol(gg), b=nrow(gg))


grid.draw(gg)
grid.draw(tableGrob(df, rows=NULL))




df <- data.frame(g[,1],g[,2])
colnames(df) <- c("Assembly name","Blast hit count \nvs. Noncode")


##Spalten spezifisch formatieren;

hj <- matrix(c(1, 0.5), ncol=2, nrow=nrow(df), byrow=TRUE)
x <- matrix(c(0.95, 0.5), ncol=2, nrow=nrow(df), byrow=TRUE)

tt1 <- ttheme_default(core=list(fg_params=list(hjust = as.vector(hj), 
                                               x = as.vector(x))),
                      colhead=list(fg_params=list(hjust=0.5, x=0.5)))




dfGrob <- tableGrob(df, rows = NULL, theme = tt1)




#Grid außen
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 1, r = ncol(dfGrob), b=nrow(dfGrob))

#Grid nach erster Zeile - quer
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 2, l = 1, r = ncol(dfGrob), b=nrow(dfGrob))

#Grid nach erster Spalte - horizontal
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 1, l = 2, r = ncol(dfGrob), b=nrow(dfGrob))


grid.newpage()
grid.draw(dfGrob)


#### Quast comparisson
df <- data.frame(g[,1],g[,2])
colnames(df) <- c("Assembly name","Blast hit count \nvs. Noncode")

q = Firstrow_quast[,1:10]

quast_table <- tableGrob(q[,1:12], theme=tt, widths=1)

hj <- matrix(c(1, 0.5), ncol=2, nrow=nrow(quast_table), byrow=TRUE)
x <- matrix(c(0.95, 0.5), ncol=2, nrow=nrow(quast_table), byrow=TRUE)

tt1 <- ttheme_default(core=list(fg_params=list(hjust = as.vector(hj), 
                                               x = as.vector(x))),
                      colhead=list(fg_params=list(hjust=0.5, x=0.5)))




#Grid außen
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 1, r = ncol(dfGrob), b=nrow(dfGrob))

#Grid nach erster Zeile - quer
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 2, l = 1, r = ncol(dfGrob), b=nrow(dfGrob))

#Grid nach erster Spalte - horizontal
dfGrob <- gtable_add_grob(dfGrob,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 1, l = 2, r = ncol(dfGrob), b=nrow(dfGrob))

tt3 <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE,hjust=1,x=0.02),
                      core=list(fg_params=list(hjust=0, x=0.02), padding = unit(c(3, 3), "mm"),base_size = 9, core=list(fg_fun=list(width = unit(10, "cm")))))

quast_table <- tableGrob(q[,1:12], theme=tt3)
grid.newpage()
grid.draw(quast_table)


tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold",parse=TRUE,hjust=0,x=0.02)),
                     core=list(bg_params = list(fill= "black", col = "white")))


##QUAST1

tt4 <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1, fontface = "bold", fontsize=10)),
                      rowhead=list(fg_params=list(hjust=1, x=0.98, unit(14, "mm"), fontsize=8, fontface = "bold")),
                      colhead=list(fg_params = list(fontsize=11,fontface="bold", hjust=0,x=0.02)))

q = Firstrow_quast[,1:10]
quast_table <- tableGrob(q[,1:10], theme=tt4)

quast_table$widths <- unit(rep(1/23, 23), "npc")
quast_table$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table <- gtable_add_grob(quast_table,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 2, r =11, b=nrow(quast_table))

#Grid nach erster Zeile - quer
quast_table <- gtable_add_grob(quast_table,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 2, l = 2, r = 11, b=nrow(quast_table))

#Grid nach 2. und 4. Spalte
quast_table <- gtable_add_grob(quast_table,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 4, r = 5, b=nrow(quast_table))

#Grid nach 2. und 4. Spalte
quast_table <- gtable_add_grob(quast_table,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                               t = 1, l = 8, r = 9, b=nrow(quast_table))

#Grid nach achter Zeile - quer
quast_table <- gtable_add_grob(quast_table,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                               t = 10, l = 2, r = 11, b=nrow(quast_table))

grid.draw(quast_table)

##QUAST2

tt4 <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1, fontface = "bold", fontsize=10)),
                      rowhead=list(fg_params=list(hjust=1, x=0.98, unit(14, "mm"), fontsize=8, fontface = "bold")),
                      colhead=list(fg_params = list(fontsize=11,fontface="bold", hjust=0,x=0.02)))

q = `2ndrow_quast`[,1:10]
quast_table2 <- tableGrob(q[,1:10], theme=tt4)

quast_table2 = 
  rownames() = c(   "contigs (>= 0 bp)",
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

quast_table2$widths <- unit(rep(1/23, 23), "npc")
quast_table2$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table2 <- gtable_add_grob(quast_table2,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                               t = 1, l = 2, r =11, b=nrow(quast_table2))

#Grid nach erster Zeile - quer
quast_table2 <- gtable_add_grob(quast_table2,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                               t = 2, l = 2, r = 11, b=nrow(quast_table2))

#Grid nach 2. und 4. Spalte
quast_table2 <- gtable_add_grob(quast_table2,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                               t = 1, l = 4, r = 5, b=nrow(quast_table2))

#Grid nach 2. und 4. Spalte
quast_table2 <- gtable_add_grob(quast_table2,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                               t = 1, l = 8, r = 9, b=nrow(quast_table2))

#Grid nach achter Zeile - quer
quast_table2 <- gtable_add_grob(quast_table2,
                               grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                               t = 10, l = 2, r = 11, b=nrow(quast_table2))

grid.draw(quast_table2)

##Quast 3 Row

q = `3row_quast`[,1:10]
quast_table3 <- tableGrob(q[,1:10], theme=tt4)

quast_table3$widths <- unit(rep(1/23, 23), "npc")
quast_table3$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table3 <- gtable_add_grob(quast_table3,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =11, b=nrow(quast_table3))

#Grid nach erster Zeile - quer
quast_table3 <- gtable_add_grob(quast_table3,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 2, l = 2, r = 11, b=nrow(quast_table3))

#Grid nach 2. und 4. Spalte
quast_table3 <- gtable_add_grob(quast_table3,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 4, r = 5, b=nrow(quast_table3))

#Grid nach 2. und 4. Spalte
quast_table3 <- gtable_add_grob(quast_table3,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 8, r = 9, b=nrow(quast_table3))

#Grid nach achter Zeile - quer
quast_table3 <- gtable_add_grob(quast_table3,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 10, l = 2, r = 11, b=nrow(quast_table3))

grid.draw(quast_table3)

##Quast 4 Row

q = `4row_quast`[,1:10]
quast_table4 <- tableGrob(q[,1:10], theme=tt4)

quast_table4$widths <- unit(rep(1/23, 23), "npc")
quast_table4$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table4 <- gtable_add_grob(quast_table4,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =11, b=nrow(quast_table4))

#Grid nach erster Zeile - quer
quast_table4 <- gtable_add_grob(quast_table4,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 2, l = 2, r = 11, b=nrow(quast_table4))

#Grid nach 2. und 4. Spalte
quast_table4 <- gtable_add_grob(quast_table4,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 4, r = 5, b=nrow(quast_table4))

#Grid nach 2. und 4. Spalte
quast_table4 <- gtable_add_grob(quast_table4,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 8, r = 9, b=nrow(quast_table4))

#Grid nach achter Zeile - quer
quast_table4 <- gtable_add_grob(quast_table4,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 10, l = 2, r = 11, b=nrow(quast_table4))

grid.draw(quast_table4)

##Quast 5 Row

q = `5row_quast`[,1:10]
quast_table5 <- tableGrob(q[,1:10], theme=tt4)

quast_table5$widths <- unit(rep(1/23, 23), "npc")
quast_table5$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table5 <- gtable_add_grob(quast_table5,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =11, b=nrow(quast_table5))

#Grid nach erster Zeile - quer
quast_table5 <- gtable_add_grob(quast_table5,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 2, l = 2, r = 11, b=nrow(quast_table5))

#Grid nach 2. und 4. Spalte
quast_table5 <- gtable_add_grob(quast_table5,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 4, r = 5, b=nrow(quast_table5))

#Grid nach 2. und 4. Spalte
quast_table5 <- gtable_add_grob(quast_table5,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 8, r = 9, b=nrow(quast_table5))

#Grid nach achter Zeile - quer
quast_table5 <- gtable_add_grob(quast_table5,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 10, l = 2, r = 11, b=nrow(quast_table5))

grid.draw(quast_table5)


##Quast 6 Row

q = `6row_quast`[,1:8]
quast_table6 <- tableGrob(q[,1:8], theme=tt4)

quast_table6$widths <- unit(rep(1/23, 23), "npc")
quast_table6$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =9, b=nrow(quast_table6))

#Grid nach erster Zeile - quer
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 2, l = 2, r = 9, b=nrow(quast_table6))

#Grid nach 2. und 4. Spalte
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 4, r = 5, b=nrow(quast_table6))

#Grid nach 2. und 4. Spalte
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 8, r = 9, b=nrow(quast_table6))

#Grid nach achter Zeile - quer
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 10, l = 2, r = 9, b=nrow(quast_table6))

grid.draw(quast_table6)




######CNCI

CNCI = as.data.frame(CNCI_table[2:28,])

tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE),core=list(fg_params=list(hjust=1, x=0.95)))


colnames(CNCI) <- c("n(Seqs)","Noncoding \n CNCI", "%","NON\nCODE","% ")
CNCI

g_CNCI <- tableGrob(CNCI[,1:5], theme=tt)

g_CNCI <- gtable_add_grob(g_CNCI,
                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                      t = 1, l = 2, r = ncol(g_CNCI), b=nrow(g_CNCI))
g_CNCI <- gtable_add_grob(g_CNCI,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 2, r = ncol(g_CNCI),b = 1)
g_CNCI <- gtable_add_grob(g_CNCI,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                          t = 1, l = 3, r = 4,b = nrow(g_CNCI))


grid.draw(g_CNCI)

######CNCI + NONCODE

CNCI = as.data.frame(CNCI_table[2:28,])

tt <- ttheme_default(colhead=list(fg_params = list(fontface="bold"),parse=TRUE),core=list(fg_params=list(hjust=1, x=0.95)))


colnames(CNCI) <- c("n(Seqs)","Noncoding \n CNCI", "%")
CNCI

g_CNCI <- tableGrob(CNCI[,1:5], theme=tt)

g_CNCI <- gtable_add_grob(g_CNCI,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 2, r = ncol(g_CNCI), b=nrow(g_CNCI))
g_CNCI <- gtable_add_grob(g_CNCI,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                          t = 1, l = 2, r = ncol(g_CNCI),b = 1)
g_CNCI <- gtable_add_grob(g_CNCI,
                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                          t = 1, l = 3, r = 3,b = nrow(g_CNCI))


grid.draw(g_CNCI)




########Quast ABAL.trimmed.fasta

input_data = as.data.frame(read.csv("U:/OneDrive/Master/Tabellen/quastABAL.txt",sep="\t",header=TRUE))

q = as.data.frame(input_data[,2])


rownames(q) = input_data[,1]
colnames(q) = c("ABALtrimmed")

quast_table6 <- tableGrob(q, theme=tt4)


quast_table6$widths <- unit(rep(1/9, 9), "npc")
quast_table6$widths[1] <- unit(1/5, "npc")
grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =2, b=nrow(quast_table6))

#Grid nach erster Zeile - quer
quast_table6 <- gtable_add_grob(quast_table6,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =2)



grid.draw(quast_table6)

#######busco geneset table

q1 = BUSCO_geneset_table3withpercentage[,1:7]
q = as.data.frame(q1[1:9,2:7])
rownames(q) = q1[1:9,1]
colnames(q) = c("Complete\nBUSCOs","Single-copy","Duplicate",	"Fragmented",	"Missing",	"Total\nBUSCOs")

busco_gene_table <- tableGrob(q[,1:6], theme=tt4)


busco_gene_table$widths <- unit(rep(1/12, 12), "npc")
busco_gene_table$widths[1] <- unit(1/5, "npc")#

grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
busco_gene_table <- gtable_add_grob(busco_gene_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 7, r =2, b=nrow(busco_gene_table))

#Grid nach erster Zeile - quer
busco_gene_table <- gtable_add_grob(busco_gene_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =7)

busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 1, l = 5, r =6, b=nrow(busco_gene_table))

busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 1, l = 3, r =5, b=nrow(busco_gene_table))


grid.draw(busco_gene_table)



######BUSCO assemblies table

q1 = BUSCO_assembly_table[,1:7]
q = as.data.frame(q1[2:nrow(q1),2:7])
rownames(q) = q1[2:nrow(q1),1]

colnames(q) = c("Complete\nBUSCOs","Single-copy","Duplicate",	"Fragmented",	"Missing",	"Total\nBUSCOs")

busco_gene_table <- tableGrob(q[,1:5], theme=tt4)


busco_gene_table$widths <- unit(c(0.2, 0.08, 0.08, 0.062, 0.08, 0.06), "npc")
busco_gene_table$widths[1] <- unit(1/5, "npc")







grid.newpage()

pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                    t = 1, l = 6, r =2, b=nrow(busco_gene_table))

#Grid nach erster Zeile - quer
busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                    t = 1, l = 2, r =6)

busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 1, l = 5, r =5, b=nrow(busco_gene_table))

busco_gene_table <- gtable_add_grob(busco_gene_table,
                                    grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                    t = 1, l = 3, r =5, b=nrow(busco_gene_table))


grid.draw(busco_gene_table)

#######average median contig length

df_average_median = as.data.frame(read.table("U:/OneDrive/Master/Tabellen/average_mean_contig_length.txt",sep="\t",header=TRUE))
colnames(df_average_median) = c("Assembly \n Name","Min. \n contig \n length","Max. \n contig \n length","Average \n contig \n length","Median \n contig \n length")

contig_table <- tableGrob(df_average_median, theme=tt4)


contig_table$widths <- unit(rep(1/12, 6), "npc")

grid.newpage()

#pushViewport(viewport(just="left",x=0.1,width=1,height = 0.9))

#Grid außen
contig_table <- gtable_add_grob(contig_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                                t = 1, l = 2, r =6, b=nrow(contig_table))

#Grid nach erster Zeile - quer
contig_table <- gtable_add_grob(contig_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                t = 1, l = 2, r = 6, b=1)

#Grid nach erster Zeile - längs
contig_table <- gtable_add_grob(contig_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                t = 1, l = 3, r = 3, b=nrow(contig_table))
#Grid nach erster Zeile - längs
contig_table <- gtable_add_grob(contig_table,
                                grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                t = 1, l = 5, r = 5, b=nrow(contig_table))


grid.draw(contig_table)

##Mean Read length Calc
read_lengths = read.table("U:/OneDrive/Master/Tabellen/read_lenghts.txt",sep="\t")
read_lengths2 = read_lengths[,1]*read_lengths[,2]
read_lengths3 = cbind(read_lengths,read_lengths2)
mean_read_length = sum(read_lengths3[,3])/sum(read_lengths3[,1]); mean_read_length


##Coverage, genome length and breadth for 25kmer assemblies

#average per base coverage * (total number of covered bases divided by reference genome length)
trinity_genome_length = 71432699
trinity_breadth = 86.8001 * (60031360 / 71432699); trinity_breadth
#72.94598

velvet_genome_length = 74392846
velvet_breadth = 102.942 * (66002274 / 74392846); velvet_breadth
#91.33144

soap_genome_length = 46434879
soap_breadth = 79.0475 * (35358152 / 46434879); soap_breadth
#60.19125

##average genome coverage
#The average coverage for a whole genome can be calculated from the length of the original genome (G), the number of reads (N), and the average read length (L)
#as N * L / G
n_reads = (18128695)
mean_read_length = sum(read_lengths3[,3])/sum(read_lengths3[,1]); mean_read_length
trinity_genome_length = 71432699

trinity_genome_coverage = n_reads * mean_read_length / trinity_genome_length; trinity_genome_coverage
velvet_genome_coverage = n_reads * mean_read_length / velvet_genome_length; velvet_genome_coverage
soap_genome_coverage = n_reads * mean_read_length / soap_genome_length; soap_genome_coverage

##Coverage for clustered assemblies:
trinity_clustered_length = 58049616
velvet_clustered_length = 52236751
soap_clustered_length = 44710937


trinity_clustered_coverage = n_reads * mean_read_length / trinity_clustered_length; trinity_clustered_coverage
velvet_clustered_coverage = n_reads * mean_read_length / velvet_clustered_length; velvet_clustered_coverage
soap_clustered_coverage = n_reads * mean_read_length / soap_clustered_length; soap_clustered_coverage

######Plot Venndiagram for Database Hits;

library(readxl)
Trinity25_Trinotate_report_hits_venn <- read_excel("U:/OneDrive/Master/Tabellen/Trinity25_Trinotate_report_hits_venn.xlsx")
trinity_venn = Trinity25_Trinotate_report_hits_venn

#upload library
library(VennDiagram)

#Aim is to plot populations of Trinotate Report (Databases are populations)
trin_blastp = na.omit(as.character(trinity_venn$Blastp))
trin_pfam = na.omit(as.character(trinity_venn$PFAM))
trin_signalp = na.omit(as.character(trinity_venn$SignalP))
trin_tmhmm = na.omit(as.character(trinity_venn$TmHMM))
trin_goblast = na.omit(as.character(trinity_venn$`GO blast`))



#The goal of the Venn Diagram is to count how many words are common between SNP_pop_1 and SNP_pop_2, between SNP_pop_1 and SNP_pop_3 and so on...
#The venn.diagram function do it automatically and draw it! (you will get a png file in your current working directory)

venn.diagram(
  list(trin_blastp, trin_pfam, trin_signalp,trin_tmhmm,trin_goblast),
  category.names = c("Uniprot","Pfam", "signalP","TmHMM","GO Blast"),
  filename = '#14_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  #fill = c('yellow', 'purple', 'green',"magenta", "dark green"),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  #cat.fontface = "bold",
  #cat.default.pos = "outer",
  #cat.pos = c(-27, 27, -54, 54, 135),
  #cat.dist = c(0.055, 0.055, 0.35, 0.35, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

library(VennDiagram)
library(gridExtra)
library(grid)
library(lattice)
library(grid)
temp1 <- venn.diagram(list(trin_blastp, trin_signalp, trin_tmhmm),
                     fill = c('purple', "magenta","dark green"), 
                     alpha = c(0.5, 0.5,0.5), 
                     category.names = c("Uniprot","signalP","TmHMM"),
                     cex = 1,cat.fontface = 4,
                     cat.cex = 0.6, fontface = "bold",lty =2, fontfamily = "sans", filename = NULL, lwd=2, rotation = 1
                    
                     )


temp2 <- venn.diagram(list(trin_blastp, trin_pfam, trin_goblast, trin_kegg),
                     fill = c('purple', "magenta","dark green","red"), 
                     alpha = c(0.5, 0.5,0.5,0.5), 
                     category.names = c("Uniprot","Pfam","GO Blast","Kegg"),
                     cex = 1,cat.fontface = 4,
                     cat.cex = 0.6, fontface = "bold",lty =2, fontfamily = "sans", filename = NULL, lwd=2
                     
)
grid.draw(temp2)
grid.draw(temp1)
grid.arrange(temp1,temp2, ncol=2)

venn1 = venn.diagram(list(trin_blastp, trin_pfam, trin_goblast, trin_tmhmm),
                     fill = c('purple', "magenta","dark green","red"), 
                     alpha = c(0.5, 0.5,0.5,0.5), 
                     category.names = c("1","2","3","4"),
                     filename = NULL,
                     rotation = 1
                     )
grid.draw(venn1)



# start new page
plot.new() 

pdf("testpdf", width = 14, height = 7)
# setup layout
gl <- grid.layout(nrow=1, ncol=2)
# grid.show.layout(gl)

# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 

# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)

# start new base graphics in first viewport
par(new=TRUE, fig=gridFIG())

grid.draw(temp2)

# done with the first viewport
popViewport()

# move to the next viewport
pushViewport(vp.2)

grid.draw(temp1)

# done with this viewport
popViewport(1)

dev.off()



####Venn for AugustusIN
library(readxl)
augIN_venn <- read_excel("U:/OneDrive/Master/Trinotate_reports/Genes/reportxls_withFA/Reports/augIN.proteinID.blast.pfam.egg.kegg.go.signal.tmhmm.xlsx")


#upload library
library(VennDiagram)

#Aim is to plot populations of Trinotate Report (Databases are populations)
augIN_blastp = na.omit(as.character(augIN_venn$Blastp))
augIN_pfam = na.omit(as.character(augIN_venn$PFAM))
augIN_signalp = na.omit(as.character(augIN_venn$SignalP))
augIN_tmhmm = na.omit(as.character(augIN_venn$TmHMM))
augIN_goblast = na.omit(as.character(augIN_venn$`GO Blast`))
augIN_gopfam = na.omit(as.character(augIN_venn$`GO PFAM`))
augIN_egg = na.omit(as.character(augIN_venn$`EGGNOG`))
augIN_kegg = na.omit(as.character(augIN_venn$`KEGG`))
temp3 <- venn.diagram(list(augIN_blastp, augIN_pfam, augIN_goblast, augIN_kegg),
                      fill = c('purple', "magenta","dark green","red"), 
                      alpha = c(0.5, 0.5,0.5,0.5), 
                      category.names = c("Uniprot","Pfam","GO Blast","Kegg"),
                      cex = 1,cat.fontface = 4,
                      cat.cex = 0.6, fontface = "bold",lty =2, fontfamily = "sans", filename = NULL, lwd=2
                      
)
grid.draw(temp3)

temp3 <- venn.diagram(list(augIN_blastp, augIN_pfam, augIN_goblast, augIN_kegg,augIN_egg),
                      fill = c('purple', "magenta","dark green","red","yellow"), 
                      alpha = c(0.5, 0.5,0.5,0.5,0.3), 
                      category.names = c("Uniprot","Pfam","GO Blast","Kegg","Eggnog"),
                      cex = 1,cat.fontface = 4,
                      cat.cex = 1.5, fontface = "bold",lty =1, fontfamily = "serif", filename = NULL, lwd=2
                      
)
grid.draw(temp3)

temp4 <- venn.diagram(list(augIN_blastp, augIN_pfam, augIN_goblast, augIN_kegg),
                      fill = c('purple', "magenta","dark green","red"), 
                      alpha = c(0.5, 0.5,0.5,0.5), 
                      category.names = c("Uniprot","Pfam","GO Blast","Kegg"),
                      cex = 1,cat.fontface = 4,
                      cat.cex = 1, fontface = "bold",lty =1, fontfamily = "serif", filename = NULL, lwd=2
                      
)
grid.draw(temp4)

