## se carga la libreria

library(qtl)

##hyper 
hyper <- sim.geno(supermap, n.draws=16, error.prob=0.001)

supermap <- read.csv("mapa_Ana.csv")

# PESO 
qtl_peso_b<-scanone(supermap,pheno.col = 5)
plot(qtl_peso_b, chr =1:12, main ="Peso del Fruto")
operm<-scanone(supermap, method="hk", pheno.col = 5,n.perm = 100, perm.Xsp = T)
add.threshold(qtl_peso_b, perms = operm, alpha = 0.05)
which.max(qtl_peso_b$lod)
bayesint(qtl_peso_b, 2)
effectplot(hyper, mname1="chr2_157682163", 
           pheno.col = 5, main = "Grafica de Efecto en Peso del Fruto")


# LARGO
qtl_largo_b<-scanone(supermap,pheno.col = 6)
plot(qtl_largo_b, chr =1:12, main ="Largo del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 6,n.perm = 100, perm.Xsp = T)
add.threshold(qtl_largo_b, perms = operm, alpha = 0.05)
which.max(qtl_largo_b$lod)
bayesint(qtl_largo_b, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")


### ANCHO 
qtl_ancho_b<-scanone(supermap,pheno.col = 7)
plot(qtl_ancho_b, chr =1:12, main ="Ancho del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 7,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_ancho_b, perms = operm, alpha = 0.05)
which.max(qtl_ancho_b$lod)
bayesint(qtl_ancho_b, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 7, main = "Grafica de Efecto en Ancho del Fruto")


# FORMA DEL FRUTO 
qtl_tamano_b<-scanone(supermap,pheno.col = 8)
plot(qtl_tamano_b, chr =1:12, main ="Tamaño del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 8,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_tamano_b, perms = operm, alpha = 0.05)
which.max(qtl_tamano_b$lod)
bayesint(qtl_tamano_b, 2)
effectplot(hyper, mname1="chr2_155597316", pheno.col = 8, 
           main = "Grafica de efecto del tamaño del Fruto")
bayesint(qtl_tamano_b, 12)
effectplot(hyper, mname1="chr12_222335343", pheno.col = 8, 
           main = "Grafica de efecto del tamaño del Fruto")


# DEHISCENCIA 
qtl_dehis_b<-scanone(supermap,pheno.col = 9)
plot(qtl_dehis_b, chr =1:12, main ="Dehiscencia del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 9,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_dehis_b, perms = operm, alpha = 0.05)
which.max(qtl_dehis_b$lod)
bayesint(qtl_dehis_b, 10)
effectplot(hyper, mname1="chr10_228669210", pheno.col = 9, 
           main = "Grafica de Efecto en Dehiscencia del Fruto")


# AREA BLUPS
qtl_area_b<-scanone(supermap,pheno.col = 10)
plot(qtl_area_b, chr =1:12, main ="Area del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 10,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_area_b, perms = operm, alpha = 0.05)
which.max(qtl_area_b$lod)
bayesint(qtl_area_b, 2)
effectplot(hyper, mname1="chr2_161091689", pheno.col = 10, 
           main = "Grafica de Efecto en Area del Fruto")


#len_2
qtl_len<-scanone(supermap,pheno.col = 13)
plot(qtl_len, chr =1:12, main ="Largo del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 13,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_len, perms = operm, alpha = 0.05)
which.max(qtl_len$lod)
bayesint(qtl_len, 2)
effectplot(hyper, mname1="chr2_156945939", 
           pheno.col = 13, main = "Grafica de Efecto en Largo del Fruto")


#wid_2
qtl_wid<-scanone(supermap,pheno.col = 14)
plot(qtl_wid, chr =1:12, main ="Ancho del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 14,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_wid, perms = operm, alpha = 0.05)
which.max(qtl_wid$lod)
bayesint(qtl_wid, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")


#area_2
qtl_area<-scanone(supermap,pheno.col = 15)
plot(qtl_area, chr =1:12, main ="Area del fruto")
operm<-scanone(supermap, method="hk", pheno.col = 15,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_area, perms = operm, alpha = 0.05)
which.max(qtl_area$lod)
bayesint(qtl_area, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")


#lw_2
qtl_lw<-scanone(supermap,pheno.col = 16)
plot(qtl_lw, chr =1:12, main ="indice de forma fruto")
operm<-scanone(supermap, method="hk", pheno.col = 16,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_lw, perms = operm, alpha = 0.05)
which.max(qtl_lw$lod)
bayesint(qtl_lw, 12)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")


#perimetro_2
qtl_per<-scanone(supermap,pheno.col = 17)
plot(qtl_per, chr =1:12, main ="Perímetro del Fruto")
operm<-scanone(supermap, method="hk", pheno.col = 17,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_per, perms = operm, alpha = 0.05)
which.max(qtl_per$lod)
bayesint(qtl_per, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en perímetro del Fruto")
qtl_dehis[782,]
qtl_dehis[c(781:783),]

#r_mean_2
qtl_r_mean<-scanone(supermap,pheno.col = 18)
plot(qtl_r_mean, chr =1:12, main ="R_mean")
operm<-scanone(supermap, method="hk", pheno.col = 18,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_r_mean, perms = operm, alpha = 0.05)
which.max(qtl_r_mean$lod)
bayesint(qtl_r_mean, 1)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")
qtl_dehis[782,]
qtl_dehis[c(781:783),]

#b_mean_2
qtl_b_mean<-scanone(supermap,pheno.col = 19)
plot(qtl_b_mean, chr =1:12, main ="B_mean")
operm<-scanone(supermap, method="hk", pheno.col = 19,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_b_mean, perms = operm, alpha = 0.05)
which.max(qtl_b_mean$lod)
bayesint(qtl_b_mean, 3)
effectplot(hyper, mname1="chr3_202213310", 
           pheno.col = 19, main = "Grafica de Efecto en Coloración azul Fruto")
qtl_b_mean[286,]
qtl_dehis[c(781:783),]

#g_mean_2
qtl_g_mean<-scanone(supermap,pheno.col = 20)
plot(qtl_g_mean, chr =1:12, main ="G_mean")
operm<-scanone(supermap, method="hk", pheno.col = 20,n.perm = 1000, perm.Xsp = T)
add.threshold(qtl_g_mean, perms = operm, alpha = 0.05)
which.max(qtl_g_mean$lod)
bayesint(qtl_g_mean, 2)
effectplot(hyper, mname1="chr2_161091689", 
           pheno.col = 6, main = "Grafica de Efecto en Largo del Fruto")
qtl_dehis[782,]
qtl_dehis[c(781:783),]