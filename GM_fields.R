# creates p vector field data 
# each row 4 entries: (q1,q2,p1,p2)

setwd(paste(directory,"/Results GM",sep=""))

q1 <- seq(min(qt[,2]*0.99),max(qt[,2]*1.01),by=q0_OGS[r-1]/500) 

q2 <- seq(min(qt[,3]*0.99),max(qt[,3]*1.01),by=q0_OGS[r]/500) 

Eq <- q1%*%t(q2)*0

for (i in 1:dim(Eq)[1]) for (j in 1:dim(Eq)[2]) Eq[i,j] <- varphi(c(1, q1[i], q2[j] ),nt)

#creates points q1xq2 grid where fields are calculated
q12 <- as.matrix(expand.grid(q1,q2))

muE.q12 <- 0*q12[,1]
for (i in 1:dim(q12)[1]) muE.q12[i] <- varphi(c(1,q12[i,]),0)

varphi.field <- cbind(q12,muE.q12)

varphi.field <- rbind(c("x","y","muE"),varphi.field)

length.p <- dim(q12)[1]

p.q12 <- 0*q12
for (i in 1:length.p) p.q12[i,] <- A_i(c(1,q12[i,]),0)

A.field <- cbind(q12,p.q12)

A.field.scaled <- cbind(A.field[,1:2], (A.field[,1:2] + 0.002*A.field[,3:4]))

#pdf(paste(modelname,"Afield.pdf",sep=""),title="Fields",width=5,height=5)

png(paste(modelname,"Afield.png",sep=""),width=5,height=5,res=600, units = "in")
par(mfrow=c(1,1))

#contourplot(Eq,cuts=100, xlab="q_2",ylab="q_3",row.values =q1, column.values =q2, xlim=c(min(q1),max(q1)), 
#            ylim=c(min(q2),max(q2)), labels = F , asp=1 )
#par(new=TRUE)

#plot(q0_OGS[2],q0_OGS[3],col=4, xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)) , pch=19 , asp=1)
#par(new=TRUE)
#plot(qt[,2],qt[,3],col=2, xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)) , pch=19 , asp=1, axes=F)

plot(q0_OGS[2],q0_OGS[3], xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)) , col=4, pch=19,xlab = bquote(q[1]),ylab=bquote(q[2]))
par(new=TRUE)
plot(qt[,2],qt[,3],col=2, xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)) , pch=19 , axes=F,xlab = bquote(q[1]),ylab=bquote(q[2]))
for (i in 1:length.p) arrows(A.field.scaled[i,1],A.field.scaled[i,2],A.field.scaled[i,3],A.field.scaled[i,4],length = 0.05, 
                             xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)))

dev.off()

# Field E ######################################################################

q1E <- seq(min(qt[,2]*0.999),max(qt[,2]*1.001),by=q0_OGS[r-1]/2000) 

q2E <- seq(min(qt[,3]*0.999),max(qt[,3]*1.001),by=q0_OGS[r]/2000) 

q12E <- as.matrix(expand.grid(q1E,q2E))

length.E <- dim(q12E)[1]

E.q12 <- 0*q12E
for (i in 1:length.E) E.q12[i,] <- E_i(c(1,q12E[i,]),0)

E.field <- cbind(q12E,E.q12)

E.field.scaled <- cbind(E.field[,1:2], (E.field[,1:2] + 0.004*E.field[,3:4]))

#pdf(paste(modelname,"Efield.pdf",sep=""),title="Fields",width=5,height=5)

png(paste(modelname,"Efield.png",sep=""),width=5,height=5,res=600, units = "in")
par(mfrow=c(1,1))

plot(q0_OGS[2],q0_OGS[3], xlim=c(min(q1E),max(q1E)), ylim=c(min(q2E),max(q2E)) , col=4, pch=19,xlab = bquote(q[1]),ylab=bquote(q[2]))
par(new=TRUE)
plot(qt[,2],qt[,3],col=2, xlim=c(min(q1E),max(q1E)), ylim=c(min(q2E),max(q2E)) , pch=19, axes=F,xlab = bquote(q[1]),ylab=bquote(q[2]))
for (i in 1:length.E) arrows(E.field.scaled[i,1],E.field.scaled[i,2],E.field.scaled[i,3],E.field.scaled[i,4],length = 0.05, xlim=c(min(q1E),max(q1E)), ylim=c(min(q2E),max(q2E)) )

dev.off()

A.field <- rbind(c("x","y","p1","p2"),A.field.scaled)

# export muE field data
write.table(varphi.field, file = paste(modelname,"_muEfield", ".dat",sep=""),col.names = F, row.names = F, quote=F)

# export p field data
write.table(A.field, file = paste(modelname,"_pfield", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export qt x data
write.table(qt[,-1], file = paste(modelname,"_q", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export qt x mu data
write.table(cbind(qt[,-1],mu_opt), file = paste(modelname,"_qmu", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export qt x mu data
write.table(cbind(qt[,-1],1), file = paste(modelname,"_qmu1", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export qt x mu data
write.table(cbind(qt[,-1],varphi0), file = paste(modelname,"_qmuE", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export qt x mu data
write.table(cbind(qt[,-1],1.12), file = paste(modelname,"_qmu12", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# initial time point
int <- 390

# averages stating at qt[370,]
aqt2 <- qt[int:nt,]
for (t in int:nt) aqt2[t-int + 1,] <- colSums(qt[(int-1):t,])/(t-int + 2)

amu2 <- 0
for (t in int:nt) amu2[t-int +1] <- sum(mu_opt[(int -1):t])/(t-int +2)

# export average qt x mu data
write.table(cbind(aqt2[,-1],amu2), file = paste(modelname,"_aqmu", ".dat",sep=""), col.names = F, row.names = F, quote=F)

# export average qt x mu data
write.table(cbind(aqt2[,-1],1), file = paste(modelname,"_aqmu1", ".dat",sep=""), col.names = F, row.names = F, quote=F)

setwd(directory)