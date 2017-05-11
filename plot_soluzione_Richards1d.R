#PER PLOTTARE I RISULTATI DI RICHARDS 1D risolto con nesed newton method
# vedi articolo di Casulli per il confronto della soluzione

setwd('')

#readVector<-read.table('Psi initial condition.txt', skip=3, sep = "", dec=".")
#space = readVector[1:nrow(readVector),2]
#initialCondition = readVector[1:nrow(readVector),1]

readVector<-read.table('Psi_35.txt', skip=3, sep = "", dec=".")
space = readVector[1:nrow(readVector),2]
psi1 = readVector[1:nrow(readVector),1]

readVector<-read.table('Psi_155.txt', skip=3, sep = "", dec=".")
psi2 = readVector[1:nrow(readVector),1]

readVector<-read.table('Psi_300.txt', skip=3, sep = "", dec=".")
psi3 = readVector[1:nrow(readVector),1]

y_min = min( min(psi1),min(psi2),min(psi3) )
y_max = max( max(psi1),max(psi2),max(psi3) )
#plot(space,initialCondition, type = "p", cex=0.6, xlab="Depth [m]", ylab="Pressure [m]", ylim=c(y_min, y_max), col="yellow")

plot(space,psi1, type = "p", cex=0.6, xlab="Depth [m]", ylab="Pressure [m]", ylim=c(y_min, y_max), col="lightgrey")
points(space,psi2, type = "p", col="darkgrey")
points(space,psi3, type = "p", col="black")
legend(1.25, -1.5, legend=c("Time 35000 s", "Time 155000 s", "Time 300000 s"),
       col=c("lightgrey", "darkgrey", "black"), pch=1, cex=0.93)
box()

