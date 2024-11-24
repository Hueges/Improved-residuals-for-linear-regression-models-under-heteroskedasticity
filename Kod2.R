library(nortest)
library(car)
library( tseries ) #biblioteka za Dickey-Fuller test
podaci <- read.csv(file.choose(), header=TRUE)
podaci <- podaci[1:50,]  #zelimo manju bazu

attach(podaci)
#gradimo model zavisnosti duzine kocenja od brzine i kvadrata brzine 
model <- lm(api00 ~ meals + full + enroll + I(full^2) , data=podaci)
summary(model)

n <- length(api00)

rez_scaled <- scale(residuals(model))

#Primetimo da obicni reziduali otprilike prate normalnu rasodelu
qqPlot(rstandard(model), ylab="Residuals")

#homoskedasticnost ne deluje ispunjeno
plot(fitted(model), scale(residuals(model)), ylab="Std residuals", xlab="fitted")
#Potvdjujemo prethodno i na ovom plotu
plot(rstandard(model)^2 ~ fitted(model))

#Resavamo homeskedasticnost pravljenjem novih PCA reziduala
# Pravimo dizajn i Hat matricu
X <- model.matrix( ~ meals + full + enroll + I(full^2) , data=podaci)
XtXi <- solve(t(X) %*% X)
H <- X %*% XtXi %*% t(X)

#Pravimo prvo homoskedasticne PCA reziduale da vidimo kako se ponasaju
#PCA se definisu kao sopstveni vektori (I-H) * residuali (homoskedasticni).
PCA_homosk_rez <- t((eigen(diag(n)-H))$vectors)%*%(residuals(model))
#vidimo da su poslednjih 5 nula, kao sto treba da bude.
PCA_homosk_rez

y <- PCA_homosk_rez * PCA_homosk_rez #vektor kvadriranih PCA reziduala

sigma <- vector()
for (i in 1:n){
  sigma[i]=(sum(y)-y[i])/(n-(qr(X)$rank)-1) 
} # ocena disperzije gresaka

PCA_homosk_st_rez <- vector() #vektor za standardizovane PCA reziduale.

#standardizacija PCA reziduala
for(i in 1:n){
  PCA_homosk_st_rez[i] <- PCA_homosk_rez[i]/sqrt(sigma[i]) 
} 

PCA_homosk_st_rez

#QQplot u odnosu na kvantile Studentove raspodele sa n-p-1 stepeni slobode.
#PCA reziduali i ne prate bas studentovu raspodelu, sugerisu vrlo debeo rep.
qqplot(qt(ppoints(PCA_homosk_st_rez), df = n-(qr(X)$rank)-1), PCA_homosk_st_rez ,
       xlab = "Teoretski kvantili t raspodele", ylim=c(-3,3))
qqline(PCA_homosk_st_rez)

# Interval poverenja Studentove raspodele
confidence_interval <- function(vector, interval) {
  vec_sd <- sd(vector)
  
  n <- length(vector)
  
  vec_mean <- mean(vector)
  
  error <- qt((interval + 1)/2, df = n-(qr(X)$rank)-1) * vec_sd / sqrt(n)
  
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
#0 pripada intervalu, te ne sumnjamo na pretpostavku ocekivanja 0
confidence_interval(PCA_homosk_st_rez, 0.95)


#Homoskedasticni PCA reziduali cvrsto potvrduju da je model heteroskedastican.
plot(c(1:45), PCA_homosk_st_rez[1:45], ylab = "PCA homosk. residuals", xlab = "index")
abline(h=0)

#Pravimo heteroskedasticne PCA reziduale.
#Matrice E_i stavimo u listu E, 
E <- list(E0=diag(n), E1=(n/(n-5))*diag(n),
          E2=diag(1/(1-diag(H)), nrow=n, ncol = n),
          E3=diag(1/((1-diag(H))^2), nrow=n, ncol = n))

KV_rez <- diag((residuals(model))^2) #dijagonalna matrica kvadratnih reziduala.

#Sigma_i kreiramo i stavljamo u listu Omega.
Sigma <- list(Sigma_0=KV_rez, Sigma_1=E$E1%*%KV_rez,
              Sigma_2=E$E2%*%KV_rez, Sigma_3=E$E3%*%KV_rez)

#Qi sopstvene vektore stavljamo u listu Q.
Q <- list(Q0=eigen((diag(n)-H)%*%Sigma$Sigma_0),
          Q1=eigen((diag(n)-H)%*%Sigma$Sigma_1),
          Q2=eigen((diag(n)-H)%*%Sigma$Sigma_2),
          Q3=eigen((diag(n)-H)%*%Sigma$Sigma_3))

#Najzad formiramo heteroskedasticne PCA reziduale.
PCA_heter_rez <- list(R0=t(Q$Q0$vectors)%*%(residuals(model)),
                      R1=t(Q$Q1$vectors)%*%(residuals(model)),
                      R2=t(Q$Q2$vectors)%*%(residuals(model)),
                      R3=t(Q$Q3$vectors)%*%(residuals(model)))

#vrsimo standardizaciju deljenjem odgovarajucim sopstvenim vrednostima.
PCA_heter_st_rez <- list(R0=(PCA_heter_rez$R0)/(sqrt(Q$Q0$values)),
                         R1=(PCA_heter_rez$R1)/(sqrt(Q$Q1$values)),
                         R2=(PCA_heter_rez$R2)/(sqrt(Q$Q2$values)),
                         R3=(PCA_heter_rez$R3)/(sqrt(Q$Q3$values)))

#crtamo grafike heteroskedasticnih PCA reziduala sa razlicitim ocenama kovarijacione matrice
par(mfrow=c(2,2))
plot(c(1:45), PCA_heter_st_rez$R0[1:45], xlab="index", ylab="PCA residuals sa omega0")
plot(c(1:45), PCA_heter_st_rez$R1[1:45], xlab="index", ylab="PCA residuals sa omega1")
plot(c(1:45), PCA_heter_st_rez$R2[1:45], xlab="index", ylab="PCA residuals sa omega2")
plot(c(1:45), PCA_heter_st_rez$R3[1:45], xlab="index", ylab="PCA residuals sa omega3")
#Vidimo da su prva 3 grafika vrlo lepo rasprsena pod pretpostavkom heteroskedasticnosti
#Cetvrti grafik je malo slabije rasprsen tako da treba izabrati neku od ostalih ocena.