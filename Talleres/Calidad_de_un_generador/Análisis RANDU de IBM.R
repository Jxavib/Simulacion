# ==================================================================
# Análisis del generador congruencial RANDU de IBM
# ==================================================================
source("Funciones/RANDC.r")   # Cargar RANDC.r
initRANDC(semilla = 321,a =  5,c =  1,m =  512)
initRANDC(54321)    # Fijar semilla para reproductibilidad
nsim <- 9999
u <- RANDCN(nsim)  # Generar

# ------------------------------------------------------------------
# EJEMPLO problemas generador RANDU de IBM
# ------------------------------------------------------------------
require(rgl)
y <- matrix(u, ncol = 3, byrow = TRUE)
plot3d(y)           # Rotar para ver los hiperplanos


# ------------------------------------------------------------------
# Bondad de ajuste
# ------------------------------------------------------------------

# Histograma
hist(u, freq = FALSE)
curve(dunif(x, 0, 1), add = TRUE)   # abline(h=1)

# Distribución empírica
curve(ecdf(u)(x), type = "s", lwd = 2)
curve(punif(x, 0, 1), add = TRUE)

# Test chi-cuadrado
source("Funciones/Test Chi-cuadrado continua.r")    # Cargar chisq.test.cont
chisq.test.cont(u, distribution="unif", nclasses=100, output= FALSE, nestpar=0, min=0, max=1)

# Test de Kolmogorov-Smirnov
ks.test(u, "punif", 0, 1)


# ------------------------------------------------------------------
# Aleatoriedad
# ------------------------------------------------------------------

# Gráfico secuencial
plot(as.ts(u))    # plot(u, type = 'l')

# Gráfico de dispersion retardado
plot(u[-nsim], u[-1], xlab="u_t", ylab="u_t+1", pch=21, bg="white")

#plot(matrix(u,ncol=2,byrow=T))  # Alternativa
# Gráfico de autocorrelaciones
#acf(u) #correlaciones
#pacf(u) #correlaciones parciales

# Test de rachas
library (tseries)
runs.test(as.factor(u > median(u)))

# Test de Ljung-Box
Box.test(u, lag = 10, type = "Ljung")
Box.test(u, lag = 10*log10(length(u))-1, type = "Ljung")
# ------------------------------------------------------------------
# Repetición de contrastes
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Valores iniciales
#initRANDC(semilla = 321,a =  5,c =  1,m =  512)
initRANDC(54321)  #Fijar semilla para reproductibilidad
# set.seed(54321)
n <- 500
nsim <- 1000
estadistico <- numeric(nsim)
pvalor <- numeric(nsim)


# ------------------------------------------------------------------
# Realizar contrastes
for(isim in 1:nsim) {
  u <- RANDCN(n)  #Generar
  #u <- runif(n)
  tmp <- chisq.test.cont(u, distribution="unif", nclasses=100, output= FALSE, nestpar=0, min=0, max=1)
  estadistico[isim] <- tmp$statistic
  pvalor[isim] <- tmp$p.value
}


# ------------------------------------------------------------------
# Proporción de rechazos
cat("\nProporción de rechazos al 1% =", sum(pvalor < 0.01)/nsim, "\n")
cat("Proporción de rechazos al 5% =", sum(pvalor < 0.05)/nsim, "\n")
cat("Proporción de rechazos al 10% =", sum(pvalor < 0.1)/nsim, "\n")


# ------------------------------------------------------------------
# Análisis del estadístico contraste

# Histograma
hist(estadistico, freq=FALSE)
curve(dchisq(x,99), add=TRUE)

# Test ji-cuadrado
chisq.test.cont(estadistico, distribution="chisq", nclasses=20, nestpar=0, df=99)

# Test de Kolmogorov-Smirnov
ks.test(estadistico, "pchisq", df=99,alternative = )


# ------------------------------------------------------------------
# Análisis de los p-valores

# Histograma
hist(pvalor, freq=FALSE)
curve(dunif(x,0,1), add=TRUE)   #abline(h=1)

# Test ji-cuadrado
chisq.test.cont(pvalor, distribution="unif", nclasses=20, nestpar=0, min=0, max=1)

# Test de Kolmogorov-Smirnov
ks.test(pvalor, "punif",  min=0, max=1)
