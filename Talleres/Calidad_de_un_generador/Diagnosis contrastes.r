# ------------------------------------------------------------------
# Repetición de contrastes
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Valores iniciales
set.seed(1) # Fijar semilla para reproductibilidad
n <- 50
nsim <- 100
estadistico <- numeric(nsim)
pvalor <- numeric(nsim)


# ------------------------------------------------------------------
# Realizar contrastes
for(isim in 1:nsim) {
  u <- runif(n)
  tmp <- ks.test(u, "punif",  min=0, max=1)
  estadistico[isim] <- tmp$statistic
  pvalor[isim] <- tmp$p.value
}
#Probar a cambiar por u <- rnorm(n) y tmp <- ks.test(u, "pnorm",  mean = mean(u), sd = sd(u))
# ------------------------------------------------------------------
# Proporción de rechazos
cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")


# ------------------------------------------------------------------
# Análisis del estadístico contraste

# Histograma
hist(estadistico, freq=FALSE)

#...

# ------------------------------------------------------------------
# Análisis de los p-valores
# (si todo 'correcto' con distribución uniforme)

# Histograma
hist(pvalor, freq=FALSE)
curve(dunif(x,0,1), add=TRUE)   #abline(h=1)

# Test de Kolmogorov-Smirnov
ks.test(pvalor, "punif",  min=0, max=1)

