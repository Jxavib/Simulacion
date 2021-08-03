# ------------------------------------------------------------------
# Ejercicio 5
# ------------------------------------------------------------------
# MC clásico
nsim <- 1000
lambda <- 0.5
set.seed(1)
u <- runif(nsim%/%2)
x <- - log(u) / lambda

mean(x) # Aprox por MC  # valor teor 1/lambda = 2
vx<-var(x)  # medida precisión

xa <- - log(1-u) / lambda
mean(xa) # Aprox por MC  # valor teor 1/lambda = 2
vxa<-var(xa)  # medida precisión
# NOTA: Varianza supoñendo independencia

corr <- cor(x,xa)
xt <- c(x,xa)
mean(xt) # Aprox por MC  # valor teor 1/lambda = 2
vxt <- var(xa)*(1 + corr) # Estimación varianza supoñendo dependencia

# Estimación del porcentaje de reducción en la varianza
100*(var(x)-vxt )/var(x)
