# ==================================================================
# Práctica 7: Aplicaciones de la simulación II
# ==================================================================

# ==================================================================
# Integración Monte Carlo
# Error O(n^-(1/2)) (Independiente del nº de dimensiones)

# ------------------------------------------------------------------
# Ejercicio 1
# ------------------------------------------------------------------

# -------------------------------
mc.integral <- function(fun, a, b, n) {
    # Integración Monte Carlo de fnt entre a y b utilizando una muestra de tamaño n
    # 
    # fnt es una función de una sola variable Se asume a < b y n entero positivo
    # -------------------------------
    x <- runif(n, a, b)
    fx <- sapply(x, fun)
    return(mean(fx) * (b - a))
}


fun <- function(x) return(4 * x^3)
set.seed(54321)
mc.integral(fun, 0, 1, 20)
mc.integral(fun, 0, 1, 100)

# No es eficiente:
S <- function(n) mc.integral(fun, 0, 1, n)
n.vec <- seq(10, 2000, by = 10)
S.vec <- sapply(n.vec, S)
plot(n.vec, S.vec, type = "l", xlab = "n", ylab = "integral")
abline(h = 1, lty = 2)


# Alternativa para representar gráficamente
mc.integral <- function(fun, a, b, n, plot = TRUE) {
    fx <- sapply(runif(n, a, b), fun) * (b - a)
    if (plot) {
        estint <- cumsum(fx)/(1:n)
        esterr <- sqrt(cumsum((fx - estint)^2))/(1:n)
        plot(estint, ylab = "Media y rango de error", type = "l", lwd = 2, ylim = mean(fx) + 
            2 * c(-esterr[1], esterr[1]), xlab = "Iteraciones")
        lines(estint + 2 * esterr, col = "darkgray", lwd = 2)
        lines(estint - 2 * esterr, col = "darkgray", lwd = 2)
        return(list(valor = estint[n], error = 2 * esterr[n]))
    } else return(list(valor = mean(fx), error = 2 * sd(fx)/sqrt(n)))
}

set.seed(54321)
mc.integral(fun, 0, 1, 5000)
abline(h = 1, lty = 2)

set.seed(54321)
mc.integral(fun, 0, 1, 5000, plot = FALSE)


# ------------------------------------------------------------------
# Ejercicio 2
# ------------------------------------------------------------------

curve(dnorm(x), 4.5, 6, ylab = "dnorm(x) y dexp(x-4.5)*k")
abline(v = 4.5)
abline(h = 0)
curve(dexp(x - 4.5) * dnorm(4.5), add = TRUE, lty = 2)  # Reescalado para comparación

set.seed(54321)
nsim <- 10^3
y <- rexp(nsim) + 4.5
w <- dnorm(y)/dexp(y - 4.5)
plot(cumsum(w)/1:nsim, type = "l", ylab = "Aproximación", xlab = "Iteraciones")
abline(h = pnorm(-4.5), lty = 2)

pnorm(-4.5)
mean(w)
sqrt(var(w)/nsim)

est <- mean(rnorm(nsim) > 4.5)
est
sqrt(est * (1 - est)/nsim)


# ------------------------------------------------------------------ 
# Ejercicio 3
# ------------------------------------------------------------------

nsim <- 10^5
set.seed(4321)
y <- rnorm(nsim)
w <- dcauchy(y)/dnorm(y)

plot(cumsum(w * (y > 2) * (y < 6))/1:nsim, type = "l", ylab = "Aproximación", xlab = "Iteraciones")
abline(h = pcauchy(6) - pcauchy(2), lty = 2)
# Mala elección de la densidad auxiliar

boxplot(nsim * w/sum(w))  #Gráfico pesos de cada observación (la suma total es nsim)


# ------------------------------------------------------------------ 
# Ejercicio 4
# ------------------------------------------------------------------ 
# Simular normal a partir de Cauchy (Sampling Importance Resampling)

nsim <- 10^3
nsim2 <- 10^5
set.seed(4321)
y <- rcauchy(nsim2)
w <- dnorm(y)/dcauchy(y)
rx <- sample(y, nsim, prob = w/sum(w))
hist(rx, freq = FALSE)
curve(dnorm, add = TRUE)


# ------------------------------------------------------------------ 
# Ejercicio 5
# ------------------------------------------------------------------ 
# Optimización Monte Carlo

# ------------------------------- 
# Apartado a 
# Optimización numérica

# Muestra
set.seed(12345)
data <- sample(rbind(rnorm(10^2), 2.5 + rnorm(3 * 10^2)))

# Logaritmo (negativo) de la función de verosimilitud
like <- function(mu) {
    -sum(log((0.25 * dnorm(data - mu[1]) + 0.75 * dnorm(data - mu[2]))))
}

# Representar la superficie del logaritmo de la verosimilitud
mu1 <- mu2 <- seq(-2, 5, le = 250)
lli <- matrix(0, nco = 250, nro = 250)

for (i in 1:250) for (j in 1:250) lli[i, j] <- like(c(mu1[i], mu2[j]))

par(mar = c(4, 4, 1, 1))
image(mu1, mu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mu1, mu2, -lli, nle = 100, add = TRUE)

# Valores iniciales
starts <- matrix(c(1, 1, -1, -1, 4, -1, 4, 2.5, 1, 4, -1, 4.5), nrow = 2)
points(t(starts), col = "blue", pch = 19)

# Minimización numérica con nlm
for (j in 1:ncol(starts)) {
    sta <- starts[, j]
    mmu <- sta
    for (i in 1:(nlm(like, sta)$it)) mmu <- rbind(mmu, nlm(like, sta, iterlim = i)$est)
    lines(mmu, lwd = 2, col = "blue")
    points(mmu[i + 1, 1], mmu[i + 1, 2], pch = 19)
    print(c(mmu[i + 1, ], like(mmu[i + 1, ])))
}


# ------------------------------- 
# Apartado b 
# Temple simulado

# Representar la superficie del logaritmo de la verosimilitud
image(mu1, mu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mu1, mu2, -lli, nle = 100, add = TRUE)

# Valores iniciales
starts <- matrix(c(1, 1, -1, -1, 4, -1, 4, 2.5, 1, 4, -1, 4.5), nrow = 2)

SA2 <- function(fun, xini, minx = -2, maxx = 5, tolerance = 10^(-4), factor = 1) {
    temp <- scale <- iter <- dif <- 1
    the <- matrix(xini, ncol = 2)
    curfun <- hval <- fun(xini)
    while (dif > tolerance) {
        prop <- the[iter, ] + rnorm(2) * scale[iter]
        if ((min(prop) < minx) || (max(prop) > maxx) || (temp[iter] * log(runif(1)) > 
            -fun(prop) + curfun)) 
            prop <- the[iter, ]
        curfun <- fun(prop)
        hval <- c(hval, curfun)
        the <- rbind(the, prop)
        iter <- iter + 1
        temp <- c(temp, 1/log(iter + 1))  # Actualizar la temperatura
        # Se controla el nº de perturbaciones aceptadas
        ace <- length(unique(the[(iter/2):iter, 1]))
        if (ace == 1) 
            # si es muy pequeño se disminuye la escala de la perturbación
        factor <- factor/10
        if (2 * ace > iter) 
            # si es muy grande se aumenta
        factor <- factor * 10
        scale <- c(scale, max(2, factor * sqrt(temp[iter])))  # Actualizar la escala de la perturbación
        dif <- (iter < 100) + (ace < 2) + (max(hval) - max(hval[1:(iter/2)]))
    }
    list(theta = the, val = hval, ite = iter)
}

set.seed(12345)
points(t(starts), col = "blue", pch = 19)

for (j in 1:ncol(starts)) {
    sar <- SA2(like, starts[, j])
    lines(sar$the[, 1], sar$the[, 2], lwd = 2, col = "blue")
    points(sar$the[sar$it, 1], sar$the[sar$it, 2], pch = 19)
}


# ------------------------------- 
# Apartado c 
# Algoritmo genético

require(DEoptim)

# Representar la superficie del logaritmo de la verosimilitud
image(mu1, mu2, -lli, xlab = expression(mu[1]), ylab = expression(mu[2]))
contour(mu1, mu2, -lli, nle = 100, add = TRUE)

# No requiere valores iniciales (los genera al azar en el rango)

# Optimización con DEoptim
lower <- c(-2, -2)
upper <- c(5, 5)
set.seed(12345)
der <- DEoptim(like, lower, upper, DEoptim.control(NP = 20, itermax = 100))

lines(der$member$bestmemit[, 1], der$member$bestmemit[, 2], lwd = 2, col = "blue")
points(der$optim$bestmem[1], der$optim$bestmem[2], pch = 19)
 
