#===============================================================================
# Test chi cuadrado distribución continua
#===============================================================================

#-------------------------------------------------------------------------------
# chisq.test.cont(x, distribution, nclasses, output, nestpar,...)
#-------------------------------------------------------------------------------
# Realiza el test ji-cuadrado de bondad de ajuste para una distribución continua
# discretizando en intervalos equiprobables.
# Parámetros:
#   distribution = "norm","unif",etc
#   nclasses = floor(length(x)/5)
#   output = TRUE
#   nestpar = 0= nº de parámetros estimados
#   ... = parámetros distribución
# Ejemplo:
#   chisq.test.cont(x, distribution="norm", nestpar=2, mean=mean(x), sd=sqrt((nx-1)/nx)*sd(x))
#-------------------------------------------------------------------------------
chisq.test.cont <- function(x, distribution = "norm", nclasses = floor(length(x)/5), 
    output = TRUE, nestpar = 0, ...) {
    # Funciones distribución
    q.distrib <- eval(parse(text = paste("q", distribution, sep = "")))
    d.distrib <- eval(parse(text = paste("d", distribution, sep = "")))
    
    # Puntos de corte
    q <- q.distrib((1:(nclasses - 1))/nclasses, ...)
    tol <- sqrt(.Machine$double.eps)
    xbreaks <- c(min(x) - tol, q, max(x) + tol)
    
    # Gráficos y frecuencias
    if (output) {
        xhist <- hist(x, breaks = xbreaks, freq = FALSE, lty = 2, border = "grey50")
        curve(d.distrib(x, ...), add = TRUE)
    } else {
        xhist <- hist(x, breaks = xbreaks, plot = FALSE)
    }
    
    # Cálculo estadístico y p-valor
    O <- xhist$counts  # Equivalente a table(cut(x, xbreaks)) pero más eficiente
    E <- length(x)/nclasses
    DNAME <- deparse(substitute(x))
    METHOD <- "Pearson's Chi-squared test"
    STATISTIC <- sum((O - E)^2/E)
    names(STATISTIC) <- "X-squared"
    PARAMETER <- nclasses - nestpar - 1
    names(PARAMETER) <- "df"
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    
    # Preparar resultados
    classes <- format(xbreaks)
    classes <- paste("(", classes[-(nclasses + 1)], ",", classes[-1], "]", 
        sep = "")
    RESULTS <- list(classes = classes, observed = O, expected = E, residuals = (O - 
        E)/sqrt(E))
    if (output) {
        cat("\nPearson's Chi-squared test table\n")
        print(as.data.frame(RESULTS))
    }
    if (any(E < 5)) 
        warning("Chi-squared approximation may be incorrect")
    structure(c(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, 
        method = METHOD, data.name = DNAME), RESULTS), class = "htest")
}
