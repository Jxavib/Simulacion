# ==================================================================
# Evaluaci?n del tiempo de CPU y rendimiento
# ==================================================================

# Funci?n de prueba
funtest <- function(n) median(runif(n))

# ------------------------------------------------------------------
# Utilidades de R
# ------------------------------------------------------------------

# -------------------------------
# proc.time() permite obtener tiempo de computaci?n real y de CPU

tini <- proc.time()
# C?digo a evaluar
funtest(100000)
# Fin c?digo a evaluar
tiempo <- proc.time() - tini
tiempo

# -------------------------------
# system.time(expresi?n) permite obtener tiempo de computaci?n real y de CPU de expresi?n

system.time(funtest(100000))

# Un paquete de inter?s puede ser microbenchmark
# microbenchmark: Accurate Timing Functions
     
# -------------------------------
# Rprof(fichero): permite evaluar el rendimiento muestreando la pila en intervalos para determinar en que funciones se emplea el tiempo de cpu
# Rprof(NULL): desactiva el muestreo
# summaryRprof(fichero): muestra los resultados

Rprof("rendimiento", interval = 0.02, memory.profiling=TRUE)
funtest(100000)
Rprof(NULL)
summaryRprof("rendimiento")

# Ver ?Tidying and profiling R code? en ?Writing R Extensions?


# ------------------------------------------------------------------
# Utilidades tiempo de CPU
# ------------------------------------------------------------------
# Llamar a CPUtimeini() donde se quiere empezar a contar.
# Llamar a CPUtimeprint() para imprimir tiempo total y tiempo desde la ?ltima llamada a una de estas funciones. 

# -------------------------------
CPUtimeini <- function() {
# -------------------------------
  .tiempo.ini <<- proc.time()
  .tiempo.last <<- .tiempo.ini
}

# -------------------------------
CPUtimeprint <- function() {
# -------------------------------
  tmp <- proc.time()
  cat("\nTiempo ?ltima operaci?n:\n")
  print(tmp-.tiempo.last)
  cat("Tiempo total operaci?n:\n")
  print(tmp-.tiempo.ini)
  .tiempo.last <<- tmp
}


# ------------------------------------------------------------------
# Ejemplo
# ------------------------------------------------------------------

CPUtimeini()
funtest(100000)
CPUtimeprint()
funtest(100000)
CPUtimeprint()