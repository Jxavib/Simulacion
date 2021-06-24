# Para instalar los paquetes empleados basta con ejecutar en la consola de R:

pkgs <- c('Rcmdr','MASS','DEoptim','nortest','tseries','geoR','copula')
# install.packages(pkgs, dependencies=TRUE)
install.packages(setdiff(pkgs, installed.packages()[,"Package"]), dependencies = TRUE)

