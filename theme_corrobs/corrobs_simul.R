## Aim: simulate data for the "correlated observations" theme
## Copyright (C) 2015-2016 CIRAD, INRA
## License: GPL-3+
## Persons: Marie Denis [cre,aut], Timothee Flutre [ctb]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

corrobs.simul.name <- "corrobs_simul"
corrobs.simul.version <- "0.2.2" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
corrobsSimulHelp <- function(){
  txt <- paste0("`", corrobs.simul.name, "' simulates data for the \"correlated observations\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --model\ttype of model for simulation (default=VAR1/DLM1)\n")
  txt <- paste0(txt, "      --ninds\tnumber of individuals (default=100)\n")
  txt <- paste0(txt, "      --ntimes\tnumber of observed times (default=10)\n")
  txt <- paste0(txt, "      --su2\tresidual variance component for process (default=5)\n")
  txt <- paste0(txt, "      --sy2\tresidual variance component for observation (default=1)\n")
  txt <- paste0(txt, "      --rho\tcorrelation parameter (default=0.8)\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <marie.denis@cirad.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
corrobsSimulVersion <- function(){
  txt <- paste0(corrobs.simul.name, " ", corrobs.simul.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 CIRAD, INRA.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Marie Denis [cre,aut], Timothee Flutre [ctb].")
  write(txt, stdout())
}

##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Timothee Flutre
corrobsSimulParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      corrobsSimulHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      corrobsSimulVersion()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-v" || prog.args[i] == "--verbose"){
      prog.opts$verbose <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--model"){
      prog.opts$model <- prog.args[i+1]
      i <- i + 1
    }
    else if(prog.args[i] == "--ninds"){
      prog.opts$nb.inds <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--ntimes"){
      prog.opts$nb.times <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--su2"){
      prog.opts$sigma.u2 <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--sy2"){
      prog.opts$sigma.y2 <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--rho"){
      prog.opts$rho <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else{
      write(paste0(corrobs.simul.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      corrobsSimulHelp()
      quit("no", status=1)
    }
  }

  return(prog.opts)
}

##' Check program options
##'
##'
##' @param prog.opts named list
##' @return named list
##' @author Marie Denis, Timothee Flutre
corrobsSimulCheckOptions <- function(prog.opts){
  if(! prog.opts$model %in% c("VAR1", "DLM1")){
    msg <- paste0("ERROR: unknown --model ", prog.opts$model)
    write(msg, stderr())
    quit(save="no", status=1)
  }

  return(prog.opts)
}

##' VAR1
##'
##' Simulate the observation from a VAR1 model
##' @param N integer, number of individuals (default = 100)
##' @param T integer, number of observation times (default = 10)
##' @param su2 variance parameter (default = 5)
##' @param rho, correlation parameter (default =0.8)
##' @return return correlated observations for ninds individuals
##'
##' @author Marie Denis
simul.VAR1 <- function(N = 100, T=10, rho=0.8, su2=5 ){
  Y <- NULL
  for (j in 1:N){
    tmp <- arima.sim(n=T,model=list(ar=rho), sd = sqrt(su2))
    Y <- cbind(Y,tmp)
  }
  colnames(Y) = paste("Ind", 1:N)
  rownames(Y) = paste("T",1:T)
  y <- matrix(as.vector(Y),ncol=1)
  data <- data.frame(y=y, I= factor(rep(1:N,each=(T))))
  return(list(data = data, Y = Y))
}

##' DLM1
##'
##' Simulate the observation from a DLM1 model
##' @param T integer, number of observation times (default = 10)
##' @param su2 variance parameter for latent process (default = 5)
##' @param sy2 variance parameter for observations (default = 1)
##' @param rho, correlation parameter (default =0.8)
##' @return return observations and states from a first order DLM
##'
##' @author Marie Denis
simul.DLM1 <- function(T=10, rho=0.8, su2=5, sy2=1){
  w  <- rnorm(T,0,sqrt(su2))
  v  <- rnorm(T,0,sqrt(sy2))
  x  = y = rep(0,T)
  x[1] = rnorm(1,0,sqrt(su2/(1-rho^2)))
  y[1] = x[1] + v[1]
  for (t in 2:T){
    x[t] = rho*x[t-1] + w[t]
    y[t] = x[t]   + v[t]
  }
  data <- data.frame(y=y)
  return(list(data = data))
}

##' Perform the simulation(s)
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Marie Denis, Timothee Flutre
corrobsSimulRun <- function(prog.opts, simul.file, seed){
  set.seed(seed)
  if(prog.opts$verbose > 0)
    write(paste0("simulate ", prog.opts$model, " ..."), stdout())

  if(prog.opts$model == "VAR1"){
    sim <- simul.VAR1(N = prog.opts$nb.inds,
                      T = prog.opts$nb.times,
                      rho = prog.opts$rho,
                      su2 = prog.opts$sigma.u2)
  } else{
    sim <- simul.DLM1(T = prog.opts$nb.times,
                      rho = prog.opts$rho,
                      sy2 = prog.opts$sigma.y2,
                      su2 = prog.opts$sigma.u2)
  }

  if(prog.opts$verbose > 0)
    write("save into file ...", stdout())
  save(seed, sim, file=simul.file)
}

##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Thimothee Flutre
corrobsSimulMain <- function(prog.args, simul.file, seed){
  prog.opts <- list(verbose=1,
                    nb.inds=100,
                    nb.times=10,
                    model="VAR1",
                    sigma.u2=5,
                    sigma.y2=1,
                    rho=0.8)

  prog.opts <- corrobsSimulParseArgs(prog.args, prog.opts)

  prog.opts <- corrobsSimulCheckOptions(prog.opts)

  corrobsSimulRun(prog.opts, simul.file, seed)
}
