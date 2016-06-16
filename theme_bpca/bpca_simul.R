# Aim: wrapper to simulate data for the "bpca" theme
# Copyright (C) 2015-2016 INRA
# License: GPL-3+
# Persons: Gabrielle Weinrott [cre,aut]
# Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

bpca.simul.name <- "bpca_simul"
bpca.simul.version <- "0.1.0" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
bpcaSimulHelp <- function(){
  txt <- paste0("`", bpca.simul.name, "' simulates data for the \"Bayesian Principal Component Analysis\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --ninds\tnumber of individuals (default=10)\n")
  txt <- paste0(txt, "      --nobs\tnumber of observations (default=150)\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <gabrielle.weinrott@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
bpcaSimulVersion <- function(){
  txt <- paste0(bpca.simul.name, " ", bpca.simul.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License AGPLv3+: GNU AGPL version 3 or later <http://gnu.org/licenses/agpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Gabrielle Weinrott [cre,aut].")
  write(txt, stdout())
}


##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Gabrielle Weinrott
bpcaSimulParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      bpcaSimulHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      bpcaSimulVersion()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-v" || prog.args[i] == "--verbose"){
      prog.opts$verbose <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--ninds"){
      prog.opts$nb.inds <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--nobs"){
      prog.opts$nb.obs <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else{
      write(paste0(bpca.simul.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      bpcaSimulHelp()
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
##' @author Gabrielle Weinrott
bpcaSimulCheckOptions <- function(prog.opts){
  suppressPackageStartupMessages(library(MASS))
  return(prog.opts)
}


##' Perform the simulation(s)
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Gabrielle Weinrott
bpcaSimulRun <- function(prog.opts, simul.file, seed){
  set.seed(seed)

  if(prog.opts$verbose > 0)
    write("simulate dataset... ", stdout())

    data <- simul.OU(prog.opts$nb.inds, prog.opts$nb.obs)

  if(prog.opts$verbose > 0)
    write("save into file ...", stdout())

  save(data, file=simul.file)
}


##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Gabrielle Weinrott
bpcaSimulMain <- function(prog.args, simul.file, seed){
  prog.opts <- list(verbose=1,
                    nb.inds = 10,
                    nb.obs = 150)

  prog.opts <- bpcaSimulParseArgs(prog.args, prog.opts)

  prog.opts <- bpcaSimulCheckOptions(prog.opts)

  bpcaSimulRun(prog.opts, simul.file, seed)
}

##' Function that simulates a dataset following an Ornstein-Uhlenbeck process
##'
##' @param ninds the number of curves to simulate
##' @param nobs the number of observation times for each curve
##'
##' @return data a simulated dataset
##'
##' @author Gabrielle Weinrott
##'
simul.OU <- function(ninds, nobs){

  Posdef <- function (nobs, ev = runif(nobs, 0, ninds))
  {
    Z <- matrix(ncol=nobs, rnorm(nobs^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }

  Sigma <- Posdef(ninds)
  mu <- rep(1, ninds)
  data <- mvrnorm(nobs, mu, Sigma)
  return(data)
}
