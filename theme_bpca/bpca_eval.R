## Aim: evaluate inference for the "bpca" theme
## Copyright (C) 2015-2016 INRA
## License: AGPL-3+
## Persons: Gabrielle Weinrott [aut]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

bpca.eval.name <- "bpca_eval"
bpca.eval.version <- "0.1.0" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
bpcaEvalHelp <- function(){
  txt <- paste0("`", bpca.eval.name, "' evaluates for the \"bpca\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --pkg\tpackage used to fit the model\n")
  txt <- paste0(txt, "\t\tamong: rjags, rstan\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <gabrielle.weinrott@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
bpcaEvalVersion <- function(){
  txt <- paste0(bpca.eval.name, " ", bpca.eval.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License AGPLv3+: GNU AGPL version 3 or later <http://gnu.org/licenses/agpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Gabrielle Weinrott [aut].")
  write(txt, stdout())
}

##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Gabrielle Weinrott
bpcaEvalParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      bpcaEvalHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      bpcaEvalVersion()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-v" || prog.args[i] == "--verbose"){
      prog.opts$verbose <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--pkg"){
      prog.opts$pkg <- prog.args[i+1]
      i <- i + 1
    }
    else{
      write(paste0(bpca.eval.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      bpcaEvalHelp()
      quit("no", status=1)
    }
  }

  return(prog.opts)
}

##' ##' Check program options
##'
##'
##' @param prog.opts named list
##' @return nothing
##' @author Gabrielle Weinrott
bpcaEvalCheckOptions <- function(prog.opts){
  if(is.null(prog.opts$pkg)){
    write("ERROR: missing compulsory option --pkg", stderr())
    quit("no", status=1)
  }
  if(! prog.opts$pkg %in% c("rstan", "rjags")){
    write(paste0("ERROR: unknown package '", prog.opts$pkg, "'"), stderr())
    quit("no", status=1)
  }

  suppressPackageStartupMessages(library(tools)) # for file_path_sans_ext()
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(coda))
  suppressPackageStartupMessages(library(prog.opts$pkg, character.only=TRUE))
  suppressPackageStartupMessages(library(mcmcse))
}

##' Perform the evaluation(s)
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param infer.file character
##' @param eval.file character
##' @return nothing
##' @author Gabrielle Weinrott
bpcaEvalRun <- function(prog.opts, simul.file, infer.file, eval.file){

  load(simul.file)
  sid <- strsplit(x=file_path_sans_ext(basename(simul.file)),
                  split="_")[[1]][3]
  load(infer.file)
  infer.dir <- dirname(infer.file)

  if(prog.opts$verbose > 0)
    write(paste0("evaluate ", prog.opts$pkg, " ..."), stdout())

  if(prog.opts$pkg == "rjags"){
    ## print(summary(fit))

    write(sprintf("niters=%.2f", nrow(fit[[1]])), stdout())

    ess <- effectiveSize(fit)

    var.names <- dimnames(fit[[1]])[[2]]
    n.params <- length(var.names)
    solve.f <- abs((1 - sqrt(1 + 4*n.params))/2)

    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste("ess(", var.names[i], ")=", round(ess[var.names[i]], 2), sep = "")),
            stdout())
    }

    for(i in 1:(solve.f^2)){
      write(sprintf(paste("ess(", var.names[i], ")=", round(ess[var.names[i]], 2), sep = "")),
              stdout())
    }

    ac <- autocorr.diag(fit)

    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste("autocorr ", rownames(ac)[2], " : ", var.names[i], "=",
                          round(ac[2, var.names[i]], 2), sep = "")),
            stdout())
    }

    for(i in 1:(solve.f^2)){
      write(sprintf(paste("autocorr ", rownames(ac)[2], " : ", var.names[i], "=",
                          round(ac[2, var.names[i]], 2), sep = "")),
            stdout())
    }

    # batch means
    mcse.post <- c()
    mcse.se <- c()
    for(i in 1:ncol(fit[[1]])){
      fleg <- mcse(fit[[1]][ , i])
      mcse.post[i] <- fleg$est
      mcse.se[i] <- fleg$se
    }

  # mcse.mat(fit[[1]])

    write("posterior estimates for mu...", stdout())
    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste(var.names[i], " = ",
                          round(mcse.post[i], 2),
                          " +/- ",
                          round(mcse.se[i], 3), sep = "")),
            stdout())
    }

    write("posterior estimates for Sigma...", stdout())
    for(i in 1:(solve.f^2)){
      write(sprintf(paste(var.names[i], " = ",
                          round(mcse.post[i], 2),
                          " +/- ",
                          round(mcse.se[i], 3), sep = "")),
            stdout())
    }
  }

  if(prog.opts$pkg == "rstan"){
    ## monitor() requires >1 chain, thus use the coda package instead

    coda.obj <- function(fit){
      codafit <- mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
      return(codafit)
    }

    fit <- coda.obj(fit)

    write(sprintf("niters=%.2f", nrow(fit[[1]])), stdout())

    ess <- effectiveSize(fit)

    var.names <- dimnames(fit[[1]])[[2]]
    n.params <- length(var.names)
    solve.f <- abs((1 - sqrt(1 + 4*n.params))/2)

    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste("ess(", var.names[i], ")=", round(ess[var.names[i]], 2), sep = "")),
            stdout())
    }

    for(i in 1:(solve.f^2)){
      write(sprintf(paste("ess(", var.names[i], ")=", round(ess[var.names[i]], 2), sep = "")),
            stdout())
    }

    ac <- autocorr.diag(fit)

    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste("autocorr ", rownames(ac)[2], " : ", var.names[i], "=",
                          round(ac[2, var.names[i]], 2), sep = "")),
            stdout())
    }

    for(i in 1:(solve.f^2)){
      write(sprintf(paste("autocorr ", rownames(ac)[2], " : ", var.names[i], "=",
                          round(ac[2, var.names[i]], 2), sep = "")),
            stdout())
    }

    # batch means
    mcse.post <- c()
    mcse.se <- c()
    for(i in 1:ncol(fit[[1]])){
      fleg <- mcse(fit[[1]][ , i])
      mcse.post[i] <- fleg$est
      mcse.se[i] <- fleg$se
    }

    # mcse.mat(fit[[1]])

    write("posterior estimates for mu...", stdout())
    for(i in (solve.f^2 + 1):length(var.names)){
      write(sprintf(paste(var.names[i], " = ",
                          round(mcse.post[i], 2),
                          " +/- ",
                          round(mcse.se[i], 3), sep = "")),
            stdout())
    }

    write("posterior estimates for Sigma...", stdout())
    for(i in 1:(solve.f^2)){
      write(sprintf(paste(var.names[i], " = ",
                          round(mcse.post[i], 2),
                          " +/- ",
                          round(mcse.se[i], 3), sep = "")),
            stdout())
    }
  }
  }

##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param infer.file character
##' @param eval.file character
##' @return nothing
##' @author Gabrielle Weinrott
bpcaEvalMain <- function(prog.args, simul.file, infer.file, eval.file){
  prog.opts <- list(verbose=1, pkg=NULL)

  prog.opts <- bpcaEvalParseArgs(prog.args, prog.opts)

  bpcaEvalCheckOptions(prog.opts)

  bpcaEvalRun(prog.opts, simul.file, infer.file, eval.file)
  }
