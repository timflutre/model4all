## Aim: perform inference for the "bpca" theme
## Copyright (C) 2015-2016 INRA
## License: AGPL-3+
## Persons: Gabrielle Weinrott [cre,aut], Timothee Flutre [ctb]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

bpca.infer.name <- "bpca_infer"
bpca.infer.version <- "0.1.0" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
bpcaInferHelp <- function(){
  txt <- paste0("`", bpca.infer.name, "' performs inference for the \"latent factor\" theme.\n")
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
bpcaInferVersion <- function(){
  txt <- paste0(bpca.infer.name, " ", bpca.infer.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License AGPLv3+: GNU AGPL version 3 or later <http://gnu.org/licenses/agpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Gabrielle Weinrott [cre,aut], Timothee Flutre [ctb].")
  write(txt, stdout())
}

##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Gabrielle Weinrott
bpcaInferParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      bpcaInferHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      bpcaInferVersion()
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
      write(paste0(bpca.infer.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      bpcaInferHelp()
      quit("no", status=1)
    }
  }

  return(prog.opts)
}

##' Check program options
##'
##'
##' @param prog.opts named list
##' @return nothing
##' @author Gabrielle Weinrott
bpcaInferCheckOptions <- function(prog.opts){
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
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(prog.opts$pkg, character.only=TRUE))
  suppressPackageStartupMessages(library(MASS))
  suppressPackageStartupMessages(library(coda))

  return(prog.opts)
}

##'
##'
##' @param jags.file
##' @return nothing
##' @author Gabrielle Weinrott
bpcaInfer_rjags_writeModel <- function(jags.file){
  ## meta-data
  model.code <- "# model: bpca
# author: Gabrielle Weinrott (INRA)"
  model.code <- paste0(model.code, "
# written: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")

  ## block: model
  model.code <- paste0(model.code, "
model {
  ## priors
  mu[1:P] ~ dmnorm(mu.prior[], mu.prior.prec[,])
  prec[1:P,1:P] ~ dwish(covmat.prior[,], covmat.prior.DF)
  cov[1:P,1:P] <- inverse(prec[,])

  ## likelihood
  for (i in 1:N){
    Y[i,1:P] ~ dmnorm(mu[], prec[,])
  }
")

  model.code <- paste0(model.code, "}")

  cat(model.code, file = jags.file)
}

##' Fit a Bayesian PCA JAGS model to a data set
##'
##'
##' @param data matrix of simulated data
##' @param nchains number of chains
##' @param niter number of iterations
##' @param burnin number of burnin iterations
##' @param thin
##' @return fit, a BUGS object
##' @author Gabrielle Weinrott
bpcaInfer_rjags <- function(infer.dir,
                            fit.data,
                            nchains,
                            niter,
                            burnin,
                            thin){
  fit <-  NULL

  N = nrow(fit.data)
  P = ncol(fit.data)

  covmat.prior = as.matrix(diag(1/1000, P))
  covmat.prior.DF = P
  mu.prior.cov = as.matrix(diag(1000, P))
  mu.prior = rep(0, P)

  mu.prior.prec = ginv(mu.prior.cov)

  jags.file <- paste0(infer.dir, "/bpca.jags")
  bpcaInfer_rjags_writeModel(jags.file)

  file.sm <- paste0(infer.dir, "/bpca_rjags_sm.RData")

  jags <- jags.model(file = jags.file,
                     data = list(Y = as.matrix(fit.data),
                                 N = N,
                                 P = P,
                                 covmat.prior = covmat.prior,
                                 mu.prior = mu.prior,
                                 covmat.prior.DF = covmat.prior.DF,
                                 mu.prior.prec = mu.prior.prec),
	                   inits = NULL,
	                   n.chains = nchains,
	                   n.adapt = burnin)

  fit <- coda.samples(model = jags,
		      variable.names = c("cov", "mu"),
		      n.iter = niter,
		      thin = thin)

  return(fit)
}

##'
##'
##' @param stan.file
##' @return nothing
##' @author Gabrielle Weinrott
bpcaInfer_rstan_writeModel <- function(stan.file){
  ## meta-data
  model.code <- "# model: bpca
# copyright: INRA
# license: AGPL-3
# persons: Gabrielle Weinrott [cre,aut]"
  model.code <- paste0(model.code, "
# date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")

  ## block: data
  model.code <- paste0(model.code, "
data {
  int<lower=1> N;
  int<lower=1> P;
  cov_matrix[P] covmatprior;
  vector[P] dat[N];
}

parameters {
  cov_matrix[P] prec;
  vector[P] mu;
}

model {
  mu ~ normal(0, 1000);
  prec ~ inv_wishart(P, covmatprior);
  for (n in 1:N){
    dat[n] ~ multi_normal(mu, prec);
  }
}")

  cat(model.code, file=stan.file)
}

##' Fit a Bayesian Latent Factor STAN model to a data set
##'
##'
##' @param data matrix of simulated data (projected onto the histogram basis)
##' @param nchains number of chains
##' @param niter number of iterations
##' @param burnin number of burnin iterations
##' @param thin
##' @return stanfit, a STAN object
##' @author Gabrielle Weinrott
bpcaInfer_rstan <- function(infer.dir, fit.data,
                            nchains, niter, burnin, thin){
  fit <-  NULL

  if(is.null(fit.data)){
    stop("No data specified!")
  }

  rstan_options(auto_write = TRUE)
  ## options(mc.cores = parallel::detectCores())

  N <- nrow(fit.data)
  P <- ncol(fit.data)

  mu.prior.prec = diag(1000, P)
  covmat.prior = diag(1000, P)

  stan_data <- list(P=P, N=N, dat = fit.data, covmatprior = covmat.prior)

  stan.file <- paste0(infer.dir, "/bpca.stan")

  bpcaInfer_rstan_writeModel(stan.file)

  fit <- stan(file = stan.file,
              dat = stan_data,
              iter = niter,
              warmup = burnin,
              thin = thin,
              chains = nchains)

  return(fit)
}

##' Main function to fit the Bayesian Latent Factor model
##'
##' @param prog.opts a list of parameters for the modelFit function
##' @param out.file the path to the file that will be saved
##' @param verbose binary variable (default = 1)
##' @param obj.version output version (eg. stanfit1)
##' @return nothing
##' @author Gabrielle Weinrott
bpcaInferRun <- function(prog.opts, simul.file, infer.dir, infer.file){
  load(simul.file)
  sid <- strsplit(x=file_path_sans_ext(basename(simul.file)),
                  split="_")[[1]][3]

  pkg.ver <- packageVersion(prog.opts$pkg)
  if(prog.opts$verbose > 0)
    write(paste0("infer with ", prog.opts$pkg, " (v", pkg.ver,
                 ")..."),
          stdout())

  fit <- NULL

  if(prog.opts$pkg == "rstan"){
    st <- system.time(fit <- bpcaInfer_rstan(
                          infer.dir = infer.dir,
                          fit.data = data,
                          nchains = prog.opts$nb.chains,
                          niter = prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin)
                      )
  }

  if(prog.opts$pkg == "rjags"){
    st <- system.time(fit <- bpcaInfer_rjags(
                          infer.dir = infer.dir,
                          fit.data = data,
                          nchains = prog.opts$nchains,
                          niter = prog.opts$niter,
                          burnin = prog.opts$burnin,
                          thin = prog.opts$thin)
                      )
  }

  if(prog.opts$verbose > 0)
    write("save into file ...", stdout())
  save(fit, pkg.ver, st, file=paste0(infer.dir, "/", infer.file))
}

##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param infer.dir character
##' @param infer.file character
##' @return nothing
##' @author Gabrielle Weinrott
bpcaInferMain <- function(prog.args, simul.file, infer.dir, infer.file,
                          nb.chains=2, nb.iters=10^3, burnin=10^2, thin=5){
  prog.opts <- list(verbose=1,
                    pkg=NULL,
                    nb.chains=nb.chains,
                    nb.iters=nb.iters,
                    burnin=burnin,
                    thin=thin)

  prog.opts <- bpcaInferParseArgs(prog.args, prog.opts)

  prog.opts <- bpcaInferCheckOptions(prog.opts)

  bpcaInferRun(prog.opts, simul.file, infer.dir, infer.file)
}
