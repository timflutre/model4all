## Aim: perform inference for the "correlated observations" theme
## Copyright (C) 2015-2016 CIRAD, INRA
## License: GPL-3+
## Persons: Marie Denis [cre,aut], Timothee Flutre [ctb]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

corrobs.infer.name <- "corrobs_infer"
corrobs.infer.version <- "0.3.0" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
corrobsInferHelp <- function(){
  txt <- paste0("`", corrobs.infer.name, "' performs inference for the \"correlated observations\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --pkg\tpackage used to fit the model\n")
  txt <- paste0(txt, "\t\tamong: INLA, nlme, rstan\n")
  txt <- paste0(txt, "      --model\ttype of model for fitting (default=VAR1/DLM1)\n")
  txt <- paste0(txt, "      --method\ttype of estimation for nlme (default=REML)\n")
  txt <- paste0(txt, "      --co\tcompile only (for rstan)\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <marie.denis@cirad.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
corrobsInferVersion <- function(){
  txt <- paste0(corrobs.infer.name, " ", corrobs.infer.version, "\n")
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
##' @author Timothee Flutre, Marie Denis
corrobsInferParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      corrobsInferHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      corrobsInferVersion()
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
    else if(prog.args[i] == "--model"){
      prog.opts$model <- prog.args[i+1]
      i <- i + 1
    }
    else if(prog.args[i] == "--meth"){
      prog.opts$method <- prog.args[i+1]
      i <- i + 1
    }
    else if(prog.args[i] == "--co"){
      prog.opts$compile.only <- TRUE
    }
    else{
      write(paste0(corrobs.infer.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      corrobsInferHelp()
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
##' @author Timothee Flutre, Marie Denis
corrobsInferCheckOptions <- function(prog.opts){
  if(! prog.opts$model %in% c("VAR1", "DLM1")){
    msg <- paste0("ERROR: unknown --model ", prog.opts$model)
    write(msg, stderr())
    quit(save="no", status=1)
  }
  if(is.null(prog.opts$pkg)){
    write("ERROR: missing compulsory option --pkg", stderr())
    quit("no", status=1)
  }
  if(! prog.opts$pkg %in% c("INLA", "nlme", "rstan")){
    write(paste0("ERROR: unknown package '", prog.opts$pkg, "'"), stderr())
    quit("no", status=1)
  }

  suppressPackageStartupMessages(library(tools)) # for file_path_sans_ext()
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(prog.opts$pkg, character.only=TRUE))

  return(prog.opts)
}

##' Fit corrobs with nlme
##'
##'
##' @param data
##' @param method (by default REML)
##' @return list
##' @author Marie Denis
corrobsInfer_nlme <- function(data, method){
  gls.res = gls(y~1,
                data = data,
                correlation = corAR1(form = ~ 1|I ),
                method = method)
  return(gls.res)
}

##' Fit corrobs with INLA
##'
##' integrated nested Laplace approximation
##' http://www.r-inla.org/comments-1?place=msg%2Fr-inla-discussion-group%2FiDmuCF6dp6I%2F8KLAikmayPMJ
##' @param data
##' @param verbose
##' @return list
##' @author Marie Denis
corrobsInfer_INLA <- function(data, model = "VAR1", verbose=0){
  fit <- NULL

  if (model == "VAR1"){
    N <- length(table(data$I))
    T <- table(data$I)[1]
    zeta <- rep(1:T, N)
    id <- rep(1:N, each=T)
    fit <- inla(formula= y ~ f(zeta, model="ar1", replicate=id),
                family="gaussian",
                data=data.frame(y=data$y, zeta, id),
                control.compute=list(dic=TRUE),
                control.family=list(initial=10, fixed=TRUE),
                verbose=ifelse(verbose > 0, TRUE, FALSE),
                silent=TRUE)
  } else{
    id <- 1:length(data$y)
    formula <- y ~ f(id, model="ar1") - 1
    fit <- inla(formula, family="gaussian",
                data=data.frame(y=data$y, id),
                control.predictor=list(compute=TRUE),
                control.family = list(hyper=list(prec=list(prior= "gaussian",
                                                           param=c(0,1)))))
  }

  return(fit)
}

##'
##'
##' @param stan.file
##' @return nothing
##' @author Timothee Flutre, Marie Denis
corrobsInfer_rstan_writeModel <- function(stan.file, model){
  ## meta-data
  model.code <- "# model: quantgen
# copyright: INRA
# license: GPL-3
# persons: Timothée Flutre [cre,aut]"
  model.code <- paste0(model.code, "
# date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
  if ( model == "VAR1"){
  ## block: data
  model.code <- paste0(model.code, "
data {
    int<lower=0> T; //num of observations per individual
    int N;          //num of individuals
    matrix[T,N] y;  //the matrix of observations")
  model.code <- paste0(model.code, "
")
  model.code <- paste0(model.code, "}\n")

## block: parameters
  model.code <- paste0(model.code, "
parameters {
    real alpha;
    real rho;
     real<lower=0> sigma;")
model.code <- paste0(model.code, "
}
")
## block: model
  model.code <- paste0(model.code, "
model {")
  model.code <- paste0(model.code, "
              alpha ~ normal(0,1000);
              rho ~ normal(0,1000);
              sigma ~ inv_gamma(0.001,0.001);
              for (n in 1:N)
                for (t in 2:T)
                  y[t,n] ~ normal(alpha + rho * y[t-1,n], sigma);
")
}else{
  ## block: data
  model.code <- paste0(model.code, "
data {
                       int<lower=0> T; //num of observations per individual
                       vector[T] y;  //the matrix of observations")
  model.code <- paste0(model.code, "
                       ")
  model.code <- paste0(model.code, "}\n")

## block: parameters
  model.code <- paste0(model.code, "
                       parameters {
                        real rho;
                        vector[T] x;
                        real<lower=0> sigma_y;
                        real<lower=0> sigma_x;")

model.code <- paste0(model.code, "
                       }
                     ")
## block: model
  model.code <- paste0(model.code, "
                       model {")
  model.code <- paste0(model.code, "
                  rho ~ normal(0,5);
                  sigma_y ~ inv_gamma(0.001,0.001);
                  sigma_x ~ inv_gamma(0.001,0.001);
                  for (i in 2:T)
                    x[i] ~ normal(rho*x[i-1], sigma_x);
                  y ~ normal(x, sigma_y);
                       ")

}
model.code <- paste0(model.code, "}")
cat(model.code, file=stan.file)
}

##' Fit corrobs with rstan
##'
##' HMC
##' @param infer.dir
##' @param compile.only
##' @param y
##' @param data
##' @param nb.iters
##' @param burnin
##' @param thin
##' @param verbose
##' @return list
##' @author Timothee Flutre, Marie Denis
corrobsInfer_rstan <- function(infer.dir, compile.only, model,
                               data, nb.iters, burnin, thin, verbose=0){
  fit <- NULL

  stan.file <- paste0(infer.dir, "/corrobs.stan")
  corrobsInfer_rstan_writeModel(stan.file, model)

  file.sm <- paste0(infer.dir, "/corrobs_rstan_sm.RData")
  if(compile.only){
    write(paste0("compile with rstan ..."), stdout())
    rt <- stanc(file=stan.file, model_name="corrobs")
    sm <- stan_model(stanc_ret=rt,
                     verbose=ifelse(verbose > 0, TRUE, FALSE))
    save(sm, file=file.sm)
  } else{
    if(model == "VAR1"){
      T <- nrow(data$Y)
      N <- ncol(data$Y)
      ldat <- list(T=T,
                   N=N,
                   y=data$Y)
    } else{
      ldat <- list(T=length(data$data$y),
                   y=data$data$y)
    }
    if(file.exists(file.sm)){
      load(file.sm)
      fit <- sampling(sm,
                      data=ldat,
                      iter=nb.iters + burnin,
                      warmup=burnin,
                      thin=thin,
                      chains=1)
    } else
      fit <- stan(file=stan.file,
                  data=ldat,
                  iter=nb.iters + burnin,
                  warmup=burnin,
                  thin=thin,
                  chains=1)
  }

  return(fit)
}

##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param infer.dir character
##' @param infer.file character
##' @return nothing
##' @author Timothee Flutre, Marie Denis
corrobsInferRun <- function(prog.opts, simul.file, infer.dir, infer.file){
  load(simul.file)
  sid <- strsplit(x=file_path_sans_ext(basename(simul.file)),
                  split="_")[[1]][3]

  pkg.ver <- packageVersion(prog.opts$pkg)

  fit <- NULL

  if(prog.opts$pkg == "nlme"){
    st <- system.time(fit <- corrobsInfer_nlme(
                          data= sim$data,
                          method=prog.opts$method))
  }

  if(prog.opts$pkg == "INLA"){
    st <- system.time(fit <- corrobsInfer_INLA(
                          data=sim$data,
                          model=prog.opts$model,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "rstan"){
    st <- system.time(fit <- corrobsInfer_rstan(
                          infer.dir=infer.dir,
                          compile.only=prog.opts$compile.only,
                          model = prog.opts$model,
                          data=sim,
                          nb.iters=prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin,
                          verbose=prog.opts$verbose-1))
  }

  if(! prog.opts$co){
    if(prog.opts$verbose > 0)
      write("save into file ...", stdout())
    save(fit, pkg.ver, st, file=paste0(infer.dir, "/", infer.file))
  }
}

##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param infer.dir character
##' @param infer.file character
##' @return nothing
##' @author Timothee Flutre, Marie Denis
corrobsInferMain <- function(prog.args, simul.file, infer.dir, infer.file,
                             nb.chains=2, nb.iters=10^3, burnin=10^2, thin=5){
  prog.opts <- list(verbose=1,
                    pkg=NULL,
                    nb.chains=nb.chains,
                    nb.iters=nb.iters,
                    burnin=burnin,
                    thin=thin,
                    model="VAR1",
                    method="REML",
                    compile.only=FALSE)

  prog.opts <- corrobsInferParseArgs(prog.args, prog.opts)

  prog.opts <- corrobsInferCheckOptions(prog.opts)

  corrobsInferRun(prog.opts, simul.file, infer.dir, infer.file)
}
