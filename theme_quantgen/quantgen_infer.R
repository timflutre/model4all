## Aim: perform inference for the "quantitative genetics" theme
## Copyright (C) 2015-2016 INRA
## License: GPL-3+
## Persons: Timothée Flutre [cre,aut]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

quantgen.infer.name <- "quantgen_infer"
quantgen.infer.version <- "0.5.2" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
quantgenInferHelp <- function(){
  txt <- paste0("`", quantgen.infer.name, "' performs inference for the \"quantitative genetics\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --pkg\tpackage used to fit the model\n")
  txt <- paste0(txt, "\t\tamong: rrBLUP, lme4, MCMCglmm, rjags, rstan, INLA, BGLR, np\n")
  txt <- paste0(txt, "      --dom\tuse dominant relationships\n")
  txt <- paste0(txt, "\t\tby default, only additive relationships are used\n")
  txt <- paste0(txt, "      --mark\texplicitly use markers' genotypes for inference\n")
  txt <- paste0(txt, "      --priormk\tprior on the markers' effects\n")
  txt <- paste0(txt, "\t\tBGLR: BRR/BayesA/BL/BayesB/BayesC\n")
  txt <- paste0(txt, "      --co\tcompile only (for rstan)\n")
  txt <- paste0(txt, "      --errst\tuse Student's t for errors (for rstan)\n")
  txt <- paste0(txt, "      --dbwsel\tuse default (i.e. slow) bandwith selection (for np)\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <timothee.flutre@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
quantgenInferVersion <- function(){
  txt <- paste0(quantgen.infer.name, " ", quantgen.infer.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothée Flutre [cre,aut].")
  write(txt, stdout())
}

##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Timothee Flutre
quantgenInferParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      quantgenInferHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      quantgenInferVersion()
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
    else if(prog.args[i] == "--dom"){
      prog.opts$use.dominant <- TRUE
      i <- i + 1
    }
    else if(prog.args[i] == "--mark"){
      prog.opts$use.markers <- TRUE
    }
    else if(prog.args[i] == "--priormk"){
      prog.opts$prior.mark <- prog.args[i+1]
      i <- i + 1
    }
    else if(prog.args[i] == "--co"){
      prog.opts$compile.only <- TRUE
    }
    else if(prog.args[i] == "--errst"){
      prog.opts$errors.Student <- TRUE
    }
    else if(prog.args[i] == "--dbwsel"){
      prog.opts$default.bw.sel <- TRUE
    }
    else{
      write(paste0(quantgen.infer.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      quantgenInferHelp()
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
##' @author Timothee Flutre
quantgenInferCheckOptions <- function(prog.opts){
  if(is.null(prog.opts$pkg)){
    write("ERROR: missing compulsory option --pkg", stderr())
    quit("no", status=1)
  }
  if(! prog.opts$pkg %in% c("rrBLUP", "rstan", "lme4", "MCMCglmm", "INLA",
                            "R2OpenBUGS", "rjags", "BGLR", "np")){
    write(paste0("ERROR: unknown package '", prog.opts$pkg, "'"), stderr())
    quit("no", status=1)
  }
  if(prog.opts$use.markers){
    if(! prog.opts$pkg %in% c("BGLR", "rrBLUP", "rstan", "np")){
      write(paste0("ERROR: package '", prog.opts$pkg,
                   "' can't explicitly handle markers' genotypes (yet)"),
            stderr())
      quit("no", status=1)
    }
  }
  if(prog.opts$pkg == "BGLR" & is.null(prog.opts$prior.mark)){
    prog.opts$prior.mark <- "BRR"
  }

  suppressPackageStartupMessages(library(tools)) # for file_path_sans_ext()
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(prog.opts$pkg, character.only=TRUE))

  return(prog.opts)
}

##' Fit quantgen with BGLR
##'
##'
##' @param y
##' @param W
##' @param Z
##' @param X
##' @param prior.mark
##' @param nb.iters
##' @param burnin
##' @param thin
##' @param saveAt
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_BGLR <- function(y, W, Z, X, prior.mark,
                               nb.iters, burnin, thin,
                               saveAt, verbose=0){
  fit <- BGLR(y=y,
              response_type="gaussian",
              ETA=list(list(X=W[,-1], model="FIXED"),
                       list(X=Z %*% X, model=prior.mark)),
              nIter=nb.iters,
              burnIn=burnin,
              thin=thin,
              saveAt=saveAt,
              verbose=ifelse(verbose > 0, TRUE, FALSE))
  return(fit)
}

##' Fit quantgen with INLA
##'
##' integrated nested Laplace approximation
##' http://www.r-inla.org/comments-1?place=msg%2Fr-inla-discussion-group%2FiDmuCF6dp6I%2F8KLAikmayPMJ
##' @param data
##' @param A
##' @param D
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_INLA <- function(data, A, D, verbose=0){
  fit <- NULL

  tmp <- data
  colnames(tmp)[colnames(tmp) == "geno"] <- "geno.add"

  if(is.null(D)){
    fit <- rutilstimflutre::inlaAM(data=tmp,
                                   relmat=list(geno.add=A),
                                   verbose=ifelse(verbose > 0, TRUE, FALSE),
                                   silent=TRUE)
  } else{
    tmp$geno.dom <- tmp$geno.add
    fit <- rutilstimflutre::inlaAM(data=tmp,
                                   relmat=list(geno.add=A, geno.dom=D),
                                   verbose=ifelse(verbose > 0, TRUE, FALSE),
                                   silent=TRUE)
  }

  return(fit)
}

##' Fit quantgen with lme4
##'
##' REML
##' @param data
##' @param A
##' @param D
##' @param meth.ci
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_lme4 <- function(data, A, D=NULL, meth.ci="profile", verbose=0){
  fit <- NULL

  tmp <- data
  colnames(tmp)[colnames(tmp) == "geno"] <- "geno.add"

  if(is.null(D)){
    fit <- rutilstimflutre::lmerAM(
        formula=response1 ~ 1 + year + (1|geno.add),
        data=tmp, relmat=list(geno.add=A),
        ci.meth=meth.ci,
        verbose=verbose)
  } else{
    tmp$geno.dom <- tmp$geno.add
    fit <- rutilstimflutre::lmerAM(
        formula=response1 ~ 1 + year + (1|geno.add) + (1|geno.dom),
        data=tmp, relmat=list(geno.add=A, geno.dom=D),
        ci.meth=meth.ci,
        verbose=verbose)
  }

  return(fit)
}

##' Fit quantgen with MCMCglmm
##'
##' MCMC (mostly Gibbs sampling, slice sampling possible for binary responses)
##' priors: Normal for fixed effects, inv-Wishart for (co)variances, scaled non-central F possible, improper possible
##' aim: store 1,000-2,000 iterations and autocorr(1) less than 0.1
##' @param data
##' @param A
##' @param D
##' @param nb.iters
##' @param burnin
##' @param thin
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_MCMCglmm <- function(data, A, D, nb.iters, burnin, thin,
                                   verbose=0){
  Q <- nlevels(data$year)

  tmp <- data
  colnames(tmp)[colnames(tmp) == "geno"] <- "geno.add"

  if(is.null(D)){
    fit <- MCMCglmm(fixed=response1 ~ 1 + year,
                    random=~geno.add,
                    rcov=~units,
                    family="gaussian",
                    data=tmp,
                    start=NULL,
                    prior=list(B=list(mu=c(mean(tmp$response1), rep(0,Q-1)),
                                      V=10^8 * diag(Q)),
                               ## var ~ IG(0.001,0.001)
                               G=list(G1=list(V=1, nu=0.002)),
                               R=list(V=1, nu=0.002)),
                    ## ginverse=list(geno.add=Matrix(solve(A), sparse=TRUE)),
                    nitt=nb.iters,
                    burnin=burnin,
                    thin=thin,
                    pr=TRUE,
                    DIC=TRUE,
                    verbose=ifelse(verbose > 0, TRUE, FALSE))
  } else{
    tmp$geno.dom <- tmp$geno.add
    fit <- MCMCglmm(fixed=response1 ~ 1 + year,
                    random=~geno.add + geno.dom,
                    rcov=~units,
                    family="gaussian",
                    data=tmp,
                    start=NULL,
                    prior=list(B=list(mu=c(mean(tmp$response1), rep(0,Q-1)),
                                      V=10^8 * diag(Q)),
                               ## var ~ IG(0.001,0.001)
                               ## G=list(G1=list(V=1, nu=0.002)),
                               ## stdev ~ half-Cauchy(loc=0, scale=5)
                               G=list(G1=list(V=1, nu=1, alpha.mu=0,
                                              alpha.V=5^2),
                                      G2=list(V=1, nu=1, alpha.mu=0,
                                              alpha.V=5^2)),
                               R=list(V=1, nu=0.002)),
                    ginverse=list(geno.add=Matrix(solve(A), sparse=TRUE),
                                  geno.dom=Matrix(solve(D), sparse=TRUE)),
                    nitt=nb.iters,
                    burnin=burnin,
                    thin=thin,
                    pr=TRUE,
                    DIC=TRUE,
                    verbose=ifelse(verbose > 0, TRUE, FALSE))
  }

  return(fit)
}

##' Fit quantgen with np
##'
##'
##' @param infer.dir
##' @param y
##' @param Z
##' @param X
##' @param default.bw.sel
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_np <- function(infer.dir, y, Z, X, default.bw.sel,
                             verbose=0){
  fit <- NULL

  bw <- npregbw(ydat=y[,1], xdat=Z %*% X,
                regtype="lc", # local-constant estimator (Nadaraya-Watson)
                bwmethod="cv.ls",
                ftol=ifelse(default.bw.sel, 10 * sqrt(.Machine$double.eps), 0.01),
                tol=ifelse(default.bw.sel, 10^4 * sqrt(.Machine$double.eps), 0.01),
                nmulti=ifelse(default.bw.sel, min(5,ncol(X)), 3))
  print(str(bw))

  pred <- npreg(bws=bw, gradients=TRUE, residuals=TRUE)

  fit <- list(bw=bw, pred=pred)
  return(fit)
}

##' Fit quantgen with rjags
##'
##' MCMC (Gibbs)
##' @param infer.dir
##' @param data
##' @param W
##' @param Z
##' @param A
##' @param D
##' @param nb.chains
##' @param nb.iters
##' @param burnin
##' @param thin
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_rjags <- function(infer.dir, data, W, Z, A, D,
                                nb.chains, nb.iters, burnin, thin,
                                verbose=0){
  fit <- NULL

  tmp <- data
  colnames(tmp)[colnames(tmp) == "geno"] <- "geno.add"

  if(is.null(D)){
    fit <- rutilstimflutre::jagsAM(data=tmp,
                                   relmat=list(geno.add=A),
                                   nb.chains=nb.chains, burnin=burnin,
                                   nb.iters=nb.iters, thin=thin,
                                   verbose=verbose)
  } else{
    tmp$geno.dom <- tmp$geno.add
    fit <- rutilstimflutre::jagsAM(data=tmp,
                                   relmat=list(geno.add=A, geno.dom=D),
                                   nb.chains=nb.chains, burnin=burnin,
                                   nb.iters=nb.iters, thin=thin,
                                   verbose=verbose)
  }

  return(fit)
}

##' Fit quantgen with rrBLUP
##'
##' REML
##' @param y
##' @param W
##' @param Z
##' @param method
##' @param use.markers
##' @param X
##' @param A
##' @return list
##' @author Timothee Flutre
quantgenInfer_rrBLUP <- function(model, method="REML", use.markers=FALSE,
                                 X=NULL, A=NULL){
  fit <- NULL
  if(use.markers){
    fit <- mixed.solve(y=model$Y[,1],
                       Z=model$Z %*% X,
                       X=model$W,
                       method=method, SE=TRUE, return.Hinv=TRUE)
  } else{
    fit <- mixed.solve(y=model$Y[,1],
                       Z=model$Z,
                       K=A,
                       X=model$W,
                       method=method, SE=TRUE, return.Hinv=TRUE)
  }
  return(fit)
}

##' Fit quantgen with rstan
##'
##' HMC
##' @param infer.dir
##' @param compile.only
##' @param data
##' @param A
##' @param errors.Student
##' @param nb.chains
##' @param nb.iters
##' @param burnin
##' @param thin
##' @param verbose
##' @return list
##' @author Timothee Flutre
quantgenInfer_rstan <- function(infer.dir, compile.only,
                                data, A=NULL, errors.Student,
                                nb.chains, nb.iters, burnin, thin, verbose=0){
  fit <- NULL

  tmp <- data
  colnames(tmp)[colnames(tmp) == "geno"] <- "geno.add"

  fit <- rutilstimflutre::stanAM(data=tmp,
                                 relmat=list(geno.add=A),
                                 errors.Student=errors.Student,
                                 nb.chains=nb.chains, burnin=burnin,
                                 nb.iters=nb.iters, thin=thin,
                                 out.dir=infer.dir,
                                 task.id="model4all-quantgen-infer-rstan",
                                 compile.only=compile.only,
                                 rm.stan.file=FALSE,
                                 rm.sm.file=FALSE,
                                 verbose=verbose)

  return(fit)
}

##' Fit quantgen with R2OpenBUGS
##'
##' MCMC (Gibbs)
##' @param data
##' @param W
##' @param Z
##' @param A
##' @param nb.chains
##' @param nb.iters
##' @param burnin
##' @param thin
##' @return list
##' @author Timothee Flutre
quantgenInfer_R2OpenBUGS <- function(data, W, Z, A, nb.chains, nb.iters, burnin,
                                     thin){
  N <- length(data$response1)
  Q <- nlevels(data$year)
  I <- nlevels(data$ind)
  fit <- bugs(data=list(N=N, P=P, Q=Q, W=W, Z=Z,
                        Ainv=solve(A), y=data$response1,
                        Id=diag(I), muu=rep(0,I)),
              inits=NULL,
              parameters.to.save=c("alpha", "u", "sigmau2", "sigma2"),
              model.file="quantgen_R2OpenBUGS.txt",
              working.directory=getwd(),
              n.chains=nb.chains,
              n.iter=nb.iters,
              n.burnin=burnin,
              n.thin=thin,
              DIC=TRUE)
  return(fit)
}

##' Perform inference
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param infer.dir character
##' @param infer.file character
##' @return nothing
##' @author Timothée Flutre
quantgenInferRun <- function(prog.opts, simul.file, infer.dir, infer.file){
  load(simul.file)
  sid <- strsplit(x=file_path_sans_ext(basename(simul.file)),
                  split="_")[[1]][3]
  if(exists("D"))
    if(! prog.opts$use.dominant)
      D <- NULL

  pkg.ver <- packageVersion(prog.opts$pkg)
  if(prog.opts$verbose > 0 && ! prog.opts$co)
    write(paste0("infer with ", prog.opts$pkg, " (v", pkg.ver,
                 ", markers=", prog.opts$use.markers, ") ..."),
          stdout())

  fit <- NULL

  if(prog.opts$pkg == "BGLR"){
    st <- system.time(fit <- quantgenInfer_BGLR(
                          y=model$Y[,1],
                          W=model$W,
                          Z=model$Z,
                          X=genomes$genos,
                          prior.mark=prog.opts$prior.mark,
                          nb.iters=prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin,
                          saveAt=paste0(infer.dir, "/sid_", sid, "_"),
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "INLA"){
    st <- system.time(fit <- quantgenInfer_INLA(
                          data=model$data,
                          A=A,
                          D=D,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "lme4"){
    st <- system.time(fit <- quantgenInfer_lme4(
                          data=model$data,
                          A=A,
                          D=D,
                          meth.ci="profile"))
  }

  if(prog.opts$pkg == "MCMCglmm"){
    st <- system.time(fit <- quantgenInfer_MCMCglmm(
                          data=model$data,
                          A=A,
                          D=D,
                          nb.iters=prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "np"){
    st <- system.time(fit <- quantgenInfer_np(
                          y=model$Y[,1],
                          Z=model$Z,
                          X=genomes$genos,
                          default.bw.sel=prog.opts$default.bw.sel,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "rjags"){
    st <- system.time(fit <- quantgenInfer_rjags(
                          infer.dir=infer.dir,
                          data=model$data,
                          W=model$W,
                          Z=model$Z,
                          A=A,
                          D=D,
                          nb.chains=prog.opts$nb.chains,
                          nb.iters=prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "rrBLUP"){
    st <- system.time(fit <- quantgenInfer_rrBLUP(
                          model=model,
                          use.markers=prog.opts$use.markers,
                          X=genomes$genos,
                          A=A))
  }

  if(prog.opts$pkg == "rstan"){
    st <- system.time(fit <- quantgenInfer_rstan(
                          infer.dir=infer.dir,
                          compile.only=prog.opts$compile.only,
                          data=model$data,
                          A=A,
                          errors.Student=prog.opts$errors.Student,
                          nb.chains=prog.opts$nb.chains,
                          nb.iters=prog.opts$nb.iters,
                          burnin=prog.opts$burnin,
                          thin=prog.opts$thin,
                          verbose=prog.opts$verbose-1))
  }

  if(prog.opts$pkg == "R2OpenBUGS"){
    st <- system.time(fit <- quantgenInfer_R2OpenBUGS(
                          data=model$data,
                          W=model$W,
                          Z=model$Z,
                          A=A,
                          nb.chains=prog.opts$nb.chains,
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
##' @author Timothee Flutre
quantgenInferMain <- function(prog.args, simul.file, infer.dir, infer.file,
                              nb.chains=2, nb.iters=10^3, burnin=10^2, thin=5){
  prog.opts <- list(verbose=1,
                    pkg=NULL,
                    nb.chains=nb.chains,
                    nb.iters=nb.iters,
                    burnin=burnin,
                    thin=thin,
                    use.dominant=FALSE,
                    use.markers=FALSE,
                    prior.mark=NULL,
                    compile.only=FALSE,
                    errors.Student=FALSE,
                    default.bw.sel=FALSE)

  prog.opts <- quantgenInferParseArgs(prog.args, prog.opts)

  prog.opts <- quantgenInferCheckOptions(prog.opts)

  quantgenInferRun(prog.opts, simul.file, infer.dir, infer.file)
}
