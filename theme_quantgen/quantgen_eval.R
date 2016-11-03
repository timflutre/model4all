## Aim: evaluate inference for the "quantitative genetics" theme
## Copyright (C) 2015-2016 INRA
## License: GPL-3+
## Persons: Timothee Flutre [cre,aut]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

quantgen.eval.name <- "quantgen_eval"
quantgen.eval.version <- "0.3.0" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
quantgenEvalHelp <- function(){
  txt <- paste0("`", quantgen.eval.name, "' evaluates for the \"quantitative genetics\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --pkg\tpackage used to fit the model\n")
  txt <- paste0(txt, "      --mark\tmarkers' genotypes explicitly used for inference\n")
  txt <- paste0(txt, "\t\tamong: BGLR, INLA, lme4, MCMCglmm, rjags, rrBLUP, rstan, R2OpenBUGS\n")
  txt <- paste0(txt, "      --priormk\tprior on the markers' effects\n")
  txt <- paste0(txt, "\t\tBGLR: BRR/BayesA/BL/BayesB/BayesC\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <timothee.flutre@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
quantgenEvalVersion <- function(){
  txt <- paste0(quantgen.eval.name, " ", quantgen.eval.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothee Flutre [cre,aut].")
  write(txt, stdout())
}

##' Parse program arguments
##'
##' Allow short and long options
##' @param prog.args character vector
##' @param prog.opts named list of program options with default values
##' @return named list
##' @author Timothee Flutre
quantgenEvalParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      quantgenEvalHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      quantgenEvalVersion()
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
    else if(prog.args[i] == "--mark"){
      prog.opts$use.markers.infer <- TRUE
    }
    else if(prog.args[i] == "--priormk"){
      prog.opts$prior.mark <- prog.args[i+1]
      i <- i + 1
    }
    else{
      write(paste0(quantgen.eval.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      quantgenEvalHelp()
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
##' @author Timothee Flutre
quantgenEvalCheckOptions <- function(prog.opts){
  if(is.null(prog.opts$pkg)){
    write("ERROR: missing compulsory option --pkg", stderr())
    quit("no", status=1)
  }
  if(! prog.opts$pkg %in% c("rrBLUP", "rstan", "lme4", "MCMCglmm", "INLA", "rjags", "BGLR")){
    write(paste0("ERROR: unknown package '", prog.opts$pkg, "'"), stderr())
    quit("no", status=1)
  }

  suppressPackageStartupMessages(library(tools)) # for file_path_sans_ext()
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(coda))
  suppressPackageStartupMessages(library(prog.opts$pkg, character.only=TRUE))
}

##'
##'
##' https://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation
##' @param error
##' @param mean.y
##' @return numeric(1)
##' @author Timothee Flutre
rmse <- function(error, mean.y=NULL){
  out <- sqrt(mean(error^2))
  if(! is.null(mean.y))
    out <- out / mean.y
  return(out)
}

##' Perform the evaluation(s)
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param infer.file character
##' @param eval.file character
##' @return nothing
##' @author Timothee Flutre
quantgenEvalRun <- function(prog.opts, simul.file, infer.file, eval.file){
  load(simul.file)
  sid <- strsplit(x=file_path_sans_ext(basename(simul.file)),
                  split="_")[[1]][3]
  load(infer.file)
  infer.dir <- dirname(infer.file)

  use.markers.simul <- "betat" %in% names(model)

  pkg.ver <- packageVersion(prog.opts$pkg)
  if(prog.opts$verbose > 0)
    write(paste0("evaluate ", prog.opts$pkg, " (v", pkg.ver, ") ..."), stdout())

  write(sprintf("time: user=%.3f system=%.3f elapsed=%.3f",
                st[1], st[2], st[3]),
        stdout())

  if(prog.opts$pkg == "BGLR"){
    tmp <- data.frame(mu=scan(paste0(infer.dir, "/sid_", sid,
                                     "_mu.dat"), quiet=TRUE))
    tmp <- cbind(tmp, read.table(paste0(infer.dir, "/sid_", sid,
                                        "_ETA_1_b.dat"),
                                 header=TRUE))
    tmp <- cbind(tmp, sigma2=scan(paste0(infer.dir, "/sid_", sid,
                                         "_varE.dat"), quiet=TRUE))
    if(prog.opts$prior.mark == "BRR"){
      tmp <- cbind(tmp, sigmau2=scan(paste0(infer.dir, "/sid_", sid,
                                            "_ETA_2_varB.dat"), quiet=TRUE))
    } else if(prog.opts$prior.mark == "bayesC")
      tmp <- cbind(tmp, sigmau2=scan(paste0(infer.dir, "/sid_", sid,
                                            "_ETA_2_parBayesC.dat"), quiet=TRUE))
    tmp <- as.mcmc(tmp)

    message("convergence diagnostics:")
    write(sprintf("nb.iters=%i", nrow(tmp)), stdout())

    ess <- effectiveSize(tmp)
    write(sprintf("ess(c[1])=%.2f ess(c[2])=%.2f ess(c[3])=%.2f",
                  ess[1], ess[2], ess[3]),
          stdout())
    write(sprintf("ess(V.E)=%.2f",
                  ess["sigma2"]),
          stdout())

    ac <- autocorr.diag(tmp)
    write(sprintf("autocorr %s: c_1=%.3f c_2=%.3f c_3=%.3f",
                  rownames(ac)[2], ac[2,1], ac[2,2], ac[2,3]),
          stdout())
    write(sprintf("autocorr %s: V.E=%.3f",
                  rownames(ac)[2], ac[2,"sigma2"]),
          stdout())

    fmt <- paste("c=[%.2f %.2f %.2f]",
                 "post.mean=[%.2f %.2f %.2f]",
                 "post.sd=[%.3f %.3f %.3f]")
    write(sprintf(fmt,
                  model$c[1], model$c[2], model$c[3],
                  fit$mu,
                  fit$ETA[[1]]$b[1],
                  fit$ETA[[1]]$b[2],
                  fit$SD.mu,
                  fit$ETA[[1]]$SD.b[1],
                  fit$ETA[[1]]$SD.b[2]),
          stdout())
    fmt <- paste("V.E=%.2f",
                 "post.mean=%.2f",
                 "post.sd=%.3f")
    write(sprintf(fmt,
                  model$sigma2,
                  fit$varE,
                  fit$SD.varE),
          stdout())
    ## write(sprintf("sigma_betat^2=%.4f hat{sigma}_betat^2=%.4f",
    ##               model$sigma.a2,
    ##               fit$ETA[[2]]$varB),
    ##       stdout())
    write(sprintf("Cor(a,hat{a}): S=%.2f P=%.2f",
                  cor(model$a, fit$ETA[[2]]$b, method="spearman"),
                  cor(model$a, fit$ETA[[2]]$b, method="pearson")),
          stdout())
    true.bv <- genomes$genos %*% model$a
    infer.bv <- genomes$genos %*% fit$ETA[[2]]$b
    write(sprintf("Cor(X a, X a): S=%.2f P=%.2f",
                  cor(true.bv, infer.bv, method="spearman"),
                  cor(true.bv, infer.bv, method="pearson")),
          stdout())
  }

  if(prog.opts$pkg == "INLA"){
    message("fixed effects:")
    ## print(summary(fit))
    ## print(fit$summary.fixed)
    hpdi.c1 <- inla.hpdmarginal(0.95, fit$marginals.fixed[["(Intercept)"]])
    hpdi.c2 <- inla.hpdmarginal(0.95, fit$marginals.fixed[["year2011"]])
    hpdi.c3 <- inla.hpdmarginal(0.95, fit$marginals.fixed[["year2012"]])
    fmt <- paste0("c=[%.2f %.2f %.2f]",
                  "\npost.mean=[%.2f %.2f %.2f]",
                  "\npost.sd=[%.3f %.3f %.3f]",
                  "\nHPDi_95(post{c_1})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_2})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_3})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$C[1,1], model$C[2,1], model$C[3,1],
                  fit$summary.fixed[1, "mean"],
                  fit$summary.fixed[2, "mean"],
                  fit$summary.fixed[3, "mean"],
                  fit$summary.fixed[1, "sd"],
                  fit$summary.fixed[2, "sd"],
                  fit$summary.fixed[3, "sd"],
                  hpdi.c1[1,"low"], hpdi.c1[1,"high"],
                  hpdi.c2[1,"low"], hpdi.c2[1,"high"],
                  hpdi.c3[1,"low"], hpdi.c3[1,"high"]),
          stdout())

    message("random effects:")
    ## print(fit$summary.hyperpar)
    ## print(fit$marginals.hyperpar)
    hpdi.V.E <- inla.hpdmarginal(0.95, fit$marginals.hyperpar[["Precision for the Gaussian observations"]])
    hpdi.V.G.A <- inla.hpdmarginal(0.95, fit$marginals.hyperpar[["Precision for geno.add"]])
    fmt <- paste0("V.E=%.2f post.mean=%.2f post.sd=%.3f",
                  "\nHPDi_95(post{V.E})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.E,
                  1/fit$summary.hyperpar["Precision for the Gaussian observations", "mean"],
                  fit$summary.hyperpar["Precision for the Gaussian observations", "sd"],
                  1/hpdi.V.E[1,"high"], 1/hpdi.V.E[1,"low"]),
          stdout())
    fmt <- paste0("V.G.A=%.2f post.mean=%.2f post.sd=%.3f",
                  "\nHPDi_95(post{V.G.A})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.G.A,
                  1/fit$summary.hyperpar["Precision for geno.add", "mean"],
                  fit$summary.hyperpar["Precision for geno.add", "sd"],
                  1/hpdi.V.G.A[1,"high"], 1/hpdi.V.G.A[1,"low"]),
          stdout())
    ## print(head(fit$summary.random[["geno.add"]]))
    write(sprintf("Cor(G.A,post.mean): S=%.2f P=%.2f",
                  cor(model$G.A[,1], fit$summary.random[["geno.add"]][,"mean"],
                      method="spearman"),
                  cor(model$G.A[,1], fit$summary.random[["geno.add"]][,"mean"],
                      method="pearson")),
          stdout())
  }

  if(prog.opts$pkg == "lme4"){
    message("fixed effects:")
    fix <- fixef(fit$merMod)
    txt <- sprintf(paste("c=[%.2f %.2f %.2f]",
                         "\nhat{c}=[%.2f %.2f %.2f]"),
                   model$C[1,1], model$C[2,1], model$C[3,1],
                   fix[1], fix[2], fix[3])
    write(txt, stdout())

    se.fix <- coefficients(summary(fit$merMod))[, "Std. Error"]
    txt <- sprintf(paste("se(hat{c})=[%.2f %.2f %.2f]",
                         "\nci_95(hat{c_1})=[%.2f %.2f]",
                         "\nci_95(hat{c_2})=[%.2f %.2f]",
                         "\nci_95(hat{c_3})=[%.2f %.2f]"),
                   se.fix[1], se.fix[2], se.fix[3],
                   fit$ci["(Intercept)", "2.5 %"], fit$ci["(Intercept)", "97.5 %"],
                   fit$ci["year2011", "2.5 %"], fit$ci["year2011", "97.5 %"],
                   fit$ci["year2012", "2.5 %"], fit$ci["year2012", "97.5 %"])
    write(txt, stdout())

    message("random effects:")
    tmp <- as.data.frame(VarCorr(fit$merMod))
    write(sprintf("V.E^2=%.2f hat{V.E}^2=%.2f ci_95=[%.2f %.2f]",
                  model$V.E,
                  tmp[tmp$grp == "Residual", "vcov"],
                  fit$ci["sigma", "2.5 %"]^2, fit$ci["sigma", "97.5 %"]^2),
          stdout())
    write(sprintf("V.G.A^2=%.2f hat{V.G.A}_u^2=%.2f ci_95=[%.2f %.2f]",
                  model$V.G.A,
                  tmp[tmp$grp == "geno.add", "vcov"],
                  fit$ci["sd_(Intercept)|geno.add", "2.5 %"]^2,
                  fit$ci["sd_(Intercept)|geno.add", "97.5 %"]^2),
          stdout())
    write(sprintf("Cor(G.A,hat{G.A}): S=%.2f P=%.2f",
                  cor(model$G.A[,1], ranef(fit$merMod)$geno.add, method="spearman"),
                  cor(model$G.A[,1], ranef(fit$merMod)$geno.add, method="pearson")),
          stdout())
  }

  if(prog.opts$pkg == "MCMCglmm"){
    message("convergence diagnostics:")
    write(sprintf("nb.iters=%i", nrow(fit$Sol)), stdout())

    ess <- effectiveSize(fit$Sol)
    write(sprintf("ess(c[1])=%.2f ess(c[2])=%.2f ess(c[3])=%.2f",
                  ess[1], ess[2], ess[3]),
          stdout())
    ess <- effectiveSize(fit$VCV)
    write(sprintf("ess(V.E)=%.2f ess(V.G.A)=%.2f",
                  ess[2], ess[1]),
          stdout())

    ac <- autocorr.diag(fit$Sol[,1:3])
    write(sprintf("autocorr %s: c_1=%.3f c_2=%.3f c_3=%.3f",
                  rownames(ac)[2], ac[2,1], ac[2,2], ac[2,3]),
          stdout())
    ac <- autocorr.diag(fit$VCV)
    write(sprintf("autocorr %s: V.E=%.3f V.G.A=%.3f",
                  rownames(ac)[2], ac[2,"units"], ac[2, "geno.add"]),
          stdout())

    message("fixed effects:")
    ## print(str(fit$Sol))
    hpdi <- coda::HPDinterval(fit$Sol, prob=0.95)
    fmt <- paste0("c=[%.2f %.2f %.2f]",
                 "\npost.mean=[%.2f %.2f %.2f]",
                 "\npost.sd=[%.2f %.2f %.2f]",
                 "\nHPDi_95(post{c_1})=[%.2f %.2f]",
                 "\nHPDi_95(post{c_2})=[%.2f %.2f]",
                 "\nHPDi_95(post{c_3})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$C[1,1], model$C[2,1], model$C[3,1],
                  mean(fit$Sol[,1]), mean(fit$Sol[,2]), mean(fit$Sol[,3]),
                  sd(fit$Sol[,1]), sd(fit$Sol[,2]), sd(fit$Sol[,3]),
                  hpdi["(Intercept)", "lower"], hpdi["(Intercept)", "upper"],
                  hpdi["year2011", "lower"], hpdi["year2011", "upper"],
                  hpdi["year2012", "lower"], hpdi["year2012", "upper"]),
          stdout())

    message("random effects:")
    ## print(str(fit$VCV))
    hpdi <- coda::HPDinterval(fit$VCV, prob=0.95)
    fmt <- paste0("V.E=%.2f post.mean=%.2f post.sd=%.2f",
                  "\nHPDi_95(post{V.E})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.E,
                  mean(fit$VCV[,"units"]),
                  sd(fit$VCV[,"units"]),
                  hpdi["units", "lower"], hpdi["units", "upper"]),
          stdout())
    fmt <- paste0("V.G.A=%.2f post.mean=%.2f post.sd=%.2f",
                  "\nHPDi_95(post{V.G.A})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.G.A,
                  mean(fit$VCV[,"geno.add"]),
                  sd(fit$VCV[,"geno.add"]),
                  hpdi["geno.add", "lower"], hpdi["geno.add", "upper"]),
          stdout())

    write(sprintf("Cor(G.A,post.mean): S=%.2f P=%.2f",
                  cor(model$G.A[,1], colMeans(fit$Sol[,grepl("ind",colnames(fit$Sol))]),
                      method="spearman"),
                  cor(model$G.A[,1], colMeans(fit$Sol[,grepl("ind",colnames(fit$Sol))]),
                      method="spearman")),
          stdout())
  }

  if(prog.opts$pkg == "rjags"){
    ## print(summary(fit))

    message("convergence diagnostics:")
    write(sprintf("niters=%o", nrow(fit[[1]])), stdout())

    ess <- effectiveSize(fit)
    write(sprintf("ess(c[1])=%.2f ess(c[2])=%.2f ess(c[3])=%.2f",
                  ess["c[1]"], ess["c[2]"], ess["c[3]"]),
          stdout())
    write(sprintf("ess(V.E)=%.2f ess(V.G.A)=%.2f",
                  ess["V.E"], ess["sigma.A2"]),
          stdout())

    ac <- autocorr.diag(fit)
    write(sprintf("autocorr %s: c_1=%.3f c_2=%.3f c_3=%.3f",
                  rownames(ac)[2], ac[2,"c[1]"], ac[2,"c[2]"], ac[2,"c[3]"]),
          stdout())
    write(sprintf("autocorr %s: V.E=%.3f V.G.A=%.3f",
                  rownames(ac)[2], ac[2,"V.E"], ac[2, "sigma.A2"]),
          stdout())

    message("fixed effects:")
    hpdi <- coda::HPDinterval(fit[[1]], prob=0.95)
    fmt <- paste0("c=[%.2f %.2f %.2f]",
                  "\npost.mean=[%.2f %.2f %.2f]",
                  "\npost.sd=[%.3f %.3f %.3f]",
                  "\nHPDi_95(post{c_1})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_2})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_3})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$C[1,1], model$C[2,1], model$C[3,1],
                  mean(fit[[1]][,"c[1]"]),
                  mean(fit[[1]][,"c[2]"]),
                  mean(fit[[1]][,"c[3]"]),
                  sd(fit[[1]][,"c[1]"]),
                  sd(fit[[1]][,"c[2]"]),
                  sd(fit[[1]][,"c[3]"]),
                  hpdi["c[1]", "lower"], hpdi["c[1]", "upper"],
                  hpdi["c[2]", "lower"], hpdi["c[2]", "upper"],
                  hpdi["c[3]", "lower"], hpdi["c[3]", "upper"]),
          stdout())

    message("random effects:")
    fmt <- paste0("V.E=%.2f post.mean=%.2f post.sd=%.3f",
                  "\nHPDi_95(post{V.E})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.E,
                  mean(fit[[1]][,"V.E"]),
                  sd(fit[[1]][,"V.E"]),
                  hpdi["V.E", "lower"], hpdi["V.E", "upper"]),
          stdout())

    fmt <- paste("V.G.A=%.2f post.mean=%.2f post.sd=%.3f",
                 "\nHPDi_95(post{V.G.A})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$V.G.A,
                  mean(fit[[1]][,"sigma.A2"]),
                  sd(fit[[1]][,"sigma.A2"]),
                  hpdi["sigma.A2", "lower"], hpdi["sigma.A2", "upper"]),
          stdout())

    idx <- grep("g.A\\[", colnames(fit[[1]]))
    write(sprintf("Cor(G.A,post.mean): S=%.2f P=%.2f",
                  cor(model$G.A[,1], colMeans(fit[[1]][,idx]),
                      method="spearman"),
                  cor(model$G.A[,1], colMeans(fit[[1]][,idx]),
                      method="pearson")),
          stdout())
  }

  if(prog.opts$pkg == "rrBLUP"){
    message("fixed effects:")
    write(sprintf("c=[%.2f %.2f %.2f] hat{c}=[%.2f %.2f %.2f]",
                  model$C[1,1], model$C[2,1], model$C[3,1],
                  fit$beta[1], fit$beta[2], fit$beta[3]),
          stdout())

    message("random effects:")
    write(sprintf("V.E=%.2f hat{V.E}^2=%.2f",
                  model$V.E,
                  fit$Ve),
          stdout())
    if(use.markers.simul){
      true.bv <- (genomes$genos %*% model$betat + model$u)[,1]
      if(prog.opts$use.markers.infer){
        write(sprintf("sigma_betat^2=%.4f hat{sigma}_betat^2=%.4f",
                      model$sigma.betat2,
                      fit$Vu),
              stdout())
        write(sprintf("Cor(betat,hat{betat}): S=%.2f P=%.2f",
                      cor(model$betat, fit$u, method="spearman"),
                      cor(model$betat, fit$u, method="pearson")),
              stdout())
        infer.bv <- genomes$genos %*% fit$u
        write(sprintf("Cor(X bt + u, X bt): S=%.2f P=%.2f",
                      cor(true.bv, infer.bv, method="spearman"),
                      cor(true.bv, infer.bv, method="pearson")),
              stdout())
        write(sprintf("CV-RMSE(X bt + u)=%.3f",
                      rmse(infer.bv - true.bv, mean(model$y))),
              stdout())
      } else{
        write(sprintf("Var(X bt + u)=%.2f hat{sigma}_u^2=%.2f",
                      var(true.bv),
                      fit$Vu),
              stdout())
        write(sprintf("Cor(X bt + u, hat{u}): S=%.2f P=%.2f",
                      cor(true.bv, fit$u, method="spearman"),
                      cor(true.bv, fit$u, method="pearson")),
              stdout())
        write(sprintf("CV-RMSE(X bt + u)=%.3f",
                      rmse(fit$u - true.bv, mean(model$Y[,1]))),
              stdout())
      }
    } else{
      write(sprintf("V.G.A=%.2f hat{V.G.A}=%.2f",
                    model$V.G.A,
                    fit$Vu),
            stdout())
      write(sprintf("Cor(G.A,hat{G.A}): S=%.2f P=%.2f",
                    cor(model$G.A[,1], fit$u, method="spearman"),
                    cor(model$G.A[,1], fit$u, method="pearson")),
            stdout())
    }
  }

  if(prog.opts$pkg == "rstan"){
    ## rstan::monitor() requires >1 chain, thus use the coda package instead
    tmp <- rstan::As.mcmc.list(fit)

    message("convergence diagnostics:")
    write(sprintf("niters=%i", nrow(tmp[[1]])), stdout())

    ess <- effectiveSize(tmp)
    write(sprintf("ess(c[1])=%.2f ess(c[2])=%.2f ess(c[3])=%.2f",
                  ess["c[1]"], ess["c[2]"], ess["c[3]"]),
          stdout())
    write(sprintf("ess(sigma.E)=%.2f ess(sigma.g.A)=%.2f",
                  ess["sigma_E"], ess["sigma_g_A"]),
          stdout())

    ac <- autocorr.diag(tmp)
    write(sprintf("autocorr %s: c_1=%.3f c_2=%.3f c_3=%.3f",
                  rownames(ac)[2], ac[2,"c[1]"], ac[2,"c[2]"], ac[2,"c[3]"]),
          stdout())
    write(sprintf("autocorr %s: sigma.E=%.3f sigma.g.A=%.3f",
                  rownames(ac)[2], ac[2,"sigma_E"], ac[2, "sigma_g_A"]),
          stdout())

    message("fixed effects:")
    hpdi <- coda::HPDinterval(tmp[[1]], prob=0.95)
    fmt <- paste0("c=[%.2f %.2f %.2f]",
                  "\npost.mean=[%.2f %.2f %.2f]",
                  "\npost.sd=[%.3f %.3f %.3f]",
                  "\nHPDi_95(post{c_1})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_2})=[%.2f %.2f]",
                  "\nHPDi_95(post{c_3})=[%.2f %.2f]")
    write(sprintf(fmt,
                  model$C[1,1], model$C[2,1], model$C[3,1],
                  mean(tmp[[1]][,"c[1]"]),
                  mean(tmp[[1]][,"c[2]"]),
                  mean(tmp[[1]][,"c[3]"]),
                  sd(tmp[[1]][,"c[1]"]),
                  sd(tmp[[1]][,"c[2]"]),
                  sd(tmp[[1]][,"c[3]"]),
                  hpdi["c[1]", "lower"], hpdi["c[1]", "upper"],
                  hpdi["c[2]", "lower"], hpdi["c[2]", "upper"],
                  hpdi["c[3]", "lower"], hpdi["c[3]", "upper"]),
          stdout())

    message("random effects:")
    fmt <- paste0("sigma_E=%.2f post.mean=%.2f post.sd=%.2f",
                  "\nHPDi_95(post{sigma_E})=[%.2f %.2f]")
    write(sprintf(fmt,
                  sqrt(model$V.E),
                  ## mean(extract(fit, "sigma_E")[[1]]^2)),
                  mean(tmp[[1]][,"sigma_E"]),
                  sd(tmp[[1]][,"sigma_g_A"]),
                  hpdi["sigma_E", "lower"], hpdi["sigma_E", "upper"]),
          stdout())
    fmt <- paste0("sigma.g.A=%.2f post.mean=%.2f post.sd=%.2f",
                  "\nHPDi_95(post{sigma.G.A})=[%.2f %.2f]")
    write(sprintf(fmt,
                  sqrt(model$V.G.A),
                  ## mean(extract(fit, "sigma_g_A")[[1]]^2)),
                  mean(tmp[[1]][,"sigma_g_A"]),
                  sd(tmp[[1]][,"sigma_g_A"]),
                  hpdi["sigma_g_A", "lower"], hpdi["sigma_g_A", "upper"]),
          stdout())
    idx <- grep("g_A\\[", colnames(tmp[[1]]))
    write(sprintf("Cor(G.A,post.mean): S=%.2f P=%.2f",
                  cor(model$G.A[,1],
                      colMeans(tmp[[1]][,idx]),
                      method="spearman"),
                  cor(model$G.A[,1],
                      colMeans(tmp[[1]][,idx]),
                      method="pearson")),
          stdout())
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
##' @author Timothee Flutre
quantgenEvalMain <- function(prog.args, simul.file, infer.file, eval.file){
  prog.opts <- list(verbose=1,
                    pkg=NULL,
                    use.markers.infer=FALSE,
                    prior.mark=NULL)

  prog.opts <- quantgenEvalParseArgs(prog.args, prog.opts)

  quantgenEvalCheckOptions(prog.opts)

  quantgenEvalRun(prog.opts, simul.file, infer.file, eval.file)
}
