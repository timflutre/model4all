## Aim: simulate data for the "quantitative genetics" theme
## Copyright (C) 2015-2016 INRA
## License: GPL-3+
## Persons: Timothee Flutre [cre,aut]
## Versioning: https://mulcyber.toulouse.inra.fr/projects/comp-fit-gmrf/

quantgen.simul.name <- "quantgen_simul"
quantgen.simul.version <- "0.4.1" # http://semver.org/

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
quantgenSimulHelp <- function(){
  txt <- paste0("`", quantgen.simul.name, "' simulates data for the \"quantitative genetics\" theme.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --ntraits\tnumber of traits\n")
  txt <- paste0(txt, "\t\tdefault=1\n")
  txt <- paste0(txt, "      --ngenos\tnumber of genotypes\n")
  txt <- paste0(txt, "\t\tdefault=100\n")
  txt <- paste0(txt, "      --Ne\teffective population size\n")
  txt <- paste0(txt, "\t\tdefault=10^4\n")
  txt <- paste0(txt, "      --cL\tchromosome length in bp (same for all 20 of them)\n")
  txt <- paste0(txt, "\t\tdefault=10^5\n")
  txt <- paste0(txt, "      --mu\tneutral mutation rate in events/base/generation\n")
  txt <- paste0(txt, "\t\tdefault=10^-8\n")
  txt <- paste0(txt, "      --c\tneutral recombination rate in events/base/generation\n")
  txt <- paste0(txt, "\t\tdefault=10^-8\n")
  txt <- paste0(txt, "      --nyears\tnumber of years\n")
  txt <- paste0(txt, "\t\tdefault=3\n")
  txt <- paste0(txt, "      --VGA\tvariance of the breeding values\n")
  txt <- paste0(txt, "\t\tdefault=15; requires --ntraits 1\n")
  txt <- paste0(txt, "\t\tset --VGA NULL and use --nuWA when --ntraits > 1\n")
  txt <- paste0(txt, "      --nuWA\tdegrees of freedom of the Wishart prior for V_{G_A}\n")
  txt <- paste0(txt, "\t\tdefault=5; used only if --ntraits > 1\n")
  txt <- paste0(txt, "      --dom\tuse dominant relationships\n")
  txt <- paste0(txt, "\t\tby default, only additive relationships are used\n")
  txt <- paste0(txt, "      --VGD\tvariance of the dominant deviations\n")
  txt <- paste0(txt, "\t\tdefault=3 if --dom; requires --ntraits 1\n")
  txt <- paste0(txt, "\t\tset --VGD NULL and use --nuWD when --ntraits > 1\n")
  txt <- paste0(txt, "      --nuWD\tdegrees of freedom of the Wishart prior for V_{G_D}\n")
  txt <- paste0(txt, "\t\tdefault=5; used only if --ntraits > 1\n")
  txt <- paste0(txt, "      --VE\tvariance of the errors\n")
  txt <- paste0(txt, "\t\tdefault=5; requires --ntraits 1; ignored if --errdf not Inf\n")
  txt <- paste0(txt, "\t\tset --VE NULL and use --nuWE when --ntraits > 1\n")
  txt <- paste0(txt, "      --errdf\tdegrees of freedom of the errors' Student's t-distribution\n")
  txt <- paste0(txt, "\t\tdefault=Inf (Normal); e.g. 3\n")
  txt <- paste0(txt, "\t\tignored if --ntraits > 1\n")
  txt <- paste0(txt, "      --nuWE\tdegrees of freedom of the Wishart prior for V_E\n")
  txt <- paste0(txt, "\t\tdefault=5; used only if --ntraits > 1\n")
  txt <- paste0(txt, "      --na\tpercentage of missing phenotypes, at random (default=0)\n")
  txt <- paste0(txt, "      --mark\texplicitly use markers' genotypes when simulating phenotypes\n")
  txt <- paste0(txt, "\t\totherwise (i.e. by default), genetic relationships matrices are used\n")
  txt <- paste0(txt, "\t\tfor the moment, only works with --ntraits 1\n")
  txt <- paste0(txt, "      --pi\tproportion of SNP effects that are non-zero\n")
  txt <- paste0(txt, "\t\tonly used with --mark\n")
  txt <- paste0(txt, "\t\tdefault: log(pi) ~ Unif(log(1/P), log(1))\n")
  txt <- paste0(txt, "      --pveA\tproportion of phenotypic variance explained by SNPs\n")
  txt <- paste0(txt, "\t\twith non-zero effect\n")
  txt <- paste0(txt, "\t\trequires --VGA not NULL to deduce VE\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <timothee.flutre@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
quantgenSimulVersion <- function(){
  txt <- paste0(quantgen.simul.name, " ", quantgen.simul.version, "\n")
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
quantgenSimulParseArgs <- function(prog.args, prog.opts){
  i <- 0
  while(i < length(prog.args)){
    i <- i + 1
    if(prog.args[i] == "-h" || prog.args[i] == "--help"){
      quantgenSimulHelp()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-V" || prog.args[i] == "--version"){
      quantgenSimulVersion()
      quit("no", status=0)
    }
    else if(prog.args[i] == "-v" || prog.args[i] == "--verbose"){
      prog.opts$verbose <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--ntraits"){
      prog.opts$nb.traits <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--ngenos"){
      prog.opts$nb.genos <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--Ne"){
      prog.opts$Ne <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--cL"){
      prog.opts$chrom.len <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--mu"){
      prog.opts$mu <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--c"){
      prog.opts$c <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--nyears"){
      prog.opts$nb.years <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--VGA"){
      if(prog.args[i+1] == "NULL"){
        prog.opts$V.G.A <- NULL
      } else
        prog.opts$V.G.A <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--nuWA"){
      prog.opts$nu.G.A <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--dom"){
      prog.opts$use.dominant <- TRUE
    }
    else if(prog.args[i] == "--VGD"){
      if(prog.args[i+1] == "NULL"){
        prog.opts$V.G.D <- NULL
      } else
        prog.opts$V.G.D <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--nuWD"){
      prog.opts$nu.G.D <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--VE"){
      if(prog.args[i+1] == "NULL"){
        prog.opts$V.E <- NULL
      } else
        prog.opts$V.E <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--errdf"){
      prog.opts$err.df <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--nuWE"){
      prog.opts$nu.E <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--na"){
      prog.opts$perc.NA <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--mark"){
      prog.opts$use.markers <- TRUE
    }
    else if(prog.args[i] == "--pi"){
      prog.opts$pi <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else if(prog.args[i] == "--pveA"){
      prog.opts$pve.A <- as.numeric(prog.args[i+1])
      i <- i + 1
    }
    else{
      write(paste0(quantgen.simul.name, ": invalid option -- ", prog.args[i], "\n"), stderr())
      quantgenSimulHelp()
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
##' @author Timothee Flutre
quantgenSimulCheckOptions <- function(prog.opts){
  suppressPackageStartupMessages(library(rutilstimflutre)) # https://github.com/timflutre/rutilstimflutre
  if(packageVersion("rutilstimflutre") < "0.60.0"){
    msg <- "ERROR: version of 'rutilstimflutre' should be >= 0.60.0"
    write(msg, stderr())
    quit(save="no", status=1)
  }
  if(prog.opts$use.markers & prog.opts$nb.traits > 1){
    msg <- "ERROR: --mark should go with --ntraits 1"
    write(msg, stderr())
    quit(save="no", status=1)
  }
  return(prog.opts)
}

##' Perform the simulation(s)
##'
##'
##' @param prog.opts named list
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Timothee Flutre
quantgenSimulRun <- function(prog.opts, simul.file, seed){
  set.seed(seed)

  if(prog.opts$verbose > 0){
    msg <- paste0("pkg 'rutilstimflutre' v", packageVersion("rutilstimflutre"))
    write(msg, stdout())
  }

  genomes <- NULL
  A <- NULL
  model <- NULL

  if(prog.opts$verbose > 0){
    msg <- paste0("simulate ", prog.opts$nb.genos, " genomes ...")
    write(msg, stdout())
  }
  genomes <- simulCoalescent(nb.inds=prog.opts$nb.genos,
                                 pop.mut.rate=4 * prog.opts$Ne *
                                   prog.opts$mu * prog.opts$chrom.len,
                                 pop.recomb.rate=4 * prog.opts$Ne *
                                   prog.opts$c * prog.opts$chrom.len,
                                 chrom.len=prog.opts$chrom.len,
                                 verbose=prog.opts$verbose - 1)
  if(prog.opts$verbose > 0)
    write(paste0(ncol(genomes$genos), " markers generated ..."), stdout())

  if(prog.opts$use.markers){
    if(prog.opts$verbose > 0)
      write("simulate phenotypes (use markers explicitly) ...", stdout())

    model <- simulBvsr(Q=prog.opts$nb.years,
                       X=genomes$genos,
                       pi=prog.opts$pi,
                       pve.A=prog.opts$pve.A,
                       sigma.a2=prog.opts$V.G.A,
                       err.df=prog.opts$err.df,
                       perc.NA=prog.opts$perc.NA)
    A <- model$A
  } else{
    if(prog.opts$verbose > 0)
      write("estimate genetic relationships ...", stdout())
    A <- estimGenRel(X=genomes$genos,
                     relationships="additive", method="vanraden1",
                     verbose=prog.opts$verbose - 1)
    D <- NULL
    if(prog.opts$use.dominant)
      D <- estimGenRel(X=genomes$genos,
                       relationships="dominant", method="vitezica",
                       verbose=prog.opts$verbose - 1)

    if(prog.opts$verbose > 0)
      write(paste0("simulate phenotypes (use relationships, ",
                   prog.opts$nb.traits, " trait",
                   ifelse(prog.opts$nb.traits > 1, "s", ""), ") ..."),
            stdout())
    model <- simulAnimalModel(T=prog.opts$nb.traits,
                              Q=prog.opts$nb.years,
                              A=A,
                              V.G.A=prog.opts$V.G.A,
                              nu.G.A=prog.opts$nu.G.A,
                              D=D,
                              V.G.D=prog.opts$V.G.D,
                              V.E=prog.opts$V.E,
                              err.df=prog.opts$err.df,
                              perc.NA=prog.opts$perc.NA)
  }

  if(prog.opts$verbose > 0)
    write("save into file ...", stdout())
  save(seed, genomes, A, D, model, file=simul.file)
}

##' Entry point of the program
##'
##'
##' @param prog.args character
##' @param simul.file character
##' @param seed integer
##' @return nothing
##' @author Timothee Flutre
quantgenSimulMain <- function(prog.args, simul.file, seed){
  prog.opts <- list(verbose=1,
                    nb.traits=1,
                    nb.genos=100,
                    Ne=10^4,
                    chrom.len=10^5,
                    mu=10^(-8),
                    c=10^(-8),
                    nb.years=3,
                    scale.halfCauchy=NULL,
                    V.G.A=15,
                    nu.G.A=5,
                    use.dominant=FALSE,
                    V.G.D=3,
                    nu.G.D=5,
                    V.E=5,
                    err.df=Inf,
                    nu.E=5,
                    perc.NA=0,
                    use.markers=FALSE,
                    pi=NULL,
                    pve.A=NULL)

  prog.opts <- quantgenSimulParseArgs(prog.args, prog.opts)

  prog.opts <- quantgenSimulCheckOptions(prog.opts)

  quantgenSimulRun(prog.opts, simul.file, seed)
}
