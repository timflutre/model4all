Project "model4all"
===================

This directory contains statistical models and computer code as part of a **collaborative, pedagogical initiative for statistical modelling based on simulation, between statisticians, and plant geneticists and eco-physiologists**.

The coordinators are [Timothée Flutre](http://openwetware.org/wiki/User:Timothee_Flutre) ([INRA](https://en.wikipedia.org/wiki/Institut_national_de_la_recherche_agronomique)) and [Marie Denis](https://www.researchgate.net/profile/Marie_Denis2) ([CIRAD](https://en.wikipedia.org/wiki/Centre_de_coop%C3%A9ration_internationale_en_recherche_agronomique_pour_le_d%C3%A9veloppement)), both working at [AGAP](http://umr-agap.cirad.fr/en) (a joint research unit in Montpellier, France).
The project hence is currently funded, and the copyright owned, by INRA and CIRAD.
However, the goal is to actively **promote collaboration "in the open"**, a much-needed endeavour in face of global challenges.
Consequently, all contributions in terms of documents are under the [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/) license, and all contributions in terms of code are under the [GNU AGPL v3](https://www.gnu.org/licenses/agpl.html).
Note however that using system calls to assess external, potentially private programs, is possible.

The content of this directory is versioned using [git](https://en.wikipedia.org/wiki/Git_(software)), the central repository being hosted [**here**](https://github.com/timflutre/model4all), on [GitHub](https://en.wikipedia.org/wiki/GitHub).
If you are not familiar with these, please, do read the official [git book](https://www.git-scm.com/book/en/v2) (at the very least the first three chapters), and don't restrain yourself from reading [GitHub's help](https://help.github.com/).


Overview
--------

A **theme** gathers a set of simulation and inferential models which are used to explore a common research question.
For instance, the theme `quantgen` is concerned with quantitative genetics, the theme `corrobs` deals with correlated observations, etc.
See the various `theme_<...>` sub-directories.

The general approach is decomposed into **three steps**:

1. **simulate** data according to a given model;

2. **infer** with one or several statistical models (potentially different from those used at step 1), via one or several inferential methods, implemented in one or several softwares;

3. **evaluate** model fit, parameters' estimation, prediction accuracy, computation time, memory requirements, etc, for each inference performed at step 2.

Each step is performed by a **generic program**: `model4all_simul`, `model4all_infer` and `model4all_eval`.
Each generic program then calls or executes a piece of code specific to the theme of interest.

More **details** can be found in `doc/README_model4all.html`.

A presentation describing the origin and motivation of this initiative is also [available](http://prodinra.inra.fr/?locale=en#!ConsultNotice:359165) in french.


Installation
------------

Currently, the generic part of the project is based on the [R](https://en.wikipedia.org/wiki/R_(programming_language)) software environment (but any external program can be executed via a system call).
At least version 3.0 is necessary, along with the [rmarkdown](https://cran.r-project.org/package=rmarkdown) package.
Depending on the theme, various other software components may be required, such as R packages, external programs, etc.
Please refer to each particular theme.
Note that as the total number of tools can be large, it is judicious to use `model4all` on a shared computer (e.g. server, cluster) so that tools are installed once and available to all users.

On a computer with GNU/Linux or Mac OS, the content of the git repository can be downloaded and installed with the following command-lines executed in a terminal:

```
wget -O model4all-master https://github.com/timflutre/model4all/archive/master.zip
unzip model4all-master.zip
cd model4all-master
make
make check
make install
```

By default, the three generic programs (see above) are installed in a directory `bin` in `${HOME}`.
If the `bin` directory does not exist, it will be created.
To install elsewhere, use the option `INSTALL`:

```
sudo make install INSTALL="/usr/local"
```

To check if the generic programs are properly installed and available in your PATH, execute the following command-line:

```
which model4all_simul
model4all_simul -h
```

The **main entry point** to the project is the file `doc/README_model4all.Rmd`.
Using the `rmarkdown` package, it can be converted into `html` or `pdf` via the following command-lines:

```
make html
make pdf
```

To convert into `pdf`, note that a valid [LaTeX](https://en.wikipedia.org/wiki/LaTeX) distribution, such as [TeX Live](https://en.wikipedia.org/wiki/TeX_Live), should be available on your computer.


Contribution
------------

There are so many ways to improve this initiative that it's hard to know where to start.
Nevertheless, here is a try:
- fix a typo in the documentation;
- fix a programming error in the computer code;
- improve an existing theme (new simulation, inference method, software, evaluation criterion, etc);
- add a whole new theme;
- etc

Concretely, you can contribute by:
- reporting [issues](https://github.com/timflutre/model4all/issues) online;
- proposing [pull requests](https://github.com/timflutre/model4all/pulls) (to the "dev" branch, not "master"; read more on [branching](https://www.git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell)).
