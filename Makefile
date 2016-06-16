INSTALL=${HOME}
VERSION=0.1.0

GENERIC_EXE=src/model4all_simul src/model4all_infer src/model4all_eval

all:
	@echo "no compilation required (yet)"

check:
	@echo "no test implemented (yet)"

install: $(GENERIC_EXE)
	mkdir -p ${INSTALL}/bin
	cp $(GENERIC_EXE) ${INSTALL}/bin

html: doc/README_model4all.Rmd
	cd doc/; echo 'library(rmarkdown); render("README_model4all.Rmd", "html_document")' | R --vanilla --quiet

pdf: doc/README_model4all.Rmd
	cd doc/; echo 'library(rmarkdown); render("README_model4all.Rmd", "pdf_document")' | R --vanilla --quiet

dist:
	rm -rf model4all-${VERSION}
	rm -f model4all-${VERSION}.tar.gz
	mkdir -p model4all-${VERSION}
	cp Makefile README.md model4all-${VERSION}/
	mkdir -p model4all-${VERSION}/doc
	cp doc/README_model4all.Rmd model4all-${VERSION}/doc
	mkdir -p model4all-${VERSION}/src
	cp src/model4all_simul model4all-${VERSION}/src
	cp src/model4all_infer model4all-${VERSION}/src
	cp src/model4all_eval model4all-${VERSION}/src
	for theme in quantgen corrobs bpca; do \
	 	mkdir -p model4all-${VERSION}/theme_$$theme ; \
	 	cp theme_$$theme/README_theme-$$theme.Rmd model4all-${VERSION}/theme_$$theme; \
	 	cp theme_$$theme/Makefile model4all-${VERSION}/theme_$$theme; \
	 	cp theme_$$theme/$${theme}_simul.R model4all-${VERSION}/theme_$$theme; \
	 	cp theme_$$theme/$${theme}_infer.R model4all-${VERSION}/theme_$$theme; \
	 	cp theme_$$theme/$${theme}_eval.R model4all-${VERSION}/theme_$$theme; \
	done
	tar -czf model4all-${VERSION}.tar.gz model4all-${VERSION}/
	rm -rf model4all-${VERSION}

uninstall:
	rm ${INSTALL}/bin/$(GENERIC_EXE)
