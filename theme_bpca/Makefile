README_theme-bpca.pdf: README_theme-bpca.Rmd
	echo 'library(rmarkdown); render("README_theme-bpca.Rmd", "pdf_document")' | R --vanilla --quiet

README_theme-bpca.html: README_theme-bpca.Rmd
	echo 'library(rmarkdown); render("README_theme-bpca.Rmd", "html_document")' | R --vanilla --quiet

clean:
	rm -f README_theme-bpca.md \
		README_theme-bpca.html \
		README_theme-bpca.pdf \
		*~
	rm -rf README_theme-bpca_cache/
	rm -rf README_theme-bpca_files/
