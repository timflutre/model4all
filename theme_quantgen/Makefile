README_theme-quantgen.pdf: README_theme-quantgen.Rmd
	echo 'library(rmarkdown); render("README_theme-quantgen.Rmd", "pdf_document")' | R --vanilla --quiet

README_theme-quantgen.html: README_theme-quantgen.Rmd
	echo 'library(rmarkdown); render("README_theme-quantgen.Rmd", "html_document")' | R --vanilla --quiet

clean:
	rm -f README_theme-quantgen.md \
		README_theme-quantgen.html \
		README_theme-quantgen.pdf \
		*~
	rm -rf README_theme-quantgen_cache/
	rm -rf README_theme-quantgen_files/
