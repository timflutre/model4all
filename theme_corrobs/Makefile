README_theme-corrobs.pdf: README_theme-corrobs.Rmd
	echo 'library(rmarkdown); render("README_theme-corrobs.Rmd", "pdf_document")' | R --vanilla --quiet

README_theme-corrobs.html: README_theme-corrobs.Rmd
	echo 'library(rmarkdown); render("README_theme-corrobs.Rmd", "html_document")' | R --vanilla --quiet

clean:
	rm -f README_theme-corrobs.md \
		README_theme-corrobs.html \
		README_theme-corrobs.pdf \
		*~
	rm -rf README_theme-corrobs_cache/
	rm -rf README_theme-corrobs_files/
