.PHONY: document test check install format lint site clean

document:
	Rscript -e "devtools::document()"

test:
	Rscript -e "devtools::test()"

check:
	Rscript -e "devtools::check()"

install:
	Rscript -e "devtools::install()"

format:
	air format .

lint:
	Rscript -e "lintr::lint_package()"

site:
	Rscript -e "altdoc::render_docs()"

site-preview:
	Rscript -e "altdoc::preview_docs()"

clean:
	rm -rf docs/ man/ causatr.Rcheck