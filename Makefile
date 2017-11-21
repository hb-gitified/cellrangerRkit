all:	rkit vignette

rkit:
	mkdir -p build
	R CMD build --no-build-vignettes .

# NOTE: You have to run this TWICE to update the table of contents.
vignette:
	Rscript scripts/run_knitr.R

clean:
	rm -rf build
