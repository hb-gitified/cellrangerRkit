all:	rkit vignette

rkit:
	mkdir -p build
	R CMD build --no-build-vignettes .

vignette:
	Rscript scripts/run_knitr.R

clean:
	rm -rf build
