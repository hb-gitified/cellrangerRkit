#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

# compile document with knitr
library('knitr')
setwd(dir='vignettes')

cat('Running Sweave2knitr...\n')
Sweave2knitr("cellrangerrkit-PBMC-vignette.Rnw")

cat('Running knit...\n')
knit(file.path(getwd(), 'cellrangerrkit-PBMC-vignette-knitr.Rnw'))

cat ('Running pdflatex...\n')
system("pdflatex cellrangerrkit-PBMC-vignette-knitr.tex")
