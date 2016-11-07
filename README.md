# 10X Genomics Cell Ranger R kit

## Description
R package that accompanies the 10X Genomics Cell Ranger pipeline. It supports manipulation and visualization of output from the pipeline as well as re-running of secondary analyses.

## Development
Make sure you have R packages devtools and roxygen2 installed and loaded:

```R
install.packages( c('devtools', 'roxygen2') )
library(devtools)
library(roxygen2)
```

To install the current version from master:

```R
install_github( '10XDev/cellrangerRkit', user = 'github_user',
    auth_token = 'some_auth_token' )
```

The auth token can be generated using the default settings in
https://github.com/settings/applications under 'Personal access tokens'. Make
sure to save the token, otherwise you need to generate another one.

## Workflow
Since you've already installed devtools, you should consider using it to help
develop the package. Note that all `'path/to/working/copy'` strings can be
replaced by the empty string if you're inside the root of the package directory
(which I strongly recommend). There is a developer mode you can enter by
typing:

```R
dev_mode()
```
and exit similarly. What this does is create a sandboxed development
environment for your packages.

Throughout development, you probably will need to load updates regularly:

```R
install('path/to/working/copy')
```

If you're trying to access a new function outside of the package (in
userspace land) make sure to add the `@export` decorator and generate
documentation for it using roxygen2:

```R
document('path/to/working/copy')
```

Once you dig the changes and you think it's stable, drop out of developer mode
and install it into your main package set.

### Testing

We are using the testthat package for unit testing. All tests live in
`inst/tests/`. The directory `tests/` simply contains a stub to test all things
that live in `inst/tests/`. Using devtools you can call all tests by running:

```R
test('path/to/working/copy')
```

### Building for distribution

devtools::document()
devtools::test()
devtools::build()
devtools::build_vignettes()


## Uninstalling
```R
remove.packages('cellrangerRkit')
```

## More reading
* Some basics on making/installing packages:
  http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
  * More details on philosophy and workflow: http://r-pkgs.had.co.nz/