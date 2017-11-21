#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
check_package_version = function(pkg_name, min_version) {
    package_version = packageVersion(pkg_name)
    if(package_version < min_version) stop(sprintf("Version %s of %s found, need at least %s", package_version, pkg_name, min_version))
}

message('\nCell Ranger R Kit Installer')
message('Copyright (c) 2016 10x Genomics, Inc. All rights reserved.')
message('----------------------------------------------------------')
eula <- readline("By continuing, you are agreeing to the EULA at: http://support.10xgenomics.com/license\n[ hit enter to continue ]\n")

# Try to install rhdf5 correctly for all versions of R
source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates=TRUE)

biocLite("rhdf5", siteRepos="http://bioconductor.org/packages/3.3/bioc", suppressUpdates=TRUE)
check_package_version("rhdf5", 2.14)

packages = c(
  "Matrix",
  "ggplot2",
  "bit64",
  "data.table",
  "Rtsne",
  "pheatmap",
  "irlba",
  "Rmisc"
)

for(package in packages) {
  install.packages(package)
}

install.packages("http://s3-us-west-2.amazonaws.com/10x.files/code/cellrangerRkit-2.0.0.tar.gz", repos=NULL)
