install_if_not_installed <- function(pkg, install_function = install.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (identical(install_function, install.packages)) {
      install.packages(pkg)
    } else {
      install_function(pkg)
    }
  }
}


library(VariantAnnotation)
install.packages(
  'gwasglue',
  repos = c(
    'https://mrcieu.r-universe.dev',
    'https://cloud.r-project.org'
  )
)

BiocManager::version()
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/Bioconductor")
BiocManager::install("VariantAnnotation")
BiocManager::install('VariantAnnotation')
BiocManager::install("VariantAnnotation")

devtools::install_github("Bioconductor/VariantAnnotation")

BiocManager::install("VariantAnnotation", version = "3.19")









