#' Check for missing packages
#'
#' @noRd
check_package <- function(package, repo = "CRAN", git_repo = ""){
    
    if (repo == "CRAN"){
        install_function <- "install.packages('"
    } else if (repo == "github") {
        install_function <- paste0("remotes::install_github('", git_repo, "/")
    } else if (repo == "bioc"){
        install_function <- "BiocManager::install('"
    } else{
        install_function <- "Unknown repository.. "
    }
    
    if(!requireNamespace(package, quietly = TRUE)){
        stop(
            "Package ", package," is not installed.\n",
            "Please run: ", install_function, package, "')")
    }
    requireNamespace(package, quietly = TRUE)
}

set_apply <- function(mc.cores, pb = TRUE) {
    if(mc.cores == 1) {
        APPLY <- lapply
    } else {
        if(pb) {
            cyCombine:::missing_package(package = "pbmcapply")
            APPLY <- pbmcapply::pbmclapply
        } else {
            APPLY <- parallel::mclapply
        }
        formals(APPLY)$mc.cores <- mc.cores
    }
    return(APPLY)
}