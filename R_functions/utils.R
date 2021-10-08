
check_directory <- function(dir_name, recursive=TRUE){
  ## checks if directory exists, if not, creates it
  if (!dir.exists(dir_name)){
    dir.create(dir_name, recursive=recursive)
    warning(paste0("The following given directory did not exist and was created by 'check_directory': ", dir_name))
  }
}

safe_save_rds <- function(object, file_path, recursive=TRUE){
  ## object character -> ___
  ## save object to file_path. If dirname(file_path) does not exists, create it
  
  dir_name <- dirname(file_path)
  check_directory(dir_name, recursive=recursive)
  
  saveRDS(object, file = file_path)
  
}

last_elem <- function(x){
  # Fastest way to get the last element of a vector x in O(1)
  # (better than Rcpp::mylast(x), tail(x, n=1), dplyr::last(x), x[end(x)[1]]], and rev(x)[1])
  x[length(x)]
}

