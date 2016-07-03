
## give the lower and upper of region of uncertainty and the output is the lower and upper for each parameter in a list
make_vertices <- function(lower, upper){
  par_list <- vector("list",length(lower))
  for(i in 1:length(lower))
    par_list[[i]] <- c(lower[i], upper[i])

  vertices <- matrix(unlist(expand.grid(par_list)), nrow = 2^length(upper))
  return(vertices)
}




