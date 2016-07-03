#'  Neighborhood based crowding Differential Evolution (NCDE)  for multimodal optimization
#'
#'
#'  A niching differential evolution (DE) algorithm with neighborhood mutation strategy to solve multimodal optimization problems.
#'
#'
#'
#'
#' @param fn objective function that should be maximized. Should not return \code{NaN}.
#' @param lower lower bound.
#' @param upper upper bound.
#' @param control control parameters for the algorithm. See "Details".
#' @param ... extra arguments are passed to \code{fn}.
#'
#' @references
#' Qu, B. Y., Suganthan, P. N., & Liang, J. J. (2012). Differential evolution with neighborhood mutation for multimodal optimization. IEEE transactions on evolutionary computation, 16(5), 601-614.\cr
#'  Based on a MATLAB code that can be found in Suganthan's home page.
#'
#' @note
#' The function is maximized.\cr
#'
#'
#' Check whether \code{fn} does not return any \code{NaN}.\cr
#'
#' The only stopping rule is the number of iterations.\cr
#'
#' The implemented algorithm here slightly differs from the original
#'  MATLAB code in a way that the off-springs for all \code{NP} members calculated by vectorization (\code{sapply})
#'   and not with \code{for} loop as the original MATLAB code.
#' @details
#'
#' @details
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{iter}}{number of iterations. Defaults to \code{200}.}
#'   \item{\code{SF}}{F is a scale factor in \eqn{[0, 2]} used for scaling the differential vector. Defaults to \code{0.9}.}
#'   \item{\code{CR}}{crossover probability from \eqn{[0, 1]}. Defaults to \code{0.2}.}
#'   \item{\code{NP}}{Population size. Defaults to \code{10 * length(lower)}.}

#'   \item{\code{seed}}{random seed.}
#'   \item{\code{hybrid}}{logical; if true, before quiting the algorithm,  an \code{L-BFGS-B}
#'    search with the provided position as initial guess is done to improve the accuracy of the results.
#'      Defaults to \code{TRUE}. Note that no attempt is done to control the maximal number of
#'       function evaluations within the local search step (this can be done separately through \code{hybrid.control})
#'       }
#'    \item{\code{hybrid_control}}{List with any additional control parameters to pass on to \code{\link{stats}{optim}}
#'     when using \code{L-BFGS-B} for the local search. Defaults to \code{NULL}.}
#' }
#'
#' The DE strategy is 'DE/rand/1'.\cr
#'
#' @return
#'  a list contains:
#' \describe{
#'   \item{\code{pop}}{a matrix; position of the population.}
#'   \item{\code{popval}}{a vector; corresponding population values.}
#'   \item{\code{nfeval}}{number of function evaluations.}
#'   \item{\code{maxima}}{position of particles after using \code{L-BFGS-B} when \code{control$hybrid = TRUE}. It should be local maxima. If \code{control$hybrid == FALSE}, then it is \code{NA}.}
#'   \item{\code{maximaval}}{fitness values of \code{maxima}.}
#' }
#'
#'
#' The concept of
#' niching is inspired by the way organisms evolve in nature.
#' When integrated with EAs, niching involves the formation
#' of subpopulations within a population. Each subpopulation
#' aims to locate one optimal solution and together the whole
#' population is expected to locate multiple optimal solutions in a single run (Qu et al., 2012).\cr
#'
#' @examples
#'  ####################################################################################
#'  ## Two-Peak Trap:  global maximum on x = 20 and local maximum on x = 0
#' two_peak <- function(x)
#'    y <- (160/15) * (15 - x) * (x < 15) + 40 * (x - 15) * (x >= 15)
#'
#'  ncde(two_peak, 0, 20, control = list(seed = 66, NP = 50))
#'
#'  # without local search
#'  ncde(two_peak, 0, 20, control = list(seed = 66, NP = 50))
#'
#'  # without refining and local search
#'  ncde(two_peak, 0, 20, control = list(seed = 66, NP = 50,  hybrid = FALSE))
#'
#'  ####################################################################################
#' # Decreasing Maxima: one global on x = 0.1 and four local maxima
#' dmaxima <- function(x)
#'   y <- exp(-2 * log(2) * ((x - 0.1)/0.8)^2) * (sin(5 * pi * x))^6
#'
#' res <- ncde(dmaxima, 0, 1, control = list(seed = 66, NP = 10))
#' unique(round(res$maxima, 5)) ## ferpsols can not find the local maxima
#'
#' ## plot
#' x <- seq(0, 1, length.out = 400)
#' plot(x,  dmaxima(x), type = "l")
#'
#'  ####################################################################################
#' # Himmelblau's function: four global optima on
#' # x = c(-2.80512, 3.13131), c(3.00000, 2.00000), c(-3.77931, -3.28319) and c(3.58443, -1.84813)
#' Himmelblau <- function(x){
#'   y <- - (x[1]^2 + x[2] - 11)^2 - (x[1] + x[2]^2 - 7)^2
#'   return(y)
#' }
#'
#'
#' res <- ncde(Himmelblau, c(-6, -6), c(6, 6), control = list(seed = 66, NP = 50))
#' unique(round(res$maxima, 5))
#'
#'
#' Himmelblau_plot <- function(x, y)
#'   Himmelblau(x = c(x, y))
#' Himmelblau_plot <- Vectorize(Himmelblau_plot)
#' x <- y <- seq(-6, 6, length.out = 100)
#' persp(x, y, z = outer(X = x, Y = y, FUN = Himmelblau_plot))
#'
#'  ####################################################################################
#' # Six-Hump Camel Back: two global and two local maxima
#' Six_Hump <- function(x){
#'   factor1 <- (4 - 2.1 * (x[1]^2) + (x[1]^4)/3) * (x[1]^2) + x[1] * x[2]
#'   factor2 <- (-4 + 4 * (x[2]^2)) * (x[2]^2)
#'   y <- -4 * (factor1 + factor2)
#'   return(y)
#' }
#' res <- ncde(Six_Hump, c(-1.9, -1.1), c(1.9, 1.1), control = list(seed = 66, NP = 10))
#' unique(round(res$maxima, 5))
#'
#'  ####################################################################################
#'  ## 2D Inverted Shubert function :
#'  # The global minima: 18 global minima  f(x*) = -186.7309.
#'  # the local maxima: sevral
#'
#'  Shubert <- function(x){
#'     j <- 1:5
#'     out <- -(sum(j * cos((j + 1) * x[1] + j)) * sum(j * cos((j + 1) * x[2] + j)))
#'     return(out)
#'  }
#'
#'  res <- ncde(Shubert, rep(-10, 2), rep(10, 2), control = list(seed = 66, NP = 70))
#'  unique(round(res$maxima, 5)) ## all global maxima but not local maxima
#'
#'  ## plotting
#' Shubert_plot <- function(x, y)
#'   Shubert(x = c(x, y))
#' Shubert_plot <- Vectorize(Shubert_plot)
#'  y <- x <- seq(-10, 10, length.out = 40)
#'  persp(x, y, z = outer(X = x, Y = y, FUN =Shubert_plot))
#' @export


ncde <-  function(fn, lower, upper, control = list(),...){

  if (missing(fn))
    stop("'fn' is missing")
  if (missing(fn))
    stop("'lower' is missing")
  if (missing(fn))
    stop("'upper' is missing")

  fn1 <- function(par)
    fn(par, ...)
  #fn1 <- fn # delet later



  # differential weighting factor from interval [0,2]
  if (is.null(control$SF))
    control$SF <- 0.9
  ### algorithms parameter
  #if (is.null(control$st))
  control$st <- 2 ## strategy in DE algorithm
  # crossover probability from interval [0,1]
  if (is.null(control$CR))
    control$CR <- 0.2
  if (is.null(control$NP)) ## humber of intial solutions
    control$NP <- 10 * length(lower)
  NP <- control$NP


  ### comomon param
  if (is.null(control$iter))
    control$iter <- 200
  maxiter <- control$iter
  if(!is.null(control$seed))
    set.seed(control$seed)
  if (is.null(control$hybrid))
    control$hybrid <- TRUE
  if (is.null(control$hybrid_control))
    control$hybrid_control <- NULL

  ##############################################################################
  ##initializing
  ## creating the lower matrix, required for vectorization
  lowermat <- matrix(rep(lower, NP), nrow = NP, byrow = TRUE)
  uppermat <- matrix(rep(upper, NP), nrow = NP, byrow = TRUE)
  D <- length(lower) ##dimension
  ## initialazing the population
  pop <- lowermat + (uppermat - lowermat)* matrix(runif(D * NP), ncol = D, nrow = control$NP)
  # adding the vertices
  vertices <- make_vertices(lower = lower, upper = upper)
  pop[1:dim(vertices)[1], ] <- vertices
  ## initial values
  popval <- apply(pop, 1, fn1)
  nfeval <- NP # umber of function evaluation up to now
  ibest <- which.max(popval) # best solution in population
  ##############################################################################


  for(ii in 2:maxiter){

    ############################################################################
    # neighborhood size
    ## 'm' in table 8
    # In neighborhood mutation, there is only one parameter
    # m which is the neighborhood size. This parameter controls
    # how many individuals are selected in each subpopulation
    # Generally, m should be chosen between 1/20 of the population
    # size and 1/5 of the population size. It can also be dynamically
    # set, from a relatively large value to a small value. Different
    # Fig. 1. Illustration of neighborhood mutation.
    # from other niching parameters, neighborhood size is easy to
    # choose as it can be made proportional to the population size.
    # Moreover, the performance of algorithm is not sensitive to the
    # change of the neighborhood size as evidenced in Table XXI (NCDE).

    if (NP <= 200)
      neigh_size <- 5 + 5 * ((maxiter - ii)/maxiter) else
        neigh_size <- 20 + 30 * ((maxiter- ii)/maxiter)
      # if neigh_size is less than 5 it produces an error in DE when we want to sample pm_id
      neigh_size <- floor(neigh_size)
      ## neigh_size is the number of elemenst in subpoulation
      ##########################################################################

      ##########################################################################
      # computing the euclidean distance
      euc_dist  <- as.matrix(dist(pop, method = "euclidean"))
      ##########################################################################

      ###########################################################################
      ## producing the off spring u_i
      ## producing the ofsprings u_i in each subpopulation_i that is pop[order(euc_dist[j, ]), ][1:neigh_size,]
      ui <- t(sapply(1:NP, function(j)DE(newpop = pop[order(euc_dist[j, ]), ,drop = FALSE][1:neigh_size, , drop = FALSE],
                                         bestpop = pop[j, , drop = FALSE],
                                         target = pop[j, , drop = FALSE],
                                         st = control$st, ## strategy
                                         SF = control$SF,
                                         CR = control$CR,
                                         lower = lower,
                                         upper = upper,
                                         dimension = D)))
      if (D == 1)
        ui <- t(ui)
      ##########################################################################

      ##########################################################################
      # compute the cost function for u_i
      tempval <- apply(ui, MARGIN = 1, fn1)
      nfeval <- nfeval + NP ## update number of function evaluations.
      ##########################################################################

      ##########################################################################
      for(i in 1:NP){

        ########################################################################
        ## randomly perumted the population
        pp <- sample(1:NP, size = NP, replace = FALSE)
        qq <- pop[pp, , drop = FALSE]
        ## calulate the distance of each ui to qq and find the smallest one
        min_id <- which.min(sqrt(apply(sweep(qq, 2, ui[i, , drop = FALSE])^2, 1, sum)))

        if (popval[pp[min_id]]< tempval[i]){
          pop[pp[min_id], ] <- ui[i, , drop = FALSE]
          popval[pp[min_id]] <- tempval[i]
        }
      }
      ##########################################################################

      ##########################################################################
      #### vectorization of comparing the ui with population


      ##########################################################################

      # update index of the global best
      ibest <- which.max(popval)

  }


  ##############################################################################
  # sorting the pop and its corresponding values
  pop <- pop[order(popval, decreasing = TRUE), , drop = FALSE]
  popval <- sort(popval, decreasing = TRUE)
  ##############################################################################

  if (control$hybrid){
    fn2 <- function(par)
      -fn1(par) ## because optim works with min
    #
    temp_optim <- t(sapply(1:NP, FUN =  function(j)optim(par = pop[j, ], fn = fn2, method = "L-BFGS-B", lower = lower, upper = upper,
                                                         control = control$hybrid_control)))
    temp_optim <- t(sapply(1:NP, function(j) c(temp_optim[j, ]$par, temp_optim[j, ]$value)))
    #temp_optim <- round(temp_optim, control$digits)
    #temp_optim <- unique(temp_optim)
    maxima <- temp_optim[, -(D+1)]
    maximaval <- -temp_optim[, (D+1)]
  }else{
    maxima <- maximaval <- NA
  }

  return(list(pop = pop, popval = popval, nfeval = nfeval, maxima = maxima, maximaval = maximaval))
}


###########################################################################################
###########################################################################################

# fn <- function(xx)
# {
#
#
#   factor1=(4-2.1*(xx[1]^2)+(xx[1]^4)/3)*(xx[1]^2)+xx[1]*xx[2]
#   factor2=(-4+4*(xx[2]^2))*(xx[2]^2)
#   y=-4*(factor1+factor2)
#
#   return(y)
# }
#
# lower <- c(-1.9, -1.1)
# upper <- c(1.9, 1.1)
#
#
# #two maxima
# fn(c(0.089842, -0.712656))
# fn(c(-0.089842, 0.712656))
# #
# #
# # control <- list()
#
# test <- ncde(fn, lower, upper, control = list(NP = 200, maxiter = 50))
# max(test$val)
# #
# time1 <- proc.time()
# Rprof(filename = "ncde_rproof.out", line.profiling = T)
# test <- ncde(fn, lower, upper, control = list(seed = 900, NP = 100))
# Rprof(NULL)
# summaryRprof("ncde_rproof.out", lines = "show")
# time1 <- proc.time() - time1


