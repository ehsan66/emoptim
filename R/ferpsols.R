
#' Fitness-Euclidean distance ratio particle swarm optimization (FER-PSO) with local search for multimodal optimization
#'
#'
#' In  FER-PSO, Fitness and Euclidean distance Ratio (FER)
#' is calculated based on the fitness difference and
#' the Euclidean distance between a particle's personal best and
#' other personal bests of the particles in the population. The
#' key advantage is that FER-PSO removes the need of prespecifying
#' niching parameters that are commonly required in
#' existing niching evolutiobary algorithms for multimodal optimization.\cr
#' Furthermore,  a local search technique has been used to enhance the ability to
#' locate most global or local optima.
#'
#' @param fn objective function that should be maximized. Should not return \code{NaN}.
#' @param lower lower bound.
#' @param upper upper bound.
#' @param control control parameters for the algorithm. See "Details".
#' @param ... extra arguments are passed to \code{fn}.
#'
#' @references
#'  Qu, B. Y., Liang, J. J., & Suganthan, P. N. (2012). Niching particle swarm optimization with local search for multi-modal optimization. Information Sciences, 197, 131-143.\cr
#'
#'  Li, X. (2007, July). A multimodal particle swarm optimizer based on fitness Euclidean-distance ratio. In Proceedings of the 9th annual conference on Genetic and evolutionary computation (pp. 78-85). ACM. \cr
#'
#'  Based on a MATLAB code that can be found in Suganthan's home page.
#'
#' @note
#' The function is maximized.\cr
#'
#' Check whether \code{fn} does not return any \code{NaN}.\cr
#'
#' The only stopping rule is the number of iterations.
#'
#' @details
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{\code{swarm}}{swarm size. Defaults to \code{50}.}
#'   \item{\code{iter}}{number of iterations. Defaults to \code{200}.}
#'   \item{\code{w}}{inertia weight. Defaults to \code{0.729843788}.}
#'   \item{\code{c1}}{acceleration factor. Defaults to \code{2.05}.}
#'   \item{\code{c2}}{acceleration factor. Defaults to \code{2.05}.}
#'   \item{\code{local}}{logical; local search should be performed? Defaults to \code{TRUE}}
#'   \item{\code{vectorize_local}}{logical; vectorization for local search? Defaults to \code{TRUE}.}
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
#' @return
#'  a list contains:
#' \describe{
#'   \item{\code{pbest}}{a matrix; position of the particles.}
#'   \item{\code{pbestval}}{a vector; corresponding fitness value of each particle.}
#'   \item{\code{nfeval}}{number of function evaluations.}
#'   \item{\code{maxima}}{position of particles after using \code{L-BFGS-B}. It should be local maxima. If \code{control$ == FALSE}, then it is \code{NA}.}
#'   \item{\code{maximaval}}{fitness values of \code{maxima}.}
#' }
#' @examples
#'####################################################################################
#'## Two-Peak Trap:  global maximum on x = 20 and local maximum on x = 0
#'two_peak <- function(x)
#'  y <- (160/15) * (15 - x) * (x < 15) + 40 * (x - 15) * (x >= 15)
#'
#'ferpsols(two_peak, 0, 20, control = list(seed = 66, swarm = 50))
#'
#'# without local search
#'ferpsols(two_peak, 0, 20, control = list(seed = 66, swarm = 50))
#'
#'# without refining and local search
#'ferpsols(two_peak, 0, 20, control = list(seed = 66, swarm = 50,  hybrid = FALSE))
#'
#'  ####################################################################################
#'# Decreasing Maxima: one global on x = 0.1 and four local maxima
#'dmaxima <- function(x)
#'  y <- exp(-2 * log(2) * ((x - 0.1)/0.8)^2) * (sin(5 * pi * x))^6
#'
#'res <- ferpsols(dmaxima, 0, 1, control = list(seed = 66, swarm = 100))
#'unique(round(res$maxima, 5)) ## ferpsols can not find the local maxima
#'
#'## plot
#'x <- seq(0, 1, length.out = 400)
#'plot(x,  dmaxima(x), type = "l")
#'
#'####################################################################################
#'# Himmelblau's function: four global optima on
#'# x = c(-2.80512, 3.13131), c(3.00000, 2.00000), c(-3.77931, -3.28319) and c(3.58443, -1.84813)
#'Himmelblau <- function(x){
#'  y <- - (x[1]^2 + x[2] - 11)^2 - (x[1] + x[2]^2 - 7)^2
#'  return(y)
#'}
#'
#'
#'res <- ferpsols(Himmelblau, c(-6, -6), c(6, 6), control = list(seed = 66, swarm = 50))
#'unique(round(res$maxima, 5))
#'
#'Himmelblau_plot <- function(x, y)
#'  Himmelblau(x = c(x, y))
#'Himmelblau_plot <- Vectorize(Himmelblau_plot)
#'x <- y <- seq(-6, 6, length.out = 100)
#'persp(x, y, z = outer(X = x, Y = y, FUN = Himmelblau_plot))
#'
#'####################################################################################
#'# Six-Hump Camel Back: two global and two local maxima
#'Six_Hump <- function(x){
#'  factor1 <- (4 - 2.1 * (x[1]^2) + (x[1]^4)/3) * (x[1]^2) + x[1] * x[2]
#'  factor2 <- (-4 + 4 * (x[2]^2)) * (x[2]^2)
#'  y <- -4 * (factor1 + factor2)
#'  return(y)
#'}
#'res <- ferpsols(Six_Hump, c(-1.9, -1.1), c(1.9, 1.1), control = list(seed = 66, swarm = 200))
#'unique(round(res$maxima, 5)) ## can not find the local maxima
#'
#'####################################################################################
#'## 2D Inverted Shubert function :
#'# The global minima: 18 global minima  f(x*) = -186.7309.
#'# the local maxima: sevral
#'
#'Shubert <- function(x){
#'    j <- 1:5
#'    out <- -(sum(j * cos((j + 1) * x[1] + j)) * sum(j * cos((j + 1) * x[2] + j)))
#'    return(out)
#'}
#'
#'res <- ferpsols(Shubert, rep(-10, 2), rep(10, 2), control = list(seed = 66, swarm = 200))
#'unique(round(res$maxima, 5)) ## only the global maxima are found
#'
#'## plotting
#'Shubert_plot <- function(x, y)
#'  Shubert(x = c(x, y))
#'Shubert_plot <- Vectorize(Shubert_plot)
#'y <- x <- seq(-10, 10, length.out = 40)
#'persp(x, y, z = outer(X = x, Y = y, FUN =Shubert_plot))
#' @importFrom lhs randomLHS augmentLHS
#' @importFrom stats dist optim runif
#' @export


ferpsols <-  function(fn, lower, upper, control = list(), ...){


  if (missing(fn))
    stop("'fn' is missing")
  if (missing(fn))
    stop("'lower' is missing")
  if (missing(fn))
    stop("'upper' is missing")


  fn1 <- function(par)
    fn(par, ...)

  #swarm size
  if (is.null(control$swarm))
    control$swarm <- 50
  if (is.null(control$w)) ## inertia weight
    control$w<- 0.729843788
  if (is.null(control$c1))   ##accelaration factor accroding to equation 4
    control$c1 <- 2.05
  if (is.null(control$c2))
    control$c2 <- 2.05
  if (is.null(control$local))
    control$local <- TRUE
  if (is.null(control$vectorize_local))
    control$vectorize_local <- TRUE

  ## common param
  if (is.null(control$iter))
    control$iter <- 200
  if (!is.null(control$seed))
    set.seed(control$seed)
  if (is.null(control$hybrid))
    control$hybrid <- TRUE

  if (is.null(control$hybrid_control))
    control$hybrid_control <- NULL




  # parameters fixed during the run
  # max velocity and minus min veloicty
  maxvel <- 0.5 * (upper-lower)
  swarm <- control$swarm


  ####################################################################
  ## max velocity as a matrix, needed for vectorization
  maxvelmat <- matrix(rep(maxvel, swarm), nrow = swarm, byrow = TRUE)


  #lower matrix and upper matrix
  lowermat <- matrix(rep(lower, swarm), nrow = swarm, byrow = TRUE)
  uppermat <- matrix(rep(upper, swarm), nrow = swarm, byrow = TRUE)
  #####################################################################

  ##size of the search space Eq. 7. sos is  ||s|| in the numerator of alpha Li (2007)
  sos <- sqrt(sum((upper - lower)^2))
  D <- length(lower) ##dimension


  #####################################################################
  ## initializing
  ### we want to add the vertices to the inital positions later
  vertices <- make_vertices(lower = lower, upper = upper)
  A <- lhs::randomLHS(floor((swarm - dim(vertices)[1])/2), length(lower))
  pos <- lhs::augmentLHS(A, ceiling((swarm - dim(vertices)[1])/2))

  vrunif <- Vectorize(runif)
  #   pos <- vrunif (swarm, lower, upper)
 #pos[1:dim(vertices)[1], ] <- vertices
  pos <- rbind(pos, vertices)

  # fitness
  fit <- apply(pos, 1, fn1)
  ##initializing the fitness count

  nfeval <- swarm
  pbest <- pos
  pbestval <- fit
  #initialize the velocity of the particles
  vel <- lowermat + (2 * uppermat * vrunif(swarm, rep(0, D), rep(1, D)))
  #####################################################################


  for(iter in 2:(control$iter)){

    # global best fitness value and id
    gbestval <- max(fit)
    gbestid <- which.max(fit)
    gworstval <- min(fit)
    gworstid <- which.min(fit)

    # compute the scaling factor (alpha in equation 7)
    sf <- sos/(gbestval-gworstval)

    # inital the nbest position
    nbest <- pbest

    # inital the nbest value
    nbestval <- pbestval


    # update the nbest and nbestval
    ##caluclate the euclidean distance for all pbest
    eucdist <- dist(pbest, method = "euclidean")
    ### update_nbest update the nbest and nbest value
    nbesttemp <- update_nbest(eucdist = as.matrix(eucdist) , pbest = pbest,
                              nbest = nbest, pbestval = pbestval, nbestval = nbestval, sf = sf)


    nbest <- nbesttemp$nbest
    nbestval <- nbesttemp$nbestval


    ## update velocity
    temp1 <- control$c1 * vrunif(swarm, rep(0, D), rep(1, D)) * (pbest - pos) + control$c2 * vrunif(swarm, rep(0, D), rep(1, D)) * (nbest - pos)
    vel   <- control$w * vel + temp1
    temp1 <- NA ## cleaning


    #######################################################################################################
    # velocity
    ##check the upper bound
    vel <- (vel > maxvelmat) * maxvelmat + (vel <= maxvelmat) * vel
    ##check the lower bound
    vel <- (vel < -maxvelmat) * -maxvelmat + (vel >= -maxvelmat) * vel
    ########################################################################################################


    ########################################################################################################
    ##update position
    pos <- pos + vel
    ## correct out of bounds value
    pos <-  ((pos >= lowermat) & (pos <= uppermat)) * pos + ##if inside the box return pos
      (pos < lowermat) * (lowermat + 0.25 * (uppermat - lowermat) * vrunif(swarm, rep(0, D), rep(1, D))) + #lower bound
      (pos > uppermat) * (uppermat - 0.25 * (uppermat - lowermat) * vrunif(swarm, rep(0, D), rep(1, D))) ##upper bound
    ########################################################################################################

    ########################################################################################################
    # computing the fitness values and updating the pbes and pbestval
    fit <- apply(pos, 1, fn1)
    nfeval <- swarm + nfeval
    temp2 <- (pbestval > fit)
    pbestval <- temp2 * pbestval + (1 - temp2) * fit
    pbest <- temp2 * pbest + (1 - temp2) * pos
    ########################################################################################################


    ########################################################################################################
    ## local search
    if (control$local && control$vectorize_local){

      temp_euc <-  as.matrix(dist(pbest, method = "euclidean"))
      temp_near_id <- sapply(1:swarm, FUN = function(j)which(order(temp_euc[j, ]) == 2))


      near <- pbest - pbest[temp_near_id , , drop = FALSE] * (pbestval >= pbestval[temp_near_id]) +
        (pbest[temp_near_id , , drop = FALSE] - pbest) * (pbestval < pbestval[temp_near_id])

      temp_pbest <- pbest + 2 * matrix(runif(D * swarm), nrow = swarm, ncol = D) * near


      temp_pbest <- ((temp_pbest >= lowermat) & (temp_pbest <= uppermat)) * temp_pbest +
        (temp_pbest < lowermat) * (lowermat + 0.25 * (uppermat - lowermat) * vrunif(swarm, rep(0, D), rep(1, D))) + #lower bound
        (temp_pbest > uppermat) * (uppermat - 0.25 * (uppermat - lowermat) * vrunif(swarm, rep(0, D), rep(1, D)))


      temp_pbestval <- apply(temp_pbest, 1, fn1)
      nfeval <- nfeval + swarm ## updating the number of fitness evaluations
      temp_improve <- temp_pbestval >  pbestval
      ## updating pbest and pbestval if there is improvement
      pbest <- pbest * (1-temp_improve) + temp_pbest * temp_improve
      pbestval <- pbestval * (1-temp_improve) + temp_pbestval * temp_improve
    }


    if (control$local && !control$vectorize_local){
      ##local search
      for (k in 1:swarm){
        temp_euc <- order(sqrt(apply(sweep(pbest, 2, pbest[k, , drop = FALSE])^2, 1, sum)))
        temp_near_id <- which(temp_euc == 2)

        ## producing a point near the pbest[k]
        near <- pbest[k, , drop = FALSE] - pbest[temp_near_id , , drop = FALSE] * (pbestval[k] >= pbestval[temp_near_id]) +
          (pbest[temp_near_id , , drop = FALSE] - pbest[k, , drop = FALSE]) * (pbestval[k] < pbestval[temp_near_id])

        temp_pbest <- pbest[k, , drop = FALSE] + 2 * runif(D) * near
        temp_pbest  <- ((temp_pbest >= lower) & (temp_pbest <= upper)) * temp_pbest +
          (temp_pbest < lower) * (lower + 0.25 * (upper - lower)*runif(D)) + #lower bound
          (temp_pbest > upper) * (upper - 0.25 * (upper - lower)*runif(D)) # upper bound

        temp_pbestval <- fn1(temp_pbest)
        nfeval <- nfeval + 1
        temp_improve <- temp_pbestval >  pbestval[k]
        ## updating the number of fitness evaluations
        pbest[k, ] <- pbest[k, , drop = FALSE] * (1-temp_improve) + temp_pbest * temp_improve
        pbestval[k] <- pbestval[k] * (1-temp_improve) + temp_pbestval * temp_improve

      }
    }
  } ## end of local search
  ########################################################################################################

  pbest <- pbest[order(pbestval, decreasing = TRUE), , drop = FALSE]
  pbestval <- sort(pbestval, decreasing = TRUE)


  ########################################################################################################
  if (control$hybrid){
    ##### now we imrove all the particles with optim
    fn2 <- function(par)
      -fn1(par) ## because optim works with min

    temp_optim <- t(sapply(1:swarm, FUN =  function(j)optim(par = pbest[j, ], fn = fn2, method = "L-BFGS-B",
                                                            lower = lower, upper = upper,
                                                            control= control$hybrid_control)))
    temp_optim <- t(sapply(1:swarm, function(j) c(temp_optim[j, ]$par, temp_optim[j, ]$value)))
    #temp_optim <- round(temp_optim, control$digits)
    #temp_optim <- unique(temp_optim)
    maxima <- temp_optim[, -(D+1)]
    maximaval <- -temp_optim[, (D+1)]
  } else
    maxima <- maximaval <- NA

  ########################################################################################################

  return(list(pbest = pbest, pbestval = pbestval, nfeval = nfeval, maxima = maxima, maximaval = maximaval))
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
#
#
# control = list()
#
# time1 <- proc.time()
# Rprof(filename = "ferpso_rproof.out", line.profiling = T)
# test <- ferpso(fn, lower, upper, control = list(local = F, seed = 900, swarm = 500))
# Rprof(NULL)
# summaryRprof("ferpso_rproof.out", lines = "show")
# time1 <- proc.time() - time1

###########################################################################################
###########################################################################################

# update_nbest(eucdist = eucdist,
#              pbest = pbest, nbest = nbest, pbestval = pbestval, nbestval = nbestval, sf = sf)


# Two-Peak Trap
#ferpso(fn1, 0, 20, control = list(local = TRUE, seed = 66, swarm = 50))

