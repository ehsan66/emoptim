# upate the new position based on DE
DE <- function(target, newpop, bestpop, st = 2, SF = .9, CR = .1, lower, upper, dimension){
  # target: target vector
  # pop: population
  # bestpop: the best pop
  # st: strategy
  # CR: control parameter of DE that decides in a comparison with random number. [0, 1]
  # SF: differential weighting factor from interval [0,2]
  # return a vector of trial


  pm_id <- sample(1:dim(newpop)[1], size = 5, replace = FALSE)
  pm1 <- newpop[pm_id[1], ]
  pm2 <- newpop[pm_id[2], ]
  pm3 <- newpop[pm_id[3], ]
  pm4 <- newpop[pm_id[4], ]
  pm5 <- newpop[pm_id[5], ]

  mui <- runif(dimension) < CR
  # randomly choosen integer from 1:dim(newpop)[1]

  # if all(mui) == FALSE the we have ui = target, to avoid this we set one of the mui by random to 1
  if(!all(mui)){
    mui[sample(1:dimension, 1)] <- TRUE
  }



  if (st == 1){                # DE/best/1   6
    ui = bestpop + SF * (pm1 - pm2)      # differential variation
    ui = (target * !mui) + (ui * mui)    # binomial crossover
  }
  if (st == 2){                           # DE/rand/1
    ui = pm3 + SF * (pm1 - pm2)       # differential variation
    ui = (target * !mui) + (ui * mui)          # crossover
  }

  if (st == 3)
    stop("The strategy 3 is not implemented. Please contact the maintainer.")
  # if (st == 3){                  # DE/rand-to-best/1    8
  #   ui = popold + SF * (bestpop - popold) + SF * (pm1 - pm2)
  #   ui = (target * !mui) + (ui * mui)     # crossover
  # }
  if (st == 4){                  # DE/best/2           9
    ui = bestpop + SF * (pm1 - pm2 + pm3 - pm4)  # differential variation
    ui = (target * !mui) + (ui * mui)           # crossover
  }
  if (st == 5){                  # DE/rand/2           10
    ui = pm5 + SF * (pm1 - pm2 + pm3 - pm4)  # differential variation
    ui = (target * !mui) + (ui * mui)            # crossover
  }

  ## correct out of bounds values
  # if (runif(1) > 0.5){
  #   ui <- (ui < lower) * lower+ (ui >= lower) * ui
  #   ui <- (ui > upper) * upper + (ui <= upper) * ui
  # }else{
  ui <- (ui < lower) * (lower + runif(dimension) * (upper - lower))+ (ui >= lower) * ui
  ui <- (ui > upper) * (lower + runif(dimension) * (upper - lower))+ (ui <= upper) * ui
  # }

  return(ui)
}
