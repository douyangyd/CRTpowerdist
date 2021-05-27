### package needed
### Run helper function file first
require(arrangements)
require(doParallel)
require(parallel)
require(lme4)
require(ggplot2)
require(doRNG)
require(compiler)

### Helper 1: Convert allocations to treatment indicator and calculate
### weights
rand <- function(I, J, allocs, n) {
  ord <- c()
  for (i in 1:length(I)) {
    deno <- c()
    for (j in 1:J) {
      temp <- c()
      temp <- rep(j, allocs[i, j, n])
      deno <- c(deno, temp)
    }
    ord <- c(ord, deno)
  }
  return(ord)
}

multi <- function(M, N, allocs, k) {
  # k is the kth unique allocation you want to evaluate
  temp <- c()
  final <- c()
  for (i in 1:length(M)) {
    for (j in 1:length(N)) {
      temp[j] <- factorial(allocs[i, j, k])
      deno <- prod(temp)
    }
    final[i] <- factorial(M[i])/deno
  }
  return(prod(final))
}


### Helper 2: function to generate all unique allocation when length(I)
### >3
alloc1 <- function(I.now, period = 1, stem = vector(), store.alloc = F) {
  if (period == J) {
    complete.alloc <- cbind(stem, I.now)
    if (all(apply(complete.alloc, 1, sum) == I)) {
      count <<- count + 1
      # print(paste('count = ', count))
      if (store.alloc) {
        allocs <<- c(allocs, complete.alloc)
      }
    }
  } else {
    min.type <- pmax(0, N[period] - sum(I.now) + I.now)
    max.type <- pmin(I.now, N[period])
    # print(paste('min.type = ', min.type)) print(paste('max.type = ',
    # max.type))
    n.per.type <- max.type - min.type + 1
    n.alloc <- vector()
    # print(paste('n.per.type', paste(n.per.type))) print(prod(n.per.type))
    for (j in 0:(prod(n.per.type) - 1)) {
      # print(paste('j =', j))
      for (i in 1:M) {
        if (i == 1)
          n.alloc[i] <- j%/%prod(n.per.type[(i + 1):M]) + min.type[i] else if (i == M)
            n.alloc[i] <- j%%(n.per.type[M]) + min.type[i] else n.alloc[i] <- ((j%/%prod(n.per.type[(i + 1):M]))%%n.per.type[i]) +
                min.type[i]
      }
      # print(n.alloc)
      if (sum(n.alloc) == N[period])
        alloc1(I.now = I.now - n.alloc, period = period + 1, stem = cbind(stem,
                                                                          n.alloc), store.alloc = store.alloc)
    }
  }
}


### Helper 3: funnction to generate all unique possible allocations.Each
### rwo represents a unique allocation. Column number are cluster id
### (arranged as (S,M,L) or (S,L)). Entries are at which period the
### cluster transit.

all_allocs <- function(I, J, P) {
  if (length(unique(P)) != 1) {
    count <<- 0
    allocs <<- vector()
    M <- length(I)
    N <- P
    environment(alloc1) <- environment()
    alloc1(I.now = I, store.alloc = T)
    dim(allocs) <- c(length(I), J, count)
  } else {
    if (length(I) == 1) {
      allocs <- c(rep(1:J, each = I/J))
    }
    if (length(I) == 2) {
      X0 <- permutations(c(0:min(I[2], sum(I)/J)), k = (J - 1), replace = T)
      XL <- cbind(X0, apply(X0, 1, function(x) {
        ifelse(sum(x) <= I[2] && sum(x) >= I[2] - sum(I)/J, I[2] -
                 sum(x), NA)
      }))
      XL <- XL[!rowSums(!is.finite(XL)), ]
      XS <- matrix(sum(I)/J, nrow = dim(XL)[1], ncol = dim(XL)[2]) -
        XL
      allocs <- list()
      for (i in 1:dim(XL)[1]) {
        allocs[[i]] <- rbind(XS[i, ], XL[i, ])
      }

      allocs <- array(unlist(allocs), dim = c(length(I), J, dim(XL)[1]))
    }
    if (length(I) == 3) {
      X0 <- permutations(c(0:min(I[3], sum(I)/J)), k = (J - 1), replace = T)
      X1 <- cbind(X0, apply(X0, 1, function(x) {
        ifelse(sum(x) <= I[3] && sum(x) >= I[3] - sum(I)/J, I[3] -
                 sum(x), NA)
      }))
      X1 <- X1[!rowSums(!is.finite(X1)), ]  #X1 gives all the ways the large clusters can be assigned

      X2 <- permutations(c(0:min(I[2], sum(I)/J)), k = (J - 1), replace = T)
      X3 <- cbind(X2, apply(X2, 1, function(x) {
        ifelse(sum(x) <= I[2] && sum(x) >= I[2] - sum(I)/J, I[2] -
                 sum(x), NA)
      }))
      X3 <- X3[!rowSums(!is.finite(X3)), ]  #X3 gives all the ways the medium clusters can be assigned

      X6 <- matrix(0, ncol = 2 * J)
      for (i in 1:dim(X1)[1]) {
        X4 <- matrix(X1[i, ], nrow = dim(X3)[1], ncol = J, byrow = T)
        X5 <- cbind(X4, X3, apply(X3 + X4, 1, function(x) {
          ifelse(max(x) <= sum(I)/J && sum(x) == (sum(I) - I[1]),
                 sum(x), NA)
        }))
        X5 <- X5[!rowSums(!is.finite(X5)), ]
        X6 <- rbind(X6, X5[, 1:(2 * J)])  #X6 is a dummy matrix that contains all the possible unique allocations; columns 1 to p for large clusters, and columns (p+1) to 2p for medium clusters
      }

      X6 <- X6[2:(dim(X6)[1]), ]
      XS <- matrix(sum(I)/J, ncol = J, nrow = dim(X6)[1]) - X6[,
                                                               1:J] - X6[, (J + 1):(2 * J)]  #allocations for small clusters
      XL <- X6[, 1:J]  #allocations for large clusters
      XM <- X6[, (J + 1):(2 * J)]  #allocations for medium clusters
      #### ========================================================


      allocs <- list()
      for (i in 1:dim(X6)[1]) {
        allocs[[i]] <- rbind(XS[i, ], XM[i, ], XL[i, ])
      }
      allocs <- array(unlist(allocs), dim = c(length(I), J, dim(X6)[1]))
    }
    if (length(I) > 3) {
      count <<- 0
      allocs <<- vector()
      M <- length(I)
      N <- rep(sum(I)/J, J)
      environment(alloc1) <- environment()
      alloc1(I.now = I, store.alloc = T)
      dim(allocs) <- c(length(I), J, count)
    }

  }

  if (length(I) == 1) {
    rand.order <- t(as.matrix(allocs))
    w <- 1
  } else {
    rand.order <- matrix(0, ncol = sum(I), nrow = dim(allocs)[3])
    for (i in 1:dim(allocs)[3]) {
      rand.order[i,] <- rand(I, J, allocs, i)
    }

    multiplicity <- c()
    for (i in 1:dim(allocs)[3]) {
      multiplicity[i] <- multi(I, J, allocs, i)
    }
    w <- multiplicity/sum(multiplicity)  ## weights
  }
  return(list(order = rand.order, weights = w))
}


### Helper 3 a function from hfunk package
freq <- function(x, elts) {
  res <- sapply(elts, function(el) sum(x == el))
  names(res) <- elts
  return(res)
}


### Helper 4 for generating a list of design matrices
DesignMatrix <- function(I, J, P, K, S, factor.time = FALSE, user.allocs = NULL, design, strat) {
  # I = vector of cluster types J = # of steps (excluding baseline) K =
  # vector of cluster sizes factor.time = logical; for whether time
  # should be considered a factor user.allocs = NULL to generate design
  # matrices for all unique allocations; otherwise only for the specified
  # allocation
  # S = indicator for stratum
  # design = pcrt or sw
  # strat = logical for whether there is stratification



  if (!is.null(user.allocs)) {
    allocs <- user.allocs
    weights = NULL
  } else {
    if (strat){
      all_allocs_out = all_allocs_strat(I=I, P=P, S=S, K=K)
    } else {
      all_allocs_out <- all_allocs(I, J, P)
    }
    weights = all_allocs_out$weights
    allocs = all_allocs_out$order
  }

  if (!(class(allocs) == "matrix" || class(allocs) == "data.frame")) allocs <- matrix(allocs, nrow = 1)

  if (design == "pcrt"){

    if(!all(K >= 1)){warning("Some clusters have size < 1 and may be empty with no observations", immediate. = T)}
    size <- K
    #size <- c()
    #for (i in 1:length(K)){
    #  size[i] <- sample(c(floor(K[i]),floor(K[i])+1),
    #                    prob= c(1-(K[i]-floor(K[i])),
    #                            (K[i]-floor(K[i]))),
    #                    1, replace=T)
    #}


    XMat <- lapply(1:nrow(allocs), FUN = function(ii) {
      alloc_i <- (allocs[ii, ])

      Tx <- unlist(sapply(1:sum(I), function(kk) {
        rep(ifelse(alloc_i[kk] == 1, 0, 1), size[kk])
      }, simplify = F))

      cbind(Tx, Intercept = 1)
    })


  }


  if (design == "sw"){
    #size = K * (J+1)
    size = K
    if(!all(K >= 1)){warning("Some clusters have size < 1 and may be empty with no observations", immediate. = T)}

    #for (i in 1:length(K)){
    #  size[i] <- sample(c(floor(K[i]),floor(K[i])+1),
    #                    prob= c(1-(K[i]-floor(K[i])),
    #                            (K[i]-floor(K[i]))),
    #                    1, replace=T)
    #}

      #Period <- unlist(lapply(1:length(K), FUN= function(jj){
      #  rep(0:J,each=K[jj],replace = T)}))

    Period <- unlist(lapply(1:length(K), FUN= function(jj){
        1:K[jj] %% (J+1)}))

    #if(all(floor(K / (J + 1)) == K / (J + 1))){
    #  Period <- unlist(lapply(1:length(size), FUN= function(jj){
    #    rep(0:J,each=size[jj]/(J+1),replace = T)}))
    #} else {
    #  Period <- unlist(lapply(1:length(size), FUN= function(jj){
    #    sample(0:J,size[jj],replace = T)}))
    #}

    # Cluster <- unlist(sapply(1:sum(I), function(kk) {
    #   rep(kk, size[kk])
    # }, simplify = F))

    XMat <- lapply(1:nrow(allocs), FUN = function(ii) {
      alloc_i <- (allocs[ii, ])

      Tx <- ifelse(unlist(sapply(1:sum(I), function(kk) {
        rep(alloc_i[kk], size[kk])
      }, simplify = F)) > Period, 0, 1)

      if (factor.time) {
        Period <- as.factor(Period)
        Period <- model.matrix(~. + 0, as.data.frame(Period))
        cbind(Tx, Intercept = 1, Period[, -1])
      } else cbind(Tx, Intercept = 1, Period)
    })
  }
  return(list(XMat=XMat, size=size, weights=weights))
}





# function to generate unique allocations for stratified randomization
all_allocs_strat <- function(I, P, S, K) {
  # identify levels of stratification
  strat_name = levels(S)

  # generate list of number of cluster with different sizes
  param_strat_list <- lapply(strat_name, function(name_i) {
    K_strat = K[which(S == name_i)]
    I_strat = sapply(unique(K), function(kk){
      length(which(K_strat == kk))
    })

    if(!all(I_strat != 0)) {
      if(length(I_strat) == 2) I_strat = c(I_strat, 0, 0)
      if(length(I_strat) == 3) I_strat = c(I_strat, 0)
    }
    P_strat = P[[which(names(P) == name_i)]]
    id_strat = which(S == name_i)
    list(I_strat = I_strat, P_strat = P_strat, id_strat = id_strat )
  })


  J = length(P[[1]])

  # generate allocations within each stratum
  alloc_strat_list <- lapply(param_strat_list, function(param_strat) {
    all_allocs(I = param_strat$I_strat, J = J, P = param_strat$P_strat)$order

  })

  # all possible combinations of unique allocations from each stratum
  n.alloc = list()
  for (ii in 1:length(alloc_strat_list)) {n.alloc[[ii]] = 1:nrow(alloc_strat_list[[ii]])}
  combi = expand.grid(n.alloc)

  # convert allocations to array form
  array_list = mapply(alloc_strat_list, param_strat_list, FUN = function(alloc_strat, param_strat){
    nn = nrow(alloc_strat)
    if(length(param_strat$I_strat) != length(I)) param_strat$I_strat = param_strat$I_strat[1:length(I)]
    clus.id = rep(1:length(param_strat$I_strat), times = param_strat$I_strat)

    # list of allocations grouped by cluster type
    alloc_byclus <- lapply(1:length(param_strat$I_strat), function(ii) {
      alloc_ii <- matrix(alloc_strat[, which(clus.id == ii)], nrow = nrow(alloc_strat))
    })

    # convert back to array type
    alloc_array <- array(dim = c(length(param_strat$I_strat), J, nn))

    for (ii in 1:nn) {
      alloc_array[, , ii] <- t(sapply(alloc_byclus, function(alloc_byclus) {
        table(factor(alloc_byclus[ii, ], levels = 1:J))
      }))
    }
    alloc_array
  }, SIMPLIFY = F)

  # combine across strata
  alloc_merge =lapply(1:nrow(combi), function(ii){
    temp = list()
    for (jj in 1:ncol(combi)){
      temp[[jj]] = array_list[[jj]][, , combi[ii, jj]]
    }
    Reduce("+", temp)

  })

  # only keep unique allocations
  alloc_merge = unique(alloc_merge)

  # convert back to matrix form
  alloc = sapply(1:length(alloc_merge), function(n){
    ord <- c()
    for (i in 1:length(I)) {
      deno <- c()
      for (j in 1:J) {
        temp <- c()
        temp <- rep(j, alloc_merge[[n]][i,j])
        deno <- c(deno, temp)
      }
      ord <- c(ord, deno)
    }
    ord
  })

  alloc = t(alloc)

  # weight calculations...
  # not sure about this part, since multiplicity under stratification
  # may not be identical to that without stratification
  temp <- array(unlist(alloc_merge), dim = c(length(I), J, nrow(alloc)))
  multiplicity <- c()
  for (i in 1:dim(temp)[3]) {
    multiplicity[i] <- multi(I, J, temp, i)
  }
  weight <- multiplicity/sum(multiplicity)  ## weights
  return(list(order = alloc, weights = weight))
}

###
rand <- compiler::cmpfun(rand)
multi <- compiler::cmpfun(multi)
alloc1 <- compiler::cmpfun(alloc1)
all_allocs <- compiler::cmpfun(all_allocs)
freq <- compiler::cmpfun(freq)
DesignMatrix <- compiler::cmpfun(DesignMatrix)
all_allocs_strat <- compiler::cmpfun(all_allocs_strat)

##-------------------------------------------------------
## Functions used to do power calculations and construct
## power distribution. These functions are wrapped in main
## functions
##-------------------------------------------------------

### Simulation based calculation

sim.ap <- function(I, P, K, mu0, Tx.effect, Time.effect = NULL, factor.time = FALSE,
                   design, user.allocs, n.sims, rho = NULL, sigma.a = NULL, sigma.e = NULL,
                   sig.level = 0.05, family = "gaussian") {
  ### I: number of clusters of each type (categorized by cluster sizes) J:
  ### number of transitioning period K: cluster size (Length of I must be
  ### the same as lenght of K) mu0: Baseline mean/rate (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes) Tx.effect:
  ### treatment effect (Linear scale for gaussian outcome, natural scale
  ### for binary/count outcomes) Time.effect: time effect (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes)
  ### factor.time: T/F whether treat time as continous or categorical
  ### variable design: sw or pcrt user.allocs: The randomization order(s)
  ### to evaluate. Entered as matrix. Each row represent a randomization
  ### order the user try to evaluate. The entries are when the cluster
  ### transits from control to treatment group.  n.sims: the number of
  ### trials the user want to simulate to calculation power rho: ICC
  ### sigma.a: between-cluster variance sigma.e: within-cluster variance
  ### (Only two of rho,sigma.a and sigma.e should be provided) sig.level:
  ### significance level (Default=0.05) family: one of gaussian, binomial
  ### or poisson (case sensitive)

  ### check input vadility
  if (sum(I) != length(K)) {
    stop("Number of type of clusters must match the number of cluster sizes")
  }
  if (family != "gaussian" & family != "binomial" & family != "poisson") {
    stop("Family must be one of gaussian, binomial or poisson")
  }
  ### sigma and rho conversion-----
  if (family == "gaussian") {
    if(sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) !=1)
      stop("Exactly one of rho, sigma.a, sigma.e must be null")

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(sigma.e)) { sigma.e <- sqrt(sigma.a^2 * (1 - rho)/rho) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "binomial"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Exactly one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline risk and treatment effect instead of the input value")}

    mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))

    sigma.e <- sqrt((mu0 * (1 - mu0) + mu1 * (1 - mu1))/2)

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "poisson"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Only one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline rate and treatment effect instead of the input value")}
    mu1 <- Tx.effect * mu0

    sigma.e <- sqrt((mu0 + mu1)/2)

    if (is.null(sigma.a)) {
      sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho))}
    if (is.null(rho)) {
      rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }

  }
  ### -----
  if (design != "sw" & design != "pcrt") {
    stop("Design should be one of sw or pcrt")
  }
  if (design == "pcrt" & length(P)>2) {
    stop("Parallel CRT design can only have two groups")
  }


  if (length(sig.level) > nrow(user.allocs)) {
    stop("You entered too many sig.level")
  }

  if (length(Time.effect)>1 & length(Time.effect) != (length(P))) {
    stop("length of Time.effect should match the number of transition steps")
  }

  sample.allocs <- user.allocs
  if (length(sig.level) == 1){
    sig.level <- rep(sig.level,nrow(sample.allocs))
  } else {sig.level <- sig.level}
  ### Parallel CRT

  if (design == "pcrt") {
    TTC <- NULL
    J <- 2
    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list )
        ### generate conuterfactural
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        y_trt <- mu0 + clus.mean + Tx.effect  + error.y
        y_ctr <- mu0 + clus.mean + error.y
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            m <-  lmer(Y~ Treatment +(1| Cluster),  control=lmerControl(check.conv.grad = .makeCC("ignore", tol = 1e-3, relTol = NULL)), data = analytic.data)
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j]=2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"=pval, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$pval <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- list("attained.power.sim"=pval,"TGI"=abs(con-int),"SE"=SE,"trt.estimate"=coef)
    }

    ### Binary outcomes
    if (family == "binomial") {
      mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean +
          b1
        y0 <- b0 + clus.mean
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        ### generate conuterfactural
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "binomial",
                             analytic.data, control=glmerControl( check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- list("attained.power.sim"=pval,"TGI"=abs(con-int),"SE"=SE,"trt.estimate"=coef)
    }

    ### count outcomes
    if (family == "poisson") {
      mu1 <- mu0 * Tx.effect
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        ### generate conuterfactural
        y1 <- b0 + clus.mean + b1
        y0 <- b0 + clus.mean
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "poisson",
                             analytic.data, control= glmerControl( check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),)
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- list("attained.power.sim"=pval,"TGI"=abs(con-int),"SE"=SE,"trt.estimate"=coef)

    }
  }

  ### SW design
  if (design == "sw") {
    J <- length(P)
    K <- K
    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list)
        if (is.null(Time.effect)) {
          Time.effect <- 0
        } else Time.effect <- Time.effect
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        ### generate conuterfactural
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        if (length(Time.effect)>1){
          Time.effect <- c(0,Time.effect) # add 0 to baseline
          y_trt <- mu0 + clus.mean + Tx.effect + Time.effect[period+1] + error.y
          y_ctr <- mu0 + clus.mean + Time.effect[period+1] + error.y
        } else {
          y_trt <- mu0 + clus.mean + Tx.effect + period * Time.effect + error.y
          y_ctr <- mu0 + clus.mean + period * Time.effect + error.y}
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          cor <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <-  lmer(Y~ factor(Period) + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            } else {
              m <-  lmer(Y~ Period + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            }
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j] <- 2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"= pval,"con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$pval <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- list("attained.power.sim"=pval,"SE"=SE,"trt.estimate"=coef,"TGI"=abs(con-int),"TTC"=timecor)
    }

    ### Binary outcomes
    if (family == "binomial") {
      mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list )
        clus.id <- rep(1:sum(I), clus.size.list )
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        ### generate conuterfactural
        if (length(Time.effect)>1){
          b2 <- c(0, b2) # add 0 to baseline
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + b2 * period
          y0 <- b0 + clus.mean + b2 * period
        }
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims, ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "binomial", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "binomial", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                          check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- list("attained.power.sim"=pval,"SE"=SE,"trt.estimate"=coef,"TGI"=abs(con-int),"TTC"=timecor)
    }

    ### count outcomes
    if (family == "poisson") {
      mu1 <- mu0 * Tx.effect
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list)
        clus.id <- rep(1:sum(I), clus.size.list )
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        ### generate conuterfactural
        if (length(Time.effect)>1){
          b2 <- c(0, b2) # add 0 to baseline
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + period * b2
          y0 <- b0 + clus.mean + period * b2
        }
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "poisson", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "poisson", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- list("attained.power.sim"=pval,"SE"=SE,"trt.estimate"=coef,"TGI"=abs(con-int),"TTC"=timecor)
    }
  }
  return(pred.data)
  #return(list(results = pred.data, input = list(I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, factor.time = factor.time,
  #                                              design = design, n.sims = n.sims, rho = rho,
  #                                              sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig.level,
  #                                              family = family)))
}

sim.pd <- function(I, P, K, mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL, factor.time = FALSE,
                   design, gen.all = TRUE, n.allocs = NULL, n.sims, rho = NULL,
                   sigma.a = NULL, sigma.e = NULL, print.allocs = F, plot = FALSE, sig.level = 0.05,
                   family = "gaussian") {
  ### I: number of clusters of each type (categorized by cluster sizes) P:
  ### number of transiting at each period, length of P will be the total
  ### number of transition period K: cluster size (Length of I must be the
  ### same as lenght of K) mu0: Baseline mean/rate (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes) Tx.effect:
  ### treatment effect (Linear scale for gaussian outcome, natural scale
  ### for binary/count outcomes) Time.effect: time effect (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes)
  ### factor.time: T/F whether treat time as continous or categorical
  ### variable design: choice of pcrt or sw gen.all: list all unique
  ### allocation or not n.allocs: The number of randomization order(s) to
  ### evaluate. It could be a 'A' where all possible unique allocations
  ### will be evaluated. If a number is entered, a sample of unique
  ### allocations will be sampled and evaluate n.sims: the number of trials
  ### the user want to simulate to calculation power rho: ICC sigma.a:
  ### between-cluster variance sigma.e: within-cluster variance (Only two
  ### of rho,sigma.a and sigma.e should be provided) print.allocs: TRUE or
  ### FALSE. Whether print the allocs has been evaluated plot: T/F display
  ### histogram of power distribution or not sig.level: significance level
  ### (Default=0.05) family: one of gaussian, binomial or poisson (case
  ### sensitive)

  if ((length(unique(K) == 1) == T)) {n.allocs = 1}

  ### check input vadility
  if (sum(I) != length(K)) {
    stop("Number of type of clusters must match the number of cluster sizes")
  }
  if (family != "gaussian" & family != "binomial" & family != "poisson") {
    stop("Family must be one of gaussian, binomial or poisson")
  }
  if (sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) == 0) {
    stop("One of rho, sigma.a, sigma.e must be null")
  }
  ### sigma and rho conversion-----
  if (family == "gaussian") {
    if(sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) !=1)
      stop("Exactly one of rho, sigma.a, sigma.e must be null")

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(sigma.e)) { sigma.e <- sqrt(sigma.a^2 * (1 - rho)/rho) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "binomial"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Exactly one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline risk and treatment effect instead of the input value")}

    mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))

    sigma.e <- sqrt((mu0 * (1 - mu0) + mu1 * (1 - mu1))/2)

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "poisson"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Only one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline rate and treatment effect instead of the input value")}
    mu1 <- Tx.effect * mu0

    sigma.e <- sqrt((mu0 + mu1)/2)

    if (is.null(sigma.a)) {
      sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho))}
    if (is.null(rho)) {
      rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }

  }
  ### -----
  if (design != "sw" & design != "pcrt") {
    stop("Design should be one of sw or pcrt")
  }
  if (design == "pcrt" & length(P)>2) {
    stop("Parallel CRT design can only have two groups")
  }

  if (length(I) ==1 & length(unique(P)) != 1){
    stop("Please consider function sim.ap as the parameters you entered only have one unique allocation")
  }
  if (gen.all == FALSE & is.null(n.allocs)) {
    stop("Evaluating all allocations is impossible without listing out all allocations")
  }
  if (length(Time.effect)>1 & length(Time.effect) != (length(P))) {
    stop("length of Time.effect should match the number of transition steps")
  }

  if (design == "pcrt") {
    pred.power <- FALSE
    TTC <- NULL
    J <- 2
    if (gen.all == FALSE) {
      pred.power <- FALSE
      sample.allocs <- t(replicate(n.allocs, sample(0:1, sum(I),
                                                    prob = c(1/2, 1/2), replace = T)))
      index <- split(1:sum(I), rep(1:length(I), I), drop = TRUE)
      all <- list()
      for (n in 1:nrow(sample.allocs)) {
        all[[n]] <- t(sapply(1:length(I), FUN = function(i) {
          freq(sample.allocs[n, index[[i]]], 1:2)
        }))
      }
      temp <- array(unlist(all), dim = c(length(I), 2, nrow(sample.allocs)))
      multiplicity <- c()
      for (i in 1:dim(temp)[3]) {
        multiplicity[i] <- multi(I, 2, temp, i)
      }
      sample.weight <- multiplicity/sum(multiplicity)  ## weights
    }
    all.ordwgh <- all_allocs(I, 2, P)
    allocs <- all.ordwgh$order - 1
    wgh <- all.ordwgh$weights
    if (is.null(n.allocs)) {
      pred.power <- FALSE
      sample.allocs <- allocs
      sample.weight <- wgh
      sample.weight.c <- NULL
      print(paste0("Total number of unique allocations: ", nrow(sample.allocs)))
    } else
    {
      if ((length(unique(K) == 1) == T)) { pred.power <- FALSE} else {
        pred.power <- FALSE
      }
      index <- sample(1:nrow(allocs), n.allocs)
      sample.allocs <- allocs[index, , drop = FALSE]  #maintain matrix type when n.allocs=1
      sample.weight <- wgh[index]
      sample.weight.c <- wgh[-index]
    }
    #if (print.allocs == T) {
    #  print("Evaluated allocations: ")
    #  print(sample.allocs)
    #}

    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list )
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        y_trt <- mu0 + clus.mean + Tx.effect  + error.y
        y_ctr <- mu0 + clus.mean + error.y
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <-  lmer(Y~ Treatment +(1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)), data = analytic.data)
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j]=2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"=pval, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$pval <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }

    ### Binary outcomes
    if (family == "binomial") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean +
          b1
        y0 <- b0 + clus.mean
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "binomial",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }

    ### count outcomes
    if (family == "poisson") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean + b1
        y0 <- b0 + clus.mean
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "poisson",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }
  }

  if (design == "sw") {
    J <- length(P)
    K <- K
    if ((sum(I)/J)%%1 != 0) {
      stop("clusters cannot be evenly distributed into periods")
    }
    if (gen.all == FALSE) {
      pred.power <- FALSE
      sample.allocs <- t(replicate(2, sample(rep(1:J, each = sum(I)/J))))
      index <- split(1:sum(I), rep(1:length(I), I), drop = TRUE)
      all <- list()
      for (n in 1:nrow(sample.allocs)) {
        all[[n]] <- t(sapply(1:length(I), FUN = function(i) {
          freq(sample.allocs[n, index[[i]]], 1:J)
        }))
      }
      temp <- array(unlist(all), dim = c(length(I), J, nrow(sample.allocs)))
      multiplicity <- c()
      for (i in 1:dim(temp)[3]) {
        multiplicity[i] <- multi(I, J, temp, i)
      }
      sample.weight <- multiplicity/sum(multiplicity)  ## weights
    }
    all.ordwgh <- all_allocs(I, J, P)
    allocs <- all.ordwgh$order
    wgh <- all.ordwgh$weights
    if (is.null(n.allocs)) {
      pred.power <- FALSE
      sample.allocs <- allocs
      sample.weight <- wgh
      sample.weight.c <- NULL
      print(paste0("Total number of unique allocations: ", nrow(sample.allocs)))
    } else {
      if ((length(unique(K) == 1) == T)) { pred.power <- FALSE} else {
        pred.power <- TRUE
      }
      index <- sample(1:nrow(allocs), n.allocs)
      sample.allocs <- allocs[index, , drop = FALSE]  #maintain matrix type when n.allocs=1
      sample.weight <- wgh[index]
      unsample.allocs <- allocs[-index, , drop = FALSE]
      sample.weight.c <- wgh[-index]
    }
    #if (print.allocs == T) {
    #  print("Evaluated allocations: ")
    #  print(sample.allocs)
    #}
    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list )
        if (is.null(Time.effect)) {
          Time.effect <- 0
        } else Time.effect <- Time.effect
        clus.id <- rep(1:sum(I), clus.size.list )
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        if (length(Time.effect)>1){
          Time.effect <- c(0, Time.effect) # add 0 to baseline
          y_trt <- mu0 + clus.mean + Tx.effect + Time.effect[period+1] + error.y
          y_ctr <- mu0 + clus.mean + Time.effect[period+1] + error.y
        } else {
          y_trt <- mu0 + clus.mean + Tx.effect + period * Time.effect + error.y
          y_ctr <- mu0 + clus.mean + period * Time.effect + error.y}

        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          cor <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <-  lmer(Y~ factor(Period) + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            } else {
              m <-  lmer(Y~ Period + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            }
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j] <- 2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"= pval,"con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$pval <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)
    }

    ### Binary outcomes
    if (family == "binomial") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list )
        clus.id <- rep(1:sum(I), clus.size.list )
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        if (length(Time.effect)>1){
          b2 <- c(0, b2) # add 0 to baseline
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + b2 * period
          y0 <- b0 + clus.mean + b2 * period
        }
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "binomial", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "binomial", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)
    }

    ### count outcomes
    if (family == "poisson") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list )
        clus.id <- rep(1:sum(I), clus.size.list )
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        if (length(Time.effect)>1){
          b2 <- c(0, b2) # add 0 to baseline
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + period * b2
          y0 <- b0 + clus.mean + period * b2
        }
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "poisson", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "poisson", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
        timecor[i] <- mean(temp[[i]]$cor)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)

    }
  }
  if (pred.power == TRUE) {
      period <- unlist(lapply(1:length(K), FUN= function(i){
        1:K[i] %% (length(P)+1)}))
      Tx <- lapply(1:nrow(unsample.allocs), FUN = function(i) {
        temp.allocs <- rep(unsample.allocs[i, ], K)
      })
      pred.value <- foreach(i = 1:nrow(unsample.allocs), .combine = c) %dorng%
        {
          Tx <- ifelse(period < Tx[[i]], 0, 1)
          pred.TTC <- cor(period, Tx)
          pred.TGI <- abs(sum(K) - sum(Tx) - sum(Tx))
          pred.fit <- glm(cbind(n.sims*power,n.sims*(1-power)) ~ TTC + TGI + I(TTC^2), family = "binomial",
                          data = pred.data)
          pred.temp <- predict(pred.fit, data.frame(TTC = pred.TTC,
                                                    TGI = pred.TGI), type = "response")
        }
    if (plot == TRUE) {
      plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
                             power = c(pred.data$power, pred.value))
      print(ggplot(plotdata, aes(x = power, y = ..density.., weight = weight)) +
              geom_histogram(binwidth = 0.01, fill = "blue") + ggtitle("Power distribution"))
    }
    plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
                           power = c(pred.data$power, pred.value))
    wgt.power <- sum(plotdata$power * plotdata$weight)  # calculate the weighted power
  }
  if (pred.power == FALSE) {
    pred.value <- NULL
    wgt.power <- NULL
    if (plot == TRUE) {
      plotdata <- data.frame(weight = sample.weight/(sum(sample.weight)),
                             power = pred.data$power)
      print(ggplot(plotdata, aes(x = power, y = ..density.., weight = weight)) +
              geom_histogram(binwidth = 0.01, fill = "blue") + ggtitle("Power distribution"))
    }
    wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))
  }
  if (pred.power == TRUE){
    plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
                             power = c(pred.data$power, pred.value))}
  else {
    plotdata <- data.frame(weight = c(sample.weight),
                             power = c(pred.data$power))
    }
  plotdata1 <- plotdata[order(plotdata$power),]
  risk <- sum(plotdata1$weight[which(plotdata1$power <= pwr.thrd)])
  if (is.null(pwr.thrd)){
    risk <- NULL
  }
  list(attained.power.sim = c(pred.data[,1],as.vector(pred.value)), coefs = pred.data[,-1], weights = c(sample.weight, sample.weight.c),
       PREP.sim = wgt.power, risk.sim = risk, allocation = sample.allocs)
}

sim.strata.pd <- function(I, P, S, K, mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL, factor.time = FALSE,
                          design, gen.all = TRUE, n.allocs = NULL, n.sims, rho = NULL,
                          sigma.a = NULL, sigma.e = NULL, print.allocs = TRUE, plot = FALSE, sig.level = 0.05,
                          family = "gaussian") {

  if ((length(unique(K) == 1) == T)) {n.allocs = 1}

  ### check input vadility
  if (sum(I) != length(K)) {
    stop("Number of type of clusters must match the number of cluster sizes")
  }
  if (family != "gaussian" & family != "binomial" & family != "poisson") {
    stop("Family must be one of gaussian, binomial or poisson")
  }
  if (sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) == 0) {
    stop("One of rho, sigma.a, sigma.e must be null")
  }
  ### sigma and rho conversion-----
  if (family == "gaussian") {
    if(sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) !=1)
      stop("Exactly one of rho, sigma.a, sigma.e must be null")

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(sigma.e)) { sigma.e <- sqrt(sigma.a^2 * (1 - rho)/rho) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "binomial"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Exactly one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline risk and treatment effect instead of the input value")}

    mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))

    sigma.e <- sqrt((mu0 * (1 - mu0) + mu1 * (1 - mu1))/2)

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "poisson"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Only one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline rate and treatment effect instead of the input value")}
    mu1 <- Tx.effect * mu0

    sigma.e <- sqrt((mu0 + mu1)/2)

    if (is.null(sigma.a)) {
      sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho))}
    if (is.null(rho)) {
      rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }

  }
  ### -----
  if (design != "sw" & design != "pcrt") {
    stop("Design should be one of sw or pcrt")
  }
  if (length(I) ==1 & length(unique(P)) != 1){
    stop("Please consider function sim.ap as the parameters you entered only have one unique allocation")
  }
  if (gen.all == FALSE & is.null(n.allocs)) {
    stop("Evaluating all allocations is impossible without listing out all allocations")
  }

  J = length(P[[1]])

  if (design == "pcrt" & length(P[[1]]) > 2) {
    stop("Parallel CRT design can only have two groups")
  }

  if(sum(I) != length(S)) stop("the total number of clusters must match the length of S")

  if (length(Time.effect)>1 & length(Time.effect) != (length(P))) {
    stop("length of Time.effect should match the number of transition steps")
  }

  check_strata_sum = unlist(lapply(levels(S), function(name_i){
    length(which(S == name_i)) == sum(P[[which(names(P) == name_i)]])  }))
  if (!is.list(P)) stop("P should be a list for stratified randomization")

  if (!all(check_strata_sum)) stop("the number of clusters within each stratum identified by S must match the sum of clusters in the corresponding row in P")
  if (!all(levels(S) %in% names(P)) | !all(names(P) %in% levels(S))) stop("the strata names/indicators in P must match those in S")

  period_count = unlist(lapply(P, length))
  if (!all(period_count == mean(period_count))) stop("all objects in P must have the same length, i.e. they should specify same number of periods)")

  if (design == "pcrt") {
    TTC <- NULL
    if (gen.all == FALSE) {
      pred.power <- FALSE
      sample.allocs <- t(replicate(n.allocs, sample(0:1, sum(I),
                                                    prob = c(1/2, 1/2), replace = T)))
      index <- split(1:sum(I), rep(1:length(I), I), drop = TRUE)
      all <- list()
      for (n in 1:nrow(sample.allocs)) {
        all[[n]] <- t(sapply(1:length(I), FUN = function(i) {
          freq(sample.allocs[n, index[[i]]], 1:2)
        }))
      }
      temp <- array(unlist(all), dim = c(length(I), 2, nrow(sample.allocs)))
      multiplicity <- c()
      for (i in 1:dim(temp)[3]) {
        multiplicity[i] <- multi(I, 2, temp, i)
      }
      sample.weight <- multiplicity/sum(multiplicity)  ## weights
    }
    all.ordwgh <- all_allocs_strat(I = I, P=P, S = S, K=K)
    allocs <- all.ordwgh$order - 1
    wgh <- all.ordwgh$weights
    if (is.null(n.allocs)) {
      pred.power <- FALSE
      sample.allocs <- allocs
      sample.weight <- wgh
      sample.weight.c <- NULL
      print(paste0("Total number of unique allocations: ", nrow(sample.allocs)))
    } else
    {
      pred.power <- TRUE
      index <- sample(1:nrow(allocs), n.allocs)
      sample.allocs <- allocs[index, , drop = FALSE]  #maintain matrix type when n.allocs=1
      sample.weight <- wgh[index]
      sample.weight.c <- wgh[-index]
    }
    #if (print.allocs == T) {
    #  print("Evaluated allocations: ")
    #  print(sample.allocs)
    #}

    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list )
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        y_trt <- mu0 + clus.mean + Tx.effect  + error.y
        y_ctr <- mu0 + clus.mean + error.y
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <-  lmer(Y~ Treatment +(1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)), data = analytic.data)
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j]=2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"=pval, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$pval <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }

    ### Binary outcomes
    if (family == "binomial") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean +
          b1
        y0 <- b0 + clus.mean
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "binomial",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }

    ### count outcomes
    if (family == "poisson") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean + b1
        y0 <- b0 + clus.mean
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "poisson",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      coef <- c()
      SE <- c()
      con <- c()
      int <- c()
      timecor <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
        coef[i] <- mean(temp[[i]]$coef[,1])
        SE[i] <- mean(temp[[i]]$coef[,2])
        con[i] <- mean(temp[[i]]$con)
        int[i] <- mean(temp[[i]]$int)
      }
      pred.data <- data.frame("power"=pval,"control.size"=con,"trt.size"=int,"SE"=SE,"coefs"=coef)
      wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))  # calculate the weighted power
    }
  }
  if(design == "sw"){
  K <- K
  all.ordwgh <- all_allocs_strat(I, P, S, K)
  allocs <- all.ordwgh$order
  wgh <- all.ordwgh$weights
  if (is.null(n.allocs)) {
    pred.power <- FALSE
    sample.allocs <- allocs
    sample.weight <- wgh
    sample.weight.c <- NULL
    index <- NULL
    print(paste0("Total number of unique allocations: ", nrow(sample.allocs)))
  } else {
    pred.power <- TRUE
    index <- sample(1:nrow(allocs), n.allocs)
    sample.allocs <- allocs[index, , drop = FALSE]  #maintain matrix type when n.allocs=1
    sample.weight <- wgh[index]
    unsample.allocs <- allocs[-index, , drop = FALSE]
    sample.weight.c <- wgh[-index]
  }
  #if (print.allocs == T) {
  #  print("Evaluated allocations: ")
  #  print(sample.allocs)
  #}
  ### Gaussian outcomes
  ### Gaussian outcomes
  if (family == "gaussian") {
    data <- foreach(iterators::icount(n.sims)) %dorng% {
      clus.size.list <- K
      #clus.size.list <- c()
      #for (i in 1:length(K)){
      #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
      #                 prob= c(1-(K[i]-floor(K[i])),
      #                         (K[i]-floor(K[i]))),
      #                 1, replace=T)
      #  clus.size.list <- c(clus.size.list,temp)}
      ### convert cluster size per period as a vector
      period <- unlist(lapply(1:length(K), FUN= function(i){
        1:K[i] %% (length(P[[1]])+1)}))
      size <- sum(clus.size.list )
      if (is.null(Time.effect)) {
        Time.effect <- 0
      } else Time.effect <- Time.effect
      clus.id <- rep(1:sum(I), clus.size.list )
      clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
      error.y = rnorm(sum(clus.size.list), sd = sigma.e)
      if (length(Time.effect)>1){
        Time.effect <- c(0,Time.effect) # add 0 to baseline
        y_trt <- mu0 + clus.mean + Tx.effect + Time.effect[period+1] + error.y
        y_ctr <- mu0 + clus.mean + Time.effect[period+1] + error.y
      }else {
        y_trt <- mu0 + clus.mean + Tx.effect + period * Time.effect + error.y
        y_ctr <- mu0 + clus.mean + period * Time.effect + error.y}
      return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

    }
    temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
      {
        output <- matrix(0, nrow=n.sims,ncol=3)
        con <- c()
        int <- c()
        cor <- c()
        pval <- c()
        for (j in 1:n.sims){
          temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
          Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
          analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
          analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
          analytic.data <- as.data.frame(analytic.data)
          colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                       "Treatment")
          con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
          int[j] <- sum(Tx) # size in intervention across all periods and clusters
          cor[j] <- cor(data[[j]]$period,Tx)
          if (factor.time == TRUE) {
            m <-  lmer(Y~ factor(Period) + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                       data = analytic.data)
          } else {
            m <-  lmer(Y~ Period + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                       data = analytic.data)
          }
          coef <- coef(summary(m))
          output[j,] <- coef(summary(m))["Treatment",]
          pval[j] <- 2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
        }
        list("coef"=output, "pval"= pval,"con"= con, "int"= int, "cor"= cor)
      }
    pval <- c()
    coef <- c()
    SE <- c()
    con <- c()
    int <- c()
    timecor <- c()
    for (i in 1:length(temp)){
      pval[i] <- mean(temp[[i]]$pval <= 0.05)
      coef[i] <- mean(temp[[i]]$coef[,1])
      SE[i] <- mean(temp[[i]]$coef[,2])
      con[i] <- mean(temp[[i]]$con)
      int[i] <- mean(temp[[i]]$int)
      timecor[i] <- mean(temp[[i]]$cor)
    }
    pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)
  }

  ### Binary outcomes
  if (family == "binomial") {
    data <- foreach(iterators::icount(n.sims)) %dorng% {
      clus.size.list <- K
      #clus.size.list <- c()
      #for (i in 1:length(K)){
      #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
      #                 prob= c(1-(K[i]-floor(K[i])),
      #                         (K[i]-floor(K[i]))),
      #                 1, replace=T)
      #  clus.size.list <- c(clus.size.list,temp)}
      ### convert cluster size per period as a vector
      period <- unlist(lapply(1:length(K), FUN= function(i){
        1:K[i] %% (length(P[[1]])+1)}))
      size <- sum(clus.size.list )
      clus.id <- rep(1:sum(I), clus.size.list )
      b0 <- log(mu0/(1 - mu0))
      b1 <- log(Tx.effect)
      if (is.null(Time.effect)) {
        b2 <- 0
      } else b2 <- Time.effect
      clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
      if (length(Time.effect)>1){
        b2 <- c(0, b2) # add 0 to baseline
        y1 <- b0 + clus.mean + b1 + b2[period+1]
        y0 <- b0 + clus.mean + b2[period+1]
      } else {
        y1 <- b0 + clus.mean + b1 + b2 * period
        y0 <- b0 + clus.mean + b2 * period
      }
      p1 <- exp(y1)/(1 + exp(y1))
      p0 <- exp(y0)/(1 + exp(y0))
      y_trt <- rbinom(size, 1, p1)
      y_ctr <- rbinom(size, 1, p0)
      return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

    }
    temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
      {
        output <- matrix(0, nrow=n.sims,ncol=4)
        con <- c()
        int <- c()
        cor <- c()
        for (j in 1:n.sims){
          temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
          Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
          analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
          analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
          analytic.data <- as.data.frame(analytic.data)
          colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                       "Treatment")
          con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
          int[j] <- sum(Tx) # size in intervention across all periods and clusters
          cor[j] <- cor(data[[j]]$period,Tx)
          if (factor.time == TRUE) {
            m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                               (1 | Cluster), family = "binomial", analytic.data,
                             control = glmerControl(optimizer = "bobyqa",
                                                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
          } else {
            m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                             family = "binomial", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                        check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
          }
          output[j,] <- summary(m)$coefficients["Treatment",]
        }
        list("coef"=output, "con"= con, "int"= int, "cor"= cor)
      }
    pval <- c()
    coef <- c()
    SE <- c()
    con <- c()
    int <- c()
    timecor <- c()
    for (i in 1:length(temp)){
      pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
      coef[i] <- mean(temp[[i]]$coef[,1])
      SE[i] <- mean(temp[[i]]$coef[,2])
      con[i] <- mean(temp[[i]]$con)
      int[i] <- mean(temp[[i]]$int)
      timecor[i] <- mean(temp[[i]]$cor)
    }
    pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)
  }

  ### count outcomes
  if (family == "poisson") {
    data <- foreach(iterators::icount(n.sims)) %dorng% {
      clus.size.list <- K
      #clus.size.list <- c()
      #for (i in 1:length(K)){
      #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
      #                 prob= c(1-(K[i]-floor(K[i])),
      #                         (K[i]-floor(K[i]))),
      #                 1, replace=T)
      #  clus.size.list <- c(clus.size.list,temp)}
      ### convert cluster size per period as a vector
      period <- unlist(lapply(1:length(K), FUN= function(i){
        1:K[i] %% (length(P[[1]])+1)}))
      size <- sum(clus.size.list )
      clus.id <- rep(1:sum(I), clus.size.list )
      b0 <- log(mu0)
      b1 <- log(Tx.effect)
      if (is.null(Time.effect)) {
        b2 <- 0
      } else b2 <- Time.effect
      clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
      if (length(Time.effect)>1){
        b2 <- c(0, b2) # add 0 to baseline
        y1 <- b0 + clus.mean + b1 + b2[period+1]
        y0 <- b0 + clus.mean + b2[period+1]
      } else {
        y1 <- b0 + clus.mean + b1 + period * b2
        y0 <- b0 + clus.mean + period * b2
      }
      y_trt <- rpois(size, exp(y1))
      y_ctr <- rpois(size, exp(y0))
      return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

    }
    temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
      {
        output <- matrix(0, nrow=n.sims,ncol=4)
        con <- c()
        int <- c()
        cor <- c()
        for (j in 1:n.sims){
          temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
          Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
          analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
          analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
          analytic.data <- as.data.frame(analytic.data)
          colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                       "Treatment")
          con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
          int[j] <- sum(Tx) # size in intervention across all periods and clusters
          cor[j] <- cor(data[[j]]$period,Tx)
          if (factor.time == TRUE) {
            m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                               (1 | Cluster), family = "poisson", analytic.data,
                             control = glmerControl(optimizer = "bobyqa",
                                                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
          } else {
            m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                             family = "poisson", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                       check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
          }
          output[j,] <- summary(m)$coefficients["Treatment",]
        }
        list("coef"=output, "con"= con, "int"= int, "cor"= cor)
      }
    pval <- c()
    coef <- c()
    SE <- c()
    con <- c()
    int <- c()
    timecor <- c()
    for (i in 1:length(temp)){
      pval[i] <- mean(temp[[i]]$coef[,4] <= 0.05)
      coef[i] <- mean(temp[[i]]$coef[,1])
      SE[i] <- mean(temp[[i]]$coef[,2])
      con[i] <- mean(temp[[i]]$con)
      int[i] <- mean(temp[[i]]$int)
      timecor[i] <- mean(temp[[i]]$cor)
    }
    pred.data <- data.frame("power"=pval,"control.size"=con, "trt.size"=int, "TGI"=abs(con-int),"SE"=SE,"coefs"=coef,"TTC"=timecor)

  }
  }
  if (pred.power == TRUE) {
      period <- unlist(lapply(1:length(K), FUN= function(i){
        1:K[i] %% (length(P)+1)}))
      Tx <- lapply(1:nrow(unsample.allocs), FUN = function(i) {
        temp.allocs <- rep(unsample.allocs[i, ], K)
      })
      pred.value <- foreach(i = 1:nrow(unsample.allocs), .combine = c) %dorng%
        {
          Tx <- ifelse(period < Tx[[i]], 0, 1)
          pred.TTC <- cor(period, Tx)
          pred.TGI <- abs(sum(K) - sum(Tx))
          pred.fit <- glm(cbind(n.sims*power,n.sims*(1-power)) ~ TTC + TGI + I(TTC^2), family = "binomial",
                          data = pred.data)
          pred.temp <- predict(pred.fit, data.frame(TTC = pred.TTC,
                                                    TGI = pred.TGI), type = "response")
        }
    if (plot == TRUE) {
      plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
                             power = c(pred.data$power, pred.value))
      print(ggplot(plotdata, aes(x = power, y = ..density.., weight = weight)) +
              geom_histogram(binwidth = 0.01, fill = "blue") + ggtitle("Power distribution"))
    }
    plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
                           power = c(pred.data$power, pred.value))
    wgt.power <- sum(plotdata$power * plotdata$weight)  # calculate the weighted power
  }
  if (pred.power == FALSE) {
    pred.value <- NULL
    wgt.power <- NULL
    if (plot == TRUE) {
      plotdata <- data.frame(weight = sample.weight/(sum(sample.weight)),
                             power = pred.data$power)
      print(ggplot(plotdata, aes(x = power, y = ..density.., weight = weight)) +
              geom_histogram(binwidth = 0.01, fill = "blue") + ggtitle("Power distribution"))
    }
    plotdata <- data.frame(weight = sample.weight/(sum(sample.weight)),
                           power = pred.data$power)
    wgt.power <- sum(pred.data$power * (sample.weight/(sum(sample.weight))))
  }
  #if (pred.power == TRUE){
  #  plotdata <- data.frame(weight = c(wgh[index], wgh[-index]),
  #                           power = c(pred.data$power, pred.value))
  #  }
  #else {
  #  plotdata <- data.frame(weight = c(sample.weight),
  #                           power = c(pred.data$power))
  #}
  plotdata1 <- plotdata[order(plotdata$power),]
  risk <- sum(plotdata1$weight[which(plotdata1$power <= pwr.thrd)])
  if (is.null(pwr.thrd)){
    risk <- NULL
  }
  list(attained.power.sim = c(pred.data[,1],as.vector(pred.value)), coefs = pred.data[,-1], weights = c(sample.weight,sample.weight.c),
       PREP.sim = wgt.power, risk.sim = risk, allocation = sample.allocs)
}



siglevel <- function(I, P, K, mu0, Tx.effect, Time.effect = NULL, factor.time = FALSE,
                     design, user.allocs, n.sims, rho = NULL, sigma.a = NULL, sigma.e = NULL,
                     sig.level = 0.05, family = "gaussian") {
  ### I: number of clusters of each type (categorized by cluster sizes) J:
  ### number of transitioning period K: cluster size (Length of I must be
  ### the same as lenght of K) mu0: Baseline mean/rate (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes) Tx.effect:
  ### treatment effect (Linear scale for gaussian outcome, natural scale
  ### for binary/count outcomes) Time.effect: time effect (Linear scale for
  ### gaussian outcome, natural scale for binary/count outcomes)
  ### factor.time: T/F whether treat time as continous or categorical
  ### variable design: sw or pcrt user.allocs: The randomization order(s)
  ### to evaluate. Entered as matrix. Each row represent a randomization
  ### order the user try to evaluate. The entries are when the cluster
  ### transits from control to treatment group.  n.sims: the number of
  ### trials the user want to simulate to calculation power rho: ICC
  ### sigma.a: between-cluster variance sigma.e: within-cluster variance
  ### (Only two of rho,sigma.a and sigma.e should be provided) sig.level:
  ### significance level (Default=0.05) family: one of gaussian, binomial
  ### or poisson (case sensitive)


  ### check input vadility
  if (sum(I) != length(K)) {
    stop("Number of type of clusters must match the number of cluster sizes")
  }
  if (family != "gaussian" & family != "binomial" & family != "poisson") {
    stop("Family must be one of gaussian, binomial or poisson")
  }

  if (sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) == 0) {
    stop("One of rho, sigma.a, sigma.e must be null")
  }
  ### sigma and rho conversion-----
  if (family == "gaussian") {
    if(sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) !=1)
      stop("Exactly one of rho, sigma.a, sigma.e must be null")

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(sigma.e)) { sigma.e <- sqrt(sigma.a^2 * (1 - rho)/rho) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "binomial"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Exactly one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline risk and treatment effect instead of the input value")}

    mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))

    sigma.e <- sqrt((mu0 * (1 - mu0) + mu1 * (1 - mu1))/2)

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }
  if (family == "poisson"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Only one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline rate and treatment effect instead of the input value")}
    mu1 <- Tx.effect * mu0

    sigma.e <- sqrt((mu0 + mu1)/2)

    if (is.null(sigma.a)) {
      sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho))}
    if (is.null(rho)) {
      rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }

  }
  ### -----
  if (design != "sw" & design != "pcrt") {
    stop("Design should be one of sw or pcrt")
  }
  if (design == "pcrt" & length(P)>2) {
    stop("Parallel CRT design can only have two groups")
  }

  if (length(sig.level) > nrow(user.allocs)) {
    stop("You entered too many sig.level")
  }

  if (length(Time.effect)>1 & length(Time.effect) != (length(P))) {
    stop("length of Time.effect should match the number of transition steps")
  }

  sample.allocs <- user.allocs
  if (length(sig.level) == 1){
    sig.level <- rep(sig.level,nrow(sample.allocs))
  } else {sig.level <- sig.level}
  ### Parallel CRT
  if (design == "pcrt") {
    TTC <- NULL
    J <- 2
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list )
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        y_trt <- mu0 + clus.mean + Tx.effect  + error.y
        y_ctr <- mu0 + clus.mean + error.y
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <-  lmer(Y~ Treatment +(1| Cluster),  data = analytic.data)
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j]= 2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"=pval, "con"= con, "int"= int)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- quantile(temp[[i]]$pval, sig.level[i])
      }
    }

    ### Binary outcomes
    if (family == "binomial") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean + b1
        y0 <- b0 + clus.mean
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters

            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "binomial",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),)
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
      }
      }

    ### count outcomes
    if (family == "poisson") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        size <- sum(clus.size.list)
        period <- rep(1, size) #not used in pcrt
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        y1 <- b0 + clus.mean + b1
        y0 <- b0 + clus.mean
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- temp.allocs
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            m <- lme4::glmer(Y ~ Treatment + (1 | Cluster), family = "poisson",
                             analytic.data, control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- mean(temp[[i]]$coef[,4] <= sig.level[i])
      }
    }
  }

  ### SW design
  if (design == "sw") {
    J <- length(P)
    K <- K
    ### Gaussian outcomes
    if (family == "gaussian") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list)
        if (is.null(Time.effect)){
          Time.effect <- 0
        } else Time.effect <- Time.effect
        clus.id <- rep(1:sum(I), clus.size.list)
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        error.y = rnorm(sum(clus.size.list), sd = sigma.e)
        if (length(Time.effect)>1){
          Time.effect <- c(0,Time.effect) # add 0 to baseline
          y_trt <- mu0 + clus.mean + Tx.effect + Time.effect[period+1] + error.y
          y_ctr <- mu0 + clus.mean + Time.effect[period+1] + error.y
        } else {
          y_trt <- mu0 + clus.mean + Tx.effect + period * Time.effect + error.y
          y_ctr <- mu0 + clus.mean + period * Time.effect + error.y
        }

        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=3)
          con <- c()
          int <- c()
          cor <- c()
          pval <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <-  lmer(Y~ factor(Period) + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            } else {
              m <-  lmer(Y~ Period + Treatment + (1| Cluster), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),

                         data = analytic.data)
            }
            coef <- coef(summary(m))
            output[j,] <- coef(summary(m))["Treatment",]
            pval[j] <- 2 * pnorm(abs(coef['Treatment', 'Estimate']/coef['Treatment', 'Std. Error']), lower.tail = F)
          }
          list("coef"=output, "pval"= pval, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- quantile(temp[[i]]$pval, sig.level[i])
      }
    }

    ### Binary outcomes
    if (family == "binomial") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list )
        clus.id <- rep(1:sum(I), clus.size.list )
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        if (length(Time.effect)>1){
          b2 <- c(0, b2) # add 0 to baseline
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + b2 * period
          y0 <- b0 + clus.mean + b2 * period
        }
        p1 <- exp(y1)/(1 + exp(y1))
        p0 <- exp(y0)/(1 + exp(y0))
        y_trt <- rbinom(size, 1, p1)
        y_ctr <- rbinom(size, 1, p0)
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims, ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "binomial", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "binomial", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                          check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- quantile(temp[[i]]$coef[,4], sig.level[i])
      }
    }

    ### count outcomes
    if (family == "poisson") {
      data <- foreach(iterators::icount(n.sims)) %dorng% {
        clus.size.list <- K
        #clus.size.list <- c()
        #for (i in 1:length(K)){
        #  temp <- sample(c(floor(K[i]),floor(K[i])+1),
        #                 prob= c(1-(K[i]-floor(K[i])),
        #                         (K[i]-floor(K[i]))),
        #                 1, replace=T)
        #  clus.size.list <- c(clus.size.list,temp)}
        ### convert cluster size per period as a vector
        period <- unlist(lapply(1:length(K), FUN= function(i){
          1:K[i] %% (length(P)+1)}))
        size <- sum(clus.size.list)
        clus.id <- rep(1:sum(I), clus.size.list)
        b0 <- log(mu0)
        b1 <- log(Tx.effect)
        if (is.null(Time.effect)) {
          b2 <- 0
        } else b2 <- Time.effect
        clus.mean <- rep(rnorm(sum(I), sd = sigma.a), clus.size.list)
        if (length(Time.effect)>1){
          b2 <- c(0, b2)
          y1 <- b0 + clus.mean + b1 + b2[period+1]
          y0 <- b0 + clus.mean + b2[period+1]
        } else {
          y1 <- b0 + clus.mean + b1 + period * b2
          y0 <- b0 + clus.mean + period * b2
        }
        y_trt <- rpois(size, exp(y1))
        y_ctr <- rpois(size, exp(y0))
        return(list("clus.id"=clus.id, "clus.size.list"=clus.size.list,"period"=period,"y_trt"=y_trt,"y_ctr"=y_ctr))

      }
      temp <- foreach(i = 1:nrow(sample.allocs), .packages = "lme4") %dorng%
        {
          output <- matrix(0, nrow=n.sims,ncol=4)
          con <- c()
          int <- c()
          cor <- c()
          for (j in 1:n.sims){
            temp.allocs <- rep(sample.allocs[i,],data[[j]]$clus.size.list)
            Tx <- ifelse(data[[j]]$period < temp.allocs, 0, 1)
            analytic.data <- cbind(data[[j]]$clus.id,data[[j]]$period,data[[j]]$y_trt,data[[j]]$y_ctr,Tx)
            analytic.data <- rbind(analytic.data[which(Tx==0),][,-3],analytic.data[which(Tx==1),][,-4])
            analytic.data <- as.data.frame(analytic.data)
            colnames(analytic.data) <- c("Cluster", "Period", "Y",
                                         "Treatment")
            con[j] <- nrow(analytic.data)-sum(Tx) # size in control across all periods and clusters
            int[j] <- sum(Tx) # size in intervention across all periods and clusters
            cor[j] <- cor(data[[j]]$period,Tx)
            if (factor.time == TRUE) {
              m <- lme4::glmer(Y ~ Treatment + factor(Period) +
                                 (1 | Cluster), family = "poisson", analytic.data,
                               control = glmerControl(optimizer = "bobyqa",
                                                      check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            } else {
              m <- lme4::glmer(Y ~ Treatment + Period + (1 | Cluster),
                               family = "poisson", analytic.data, control = glmerControl(optimizer = "bobyqa",
                                                                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
            }
            output[j,] <- summary(m)$coefficients["Treatment",]
          }
          list("coef"=output, "con"= con, "int"= int, "cor"= cor)
        }
      pval <- c()
      for (i in 1:length(temp)){
        pval[i] <- quantile(temp[[i]]$coef[,4], sig.level[i])
      }
    }
  }
  return(pval)
}

### Formula based calculation
analytic.pd <- function(I, P = NULL, K = NULL, S = NULL, user.allocs = NULL, pwr.thrd = NULL, factor.time = FALSE,
                        mu0, Tx.effect, Time.effect = NULL, rho = NULL, family, design, sig.level = 0.05, plot= TRUE,
                        sigma.e = NULL, sigma.a = NULL) {
  # I = vector of cluster types P = # of clusters at each step (excluding
  # baseline) K = vector of cluster sizes factor.time = logical; for
  # whether time should be considered a factor rho = ICC mu0 = baseline
  # risk/rate for binary and count outcomes Tx.effect = treatment effect
  # for continuous outcome, odds ratio for binary and count outcomes
  # Time.effect (for SW): Time effect; linear scale for continuous, natural scale for binomial or poisson If null, default is half of treatment effect on the linear scale (0.5 * Tx.effect for continuous, exp(0.5 * log(Tx.effect)) for binomial and poisson.
  # family should be one of 'gaussian', 'binomial' or 'poisson'
  # design should be one of 'pcrt' or 'sw'
  # sigma.a: between-cluster variance
  # sigma.e: within-cluster variance
  # calculation with design 'pcrt'
  # rep = number of times analytical formulae for 'sw' type should be repeated
  # S = stratum indicator

  if (sum(I) != length(K)) {
    stop("Number of type of clusters must match the number of cluster sizes")
  }
  if (family != "gaussian" & family != "binomial" & family != "poisson") {
    stop("Family must be one of gaussian, binomial or poisson")
  }



  if (design != "sw" & design != "pcrt") {
    stop("Design should be one of sw or pcrt")
  }

  if (design == "pcrt" & if (all(is.null(Time.effect)))
    FALSE else !all(Time.effect == 0)) {
    meesage("Parallel CRT design usually has no time effect")
  }

  rep = 1

  ### old development
  #if (is.null(rep)) {
  #  if(!all(K == floor(K))) stop("'rep' should be specified for non-integer sizes")
  #  else rep = 1
  #} else if(all(K == floor(K)) & rep > 1) warning("'rep' can be set to 0 or 'NULL' for integer cluster sizes for speedier calculations", immediate. = T)

  # sanity checks for stratified randomization
  if (!is.null(S)) {
    strat = TRUE
    J = length(P[[1]])

    if (design == "pcrt" & length(P[[1]]) > 2) {
      stop("Parallel CRT design can only have two groups")
    }

     if(sum(I) != length(S)) stop("the total number of clusters must match the length of S")

     check_strata_sum = unlist(lapply(levels(S), function(name_i){
       length(which(S == name_i)) == sum(P[[which(names(P) == name_i)]])  }))
     if (!is.list(P)) stop("P should be a list for stratified randomization")

     if (!all(check_strata_sum)) stop("the number of clusters within each stratum identified by S must match the sum of clusters in the corresponding row in P")
     if (!all(levels(S) %in% names(P)) | !all(names(P) %in% levels(S))) stop("the strata names/indicators in P must match those in S")

     period_count = unlist(lapply(P, length))
     if (!all(period_count == mean(period_count))) stop("all objects in P must have the same length, i.e. they should the same number of periods)")
  } else {
    strat = FALSE
    J <- length(P)
    }

  # computing any missing rho/sigma.a/sigma.e
  if (family == "gaussian") {
    if(sum(c(is.null(rho), is.null(sigma.a), is.null(sigma.e))) !=1)
      stop("Exactly one of rho, sigma.a, sigma.e must be null")

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(sigma.e)) { sigma.e <- sqrt(sigma.a^2 * (1 - rho)/rho) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }

   if (family == "binomial"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Exactly one of rho and sigma.a must be null")}

     if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline risk and treatment effect instead of the input value")}

    mu1 <- Tx.effect * (mu0/(1 - mu0))/(1 + Tx.effect * (mu0/(1 - mu0)))

    sigma.e <- sqrt((mu0 * (1 - mu0) + mu1 * (1 - mu1))/2)

    if (is.null(sigma.a)) { sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho)) }
    if (is.null(rho)) { rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }
  }

  if (family == "poisson"){

    if (sum(c(is.null(rho), is.null(sigma.a))) !=1) {
      stop("Only one of rho and sigma.a must be null")}

    if(!is.null(sigma.e)){warning("sigma.e will be computed from the baseline rate and treatment effect instead of the input value")}
    mu1 <- Tx.effect * mu0

    sigma.e <- sqrt((mu0 + mu1)/2)

    if (is.null(sigma.a)) {
      sigma.a <- sqrt(sigma.e^2 * rho/(1 - rho))}
    if (is.null(rho)) {
      rho <- (sigma.a^2/(sigma.e^2 + sigma.a^2)) }

  }


  if (design == "pcrt"){
    # allocation-specific power:
    power.rep <- foreach(iterators::icount(rep), .packages = c("arrangements", "Matrix") , .export = c("rand", "multi", "alloc1", "freq", "DesignMatrix", "all_allocs", "all_allocs_strat"), .combine = rbind ) %dorng% {
      # generate design matrices for all unique allocations
      Design_out <- DesignMatrix(I = I, J = J, P = P, K = K, S = S, strat = strat, factor.time = factor.time, user.allocs = user.allocs, design = design)

      XMat = Design_out$XMat
      size = Design_out$size

      if (family == "gaussian") {

        # generate covariance matrix for cross-sectional design
        Ve <- diag(1, nrow = sum(size))

        Vclus <- lapply(1:sum(I), FUN = function(clus) {
          matrix(1, nrow = size[clus], ncol = size[clus])
        })
        Vclus <- bdiag(Vclus)
        VV <- as.matrix(sigma.e^2 * Ve + sigma.a^2 * Vclus)

        # power calculation
        power <- sapply(XMat, FUN = function(XMat_i) {
          VarTx <- solve(t(XMat_i) %*% solve(VV) %*% XMat_i)[1, 1]
          pnorm(Tx.effect/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
            pnorm(-Tx.effect/sqrt(VarTx) - qnorm(1 - sig.level/2),
                  lower.tail = T)
        })

      }
      if (family == "binomial") {
        # only AML has been implemented for binary outcome
        b0 <- log(mu0/(1 - mu0))
        b1 <- log(Tx.effect)

        power <- sapply(XMat, FUN = function(XMat_i) {
          vec_pi <- as.vector(1 + exp(-(XMat_i %*% c(b1, b0))))^(-1)
          clus.id <- split(1:sum(size), rep(1:sum(I), size))

          FIM <- lapply(1:sum(I), FUN = function(ii) {
            VClus <- matrix(1, nrow = size[ii], ncol = size[ii])
            WWi_inv <- solve(diag(vec_pi[clus.id[[ii]]] * (1 - vec_pi[clus.id[[ii]]])))
            VV <- WWi_inv + sigma.a^2 * VClus
            t(XMat_i[clus.id[[ii]], ]) %*% solve(VV) %*% XMat_i[clus.id[[ii]],
            ]
          })
          FIM <- Reduce("+", FIM)

          VarTx <- solve(FIM)[1, 1]

          pnorm(b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
            pnorm(-b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T)
        })
      }

      if (family == "poisson") {
        # only AML has been implemented for count outcome
        b0 <- log(mu0)

        b1 <- log(Tx.effect)

        power <- sapply(XMat, FUN = function(XMat_i) {
          vec_lambda <- as.vector(exp(XMat_i %*% c(b1, b0)))
          clus.id <- split(1:sum(size), rep(1:sum(I), size))

          FIM <- lapply(1:sum(I), FUN = function(ii) {
            VClus <- matrix(1, nrow = size[ii], ncol = size[ii])
            WWi_inv <- solve(diag(vec_lambda[clus.id[[ii]]]))
            VV <- WWi_inv + sigma.a^2 * VClus
            t(XMat_i[clus.id[[ii]], ]) %*% solve(VV) %*% XMat_i[clus.id[[ii]],
            ]
          })
          FIM <- Reduce("+", FIM)

          VarTx <- solve(FIM)[1, 1]

          pnorm(b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
            pnorm(-b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T)
        })
      }
      power
    }

    CV <- sd(K)/mean(K)

    if (family == "gaussian") {
      # Manatunga et al. (2001) Sample size estimation in cluster randomized
      # studies with varying cluster size:
      power.CV <- pnorm(sqrt(sum(I)/2 * mean(K) * Tx.effect^2/2/sigma.e^2/(1 + ((CV^2 + 1) * mean(K) - 1) * rho)) - qnorm(p = (1 - sig.level/2)))
    }

    if (family == "binomial") {

      # Kang et al. (2003) Sample size calculation for dichotomous outcomes
      # in cluster randomization trials with varying cluster size:
      power.CV <- pnorm(sqrt(sum(I)/2 * mean(K) * (mu1 - mu0)^2/(1 + ((CV^2 + 1) * mean(K) - 1) * rho)/(mu0 * (1 - mu0) + mu1 * (1 - mu1))) - qnorm(p = (1 - sig.level/2)))
    }
    if (family == "poisson") {
      # Wang et al. (2018) Sample size calculation for count outcomes in
      # cluster randomization trials with varying cluster sizes:
      power.CV <- pnorm(sqrt(sum(I)/2 * (mu1 - mu0)^2 * mean(K)/(mu1 + mu0)/(1 +  ((CV^2 + 1) * mean(K) - 1) * rho)) - qnorm(p = (1 - sig.level/2)))
    }

  }



  if (design == "sw") {
    if (is.list(P)){np <- length(P[[1]])} else {np <- length(P)}
    K_per <- K / (np+1)
    ncpp <- round(sum(I)/np) # number of cluster per period
    CV <- sd(K_per*(np+1))/mean(K_per*(np+1))
    tol <- sum(K_per*(np+1)) ### total size
    m <- mean(K_per)

    DE_C <- 1+(m*(1+CV^2)-1)*rho
    r <- (m*(1+CV^2)*rho)/DE_C
    DE_R <- 3*np*(1-r)*(1+np*r)/(((np^2)-1)*(2+np*r))
    ## Power.CV is based on Hemming et al 2020 "A tutorial on sample size calculation for multiple-period cluster randomized parallel,
    ## cross-over and stepped-wedge trials using the Shiny CRT Calculator"
    if (family == "gaussian"){
      ses <- abs(Tx.effect)/sigma.e
      power.CV <- pnorm(sqrt(np*m*ncpp*(ses^2)/(4*DE_C*DE_R))-qnorm(1-sig.level/2))
    }
    if (family == "binomial"){
      ses <- abs(mu1-mu0)/sigma.e
      power.CV <- pnorm(sqrt(np*m*ncpp*(ses^2)/(4*DE_C*DE_R))-qnorm(1-sig.level/2))
    }
    if (family == "poisson"){
      ses <- abs(mu1-mu0)/sigma.e
      power.CV <- pnorm(sqrt(np*m*ncpp*(ses^2)/(4*DE_C*DE_R))-qnorm(1-sig.level/2))
    }

      power.rep <- foreach(iterators::icount(rep), .packages = c("arrangements", "Matrix") , .export = c("rand", "multi", "alloc1", "freq", "DesignMatrix", "all_allocs", "all_allocs_strat"), .combine = rbind ) %dorng% {
        # generate design matrices for all unique allocations
        Design_out <- DesignMatrix(I = I, J = J, P = P, K = K, S = S, strat = strat, factor.time = factor.time, user.allocs = user.allocs, design = design)

        XMat = Design_out$XMat
        size = Design_out$size


        if (family == "gaussian") {
          Ve <- diag(1, nrow = sum(size))

          Vclus <- lapply(1:sum(I), FUN = function(clus) {
            matrix(1, nrow = size[clus], ncol = size[clus])
          })
          Vclus <- bdiag(Vclus)
          VV <- as.matrix(sigma.e^2 * Ve + sigma.a^2 * Vclus)

          # power calculation
          power <- sapply(XMat, FUN = function(XMat_i) {
            VarTx <- solve(t(XMat_i) %*% solve(VV) %*% XMat_i)[1, 1]
            pnorm(Tx.effect/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
              pnorm(-Tx.effect/sqrt(VarTx) - qnorm(1 - sig.level/2),
                    lower.tail = T)
          })

        }
        if (family == "binomial") {
          # only AML has been implemented for binary outcome
          b0 <- log(mu0/(1 - mu0)) # baseline risk
          b1 <- log(Tx.effect) # log(OR)

          if (is.null(Time.effect)) # if Time.effect is not specified:
            if (factor.time) b2 <- 0 * b1 * c(1:J) else # when time is categorical, default values are equivalent to continuous time ???
              b2 <- 0 * b1 # default is half of treatment effect when time is continuous

          else { # f Time.effect is specified:
            if (factor.time) {
              if (length(Time.effect) != (J)) stop("length of Time.effect should match the number of transition steps")
            } else {
              if (length(Time.effect) != 1) stop("Time.effect should have length 1 when time is a continuous covariate")
            }
            b2 <- Time.effect
          }

          power <- sapply(XMat, FUN = function(XMat_i) {
            vec_pi <- as.vector(1 + exp(-(XMat_i %*% c(b1, b0, b2))))^(-1)
            clus.id <- split(1:sum(size), rep(1:sum(I), size))

            FIM <- lapply(1:sum(I), FUN = function(ii) {
              VClus <- matrix(1, nrow = size[ii], ncol = size[ii])
              WWi_inv <- solve(diag(vec_pi[clus.id[[ii]]] * (1 - vec_pi[clus.id[[ii]]])))
              VV <- WWi_inv + sigma.a^2 * VClus
              t(XMat_i[clus.id[[ii]], ]) %*% solve(VV) %*% XMat_i[clus.id[[ii]],
                                                                  ]
            })
            FIM <- Reduce("+", FIM)

            VarTx <- solve(FIM)[1, 1]

            pnorm(b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
              pnorm(-b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T)
          })

        }
        if (family == "poisson") {


          # only AML has been implemented for count outcome

          b0 <- log(mu0) # baseline rate

          b1 <- log(Tx.effect) # log(RR)

          if (is.null(Time.effect)) # if Time.effect is not specified:
            if (factor.time) b2 <- 0 * c(1:J) else # when time is categorical, default is ???
              b2 <- 0 * b1 # default is half of treatment effect when time is continuous

          else { # f Time.effect is specified:
            if (factor.time) {
              if (length(Time.effect) != (J)) stop("length of 'Time.effect should match the number of transition steps")
            } else {
              if (length(Time.effect) != 1) stop("Time.effect should have length 1 when time is a continuous covariate")
            }
            b2 <- Time.effect
          }

          power <- sapply(XMat, FUN = function(XMat_i) {
            vec_lambda <- as.vector(exp(XMat_i %*% c(b1, b0, b2)))
            clus.id <- split(1:sum(size), rep(1:sum(I), size))

            FIM <- lapply(1:sum(I), FUN = function(ii) {
              VClus <- matrix(1, nrow = size[ii], ncol = size[ii])
              WWi_inv <- solve(diag(vec_lambda[clus.id[[ii]]]))
              VV <- WWi_inv + sigma.a^2 * VClus
              t(XMat_i[clus.id[[ii]], ]) %*% solve(VV) %*% XMat_i[clus.id[[ii]],
                                                                  ]
            })
            FIM <- Reduce("+", FIM)

            VarTx <- solve(FIM)[1, 1]

            pnorm(b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T) +
              pnorm(-b1/sqrt(VarTx) - qnorm(1 - sig.level/2), lower.tail = T)
          })
        }

        power

      }
  }

  power.rep = matrix(power.rep, nrow=1)

  power.mean = apply(power.rep, 2, mean)
  if (is.null(user.allocs)){

    if (strat){
      all_allocs_out =  all_allocs_strat(I=I, P=P, S=S, K=K)
    } else {
      all_allocs_out = all_allocs(I=I, J=J, P=P)
    }
    wgh = all_allocs_out$weights
    allocs = all_allocs_out$order

    final.data <- data.frame(weight = c(wgh),
                             power = c(power.mean))
    final.data1 <- final.data[order(final.data$power),]
    risk <- sum(final.data1$weight[which(final.data1$power <= pwr.thrd)])
    if (is.null(pwr.thrd)){
      risk <- NULL
    }
    if (plot == TRUE){
      print(ggplot(final.data, aes(x = power, y = ..density.., weight = weight)) +
              geom_histogram(binwidth = 0.01, fill = "blue") + ggtitle("Power distribution (Analytic-based)"))
    }
  } else {
    allocs <- user.allocs
    wgh <- NULL
    risk <- NULL}

  if (design == "sw") return(list(attained.power.analytic=power.mean, PREP.analytic = sum(power.mean*wgh), CV = CV, PREP.CV= power.CV, allocations = allocs, risk.analytic = risk))
  #if (design == "sw" & !is.null(user.allocs)) return(list(attained.power.analytic=power.mean, CV = CV, PREP.CV= power.CV, allocations = allocs, risk.analytic = risk))
  if (design == "pcrt")  return(list(attained.power.analytic=power.mean, PREP.analytic = sum(power.mean*wgh), CV = CV, PREP.CV= power.CV, allocations = allocs, risk.analytic = risk))
  #if (design == "pcrt" & !is.null(user.allocs))  return(list(attained.power.analytic=power.mean, CV = CV, PREP.CV= power.CV, allocations = allocs, risk.analytic = risk))

}

sim.ap <- compiler::cmpfun(sim.ap)
sim.pd <- compiler::cmpfun(sim.pd)
sim.strata.pd <- compiler::cmpfun(sim.strata.pd)
siglevel <- compiler::cmpfun(siglevel)
analytic.pd <- compiler::cmpfun(analytic.pd)


## ---------------------------------------------------
## Main functions and wrapper function to combine
## above functions
## ---------------------------------------------------



power.pd <- function(I, P, K, mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL, factor.time = TRUE,
                     design, rho = NULL, family, sig.level = 0.05,
                     sigma.e = NULL, sigma.a = NULL, method = "analytic",
                     plot = FALSE, gen.all = TRUE, n.allocs,  n.sims, seed = NULL){
  n.cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores=n.cores-1)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (method == "analytic"){
    gen.all = NULL; n.allocs = NULL;  n.sims = NULL
    res <- analytic.pd (I = I, P = P, K = K, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot = plot,
                        design = design, sig.level = sig.level, sigma.e = sigma.e, sigma.a = sigma.a)

  }
  if (method == "sim"){
    res <- sim.pd (I = I, P = P, K = K, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                   design = design, gen.all = gen.all, n.allocs = n.allocs,  n.sims = n.sims, rho = rho,
                   sigma.e = sigma.e, sigma.a = sigma.a, print.allocs = TRUE, plot = plot, sig.level = sig.level,
                   family = family)
  }
  if (method == "both"){

    res <- sim.pd (I = I, P = P, K = K, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                   design = design, gen.all = gen.all, n.allocs = n.allocs,  n.sims = n.sims, rho = rho,
                   sigma.e = sigma.e, sigma.a = sigma.a, print.allocs = TRUE, plot = plot, sig.level = sig.level,
                   family = family)

    res1 <- analytic.pd (I = I, P = P, K = K, user.allocs = res$allocation, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot = plot,
                        design = design, sig.level = sig.level, sigma.e = sigma.e, sigma.a = sigma.a)

    allocs <- res$allocation
    wgh <- res$weights
    final.data <- data.frame(weight = c(wgh),
                             power = c(res1$attained.power))
    final.data1 <- final.data[order(final.data$power),]
    risk <- sum(final.data1$weight[which(final.data1$power <= pwr.thrd)])
    if (is.null(pwr.thrd)){
      risk <- NULL
    }
    res1$risk.analytic <- risk
    res1$allocations <- NULL
    ### Redefine PREP in case only part of it has been evaulated
    res1$PREP.analytic <- sum(res1$attained.power.analytic*(wgh/sum(wgh)))
  }
  if(method =="both") {return(list(results = c(res, res1), inputs = list(I = I, P = P, K = K, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                                                                                                   design = design, gen.all = gen.all, n.allocs = n.allocs,  n.sims = n.sims, rho = rho,
                                                                                                   sigma.e = sigma.e, sigma.a = sigma.a, plot = plot, sig.level = sig.level,
                                                                                                   family = family)))} else {return(list(results = res, inputs = list(I = I, P = P, K = K, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                                                                                                                                                                      design = design, gen.all = gen.all, n.allocs = n.allocs,  n.sims = n.sims, rho = rho,
                                                                                                                                                                      sigma.e = sigma.e, sigma.a = sigma.a, plot = plot, sig.level = sig.level,
                                                                                                                                                                      family = family)))}

}



power.ap <- function(I, P , K , mu0, Tx.effect, Time.effect = NULL, user.allocs , factor.time = TRUE,
                     family, design, rho = NULL, sigma.e = NULL, sigma.a = NULL, sig.level = 0.05,
                     method = "analytic", adj.pwr=FALSE,
                     seed = NULL, n.sims = 1000){
  n.cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores=n.cores-1)
  if (adj.pwr == TRUE){
    if (method == "analytic" | method == "both"){stop("Adjust for Type I error only works when Method = 'sim'")}
    if (!is.null(seed)){
      set.seed(seed)
    }
    if (family != "gaussian"){Tx = 1} else {Tx = 0}
    temp <- siglevel(I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx, Time.effect = Time.effect,
                      factor.time = factor.time, design = design, n.sims = n.sims, rho = rho,
                      sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig.level,
                      family = family)
    sig <- temp
  } else {sig <- sig.level}
  if (method == "analytic"){
    n.sims = NULL
    if (!is.null(seed)){
      set.seed(seed)
    }
    res <- analytic.pd (I = I, P = P, K = K, user.allocs = user.allocs, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot = F,
                        design = design, sig.level = sig, sigma.e = sigma.e, sigma.a = sigma.a)
    res$PREP.CV <- NULL
    res$allocations <- NULL
    res$risk.analytic <- NULL
    res$PREP.analytic <- NULL
  }
  if (method == "sim"){
    if (!is.null(seed)){
      set.seed(seed)
    }
    res <- sim.ap (I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, factor.time = factor.time,
                   design = design, n.sims = n.sims, rho = rho,
                   sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig,
                   family = family)
  }

  if (method == "both"){
    if (!is.null(seed)){
      set.seed(seed)
    }
    res <- sim.ap (I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, factor.time = factor.time,
                   design = design, n.sims = n.sims, rho = rho,
                   sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig,
                   family = family)
    res1 <- analytic.pd (I = I, P = P, K = K, user.allocs = user.allocs, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot = F,
                        design = design, sig.level = sig, sigma.e = sigma.e, sigma.a = sigma.a)
    res1$PREP.CV <- NULL
    res1$allocations <- NULL
    res1$risk.analytic <- NULL
    res1$PREP.analytic <- NULL
  }
  if(method =="both") {return(list(results = c(res,res1), inputs = list(I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, factor.time = factor.time,
                                                                    design = design, n.sims = n.sims, rho = rho,
                                                                    sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig.level,
                                                                    family = family)))} else {return(list(results = res, inputs = list(I = I, P = P, K = K, user.allocs = user.allocs, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, factor.time = factor.time,
                                                                                                                                         design = design, n.sims = n.sims, rho = rho,
                                                                                                                                         sigma.e = sigma.e, sigma.a = sigma.a, sig.level = sig.level,
                                                                                                                                         family = family)))}
}




power.strat.pd <- function(I, P , K , S , mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL,
                          factor.time = TRUE, rho = NULL, family, design, gen.all = TRUE, n.allocs,
                          n.sims, sig.level = 0.05, plot = FALSE, seed = NULL,
                          sigma.e = NULL, sigma.a = NULL, method = "analytic"){
  n.cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores=n.cores-1)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (method == "analytic"){
    gen.all = NULL; n.allocs = NULL; n.sims = NULL
    res <- analytic.pd (I = I, P = P, K = K, S = S, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot =plot,
                        design = design, sig.level = sig.level, sigma.e = sigma.e, sigma.a = sigma.a)
    res$PREP.CV <- NULL
  }
  if (method == "sim"){
    res <- sim.strata.pd (I = I, P = P, K = K, S = S, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                          design = design, gen.all = gen.all, n.allocs = n.allocs, n.sims = n.sims, rho = rho,
                          sigma.e = sigma.e, sigma.a = sigma.a, print.allocs = TRUE, plot = plot, sig.level = sig.level,
                          family = family)
  }
  if (method == "both"){
    res <- sim.strata.pd (I = I, P = P, K = K, S = S, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                          design = design, gen.all = gen.all, n.allocs = n.allocs,  n.sims = n.sims, rho = rho,
                          sigma.e = sigma.e, sigma.a = sigma.a, print.allocs = TRUE, plot = plot, sig.level = sig.level,
                          family = family)
    res1 <- analytic.pd (I = I, P = P, K = K, S = S, user.allocs = res$allocation, pwr.thrd = pwr.thrd, factor.time = factor.time,
                        mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, rho = rho, family = family, plot = plot,
                        design = design, sig.level = sig.level, sigma.e = sigma.e, sigma.a = sigma.a  )
    allocs <- res$allocation
    wgh <- res$weights
    final.data <- data.frame(weight = c(wgh),
                             power = c(res1$attained.power))
    final.data1 <- final.data[order(final.data$power),]
    risk <- sum(final.data1$weight[which(final.data1$power <= pwr.thrd)])
    if (is.null(pwr.thrd)){
      risk <- NULL
    }
    res1$risk.analytic <- risk
    res1$PREP <- NULL
    res1$allocations <- NULL
    res1$PREP.analytic <- sum(res1$attained.power.analytic*(wgh/sum(wgh)))
  }
  if(method =="both") {return(list(results = list(simres = res, analyticres = res1), inputs = list(I = I, P = P, K = K, S = S, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                                                                                                   design = design, gen.all = gen.all, n.allocs = n.allocs, n.sims = n.sims, rho = rho,
                                                                                                   sigma.e = sigma.e, sigma.a = sigma.a, plot = plot, sig.level = sig.level,
                                                                                                   family = family)))} else {return(list(results = res, inputs = list(I = I, P = P, K = K, S = S, mu0 = mu0, Tx.effect = Tx.effect, Time.effect = Time.effect, pwr.thrd = pwr.thrd, factor.time = factor.time,
                                                                                                                                                                      design = design, gen.all = gen.all, n.allocs = n.allocs, n.sims = n.sims, rho = rho,
                                                                                                                                                                      sigma.e = sigma.e, sigma.a = sigma.a, print.allocs = TRUE, plot = plot, sig.level = sig.level,
                                                                                                                                                                      family = family)))}
}



power.ap <- compiler::cmpfun(power.ap)
power.pd <- compiler::cmpfun(power.pd)
power.strat.pd <- compiler::cmpfun(power.strat.pd)
siglevel <- compiler::cmpfun(siglevel)















