# Function to get the average number of valid trials per person after cleaning

get.K <- function(dat){
  dat <- subset(dat, accuracy == 1)
  dat <- subset(dat, rt > .2)
  dat <- subset(dat, rt < 2.5)
  dat <- subset(dat, congruency %in% 1:2)
  issue <- issue.sub(dat)
  dat <- subset(dat, !(subject %in% issue))
  I <- length(unique(dat$subject))
  K <- table(dat$subject)/2
  Kall <- round(mean(K, na.rm = T))
  return(Kall)
}


# Function to get the number of participants after cleaning

get.I <- function(dat){
  dat <- subset(dat, accuracy == 1)
  dat <- subset(dat, rt > .2)
  dat <- subset(dat, rt < 2.5)
  dat <- subset(dat, congruency %in% 1:2)
  issue <- issue.sub(dat)
  dat <- subset(dat, !(subject %in% issue))
  return(length(unique(dat$subject)))
}

# Identify participants with fewer than 10 trials per condition for exclusion

issue.sub <- function(dats){
  tmp <- table(dats$subject, dats$congruency)
  tmp <- cbind(tmp, as.numeric(rownames(tmp)))
  return(tmp[tmp[, 1] < 10 | tmp[, 2] < 10, 3])
}


# Estimate a simple linear hierarchical model

genModOneTask <- function(dats, M = 3000, b0 = .030, a0 = 2, b1 = 1, a1 = 2){
  if (mean(dats$congruency %in% 1:2)<1) stop("Conditions must be 1 and 2")
  sub <- as.integer(as.factor(dats$subject))
  I <- max(sub)
  N <- dim(dats)[1]
  keep <- (M/10 + 1) : M
  
  K <- table(dats$subject, dats$congruency)
  Kall <- rowSums(K)
  mn <- tapply(dats$rt, list(dats$subject,dats$congruency), mean)
  sd <- tapply(dats$rt, list(dats$subject,dats$congruency), sd)
  
  theta = alpha = matrix(nrow = M, ncol = I, 0)
  s2 <- 1:M
  muTheta = s2Theta = muAlpha = s2Alpha = 1:M
  theta[1,] <- rep(0, I)
  s2[1] <- .3^2
  muAlpha[1] <- .8
  muTheta[1] <- 0
  muTheta.m <- .05
  muTheta.v <- .1^2
  muAplha.m <- .8
  muAlpha.v <- 1^2
  s2Theta[1] <- .2^2
  s2Alpha[1] <- .3^2
  
  
  x <- matrix(nrow = I, ncol = 2)
  x[,1] <- rep(0,I)
  x[,2] <- rep(1,I)
  
  for (m in 2:M){
    #alpha
    c <- apply(K*(mn - x*theta[m-1,]), 1, sum) / s2[m-1] + muAlpha[m-1] / s2Alpha[m-1]
    v <- 1/(Kall/s2[m-1] + 1/s2Alpha[m-1])
    alpha[m,] <- rnorm(I, c*v, sqrt(v))
    #theta
    c <- K[,2]*(mn[,2] - alpha[m,]) / s2[m-1] + muTheta[m-1] / s2Theta[m-1]
    v <- 1 / (K[,2] / s2[m-1] + 1 / s2Theta[m-1])
    theta[m,] <- rnorm(I, c*v, sqrt(v))
    #s2
    scale <- sum((K-1) * sd^2 + K * (((mn - alpha[m,]) - x*theta[m,])^2)) / 2 + .5
    s2[m] <- rinvgamma(1, shape = (N+1)/2, scale = scale)
    #muAlpha
    v <- 1 / (I / s2Alpha[m-1] + 1 / muAlpha.v)
    c <- sum(alpha[m,]) / s2Alpha[m-1]
    muAlpha[m] <- rnorm(1, v*c, sqrt(v))
    #s2Theta
    scale <- sum((alpha[m,] - muAlpha[m])^2) / 2 + b1^2
    s2Alpha[m] <- rinvgamma(1, shape = I/2 + a1, scale = scale)
    #muTheta
    v <- 1 / (I / s2Theta[m-1] + 1 / muTheta.v)
    c <- sum(theta[m,]) / s2Theta[m-1]
    muTheta[m] <- rnorm(1, v*c, sqrt(v))
    #s2Theta
    scale <- sum((theta[m,] - muTheta[m])^2) / 2 + b0^2
    s2Theta[m] <- rinvgamma(1, shape = I/2 + a0, scale = scale)
  }
  return(list(s2Theta = s2Theta[keep], s2Alpha = s2Alpha[keep]
              , s2 = s2[keep]
              , alpha = alpha[keep,], theta = theta[keep,]))
}

# helper function to estimate split-half reliability using the splithalf package

get.splithalf <- function(dats){
  dats$congruency <- factor(dats$congruency, levels = 1:2, labels = c("congruent", "incongruent"))
  difference <- splithalf(data = dats,
                          outcome = "RT",
                          score = "difference",
                          halftype = "random",
                          permutations = 2000,
                          var.RT = "rt",
                          var.participant = "subject",
                          var.compare = "congruency",
                          compare1 = "congruent",
                          compare2 = "incongruent",
                          average = "mean",
                          plot = TRUE)
  return(c(sh = difference$final_estimates$splithalf
           , sb =  difference$final_estimates$spearmanbrown))
}


# This function generates a plot of the reliability coefficient against trial size.
# Function was kindly provided by Mahbod Mehrvarz.
# Parameters:
#   Var: A Boolean indicating the variant of the plot. Default is TRUE.
#   dots: An optional list of 'l' (trial sizes) and 'g' (gamma values) to plot as points.

makeRelCoefFig = function(Var = T, dots = NULL){
  # Define the reliability function
  rel = function(g2, L) g2 / (g2 + 2/L)
  
  # Set graphical parameters
  par(mgp=c(2.3,.7,0), mar = c(5,5,1,1) + 0.1) 
  
  # Define range for trial sizes
  n1 = 1:9
  n2 = 0:3
  
  # Define colors for the lines in the plot
  colours = c("firebrick3", "firebrick1", "goldenrod1", 
              "darkseagreen2", "mediumseagreen", "deepskyblue2",
              "deepskyblue4", "mediumpurple3", "magenta3"
  )
  
  # Conditional plotting based on the Var parameter
  if (Var == F){
    colours = colours[c(1,3,5,6,8,9)]
    gamma_leg = c(2, 1, .5, .2, .1, .05)
    gamma = gamma_leg^2
    gamma_lab = expression(gamma)
  } else {
    gamma = c(2, 1, .5, .2, .1, .05, .02, .01, .005)
    gamma_leg = gamma
    gamma_lab = expression(gamma^2)
  }
  n = c(as.vector(outer(n1, n2, function(x, y) {x*10^y})), 10000)
  xlims = c(0,4)
  
  # Calculate reliability coefficients
  R = outer(gamma, n, rel)
  
  # Define major ticks for the x-axis
  majors = 0:5
  
  # Start plotting
  plot(NA, NA, 
       type='l', 
       lwd=2, 
       axes=FALSE, 
       xlab="Trial Size (L)", 
       ylab="Reliability Coefficient", 
       ylim = c(-.01,1.1),
       xlim = xlims,
       xaxt = "n",
       yaxt = "n",
       frame.plot = F)
  
  # Add lines for each gamma value
  for(i in 1:length(gamma)) {
    lines(log10(n), R[i,], lwd=2, col=colours[i])
  }
  
  # Add labels and axes
  # mtext("Trial Size (L)", side = 1, line = 3)
  # mtext("Reliability Coefficient", side = 2, line = 3.5)
  axis(2, at = seq(0,1,length.out=6), labels =  seq(0,1,length.out=6), 
       las=1, cex.axis=1.3)
  axis(1, at=majors, labels=10^majors, cex.axis=1.3)
  axis(1, at=log10(n), labels=NA, cex.axis=1.3)
  
  
  # Add horizontal lines for reference
  abline(h=c(.7, .9), lty=2, col = "gray47")
  
  # Add legend
  legend("bottomright", 
         legend=gamma_leg, 
         fill=colours[1:length(gamma)], 
         title=gamma_lab, 
         cex = 1.2,
         bty = "n"
  )
  
  ldots <- length(dots$l)
  # Add points to the plot if 'dots' is provided
  if (!is.null(dots) && length(dots) > 0){
    cols = adjustcolor(1, .5)
    for (i in 1:ldots){
      l = dots$l[i]
      g2 = dots$g[i]^2
      points(log10(l), rel(g2,l), pch = 19, col = cols, cex = 2)
    }
  }
}

## Function to do the splithalf analysis

get.splithalf <- function(dats){
  dats$congruency <- factor(dats$congruency, levels = 1:2, labels = c("congruent", "incongruent"))
  difference <- splithalf(data = dats,
                          outcome = "RT",
                          score = "difference",
                          halftype = "random",
                          permutations = 2000,
                          var.RT = "rt",
                          var.participant = "subject",
                          var.compare = "congruency",
                          compare1 = "congruent",
                          compare2 = "incongruent",
                          average = "mean",
                          plot = FALSE)
  return(c(sh = difference$final_estimates$splithalf, sb =  difference$final_estimates$spearmanbrown))
}
