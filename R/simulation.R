#' Generate synthetic data for simulation study
#' 
#' @param dimensions 
#' @param w1
#' @param n1
#' @param n2
#' @param n3
#' @param n4
#' @param m
#' @param mu1
#' @param mu2
#' @param mu3
#' @param mu4
#' @param sigma1
#' @param sigma2
#' @param sigma3
#' @param sigma4
#' @param theta
#' @param beta
#' @param seed
#' 
#' @return An object of class \code{l1synthetic} with the following elements:
#' \item{votes}{}
#' \item{legis_data}{}
#' \item{votes_data}{}
#' \item{true_wight}{}
#' \item{true_beta}{}
#' \item{code}{}
#' 
#' @references Sooahn Shin, Yohan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @export generate.data

generate.data <- function(dimensions = 2, 
                          w1 = 0.6,
                          n1 = 40, 
                          n2 = 40, 
                          n3 = 10, 
                          n4 = 10,
                          m = 500,
                          mu1 = c(-0.25,-0.25), 
                          mu2 = c(0.25,0.25),
                          mu3 = c(0.75,-0.75), 
                          mu4 = c(-0.75,0.75), 
                          sigma1 = 0.02, 
                          sigma2 = 0.02, 
                          sigma3 = 0.01, 
                          sigma4 = 0.01, 
                          theta = 1.5, 
                          beta = 10,
                          seed = NULL) {
  if(!dimensions%in%c(1,2)) stop("'dimensions' must be either 1 or 2.\n")
  if(dimensions==1) w1 <- 1
  
  set.seed(seed)
  
  ### 1. Generate ideal points
  n <- n1+n2+n3+n4
  L <- data.frame(name = paste("Legislator",1:n), 
                  group = as.factor(c(rep("group1",n1),rep("group2",n2),rep("group3",n3),rep("group4",n4))),
                  stringsAsFactors = F)
  L[c("true_ideal_point_1d","true_ideal_point_2d")] <- rbind(rmvnorm(n1,mu1,diag(2)*sigma1),
                                                             rmvnorm(n2,mu2,diag(2)*sigma2),
                                                             rmvnorm(n3,mu3,diag(2)*sigma3),
                                                             rmvnorm(n4,mu4,diag(2)*sigma4))
  
  ### 2. Generate Yea & Nay outcome
  O <- data.frame(name = paste("Vote",1:m), 
                  true_yea_position_1d = runif(m,-theta,theta),
                  true_yea_position_2d = runif(m,-theta,theta),
                  true_nay_position_1d = runif(m,-theta,theta),
                  true_nay_position_2d = runif(m,-theta,theta))
  
  ### 3. Compute Pr(Yea) & Pr(Nay)
  Dijy <- L[rep(1:n, each = m),c("true_ideal_point_1d","true_ideal_point_2d")] - 
    do.call("rbind", rep(list(O[c("true_yea_position_1d","true_yea_position_2d")]), n))
  Dijy <- as.matrix(Dijy)%*%matrix(c(w1,0,0,(1-w1)),2,2) # weighted
  Dijy <- abs(Dijy)
  Uijy <- -rowSums(Dijy)*beta
  Dijn <- L[rep(1:n, each = m),c("true_ideal_point_1d","true_ideal_point_2d")] - 
    do.call("rbind", rep(list(O[c("true_nay_position_1d","true_nay_position_2d")]), n))
  Dijn <- as.matrix(Dijn)%*%matrix(c(w1,0,0,(1-w1)),2,2) # weighted
  Dijn <- abs(Dijn)
  Uijn <- -rowSums(Dijn)*beta
  
  Ystar <- pnorm(Uijy-Uijn,0,1)
  Ystar <- matrix(Ystar,n,m,byrow=TRUE)
  
  ### 4. Set Y
  U <- matrix(runif(m*n,0,1),n,m)
  Y <- U<Ystar
  Y <- 1*Y
  rownames(Y) <- L$name
  colnames(Y) <- O$name
  
  if(dimensions==1) {
    L <- L[,c("name","group","true_ideal_point_1d")]
    O <- O[,c("true_yea_position_1d","true_nay_position_1d")]
  }
  
  dat <- list(votes = Y, 
              legis_data = L, 
              votes_data = O, 
              true_weight = c(w1, 1-w1), 
              true_beta = beta,
              code = match.call())
  
  class(dat) <- c("l1synthetic")
  
  return(dat)
}

rmvnorm <- function(n,mu,Sigma){
  E <- matrix(rnorm(n*length(mu)),n,length(mu))
  t( t(E%*%chol(Sigma)) + c(mu) )
}

#' Plot the simulation study
#' 
#' @param synthetic_data
#' @param l1object
#' @param coord_flip
#' 
#' @return 
#' 
#' @import ggplot2
#' @references Sooahn Shin, Yohan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @export plot.simulation

plot.simulation <- function(synthetic_data,
                            l1object,
                            coord_flip = FALSE){
  if(class(synthetic_data)!="l1synthetic") stop("'synthetic_data' is not of class l1synthetic.\n")
  if(class(l1object)!="l1ideal") stop("'l1object' is not of class l1ideal.\n")
  
  dim <- length(l1object$legislators)
  w1 <- synthetic_data$true_weight[1]
  
  if(dim==1){
    ideal_df <- synthetic_data$legis_data
    ideal_df <- ideal_df[!ideal_df$name %in% l1object$legis_dropped,]
    ideal_df$order <- match(ideal_df$true_ideal_point_1d, sort(ideal_df$true_ideal_point_1d, decreasing = T)) 
    
    plot_true <- ggplot() + 
      geom_text(aes(label = name, x = true_ideal_point_1d, y = order, col = group), data = ideal_df) +
      theme_classic() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position="none") +
      labs(x="Dimension 1",y=NULL,title="True Latent Space")
    
    plot_res <- plot.l1ideal(l1object, group = ideal_df$group) + theme(legend.position="none")
    
    ideal_df$ideal_point_1d <- plot_res$layers[[1]]$data[,"ideal_point_1d"]
    
    plot_comp <- ggplot() +
      geom_point(aes(x = true_ideal_point_1d, y = ideal_point_1d, col = group, shape = group), data = ideal_df, alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, col="grey") +
      theme_classic() +
      theme(legend.position="none") +
      labs(x="true",y="estimate",title="True v. Estimate")
    
    p <- list(plot_true, plot_res, plot_comp)
    
  } else {
    ideal_df <- synthetic_data$legis_data
    ideal_df <- ideal_df[!ideal_df$name %in% l1object$legis_dropped,]
    ideal_df$true_ideal_point_1d <- ideal_df$true_ideal_point_1d*w1
    ideal_df$true_ideal_point_2d <- ideal_df$true_ideal_point_2d*(1-w1)
    
    plot_true <- ggplot() + 
      geom_point(aes(x = true_ideal_point_1d, y = true_ideal_point_2d, col = group, shape = group), data = ideal_df, alpha = 0.7) +
      theme_classic() +
      theme(legend.position="none") +
      labs(x="Dimension 1",y="Dimension 2",title="True Latent Space (weighted)")
    
    plot_res <- plot.l1ideal(l1object, group = ideal_df$group, weighted = TRUE) + theme(legend.position="none")
    
    if (!coord_flip) {
      ideal_df$ideal_point_1d <- plot_res$layers[[1]]$data[,"ideal_point_1d"]
      ideal_df$ideal_point_2d <- plot_res$layers[[1]]$data[,"ideal_point_2d"]
    } else {
      plot_res <- plot_res + coord_flip()
      
      ideal_df$ideal_point_1d <- plot_res$layers[[1]]$data[,"ideal_point_2d"]
      ideal_df$ideal_point_2d <- plot_res$layers[[1]]$data[,"ideal_point_1d"]
    }
    
    plot_comp1 <- ggplot() +
      geom_point(aes(x = true_ideal_point_1d, y = ideal_point_1d, col = group, shape = group), data = ideal_df, alpha = 0.7) +
      theme_classic() +
      theme(legend.position="none") +
      labs(x="true",y="estimate",title="True v. Estimate (Dimension 1)")
    
    plot_comp2 <- ggplot() +
      geom_point(aes(x = true_ideal_point_2d, y = ideal_point_2d, col = group, shape = group), data = ideal_df, alpha = 0.7) +
      theme_classic() +
      theme(legend.position="none") +
      labs(x="true",y="estimate",title="True v. Estimate (Dimension 2)")
    
    p <- list(plot_true, plot_res, plot_comp1, plot_comp2)
  }
  return(p)
}

#' Multiple ggplots
#' 
#' @param plotlist
#' @param cols
#' @param layout
#' 
#' @return 
#' 
#' @import grid
#' @references http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @export multiplot

multiplot <- function (..., plotlist = NULL, cols = 1, layout = NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
  }
  if (numPlots == 1) {
    print(plots[[1]])
  }
  else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                               ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                                      layout.pos.col = matchidx$col))
    }
  }
}
