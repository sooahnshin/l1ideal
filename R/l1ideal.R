#' L1 norm based multidimensional ideal point estimation
#'
#' @param rollcall A \code{rollcall} object generated by Simon Jackman's \code{pscl} package.
#' @param dimensions The number of dimensions in the latent space. This function supports either 1 or 2 number of dimensions. The default is \code{2}.
#' @param mcmc The total number of MCMC iterations. The default is \code{5000}.
#' @param thin A positive integer indicating the thining interval for the Markov chain. The default is \code{1}.
#' @param burnin The burnin interval for the Markov chain. The default is \code{0}.
#' @param minvotes Minimum number of votes required for a legislator; i.e., legislators who have voted less than \code{minvotes} are excluded from the analysis. The default is \code{20}.
#' @param lop Minimum proportion of the minority votes; i.e., votes where the minority side is smaller than \code{lop} percent are excluded from the analysis. Recommend to keep the default value \code{0}, since every vote includes valuable information.
#' @param polarity A vector specifying the row index(es) of the legislator(s) in 'votes' matrix constrained to have a positive value on each dimension; e.g., \code{c(10,20)} indicates that the 10th legislator is constrained to have a positive value on the first dimension while the 20th legislator is constrained to have a positive value on the second dimension. The length of the vector must corresponds to the number of dimensions. The default is \code{c(1,1)}.
#' @param verbose A positive integer specified to print the progress on the screen; e.g., every \code{verbose}-th iteration is printed on the screen. The default is \code{100}.
#' @param seed An integer used to set the seed. The default is \code{NULL}.
#'
#' @return An object of class \code{l1ideal} with the following elements:
#' \item{beta}{An object of class \code{mcmc} containing the posterior samples of the beta parameter.}
#' \item{weight}{An object of class \code{mcmc} containing the posterior samples of the weight of the first dimension.}
#' \item{legislators}{A list of \code{mcmc} object(s) containing the posterior samples of the legislator ideal points. Posterior samples of each dimensional coordinates are stored in a separate \code{mcmc}.}
#' \item{yea_positions}{A list of \code{mcmc} object(s) containing the posterior samples of the yea positions. Posterior samples of each dimensional coordinates are stored in a separate \code{mcmc}.}
#' \item{nay_positions}{A list of \code{mcmc} object(s) containing the posterior samples of the nay positions. Posterior samples of each dimensional coordinates are stored in a separate \code{mcmc}.}
#' \item{legis_data}{A matrix or data frame containing covariates of legislators included in the \code{rollcall} object.}
#' \item{votes_data}{A matrix or data frame containing covariates of votes included in the \code{rollcall} object.}
#' \item{legis_dropped}{A vector of excluded legislators.}
#' \item{votes_dropped}{A vector of excluded votes.}
#' \item{code}{The full code that has been implemented.}
#' \item{running_time}{The total running time of the estimation.}
#' 
#' @import coda
#' @references Sooahn Shin, Johan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @useDynLib l1ideal, .registration = TRUE
#' @export l1ideal
#' @examples
#' \dontrun{
#' dat <- generate.data(seed = 1)
#' rc <- pscl::rollcall(dat$votes,
#'                      yea = 1, nay = 0, missing = 2, notInLegis = 3,
#'                      legis.names = dat$legis_data$name, legis.data = dat$legis_data,
#'                      vote.names = dat$votes_data$name, vote.data = dat$votes_data)
#' res <- l1ideal(rc, dimensions = 2, mcmc  = 100, thin = 10, burnin = 100, minvotes = 20,
#'                lop = 0, polarity = c(86,97), verbose = 100, seed = 123)
#' p <- plot.simulation(dat, res)
#' multiplot(plotlist = p, cols = 4)
#' }

l1ideal <- function(rollcall, 
                    dimensions = 2, 
                    mcmc = 5000, 
                    thin = 1, 
                    burnin = 0, 
                    minvotes = 20, 
                    lop = 0, 
                    polarity = c(1,1), 
                    verbose = 100,
                    seed = NULL){
    
    start <- proc.time()
    
    if(class(rollcall)!="rollcall") stop("'rollcall' is not of class rollcall.\n")
    if(!dimensions%in%c(1,2)) stop("'dimensions' must be either 1 or 2.\n")
    if(mcmc<1) stop("'mcmc' cannot be negative.\n") 
    if(thin<1) stop("'thin' cannot be negative.\n")
    if(burnin<0) stop("'burnin' cannot be negative.\n")
    if((lop<0) | (lop>1)) stop("'lop' must fall between 0 and 1.\n")
    if(dimensions != length(polarity)) stop("'polarity' must be a vector of length 'dimensions'\n")
    if(verbose<1) stop("'verbose' cannot be negative.\n")
    
    set.seed(seed)
    votes <- rollcall$votes
    yea_code <- rollcall$codes$yea
    nay_code <- rollcall$codes$nay
    na_code <- c(rollcall$codes$notInLegis,rollcall$codes$missing)
    
    if(any(polarity > nrow(votes))) stop("'polarity' is incorrectly specified")
    
    legis_filtered <- apply(matrix(votes%in%c(yea_code,nay_code),nrow(votes),ncol(votes)),1,sum) >= minvotes
    
    if(!prod(legis_filtered[polarity])) stop("'polarity' fails minimum vote requirements\n")
    
    lop_denominator <- pmin(apply(matrix(votes%in%yea_code,nrow(votes),ncol(votes)),2,sum),apply(matrix(votes%in%nay_code,nrow(votes),ncol(votes)),2,sum))
    lop_numerator <- apply(matrix(votes%in%c(yea_code,nay_code),nrow(votes),ncol(votes)),2,sum)
    votes_lop <- lop_denominator/lop_numerator
    votes_lop[lop_numerator==0] <- 0
    votes_filtered <- votes_lop > lop
    
    nlegis <- sum(legis_filtered)
    nvotes <- sum(votes_filtered)
    
    legis_starts <- runif(dimensions*nlegis,min=-1,max=1)
    votes_starts <- runif(2*dimensions*nvotes,min=-1,max=1)
    
    if (!is.null(rollcall$codes)) {
        votes[votes %in% na_code] <- -1
        votes[votes %in% yea_code] <- 1
        votes[votes %in% nay_code] <- 0
    }
    
    nsamp <- mcmc + burnin
    rcmatrix <- as.vector(votes[legis_filtered,votes_filtered])
    nparam <- 1 + (dimensions*nlegis) + (dimensions*2*nvotes) + 1
    output_len <- (floor(nsamp/thin)+1)*nparam
    
    res <- .C("Cl1ideal",
              rcdata = as.integer(rcmatrix),
              legis_starts = legis_starts,
              bill_starts = votes_starts,
              output = numeric(output_len),
              thin = as.integer(thin),
              ncol = as.integer(nvotes),
              nrow = as.integer(nlegis),
              nsamples = as.integer(nsamp),
              dim = as.integer(dimensions),
              verbose = as.integer(verbose))
    
    l1_mat <- matrix(res$output, ncol = nparam, byrow=TRUE)
    l1_mat <- l1_mat[(burnin/thin+1):nrow(l1_mat),]
    l1_mcmc <- mcmc(l1_mat[-nrow(l1_mat),], start = floor(burnin/thin)*thin+thin, end = ((burnin+mcmc)%/%thin)*thin, thin = thin)
    beta <- l1_mcmc[,1]
    weight <- l1_mcmc[,nparam]
    
    legislators <- vector("list",dimensions)
    for(i in 1:dimensions){
        legislators[[i]] <- l1_mcmc[,(((i-1)*nlegis)+2):(1+(nlegis*i))]
    }
    
    rollcalls <- l1_mcmc[,(2+(nlegis*dimensions)):(1+(nlegis*dimensions)+(2*nvotes*dimensions))]
    
    all_yeas <- rollcalls[,1:(nvotes*dimensions)]
    all_nays <- rollcalls[,(1+(nvotes*dimensions)):(2*nvotes*dimensions)]
    
    rollcall_yeas <- vector("list",dimensions)
    for(i in 1:dimensions){
        rollcall_yeas[[i]] <- all_yeas[,(nvotes*(i-1)+1):(nvotes*i)]
    }
    
    rollcall_nays <- vector("list",dimensions)
    for(i in 1:dimensions){
        rollcall_nays[[i]] <- all_nays[,(nvotes*(i-1)+1):(nvotes*i)]
    }
    
    
    legis_names <- rownames(votes[legis_filtered,])
    votes_names <- colnames(votes[,votes_filtered])
    
    for (i in 1:dimensions){
        colnames(legislators[[i]]) <- legis_names
        colnames(rollcall_yeas[[i]]) <- votes_names
        colnames(rollcall_nays[[i]]) <- votes_names
    }
    
    if(dimensions == 2 & mean(weight) < 0.5) {
        weight <- 1-weight
        legislators <- list(legislators[[2]], legislators[[1]])
        rollcall_yeas <- list(rollcall_yeas[[2]], rollcall_yeas[[1]])
        rollcall_nays <- list(rollcall_nays[[2]], rollcall_nays[[1]])
    }
    
    if(summary(legislators[[1]])[[1]][rownames(votes)[polarity[1]],"Mean"] < 0) {
        legislators[[1]] <- -legislators[[1]]
        rollcall_yeas[[1]] <- -rollcall_yeas[[1]]
        rollcall_nays[[1]] <- -rollcall_nays[[1]]
    }
    
    if(dimensions == 2) {
        if (summary(legislators[[2]])[[1]][rownames(votes)[polarity[2]],"Mean"] < 0) {
            legislators[[2]] <- -legislators[[2]]
            rollcall_yeas[[2]] <- -rollcall_yeas[[2]]
            rollcall_nays[[2]] <- -rollcall_nays[[2]]
        }
    } 
    
    l1res <- list(beta = beta, 
                  weight = weight,
                  legislators = legislators,
                  yea_positions = rollcall_yeas,
                  nay_positions = rollcall_nays,
                  legis_data = rollcall$legis.data,
                  votes_data = rollcall$vote.data,
                  legis_dropped = legis_names[!legis_filtered],
                  votes_dropped = votes_names[!votes_filtered],
                  code = match.call(),
                  running_time = paste0((proc.time() - start)[3]," seconds"))
    
    class(l1res) <- c("l1ideal")
    
    cat("L1 multidimensional ideal point estimation completed.\nTime to completion:", 
        (proc.time() - start)[3], "seconds.\n\n")
    
    l1res
    
}

#' Postprocess the estimates to maximize the posterior
#'
#' @param l1object An object of class \code{l1ideal}
#' 
#' @return An object of class \code{l1ideal} where the estimates for legislator ideal points, yea positions, and nay positions are shifted with the value of \Delta^\ast to maximize the posterior.
#' @references Sooahn Shin, Johan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @export postprocess.l1ideal

postprocess.l1ideal <- function(l1object){
    if(class(l1object)!="l1ideal") stop("'l1object' is not of class l1ideal.\n")
    
    len <- dim(l1object$legislators[[1]])[1]
    
    n <- dim(l1object$legislators[[1]])[2]
    m <- dim(l1object$yea_positions[[1]])[2]
    
    dim <- length(l1object$legislators)
    
    for (s in 1:dim) {
        Delta <- -(rowSums(l1object$legislators[[s]]) + rowSums(l1object$yea_positions[[s]]) + rowSums(l1object$nay_positions[[s]]))/(n+2*m)
        l1object$legislators[[s]] <- l1object$legislators[[s]] - matrix(Delta, nrow = len, ncol = n)
        l1object$yea_positions[[s]] <- l1object$yea_positions[[s]] - matrix(Delta, nrow = len, ncol = m)
        l1object$nay_positions[[s]] <- l1object$nay_positions[[s]] - matrix(Delta, nrow = len, ncol = m)
    }
    
    return(l1object)
}

#' Summarise the l1ideal object
#'
#' @param l1object An object of class \code{l1ideal}
#' 
#' @return A data.frame including the mean, standard deviation, and credible interval of each parameter.
#' @references Sooahn Shin, Johan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @export summary.l1ideal

summary.l1ideal <- function(l1object){
    
    if(class(l1object)!="l1ideal") stop("'l1object' is not of class l1ideal.\n")
    
    dim <- length(l1object$legislators)
    
    df <- data.frame(parameter = as.vector(sapply(1:dim, function(i) c(rep(paste0("ideal_point_",i,"d"), ncol(l1object$legislators[[i]])),
                                                                       rep(paste0("yea_position_",i,"d"), ncol(l1object$yea_positions[[i]])),
                                                                       rep(paste0("nay_position_",i,"d"), ncol(l1object$nay_positions[[i]]))))),
                     name = as.vector(sapply(1:dim, function(i) c(colnames(l1object$legislators[[i]]),
                                                                  colnames(l1object$yea_positions[[i]]),
                                                                  colnames(l1object$nay_positions[[i]])))),
                     mean = as.vector(sapply(1:dim, function(i) c(colMeans(l1object$legislators[[i]]),
                                                                  colMeans(l1object$yea_positions[[i]]),
                                                                  colMeans(l1object$nay_positions[[i]])))),
                     sd = as.vector(sapply(1:dim, function(i) c(apply(l1object$legislators[[i]], 2, sd),
                                                                apply(l1object$yea_positions[[i]], 2, sd),
                                                                apply(l1object$nay_positions[[i]], 2, sd)))),
                     stringsAsFactors = F)
    df[,"2.5%"] <- as.vector(sapply(1:dim, function(i) c(apply(l1object$legislators[[i]], 2, quantile, 0.025),
                                                         apply(l1object$yea_positions[[i]], 2, quantile, 0.025),
                                                         apply(l1object$nay_positions[[i]], 2, quantile, 0.025))))
    df[,"97.5%"] <- as.vector(sapply(1:dim, function(i) c(apply(l1object$legislators[[i]], 2, quantile, 0.975),
                                                          apply(l1object$yea_positions[[i]], 2, quantile, 0.975),
                                                          apply(l1object$nay_positions[[i]], 2, quantile, 0.975))))
    
    df <- rbind(df,
                c("beta", "beta", mean(l1object$beta), sd(l1object$beta), quantile(l1object$beta, c(0.025, 0.975))))
    df <- rbind(df,
                c("weight_1d", "weight_1d", mean(l1object$weight), sd(l1object$weight), quantile(l1object$weight, c(0.025, 0.975))))
    
    df$mean <- as.numeric(df$mean)
    df$sd <- as.numeric(df$sd)
    df$`2.5%` <- as.numeric(df$`2.5%`)
    df$`97.5%` <- as.numeric(df$`97.5%`)
    
    # attr(df, "rownames") <- rownames(df)
    rownames(df) <- 1:nrow(df)
    
    return(df)
}

#' Plot the l1ideal object
#'
#' @param l1object An object of class \code{l1ideal}.
#' @param group A vetor of an attribute of legislators for grouping. Note that the order of the attributes must corresponds to the order of legislators in voting matrix, while excluding the legislators who failed minimum vote requirements specified in \code{l1ideal} function. The default is \code{NULL}.
#' @param weighted Logical. If \code{TRUE}, plot the product of each-dimensional coordinate of ideal point and dimensional weight. Otherwise, plot each-dimensional coordinate of ideal point. Only valid in two dimensional case. The default is \code{TRUE}.
#' 
#' @return A ggplot of ideal points estimated by l1 norm multidimensional ideal point model.
#' 
#' @import ggplot2
#' @references Sooahn Shin, Johan Lim, and Jong Hee Park 2019. "L1 norm Based Multidimensional Ideal Point Estimation: With Application to Roll Call Voting Data" Working Paper.
#' @export plot.l1ideal

plot.l1ideal <- function(l1object,
                         color.group = NULL,
                         shape.group = NULL,
                         weighted = TRUE) {
    
    if(class(l1object)!="l1ideal") stop("'l1object' is not of class l1ideal.\n")
    
    dim <- length(l1object$legislators)
    df <- summary.l1ideal(l1object)
    
    if(!is.null(color.group)) {
        if(length(color.group)>length(df[df$parameter=="ideal_point_1d","name"])) stop("'color.group' is not valid. Check whether it excluded the legislators who failed minimum votes requirement.\n")
        if(length(color.group)<length(df[df$parameter=="ideal_point_1d","name"])) stop("'color.group' is not valid. Check whether it included all the legislators from the analysis.\n")
    }
    
    if(!is.null(shape.group)) {
        if(length(shape.group)>length(df[df$parameter=="ideal_point_1d","name"])) stop("'shape.group' is not valid. Check whether it excluded the legislators who failed minimum votes requirement.\n")
        if(length(shape.group)<length(df[df$parameter=="ideal_point_1d","name"])) stop("'shape.group' is not valid. Check whether it included all the legislators from the analysis.\n")
    }
    
    if(dim==1) {
        ideal_df <- data.frame(ideal_point_1d = df[df$parameter=="ideal_point_1d","mean"],
                               name = df[df$parameter=="ideal_point_1d","name"],
                               stringsAsFactors = F)
        ideal_df$color.group <- color.group
        ideal_df$shape.group <- shape.group
        ideal_df$order <- match(ideal_df$ideal_point_1d, sort(ideal_df$ideal_point_1d, decreasing = T)) 
        plot <- ggplot() + 
            geom_text(aes(label = name, x = ideal_point_1d, y = order, col = color.group, shape = shape.group), data = ideal_df) +
            theme_classic() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()) +
            labs(x="Dimension 1",y=NULL,title="L1 Norm Ideal Point Estimation")
    } else {
        if(weighted) {
            ideal_df <- data.frame(ideal_point_1d = df[df$parameter=="ideal_point_1d","mean"]*df[df$parameter=="weight_1d","mean"],
                                   ideal_point_2d = df[df$parameter=="ideal_point_2d","mean"]*(1-df[df$parameter=="weight_1d","mean"]),
                                   name = df[df$parameter=="ideal_point_1d","name"],
                                   stringsAsFactors = F)
            ideal_df$color.group <- color.group
            ideal_df$shape.group <- shape.group
            plot <- ggplot() + 
                geom_point(aes(x = ideal_point_1d, y = ideal_point_2d, col = color.group, shape = shape.group), data = ideal_df, alpha = 0.7) +
                theme_classic() +
                labs(x="Dimension 1",y="Dimension 2",title="L1 Norm Ideal Point Estimation (weighted)")
        } else {
            ideal_df <- data.frame(ideal_point_1d = df[df$parameter=="ideal_point_1d","mean"],
                                   ideal_point_2d = df[df$parameter=="ideal_point_2d","mean"],
                                   name = df[df$parameter=="ideal_point_1d","name"],
                                   stringsAsFactors = F)
            ideal_df$color.group <- color.group
            ideal_df$shape.group <- shape.group
            plot <- ggplot() + 
                geom_point(aes(x = ideal_point_1d, y = ideal_point_2d, col = color.group, shape = shape.group), data = ideal_df, alpha = 0.7) +
                theme_classic() +
                labs(x="Dimension 1",y="Dimension 2",title="L1 Norm Ideal Point Estimation")
        }
    }
    return(plot)
}

