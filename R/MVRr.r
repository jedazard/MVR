##########################################################################################################################################
# MVR
##########################################################################################################################################

##########################################################################################################################################
# 1. END-USER REGULARIZATION & VARIANCE STABILIZATION FUNCTIONS 
##########################################################################################################################################


##########################################################################################################################################
################
# Usage         :
################
#                   mvr(data, block=rep(1,nrow(data)), tolog=FALSE, nc.min=1, nc.max=30, probs=seq(0, 1, 0.01), B=100, parallel=FALSE, conf=NULL, verbose=TRUE)
#
################
# Description   :
################
#                   Mean-Variance Regularization (MVR) and Variance Stabilization function by similarity statistic
#                   under sample group homoscedasticity or heteroscedasticity assumption. 
#                   Return an object of class "mvr".                 
#                   
################
# Arguments     :
################
# data          :   Numeric matrix of untransformed (raw) data, where samples are by rows and variables (to be clustered) are by columns, 
#                   or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#                   Missing values (NA), NotANumber values (NaN) or Infinite values (Inf) are not allowed.
# block         :   Vector or factor group/blocking variable. 
#                   Must Have length equal to the sample size. 
#                   All group sample sizes must be > 1. 
#                   Defaults to single group situation i.e. to the assumption of equal variance between groups.
# tolog         :   Is the data to be log2-transformed first? Optional, defaults to FALSE.
#                   Negative or null values will be changed to 1 before taking log2-transformation.
# nc.min        :   Minimum # of entertained clusters, defaults to  1
# nc.max        :   Maximum # of entertained clusters, defaults to 30
# probs         :   Numeric vector of probabilities for quantile diagnostic plots. Defaults to seq(0, 1, 0.01).
# B             :   Number of Monte Carlo replicates of the inner loop of the sim statistic function.
# parallel      :   Is parallel computing to be performed? Optional, defaults to FALSE.
# conf          :   List of parameters for cluster configuration. 
#                   Inputs for R package parallel function makeCluster() for cluster setup.
#                   Optional, defaults to NULL. See details for usage.
# verbose       :   Is the output to be verbose? (defaults to TRUE).
#
################
# Values        :
################
#               :   Object of class "mvr", containing the following:
# Xraw          :   Numeric matrix of original data
# Xmvr          :   Numeric matrix of MVR-transformed/standardized data
# centering     :   Numeric vector of centering values for standardization (cluster mean of pooled sample mean)
# scaling       :   Numeric vector of scaling values for standardization (cluster mean of pooled sample std dev)
# MVR           :   List (of size #groups) containing "MeanVarReg" values 
# block         :   Value of argumemt "block"
# tolog         :   Value of argumemt "tolog"
# nc.min        :   Value of argumemt "nc.min"
# nc.max        :   Value of argumemt "nc.max"
# probs         :   Value of argumemt "probs"
#
################

mvr <- function(data, block=rep(1,nrow(data)), tolog=FALSE, nc.min=1, nc.max=30, probs=seq(0, 1, 0.01), B=100, parallel=FALSE, conf=NULL, verbose=TRUE) {
    if (any(table(block) == 1)) 
        stop("All group sample sizes must be greater than 1!")
    p <- ncol(data)
    if (is.null(colnames(data))) colnames(data) <- paste("v", 1:p, sep="")
    block <- as.factor(block)
    lev <- levels(block)
    ng <- nlevels(block)
    def <- vector(mode="list", length=ng)
    tab <- numeric(ng)
    for (g in 1:ng) {
        def[[g]] <- which(block == lev[g])
        tab[g] <- length(def[[g]])
    }
    X.raw <- data
    if (tolog) {
        if (!is.empty(which(data <= 0))) {
            data[data <= 0] <- 1
            warning("Negative or null values have been changed to 1 before taking log2-transformation \n")
        }
        data <- log2(data)
    }
    MVR <- vector(mode="list", length=ng)
    M <- matrix(data=NA, nrow=p, ncol=ng, dimnames=list(colnames(data), lev))
    for (g in 1:ng) {
        if (verbose) cat("Computing optimal cluster configuration of group ", g, " ... \n", sep="")
        MVR[[g]] <- MeanVarReg(data=data[def[[g]], , drop=FALSE], nc.min=nc.min, nc.max=nc.max, probs=probs, B=B, parallel=parallel, conf=conf, verbose=verbose)
        M[,g] <- MVR[[g]]$membership
        if (verbose) cat("group : ", g, " done \n")
    }
    clus <- merging.cluster(M)
    vec.av <- numeric(p)
    vec.sd <- numeric(p)
    for (k in unique(clus)) {
        pma <- pmatch(x=colnames(data), table=names(clus))
        ind <- (clus == k)[pma]
        vec.av[ind] <- mean(pooled.mean(x=data[, ind, drop=FALSE], block=block))
        vec.sd[ind] <- mean(pooled.sd(x=data[, ind, drop=FALSE], block=block))
    }
    X.mvr <- scale(x=data, center=vec.av, scale=vec.sd)    
    return(structure(list(Xraw = X.raw, Xmvr = X.mvr, 
                          centering = vec.av, scaling = vec.sd, 
                          MVR = MVR, block = block, tolog = tolog, 
                          nc.min = nc.min, nc.max = nc.max, probs = probs),
                     class = "mvr"))
}
##########################################################################################################################################





##########################################################################################################################################
# 2. END-USER REGULARIZED TESTS-STATISTICS FUNCTIONS
##########################################################################################################################################


##########################################################################################################################################
################
# Usage         :
################
#                   mvrt.test(data, obj=NULL, block, tolog=FALSE, nc.min=1, nc.max=30, pval=FALSE, 
#                             replace=FALSE, n.resamp=100, parallel=FALSE, conf=NULL, verbose=TRUE)
#
################
# Description   :
################
#                   Computes Mean-Variance Regularized t-test statistic and its significance (p-value)
#                   under sample group homoscedasticity or heteroscedasticity assumption. 
#                   Return an object of class "mvrt.test".                 
#                   
################
# Arguments     :
################
# data          :   Numeric matrix of untransformed (raw) data, where samples are by rows and variables (to be clustered) are by columns, 
#                   or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#                   Missing values (NA), NotANumber values (NaN) or Infinite values (Inf) are not allowed.
# obj           :   mvr object returned by mvr(). Defaults to NULL.
# block         :   Vector or factor group/blocking variable. 
#                   Number of sample groups must be >= 2.
#                   All group sample sizes must be > 1.
# tolog         :   Is the data to be log2-transformed first? Optional, defaults to FALSE.
#                   Negative or null values will be changed to 1 before taking log2-transformation.
# nc.min        :   Minimum # of entertained clusters, defaults to  1
# nc.max        :   Maximum # of entertained clusters, defaults to 30
# pval          :   Shall p-values be computed? If not, n.resamp and replace will be ignored.
#                   If FALSE (default), t-statistic only will be computed,  
#                   If TRUE, exact (permutation test) or approximate (bootstrap test) p-values be computed.
# replace       :   Shall permutation test (default) or bootstrap test be computed?
#                   If FALSE (default), permutation test will be computed with null permutation distribution,
#                   If TRUE, bootstrap test will be computed with null bootstrap distribution.
# n.resamp      :   Number of resamplings (default=100) to compute (by permutation or bootstsrap). 
# parallel      :   Is parallel computing to be performed? Optional, defaults to FALSE.
# conf          :   List of parameters for cluster configuration. 
#                   Inputs for R package parallel function makeCluster() for cluster setup.
#                   Optional, defaults to NULL. See details for usage.
# verbose       :   Is the output to be verbose? (defaults to TRUE).
#
################
# Values        :
################
#               :   Object of class "mvrt.test", containing the following:
# statistic     :   Vector of test-statistic values
# p.value       :   Vector of p-values if requested, otherwise NULL value
#
################

mvrt.test <- function(data, obj=NULL, block, tolog=FALSE, nc.min=1, nc.max=30, pval=FALSE, replace=FALSE, n.resamp=100, parallel=FALSE, conf=NULL, verbose=TRUE) {
    
    statresamp <- function(x, block, nc.min, nc.max, B, replace, verbose) {
        n <- nrow(x) 
        p <- ncol(x) 
        block <- as.factor(block)
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        stat <- matrix(data=NA, nrow=B, ncol=p)
        b <- 1
        while (b <= B) {
            if (verbose) cat("resample: ", b, "\n")
            id <- sample(x=unlist(def), size=n, replace=replace, prob=NULL)
            if (is.valid(x=id, def=def, ng=ng))  {
                stat[b,] <- mvrt(x=x[id,], obj=NULL, lev=lev, tab=tab, ng=ng, def=def, nc.min=nc.min, nc.max=nc.max, parallel=FALSE, conf=NULL, verbose=verbose)
                b <- b + 1
            }
        }
        return(stat)
    }

    if (is.null(obj)) {
        if (any(table(block) == 1)) 
            stop("All group sample sizes must be greater than 1!")
        block <- as.factor(block)
    } else {
        block <- obj$block
        data <- obj$Xraw
        nc.min <- obj$nc.min
        nc.max <- obj$nc.max
        tolog <- obj$tolog
    }
    if (nlevels(block) <= 1) 
        stop("Number of sample groups less than 2. You must have at least two sample groups!")
    p <- ncol(data)
    if (is.null(colnames(data))) colnames(data) <- paste("v", 1:p, sep="")
    if (tolog) {
        if (!is.empty(which(data <= 0))) {
            data[data <= 0] <- 1
            warning("Negative or null values have been changed to 1 before taking log2-transformation \n")
        }
        data <- log2(data)
    }
    lev <- levels(block)
    ng <- nlevels(block)
    def <- vector(mode="list", length=ng)
    tab <- numeric(ng)
    for (g in 1:ng) {
        def[[g]] <- which(block == lev[g])
        tab[g] <- length(def[[g]])
    }
    t.reg <- mvrt(x=data, obj=obj, lev=lev, tab=tab, ng=ng, def=def, nc.min=nc.min, nc.max=nc.max, parallel=parallel, conf=conf, verbose=verbose)
    if (pval == FALSE) {
        p.value <- NULL
    } else {
        if (verbose == TRUE) cat("Computing t-statistics null distributions ...\n")
        if (!parallel) {
            stat.bo <- statresamp(x=data, block=block, nc.min=nc.min, nc.max=nc.max, B=n.resamp, replace=replace, verbose=verbose)
        } else {
            if (conf$type == "SOCK") {
                cl <- makeCluster(spec=conf$names,
                                  type=conf$type,
                                  homogeneous=conf$homo,
                                  outfile=conf$outfile,
                                  verbose=conf$verbose)            
            } else {
                cl <- makeCluster(spec=conf$cpus,
                                  type=conf$type,
                                  homogeneous=conf$homo,
                                  outfile=conf$outfile,
                                  verbose=conf$verbose)
            }
            clusterEvalQ(cl=cl, expr=library("MVR"))
            clusterSetRNGStream(cl=cl)
            stat.cl <- clusterCall(cl=cl, fun=statresamp, x=data, block=block, nc.min=nc.min, nc.max=nc.max, B=ceiling(n.resamp/conf$cpus), replace=replace, verbose=verbose)
            stopCluster(cl)
            stat.bo <- matrix(data=NA, nrow=0, ncol=p)
            for (j in 1:conf$cpus) {
                stat.bo <- rbind(stat.bo, stat.cl[[j]])
            }
            n.resamp <- conf$cpus * ceiling(n.resamp/conf$cpus)
        }
        if (verbose == TRUE) cat("Computation of p-values ...\n")
        p.value <- numeric(p)
        for (j in 1:p) {
            x <- sum(abs(stat.bo[,j]) > abs(t.reg[j]))
            if (replace == FALSE) {
                p.value[j] <- permp(x=x, n1=tab[1], n2=tab[2], nperm=n.resamp)
            } else {
                p.value[j] <- x/n.resamp
            }
        }
    }
    return(structure(list(statistic = t.reg, p.value = p.value),
                     class = "mvrt.test"))
}
##########################################################################################################################################





##########################################################################################################################################
# 3. END-USER DIAGNOSTIC PLOTS FOR QUALITY CONTROL
##########################################################################################################################################


##########################################################################################################################################
################
# Usage         :
################
#                   cluster.diagnostic(obj, title="", span=0.75, degree=2, family="gaussian", device=NULL, file="Cluster Diagnostic Plots")
#
################
# Description   :
################
#                   Function for Plotting Summary Cluster Diagnostic Plots.
#
################
# Arguments     :
################
# obj           :   Object of class "mvr" as returned by mvr().
# title         :   Title of the plot. Defaults to the empty string.
# span          :   Span parameter of the loess() function, which controls the degree of smoothing. Defaults to 0.75.
# degree        :   Degree parameter of the loess() function, which controls the degree of the polynomials to be used. Defaults to 2.
#                   (Normally 1 or 2. Degree 0 is also allowed, but see the ‘Note’ in loess {stats} package.) 
# family        :   Family used for local fittin in {"gaussian", "symmetric"} of the loess() function.
#                   If "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function. 
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Cluster Diagnostic Plots".
#
################
# Values        :
################
#
################

cluster.diagnostic <- function(obj, title="", span=0.75, degree=2, family="gaussian", device=NULL, file="Cluster Diagnostic Plots") {
    
    block <- obj$block
    lev <- levels(block)
    ng <- nlevels(block)

    clusterplot <- function(obj, ng, title, span, degree, family) {
        par(mfcol=c(3, ng), oma=c(0, 0, 3, 0), mar=c(4, 3, 3, 1), mgp=c(2, 0.5, 0), xpd=FALSE)
        
        nc.min <- obj$nc.min
        nc.max <- obj$nc.max
        probs <- obj$probs
        n <- nrow(obj$Xmvr)
        p <- ncol(obj$Xmvr)
        
        X <- matrix(data=rnorm(n=n*p, mean=0, sd=1), ncol=p)
        mu.null       <- pooled.mean(x=X, block=block)
        mu.null.quant <- quantile(x=mu.null, probs=probs)
        sd.null       <- (pooled.sd(x=X, block=block) /sqrt(n-ng)) * sqrt(rchisq(n=p, df=n-ng))
        sd.null.quant <- quantile(x=sd.null, probs=probs)

        for (g in 1:ng) {
            gap <- obj$MVR[[g]]$gap
            sde <- obj$MVR[[g]]$sde
            nc <- obj$MVR[[g]]$nc
            low <- loess(gap ~ as.numeric(nc.min:nc.max), span=span, degree=degree, family=family)$fitted
            plot(nc.min:nc.max, gap, type="n", axes=FALSE, main=paste("Group ", g, sep=""), cex.main=1,
                 xlim=range(nc.min, nc.max), ylim=range(0, max(gap + sde), min(gap - sde)),
                 xlab="Number of clusters", ylab="Sim Statistic")
            lines(nc.min:nc.max, gap, type="b", pch=3, col=1, lty=1)
            lines(low, col=2, lty=5)
            arrows(nc.min:nc.max, gap + sde, nc.min:nc.max, gap - sde, angle=90, code=3, length=0.1, col=1)
            arrows(nc, min(gap) + max(gap)/5, nc, min(gap) + max(gap)/10, length=0.05, angle=20, code=2, col=2, lwd=2)
            text(x=nc, y=min(gap) + max(gap)/5, labels=nc, pos=3, col=2)
            axis(side=1, at=seq(from=nc.min, to=nc.max, by=1), pos=0)
            axis(side=2, at=NULL, pos=nc.min)
            legend("topright", inset=0.05, legend=c("Sim Stat", "Lowess"), col=c(1,2), lty=c(1,5), cex=0.7)

            mu.quant <- obj$MVR[[g]]$mu.quant
            sd.quant <- obj$MVR[[g]]$sd.quant
            w <- which((nc.min:nc.max) <= obj$MVR[[g]]$nc)

            plot(probs, mu.null.quant, type="l", lty=1, col="green", lwd=4, cex.main=1, 
                 main=paste("Group ", g, sep=""), xlab="Percentiles", ylab="Transformed Pooled Means")
            matplot(probs, t(mu.quant), type="l", lty=4, col="black", lwd=1, add=TRUE)
            matplot(probs, t(mu.quant)[,w], type="l", lty=4, col="red", lwd=1, add=TRUE)
            segments(x0=0.5, x1=0.5, y0=min(mu.null.quant, mu.quant), y1=0, col="black", lty=2, lwd=0.5)
            segments(x0=0.5, x1=0, y0=0, y1=0, col="black", lty=2, lwd=0.5)

            plot(probs, sd.null.quant, type="l", lty=1, col="green", lwd=4, cex.main=1, 
                 main=paste("Group ", g, sep=""), xlab="Percentiles", ylab="Transformed Pooled Standard Deviations")
            matplot(probs, t(sd.quant), type="l", lty=4, col="black", lwd=1, add=TRUE)
            matplot(probs, t(sd.quant)[,w], type="l", lty=4, col="red", lwd=1, add=TRUE)
            segments(x0=0.5, x1=0.5, y0=min(sd.null.quant, sd.quant), y1=1, col="black", lty=2, lwd=0.5)
            segments(x0=0.5, x1=0, y0=1, y1=1, col="black", lty=2, lwd=0.5)
        }
        mtext(text=title, cex=1, side=3, outer=TRUE)
    }
        
    if (is.null(device)) {
        clusterplot(obj=obj, ng=ng, title=title, span=span, degree=degree, family=family)
    } else if (device == "PS") {
        postscript(file=paste(getwd(), "/", file, ".ps", sep=""), width=3.5*ng, height=11, onefile=TRUE, horizontal=FALSE)
        clusterplot(obj=obj, ng=ng, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else if (device == "PDF") {
        pdf(file=paste(getwd(), "/", file, ".pdf", sep=""), width=3.5*ng, height=11, onefile=TRUE, paper="US")
        clusterplot(obj=obj, ng=ng, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else {
        stop("Currently allowed display device are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
    
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   target.diagnostic(obj, title="", device=NULL, file="Target Moments Diagnostic Plots")
#
################
# Description   :
################
#                   Function for Plotting Summary Target Moments Diagnostic Plots.
#
################
# Arguments     :
################
# obj           :   Object of class "mvr" as returned by mvr().
# title         :   Title of the plot. Defaults to the empty string.
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Target Moments Diagnostic Plots".
#
################
# Values        :
################
#
################

target.diagnostic <- function(obj, title="", device=NULL, file="Target Moments Diagnostic Plots") {
    
    targetplot <- function(obj, title) {
        par(mfrow=c(2, 3), oma=c(0, 0, 3, 0), mar=c(4, 3, 3, 1), mgp=c(2, 0.5, 0), xpd=FALSE)
        
        block <- obj$block
        lev <- levels(block)
        ng <- nlevels(block)
        n <- nrow(obj$Xmvr)
        p <- ncol(obj$Xmvr)
        data.raw <- obj$Xraw
        data.mvr <- obj$Xmvr
        
        omu <- pooled.mean(x=data.mvr, block=block)
        osd <- pooled.sd(x=data.mvr, block=block)
        X <- matrix(data=rnorm(n=n*p, mean=0, sd=1), ncol=p)
        emu <- pooled.mean(x=X, block=block)
        esd <- (pooled.sd(x=X, block=block) / sqrt(n-ng)) * sqrt(rchisq(n=p, df=n-ng))
       
        if (p > 50) {
            i <- sample(x=1:p, size=50)
        } else {
            i <- 1:p
        }
        p.tt.mu <- t.test(x=omu[i], y=emu[i], alternative="two.sided")$p.value
        p.tt.sd <- t.test(x=osd[i], y=esd[i], alternative="two.sided")$p.value
        p.ks.mu <- ks.test(x=omu, y=emu, alternative="two.sided")$p.value
        p.ks.sd <- ks.test(x=osd, y=esd, alternative="two.sided")$p.value

        dens <- density(pooled.mean(data.raw, block=block), na.rm=TRUE)
        plot(dens, type="l", col=1, lty=1, cex.main=1,
             xlim=range(0, dens$x), ylim=range(0, dens$y),
             xlab="Means", ylab="Density", main="Untransformed")
        abline(v=0, col=2, lty=2)
        abline(v=mean(pooled.mean(data.raw, block=block), na.rm=TRUE), col=1, lty=4)
        legend("topright", inset=0.05, legend="observed", lty=4, col=1, cex=0.7)

        dens.std <- density(pooled.mean(data.mvr, block=block), na.rm=TRUE)
        plot(dens.std, type="l", col=1, lty=1, cex.main=1,
             xlim=range(0, dens.std$x), ylim=range(0, dens.std$y),
             xlab="Means", ylab="Density", main="Mean-Variance Regularization")
        abline(v=0, col=2, lty=2)
        abline(v=mean(pooled.mean(data.mvr, block=block), na.rm=TRUE), col=1, lty=4)
        legend(x="topright", inset=0.05, legend=c("expected", "observed", paste("p=", round(p.tt.mu, 4), sep="")), 
               lty=c(2,4,NA), col=c(2,1,NA), cex=0.7)
    
        qqplot(x=emu, y=omu, col=1, pch=".", cex=3, cex.main=1, 
               main="QQ plot", 
               xlab="Theoretical (Normal) Pooled Mean Quantiles", 
               ylab="Observed Transformed Pooled Mean Quantiles")
        segments(x0=0, x1=0, y0=min(omu), y1=0, col="black", lty=2, lwd=0.5)
        segments(x0=min(emu), x1=0, y0=0, y1=0, col="black", lty=2, lwd=0.5)
        abline(a=0, b=1, lty=1, lwd=0.5, col="red")
        legend(x="topleft", inset=0.05, legend=paste("p=", round(p.ks.mu, 4), sep=""), cex=0.7)

        dens <- density(pooled.sd(data.raw, block=block), na.rm=TRUE)
        plot(dens, type="l", col=1, lty=1, cex.main=1,
             xlim=range(1, dens$x), ylim=range(0, dens$y),
             xlab="Standard Deviations", ylab="Density", main="Untransformed")
        abline(v=1, col=2, lty=2)
        abline(v=mean(pooled.sd(data.raw, block=block), na.rm=TRUE), col=1, lty=4)
        legend(x="topright", inset=0.05, legend="observed", lty=4, col=1, cex=0.7)

        dens.std <- density(pooled.sd(data.mvr, block=block), na.rm=TRUE)
        plot(dens.std, type="l", col=1, lty=1, cex.main=1,
             xlim=range(1, dens.std$x), ylim=range(0, dens.std$y),
             xlab="Standard Deviations", ylab="Density", main="Mean-Variance Regularization")
        abline(v=1, col=2, lty=2)
        abline(v=mean(pooled.sd(data.mvr, block=block), na.rm=TRUE), col=1, lty=4)
        legend(x="topright", inset=0.05, legend=c("expected", "observed", paste("p=", round(p.tt.sd, 4), sep="")), 
               lty=c(2,4,NA), col=c(2,1,NA), cex=0.7)

        qqplot(x=esd, y=osd, col=1, pch=".", cex=3, cex.main=1, 
               main="QQ plot", 
               xlab="Theoretical (Chi) Pooled Std. Dev. Quantiles", 
               ylab="Observed Transformed Pooled Std. Dev. Quantiles")
        segments(x0=1, x1=1, y0=min(osd), y1=1, col="black", lty=2, lwd=0.5)
        segments(x0=min(esd), x1=1, y0=1, y1=1, col="black", lty=2, lwd=0.5)
        abline(a=0, b=1, lty=1, lwd=0.5, col="red")
        legend(x="topleft", inset=0.05, legend=paste("p=", round(p.ks.sd, 4), sep=""), cex=0.7)

        mtext(text=title, cex=1, side=3, outer=TRUE)
    }
        
    if (is.null(device)) {
        targetplot(obj=obj, title=title)
    } else if (device == "PS") {
        postscript(file=paste(getwd(), "/", file, ".ps", sep=""), width=9, height=6.5, onefile=TRUE, horizontal=FALSE)
        targetplot(obj=obj, title=title)
        dev.off()
    } else if (device == "PDF") {
        pdf(file=paste(getwd(), "/", file, ".pdf", sep=""), width=9, height=6.5, onefile=TRUE, paper="US")
        targetplot(obj=obj, title=title)
        dev.off()
    } else {
        stop("Currently allowed display device are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
    
}
##########################################################################################################################################






##########################################################################################################################################
################
# Usage         :
################
#                   stabilization.diagnostic(obj, title="", span=0.5, degree=2, family="gaussian", device=NULL, file="Stabilization Diagnostic Plots")
#
################
# Description   :
################
#                   Function for Plotting Summary Variance Stabilization Diagnostic Plots.
#
################
# Arguments     :
################
# obj           :   Object of class "mvr" as returned by mvr().
# title         :   Title of the plot. Defaults to the empty string.
# span          :   Span parameter of the loess() function, which controls the degree of smoothing. Defaults to 0.5.
# degree        :   Degree parameter of the loess() function, which controls the degree of the polynomials to be used. Defaults to 2.
#                   (Normally 1 or 2. Degree 0 is also allowed, but see the ‘Note’ in loess {stats} package.) 
# family        :   Family used for local fittin in {"gaussian", "symmetric"} of the loess() function.
#                   If "gaussian" fitting is by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function. 
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Stabilization Diagnostic Plots".
#
################
# Values        :
################
#
################

stabilization.diagnostic <- function(obj, title="", span=0.5, degree=2, family="gaussian", device=NULL, file="Stabilization Diagnostic Plots") {
    
    stabplot <- function(obj, title, span, degree, family) {
        par(mfrow=c(1, 2), oma=c(0, 0, 3, 0), mar=c(4, 3, 3, 1), mgp=c(2, 0.5, 0), xpd=FALSE)
        
        block <- obj$block
        data.raw <- obj$Xraw
        data.mvr <- obj$Xmvr

        ord <- order(pooled.mean(x=data.raw, block=block))
        psd <- pooled.sd(x=data.raw, block=block)[ord]
        low <- loess(psd ~ as.numeric(1:ncol(data.raw)), span=span, degree=degree, family=family)$fitted
        plot(x=1:ncol(data.raw), y=psd, type="p", 
             main="Variance-Mean Plot of Untransformed Data", 
             xlab="rank(pooled mean)", ylab="pooled sd", col=2, pch=".", cex=3, cex.main=0.7)
        lines(low, col=1, lty=5, lwd=2)
        
        ord <- order(pooled.mean(x=data.mvr, block=block))
        psd <- pooled.sd(x=data.mvr, block=block)[ord]
        low <- loess(psd ~ as.numeric(1:ncol(data.mvr)), span=span, degree=degree, family=family)$fitted
        plot(x=1:ncol(data.mvr), y=psd, type="p", 
             main="Variance-Mean Plot of MVR-transformed Data", 
             xlab="rank(pooled mean)", ylab="pooled sd", col=2, pch=".", cex=3, cex.main=0.7)
        lines(low, col=1, lty=5, lwd=2)
        
        mtext(text=title, cex=1, side=3, outer=TRUE)
    }
    
    if (is.null(device)) {
        stabplot(obj=obj, title=title, span=span, degree=degree, family=family)
    } else if (device == "PS") {
        postscript(file=paste(getwd(), "/", file, ".ps", sep=""), width=7, height=5, onefile=TRUE, horizontal=FALSE)
        stabplot(obj=obj, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else if (device == "PDF") {
        pdf(file=paste(getwd(), "/", file, ".pdf", sep=""), width=7, height=5, onefile=TRUE, paper="US")
        stabplot(obj=obj, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else {
        stop("Currently allowed display device are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
    
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   normalization.diagnostic(obj, title="", pal, device=NULL, file="Normalization Diagnostic Plots")
#
################
# Description   :
################
#                   Function for Plotting Summary Normalization Diagnostic Plots.
#
################
# Arguments     :
################
# obj           :   Object of class "mvr" as returned by mvr().
# title         :   Title of the plot. Defaults to the empty string.
# pal           :   Color palette.
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Normalization Diagnostic Plots".
#
################
# Values        :
################
#
################

normalization.diagnostic <- function(obj, title="", pal, device=NULL, file="Normalization Diagnostic Plots") {
    
    normplot <- function(obj, title) {
        par(mfrow=c(2, 2), oma=c(0, 0, 3, 0), mar=c(4, 3, 3, 1), mgp=c(2, 0.5, 0), xpd=FALSE)
        
        block <- obj$block
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        n <- nrow(obj$Xmvr)
        p <- ncol(obj$Xmvr)
        data.raw <- obj$Xraw
        data.mvr <- obj$Xmvr
        x.tick <- seq(from=min(1:n), to=max(1:n), by=1)
        y.tick <- seq(from=min(1:p), to=max(1:p), by=min(c(p,100)))

        boxplot(as.data.frame(t(data.raw)), col=rep(2:(ng+1), tab), 
                names=block, notch=FALSE, cex.axis=0.7, cex.main=1, las=2, log="", 
                main="Boxplot of untransformed Data", xlab="Samples", ylab="intensity")
        
        boxplot(as.data.frame(t(data.mvr)), col=rep(2:(ng+1), tab), 
                names=block, notch=FALSE, cex.axis=0.7, cex.main=1, las=2, log="", 
                main="Boxplot of MVR-transformed Data", xlab="Samples", ylab="intensity")

        image(x=1:n, y=1:p, z=data.raw[,rev(1:p)], col=pal, axes=F, asp=NA, xlab="", ylab="")
        title(main="Heatmap of Untransformed Data", cex.main=1, line=2)
        axis(side=1, at=(1:n)[x.tick], labels=colnames(data.raw)[x.tick], pos=-0.5, cex.axis=0.5, las=2)
        axis(side=2, at=(1:p)[p+1-y.tick], labels=as.character((1:p)[y.tick]), pos=(1:n)[1]-0.5, cex.axis=0.5, las=2)
        mtext(text="Variables", side=2, outer=F, cex=0.7, col=1, line=2)   
        mtext(text="Samples", side=1, outer=F, cex=0.7, col=1, line=2)
        box(which="plot")

        image(x=1:n, y=1:p, z=data.mvr[,rev(1:p)], col=pal, axes=F, asp=NA, xlab="", ylab="")
        title(main="Heatmap of MVR-transformed Data", cex.main=1, line=2)
        axis(side=1, at=(1:n)[x.tick], labels=colnames(data.mvr)[x.tick], pos=-0.5, cex.axis=0.5, las=2)
        axis(side=2, at=(1:p)[p+1-y.tick], labels=as.character((1:p)[y.tick]), pos=(1:n)[1]-0.5, cex.axis=0.5, las=2)
        mtext(text="Variables", side=2, outer=F, cex=0.7, col=1, line=2)   
        mtext(text="Samples", side=1, outer=F, cex=0.7, col=1, line=2)
        box(which="plot")

        mtext(text=title, cex=1, side=3, outer=TRUE)
    }
    
    if (is.null(device)) {
        normplot(obj=obj, title=title)
    } else if (device == "PS") {
        postscript(file=paste(getwd(), "/", file, ".ps", sep=""), width=7, height=8, onefile=TRUE, horizontal=FALSE)
        normplot(obj=obj, title=title)
        dev.off()
    } else if (device == "PDF") {
        pdf(file=paste(getwd(), "/", file, ".pdf", sep=""), width=7, height=8, onefile=TRUE, paper="US")
        normplot(obj=obj, title=title)
        dev.off()
    } else {
        stop("Currently allowed display device are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
    
}
##########################################################################################################################################





##########################################################################################################################################
# 4. INTERNAL SUBROUTINES 
# (not to be called by end-user)
##########################################################################################################################################


##########################################################################################################################################
################
# Usage         :   is.empty(x) 
################
#
################
# Description   :
################
#                   Internal function called by mvr() and mvrt.test() to represent the empty 
#                   array, matrix, or vector of zero dimension or length. 
#                   Often returned by expressions and functions whose value is undefined. 
#                   
################
# Arguments     :
################
# x             :   Array, matrix or vector of any type
#
################
# Values        : 
################
#               :   Logical. Returns TRUE if its argument is empty and FALSE otherwise.
#
################

is.empty <- function(x) {
    if (is.vector(x)) {
        if((length(x) == 0) || (x == "")) {
            return(T)
        } else {
            return(F)
        }
    } else if (is.matrix(x) || is.data.frame(x)) {
        return( ((nrow(x) == 0) || (ncol(x) == 0)) )
    } else {
        return( ((length(x) == 0) || (x == "")) )
    }
}
##########################################################################################################################################





########################################################################################################################################## 
################
# Usage         :   is.valid(x, def, ng) 
################
#
################
# Description   :
################
#                   Internal function called by mvrt.test() to validate the bootstrap set of indices if a bootstrap test 
#                   is undertaken in the mvrt.test() function.
#                   Check that indices are unique per sample group. 
#                   
################
# Arguments     :
################
# x             :   Vector of bootstrap set of indices 
# ng            :   Number of sample groups in original data
# def           :   List of samples indices per sample group in original data
#
################
# Values        : 
################
# ok            :   Logical, returns TRUE if valid and FALSE otherwise.
#
################

is.valid <- function(x, def, ng) {
    ok <- TRUE
    for (g in 1:ng) {
        ok <- (ok && (length(unique(x[def[[g]]])) > 1))
    }
    return(ok)
}        
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   pooled.sd(x, block)   
#
################
# Description   :
################
#                   Pooled group sample standard deviation                  
#                   
################
# Arguments     :
################
# x             :   Numeric vector or matrix with variables in columns
# block         :   Vector or factor grouping/blocking variable. 
#                   Must Have length equal to the sample size. 
#                   All group sample sizes must be > 1.
#
################
# Values        :
################
#               :   Vector of size the dimensionality where each entry is the variable-wise pooled group sample standard deviation
#
################

pooled.sd <- function(x, block) {
    if (!is.numeric(x) && !is(x, "vector") && !is(x, "matrix") && !is(x, "data.frame")) 
        stop("The object should be of class numeric, vector, matrix or data.frame! \n")
    if (is(x, "data.frame")) {
        x <- as.matrix(x)
        mode(x) <- "numeric"
    }
    if (any(table(block) == 1))
        stop("All group sample sizes must be greater than 1! \n")
    if (length(table(block)) == 1) {
        if (!is.matrix(x)) {
            return(sd(x, na.rm=TRUE))
        } else {
            return(apply(x, 2, sd, na.rm=TRUE))
        }
    } else {
        block <- as.factor(block)
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        psd <- function(x, tab, n, ng) { sqrt(sum((tab-1)*x^2)/(n - ng)) }
        if (!is.matrix(x)) {
            n <- length(x)
            sdev.gp <- numeric(ng)
            for (g in 1:ng) {
                sdev.gp[g] <- sd(x[def[[g]]], na.rm=TRUE)
            }        
            return(psd(x=sdev.gp, tab=tab, n=n, ng=ng))
        } else {
            n <- nrow(x)
            p <- ncol(x)        
            sdev.gp <- matrix(data=NA, nrow=ng, ncol=p)
            for (g in 1:ng) {
                sdev.gp[g,] <- apply(x[def[[g]], , drop=FALSE], 2, sd, na.rm=TRUE)
            }        
            return(apply(sdev.gp, 2, function(x) psd(x, tab, n, ng)))
        }
    }
}
##########################################################################################################################################



  

##########################################################################################################################################
################
# Usage         :
################
#                   pooled.mean(x, block)
#
################
# Description   :
################
#                   Pooled group sample mean                 
#                   
################
# Arguments     :
################
# x             :   Numeric vector or matrix with variables in columns
# block         :   Vector or factor grouping/blocking variable. 
#                   Must Have length equal to the sample size.
#                   All group sample sizes must be > 1.
#
################
# Values        :
################
#               :   Vector of size the dimensionality where each entry is the variable-wise pooled group sample mean
#
################

pooled.mean <- function(x, block) {
    if (!is.numeric(x) && !is(x, "vector") && !is(x, "matrix") && !is(x, "data.frame")) 
        stop("The object should be of class numeric, vector, matrix or data.frame! \n")
    if (is(x, "data.frame")) {
        x <- as.matrix(x)
        mode(x) <- "numeric"
    }
    if (length(table(block)) == 1) {
        if (!is.matrix(x)) {
            return(mean(x, na.rm=TRUE))
        } else {
            return(apply(x, 2, mean, na.rm=TRUE))
        }
    } else {
        block <- as.factor(block)
        lev <- levels(block)
        ng <- nlevels(block)
        def <- vector(mode="list", length=ng)
        tab <- numeric(ng)
        for (g in 1:ng) {
            def[[g]] <- which(block == lev[g])
            tab[g] <- length(def[[g]])
        }
        pmean <- function(x, tab, n) { sum(tab*x)/n }
        if (!is.matrix(x)) {
            n <- length(x)
            mean.gp <- numeric(ng)
            for (g in 1:ng) {
                mean.gp[g] <- mean(x[def[[g]]], na.rm=TRUE)
            }        
            return(pmean(x=mean.gp, tab=tab, n=n))
        } else {
            n <- nrow(x)
            p <- ncol(x)        
            mean.gp <- matrix(data=NA, nrow=ng, ncol=p)
            for (g in 1:ng) {
                mean.gp[g,] <- apply(x[def[[g]], , drop=FALSE], 2, mean, na.rm=TRUE)
            }        
            return(apply(mean.gp, 2, function(x) pmean(x, tab, n)))
        }
    }
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   MeanVarReg(data, nc.min, nc.max, probs, B, parallel, conf, verbose)  
#
################
# Description   :
################
#                   Core subroutine of mvr() and mvrt.test() functions.
#                   Returns optimal cluster configuration and quantiles of means and standard deviations 
#                   for each cluter configuration when needed.
#
################
# Arguments     :
################
# data          :   Numeric matrix of untransformed (raw) data, where samples are by rows and variables (to be clustered) are by columns.
# nc.min        :   Minimum number of clusters
# nc.max        :   Maximum number of clusters
# probs         :   Numeric vector of probabilities for quantile diagnostic plots.
# B             :   Number of Monte Carlo replicates of the inner loop of the sim statistic function
# parallel      :   Is parallel computing to be performed? Optional, defaults to FALSE.
# conf          :   List of parameters for cluster configuration. 
#                   Inputs for R package parallel function makeCluster() for cluster setup.
#                   Optional, defaults to NULL. See details for usage.
# verbose       :   Is the output to be verbose?
#
################
# Values        :
################
# membership    :   Vector of cluster membership for each variable
# nc            :   Number of clusters found in optimal cluster configuration
# gap           :   Similarity statistic values
# sde           :   Standard errors of the similarity statistic values
# mu.quant      :   Numeric matrix (nc.max - nc.min + 1) x (length(probs)) of quantiles of means.
# sd.quant      :   Numeric matrix (nc.max - nc.min + 1) x (length(probs)) of quantiles of standard deviations.
#
################

MeanVarReg <- function(data, nc.min, nc.max, probs, B, parallel, conf, verbose) {
    data <- t(data)                             # transpose matrix for compatibility with clustering functions
    p <- nrow(data)                             # variables are to be clustered, so they are now by rows !
    n <- ncol(data)
    aver <- apply(data, 1, mean)
    sdev <- apply(data, 1, sd)
    mu.std <- matrix(data=NA, nrow=nc.max - nc.min + 1, ncol=p)
    sd.std <- matrix(data=NA, nrow=nc.max - nc.min + 1, ncol=p)
    gap <- rep(NA, nc.max - nc.min + 1)
    sde <- rep(NA, nc.max - nc.min + 1)
    for (k in nc.min:nc.max) {
        clus <- km.clustering(data=cbind(aver, sdev), k=k, ns=100, maxiter=1000)
        vec.av <- numeric(p)
        vec.sd <- numeric(p)
        for (i in 1:k) {
            wi <- which(clus$membership == i)
            vec.av[wi] <- clus$centers[i,1] 
            vec.sd[wi] <- clus$centers[i,2]
        }
        data.std <- as.matrix((data - (vec.av %*% t(rep(1, n)))) / (vec.sd %*% t(rep(1, n))))
        gp <- sim.dis(data=data.std, k=k, B=B, parallel=parallel, conf=conf, verbose=verbose)
        gap[k] <- gp$gapk
        sde[k] <- gp$sk
        mu.std[k,] <- apply(data.std, 1, mean)
        sd.std[k,] <- apply(data.std, 1, sd)
        if (verbose) cat("number of clusters: ", k, " done\n")
    }
    nc <- which.min(gap)
    nc <- max(which(gap <= gap[nc] + sde[nc]))
    clus <- km.clustering(data=cbind(aver, sdev), k=nc, ns=100, maxiter=1000)
    if (!is.null(probs)) {
        mu.quant <- apply(mu.std, 1, function(x) quantile(x=x, probs=probs))
        sd.quant <- apply(sd.std, 1, function(x) quantile(x=x, probs=probs))
    } else {
        mu.quant <- numeric(0)
        sd.quant <- numeric(0)
    }
    return(list(membership=clus$membership, nc=nc, gap=gap, sde=sde, mu.quant=t(mu.quant), sd.quant=t(sd.quant)))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   km.clustering(data, k, ns, maxiter)
#
################
# Description   :
################
#                   Wrapper Subroutine Around C Subroutine for 'K-means' Clustering Algorithm.
#                   Perform "K-means" clustering on a data matrix.               
#                   Internal subroutine of sim.dis(), and MeanVarReg().
#                   
################
# Arguments     :
################
# data          :   Numeric matrix of data where variables are by columns and samples to cluster are by rows.
# k             :   Integer of number of centroids to start with.
# ns            :   How many random seedings should be done.
# maxiter       :   Integer of maximum number of iterations in 'K-means' clustering algorithm (convergence criteria).
#
################
# Values        :
################
# lWk           :   Log-transformed within cluster dispersion statistic
# centers       :   Cluster centers
# membership    :   Cluster membership of each observation
# obj           :   object of class "kmeans" as returned by kmeans() 
#
################

km.clustering <- function(data, k, ns, maxiter) {
    data <- as.matrix(data)
    cn <- unique(data)
    n <- nrow(data)
    p <- ncol(data)
    nn <- nrow(cn)
    
    Z <- .C("MVR_km_clustering",  as.double(data),
                                  as.double(cn),
                                
                                  centers=double(k*p),
                                  cl=integer(n),
                                  nc=integer(k),
                                  wss=double(k),
                                  tot.wss=as.double(0.0),
                                  err=as.integer(0),
                                
                                  as.integer(n),
                                  as.integer(nn),
                                  as.integer(p),
                                  as.integer(k),
                                  as.integer(ns),
                                  as.integer(maxiter),
                                  
                                  NAOK=FALSE, 
                                  DUP=TRUE, 
                                  PACKAGE="MVR")
    
    if (Z$err > 0) {
        warning("K-means clustering did not converge in ", maxiter, " iterations.", call.=FALSE)
    }
    
    cent <- matrix(Z$centers, k)
    dimnames(cent) <- list(1L:k, dimnames(data)[[2L]])
    
    memb <- Z$cl
    if(!is.null(rn <- rownames(data)))
        names(memb) <- rn
    
    totss <- sum(scale(data, scale=FALSE)^2)
    
    obj = structure(list(cluster = memb, centers = cent, totss = totss,
                         withinss = Z$wss, tot.withinss = Z$tot.wss,
                         betweenss = totss - Z$tot.wss, size = Z$nc),
                    class = "kmeans")
    
    return(list(lWk=log(Z$tot.wss), centers=cent, membership=memb, obj=obj))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   withinsumsq(n, p, B, k)      
#
################
# Description   :
################
#                   Within-Cluster Sum of Squares Distances Subroutine.
#                   Wrapping subroutine around internal C subroutine for returning Monte-Carlo replicates 
#                   of within-cluster sum of squares distances for a given number of clusters and
#                   under a standard Gaussian reference distribution.
#                   Internal subroutine of sim.dis().
#                                          
################
# Arguments     :
################
# n             :   Number of rows of data matrix (where points to cluster are by rows (usually samples)).
# p             :   Number of coluns of data matrix (where variables are by column).
# k             :   Fixed number of clusters
# B             :   Number of Monte Carlo replicates
#
################
# Values        : 
################
# lWk.mc        :   Numeric B-vector of Monte Carlo replicates of within-cluster sum of squares.
#
################

withinsumsq <- function(n, p, B, k) {
    W <- .C("MVR_withinsumsq", as.integer(n),
                               as.integer(p),
                               as.integer(k),
                               as.integer(B),
                               lWk=double(B),
                               nstart=as.integer(10),
                               maxiter=as.integer(1000),
                               err=as.integer(0),
                               
                               NAOK=FALSE, 
                               DUP=TRUE, 
                               PACKAGE="MVR")
    
    if (W$err > 0) {
        warning("Some k-means clusterings for simulated data did not converge in 1000 iterations.", call.=FALSE)
    }
    return(W$lWk)
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   sim.dis(data, k, B, parallel, conf, verbose)      
#
################
# Description   :
################
#                   Internal similarity statistic function called by MeanVarReg() 
#                   for estimating similarity statistic with sampling variability 
#                   under a standard Gaussian reference distibution. 
#                                          
################
# Arguments     :
################
# data          :   Standardized numeric matrix of data, where points to cluster are by rows (usually samples), 
#                   or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#                   Missing values are not allowed. NaN or Inf are neither allowed.
# k             :   Fixed number of clusters
# B             :   Number of Monte Carlo replicates of the inner loop of the gap statistic function.
# parallel      :   Is parallel computing to be performed? Optional, defaults to FALSE.
# conf          :   List of parameters for cluster configuration. 
#                   Inputs for R package parallel function makeCluster() for cluster setup.
#                   Optional, defaults to NULL. See details for usage.
# verbose       :   Is the output to be verbose?
#
################
# Values        : 
################
# gapk          :   Similarity statistic value
# sk            :   Standard error of the similarity statistic value
#
################

sim.dis <- function(data, k, B, parallel, conf, verbose) {
    if (nrow(data) < ncol(data)) data <- t(data)
    lWk <- km.clustering(data=data, k=k, ns=100, maxiter=1000)$lWk
    n <- nrow(data)
    p <- ncol(data)
    if (!parallel) {
        lWk.mc <- withinsumsq(n=n, p=p, B=B, k=k)
    } else {
        if (conf$type == "SOCK") {
            cl <- makeCluster(spec=conf$names,
                              type=conf$type,
                              homogeneous=conf$homo,
                              outfile=conf$outfile,
                              verbose=conf$verbose)            
        } else {
            cl <- makeCluster(spec=conf$cpus,
                              type=conf$type,
                              homogeneous=conf$homo,
                              outfile=conf$outfile,
                              verbose=conf$verbose)
        }
        clusterEvalQ(cl=cl, expr=library("MVR"))
        clusterSetRNGStream(cl=cl)
        lWk.cl <- clusterCall(cl=cl, fun=withinsumsq, n=n, p=p, B=ceiling(B/conf$cpus), k=k)
        stopCluster(cl)
        lWk.mc <- numeric(length=0)
        for (j in 1:conf$cpus) {
            lWk.mc <- c(lWk.mc, lWk.cl[[j]])
        }
    }
    lWk.mc.bar <- sum(lWk.mc)/B
    gapk <- abs(lWk.mc.bar - lWk)
    sdk <- sqrt(sum((lWk.mc - lWk.mc.bar)^2)/B)
    return(list(gapk=gapk, sk=sdk * sqrt(1 + (1/B))))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   merging.cluster(M)
#
################
# Description   :
################
#                   Internal merging function of variable cluster configurations determined for each group.
#                   Called by mvr(). Takes as argument the cluster membership matrix for each variable by group.
#
################
# Arguments     :
################
# M             :   Cluster membership matrix for each variable by group, where variables are by rows, and groups by columns 
#                   Returned by MeanVarReg()$membership.
#
################
# Values        :
################
# clus          :   Merged clustering configuration
#
################

merging.cluster <- function(M) {
    ord <- order(M[, 1, drop=TRUE], decreasing=FALSE)
    N <- M[ord, , drop=FALSE]
    clus <- N[, 1, drop=TRUE]
    if (ncol(M) >= 2) {
        for (g in 2:ncol(M)) {
            memb <- clus
            umemb <- unique(memb)
            clus <- numeric(0)
            i <- 0
            pma <- pmatch(x=names(memb), table=names(N[,g]))
            x <- N[pma,g]
            for (k in umemb) {
                y <- sort(x[memb == k])
                u <- unique(y)
                l <- length(u)
                z <- y
                for (j in 1:l) z[y == u[j]] <- j+i         
                i <- i+l
                clus <- c(clus, z)
            }
        }
    }
    return(clus)
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   mvrt(x, obj, lev, tab, ng, def, nc.min, nc.max, parallel, conf, verbose)    
#
################
# Description   :
################
#                   Internal function called by mvrt.test() for computing the regularized test-statistic
#
################
# Arguments     :
################
# x             :   Input matrix with variables in columns
# obj           :   Object of class "mvr" as returned by mvr().
# lev           :   Levels of the group/blocking factor used in mvrt.test() 
# tab           :   Number of samples per group
# ng            :   Number of sample groups
# def           :   List of samples indices per sample group
# nc.min        :   Minimum number of clusters
# nc.max        :   Maximum number of clusters
# parallel      :   Is parallel computing to be performed? Optional, defaults to FALSE.
# conf          :   List of parameters for cluster configuration. 
#                   Inputs for R package parallel function makeCluster() for cluster setup.
#                   Optional, defaults to NULL. See details for usage.
# verbose       :   Is the output to be verbose?
#
################
# Values        :
################
# t.reg         :   Vector of test-statistic values
#
################

mvrt <- function(x, obj, lev, tab, ng, def, nc.min, nc.max, parallel, conf, verbose) {
    p <- ncol(x)
    aver <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    sdev <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    for (g in 1:ng) {
        aver[g,] <- apply(x[def[[g]], , drop=FALSE], 2, mean, na.rm=TRUE)
        sdev[g,] <- apply(x[def[[g]], , drop=FALSE], 2, sd, na.rm=TRUE)
    }
    M <- matrix(data=NA, nrow=p, ncol=ng, dimnames=list(colnames(x), lev))
    if (is.null(obj)) {
        for (g in 1:ng) {
            if (verbose) cat("Computing optimal cluster configuration of group ", g, " ... \n", sep="")
            M[,g] <- MeanVarReg(data=x[def[[g]], , drop=FALSE], nc.min=nc.min, nc.max=nc.max, probs=NULL, B=100, parallel=parallel, conf=conf, verbose=verbose)$membership
            if (verbose) cat("group : ", g, " done\n")
        }
    } else {
        if (verbose) cat("Retrieving clustering information from 'mvr' object ... \n")
        for (g in 1:ng) {
            M[,g] <- obj$MVR[[g]]$membership
            if (verbose) cat("group : ", g, " done\n")
        }
    }
    clus <- merging.cluster(M)
    aver.reg <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    sdev.reg <- matrix(data=NA, nrow=ng, ncol=p, dimnames=list(lev, colnames(x)))
    for (k in unique(clus)) {
        pma <- pmatch(x=colnames(x), table=names(clus))
        ind <- (clus == k)[pma]
        for (g in 1:ng) {
            aver.reg[g,ind] <- mean(aver[g,ind]) 
            sdev.reg[g,ind] <- sqrt(mean(sdev[g,ind]^2))
        }
    }
    t.reg <- (aver.reg[1,,drop=TRUE] - aver.reg[2,,drop=TRUE]) / sqrt(sdev.reg[1,,drop=TRUE]^2/tab[1] + sdev.reg[2,,drop=TRUE]^2/tab[2])
    return(t.reg)
}
##########################################################################################################################################


