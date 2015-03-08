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
#                   cluster.diagnostic(obj,
#                                      span=0.75,
#                                      degree=2,
#                                      family="gaussian",
#                                      title="",
#                                      device=NULL,
#                                      file="Cluster Diagnostic Plots",
#                                      path=getwd()
#                                      horizontal=FALSE,
#                                      width=8.5,
#                                      height=11, ...)
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
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (standard output screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Cluster Diagnostic Plots".
# path          :   Absolute path (without final (back)slash separator). Defaults to working directory path.
# horizontal    :   Orientation of the printed image, a logical. Defaults to FALSE, that is potrait orientation.
# width         :   Width of the graphics region in inches. Defaults to 8.5.
# height        :   Height of the graphics region in inches. Defaults to 11.
# ...           :   Generic arguments passed to other plotting functions.
#
################
# Values        :
################
#
################

cluster.diagnostic <- function(obj,
                               span=0.75,
                               degree=2,
                               family="gaussian",
                               title="",
                               device=NULL,
                               file="Cluster Diagnostic Plots",
                               path=getwd(),
                               horizontal=FALSE,
                               width=8.5,
                               height=11, ...) {

    block <- obj$block
    lev <- levels(block)
    ng <- nlevels(block)

    clusterplot <- function(obj, ng, title, span, degree, family, ...) {
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
        dev.new(width=width, height=height, title="Cluster Diagnostic Plots", noRStudioGD = TRUE)        
        clusterplot(obj=obj, ng=ng, title=title, span=span, degree=degree, family=family)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        clusterplot(obj=obj, ng=ng, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
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
#                   target.diagnostic(obj,
#                                     title="",
#                                     device=NULL,
#                                     file="Target Moments Diagnostic Plots",
#                                     path=getwd(),
#                                     horizontal=FALSE,
#                                     width=8.5,
#                                     height=6.5, ...)
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
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (standard output screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Target Moments Diagnostic Plots".
# path          :   Absolute path (without final (back)slash separator). Defaults to working directory path.
# horizontal    :   Orientation of the printed image, a logical. Defaults to FALSE, that is potrait orientation.
# width         :   Width of the graphics region in inches. Defaults to 8.5.
# height        :   Height of the graphics region in inches. Defaults to 6.5.
# ...           :   Generic arguments passed to other plotting functions.
#
################
# Values        :
################
#
################

target.diagnostic <- function(obj,
                              title="",
                              device=NULL,
                              file="Target Moments Diagnostic Plots",
                              path=getwd(),
                              horizontal=FALSE,
                              width=8.5,
                              height=6.5, ...) {

    targetplot <- function(obj, title, ...) {
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
        dev.new(width=width, height=height, title="Target Moments Diagnostic Plots", noRStudioGD = TRUE)        
        targetplot(obj=obj, title=title)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        targetplot(obj=obj, title=title)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
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
#                   stabilization.diagnostic(obj,
#                                            span=0.5,
#                                            degree=2,
#                                            family="gaussian",
#                                            title="",
#                                            device=NULL,
#                                            file="Stabilization Diagnostic Plots",
#                                            path=getwd(),
#                                            horizontal=FALSE,
#                                            width=7,
#                                            height=5, ...)
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
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (standard output screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Stabilization Diagnostic Plots".
# path          :   Absolute path (without final (back)slash separator). Defaults to working directory path.
# horizontal    :   Orientation of the printed image, a logical. Defaults to FALSE, that is potrait orientation.
# width         :   Width of the graphics region in inches. Defaults to 7.
# height        :   Height of the graphics region in inches. Defaults to 5.
# ...           :   Generic arguments passed to other plotting functions.
#
################
# Values        :
################
#
################

stabilization.diagnostic <- function(obj,
                                     span=0.5,
                                     degree=2,
                                     family="gaussian",
                                     title="",
                                     device=NULL,
                                     file="Stabilization Diagnostic Plots",
                                     path=getwd(),
                                     horizontal=FALSE,
                                     width=7,
                                     height=5, ...) {

    stabplot <- function(obj, title, span, degree, family, ...) {
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
        dev.new(width=width, height=height, title="Stabilization Diagnostic Plots", noRStudioGD = TRUE)        
        stabplot(obj=obj, title=title, span=span, degree=degree, family=family)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        stabplot(obj=obj, title=title, span=span, degree=degree, family=family)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
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
#                   normalization.diagnostic(obj,
#                                            pal,
#                                            title="",
#                                            device=NULL,
#                                            file="Normalization Diagnostic Plots",
#                                            path=getwd(),
#                                            horizontal=FALSE,
#                                            width=7,
#                                            height=8, ...)
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
# device        :   Display device in {NULL, "PS", "PDF"}. Defaults to NULL (standard output screen).
#                   Currently implemented display device are "PS" (Postscript) or "PDF" (Portable Document Format).
# file          :   File name for outputting display device. Defaults to "Normalization Diagnostic Plots".
# path          :   Absolute path (without final (back)slash separator). Defaults to working directory path.
# horizontal    :   Orientation of the printed image, a logical. Defaults to FALSE, that is potrait orientation.
# width         :   Width of the graphics region in inches. Defaults to 7.
# height        :   Height of the graphics region in inches. Defaults to 8.
# ...           :   Generic arguments passed to other plotting functions.
#
################
# Values        :
################
#
################

normalization.diagnostic <- function(obj,
                                     pal,
                                     title="",
                                     device=NULL,
                                     file="Normalization Diagnostic Plots",
                                     path=getwd(),
                                     horizontal=FALSE,
                                     width=7,
                                     height=8, ...) {

    normplot <- function(obj, title, ...) {
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
        dev.new(width=width, height=height, title="Normalization Diagnostic Plots", noRStudioGD = TRUE)        
        normplot(obj=obj, title=title)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        normplot(obj=obj, title=title)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        normplot(obj=obj, title=title)
        dev.off()
    } else {
        stop("Currently allowed display device are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }

}
##########################################################################################################################################




##########################################################################################################################################
# 4. OTHER END-USER FUNCTIONS
##########################################################################################################################################

##########################################################################################################################################
#################
#Usage         :
################
#                   MVR.news(...) 
#
################
# Description   :
################
#                   Function to display the log file of updates of the MVR package.
#
################
# Arguments     :
################
# ...               Further arguments passed to or from other methods.
#
################
# Values        :
################
#                   None.
#
##########################################################################################################################################

MVR.news <- function(...) {
    newsfile <- file.path(system.file(package="MVR"), "NEWS")
    file.show(newsfile)
}
##########################################################################################################################################

