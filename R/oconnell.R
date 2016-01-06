oconnell <- function(X,score=NULL, weights=c("unweighted","linear","quadratic"), i=NULL) {
    ## check arguments
    weights <- match.arg(weights)
    stopifnot(inherits(X,"data.frame") || inherits(X,"matrix"))
    stopifnot(ncol(X)>1 && nrow(X)>1)
    ## if i is specified, then it takes precedence over weights, else:
    if (is.null(i)) i <- switch(weights,unweighted=1,linear=2,quadratic=3)
    matX <- as.matrix(X)
    ## check for NAs, 
    if (any(is.na(matX)))
        stop("NAs in matrix: missing values or columns of different types?")
    score.labels <- sort(unique(as.vector(matX)))
    if (is.null(score)) score <- as.numeric(score.labels)
    if (any(is.na(score)))
        stop("Missing or non-numeric scores. Possibly specify the score argument")
    nrater <- ncol(X); nsubj <- nrow(X)
    raterNames <- colnames(X)
    itemNames <- rownames(X)
    ## NOTE: resp is X -> factor -> integer -> matrix -> transpose
    resp <- t(matrix(as.integer(as.factor(matX)),nrow(X),ncol(X)))
    nscore <- length(score)
    out <- .Fortran("oconnell", nrater=nrater, nsubj=nsubj, nscore=nscore, i=as.integer(i),
                    resp=resp, s1=rep(0,nsubj), s2=rep(0,nsubj),
                    p1=matrix(0,nrater,nscore),p2=rep(0,nscore),w1=matrix(0,nrater,nscore),
                    w2=matrix(0,nrater,nscore),
                    score=as.double(score),delta=matrix(0,nrater,nrater),
                    d=rep(0,nsubj),expd1=0,expd2=0,dbar=0,vars1=0,var0s1=0,vsav1=0,v0sav1=0,
                    vars2=0,var0s2=0,vsav2=0,v0sav2=0,sav1=0,sav2=0)
    colnames(out$resp) <- names(out$s1) <- names(out$s2) <- names(out$d) <- itemNames
    rownames(out$resp) <- rownames(out$p1) <- rownames(out$w1) <- rownames(out$w2) <- rownames(out$delta) <- colnames(out$delta) <- raterNames
    colnames(out$p1) <- names(out$p2) <- colnames(out$w1) <- names(score) <- score.labels
    out$p0sav1 <- 2*pnorm(abs(out$sav1/sqrt(out$v0sav1)-1),lower.tail=FALSE)
    out$p0sav2 <- 2*pnorm(abs(out$sav2/sqrt(out$v0sav2)-1),lower.tail=FALSE)
    out$X <- X
    ## out$schouten <- schouten(X,weights)
    out$call <- sys.call()
    class(out) <- "oconnell"
    out
}

summary.oconnell <- function(obj, ci.transform=c("logit","identity"), ci.p=0.95, ...) {
    ci.transform <- match.arg(ci.transform)
    expit <- function(x) 1/(1+exp(-x))
    logit <- function(x) log(x/(1-x))
    alpha <- 1-ci.p
    transf <- switch(ci.transform,
                     identity=function(mean,se) mean+qnorm(c(alpha/2,1-alpha/2))*se,
                     logit=function(mean,se)
                     expit(logit(mean)+qnorm(c(alpha/2,1-alpha/2))*se/mean/(1-mean)))
    obj$ci1 <- transf(obj$sav1,sqrt(obj$vsav1))
    obj$ci2 <- transf(obj$sav2,sqrt(obj$vsav2))
    structure(obj, class="summary.oconnell")
}

print.oconnell <- function(obj, ...) {
    weight.labels <- c("unweighted","linear weights","quadratic weights")
    if (!inherits(obj,"summary.oconnell")) obj <- summary(obj,...)
    cat(sprintf("O'Connell-Dobson estimator (%s)\n\nSav(hetero):\t%f (se: %f; 95%% CI: %f, %f)\nSav(homoge):\t%f (se: %f; 95%% CI: %f, %f)\n",
                weight.labels[obj$i],
                obj$sav1,
                sqrt(obj$vsav1),
                obj$ci1[1],
                obj$ci1[2],
                obj$sav2,
                sqrt(obj$vsav2),
                obj$ci2[1],
                obj$ci2[2],
                ...))
    cat(sprintf("Pr(Overall agreement due to chance | hetero):\t%s\n",
                format.pval(obj$p0sav1), ...))
    cat(sprintf("Pr(Overall agreement due to chance | homoge):\t%s\n",
                format.pval(obj$p0sav2), ...))
    invisible(obj)
}

print.summary.oconnell <- function(obj, ...) {
    print.oconnell(obj,...)
    cat("\nObserved marginal distributions for categories:\n\n")
    print(obj$p2)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(obj$p1)
    cat("\nAgreement statistics S_i for the individual items:\n\n")
    print(obj$s1)
    ## cat("\nAgreement statistics for each rater:\n\n")
    ## .print.summary.schouten.ka(summary(obj$schouten))
    invisible(obj)
}

plot.oconnell <- function(obj, type=c("p1","kappa by rater"), xlab=NULL, ylab=NULL, main=NULL, ...) {
    type <- match.arg(type)
    if (type=="p1") {
        if (is.null(xlab)) xlab <- "Rater"
        if (is.null(ylab)) ylab <- "Probability"
        if (is.null(main)) main <- ""
        graphics:::plot.table(obj$p1, xlab=xlab, ylab=ylab, main=main, ...)
    }
    ## if (type=="kappa by rater") plot.schouten(obj$schouten,xlab=xlab,ylab=ylab,main=main,...)
}
