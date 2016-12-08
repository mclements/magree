magree <- function(X, weights=c("unweighted","linear","quadratic"), score=NULL) {
    weights <- match.arg(weights)
    structure(list(oconnell = magree::oconnell(X,weights,score=score),
                   schouten = magree::schouten(X,weights,score=score),
                   call = sys.call()),
              class="magree")
}

print.magree <- function(x, ...)
    print.oconnell(x$oconnell, ...)

summary.magree <- function(object, ...) {
    structure(list(oconnell = summary(object$oconnell),
                   schouten = summary(object$schouten)),
              class="summary.magree")
}

print.summary.magree <- function(x, ...) {
    print.oconnell(x$oconnell, ...)
    cat("\nObserved marginal distributions for categories:\n\n")
    print(x$oconnell$p2)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(x$oconnell$p1)
    cat("\nAgreement statistics S_i for each subject:\n\n")
    print(x$oconnell$s1)
    cat("\nAgreement statistics for each observer:\n\n")
    .print.summary.schouten.ka(x$schouten)
    invisible(x)
}

plot.magree <- function(x, type=c("p1","kappa by observer"), xlab=NULL, ylab=NULL, main=NULL, ...) {
    type <- match.arg(type)
    if (type=="p1") 
        plot.oconnell(x$oconnell, xlab=xlab, ylab=ylab, main=main, ...)
    if (type=="kappa by observer")
        plot.schouten(x$schouten,xlab=xlab,ylab=ylab,main=main,...)
}


oconnell <- function(X, weights=c("unweighted","linear","quadratic"), i=NULL, score=NULL) {
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
    names(score) <- score.labels
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
    out$score <- score
    ## out$schouten <- schouten(X,weights)
    out$call <- sys.call()
    class(out) <- "oconnell"
    out
}

summary.oconnell <- function(object, ci.transform=c("logit","identity"), ci.p=0.95, ...) {
    ci.transform <- match.arg(ci.transform)
    expit <- function(x) 1/(1+exp(-x))
    logit <- function(x) log(x/(1-x))
    alpha <- 1-ci.p
    transf <- switch(ci.transform,
                     identity=function(mean,se) mean+qnorm(c(alpha/2,1-alpha/2))*se,
                     logit=function(mean,se)
                     expit(logit(mean)+qnorm(c(alpha/2,1-alpha/2))*se/mean/(1-mean)))
    object$ci1 <- transf(object$sav1,sqrt(object$vsav1))
    object$ci2 <- transf(object$sav2,sqrt(object$vsav2))
    structure(object, class="summary.oconnell")
}

print.oconnell <- function(x, ...) {
    weight.labels <- c("unweighted","linear weights","quadratic weights")
    if (!inherits(x,"summary.oconnell")) x <- summary(x,...)
    cat(sprintf("O'Connell-Dobson-Schouten estimator (%s)\n\nSav(hetero):\t%f (se: %f; 95%% CI: %f, %f)\nSav(homoge):\t%f (se: %f; 95%% CI: %f, %f)\n",
                weight.labels[x$i],
                x$sav1,
                sqrt(x$vsav1),
                x$ci1[1],
                x$ci1[2],
                x$sav2,
                sqrt(x$vsav2),
                x$ci2[1],
                x$ci2[2],
                ...))
    cat(sprintf("Pr(Overall agreement due to chance | hetero):\t%s\n",
                format.pval(x$p0sav1), ...))
    cat(sprintf("Pr(Overall agreement due to chance | homoge):\t%s\n",
                format.pval(x$p0sav2), ...))
    invisible(x)
}

print.summary.oconnell <- function(x, ...) {
    print.oconnell(x,...)
    cat("\nObserved marginal distributions for categories:\n\n")
    print(x$p2)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(x$p1)
    cat("\nAgreement statistics S_i for each subject:\n\n")
    print(x$s1)
    ## cat("\nAgreement statistics for each observer:\n\n")
    ## .print.summary.schouten.ka(summary(x$schouten))
    invisible(x)
}

plot.oconnell <- function(x, type=c("p1"), xlab=NULL, ylab=NULL, main=NULL, ...) {
    type <- match.arg(type)
    if (type=="p1") {
        if (is.null(xlab)) xlab <- "Observer"
        if (is.null(ylab)) ylab <- "Probability"
        if (is.null(main)) main <- ""
        plot(as.table(x$p1), xlab=xlab, ylab=ylab, main=main, ...)
    }
}
