schouten <- function(X,weights=c("unweighted","linear","quadratic","user"),w=NULL) {
    stopifnot(inherits(X,"data.frame") || inherits(X,"matrix"))
    stopifnot(ncol(X)>1 && nrow(X)>1)
    N <- nrow(X)
    M <- ncol(X)
    matX <- as.matrix(X)
    ## check for NAs, 
    if (any(is.na(matX)))
        stop("NAs in matrix: missing values or columns of different types?")
    score.labels <- sort(unique(as.vector(matX)))
    L <- length(unique(unlist(as.vector(matX))))
    if (is.null(w)) {
        weights <- match.arg(weights)
        w <- switch(weights,
                    unweighted=diag(L),
                    linear=outer(1:L,1:L,function(x,y) 1.0-abs(x-y)/(L-1)),
                    quadratic=outer(1:L,1:L,function(x,y) 1.0-(x-y)^2/(L-1)^2))
    } else weights <- "user"
    data2 <- matrix(as.integer(as.factor(matX)),N,M)
    raterNames <- colnames(X)
    itemNames <- rownames(X)
    out <- .Fortran("schouten",
                    N=as.integer(N),
                    M=as.integer(M),
                    L=as.integer(L),
                    data=data2,
                    w=w,
                    kab=matrix(0.0,M,M),
                    ka=rep(0.0,M),
                    kappa=0.0,
                    pab=array(0.0,c(M,M,L,L)),
                    pa=array(0.0,c(M,L,L)),
                    p=array(0.0,c(L,L)),
                    ma=array(0.0,c(M,L)),
                    qab=array(0.0,c(M,M,L,L)),
                    qa=array(0.0,c(M,L,L)),
                    q=array(0.0,c(L,L)),
                    oab=array(0.0,c(M,M)),
                    eab=array(0.0,c(M,M)),
                    oa=rep(0.0,M),
                    ea=rep(0.0,M),o=0.0,e=0.0,
                    wa=array(0.0,c(M,L)),
                    varkab=array(0.0,c(M,M)),
                    varka=rep(0.0,M),
                    vark=0.0,
                    covkka=rep(0.0,M),
                    chi=rep(0.0,M),
                    var0kab=array(0.0,c(M,M)),
                    var0ka=rep(0.0,M),
                    var0k=0.0,
                    wab=array(0.0,c(M,M)))
    out$pchi <- pchisq(out$chi,df=1,lower.tail=FALSE)
    dimnames(out$pab) <- dimnames(out$qab) <- list(raterNames,raterNames,score.labels,score.labels)
    dimnames(out$pa) <- list(raterNames,score.labels,score.labels)
    rownames(out$data) <- itemNames
    colnames(out$data) <- rownames(out$kab) <- colnames(out$kab) <- rownames(out$varkab) <- colnames(out$varkab) <-
        rownames(out$var0kab) <- colnames(out$var0kab) <-
            names(out$varka) <- names(out$var0ka) <- names(out$ka) <- raterNames
    ## colnames(out$p1) <- names(out$p2) <- colnames(out$w1) <- names(score) <- score.labels
    out$p0 <- 2*pnorm(abs(out$kappa/sqrt(out$var0k)-1),lower.tail=FALSE)
    out$p0a <- 2*pnorm(abs(out$ka/sqrt(out$var0ka)-1),lower.tail=FALSE)
    out$weights <- weights
    out$X <- X
    out$raterNames <- raterNames
    out$itemNames <- itemNames
    out$call <- sys.call()
    structure(out,class="schouten")
}


summary.schouten <- function(obj, ci.transform=c("logit","identity"), ci.p=0.95, ...) {
    ci.transform <- match.arg(ci.transform)
    expit <- function(x) 1/(1+exp(-x))
    logit <- function(x) log(x/(1-x))
    alpha <- 1-ci.p
    transf <- switch(ci.transform,
                     identity=function(mean,se) mean+qnorm(c(alpha/2,1-alpha/2))*se,
                     logit=function(mean,se)
                     expit(logit(mean)+qnorm(c(alpha/2,1-alpha/2))*se/mean/(1-mean)))
    obj$ci <- transf(obj$kappa,sqrt(obj$vark))
    obj$cia <- t(mapply(transf,obj$ka,sqrt(obj$varka)))
    rownames(obj$cia) <- names(obj$ka)
    ## obj$ciab <- mapply(transf,obj$kab,sqrt(obj$varkab))
    structure(obj, class="summary.schouten")
}

print.schouten <- function(obj, ...) {
    weight.labels <- switch(obj$weights,
                            unweighted="unweighted",
                            linear="linear weights",
                            quadratic="quadratic weights",
                            user="user-defined weights")
    if (!inherits(obj,"summary.schouten")) obj <- summary(obj,...)
    cat(sprintf("Schouten estimator (%s)\n\nAverage kappa:\t%f (se: %f; 95%% CI: %f, %f)\n",
                weight.labels,
                obj$kappa,
                sqrt(obj$vark),
                obj$ci[1],
                obj$ci[2],
                ...))
    cat(sprintf("Pr(Overall agreement due to chance):\t%s\n",
                format.pval(obj$p0), ...))
    invisible(obj)
}

.print.summary.schouten.ka <- function(obj) {
    ka <- cbind(Kappa=obj$ka,`[Lower,`=obj$cia[,1],`Upper]`=obj$cia[,2],`Pr(kappa_av=kappa_rater)`=obj$pchi)
    stats::printCoefmat(ka,eps.Pvalue=1e-5)
    invisible(obj)
}
print.summary.schouten <- function(obj, ...) {
    print.schouten(obj,...)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(apply(obj$pa,1:2,sum))
    cat("\nAgreement statistics for each rater:\n\n")
    .print.summary.schouten.ka(obj)
    invisible(obj)
}

plot.schouten <- function(obj, type=c("kappa by rater"), xlab=NULL, ylab=NULL, main=NULL, xdelta=0.1, axes=TRUE, ...) {
    type <- match.arg(type)
    ## if (type=="p1") {
    ##     if (is.null(xlab)) xlab <- "Rater"
    ##     if (is.null(ylab)) ylab <- "Probability"
    ##     if (is.null(main)) main <- ""
    ##     graphics:::plot.table(obj$p1, xlab=xlab, ylab=ylab, main=main, ...)
    ## }
    if (type=="kappa by rater") {
        su <- summary(obj)
        M <- length(obj$ka)
        if (is.null(xlab)) xlab <- "Rater"
        if (is.null(ylab)) ylab <- "Kappa"
        matplot(su$cia,type="n",ylab=ylab,xlab=xlab,axes=FALSE,...)
        if (axes) {
            axis(1,at=1:M,labels=names(obj$ka))
            axis(2)
            box()
        }
        lower <- su$cia[,1]
        upper <- su$cia[,2]
        segments(1:M,lower,1:M,upper)
        segments((1:M)-xdelta,lower,(1:M)+xdelta,lower)
        segments((1:M)-xdelta,upper,(1:M)+xdelta,upper)
        points(1:M,su$ka,pch=22,bg=1,cex=1.5)
    }
}
