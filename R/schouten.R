schouten <- function(X,weights=c("unweighted","linear","quadratic","user"),w=NULL,score=NULL) {
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
    if (is.null(score)) score <- 1:L
    if (is.null(w)) {
        weights <- match.arg(weights)
        w <- outer(score,score,
                   switch(weights,
                    unweighted=function(x,y) (x==y)+0,
                    linear=function(x,y) 1.0-abs(x-y)/(max(x)-min(x)),
                    quadratic=function(x,y) 1.0-(x-y)^2/(max(x)-min(x))^2))
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


summary.schouten <- function(object, ci.transform=c("logit","identity"), ci.p=0.95, ...) {
    ci.transform <- match.arg(ci.transform)
    expit <- function(x) 1/(1+exp(-x))
    logit <- function(x) log(x/(1-x))
    alpha <- 1-ci.p
    transf <- switch(ci.transform,
                     identity=function(mean,se) mean+qnorm(c(alpha/2,1-alpha/2))*se,
                     logit=function(mean,se)
                     expit(logit(mean)+qnorm(c(alpha/2,1-alpha/2))*se/mean/(1-mean)))
    object$ci <- transf(object$kappa,sqrt(object$vark))
    object$cia <- t(mapply(transf,object$ka,sqrt(object$varka)))
    rownames(object$cia) <- names(object$ka)
    ## object$ciab <- mapply(transf,object$kab,sqrt(object$varkab))
    structure(object, class="summary.schouten")
}

print.schouten <- function(x, ...) {
    weight.labels <- switch(x$weights,
                            unweighted="unweighted",
                            linear="linear weights",
                            quadratic="quadratic weights",
                            user="user-defined weights")
    if (!inherits(x,"summary.schouten")) x <- summary(x,...)
    cat(sprintf("O'Connell-Dobson-Schouten estimator (%s)\n\nAverage kappa:\t%f (se: %f; 95%% CI: %f, %f)\n",
                weight.labels,
                x$kappa,
                sqrt(x$vark),
                x$ci[1],
                x$ci[2],
                ...))
    cat(sprintf("Pr(Overall agreement due to chance):\t%s\n",
                format.pval(x$p0), ...))
    invisible(x)
}

.print.summary.schouten.ka <- function(x) {
    ka <- cbind(Kappa=x$ka,`[Lower,`=x$cia[,1],`Upper]`=x$cia[,2],`Pr(kappa_av=kappa_observer)`=x$pchi)
    stats::printCoefmat(ka,eps.Pvalue=1e-5)
    invisible(x)
}
print.summary.schouten <- function(x, ...) {
    print.schouten(x,...)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(apply(x$pa,1:2,sum))
    cat("\nAgreement statistics for each observer:\n\n")
    .print.summary.schouten.ka(x)
    invisible(x)
}

plot.schouten <- function(x, type=c("kappa by observer"), xlab=NULL, ylab=NULL, main=NULL, xdelta=0.1, axes=TRUE, ...) {
    type <- match.arg(type)
    ## if (type=="p1") {
    ##     if (is.null(xlab)) xlab <- "Observer"
    ##     if (is.null(ylab)) ylab <- "Probability"
    ##     if (is.null(main)) main <- ""
    ##     graphics:::plot.table(x$p1, xlab=xlab, ylab=ylab, main=main, ...)
    ## }
    if (type=="kappa by observer") {
        su <- summary(x)
        M <- length(x$ka)
        if (is.null(xlab)) xlab <- "Observer"
        if (is.null(ylab)) ylab <- "Kappa"
        matplot(su$cia,type="n",ylab=ylab,xlab=xlab,axes=FALSE,...)
        if (axes) {
            axis(1,at=1:M,labels=names(x$ka))
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
