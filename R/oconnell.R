oconnell <-
function(X,score=NULL, i=3) {
    nrater <- ncol(X); nsubj <- nrow(X)
    raterNames <- colnames(X)
    itemNames <- rownames(X)
    X <- as.matrix(X)
    score.labels <- sort(unique(as.vector(X)))
    if (is.null(score)) score <- score.labels
    X <- matrix(as.integer(as.factor(X)),nrow(X),ncol(X))
    nscore <- length(score)
    out <- .Fortran("oconnell", nrater=nrater, nsubj=nsubj, nscore=nscore, i=as.integer(i),
                    resp=t(X), s1=rep(0,nsubj), s2=rep(0,nsubj),
                    p1=matrix(0,nrater,nscore),p2=rep(0,nscore),w1=matrix(0,nrater,nscore),
                    w2=matrix(0,nrater,nscore),
                    score=as.double(score),delta=matrix(0,nrater,nrater),
                    d=rep(0,nsubj),expd1=0,expd2=0,dbar=0,vars1=0,var0s1=0,vsav1=0,v0sav1=0,
                    vars2=0,var0s2=0,vsav2=0,v0sav2=0,sav1=0,sav2=0)
    colnames(out$resp) <- names(out$s1) <- names(out$s2) <- names(out$d) <- itemNames
    rownames(out$resp) <- rownames(out$p1) <- rownames(out$w1) <- rownames(out$w2) <- rownames(out$delta) <- colnames(out$delta) <- raterNames
    colnames(out$p1) <- names(out$p2) <- colnames(out$w1) <- names(score) <- score.labels
    out$X <- X
    out$resp <- NULL
    class(out) <- "oconnell"
    out
}

print.oconnell <- function(obj, ...) {
    cat(sprintf("O'Connell-Dobson estimator\n\nSav(hetero):\t%f (se=%f)\nSav(homoge):\t%f (se=%f)\n",
            obj$sav1,
            sqrt(obj$vsav1),
            obj$sav2,
            sqrt(obj$vsav2),...))
    cat(sprintf("Pr(Overall agreement due to chance | hetero):\t%g\n",pnorm(abs(obj$sav1/sqrt(obj$v0sav1)-1),lower.tail=FALSE)), ...)
    cat(sprintf("Pr(Overall agreement due to chance | homoge):\t%g\n",pnorm(abs(obj$sav2/sqrt(obj$v0sav2)-1),lower.tail=FALSE)), ...)
}

summary.oconnell <- function(obj) structure(obj, class="summary.oconnell") # it's all there already
print.summary.oconnell <- function(obj, ...) {
    print.oconnell(obj,...)
    cat("\nObserved marginal distributions for categories:\n\n")
    print(obj$p2)
    cat("\nObserved marginal distributions for categories by observer:\n\n")
    print(obj$p1)
    cat("\nAgreement statistics S_i for the individual items:\n\n")
    print(obj$s1)
}


## Landis and Kock (1977) - output this to ./data/?
landis <- structure(c(4L, 1L, 3L, 4L, 3L, 2L, 1L, 3L, 2L, 1L, 5L, 1L, 3L, 
2L, 4L, 3L, 2L, 2L, 2L, 1L, 4L, 1L, 1L, 2L, 4L, 3L, 3L, 1L, 4L, 
3L, 1L, 3L, 2L, 3L, 5L, 2L, 3L, 3L, 5L, 5L, 3L, 1L, 2L, 4L, 3L, 
3L, 2L, 3L, 4L, 3L, 3L, 2L, 2L, 1L, 3L, 1L, 1L, 4L, 1L, 2L, 4L, 
3L, 1L, 2L, 3L, 1L, 4L, 3L, 3L, 4L, 1L, 2L, 2L, 2L, 4L, 1L, 4L, 
5L, 2L, 4L, 3L, 4L, 4L, 2L, 3L, 3L, 4L, 3L, 1L, 3L, 4L, 4L, 1L, 
3L, 4L, 3L, 1L, 2L, 3L, 2L, 3L, 3L, 2L, 1L, 3L, 3L, 2L, 3L, 1L, 
3L, 3L, 1L, 1L, 2L, 5L, 4L, 1L, 2L, 3L, 1L, 3L, 3L, 3L, 1L, 1L, 
3L, 2L, 1L, 5L, 1L, 3L, 2L, 3L, 3L, 3L, 1L, 3L, 1L, 3L, 1L, 1L, 
1L, 4L, 3L, 3L, 1L, 3L, 3L, 1L, 3L, 2L, 3L, 3L, 1L, 3L, 3L, 5L, 
3L, 2L, 1L, 3L, 4L, 3L, 2L, 3L, 3L, 3L, 3L, 3L, 2L, 3L, 1L, 3L, 
1L, 3L, 3L, 3L, 3L, 3L, 3L, 1L, 3L, 3L, 1L, 3L, 3L, 3L, 3L, 2L, 
2L, 3L, 1L, 4L, 1L, 4L, 5L, 3L, 4L, 3L, 3L, 2L, 3L, 3L, 3L, 4L, 
3L, 1L, 3L, 3L, 3L, 2L, 3L, 4L, 3L, 1L, 3L, 3L, 3L, 3L, 3L, 2L, 
1L, 3L, 3L, 3L, 3L, 1L, 3L, 3L, 1L, 1L, 2L, 3L, 3L, 1L, 3L, 4L, 
1L, 3L, 3L, 3L, 2L, 1L, 2L, 2L, 1L, 5L, 1L, 3L, 2L, 3L, 2L, 2L, 
2L, 2L, 2L, 3L, 2L, 1L, 2L, 4L, 3L, 3L, 1L, 3L, 3L, 1L, 3L, 2L, 
2L, 3L, 1L, 2L, 3L, 5L, 3L, 2L, 1L, 1L, 4L, 3L, 2L, 2L, 3L, 3L, 
2L, 3L, 2L, 2L, 1L, 3L, 2L, 2L, 3L, 2L, 2L, 3L, 3L, 1L, 2L, 2L, 
1L, 3L, 3L, 3L, 1L, 1L, 1L, 2L, 1L, 3L, 1L, 3L, 1L, 2L, 4L, 2L, 
3L, 3L, 2L, 3L, 2L, 3L, 2L, 2L, 3L, 1L, 3L, 2L, 3L, 3L, 2L, 1L, 
2L, 3L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 3L, 2L, 1L, 1L, 
1L, 4L, 4L, 1L, 1L, 2L, 1L, 3L, 4L, 3L, 1L, 1L, 3L, 2L, 1L, 4L, 
1L, 2L, 1L, 2L, 3L, 2L, 1L, 2L, 1L, 4L, 1L, 1L, 2L, 2L, 2L, 3L, 
1L, 3L, 3L, 1L, 2L, 2L, 2L, 3L, 1L, 2L, 3L, 5L, 2L, 2L, 1L, 2L, 
3L, 2L, 2L, 2L, 4L, 3L, 2L, 3L, 1L, 2L, 1L, 3L, 1L, 1L, 3L, 2L, 
2L, 3L, 4L, 1L, 2L, 3L, 1L, 3L, 2L, 3L, 3L, 1L, 2L, 1L, 2L, 2L, 
1L, 3L, 4L, 2L, 2L, 3L, 3L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 
4L, 1L, 2L, 4L, 2L, 1L, 2L, 3L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 
2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 3L, 1L, 3L, 3L, 3L, 
1L, 2L, 2L, 3L, 2L, 5L, 2L, 3L, 1L, 3L, 3L, 3L, 2L, 2L, 1L, 3L, 
2L, 1L, 2L, 4L, 3L, 3L, 1L, 3L, 3L, 1L, 3L, 3L, 3L, 4L, 2L, 3L, 
3L, 5L, 3L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 3L, 3L, 4L, 3L, 2L, 3L, 
1L, 3L, 1L, 2L, 3L, 2L, 3L, 3L, 3L, 1L, 3L, 3L, 1L, 3L, 3L, 3L, 
3L, 1L, 2L, 3L, 1L, 4L, 1L, 4L, 5L, 2L, 5L, 3L, 3L, 3L, 4L, 4L, 
3L, 4L, 3L, 2L, 4L, 2L, 4L, 2L, 4L, 4L, 3L, 1L, 4L, 3L, 3L, 3L, 
3L, 2L, 2L, 2L, 2L, 2L, 3L, 2L, 3L, 3L, 2L, 1L, 2L, 3L, 4L, 2L, 
2L, 3L, 1L, 3L, 3L, 3L, 1L, 1L, 2L, 1L, 1L, 5L, 1L, 3L, 1L, 2L, 
3L, 2L, 1L, 1L, 1L, 3L, 1L, 1L, 1L, 3L, 2L, 2L, 1L, 2L, 3L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L, 5L, 2L, 1L, 1L, 1L, 3L, 2L, 1L, 2L, 
2L, 5L, 2L, 2L, 2L, 1L, 1L, 3L, 1L, 1L, 2L, 1L, 2L, 3L, 2L, 1L, 
2L, 1L, 1L, 3L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 3L, 5L, 1L, 
1L, 3L, 3L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 3L, 1L, 3L, 1L, 2L, 3L, 
3L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 
1L, 1L, 1L, 4L, 1L, 1L, 1L, 3L, 1L, 3L, 3L, 3L, 1L, 1L, 3L, 2L, 
1L, 5L, 1L, 3L, 2L, 3L, 3L, 3L, 1L, 3L, 1L, 3L, 1L, 1L, 2L, 3L, 
3L, 3L, 1L, 3L, 3L, 1L, 3L, 2L, 3L, 3L, 1L, 3L, 3L, 5L, 3L, 2L, 
1L, 3L, 3L, 3L, 1L, 2L, 3L, 3L, 3L, 3L, 2L, 3L, 1L, 3L, 1L, 1L, 
3L, 2L, 3L, 3L, 4L, 1L, 2L, 3L, 1L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 
1L, 3L, 1L, 3L, 4L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 2L, 3L, 2L, 1L, 
3L, 2L, 3L, 2L, 3L, 4L, 3L, 1L, 2L, 3L, 1L, 3L, 3L, 1L, 1L, 3L, 
2L, 1L, 3L, 1L, 3L, 3L, 1L, 1L, 2L, 3L, 3L, 1L, 2L), .Dim = c(118L, 
7L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", 
"9", "10", "11", "12", "13", "15", "16", "17", "18", "19", "22", 
"23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", 
"34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", 
"45", "46", "47", "48", "49", "51", "52", "53", "54", "55", "56", 
"57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", 
"68", "69", "70", "71", "72", "73", "74", "76", "77", "78", "79", 
"80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", 
"91", "92", "93", "94", "95", "96", "98", "99", "100", "101", 
"102", "103", "104", "105", "106", "107", "108", "110", "111", 
"112", "113", "114", "115", "116", "117", "118", "119", "120", 
"121", "122", "123", "124", "126"), c("A", "B", "C", "D", "E", 
"F", "G")))



