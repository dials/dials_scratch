# test convolution of kernels
kconv <- function (spans=c(2,2))
{
    x = kernel('modified.daniell', spans[1]%/%2)
    k = kernel('modified.daniell', spans[2]%/%2)
    if (!is.tskernel(x))
        stop("'x' is not a kernel")
    if (!is.tskernel(k))
        stop("'k' is not a kernel")
    n <- k$m
    cat("n is:", n, "\n")
    xx <- c(rep_len(0, n), x[-x$m:x$m], rep_len(0, n))
    cat("xx is:", xx, "\n")
    cat("k$coef is:", k$coef, "\n")
    coef <- kernapply(xx, k, circular = TRUE)
    cat("full coef is:", coef, "\n")
    m <- length(coef)%/%2L
    cat("m:", m, "\n")
    cat("sliced coef is:", coef[(m + 1L):length(coef)], "\n")
    kernel(coef[(m + 1L):length(coef)], m, paste0("Composite(",
        attr(x, "name"), ",", attr(k, "name"), ")"))
}

