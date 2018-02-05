kernapply2 <- function (x, ...)
{
    UseMethod("kernapply2")
}

kernapply2.vector <- function (x, k, circular = FALSE, ...)
{
    cat('kernapply.vector\n')
    if (!is.vector(x)) stop ("'x' is not a vector")
    if (!is.tskernel(k)) stop ("'k' is not a kernel")
    m <- k$m
    if (length(x) <= 2L*m)
        stop ("'x' is shorter than kernel 'k'")
    if (m == 0L)
        return (x)
    else
    {
        n <- length(x)
        w <- c(k[0L:m], rep_len(0,n-2L*m-1L), k[-m:-1L])
        y <- fft(fft(x)*fft(w), inverse = TRUE)/n
        if (is.numeric(x)) y <- Re(y)
        if (circular)
            return (y)
        else
            return (y[(1L+m):(n-m)])
    }
}

kernapply2.default <- function (x, k, circular = FALSE, ...)
{
    cat('kernapply.default\n')
    if (is.vector(x))
        return (kernapply.vector(x, k, circular=circular))
    else if (is.matrix(x))
        return (apply(x, MARGIN=2, FUN=kernapply, k, circular=circular))
    else
        stop ("'kernapply' is not available for object 'x'")
}

kernapply2.ts <- function (x, k, circular = FALSE, ...)
{
    cat('kernapply.ts\n')
    if (!is.matrix(x))
        y <- kernapply.vector(as.vector(x), k, circular=circular)
    else
        y <- apply(x, MARGIN=2L, FUN=kernapply, k, circular=circular)
    ts (y, end=end(x), frequency=frequency(x))
}

kernapply2.tskernel <- function (x, k, ...)
{
    cat('kernapply.tskernel\n')
    if (!is.tskernel(x))
        stop ("'x' is not a kernel")
    if (!is.tskernel(k))
        stop ("'k' is not a kernel")
    n <- k$m
    xx <- c(rep_len(0,n), x[-x$m:x$m], rep_len(0,n))
    coef <- kernapply(xx, k, circular = TRUE)
    m <- length(coef) %/% 2L
    kernel(coef[(m+1L):length(coef)], m,
           paste0("Composite(", attr(x, "name"), ",", attr(k, "name"), ")"))
}
