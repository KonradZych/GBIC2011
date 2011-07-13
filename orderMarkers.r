library(qtl)
orderMarkers <- function (cross, chr, window = 7, use.ripple = TRUE, error.prob = 1e-04, 
    map.function = c("haldane", "kosambi", "c-f", "morgan"), 
    maxit = 4000, tol = 1e-04, sex.sp = TRUE, verbose = FALSE) 
{
    map.function <- match.arg(map.function)
	cat("cat line 1\n")
    if (!missing(chr)){ 
        chr <- matchchr(chr, names(cross$geno))
    }else{
		chr <- names(cross$geno)
	}
	cat("cat line 2\n")
    n.mar <- nmar(cross)
    if (verbose > 1) 
        verbose.sub <- TRUE
    else verbose.sub <- FALSE
    for (i in chr) {
        if (verbose && length(chr) > 1) 
            cat(" - Chr", i, "\n")
        if (n.mar[i] > 2) {
            neworder <- orderMarkers.sub(cross, i, window = window, 
                use.ripple = use.ripple, verbose = verbose.sub)
            cross <- switch.order(cross, i, neworder, error.prob = error.prob, 
                map.function = map.function, maxit = maxit, tol = tol, 
                sex.sp = sex.sp)
        }
    }
    cross
}

orderMarkers.sub <- function (cross, chr, window = 7, use.ripple = TRUE, verbose = FALSE) 
{
    if (missing(chr)) 
        chr <- names(cross$geno)[1]
    if (length(chr) > 1) {
        if (length(grep("^-", chr)) > 0) 
            stop("Need to give a single chromosome name.")
        warning("Need to give a single chromosome name; using just the first")
        chr <- chr[1]
    }
    if (length(matchchr(chr, names(cross$geno))) > 1) 
        stop("Chr ", chr, " not found.")
    cross <- subset(cross, chr = chr)
    names(cross$geno)[1] <- "1"
    n.mar <- totmar(cross)
    if (n.mar < 3) 
        return(1:n.mar)
    if (use.ripple && n.mar <= window) {
        rip <- summary(ripple(cross, chr = 1, window = window, 
            verbose = FALSE))
        nxo <- rip[1:2, ncol(rip)]
        if (nxo[1] <= nxo[2]) 
            return(1:n.mar)
        else return(rip[2, 1:(ncol(rip) - 1)])
    }
    nt <- ntyped(cross, "mar")
    themar <- order(nt, decreasing = TRUE)
    marnam <- markernames(cross)
    if (n.mar > 3) {
        for (i in 4:n.mar) cross <- movemarker(cross, marnam[themar[i]], 
            2)
    }
	cat("cat line 3\n")
    makeorders <- function(n) {
        orders <- matrix(ncol = n, nrow = n)
        for (k in 1:n) {
            orders[k, n - k + 1] <- n
            orders[k, -(n - k + 1)] <- 1:(n - 1)
        }
        orders
    }
    simpleswitch <- function(cross, neworder) {
        cross$geno[[1]]$data <- cross$geno[[1]]$data[, neworder]
        if (is.matrix(cross$geno[[1]]$map)) 
            cross$geno[[1]]$map <- cross$geno[[1]]$map[, neworder]
        else cross$geno[[1]]$map <- cross$geno[[1]]$map[neworder]
		cross
    }
    if (verbose) 
        cat(" --- Adding marker 3 of", n.mar, "\n")
    orders <- makeorders(3)
    nxo <- rep(NA, nrow(orders))
    nxo[1] <- sum(countXO(cross, 1))
    temp <- cross
	cat("cat line 4\n")
    for (kk in 2:nrow(orders)) {
        temp$geno[[1]]$data <- temp$geno[[1]]$data[, orders[kk, 
            ]]
        nxo[kk] <- sum(countXO(cross, 1))
		cat("cat line 5\n")
    }
	cat("cat line 6\n")
    wh <- which(nxo == min(nxo))
    if (length(wh) > 1) 
        wh <- sample(wh, 1)
    if (wh > 1) 
        cross <- simpleswitch(cross, orders[wh, ])
    if (n.mar > 3) {
		cat("cat line 7\n")
        for (k in 4:n.mar) {
			cat("cat line 8\n")
            if (verbose) 
                cat(" --- Adding marker", k, "of", n.mar, "\n")
            cross <- movemarker(cross, marnam[themar[k]], 1)
            orders <- makeorders(k)
            nxo <- rep(NA, nrow(orders))
            nxo[1] <- sum(countXO(cross, 1))
            temp <- cross
            for (kk in 2:nrow(orders)) {
                temp$geno[[1]]$data <- cross$geno[[1]]$data[, 
                  orders[kk, ]]
                nxo[kk] <- sum(countXO(temp, 1))
            }
            wh <- which(nxo == min(nxo))
            if (length(wh) > 1) 
                wh <- sample(wh, 1)
            if (wh > 1) 
                cross <- simpleswitch(cross, orders[wh, ])
        }
    }
    if (use.ripple) {
		cat("cat line 10\n")
        dif <- -8
        while (dif < 0) {
			clean_up(TRUE)
            rip <- summary(ripple(cross, chr = 1, window = window, 
                verbose = FALSE))
            dif <- diff(rip[1:2, ncol(rip)])
            if (dif < 0) 
                cross <- simpleswitch(cross, rip[2, 1:(ncol(rip) - 
                  1)])
				  
        }
		cat("cat line 11\n")
    }
    match(colnames(cross$geno[[1]]$data), marnam)
}

ripple <- function (cross, chr, window = 4, method = c("countxo", "likelihood"), 
    error.prob = 1e-04, map.function = c("haldane", "kosambi", 
        "c-f", "morgan"), maxit = 4000, tol = 1e-06, sex.sp = TRUE, 
    verbose = TRUE, n.cluster = 1) 
{
    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    if (missing(chr)) {
        chr <- names(cross$geno)[1]
        warning("chr argument not provided; assuming you want chr ", 
            chr)
    }
    else {
        if (length(chr) > 1) 
            stop("ripple only works for one chromosome at a time.")
        if (!testchr(chr, names(cross$geno))) 
            stop("Chr ", chr, "not found.")
    }
    cross <- subset(cross, chr = chr)
    chr.name <- names(cross$geno)[1]
    if (nmar(cross)[1] < 3) {
        warning("Less than three markers.")
        return(NULL)
    }
    if (error.prob < 1e-50) 
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    if (window < 2) {
        warning("The window argument must be > 1; using window=2.")
        window <- 2
    }
    window <- round(window)
    method <- match.arg(method)
    map.function <- match.arg(map.function)
    n.mar <- totmar(cross)
    if (n.mar <= window) 
        orders <- ripple.perm2(n.mar)
    else {
        temp <- ripple.perm1(window)
        n <- nrow(temp)
        orders <- cbind(temp, matrix(rep((window + 1):n.mar, n), byrow = TRUE, ncol = n.mar - window))
		write.table(orders,file="tempordering.txt",sep="\t")
		cat("cat line 11\n")
        for (i in 2:(n.mar - window + 1)) {
			clean_up(TRUE)
			cat("cat line 11,5",i,"\n")
            left <- matrix(rep(1:(i - 1), n), byrow = TRUE, ncol = i - 
                1)
            if (i < n.mar - window + 1) 
                right <- matrix(rep((i + window):n.mar, n), byrow = TRUE, 
                  ncol = n.mar - window - i + 1)
            else right <- NULL
			write.table(cbind(left, temp + i - 1, right),file="tempordering.txt",sep="\t",append=T)
            #orders <- rbind(orders, cbind(left, temp + i - 1, right))
			left <- NULL
			right <- NULL
        }
		cat("cat line 12\n")
        #orders <- as.numeric(unlist(strsplit(unique(apply(orders, 1, paste, collapse = ":")), ":")))
        #orders <- matrix(orders, ncol = n.mar, byrow = TRUE)
		orders <- read.table("tempordering.txt")
    }
    n.orders <- nrow(orders)
    if (n.orders > 49) 
        print.by <- 10
    else if (n.orders > 14) 
        print.by <- 5
    else print.by <- 2
    if (method == "likelihood") {
        loglik <- 1:n.orders
        chrlen <- 1:n.orders
        m <- seq(0, by = 5, length = n.mar)
        temcross <- cross
        if (is.matrix(cross$geno[[1]]$map)) 
            temcross$geno[[1]]$map <- rbind(m, m)
        else temcross$geno[[1]]$map <- m
        if (verbose) 
            cat("  ", n.orders, "total orders\n")
		cat("cat line 12\n")
        if (n.cluster > 1 && suppressWarnings(require(snow, quietly = TRUE))) {
            if (n.orders <= n.cluster) 
                n.cluster <- n.orders
            cl <- makeCluster(n.cluster)
            clusterStopped <- FALSE
            on.exit(if (!clusterStopped) stopCluster(cl))
            if (verbose) 
                cat("   Running in", n.cluster, "clusters\n")
            clusterEvalQ(cl, require(qtl, quietly = TRUE))
            whclust <- sort(rep(1:n.cluster, ceiling(n.orders/n.cluster))[1:n.orders])
            order.split <- vector("list", n.cluster)
            for (i in 1:n.cluster) order.split[[i]] <- orders[whclust == 
                i, , drop = FALSE]
            result <- parLapply(cl, order.split, rippleSnowLik, 
                cross = temcross, error.prob = error.prob, map.function = map.function, 
                maxit = maxit, tol = tol, sex.sp = sex.sp)
            loglik <- unlist(lapply(result, function(a) a$loglik))
            chrlen <- unlist(lapply(result, function(a) a$chrlen))
        } else {
            for (i in 1:n.orders) {
                if (verbose && (i%/%print.by) * print.by == i) 
                  cat("    --Order", i, "\n")
                temcross$geno[[1]]$data <- cross$geno[[1]]$data[, 
                  orders[i, ]]
                newmap <- est.map(temcross, error.prob = error.prob, 
                  map.function = map.function, m = 0, p = 0, 
                  maxit = maxit, tol = tol, sex.sp = sex.sp, 
                  verbose = FALSE)
                loglik[i] <- attr(newmap[[1]], "loglik")
                chrlen[i] <- diff(range(newmap[[1]]))
            }
        }
		cat("cat line 13\n")
        loglik <- (loglik - loglik[1])/log(10)
        o <- order(loglik[-1], decreasing = TRUE) + 1
        orders <- cbind(orders, LOD = loglik, chrlen)[c(1, o), 
            ]
    }
    else {
        type <- class(cross)[1]
        if (type == "f2") {
            if (class(cross$geno[[1]]) == "A") 
                func <- "R_ripple_f2"
            else func <- "R_ripple_bc"
        }
        else if (type == "bc" || type == "riself" || type == 
            "risib" || type == "dh") 
            func <- "R_ripple_bc"
        else if (type == "4way") 
            func <- "R_ripple_4way"
        else if (type == "ri4self" || type == "ri8self" || type == 
            "ri4sib" || type == "ri8sib") 
            func <- "R_ripple_ril48"
        else stop("ripple not available for cross ", type)
        genodat <- cross$geno[[1]]$data
        genodat[is.na(genodat)] <- 0
        n.ind <- nind(cross)
		cat("cat line 12\n")
        if (verbose) 
            cat("  ", n.orders, "total orders\n")
        if (n.cluster > 1 && suppressWarnings(require(snow, quietly = TRUE))) {
            if (n.orders <= n.cluster) 
                n.cluster <- n.orders
            cl <- makeCluster(n.cluster)
            clusterStopped <- FALSE
            on.exit(if (!clusterStopped) stopCluster(cl))
            if (verbose) 
                cat("   Running in", n.cluster, "clusters\n")
            clusterEvalQ(cl, require(qtl, quietly = TRUE))
            whclust <- sort(rep(1:n.cluster, ceiling(n.orders/n.cluster))[1:n.orders])
            order.split <- vector("list", n.cluster)
            for (i in 1:n.cluster) order.split[[i]] <- orders[whclust == 
                i, , drop = FALSE]
            oblxo <- unlist(parLapply(cl, order.split, rippleSnowCountxo, 
                genodat = genodat, func = func))
            stopCluster(cl)
            clusterStopped <- TRUE
        }
        else {
            z <- .C(func, as.integer(n.ind), as.integer(n.mar), 
                as.integer(genodat), as.integer(n.orders), as.integer(orders - 
                  1), oblxo = as.integer(rep(0, n.orders)), as.integer(print.by), 
                PACKAGE = "qtl")
            oblxo <- z$oblxo
        }
        o <- order(oblxo[-1]) + 1
        orders <- cbind(orders, obligXO = oblxo)[c(1, o), ]
		cat("cat line 13\n")
    }
    rownames(orders) <- c("Initial", paste(1:(nrow(orders) - 
        1)))
    class(orders) <- c("ripple", "matrix")
    attr(orders, "chr") <- chr.name
    attr(orders, "window") <- window
    attr(orders, "error.prob") <- error.prob
    attr(orders, "method") <- method
	cat("cat line 14\n")
    orders[, 1:n.mar] <- t(apply(orders[, 1:n.mar, drop = FALSE], 
        1, function(a) {
            n <- length(a)
            if ((1:n)[a == 1] > (1:n)[a == n]) 
                return(rev(a))
            else return(a)
        }))
    orders
}

clean_up <- function(verbose=FALSE){
  p_usage <- gc()[2, 3]
  s_usage <- p_usage
  n_usage <- gc()[2, 3]
  while (n_usage < p_usage) {
    p_usage = n_usage
    n_usage <- gc()[2, 3]
  }
  if (verbose) cat("clean_up: ", s_usage, " ", n_usage, "\n")
}


setwd("C:/xbinary - Konrad")
cross <- read.cross(file="mycross2.csv",format="csvr",genotypes=c("AA","AB"))
cross <- orderMarkers(cross,chr=1,window=25)