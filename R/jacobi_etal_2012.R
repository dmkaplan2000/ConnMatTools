######################################################################
# Algorithm for defining subclusters of highly-connected sites from
# the following publication:
#
# Jacobi, M. N., André, C., Döös, K., and Jonsson,
# P. R. 2012. Identification of subpopulations from connectivity
# matrices. Ecography, 35: 1004–1016.
#
# NOTE: In this paper, the connectivity matrix is oriented so that
# C_ij is dispersal from i to j, whereas in this code, the
# connectivity matrix is oriented so that C_ij is dispersal from j to
# i.  This choice of orientation is arbitrary, but one must always be
# consistent.
######################################################################

# The slit function takes a vector of sites and returns a list of one
# or two vectors of sites. It also accepts an input parameter defining
# how many times to try the spilt (the best result from the tries is
# returned).
splitCM <- function(indices,pp,tries,beta,
                    th=1e-10,alpha=0.1,maxIt=500,best=1e10) {
    # makes a submatrix of the total connectivity matrix only
    # involving the sites in the index list
    ppp = pp[indices, indices, drop=F]

    # al appears to be just 1/beta
    al = 1/beta
    
    n = dim(ppp)[2];
    s = matrix(1,nrow=n)
    
    eNoSplit = -t(s) %*% ppp %*% s + al * sum(s)^2;

    for (kk in 1:tries) {
        s = matrix( runif(n,min=-1,max=1), ncol = 1 )
        ds = s + 1
        sTot = sum(s)
        sOld = sign(s)

        for (it in 1:maxIt) {
            if (norm(s-ds,"2") <= th)
                break

            # the following three lines implements EQ. 8 in the paper
            v = ppp %*% s - al * sum(s)            
            ds = sqrt(abs(v)) * sign(v)
            s = alpha * s + (1-alpha) * ds
        }

        if (it == maxIt) warning("Reached max iterations")

        s = sign(s)
        e = -t(s) %*% ppp %*% s + al * sum(s)^2;
        # calclulates the value of the cost function

        if (e < best) {
            sBest = s
            best = e
        }
    }

    part = indices[ sBest == -1 ]
    notPart = setdiff( indices, part )

    if ( (best < eNoSplit) & (length(part)>0) & (length(notPart)>0) ) {
        return(list(part,notPart))
    }

    return(list(indices))
}

# This funtion splits the index list recursively until none of the
# subpopulations can be split further to improve the minimization
recSplitCM <- function(ta, pp, tries, beta, ...) {
    old = 0
    ret = ta
    while (length(ret) > old) {
        old = length(ret)
        ret = unlist( sapply( ret, splitCM, pp=pp, tries=tries, beta=beta,
            ..., simplify=F),
            recursive=F )
    }
    return(ret)
}

# This function tries to merge random subopoulations, checking if the
# result is a better soluton to the minimization problem
mergeCM <- function ( ta,  pp, beta ) {
    nIt = length(ta)^2
    ret = ta
    al = 1/beta
    
    for (it in 1:nIt) {
        if (length(ret) < 2) break # Must have at least 2 to merge
        
        ii = sort(sample( 1:length(ret), 2 ))
        li = c( ret[[ ii[1] ]], ret[[ ii[2] ]] )
        pTest = pp[li,li]
        
        s = matrix(1,nrow=length(li))
        eTogether = -t(s) %*% pTest %*% s + al * sum(s)^2
        s[1:length(ret[[ ii[1] ]])] = -1
        eSplit = -t(s) %*% pTest %*% s + al * sum(s)^2

        if (eTogether<eSplit) {
            ret[[ ii[1] ]] = li
            ret = ret[ -ii[2] ]
        }
    }

    return(ret)
}

# The quality function is not used to evaluate the result of the
# solution, i.e. it measures the average total leakage between the
# subpopulations
#
# This quality measure is equal to 1 - mean(RLR) of the reduced
# connectivity matrix, where RLR=relative local retention, i.e., the
# fraction of settling individuals that originated at their site of
# settlement.
qualCM <- function( ta, p ) {
    pii = matrix( 0.0, nrow=dim(p)[2], ncol=length(ta) )

    for (kk in 1:length(ta)) {
        pii[ ta[[kk]], kk ] = 1
    }

    # Note I use p instead of t(p), as was in Jacobi code.
    # This is because my matrices are oriented columns to rows
    #pt = t(pii) %*% p %*% pii %*% solve( t(pii) %*% pii )

    # Somewhat quicker method
    pt =  t(pii) %*% p %*% pii #%*% diag( 1 / sapply(ta,length) )
    # No need to normalize by number of sites in cluster of larval
    # origin because we are just going to normalize each colum
    
    ss = apply( pt, 2, sum )
    ss[ ss == 0 ] = 1
    
    return(1 - mean(diag(pt) / ss))
}

# Helper function to compute a set of beta values using formula used
# in Jacobi et al. 2012.
jacobi.et.al.2012.betas.vector <- function(n,steps=10,cycles=3/4,
                                           coeff=0.8,pwr=3.0)
  n/(1+coeff*sin(seq(0,cycles*2*pi,length.out=steps)))^pwr


# Acutal Jacobi algorithm
#
# Note that I normalize columns, not rows.  This will make "larval
# loss" uniform for all points of origin, presuming the connectivity
# matrix, p, is oriented so that settlers = p %*% eggs (i.e., each
# column represents dispersal from a specific site of origin, each row
# represents dispersal to a specific site of settlement).
jacobi.algo <- function(p, normalize.cols=TRUE,
                        make.symmetric="mean", remove.diagonal=TRUE, 
                        cycles = 2, tries=5, steps=10,
                        betas=jacobi.et.al.2012.betas.vector(dim(p)[2],steps),
                        ... ) {
    pp <- p

    # Make larval loss uniform over space
    if (normalize.cols)
        for (kk in 1:dim(pp)[2]) {
            ss = sum(pp[,kk])
            if (ss>0) pp[,kk] = pp[,kk] / ss
        }

    # Force symmetric if not already the case
    mymax = function(x) { ii = x < t(x); x[ii] = t(x)[ii]; return(x) }
    pp = switch(make.symmetric,
        mean = (pp + t(pp)) / 2.0,
        max = mymax(pp),
        stop("Bad max.symmetric string"))

    # Not sure if this is obligatory or optional
    # I believe it should be optional
    if (remove.diagonal)
        diag(pp) <- 0

    clusters = matrix(NA,nrow=dim(p)[2],ncol=cycles*length(betas))
    num.clusters = rep(NA,cycles*length(betas))
    qualities = num.clusters
    for (kk in 1:cycles) {
        print( paste( "Starting cycle",kk ) )

        # Initialize with one big cluster
        ta = list(1:dim(p)[2])

        # Loop over betas
        for (ll in 1:length(betas)) {
            beta = betas[ll]

            print( paste( "beta =", beta ) )
            
            taOld = list()
            while (!setequal(ta,taOld)) {
                taOld = ta
                ta = recSplitCM(ta, pp, tries, beta, ...)
                #if (length(unlist(ta))!=dim(pp)[2]) stop("recSplitCM error")
                ta = mergeCM(ta,  pp, beta )
                #if (length(unlist(ta))!=dim(pp)[2]) stop("mergeCM error")
            }

            # Store results
            nn = (kk-1)*length(betas)+ll
            qualities[nn] = qualCM(ta,p)
            num.clusters[nn] = length(ta)
            for (mm in 1:length(ta))
                clusters[ta[[mm]],nn] = mm
        }
    }

    # For each number of clusters, find the configuration with the
    # best quality measure
    imin <- function(x) { which( min(x) == x )[1] }
    best.clusters = by( data.frame(quality=qualities,
        index=1:length(qualities)),
        INDICES=num.clusters,
        FUN=function(d) d[imin(d$quality),] )
    
    wareturn(list(betas=rep(betas,cycles),
                clusters=clusters,qualities=qualities,
                num.clusters=num.clusters,best.clusters=best.clusters))
}

# A helper function to convert a vector of cluster identifications into
# a list appropriate for qualCM, recSplitCM, etc.
clusters.vector.to.list <- function(x) {
    xx = sort(unique(x))

    ta = list()
    for (yy in xx)
        ta[[length(ta)+1]] = which( x == yy )
    return(ta)
}
