#' Function to select optimal network of protected areas based on connectivity
#' 
#' This function finds the optimal network of protected areas based on
#' connectivity using the eigenvalue perturbation approach described
#' in Nilsson Jacobi & Jonsson (2011).
#' 
#' @param conn.mat a square connectivity matrix.
#' @param k number of calculated eigenvalues and associated
#' eigenvectors. Currently unused, since all the eigenvalues and the
#' corresponding eigenvectors of the matrix are calculated.
#' @param delta the effect of protecting site i (e.g. increase in
#' survival or fecundity in protected areas relative to unprotected
#' areas). Now a single value, in future it will be possible to
#' specify site-specific values. The perturbation theory used in the
#' construction of the algorithm assumes delta to be small
#' (e.g. delta=0.1). However, higher values give also good results.
#' @param theta the threshold of donor times recipient value that a
#' site must have to be selected.
#' @param M the maximal number of sites selected from each
#' subpopulation even if there are more sites above the threshold
#' theta
#' @param epsilon.lambda Threshold for removing complex eigenvalues.
#' @param epsilon.uv Threshold for removing eigenvectors with elements
#' of opposite signs of comparable magnitude.
#' @param only.list Logical, whether the function return only the list
#' of selected sites or also the predicted impact of each selected
#' site on the eigenvalues
#'
#' @return A list.
#'
#' @references Jacobi, M. N., and Jonsson, P. R. 2011. Optimal
#' networks of nature reserves can be found through eigenvalue
#' perturbation theory of the connectivity matrix. Ecological
#' Applications, 21: 1861–1870.
#' 
#' @author
#' Marco Andrello \email{marco.andrello@gmail.com}
#' @encoding UTF-8
#' @export
ProtectedAreaSelection <- function(d,k=dim(d)[1],delta=0.1,theta=0.05,M=20,epsilon.lambda=0.0001,epsilon.uv=0.05,only.list=T) {

    # Matrix size
    n = dim(d)[1]

    # Calculates eigenvalues and eigenvectors
    # Left and right eigenvalues must be sorted from largest to smallest; eigen() does this.
    # Currently, the calculation is very slow and can give memory problems for large matrices
    # Using arpack do not solve the problem...
    # The algorithm used by Matlab is much faster and makes use of sparse matrices; calculates the largest k eigenvalues
    # Do not know if something similar exists in R
    r <- eigen(d)
    l <- eigen(t(d))
    lambda_r 	<- r$values
    u		<- r$vectors
    v		<- l$vectors
    rm(r,l)

    # Removing complex eigenvalues and the associated eigenvectors
    index_real <- which(abs(Im(lambda_r)) < epsilon.lambda)
    lambda_r <- Re(lambda_r[index_real])
    u <- Re(u[,index_real])
    v <- Re(v[,index_real])
    k <- length(lambda_r)

    # Normalisation
    for (i in 1 : k) {
	v[,i] = v[,i] / (v[,i] %*% u[,i])
	u[,i] = u[,i] / (v[,i] %*% u[,i])
    }

    # Calculate donor times recipient value
    donorRecipient <- u * v

    # Remove vectors pairs whose donorRecipient vector has negative and positive elements of comparable magnitude

    index_pos <- array(NA,dim=k)

    for (i in 1 : k) {
        maxx <- max(donorRecipient[,i])
        minn <- min(donorRecipient[,i])

        if (maxx >= 0 && minn >= 0) {
            index_pos[i] <- 1
            next
        }

        if (maxx < 0 && minn < 0) {
            donorRecipient[,i] <- donorRecipient[,i] * -1
            index_pos[i] <- 1
            next
        }

        if ( abs(maxx/minn) > 1/epsilon.uv) {
            index_pos[i] <- 1
            next
        }

        if ( abs(maxx/minn) < epsilon.uv) {
            donorRecipient[,i] <- donorRecipient[,i] * -1
            index_pos[i] <- 1
            next
        }

        index_pos[i] <- 0

    }

    id_pos 		<- which(index_pos == 1)

    lambda_r 		<- lambda_r[id_pos]
    u 			<- u[,id_pos]
    v 			<- v[,id_pos]
    donorRecipient 	<- donorRecipient[,id_pos]
    k 			<- length(lambda_r)


    # Sorting sites according to their donorRecipient values, corresponding site index in 'sites'
    sites 		<- array(dim=c(n,k))
    valueForSites 	<- array(dim=c(n,k))
    for (i in 1 : k) {
        a <- sort(donorRecipient[,i],decreasing=T,index.return=T)
        sites[,i] <- a$ix
        valueForSites[,i] <- a$x
    }

    # Ensuring that only the top M sites are selected (and not more)
    modifiedValuesForSites 			<- valueForSites
    modifiedValuesForSites[(M+1):n,] 	<- theta-1000 # -1000 is less than theta


    # Predicting impact on the new eigenvalues
    zeroedValuesForSites <- modifiedValuesForSites
    zeroedValuesForSites[zeroedValuesForSites<theta] <- 0
    predictedValues <- lambda_r * (1 + delta * colSums(zeroedValuesForSites))

    # Creating the final list ordered by the predicted eigenvalues
    listOrder <- sort(predictedValues,decreasing=T,index.return=T)$ix
    finalList <- vector()
    for (i in 1 : k) {
        finalList <- c(finalList, sites[modifiedValuesForSites[,listOrder[i]]>theta, listOrder[i]])
    }

    finalList <- unique(finalList)

    if (only.list) {
        return(finalList)
    } else {
        return(list(finalList=finalList,predictedValues=predictedValues))
    }

}

