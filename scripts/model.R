###### model: learn elastic net model on training data 
######---------Input: trainX, trainY
######---------Return: selected features and coefficents

elastic.net <- function(trainX,trainY){	
	if(nrow(trainX)!=length(trainY)){
		stop("Number of observations is differerent")
	} 
	
	# optimize alpha---mixing parameter  
	a <- seq(0, 1, 0.1)
    search <- foreach(ai = a, .combine = rbind) %dopar% {
        cv.fit <- cv.glmnet(
			trainX,
			trainY,
			nfold = 10,
			type.measure = "mse",
			paralle = TRUE,
			alpha = ai
			)
        data.frame(
			cvm = min(cv.fit$cvm),
			lambda = cv.fit$lambda.min,
			alpha = ai
			)
	} 
    cv.opt <- search[search$cvm == min(search$cvm),] 
	
	# fit model by optimized alpha and lambda
	yfit = glmnet(
        trainX,
        trainY,
        lambda = cv.opt$lambda,
        alpha = cv.opt$alpha
		)       
	idf <- coef(yfit)
	idx <- which(idf != 0)
	selectf <- data.frame(
		features = idf@Dimnames[[1]][idx], 
		coefs = idf [idx]
	)	
}

MWAS <- function(gwas, weight, geno){
	z <- gwas %*% weight
	z.cor <- cor(geno)
	se <- sqrt(weight %*%  z.cor %*%  weight)	
	z <- z/se
	p=pnorm(abs(z),lower.tail=F)*2	
	return(c(z, p))
}