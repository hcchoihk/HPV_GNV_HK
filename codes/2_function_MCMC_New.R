# MCMC function
# main changes:
# 1) set fun_logLikelihood as an input for function_MCMC;
# 2) put all arugments as extInput_logL --set ... as the arguments for fun_logLikelihood--
# 3) allow input step size
# note: x-parameter has to be the first argument in fun_logLikelihood
#
# logPrior = function(x)
# LB, UB, minProbAccept, maxProbAccept are _vectors_ in length of the number of x-parameter


# function_MCMC_Block , allow block updating


### MCMC, single-updating
function_MCMC <- function(fun_logLikelihood, LB, UB, startingPoint, numStepsPerParameter, minProbAccept, maxProbAccept, stepSTD, logPrior, extInput_logL=NULL, reflective_update_YN = TRUE, print_updating_YN = FALSE, randomSeed=NULL)
{
	# random seed
	if (!is.null(randomSeed)){ set.seed(randomSeed)}

	numParameter = length(startingPoint);
	
	minIterBeforeRestart = 1000*numParameter;  
	
	if (length(minProbAccept)==1){
		minProbAccept = rep(minProbAccept, numParameter)
	}
	if (length(maxProbAccept)==1){
		maxProbAccept = rep(maxProbAccept, numParameter)
	}
	
	
	numStepsBetweenDisplay = max(ceiling(0.02 * numStepsPerParameter), 100);
	storageSize = min(numStepsPerParameter,1000);
	probAccept = (minProbAccept+maxProbAccept)/2;
	oldProbAccept = array(-1,c(1,numParameter));
	
	logL_wrongValue = -1e10;
	
	
	if (missing(logPrior)){
		logPrior = function(x){ return(0)}
	}
	

	if (missing(stepSTD)){
		stepSTD = rep(0.1,length(LB));
	}
	delta = 30;

	chainRecord = array(0,c(numStepsPerParameter*numParameter+1, numParameter+2)); # one row for each step; each row = +2 for value of logLikelihood and value of acceptance rate)
	chainOutcome = rep(list(NULL), numStepsPerParameter*numParameter+1);
	chainGood = -1;

	bestParameters = startingPoint;
	likelihood_out = fun_logLikelihood(bestParameters, extInput_logL);
	bestLogLikelihood = likelihood_out$logL
	bestOutcome = likelihood_out$outcome
	
	while (chainGood==-1)
	{
    	x = bestParameters;
		likelihood_out = fun_logLikelihood(x, extInput_logL);
    	xLogLikelihood = likelihood_out$logL
		xOutcome = likelihood_out$outcome;
    	xTarget = xLogLikelihood+logPrior(x);
    	
    	chainRecord[1,] = c(x, xLogLikelihood,1);
		chainOutcome[[1]] = xOutcome;
		numAccept = array(1,c(1,length(x)));
		numSteps = array(1,c(1,length(x)));
		
		stepSizeStorage = array(0,c(storageSize, numParameter)); # because the stepSTD is not changed in each new chain, generate a 1000 stepSize only to save space
    
    	for (para in 1:numParameter) # adjust step size when the mixing is not good 
    	{
  			if (probAccept[para]<minProbAccept[para])
  			{
  			  p_delta = minProbAccept[para] - probAccept[para];
  			  scalingFactor = 1/(1+p_delta*delta);
  				stepSTD[para] = scalingFactor*stepSTD[para]; 
  			}
  			if (probAccept[para]>maxProbAccept[para])
  			{
  			  p_delta = probAccept[para] - maxProbAccept[para];
  			  scalingFactor = 1+p_delta*delta;
  				stepSTD[para] = scalingFactor*stepSTD[para];
  			}
	  	}
    	
		time0 <- Sys.time();
		totalSecElapsed = 0;
		probStorage = runif(storageSize*length(x), 0, 1);
		probStorage = matrix(probStorage, ncol=length(x));
		log_probStorage = log(probStorage)
		for (para in 1:numParameter)
		{
			stepSizeStorage[,para] = rnorm(storageSize, 0, stepSTD[para]);
		}
		numTotalSteps = 1;

		paraUpdate = 1; # will update paraUpdate inside for-iter, so no need to set for-loop to paraUpdate

	
		for (iter in 2:(numStepsPerParameter*numParameter+1))
		{
			if (numSteps[paraUpdate]%%storageSize==0)
			{
				stepSizeStorage[,paraUpdate] = rnorm(storageSize, 0, stepSTD[paraUpdate]);
			}
			if (iter%%storageSize==0)
			{
				probStorage = runif(storageSize*length(x), 0, 1);
				probStorage = matrix(probStorage, ncol=length(x));
				log_probStorage = log(probStorage)
			}
        
        
			# next step in the random walk
			y = x;
			
			# Update one element (e.g. Mu or Sigma) in the parameter vector every time
			y[paraUpdate] = x[paraUpdate]+stepSizeStorage[1+(numSteps[paraUpdate]-1)%%storageSize, paraUpdate];

		if (reflective_update_YN){
			if (reflective_update_YN==1){
				# reflect once only, so use "if"
				if (y[paraUpdate]<LB[paraUpdate] || y[paraUpdate]>UB[paraUpdate])
				{
					if (y[paraUpdate]<LB[paraUpdate])
					{
						y[paraUpdate] = LB[paraUpdate]+(LB[paraUpdate]-y[paraUpdate]);
					}else
					{
						y[paraUpdate] = UB[paraUpdate]-(y[paraUpdate]-UB[paraUpdate]);
					}
				}
			} else{
				# reflect multiple times, so use "while"
				while (y[paraUpdate]<LB[paraUpdate] || y[paraUpdate]>UB[paraUpdate])
				{
					if (y[paraUpdate]<LB[paraUpdate])
					{
						y[paraUpdate] = LB[paraUpdate]+(LB[paraUpdate]-y[paraUpdate]);
					}else
					{
						y[paraUpdate] = UB[paraUpdate]-(y[paraUpdate]-UB[paraUpdate]);
					}
				}
			}
		} # reflective_update_YN
		
		
			# Compute new target likelihood
			likelihood_out = fun_logLikelihood(y, extInput_logL);
			yLogLikelihood = likelihood_out$logL;
			if (yLogLikelihood<=logL_wrongValue){
				yOutcome = xOutcome;
			} else{
				yOutcome = likelihood_out$outcome;
			}
			yTarget = yLogLikelihood+logPrior(y);
			# Metropolis-Hastings algorithm
			# if compared to log_probStorage, no need to exp every iteration
			# if (probStorage[1+((iter-1)%%storageSize)] < min(1,exp(yTarget-xTarget)))
			if (log_probStorage[1+((iter-1)%%storageSize), paraUpdate[1]] < min(0,(yTarget-xTarget)))
			{ 
				x = y;
				xLogLikelihood = yLogLikelihood;
				xOutcome = yOutcome;
				xTarget = yTarget;
				numAccept[paraUpdate] = numAccept[paraUpdate]+1;
			}
			
			if (all(xLogLikelihood>bestLogLikelihood))
			{
				bestLogLikelihood = xLogLikelihood;
				bestOutcome = xOutcome;
				bestParameters = x;
			}
			numSteps[paraUpdate] = numSteps[paraUpdate]+1;
			
			chainRecord[iter,] = c(x,xLogLikelihood,probAccept[paraUpdate]);
			chainOutcome[[iter]] = xOutcome;
        	
			probAccept[paraUpdate] = numAccept[paraUpdate]/numSteps[paraUpdate];
			
			if (probAccept[paraUpdate]>minProbAccept[paraUpdate] && probAccept[paraUpdate]<=maxProbAccept[paraUpdate])
			{
				chainGood = 1;
			} else
			{
				chainGood = -1;
			
				if (iter>minIterBeforeRestart)
				{
					if (print_updating_YN){
					
					print('');
					print('Restarting MCMC because P(Accept) is out of range...');
					cat('Parameter ', toString(paraUpdate), ': no. of steps = ', toString(numSteps[paraUpdate]), ', P(Accept) = ', toString(probAccept[paraUpdate]), ', stepSTD = ', toString(stepSTD[paraUpdate]), "\n");
					print('Outcome of the last iteration.');
					print(chainOutcome[[iter]]);

					
					Sys.sleep(1);
					break; # break for-iter
					
					} # print_updating_YN
				}
			}
			
			paraUpdate = paraUpdate+1;
			if (paraUpdate>numParameter)
			{
				paraUpdate = 1; # go back to update 1st para
			}
			
			totalSecElapsed = totalSecElapsed + difftime(Sys.time(),time0, units="sec");
			
			time0 <- Sys.time();
			
			if (paraUpdate==1 && numSteps[paraUpdate]%%numStepsBetweenDisplay==0)
			{
				if (print_updating_YN){
					cat('MCMC ', toString(round(numSteps[paraUpdate]/numStepsPerParameter*100,digits=2)),'% complete. Time elapsed = ', toString(round(totalSecElapsed/60,digits=2)),' min. Estimated remaining time = ',toString(round(totalSecElapsed/numSteps[paraUpdate]*(numStepsPerParameter-numSteps[paraUpdate])/60,digits=2)),' min.', "\n");
				} # print_updating_YN
			}
		} # end for-iter
    } # while chainGood==-1
	
	
	# -1 for removing  startingPoint
	chainLikelihood = chainRecord[-1, numParameter+1]
	chainAcceptRate = chainRecord[-1, numParameter+2]
	chainRecord = chainRecord[-1, 1:numParameter]
	chainOutcome = do.call(rbind, chainOutcome[-1]);
	
	return(list("chainRecord"=chainRecord, "probAccept"=probAccept, "stepSTD"=stepSTD, "chainLikelihood"=chainLikelihood, "chainAcceptRate"=chainAcceptRate, "chainOutcome"=chainOutcome));
} # end function_MCMC





### MCMC with block updating
function_MCMC_block <- function(fun_logLikelihood, LB, UB, startingPoint, numStepsPerParameter, minProbAccept, maxProbAccept, stepSTD, logPrior, block=NULL, extInput_logL=NULL, reflective_update_YN = TRUE, print_updating_YN = FALSE, randomSeed=NULL, startingPoint_directuse_YN=FALSE)
{
	# block, allow block updating, a list to contain the elements
	
	# random seed
	if (!is.null(randomSeed)){ set.seed(randomSeed)}
	
	
	numParameter = length(startingPoint);
	
	if (is.null(block)){
		block = as.list(1:numParameter);
	}
	numParam_block = sapply(block, length)
	num_block = length(block);

	
	
	minIterBeforeRestart = 1000*numParameter;  
	
	if (length(minProbAccept)==1){
		minProbAccept = rep(minProbAccept, numParameter)
	}
	if (length(maxProbAccept)==1){
		maxProbAccept = rep(maxProbAccept, numParameter)
	}
	
	
	numStepsBetweenDisplay = max(ceiling(0.02 * numStepsPerParameter), 100);
	storageSize = min(numStepsPerParameter,1000);
	probAccept = (minProbAccept+maxProbAccept)/2;
	oldProbAccept = array(-1,c(1,numParameter));
	
	if (missing(logPrior)){
		logPrior = function(x){ return(0)}
	}
	

	if (missing(stepSTD)){
		stepSTD = rep(0.1,length(LB));
	}
	
	delta = 30; # for adjusting the stepsize

	chainRecord = array(0,c(numStepsPerParameter*numParameter+1, numParameter+2)); # one row for each step, row=1 for initial; each row = +2 for value of logLikelihood and value of acceptance rate)
	chainOutcome = rep(list(NULL), numStepsPerParameter*numParameter+1);
	chainGood = -1;

	# Test startingPoint
	startingGood = -1;
	if (startingPoint_directuse_YN){
		bestParameters = startingPoint;
	} else{
		bestParameters = startingPoint * runif(length(startingPoint), min=0.8,max=1.2);
	}
	
	while (startingGood==-1){
	
		likelihood_out = fun_logLikelihood(bestParameters, extInput_logL);
		bestLogLikelihood = likelihood_out$logL
		bestOutcome = likelihood_out$outcome
		if (is.null(bestOutcome)){
			print("The outcome from the startingPoint is NULL. Re-sample staringPoint.");
			width_startPt = 0.2
			bestParameters = startingPoint * runif(length(startingPoint), min=1-width_startPt, max=1+width_startPt);
		} else {
			startingGood = 1;
		}
		
	} # test startingPoint
	
	if (!is.null(randomSeed)){ set.seed(randomSeed)}
	

	while (chainGood==-1)
	{

    	x = bestParameters;
		likelihood_out = fun_logLikelihood(x, extInput_logL);
    	xLogLikelihood = likelihood_out$logL
		xOutcome = likelihood_out$outcome;
    	xTarget = xLogLikelihood+logPrior(x);

		print("bestParameters")
		print(x)
		print("outcome")
		print(xOutcome)
		print("")


    	chainRecord[1,] = c(x, xLogLikelihood,1);
		chainOutcome[[1]] = xOutcome;
		numAccept = array(1,c(1,length(x)));
		numSteps = array(1,c(1,length(x)));
		
		stepSizeStorage = array(0,c(storageSize, numParameter)); # because the stepSTD is not changed in each new chain, generate a 1000 stepSize only to save space

    	if ( any(probAccept<minProbAccept | probAccept>maxProbAccept) ){
			for (para in which(probAccept<minProbAccept | probAccept>maxProbAccept)) 
			{
				if (probAccept[para]<minProbAccept[para])
				{
					p_delta = minProbAccept[para] - probAccept[para];
					scalingFactor = 1/(1+p_delta*delta);
					stepSTD[para] = scalingFactor*stepSTD[para]; 
				}
				if (probAccept[para]>maxProbAccept[para])
				{
					p_delta = probAccept[para] - maxProbAccept[para];
					scalingFactor = 1+p_delta*delta;
					stepSTD[para] = scalingFactor*stepSTD[para];
				}
			} # adjust step size when the mixing is not good
			for (iblock in 1:num_block){ # assign variability for 
				if (numParam_block[iblock]>1){
					stepSTD[block[[iblock]]] = stepSTD[block[[iblock]]] * runif(numParam_block[iblock], 0.85,1.15)
				}
			}
		}
		
    	
		time0 <- Sys.time();
		totalSecElapsed = 0;
		probStorage = runif(storageSize*length(x), 0, 1);
		probStorage = matrix(probStorage, ncol=length(x));
		log_probStorage = log(probStorage);
		for (para in 1:numParameter)
		{
			stepSizeStorage[,para] = rnorm(storageSize, 0, stepSTD[para]);
		}
		numTotalSteps = 1;

		blockUpdate = 1;

	
		iter = 2;
		# for (iter in 2:(numStepsPerParameter*numParameter+1))
		while (iter<=(numStepsPerParameter*numParameter+1))
		{
			
			paraUpdate = block[[blockUpdate]];
			paraUpdate_num = numParam_block[blockUpdate]; # number of parameters to be updated
			iter_Update = iter-1 + (1:paraUpdate_num);
		
			
			currentStep_paraUpdate = numSteps[paraUpdate[1]] # numSteps is the same for all para in the same block
			if (currentStep_paraUpdate%%storageSize==0) # draw another stepSizeStorage
			{
				stepSizeStorage[,paraUpdate] = sapply(stepSTD[paraUpdate], function(x) rnorm(storageSize, 0, x)); #rnorm(storageSize, 0, stepSTD[paraUpdate]);
			}
			if (iter%%storageSize==0) # draw new probStorage
			{
				probStorage = runif(storageSize*length(x), 0, 1);
				probStorage = matrix(probStorage, ncol=length(x));
				log_probStorage = log(probStorage)
			}
        
        
			# next step in the random walk
			y = x;
			
			# Update one block every time
			y[paraUpdate] = x[paraUpdate]+stepSizeStorage[1+(currentStep_paraUpdate-1)%%storageSize, paraUpdate];

		if (reflective_update_YN){
			while (any(y[paraUpdate]<LB[paraUpdate]) || any(y[paraUpdate]>UB[paraUpdate]))
			{
				# print( ifelse(reflective_update_YN==1, "Reflect once only", "Reflect until within LB/UB") )
				if (any(y[paraUpdate]<LB[paraUpdate]))
				{
					paraTemp = paraUpdate[which(y[paraUpdate]<LB[paraUpdate])]
					y[paraTemp] = LB[paraTemp]+(LB[paraTemp]-y[paraTemp]);
				} else if (any(y[paraUpdate]>UB[paraUpdate]))
				{
					paraTemp = paraUpdate[which(y[paraUpdate]>UB[paraUpdate])]
					y[paraTemp] = UB[paraTemp]-(y[paraTemp]-UB[paraTemp]);
				}
				if (reflective_update_YN==1){ break; } # update once only
			}
		} # reflective_update_YN			


			# Compute new target likelihood
			likelihood_out = fun_logLikelihood(y, extInput_logL);
			yLogLikelihood = likelihood_out$logL;
			yOutcome = likelihood_out$outcome;
			yTarget = yLogLikelihood+logPrior(y);
			# Metropolis-Hastings algorithm
			# if compared to log_probStorage, no need to exp every iteration
			# if (probStorage[1+((iter-1)%%storageSize)] < min(1,exp(yTarget-xTarget)))
			# print(likelihood_out)
			if (log_probStorage[1+((iter-1)%%storageSize), paraUpdate[1]] < min(0,(yTarget-xTarget)))
			{ 
				x = y;
				xLogLikelihood = yLogLikelihood;
				xOutcome = yOutcome;
				xTarget = yTarget;
				numAccept[paraUpdate] = numAccept[paraUpdate]+1;
			}
			
			if (all(xLogLikelihood>bestLogLikelihood))
			{
				bestLogLikelihood = xLogLikelihood;
				bestOutcome = xOutcome;
				bestParameters = x;
			}
			numSteps[paraUpdate] = numSteps[paraUpdate]+1;
		
			# chainRecord[iter,] = c(x,xLogLikelihood,probAccept[paraUpdate]);
			chainRecord[iter_Update,] = cbind(matrix(c(x,xLogLikelihood), ncol=numParameter+1, nrow=paraUpdate_num, byrow=TRUE), probAccept[paraUpdate]);
			if (paraUpdate_num==1){
				chainOutcome[[iter_Update]] = xOutcome
			} else{
				chainOutcome[iter_Update] = rep(list(xOutcome), paraUpdate_num);
			} 				
        	
			probAccept[paraUpdate] = numAccept[paraUpdate]/numSteps[paraUpdate];

			if (all(probAccept[paraUpdate]>minProbAccept[paraUpdate]) && all(probAccept[paraUpdate]<=maxProbAccept[paraUpdate]))
			{
				chainGood = 1;
			} else
			{
				chainGood = -1;
			
				if (iter>minIterBeforeRestart)
				{
					if (print_updating_YN){
					
					paraTemp = paraUpdate[1]
					print('');
					print('Restarting MCMC because P(Accept) is out of range...');
					cat( sprintf('Block: %d, parameter: %s\n', blockUpdate, paste(paraUpdate, collapse=",")) )
					cat( sprintf('Parameter %d : no. of steps = %d, P(Accept) = %.3f, stepSTD = %s\n', paraTemp, numSteps[paraTemp], probAccept[paraTemp], toString(stepSTD[paraTemp])));
					
					Sys.sleep(1);
					break; # break for-iter
					
					} # print_updating_YN
				}
			}
			# determine chainGood
			
			# finish update
			
			totalSecElapsed = totalSecElapsed + difftime(Sys.time(),time0, units="sec");
			
			time0 <- Sys.time();
			
			
			blockUpdate = blockUpdate + 1;
			if (blockUpdate > num_block){
				blockUpdate = 1;
			}
			
			paraTemp = block[[blockUpdate]][1];
			if (blockUpdate==1 && numSteps[paraTemp]%%numStepsBetweenDisplay==0)
			{
				if (print_updating_YN){
					cat('MCMC ', toString(round(numSteps[paraTemp]/numStepsPerParameter*100,digits=2)),'% complete. Time elapsed = ', toString(round(totalSecElapsed/60,digits=2)),' min. Estimated remaining time = ',toString(round(totalSecElapsed/numSteps[paraUpdate]*(numStepsPerParameter-numSteps[paraUpdate])/60,digits=2)[1]),' min.', "\n");
				} # print_updating_YN
			}
			
			
			
			iter = iter + paraUpdate_num
		} # end for-iter
    } # while chainGood==-1
	
	#return(chainRecord,probAccept);

	# -1 for removing startingPoint
	chainLikelihood = chainRecord[-1, numParameter+1]
	chainAcceptRate = chainRecord[-1, numParameter+2]
	chainRecord = chainRecord[-1, 1:numParameter]
	chainOutcome = do.call(rbind, chainOutcome[-1]);
	
	return(list("chainRecord"=chainRecord, "probAccept"=probAccept, "stepSTD"=stepSTD, "chainLikelihood"=chainLikelihood, "chainAcceptRate"=chainAcceptRate, 'chainOutcome'=chainOutcome));
} # end function_MCMC_block


