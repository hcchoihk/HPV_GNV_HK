## functions for fit_CxIncid_validate

## get mid-points given a vector of end-points
get_midpt = function(endpt){
	return( (head(endpt,-1) + tail(endpt,-1))/2 );
} # get_midpt

# whether a point is between another two points
is_between = function(pt2, pt1, pt3){
	return( (pt1 <= pt2 & pt2 <= pt3) | (pt1 >= pt2 & pt2 >= pt3) )
} # is_between


## log of poisson distribution
# X: observed events; lambda: rate of events based on the model
# poisspdf(X, lambda) = exp(-lambda) * lamba^X / X!
# logpoisspdf = -lambda + X*log(lambda) - log(X!)
# X! = gamma(X+1)

logpoisspdf_noC = function(X, lambda){
	# exclude the constant term (X!) for the purpose of comparing (log-)likelihood
	return( -lambda + X*log(lambda) );	
} # logpoisspdf_noC

logpoisspdf_C = function(X){
	# the constant for poisson distribution
	# log(constant) = log(1/X!) = -log(X!)	
	return( -sum(lgamma(X+1)) )
} # logpoisspdf_C



## log of binomial distribution
# n = trials, k = number of event, p = event rate 
# binompdf(k, n, p) = Comb_n_k * p^k * (1-p)^(n-k); Comb_n_k = number of combinations
# logbinompdf = k*log(p) + (n-k)*log(1-p) + log(C_n_k)

logbinompdf_noC = function(n, k, p){
	# exclude the constant term
	return( k*log(p) + (n-k)*log(1-p) )
} # logbinompdf_noC

logbinompdf_C = function(n, k){
	# the constant term only
	return( lchoose(n,k) )
} # logbinompdf_C, constant only


## log of multinomial distribution
# n = trials, k1, k2, ... = number of events, p1, p2, ... = event rates
# multinompdf(n, k's, p's) = n!/(k1! k2! ... km!) * p1^k1, * p2^k2 * ... *pm^km
# logmultinompdf = log(n) - (log(K1!) + ... + log(km!)) + k1*log(p1) + ... + km*log(pm)
logmultinompdf_noC = function(ks, ps){
	return( sum(ks * log(ps)) )
} # logmultinompdf_noC
logmultinompdf_C = function(ks){
	return( lgamma(sum(ks)+1) - sum( lgamma(ks+1) )  )
} # logmultinompdf_C, constant only


## interpolation
interpolate_outcome_agegp = function(outcome_org, agegp_org, agegp_new){
	# agegp_org, agepg_new: end-pt of age groups

	num_agegp_org = length(agegp_org);
	num_agegp_new = length(agegp_new);
	num_outcome_org = length(outcome_org);

	agegp_midpt_org = (head(agegp_org, -1)+tail(agegp_org, -1))/2;
	agegp_midpt_new = (head(agegp_new, -1)+tail(agegp_new, -1))/2;

	if (agegp_org[num_agegp_org]!=85){
		agegp_org_mod = c(agegp_org, 85);
		outcome_org_mod = c(outcome_org, tail(outcome_org,1)/10);
		agegp_midpt_org_mod = (head(agegp_org_mod,-1)+tail(agegp_org_mod,-1))/2;

		outcome_org = outcome_org_mod;
		agegp_org = agegp_org_mod;
		agegp_midpt_org = agegp_midpt_org_mod;
	}

	# the outcome may be negative if not define the bound
	# so set the outcome at last end pt
	# matlab: interp1(agegp_midpt_org, outcome_org, agegp_midpt_new, 'pchip'); 

	output = spline(x=agegp_midpt_org, y=outcome_org, xout=agegp_midpt_new)$y;
	run_pchip_YN = TRUE;
	if (run_pchip_YN & ('pracma' %in% installed.packages())>0){
		output = pracma::pchip(xi=agegp_midpt_org, yi=outcome_org, x=agegp_midpt_new); # pracma::pchip returns y only
	} else{
		# output = approx(x=agegp_midpt_org, y=outcome_org, xout=agegp_midpt_new);
		output = spline(x=agegp_midpt_org, y=outcome_org, xout=agegp_midpt_new)$y; # approx, spline return x and y
	}

	return(output)
} # interpolate_outcome_agegp


## generate the scaling factor for each age group
fun_beta_scaling_map = function(beta_scaling, scaling_age, agegp_count_num=16){
	# map fitting beta_scaling for each agegp
	# length(beta_scaling_map) = agegp_count_num
	num_beta_scaling = length(beta_scaling);
	xout_age = 1:agegp_count_num;
	
	if (missing(scaling_age)){
		if (num_beta_scaling==2){ # 2 discrete beta_scaling
			scaling_age = c(1,7);
			# scaling_age = c(1,agegp_count_num);
		} else if (num_beta_scaling==3){ # 3 discrete beta_scaling
			scaling_age = c(1, 4, 8) # at age 10, 25, 45
		} else if (num_beta_scaling==4){ # 4 beta scaling 
			scaling_age = c(1,7,11,16)
		} else if (num_beta_scaling==8){ # 8 beta scaling 
			scaling_age = union(c(seq(1,agegp_count_num, by=2)[1:7], 16))
		}
	}

	if (num_beta_scaling==2){ # 2 discrete beta_scaling
		# beta_scaling_map = unlist(mapply(function(x,y) rep(x,y), beta_scaling, c(scaling_age, agegp_count_num-scaling_age))) 
		beta_scaling_map = spline(x=scaling_age, y=beta_scaling, xout=xout_age)$y
		if (max(scaling_age)<agegp_count_num){
			beta_scaling_map[ (max(scaling_age)+1):agegp_count_num ] = beta_scaling[2]
		}
	} else if (num_beta_scaling==3){ # 3 discrete beta_scaling
		# beta_scaling_map = spline(x=scaling_age, y=beta_scaling, xout=xout_age)$y
		beta_scaling_map = pracma::pchip(xi=scaling_age, yi=beta_scaling, x=xout_age) # should be the same to spline if there are only 3 points
		if (max(scaling_age)<agegp_count_num){
			beta_scaling_map[ (max(scaling_age)+1):agegp_count_num ] = beta_scaling[num_beta_scaling]
		}
	} else if (num_beta_scaling==4){ # 4 beta scaling 
		beta_scaling_map = spline(x=scaling_age, y=beta_scaling, xout=xout_age)$y
	
	} else if (num_beta_scaling==8){ # 8 beta scaling 
		beta_scaling_map = spline(x=scaling_age, y=beta_scaling, xout=xout_age)$y

	} else if (num_beta_scaling==16){ # 16 beta scaling 
		beta_scaling_map = beta_scaling

	} else{
		stop( "generate_CxOutcome_rate. check number of beta_scaling")
	}
	
	return(beta_scaling_map)
} # fun_beta_scaling_map


# generate the disease outcome following HPV infection
generate_CxOutcome_rate = function(gamdist_par, beta_scaling, HPVincid_gp, agegp_count, pop_agegp, scaling_age=7, mat_ratioAlive){
	# HPVincid_gp: in the range 0-1

	infectCount = HPVincid_gp * pop_agegp;
	agegp_count_num = length(agegp_count);

#	n_beta_scaling = length(beta_scaling);
	beta_scaling_map = rep(beta_scaling[1], agegp_count_num);

	if (tail(beta_scaling,1) > 1){
		beta_scaling_age = tail(beta_scaling,1)
		beta_scaling = head(beta_scaling,-1)
		if (length(beta_scaling)>2){
			beta_scaling_map = as.numeric(approx(x=c(1,beta_scaling_age,agegp_count_num), y=beta_scaling, xout=1:agegp_count_num)$y)
		} else{
			beta_scaling_map = as.numeric(approx(x=c(1,beta_scaling_age), y=beta_scaling, xout=1:agegp_count_num, rule=2)$y)
		}
	} else{
		beta_scaling_map = fun_beta_scaling_map(beta_scaling)
	}


	outCxCount_new = rep(0, agegp_count_num);
	for (jj in 1:agegp_count_num){
		if (jj<agegp_count_num){
			ii_agegp_start = c(0, agegp_count[jj:(agegp_count_num-1)]);
		} else{
			ii_agegp_start = 0;
		}
		ii_agegp_end = agegp_count[jj:agegp_count_num];
		ii_agegp_start_cumsum = cumsum(ii_agegp_start);
		ii_agegp_end_cumsum = cumsum(ii_agegp_end);

		CxInc_temp = infectCount[jj] * apply(pgamma(cbind(ii_agegp_start_cumsum,ii_agegp_end_cumsum), shape=gamdist_par[1], scale=gamdist_par[2], lower.tail=TRUE),1, diff);;

		CxInc_temp = beta_scaling_map[jj] * CxInc_temp

		# consider background deaths, some may die before cancer develops
		CxInc_temp = CxInc_temp * mat_ratioAlive[jj, jj:agegp_count_num]
		
		outCxCount_new[jj:agegp_count_num] = outCxCount_new[jj:agegp_count_num] + CxInc_temp;
	}

	
	# convert to incidence rate by pop_agegp
	out_CxRate = outCxCount_new / pop_agegp;

	return(out_CxRate);

} # generate_CxOutcome_rate


## generate cancer incidence given HPV prevalence, agegp count, population per age gp
# return a list of CxInc for each type
fun_gen_CxInc_100k_byType = function(x, extInput, gen_pRC_YN=TRUE){
	# process the input parameter (x) to generate type-specific CxInc 

	# xdata
	xdata = extInput[['xdata']]

	HPVincid_data_byType = xdata[["HPVincid_data"]]

	paraest_setting = xdata[["paraest_setting"]]
	n_delay = paraest_setting["n_delay"]
	n_scaling = paraest_setting["n_scaling"]
	n_set = paraest_setting["n_set"]
	n_paramfit_noRRprog = 2*n_delay + n_scaling
	
	# extInput_9vHR
	extInput_9vHR = extInput[['extInput_9vHR']]

	set9vHR_list = extInput_9vHR[["set9vHR_list"]]
	jcol_HPVmain_list_bymain = set9vHR_list[["jcol_HPVmain_list_bymain"]]

	numHPVtype = 4;
	CxInc_model_byType = rep(list(NULL), numHPVtype)
	for (iHPV in 1:numHPVtype){
		x_iHPV = x[1:n_paramfit_noRRprog];
		xdata_iHPV = xdata;
		
		# update x, xdata per iHPV
		if (iHPV==1){
			x_iHPV = x[1:n_paramfit_noRRprog];
		} else if (iHPV%in%c(2,3)){
			x_iHPV = x[n_paramfit_noRRprog + (1:n_paramfit_noRRprog)];
		} else if (iHPV==4){
			x_iHPV = x[2*n_paramfit_noRRprog + (1:n_paramfit_noRRprog)];
		}
		
		xdata_iHPV[["HPVincid_data"]] = HPVincid_data_byType[jcol_HPVmain_list_bymain[[iHPV]]]
		
		CxInc_model_iHPV = fun_gen_CxInc_100k(x_iHPV, xdata_iHPV);
		CxInc_model_byType[[iHPV]] = CxInc_model_iHPV
	}

	
	if (gen_pRC_YN){
		sum_CxInc_model_byType = sapply(CxInc_model_byType, sum)
		pRC_model = c(HPV_16=sum_CxInc_model_byType[1], HPV_18_o9vHR=sum(sum_CxInc_model_byType[2:3]), HPV_nonV=sum_CxInc_model_byType[4])/sum(sum_CxInc_model_byType)
		
		CxInc_model_byType = list(CxInc=CxInc_model_byType,
			CxInc_pRC = pRC_model)
	}
	
	return( CxInc_model_byType )

} # fun_gen_CxInc_100k_byType


fun_gen_CxInc_100k = function(x, xdata){
	# agegp_count, number of "age" in each age group, or number of "5-year age" in each 5-year age group
	# e.g., in the age group 85-100, there are three "5-year age"

	# HPVincid_gp: in decimal numbers, not in #, i.e., range 0-1
	HPVincid_gp = xdata[["HPVincid_data"]]; 
	agegp_count = xdata[["outcome_agegp_diff_map"]];
	pop_agegp = xdata[['pop_agegp']];
	paraest_setting = xdata[['paraest_setting']];
	mat_ratioAlive = xdata[['mat_ratioAlive']]

	n_delay = paraest_setting[1];
	n_scaling = paraest_setting[2];
	n_scaling_age = paraest_setting[3];
	scaling_age = paraest_setting[4];
	
	n_gamdist = 2*n_delay;
	gamdist_par = x[1:n_gamdist];
	beta_scaling = x[(n_gamdist+1):length(x)];

	out_CxRate = generate_CxOutcome_rate(gamdist_par, beta_scaling, HPVincid_gp, agegp_count, pop_agegp, scaling_age, mat_ratioAlive);

	out_CxInc_100k = out_CxRate * 100000;

	return(out_CxInc_100k);

} # fun_gen_CxInc_100k


# fun_out_CxInc_100k_LL, returns the logLikelihood directly
# manually check lb and ub
# optim minimize fn, so set negative value for normal loglik and a very large value for x outside the range
# fmin for minimizing, set a very large value for constraint(s), the initial x0 for fmin should be valid w.r.t. the constraints
fun_out_CxInc_100k_LL = function(x, extInput) {
	lb = extInput[['lb']];
	ub = extInput[['ub']];

	wrongValue = list(logL=-1e10, outcome=NULL)
	if (any(x<lb) | any(x>ub)){
		return(wrongValue)
	}
	
	xdata = extInput[['xdata']]
	CxInc_data = extInput[['CxInc_data']];
	idx_CxInc_agegpNew_data = extInput[['idx_CxInc_agegp_data']];
	idx_CxInc_agegpNew_calibrate = extInput[['idx_CxInc_agegp_calib']];
	logL_constant = extInput[['logL_const']]
	extInput_9vHR = extInput[['extInput_9vHR']]
	
	pRC_all9vHR_inonCeCx = extInput_9vHR[["pRC_all9vHR_inonCeCx"]]

	paraest_setting = xdata[["paraest_setting"]]
	n_delay = paraest_setting["n_delay"]
	n_scaling = paraest_setting["n_scaling"]
		n_set = paraest_setting["n_set"]
	n_paramfit_noRRprog = as.numeric(2*n_delay + n_scaling)


	CxInc_model_byType = fun_gen_CxInc_100k_byType(x, extInput, gen_pRC_YN=TRUE);
	
	pRC_model = CxInc_model_byType$CxInc_pRC

	CxInc_model_allAgegp = Reduce('+', CxInc_model_byType$CxInc);
	CxInc_data_LL = CxInc_data; # CxInc for calculating log-likelihood
	if (!is.null(idx_CxInc_agegpNew_calibrate)){
		CxInc_model = CxInc_model_allAgegp[idx_CxInc_agegpNew_calibrate];
		CxInc_data_LL = CxInc_data[match(idx_CxInc_agegpNew_calibrate, idx_CxInc_agegpNew_data)]
	}



	# assume that the first few incidence should be non-decreasing
	CxInc_constraint = c(FALSE
		, any(diff(CxInc_model[1:6])<0)
		);
	if (any(CxInc_constraint)){ return(wrongValue) }
	
	
	# incidence should be non-negative
	CxInc_constraint = c(FALSE
		, any(CxInc_model<0)
		);
	if (any(CxInc_constraint)){ return(wrongValue) }
	
	# pRC_model should be non-negative
	pRC_model_constraint = c(FALSE
		, any(pRC_model<0)
		);
	if (any(pRC_model_constraint)){ return(wrongValue) }	
		


	x_set = sapply(1:n_set, simplify=FALSE, function(ii) x[n_paramfit_noRRprog*(ii-1) + (1:n_paramfit_noRRprog)])# x_set directly input for each set of iHPV

	
	maxratio_beta_scaling = 1000;
	beta_constraint = c(FALSE
		, sapply(x_set, function(xx) xx[5]/xx[3])>maxratio_beta_scaling
		, sapply(x_set, function(xx) xx[5]/xx[3])<(1/maxratio_beta_scaling)
		, sapply(x_set, function(xx) !is_between(xx[4], xx[3], xx[5]))
	)
	if (any(beta_constraint)){ 
		return(wrongValue);
	}
	
	delay_constraint = c(FALSE
		)
	if (any(delay_constraint)){ return(wrongValue) }
	
	
	# logLikelihood, (i) excluding constant and (ii) constant only
	out_temp_CxInc = sum(logpoisspdf_noC(CxInc_data_LL, CxInc_model))
	out_temp_CxInc_pRC = logmultinompdf_noC(ks=pRC_all9vHR_inonCeCx, ps=pRC_model);
	
	out_temp_woConst = out_temp_CxInc + out_temp_CxInc_pRC; # logL without constant
	out_temp_onlyConst = logL_constant + extInput_9vHR$logL_const # logL for constant only

	out_temp = out_temp_woConst + out_temp_onlyConst;# plus the logL_constant

	# using function_MCMC_New.R which requires the function to include both logL and outcome
	out_temp = list(logL=out_temp, outcome=c(CxInc_model_allAgegp, c(pRC_model, HPV_all9vHR=sum(pRC_model[c("HPV_16","HPV_18_o9vHR")]))))
	
	return(out_temp);
}; # fun_out_CxInc_100k_LL


# estimate parameters of the delay function and scaling factor
fun_estParam_HPVincid_to_Cx_byType = function(HPVincid_data, agegp_HPVincid, CxInc_data, agegp_CxInc, pop_age_1yr, pop_agegp, mat_ratioAlive, idx_CxInc_agegpNew_data=NULL, idx_CxInc_agegpNew_calibrate=NULL, runMCMC_YN=FALSE, paramInput=NULL, paraest_setting, extInput_9vHR){
# estimate the delay function and scaling factor given HPV incidence and Cx incidence
# optim, e.g., maximizing the likelihood

	# fval = [gammadist_a, gammadist_b, scaling_factor]
	# HPVincid_data, HPV incidence data; in #, i.e., range 0-100
	# agegp_HPVincid, endpt of HPV incidence data 
	# CxInc_data, cx incidence rate; in 100,000
	# agegp_CxInc, endpt of Cx incidence data; [15:5:85, 85]
	# idx_CxInc_agegpNew, index of CxInc_data in agegp_CxInc, e.g, CxInc_data maybe available for older ages only
	# pop_age_1yr, 1-yr age population
	# pop_agegp, population of age group
	# mat_ratioAlive, relative survival of older age groups when calculating the count of cancer cases, some people may die before the cancer develops
	# paramInput, input of a particular set of parameters, generate the CxInc_model if paramInput is given


	# 1. generate HPVincid with same age grouping as CxInc
	# HPVincid_data_newagegp in the range 0-100
	if (TRUE){ # identical(agegp_HPVincid, agegp_CxInc)){
		HPVincid_data_newagegp = HPVincid_data
	} else{
		HPVincid_data_newagegp = interpolate_outcome_agegp(HPVincid_data, agegp_HPVincid, agegp_CxInc);
	}


	# 2. create xdata for estimating CxInc
	# refer to fun_est_CxInc_100k
	# HPVincid_gp = xdata{1};
	# agegp_count = xdata{2};
	# pop_agegp = xdata{3};
	
	if (missing(pop_agegp)){ # generate pop_agegp from pop_age_1yr
		outcome_agegp = agegp_CxInc;
		outcome_agegp[length(outcome_agegp)] = 100; # max age
		outcome_agegp_num = length(outcome_agegp) - 1;
		outcome_agegp_diff = tail(outcome_agegp,-1) - head(outcome_agegp,-1);
		outcome_agegp_diff_map = outcome_agegp_diff/5;

		pop_agegp = c(sapply(1:(outcome_agegp_num-1), function(i) sum(pop_age_1yr[outcome_agegp[i]:(outcome_agegp[i+1]-1)])), tail(pop_age_1yr,1));
	}


	# 3. estimate the delay function and scaling factor 
	n_delay = paraest_setting[['n_delay']];
	n_scaling = paraest_setting[['n_scaling']];
	n_scaling_age = paraest_setting[['n_scaling_age']];
	scaling_age = unlist(paraest_setting[ grep("^scaling_age", names(paraest_setting), value=TRUE) ]); # paraest_setting[['scaling_age']]; # the age 
	n_RRprog = paraest_setting[['n_RRprog']];
	n_set = paraest_setting[['n_set']];
	n_paramfit = n_set *(2*n_delay + n_scaling + n_scaling_age + n_RRprog);
	
	# lb, ub, and x0 for single set
	lb = c(rep(c(1e-1,1e-1), n_delay), rep(1e-6, n_scaling), rep(1, n_scaling_age), rep(0.01, n_RRprog));
	ub = c(rep(c(20,20), n_delay), rep(0.1, n_scaling), rep(outcome_agegp_num, n_scaling_age), rep(3, n_RRprog));
	x0 = c(rep(c(2.5,2.05),n_delay), seq(0.0005,0.001, length.out=n_scaling), rep(outcome_agegp_num/2, n_scaling_age), rep(1, n_RRprog)); 
	
	# # lb, ub, and x0 for multiple set
	lb = rep(lb, n_set)
	ub = rep(ub, n_set)
	x0 = do.call(c, sapply(seq(1,0.8, length.out=n_set), simplify=FALSE, function(m) m*x0))	


	MCMC_reflective_update_YN = TRUE; # default for function_MCMC
	
	startingPoint_directuse_YN = FALSE;

	# extra input
	if (!is.null(paramInput)){
		if ("x0" %in% names(paramInput)){
			x0 = paramInput[["x0"]]
			print('x0:'); print(x0); print("")
		}
		if ("MCMC_reflective_update_YN" %in% names(paramInput)){
			MCMC_reflective_update_YN = paramInput[["MCMC_reflective_update_YN"]]
			print( sprintf('MCMC_reflective_update_YN: %s', ifelse(MCMC_reflective_update_YN, "Yes", "No")) ); print("")
		}
		if ("startingPoint_directuse_YN" %in% names(paramInput)){
			startingPoint_directuse_YN = paramInput[["startingPoint_directuse_YN"]]
		}
		if ("seed_use" %in% names(paramInput)){
			seed_use = paramInput[["seed_use"]]
		}
	} 
	if (!exists("seed_use")){ 
		seed_use = 1234
	}


	# xdata in fun_gen_CxInc_100k
	# HPVincid_gp: in the range 0-1
	xdata = list('HPVincid_data'=HPVincid_data_newagegp/100, 
		'outcome_agegp_diff_map'=outcome_agegp_diff_map, 
		'pop_agegp'=pop_agegp,
		'paraest_setting'=paraest_setting,
		'mat_ratioAlive'=mat_ratioAlive
		);

	# extInput_9vHR$logL_const = with(as.list(pRC_all9vHR_inonCeCx), sum(logbinompdf_C(n=n, k=c(n_16, n_9vHR))))
	extInput_9vHR$logL_const = with(as.list(pRC_all9vHR_inonCeCx), sum(logmultinompdf_C(ks=c(n_16, n_18_o9vHR, n_nonV))))

	extInput = list('lb'=lb, 'ub'=ub, 
		'xdata'=xdata, 'CxInc_data'=CxInc_data,
		'idx_CxInc_agegp_data'=idx_CxInc_agegpNew_data, 'idx_CxInc_agegp_calib'=idx_CxInc_agegpNew_calibrate,
		'logL_const'=logpoisspdf_C(CxInc_data), # log of the constant in poisson distribution
		'extInput_9vHR'=extInput_9vHR
		);

	# indices for the CxInc in MCMCoutput; match the specific function, e.g., fun_logLikelihood()
	idx_outcome_CxInc = 1:16;
	
	

	if (runMCMC_YN==FALSE){
		# fminsearch if not running MCMC
		fmin_out_CxInc = function(x, extInput) -fun_out_CxInc_100k_LL(x, extInput);
		if (n_scaling>1){
		# suppress Warnings of negative likelihood
		# optim_out = suppressWarnings(optim(par=x0, fn=fun_out_CxInc_100k_LL, lower=lb, upper=ub)); # optim - bounds can only be used with method L-BFGS-B (or Brent)
			out_temp = optim(par=x0, fn=fun_out_CxInc_100k_LL, extInput=extInput, method="Nelder-Mead");
			output_par = out_temp$par
			out_logL_par = -out_temp$value # recall, optim returns min value, so "reverse" for a max value
		} else{
			fmin_method = c("Nelder-Mead", "Hooke-Jeeves")[1];
			out_temp = pracma::fminsearch(x0=x0, fn=fun_out_CxInc_100k_LL, extInput=extInput, method=fmin_method)
			output_par = out_temp$xmin
			out_logL_par = -out_temp$fmin
		}
		
	} else{
		# select the approach for MCMC
		# 1, fmcmc::MCMC; 2, adaptMCMC::MCMC; 3, single-updating or block-updating
		MCMC_method = 3;

		# adaptive MCMC
		# fmcmc::MCMC, kernel-ram, with control of seed
		# adaptmcmc::MCMC, with output of log.p, set.seed before running the line
		
		
		mcmc_nsteps = 20000;
		mcmc_burnin = 1000;
		mcmc_extract = 2000;
		
		rowidx_extract = (mcmc_nsteps+mcmc_burnin)-((mcmc_extract-1):0);
		rowidx_extract = mcmc_nsteps+(1:mcmc_extract);
		
		if (MCMC_method==1 || MCMC_method==2){
			mcmc_param_multiplier = c(rep(1,2*n_delay), rep(1,n_scaling))
			fmcmc_out_CxInc_100k_LL = function(x, extInput) -fun_out_CxInc_100k_LL(x/mcmc_param_multiplier, extInput);
			x0_mcmc = x0*mcmc_param_multiplier
		} else if (MCMC_method==3){
			fmcmc_out_CxInc_100k_LL = function(x, extInput) fun_out_CxInc_100k_LL(x, extInput);
		}
		
		if (MCMC_method==1){ 
			# fmcmc::MCMC
			out_mcmc <- fmcmc::MCMC(
				initial = x0_mcmc, # lb
				fun = fmcmc_out_CxInc_100k_LL,
				nsteps  = 2*mcmc_nsteps, #+ mcmc_burnin,
				seed = seed_use,
				extInput = extInput,
				kernel  = fmcmc::kernel_ram(warmup=500, lb = lb, ub = ub, arate=0.234, eps=0.25)
				)
			
			# output_return = apply(out_mcmc[rowidx_extract,],2, function(x) quantile(x, probs=c(0.5, 0.025,0.975, 0.05,0.95, 0.1,0.9)));
			output_par = apply(out_mcmc[rowidx_extract,],2, median)
			output_par = output_par * mcmc_param_multiplier;
			out_logL_par = max(apply(out_mcmc[rowidx_extract,],1, fmcmc_out_CxInc_100k_LL, extInput=extInput))
			
			# plot( coda::as.mcmc(out_mcmc[rowidx_extract,]) )
			
		# MCMC_method==1
		} else if (MCMC_method==2){
			# adaptMCMC::MCMC
		
			adaptMCMC_scale = c(rep(c(0.8,0.6),n_delay), rep(0.1, n_scaling), rep(1, n_scaling_age))
			# adaptMCMC_scale = c(rep(c(0.05,0.05),n_delay), rep(0.05, n_scaling), rep(1, n_scaling_age))
			set.seed(seed_use) # set.seed before the lines
			# change init to x0_mcmc with multiplers
#			out_mcmc = adaptMCMC::MCMC(fmcmc_out_CxInc_100k_LL, n=mcmc_nsteps, init=x0, extInput=extInput, adapt=TRUE, acc.rate=0.234, scale=adaptMCMC_scale, showProgressBar=FALSE)
			out_mcmc = adaptMCMC::MCMC(fmcmc_out_CxInc_100k_LL, n=2*mcmc_nsteps, init=x0_mcmc, extInput=extInput, adapt=TRUE, acc.rate=0.234, scale=adaptMCMC_scale, showProgressBar=FALSE)

			# output_return = apply(out_mcmc$samples[rowidx_extract,],2, function(x) quantile(x, probs=c(0.5, 0.025,0.975, 0.05,0.95, 0.1,0.9)));
		
			out_logL_par = max(out_mcmc$log.p[rowidx_extract])
			output_par = out_mcmc$samples[which(out_mcmc$log.p==out_logL_par)[1], ]
			output_par = output_par * mcmc_param_multiplier

			out_mcmc_samples_extract = out_mcmc$samples[rowidx_extract, ]
			# out_mcmc_samples_extract = t(apply(out_mcmc_samples_extract, 1, function(x) x*mcmc_param_multiplier))
			out_mcmc_samples_extract = out_mcmc_samples_extract * t(replicate(nrow(out_mcmc_samples_extract), mcmc_param_multiplier))
			
			
			out_mcmc_logL = out_mcmc$log.p[rowidx_extract]
			
			# output_par = out_mcmc$samples[which.max(out_mcmc$log.p), ]
		
			# plot( coda::as.mcmc(out_mcmc$samples[rowidx_extract,]) )
			# plot( adaptMCMC::convert.to.coda(out_mcmc) )
			
		# MCMC_method==2
		} else if (MCMC_method==3){
			# update one-by-one, with adaptive stepSTD
			
			set.seed(seed_use) # set.seed before the lines
			
			stepSTD = c(rep(c(0.4,0.3),n_delay), rep(0.001, n_scaling), rep(0.5, n_scaling_age), rep(0.1, n_RRprog))
			minProbAccept = 0.3;
			maxProbAccept = 0.7;
			stepSTD = rep(stepSTD, n_set);
			
			singleUpdate_YN = FALSE; # single update (TRUE) or block update (FALSE)
			if (singleUpdate_YN){
				out_mcmc = function_MCMC(fun_logLikelihood=fmcmc_out_CxInc_100k_LL, LB=lb, UB=ub, startingPoint=x0, numStepsPerParameter=mcmc_burnin+mcmc_extract, minProbAccept=minProbAccept, maxProbAccept=maxProbAccept, stepSTD=stepSTD, extInput_logL=extInput)
			} else{
				
				block = rep(list(c(list(1:2), as.list(2*n_delay+(1:n_scaling)))), n_set)
				block = do.call(c, sapply(1:n_set, simplify=FALSE, function(x) lapply(block[[x]], function(y) (x-1)*(2*n_delay+n_scaling+n_RRprog) + y )))

				out_mcmc = function_MCMC_block(fun_logLikelihood=fmcmc_out_CxInc_100k_LL, LB=lb, UB=ub, startingPoint=x0, numStepsPerParameter=mcmc_burnin+mcmc_extract, minProbAccept=minProbAccept, maxProbAccept=maxProbAccept, stepSTD=stepSTD, block=block, extInput_logL=extInput, print_updating_YN=TRUE, randomSeed=seed_use, reflective_update_YN=MCMC_reflective_update_YN, startingPoint_directuse_YN=startingPoint_directuse_YN)
			}
						
			rowidx_extract = n_paramfit*(1:(mcmc_burnin+mcmc_extract)) # seq(from=n_paramfit, by=n_paramfit, length.out=mcmc_extract)
			rowidx_extract = rowidx_extract[mcmc_burnin + (1:mcmc_extract)];
			out_mcmc_samples_extract = out_mcmc$chainRecord[rowidx_extract, ]
			out_mcmc_logL = out_mcmc$chainLikelihood[rowidx_extract]
			out_mcmc_chainOutcome = out_mcmc$chainOutcome[rowidx_extract,]

			out_logL_par = max(out_mcmc_logL)
			which_max_logL = which.max(out_mcmc_logL)
			output_par = out_mcmc_samples_extract[which_max_logL,  ]
			output_par_CxInc = out_mcmc_chainOutcome[which_max_logL, idx_outcome_CxInc]
			output_par_CxInc_pRC = out_mcmc_chainOutcome[which_max_logL, -(idx_outcome_CxInc)]
			
		# MCMC_method==3
		}
		
		# browser() # plot MCMC
	} # end runMCMC_YN

	beta_scaling_map_outpar = fun_beta_scaling_map(output_par[2*n_delay + (1:n_scaling)]);
	
	out_logL_par = c(logL=out_logL_par, AIC=2*n_paramfit - 2*out_logL_par, BIC=log(length(CxInc_data))*n_paramfit - 2*out_logL_par)

	return_out = list(estParam=output_par, CxInc_estParam=output_par_CxInc, HPVincid_tofit_estParam=HPVincid_data_newagegp, logL_estParam=out_logL_par, scaling_map_estParam=beta_scaling_map_outpar, CxInc_pRC_estParam=output_par_CxInc_pRC) # CxInc_estParam=Reduce('+', fun_gen_CxInc_temp$CxInc), CxInc_pRC_estParam=fun_gen_CxInc_temp$CxInc_pRC

	if (runMCMC_YN){

		MCMCout_param = out_mcmc_samples_extract
		if (FALSE){
			MCMCout_CxInc_temp = apply(MCMCout_param,1, function(x) fun_gen_CxInc_100k_byType(x, extInput, gen_pRC_YN=TRUE))
			MCMCout_CxInc = do.call(rbind, lapply(MCMCout_CxInc_temp, function(x) Reduce("+", x$CxInc)))
			MCMCout_CxInc_pRC = do.call(rbind, lapply(MCMCout_CxInc_temp, function(x) x$CxInc_pRC))
		} else{
			# with out_mcmc$chainOutcome
			MCMCout_CxInc = out_mcmc_chainOutcome[, idx_outcome_CxInc]
			MCMCout_CxInc_pRC = out_mcmc_chainOutcome[, -(idx_outcome_CxInc)]
		}	
		
		
		MCMCout_logL = cbind(logL=out_mcmc_logL, AIC=2*n_paramfit - 2*out_mcmc_logL, BIC=log(length(CxInc_data))*n_paramfit - 2*out_mcmc_logL)

		return_out = c(return_out, MCMCoutput_estparam=list(MCMCout_param), MCMCoutput_CxInc=list(MCMCout_CxInc), MCMCoutput_logL=list(MCMCout_logL), MCMCoutput_CxInc_pRC=list(MCMCout_CxInc_pRC))
	}
	
	
	return(return_out)

} # fun_estParam_HPVincid_to_Cx


