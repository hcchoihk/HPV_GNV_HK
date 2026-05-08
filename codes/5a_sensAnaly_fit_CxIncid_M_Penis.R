## fit HPV-related cancer incidence, referring to HPV *incidence*
## refer to incidence by main HPV type
# type-specific progression


### file names of outputs
output_suffix = outputdate_string;
output_suffix = paste0(output_suffix, "_sensAnaly", "_fit_CxInc_M_Penis"); # <- change
# "F"/"M"; "Vagina_Vulva (VV)", "OPC", "Anus", "Cervix", "Penis"


fname_estoutput = sprintf("%sestoutput_%s.xlsx", output_folder_fit_CxInc, output_suffix);
fname_estoutput_mcmc = gsub("estoutput", "mcmc_estoutput", fname_estoutput)
RDataname_estoutput = gsub(".xlsx", ".RData", fname_estoutput)

# the base for fname of plots
fname_pdfplot_base = sprintf("%s%%s_%s.pdf", output_folder_fit_CxInc, output_suffix);


print( sprintf("output_suffix = %s", output_suffix) )
writeLines( "Check output_suffix.\nEnter 'c' to confirm and continue. Or, enter 'Q' to quit and modify." )
if (!run_parallel_YN){ browser() }


###

## settings

# HPV incidence setting
HPVincid_agegp_org = c(seq(10, 85, by=5), 85); # age endpt
HPVincid_agegp_new = c(seq(10, 85, by=5), 85); # grouping in CxInc_agegp
HPVincid_min = 0.001; 


# outcome (CxInc) agegp, map from individual ages to 5-year age groups
# take account of the age in HPVincid
outcome_agegp_start = seq(10, 85, by=5); # start age of each gp
outcome_agegp_num = length(outcome_agegp_start);
outcome_agegp_max = 100; # max age
outcome_agegp = c(outcome_agegp_start, outcome_agegp_max); 
outcome_agegp_diff = tail(outcome_agegp,-1) - head(outcome_agegp,-1);
outcome_agegp_diff_map = outcome_agegp_diff/5;


# match the age indices in data_CxIncRate_input and outcome_agegp
# e.g., outcome_agegp starts at age 15, cancer data starts at age 25
idx_CxInc_outcomeAgegp_data = 4:16; # CxInc_data referring to the outcome_agegp
idx_CxInc_outcomeAgegp_calibrate_basic = 4:16; # idx used for calibration, referring to the outcome_agegp
# idx_CxInc_outcomeAgegp_calibrate_basic = 4:15; # exclude the last age group (85+)


### 
fit_CxInc_youngerages_YN = FALSE;

# age group indices for HPV incidence
idx_HPVincid_outcomeAgegp = 1:16;

vec_HPV_set = 1:100;
iset_tempsave = 10;
	
# setting for diseases

CxInc_type_sex_full = list(Female = c("Vagina_Vulva", "OPC", "Anus", "Cervix"), 
	Male = c("OPC", "Penis", "Anus")
	);
sex_vec_full = c("Female","Male");


CxInc_type_sex = with(CxInc_type_sex_full, list(Female = Female, 
	Male = Male
	));
sex_vec = sex_vec_full;

	# the cancer type to fit, if fitting one cancer only in each script
	output_suffix_split = strsplit(output_suffix, "_")[[1]]
	if (any(output_suffix_split=="F")){
		sex_vec = "Female"
		if (any(output_suffix_split=="VV")){
			CxInc_type_sex = list(c("Vagina_Vulva"))
		} else if (any(output_suffix_split=="OPC")){
			CxInc_type_sex = list(c("OPC"))
		} else if (any(output_suffix_split=="Anus")){
			CxInc_type_sex = list(c("Anus"))
		} else if (any(output_suffix_split=="Cervix")){
			CxInc_type_sex = list(c("Cervix"))
		}
		names(CxInc_type_sex) = "Female";
	} else if (any(output_suffix_split=="M")){
		sex_vec = "Male"
		if (any(output_suffix_split=="Penis")){
			CxInc_type_sex = list(c("Penis"))
		} else if (any(output_suffix_split=="OPC")){
			CxInc_type_sex = list(c("OPC"))
		} else if (any(output_suffix_split=="Anus")){
			CxInc_type_sex = list(c("Anus"))
		}
		names(CxInc_type_sex) = "Male";
	}
	cat( sprintf("CxInc_type: %s, sex: %s\n", unlist(CxInc_type_sex), sex_vec) )
if (!run_parallel_YN){ browser() }

sex_vec_short = substr(sex_vec,1,1)
sex_vec_index = seq(sex_vec)

fnameout_suffix = sapply(sex_vec, function(x) paste(x, CxInc_type_sex[[x]], sep="_"))
fnameout_suffix = sapply(fnameout_suffix, function(x) gsub("Vagina_Vulva", "Vagina", x))
fnameout_suffix = sapply(fnameout_suffix, function(x) gsub("Female", "F", gsub("Male", "M", x)))

ncomb_sex_CxInc = length(fnameout_suffix); # number of comb of sex and CxInc



# ** _assumption of the relative risk__ ** for type-specific risk
# for the 4 HPV types/classes, suppose the others (HPV-18/o9vHR/nonVHR) are compared to HPV-16


## parameter estimation
n_delay = 1; # either 1 for age-independent, or same as n_scaling for separate delay functions depending on age
scaling_age = c(1, 4, 8); # knots of the age to change beta; include 1 for age group 10-14 
n_scaling = 3;
n_scaling_age = 0;
n_RRprog = 0; # relative ratio of progression for HPV-non16 vs HPV-16
n_set = 3;
paraest_setting = c(n_delay=n_delay, n_scaling=n_scaling, n_scaling_age=n_scaling_age, scaling_age=scaling_age, n_RRprog=n_RRprog, n_set=n_set);
x0 = c(rep(c(2.5, 2.5), n_delay), rep(0.005, n_scaling), rep(median(idx_HPVincid_outcomeAgegp), n_scaling_age), rep(1, n_RRprog));
x0 = rep(x0, n_set)

x0_input_otherset_YN = TRUE;
x0_directuse_YN = TRUE

runMCMC_YN = TRUE;
MCMC_reflective_update_YN = 1 # if not reflective, just let the proposal out of range


CxInc_datafit_list = rep(list(NULL), ncomb_sex_CxInc);
estparam_list = rep(list(NULL), ncomb_sex_CxInc);
CxInc_estparam_list = rep(list(NULL), ncomb_sex_CxInc);
HPVincid_tofit_list = rep(list(NULL), ncomb_sex_CxInc);
logL_estParam_list = rep(list(NULL), ncomb_sex_CxInc);
scaling_map_list = rep(list(NULL), ncomb_sex_CxInc);
CxInc_pRC_estparam_list = rep(list(NULL), ncomb_sex_CxInc);
MCMCout_estparam_list = rep(list(NULL), ncomb_sex_CxInc);
MCMCout_CxInc_list = rep(list(NULL), ncomb_sex_CxInc);
MCMCout_logL_list = rep(list(NULL), ncomb_sex_CxInc);
MCMCout_scaling_map_list = rep(list(NULL), ncomb_sex_CxInc);
MCMCout_CxInc_pRC_list = rep(list(NULL), ncomb_sex_CxInc);

	# assign names first
	names(estparam_list) = paste0('estpar_', fnameout_suffix);
	names(CxInc_estparam_list) = paste0('estCxInc_', fnameout_suffix);
	names(HPVincid_tofit_list) = paste0('estHPVincid_', fnameout_suffix);
	names(logL_estParam_list) = paste0('logL_estpar_', fnameout_suffix);
	names(scaling_map_list) = paste0('scaling_map_', fnameout_suffix);
	names(CxInc_pRC_estparam_list) = paste0('estCxInc_pRC_', fnameout_suffix);
	names(MCMCout_estparam_list) = paste0('MCMC_par_', fnameout_suffix);
	names(MCMCout_CxInc_list) = paste0('MCMC_CxInc_', fnameout_suffix);
	names(MCMCout_logL_list) = paste0('MCMC_logL_', fnameout_suffix);
	names(MCMCout_CxInc_pRC_list) = paste0('MCMC_CxInc_pRC_', fnameout_suffix);

iparamest_list = 0;
num_vec_HPV_set = length(vec_HPV_set);


for (isex in sex_vec_index){
	CxInc_type_sex_temp = CxInc_type_sex[[sex_vec[isex]]];
	num_CxInc_type_sex_temp = length(CxInc_type_sex_temp);

	if (sex_vec_short[isex]=="F"){
		HPVincid_temp = HPVincid_temp_list[["Female"]]
		pop_age_1yr_temp = pop_age_1yr_Female;
		pop_agegp_temp = pop_agegp_lifetable_Female;
		mat_ratioAlive_temp = mat_ratioAlive_Female;
	} else if (sex_vec_short[isex]=="M"){
		HPVincid_temp = HPVincid_temp_list[["Male"]]
		pop_age_1yr_temp = pop_age_1yr_Male;
		pop_agegp_temp = pop_agegp_lifetable_Male;
		mat_ratioAlive_temp = mat_ratioAlive_Male;
	}
	HPVincid_temp[HPVincid_temp==0] = HPVincid_min;
	
	if (FALSE){
		# same size of all ages
		pop_age_1yr_temp = sum(pop_age_1yr_temp)/100 * c(rep(1,84),16)
	}
	
	
	for (itype in 1:num_CxInc_type_sex_temp){
		iparamest_list = iparamest_list + 1;

		CxInc_name = CxInc_type_sex_temp[itype];
		CxInc_name_out = paste(sex_vec[isex], CxInc_name, sep='_');
		CxInc_name_short = paste(sex_vec_short[isex], CxInc_name, sep='_');
		
		CxInc_datafit = with(data_CxIncRate_input, as.numeric(data_CxIncRate_input[Cancer==CxInc_name & Gender==sex_vec_short[isex], idx_CxInc_input_mean]))

		CxInc_datafit_lower = with(data_CxIncRate_input, as.numeric(data_CxIncRate_input[Cancer==CxInc_name & Gender==sex_vec_short[isex], idx_CxInc_input_lower]))
		CxInc_datafit_upper = with(data_CxIncRate_input, as.numeric(data_CxIncRate_input[Cancer==CxInc_name & Gender==sex_vec_short[isex], idx_CxInc_input_upper]))
		
		CxInc_datafit_list[[iparamest_list]] = rbind(CxInc_datafit, CxInc_datafit_lower, CxInc_datafit_upper);
	
		idx_CxInc_outcomeAgegp_calibrate = idx_CxInc_outcomeAgegp_calibrate_basic;
	
		if (fit_CxInc_youngerages_YN){
			print( "Test CxInc younger ages." )
			CxInc_youngages = cbind(0.01,0.02,0.05);
			CxInc_datafit =  c(CxInc_youngages, CxInc_datafit);
			CxInc_datafit_list[[iparamest_list]] = cbind( pracma::repmat(CxInc_youngages, m=1,n=3), CxInc_datafit_list[[iparamest_list]])
			# update the idx to fit CxInc
			idx_CxInc_outcomeAgegp_data = 1:16; 
			idx_CxInc_outcomeAgegp_calibrate_basic = 1:16; # idx used for calibration, referring to the outcome_agegp
		}
	
		estparam_mat = matrix(0, nrow=num_vec_HPV_set, ncol=length(x0));
		CxInc_estparam_mat = matrix(0, nrow=num_vec_HPV_set, ncol=length(idx_CxInc_outcomeAgegp_data));
		HPVincid_tofit_mat = matrix(0, nrow=num_vec_HPV_set, ncol=with(HPVincid_input, numAgegp*numHPVmain));
		logL_estParam_mat = matrix(0, nrow=num_vec_HPV_set, ncol=3);
		scaling_map_mat = matrix(0, nrow=num_vec_HPV_set, ncol=HPVincid_input$numAgegp); # HPVincid_input$numAgegp*(1+n_RRprog)
		CxInc_pRC_mat_list = rep(list(NULL), num_vec_HPV_set);
		
		MCMCoutput_estparam_iset = rep(list(NULL), num_vec_HPV_set);
		MCMCoutput_CxInc_iset = rep(list(NULL), num_vec_HPV_set);
		MCMCoutput_logL_iset = rep(list(NULL), num_vec_HPV_set);
		MCMCoutput_CxInc_pRC_iset = rep(list(NULL), num_vec_HPV_set);
		

		# ** setting for fitting with 9vHR
		pRC_all9vHR_inonCeCx = df_pRC_9vHR[CxInc_name_short, ] # pRC_all9vHR[CxInc_name_short, ]
		extInput_9vHR = list("pRC_all9vHR_inonCeCx" = pRC_all9vHR_inonCeCx, "set9vHR_list"=set9vHR_list)
		
		
		time1 = Sys.time();
		print( CxInc_name_short )

		for (iset in vec_HPV_set){
			print( sprintf("iset: %d", iset) )

			paramInput = list(MCMC_reflective_update_YN=MCMC_reflective_update_YN)

			# ** set seed for random number generation
			if (exists("seed_use")){
				paramInput = append(paramInput, list(seed_use=seed_use))
			}

			if (x0_input_otherset_YN){
				paramInput = append(paramInput, list(x0=as.numeric(x0_input_otherset_read[[CxInc_name_short]][iset,])
					, startingPoint_directuse_YN = x0_directuse_YN))
			}


			HPVincid_read = as.numeric(HPVincid_temp[iset, ]);
			
			HPVincid_tofit = HPVincid_read; # read by different HPV type 
			# HPVincid_read[idx_HPVincid_outcomeAgegp];
			HPVincid_agegp_tofit = HPVincid_agegp_org[c(idx_HPVincid_outcomeAgegp, length(HPVincid_agegp_org))];
			
			# life table
			fun_output_temp = fun_estParam_HPVincid_to_Cx_byType(HPVincid_tofit, HPVincid_agegp_tofit, CxInc_datafit, HPVincid_agegp_new, pop_agegp=pop_agegp_temp, mat_ratioAlive=mat_ratioAlive_temp, idx_CxInc_agegpNew_calibrate=idx_CxInc_outcomeAgegp_calibrate, idx_CxInc_agegpNew_data=idx_CxInc_outcomeAgegp_data, runMCMC_YN=runMCMC_YN, paraest_setting=paraest_setting, extInput_9vHR=extInput_9vHR, paramInput=paramInput)

			estParam_temp =  fun_output_temp$estParam;
			CxInc_estParam_temp = fun_output_temp$CxInc_estParam
			HPVincid_tofit_temp = fun_output_temp$HPVincid_tofit_estParam
			logL_estParam_temp = fun_output_temp$logL_estParam
			scaling_map_temp = fun_output_temp$scaling_map_estParam
			CxInc_pRC_temp = fun_output_temp$CxInc_pRC_estParam
			
			estparam_mat[iset, ] = estParam_temp;
			CxInc_estparam_mat[iset, ] = CxInc_estParam_temp[idx_CxInc_outcomeAgegp_data];
			HPVincid_tofit_mat[iset, ] = HPVincid_tofit_temp;
			logL_estParam_mat[iset, ] = logL_estParam_temp;
			scaling_map_mat[iset, ] = scaling_map_temp;
			CxInc_pRC_mat_list[[iset]] = CxInc_pRC_temp;
			
			if (runMCMC_YN){
				MCMCoutput_estparam_iset[[iset]] = fun_output_temp$MCMCoutput_estparam;
				MCMCoutput_CxInc_iset[[iset]] = fun_output_temp$MCMCoutput_CxInc;
				MCMCoutput_logL_iset[[iset]] = fun_output_temp$MCMCoutput_logL;
				MCMCoutput_CxInc_pRC_iset[[iset]] = fun_output_temp$MCMCoutput_CxInc_pRC;
			}

			
			if ((iset%%iset_tempsave)==0){
				estparam_list[[iparamest_list]] = estparam_mat;
				CxInc_estparam_list[[iparamest_list]] = CxInc_estparam_mat;
				HPVincid_tofit_list[[iparamest_list]] = HPVincid_tofit_mat;
				logL_estParam_list[[iparamest_list]] = logL_estParam_mat;
				scaling_map_list[[iparamest_list]] = scaling_map_mat;
				CxInc_pRC_estparam_list[[iparamest_list]] = do.call(rbind, CxInc_pRC_mat_list);

		
				if (runMCMC_YN){
					MCMCout_estparam_list[[iparamest_list]] = do.call(rbind, MCMCoutput_estparam_iset)
					MCMCout_CxInc_list[[iparamest_list]] = do.call(rbind, MCMCoutput_CxInc_iset)
					MCMCout_logL_list[[iparamest_list]] = do.call(rbind, MCMCoutput_logL_iset)
					MCMCout_CxInc_pRC_list[[iparamest_list]] = do.call(rbind, MCMCoutput_CxInc_pRC_iset)
				}
				
				save.image(RDataname_estoutput)
				nullout = replicate(5, function(x) {gc(reset=TRUE); Sys.sleep(1)})
			} # iset_tempsave


		} # iset

		estparam_list[[iparamest_list]] = estparam_mat;
		CxInc_estparam_list[[iparamest_list]] = CxInc_estparam_mat;
		HPVincid_tofit_list[[iparamest_list]] = HPVincid_tofit_mat;
		logL_estParam_list[[iparamest_list]] = logL_estParam_mat;
		scaling_map_list[[iparamest_list]] = scaling_map_mat;
		CxInc_pRC_estparam_list[[iparamest_list]] = do.call(rbind, CxInc_pRC_mat_list);
		
		if (runMCMC_YN){
			MCMCout_estparam_list[[iparamest_list]] = do.call(rbind, MCMCoutput_estparam_iset)
			MCMCout_CxInc_list[[iparamest_list]] = do.call(rbind, MCMCoutput_CxInc_iset)
			MCMCout_logL_list[[iparamest_list]] = do.call(rbind, MCMCoutput_logL_iset)
			MCMCout_CxInc_pRC_list[[iparamest_list]] = do.call(rbind, MCMCoutput_CxInc_pRC_iset)
			
		}
		
		
		time2 = Sys.time();
		print(time2 - time1);
		
	} # itype

} # for-isex


## outputs
## plots
# density plot (and trace plot, though not exactly mcmc samples)
pdf_fname = sprintf(fname_pdfplot_base, "plot_estparam");
pdf(pdf_fname, width=3*2, height=3*length(x0))

for (ilist in 1:iparamest_list){
	estParam_mcmc_temp = estparam_list[[ilist]]
	plot(coda::as.mcmc(estParam_mcmc_temp))
	
	# an overall title of the diagnostic plots
	plotx = grconvertX(0.5, from='nfc', to="user")
	ploty = grconvertY(0.995, from='nfc', to="user")
	text(plotx, ploty, labels=names(estparam_list)[ilist], font=2, adj=c(0.5, 1), xpd=TRUE)
}
dev.off() # end traceplot


# plot of CxInc
pdf_fname = sprintf(fname_pdfplot_base, "plot_CxInc_estparam");
pdf(pdf_fname, width=5, height=5)
par(mar=c(3,4,2,1))
for (ilist in 1:iparamest_list){

	plotx_CxInc = HPVincid_agegp_new[idx_CxInc_outcomeAgegp_data];
		
	CxInc_estparam_temp = CxInc_estparam_list[[ilist]];
	CxInc_estparam_temp_median = apply(CxInc_estparam_temp, 2, median);
	
	CxInc_datafit_temp = CxInc_datafit_list[[ilist]];
	ploty_CxInc_datafit = t( CxInc_datafit_temp ); 

	ylim_plot = c(0, ceiling(max(max(CxInc_estparam_temp), max(CxInc_datafit_temp))))
	if (TRUE){
		matplot(x=plotx_CxInc, y=t(CxInc_estparam_temp), col='grey80', type='l', las=1, xlab="Age", ylab="Cancer incidence (per 100,000)", ylim=ylim_plot)
		# particular lines only
		rowidx_toplot = c(which.max(CxInc_estparam_temp[,ncol(CxInc_estparam_temp)]), which.min(CxInc_estparam_temp[,ncol(CxInc_estparam_temp)]))
		matlines(x=plotx_CxInc, y=t(CxInc_estparam_temp[rowidx_toplot,]), col='darkgreen', las=2)
	} else{
		# particular lines only
		rowidx_toplot = c(which.max(CxInc_estparam_temp[,ncol(CxInc_estparam_temp)]), which.min(CxInc_estparam_temp[,ncol(CxInc_estparam_temp)]))
		# matplot(x=plotx_CxInc, y=t(CxInc_estparam_temp[rowidx_toplot,]), col='grey80', type='l', las=1, xlab="Age", ylab="Cancer incidence (per 100,000)", ylim=ylim_plot)
	}
	lines(x=plotx_CxInc, y=CxInc_estparam_temp_median, col='orange', lwd=2);
	
	matlines(x=plotx_CxInc, y=ploty_CxInc_datafit, col='blue', lwd=c(2,1.5,1.5), lty=c(1,2,2))
	
	mtext(text=names(CxInc_estparam_list)[ilist], side=3, line=0.5)
	legend('topleft', legend=c('Empirical (10-year average)', 'Empirical (10-year lowest/highest)', 'Model estimates', 'Median of model estimates'), col=c('blue','blue', 'grey50','orange'), lwd=c(2,1.5,1.5,2), lty=c(1,2,1,1), cex=0.8)
}
dev.off() # end CxInc


# plot of CxInc - lines for CIs
pdf_fname = sprintf(fname_pdfplot_base, "plot_CxInc_estparam_CI");
pdf(pdf_fname, width=5, height=5)
par(mar=c(3,4,2,1))

plotx = seq(10, 85, by=5)
xlim_plot = range(plotx) # c(10, 85)
range_plot = 90;

for (ilist in 1:iparamest_list){

	plotx_CxInc = HPVincid_agegp_new[idx_CxInc_outcomeAgegp_data];
		
	CxInc_estparam_temp = CxInc_estparam_list[[ilist]];
	ploty_CxInc = t( apply(CxInc_estparam_temp,2, quantile, probs=c(0.5, 0.5+range_plot/100/2*c(-1,1))) )

	CxInc_datafit_temp = CxInc_datafit_list[[ilist]];
	ploty_CxInc_datafit = t( CxInc_datafit_temp ); 

	ylim_plot = c(0, ceiling(max(max(CxInc_estparam_temp), max(CxInc_datafit_temp))))

	plot(NA, xlim=xlim_plot, ylim=ylim_plot, axes=FALSE, xlab="", ylab="Cancer incidence (per 100,000)")
	mtext(side=1, text='Age (years)', line=2)
	axis(1, at=plotx)
	axis(2, las=1)
	
	matlines(x=plotx_CxInc, y=ploty_CxInc_datafit, col='blue', lwd=c(2,1.5,1.5), lty=c(1,2,2))
	matlines(x=plotx_CxInc, y=ploty_CxInc, col='orange', lwd=2, lty=c(1,2,2))

	mtext(text=names(CxInc_estparam_list)[ilist], side=3, line=0.5)

	legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot), 'Empirical (10-year average)', 'Empirical (10-year lowest/highest)'), col=c(rep('orange',2),rep('blue',2)), lty=rep(c(1,2),2), lwd=2, cex=0.8)
}
dev.off() # end CxInc


# plot of delay function
pdf_fname = sprintf(fname_pdfplot_base, "plot_delay_estparam");
pdf(pdf_fname, width=5, height=5)
par(mar=c(3,4,2,1))
for (ilist in 1:iparamest_list){
	estParam_mcmc_temp = estparam_list[[ilist]]
	
	plotx_range = c(0, 16)
	plotx = seq(plotx_range[1],plotx_range[2], by=0.1)
	ploty = apply(estParam_mcmc_temp,1, function(x) dgamma(plotx, shape=x[1], scale=x[2]))
	ploty_median = apply(ploty,1, median)
	
	matplot(x=plotx, y=ploty, type='l', xlab='x', ylab='density', col='grey80', las=1, axes=FALSE)
	lines(plotx, ploty_median, col='blue', lwd=2)
	axis(2, las=1)
	xtick_at = seq(plotx_range[1], plotx_range[2], by=1);
	xlab_at = 5*xtick_at;
	xlab_at[seq(2,length(xlab_at),by=2)] = "";
	# axis(1, at=xtick_at, labels=xlab_at, cex.axis=0.9, xpd=TRUE)
	axis(1, at=xtick_at, labels=NA)
	text(x=xtick_at, y=0, labels=xlab_at, adj=c(0.5, 3.75), xpd=TRUE)
	
	mtext(side=1, text="Time (years)", line=1.75)
	mtext(side=3, text=names(estparam_list)[ilist], line=0.5)
	legend('topright', legend=c('Estimates', 'Median of estimates'), col=c('grey50','blue'), lty=1, cex=0.85)
	
}
dev.off() # end plot of delay function


# plot of HPVincid
pdf_fname = sprintf(fname_pdfplot_base, "plot_HPVincid");
pdf(pdf_fname, width=5, height=5)
par(mar=c(3,4,2,1))
for (ilist in 1:iparamest_list){
	HPVincid_temp = HPVincid_tofit_list[[ilist]]
	HPVincid_temp = Reduce("+", lapply(set9vHR_list$jcol_HPVmain_list_bymain, function(x_jcol) HPVincid_temp[, x_jcol])) # transform to overall HPVincid
	
	plotx = head(HPVincid_agegp_new, -1);
	ploty = t( HPVincid_temp );
	ploty_median = apply(ploty,1, median)
	
	matplot(x=plotx, y=ploty, type='l', xlab='Age', ylab='HPV incidence', col='grey80', las=1)
	lines(plotx, ploty_median, col='blue', lwd=2)
	
	mtext(side=3, text=names(HPVincid_tofit_list)[ilist], line=0.5)
	legend('topright', legend=c('Estimates', 'Median of estimates'), col=c('grey50','blue'), lty=1)
}
dev.off() # end plot of HPVincid


# plot of HPVincid plus CxInc
pdf_fname = sprintf(fname_pdfplot_base, "plot_HPVincid_CxInc");
pdf(pdf_fname, width=6, height=8)
layout(rbind(1,2)) # different magnitudes, so not to put together
par(mar=c(3,4,1.5,1))
for (ilist in 1:iparamest_list){
	plotx = seq(10, 85, by=5)
	xlim_plot = range(plotx) # c(10, 85)
	range_plot = 90;
	
	# plot HPV incidence
	HPVincid_temp = HPVincid_tofit_list[[ilist]]
	HPVincid_temp = Reduce("+", lapply(set9vHR_list$jcol_HPVmain_list_bymain, function(x_jcol) HPVincid_temp[, x_jcol])) # transform to overall HPVincid
	
	plotx_HPVincid = head(HPVincid_agegp_new, -1);
	ploty_HPVincid = t( apply(HPVincid_temp,2, quantile, probs=c(0.5, 0.5+range_plot/100/2*c(-1,1))) )

	ylim_plot = c(0, ceiling(max(ploty_HPVincid)))

	matplot(x=plotx_HPVincid, y=ploty_HPVincid, type='l', xlim=xlim_plot, ylim=ylim_plot, xlab="", ylab='HPV incidence', col='darkgreen', las=1, lwd=2, lty=c(1, 2,2), axes=FALSE)
	mtext(side=1, text='Age (years)', line=2)
	axis(1, at=plotx)
	axis(2, las=1)
	legend('topright', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot)), col='darkgreen', lty=c(1,2), lwd=2, cex=0.8)

	# an overall title of the diagnostic plots
	text_overalltitle = gsub('estHPVincid_','', names(HPVincid_tofit_list)[ilist])
	mtext(side=3, text=text_overalltitle, line=0.5, font=2)
# 	plotx_title = grconvertX(0.5, from='nfc', to="user")
# 	ploty_title = grconvertY(0.95, from='nfc', to="user")
# 	text(plotx_title, ploty_title, labels=text_overalltitle, font=2, adj=c(0.5, 1), col='red', xpd=TRUE)


	# plot CxInc
	CxInc_temp = CxInc_estparam_list[[ilist]]
	plotx_CxInc = HPVincid_agegp_new[idx_CxInc_outcomeAgegp_data];
	ploty_CxInc = t( apply(CxInc_temp,2, quantile, probs=c(0.5, 0.5+range_plot/100/2*c(-1,1))) )

	CxInc_datafit_temp = CxInc_datafit_list[[ilist]];
	ploty_CxInc_datafit = t( CxInc_datafit_temp ); 
	# CxInc_datafit_plot = CxInc_datafit_temp[1,] 
	# CxInc_datafit_plot_lower = CxInc_datafit_temp[2,] 
	# CxInc_datafit_plot_upper = CxInc_datafit_temp[3,]
		
	ylim_plot = c(0, ceiling(max(max(CxInc_temp), max(ploty_CxInc_datafit))))
	
	plot(NA, xlim=xlim_plot, ylim=ylim_plot, axes=FALSE, xlab="", ylab="Cancer incidence (per 100,000)")
	mtext(side=1, text='Age (years)', line=2)
	axis(1, at=plotx)
	axis(2, las=1)

	matlines(x=plotx_CxInc, y=ploty_CxInc, col='orange', lwd=2, lty=c(1,2,2))
	matlines(x=plotx_CxInc, y=ploty_CxInc_datafit, col='blue', lwd=2, lty=c(1,2,2))
	legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot), 'Empirical (10-year average)', 'Empirical (10-year lowest/highest)'), col=c(rep('orange',2),rep('blue',2)), lty=rep(c(1,2),2), lwd=2, cex=0.8)
	
}
dev.off() # end plot of HPVincid plus CxInc


# plot of scaling_map
pdf_fname = sprintf(fname_pdfplot_base, "plot_scaling_map");
pdf(pdf_fname, width=5, height=5)
par(mar=c(3,5,2,1))
for (ilist in 1:iparamest_list){
	scaling_map_temp = scaling_map_list[[ilist]][,1:16]
	
	plotx = 1:16;
	plotx_at = seq(1, 16, by=2);
	plotx_lab = seq(10, 80, by=10);
	ploty = t( scaling_map_temp );
	ploty_median = apply(ploty,1, median)
	ylim_plot = range(pretty(c(0, max(scaling_map_temp))))
	
	matplot(x=plotx, y=ploty, type='l', xlab='Age', ylab='', col='grey80', las=1, axes=FALSE, ylim=ylim_plot)
	lines(plotx, ploty_median, col='blue', lwd=2)
	axis(1, at=plotx_at, labels=plotx_lab)
	axis(2, las=1)
	mtext(side=2, text='scaling factor', line=4);
	mtext(side=3, text=names(scaling_map_list)[ilist], line=0.5)
	legend('topright', legend=c('Estimates', 'Median of estimates'), col=c('grey50','blue'), lty=1)
}
dev.off() # end plot of scaling_map


## export estimates, modify the output after plotting
fun_estparam_stat = function(mat){
	mat = cbind(mat, blank=NA, "meandur"=mat[,1]*mat[,2], "p_below2"=pgamma(2, shape=mat[,1], scale=mat[,2], lower.tail=TRUE), "p_abv10"=1-pgamma(10, shape=mat[,1], scale=mat[,2], lower.tail=TRUE), "p_abv11"=1-pgamma(11, shape=mat[,1], scale=mat[,2], lower.tail=TRUE), "ratio"=mat[,ncol(mat)]/mat[,3])
	nstat_estparam = 5;
	if (FALSE){
		mat = cbind(mat, blank=NA, "meandur"=mat[,1]*mat[,2], "p_below2"=pgamma(2, shape=mat[,1], scale=mat[,2], lower.tail=TRUE), "p_abv10"=1-pgamma(10, shape=mat[,1], scale=mat[,2], lower.tail=TRUE), "p_abv11"=1-pgamma(11, shape=mat[,1], scale=mat[,2], lower.tail=TRUE))
		nstat_estparam = 4;
		if (n_scaling==2){
			mat = cbind(mat, "ratio"=mat[,4]/mat[,3]) 
			nstat_estparam = 5;
		}
	}
	mat_stat = apply(mat[, ncol(mat)-(nstat_estparam:1)+1],2, function(x) c(median(x), min(x), max(x), quantile(x, probs=c(0.05, 0.95))))
	colnames(mat_stat) = paste0('stat_', colnames(mat_stat))
	mat_stat = rbind(mat_stat, do.call(rbind, sapply(6:nrow(mat), simplify=FALSE, function(x) rep(NA,nstat_estparam))))
	mat = cbind(mat, mat_stat)
	return(mat)
} # fun_estparam_stat


if (length(vec_HPV_set)>=10){ 
	estparam_list = lapply(estparam_list, fun_estparam_stat)
}
openxlsx::write.xlsx(c(estparam_list, CxInc_estparam_list, HPVincid_tofit_list, logL_estParam_list, scaling_map_list, CxInc_pRC_estparam_list), file=fname_estoutput, overwrite=TRUE)

save.image(RDataname_estoutput)
