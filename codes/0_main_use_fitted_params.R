### Evaluate the cost-effectiveness of gender-neutral vaccination (GNV) of HPV vaccination

### Flow 
# 1. HPV infection and cervical cancer
# 1.1. Obtain parameter sets via calibration to empirical data of cervical cancer incidence, HPV prevalence, and transition (progression/regression) of pre-cancerour health states.
# 1.2. Use the transmission dynamic model to obtain the time-varying force-of-infection (FOI) and HPV incidence across the time horizon following HPV infection
# 1.3. Obtain the changes in cervical cancer incidence following HPV vaccination and cervical screening across time horizon 
# Note1. Outcomes of HPV vaccination on HPV infection and cervical cancer have been studied previously. The respective outcomes are saved in '.RData' file(s).
# 2. HPV-related non-cervical cancer
# 2.1. Use a three-component function to link genital HPV incidence (without HPV vaccination) to incidence non-cervical cancer. 
# Infer parameters for each cancer and for each set of genital HPV incidence generated in (1), via calibration to empirical cancer incidence.
# 2.2. Estimate the changes in incidence based on the changes in HPV incidence after HPV vaccination (which was generated in (1)).
# 3. Cost-effectiveness analysis of comparing different vaccination strategies
# 2dFOV, two-dose female-only vaccination; 2F2M, GNV of two-dose schedule for both girls and boys; 1F1M, GNV of one-dose schedule for both girls and boys; 2F1M, GNV of two-dose schedule for girls and one-dose schedule for boys


## load packages
check_installed_packages = installed.packages()
check_package_vec = c("openxlsx", "psych", "pracma", "ggplot2", "gplots", "vioplot")
check_package_YN = (check_package_vec %in% check_installed_packages)
if (any(!check_package_YN)){
	stop( sprintf("The following package(s) is/are not available:\n%s\nPlease first install the corresponding package(s).", paste(check_package_vec[which(!check_package_YN)], collapse=", ")) )
}

## load library
library(openxlsx)
library(psych) # for correlation plot for MCMC
library(pracma) # pracma::pchip interpolation

# get the date for saving the outputs 
update_outputdate = function(shortDate=TRUE) {
	ifelse(shortDate, substr(gsub("-", "", Sys.Date()),3,8), gsub("-", "", Sys.Date()))
}



## the main folder / directory to hold the codes and data
folder_main = "C:/Users/USER_NAME/Documents/NEW_FOLDER/"     # ** <-- change this


cat( sprintf("'folder_main: %s'", folder_main) )
if (!dir.exists(folder_main)){
	stop( paste0("'folder_main' does not exist. Check the path before running.\n") ) # not to create the folder automatically
}

cat( sprintf("Put the codes and data in 'codes' and 'data', respectively\n") )
cat( sprintf("Outputs will be saved in 'output_use_fitted_params'\n") )

setwd(folder_main)



### 2. Estimate parameters for modeling HPV-related non-cervical cancer.
## 2.1. Run MCMC each parameter set that is calibrated to HPV infection and cervical cancer
# This takes a few day. Suggest to multiple R sessions to run the parameter inference, using one session for one HPV-related cancer.
# An alterative method is to run for-loop in parallel, e.g., future.apply (this may require substantial computer memory)
# better to use an absolute path when running in parallel

runMCMC_nonHPVCx_YN = FALSE;  # <-- Skip inferring parameters for non-cervical cancers. Use the parameter sets that were previously fitted.
run_parallel_YN = TRUE;

seed_use = 1234;

if (run_parallel_YN){
	library(future.apply)
	if (!is(plan(), "multisession")){ # start a multisession if it does not yet exist
		num_Core = parallel::detectCores()
		plan(multisession, workers = min(6, num_Core-2)) 
	}
}

## load functions
fname_fun_fitting = "2_fun_for_fit_CxIncid_validate_byType.R";
source( paste0("codes/", fname_fun_fitting) )

fname_fun_MCMC = "2_function_MCMC_New.R"; # <-- change if needed
source( paste0("codes/", fname_fun_MCMC) )


## load data
# population data
# Lifetable data
# cancer incidence data, referring to HKCaR and other sources
# relative contribution of 9vHPV types
# HPV incidence from the calibrated model

fname_dataload = "HPVmodel_nonCeCx_parInfer_data.RData"
load( paste0("data/", fname_dataload) )


## outputs 
output_folder = "output_use_fitted_params/primary/" # "" if the same directory
outputdate_string = update_outputdate(short=TRUE)

	outputdate_string = "260430"  # <-- when running for using previous fitted parameter sets

if (!dir.exists(output_folder)){ # create folder 
	writeLines( sprintf("create output_folder: %s", output_folder) )
	dir.create(output_folder)
}

## 2.1. Parameter estimation for model calibration to non-cervical cancer incidence
# output folder for fit_CxInx
output_folder_fit_CxInc = gsub("(/){2,}", "/", paste(output_folder, sprintf("%s_fit_CxInc/", outputdate_string), sep="/"));
if (run_parallel_YN){
	output_folder_fit_CxInc = paste0(folder_main, output_folder_fit_CxInc); # an absolute path
}
if (!dir.exists(output_folder_fit_CxInc)){ # create folder 
	writeLines( sprintf("create folder for fit_CxInc: %s", output_folder_fit_CxInc) )
	dir.create(output_folder_fit_CxInc)
}

if (runMCMC_nonHPVCx_YN){
	fit_CxIncid_script_vec = sapply(c("F_VV", "F_OPC", "F_anus", "M_penis", "M_OPC", "M_anus"), function(xCx) sprintf("codes/2a_fit_CxIncid_%s.R", xCx))
	iiCx_vec_fit_CxInc = 1:length(fit_CxIncid_script_vec);
	
	if (!run_parallel_YN){
		# run one-by-one manually
		for (iiCx in iiCx_vec_fit_CxInc){
			source(fit_CxIncid_script_vec[iiCx])
		}
	} else {
		# run in parallel
		# absolute path, load global variables
		# future.seed=FALSE if the script fit_CxInc already includes the line of set.seed()
		nullout <- future_sapply(iiCx_vec_fit_CxInc, simplify=FALSE, future.globals=ls(), future.seed=FALSE, function(iiCx) source( paste0(folder_main, fit_CxIncid_script_vec[iiCx])) ) 
	}
	
	# combine MCMC results that fitted different cancer incidence
	source( "codes/2a_fit_CxIncid_mcmc_plots_bymain.R" )
	
} # if- runMCMC_nonHPVCx_YN


## 2.2. Estimate the changes in incidence of non-cervical cancer following HPV vaccination.
output_folder_proj_CxInc = gsub("(/){2,}", "/", paste(output_folder, sprintf("%s_proj_CxInc/", outputdate_string), sep="/"));
if (run_parallel_YN){
	output_folder_proj_CxInc = paste0(folder_main, output_folder_proj_CxInc); # use an absolute path
}
if (!dir.exists(output_folder_proj_CxInc)){ # create folder 
	writeLines( sprintf("create folder for proj_CxInc: %s", output_folder_proj_CxInc) )
	dir.create(output_folder_proj_CxInc)
}


runProj_HPVvacc_nonHPVCx_YN = TRUE
if (runProj_HPVvacc_nonHPVCx_YN){
	if (FALSE){ # whether to clean the space in the R session
		var_tokeep = c("source_lines","update_outputdate", "run_parallel_YN", "fname_fun_fitting", "folder_main", "output_folder", "outputdate_string", "output_folder_fit_CxInc")
		rm(list=setdiff(ls(), var_tokeep))
		source( paste0("codes/", fname_fun_fitting) )
	}
	
	proj_CxIncid_script_vec = sapply(c("GNV", "GNV_Vyr30", "2F1M", "2F1M_Vyr30", "noMaleCx", "noMaleCx_Vyr30"), function(xCx) sprintf("codes/2b_project_CxIncid_%s.R", xCx))
	iiCx_vec_proj_CxIncid = 1:length(proj_CxIncid_script_vec);
	
	if (!run_parallel_YN){
		# run one-by-one manually
		for (iiCx in iiCx_vec_proj_CxIncid){
			source(proj_CxIncid_script_vec[iiCx])
		}
	} else {
		# run in parallel
		# absolute path, load global variables
		nullout <- future_sapply(iiCx_vec_proj_CxIncid, simplify=FALSE, future.globals=ls(), future.seed=FALSE, function(iiCx) source( paste0(folder_main, proj_CxIncid_script_vec[iiCx])) ) 
	}
} # if- runProj_HPVvacc_nonHPVCx_YN




### 3. Cost-effectiveness analysis
## flow in CEA_nonCeCx_xx.R
# outcomes for cervical cancers - load results generated from C++ for the stochastic model that incorporate cervical screening and FOI following vaccination
# outcomes for non-cervical cancers - load the projected cancer incidence after vaccination
# compare cost and QALYs across different vaccination schedules
# print out numerical results in text files (txt)
# create plots in PDF files


if (TRUE){ # whether to clean the space in the R session, to save file for CEA only
	var_tokeep = c("source_lines","update_outputdate", "run_parallel_YN", "folder_main", "output_folder", "outputdate_string", "output_folder_proj_CxInc")
	rm(list=setdiff(ls(), var_tokeep))
}


# Load static data (population, lifetable, health economic parameters, data_HPVattrib, etc.)
# Population projection: %d years x %d ages"
load( "data/CEA_static_data.RData" )


# whether to consider herd protection for genital warts (sensitivity analysis)
# default is TRUE, set at FALSE for sensitivity analysis
herdProt_GWart_YN_vec = c(TRUE, FALSE);

for (herdProt_GWart_YN in herdProt_GWart_YN_vec){


	# output folder for CEA
	output_folder_CEA = gsub("(/){2,}", "/", paste(output_folder, sprintf("%s_CEA/", outputdate_string), sep="/"));

	if (herdProt_GWart_YN==FALSE){ # save results in another folder
		output_folder_CEA = gsub("(/){2,}", "/", paste(output_folder, sprintf("%s_CEA_noherdGW/", outputdate_string), sep="/"));
	}

	if (run_parallel_YN){
		output_folder_CEA = paste0(folder_main, output_folder_CEA); # better to use absolute path when running in parallel
	}
	if (!dir.exists(output_folder_CEA)){ # create folder 
		writeLines( sprintf("create folder for CEA: %s", output_folder_CEA) )
		dir.create(output_folder_CEA)
	}


	## 3.1. calculate CEA for 2F2M, 1F1M, 2F1M
	# 2F2M and 1F1M - GNV with same schedule on both schoolgirls and schoolboys
	# 2F1M - two-dose schedule for schoolgirls and one-dose schedule for schoolboys

	CEA_script_vec = sapply(c("GNV", "GNV_Vyr30", "2F1M", "2F1M_Vyr30"), function(xCx) sprintf("codes/3a_CEA_nonCeCx_%s.R", xCx))
	iiCx_vec_CEA_script_vec = 1:length(CEA_script_vec);

	if (!run_parallel_YN){
		# run one-by-one manually
		for (iiCx in iiCx_vec_CEA_script_vec){
			source(CEA_script_vec[iiCx])
		}
	} else {
		# run in parallel
		# absolute path, load global variables
		nullout <- future_sapply(iiCx_vec_CEA_script_vec, simplify=FALSE, future.globals=ls(), future.seed=FALSE, function(iiCx) source( paste0(folder_main, CEA_script_vec[iiCx])) ) 
	}

	# 3.2. compare 2F1M vs 1F1M
	source( 'codes/3b_comp_2F1M_1dGNV.R') 

} # for- herdProt_GWart_YN  


