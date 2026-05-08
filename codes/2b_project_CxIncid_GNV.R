## project HPV-related cancer incidence, based on the inferred parameters, referring to HPV *incidence*

# HPV incidence _* by cohort *_
# generate the cancer incidence for each _cohort_ (row) of HPV incidence
# one cancer for each file, each file contains nPSA sheets, each sheet contain time_horizon cancer incidence; repeat for each VaccSetting (i.e., one set of file for each VaccSetting)

# assigned PSA for HPV attribution (data_HPVattrib)
# ** Loads pre-saved data from HPVmodel_nonCeCx_project_data_2F2M_1F1M.RData **
# Run Read_data_HPVmodel_nonCeCx_project.R first to generate the RData file.


# use one folder to contain files

library(openxlsx)


# output suffix 
output_suffix = outputdate_string;

output_folder_proj_CxInc_new = gsub("(/){2,}","/", paste(output_folder_proj_CxInc, "GNV", "/", sep="/"));
dir.create(output_folder_proj_CxInc_new, showWarnings = FALSE)


# Load pre-saved data
fname_loaddata = paste0(folder_main, "data/HPVmodel_nonCeCx_project_data_GNV.RData");
load( fname_loaddata )
# Loaded variables include:
#   pop_age_1yr_Female, pop_age_1yr_Male
#   pop_agegp_lifetable_Female, pop_agegp_lifetable_Male
#   mat_ratioAlive_Female, mat_ratioAlive_Male
#   data_HPVattrib, data_VE_input
#   CxInc_type_sex_full, sex_vec, sex_vec_short, sex_vec_index
#   ncolumn_extra_cohortinfo, ncolumn_HPVincid_AgeGp, num_column_HPVincid_file
#   iHPVgp_bymain_list, jcol_HPVmain_seperate, jcol_HPVmain_list,
#   jcol_HPVmain_list_bymain, set9vHR_list, HPVincid_input_settings
#   Vcover_vec, Vcover_M_vec, Vdur_vec, VaccSetting_grid,
#   VaccSetting_base_string_vec, VaccSetting_suffix_string_list, vec_HPV_set
#   date_HPVincid_PSA, folder_PSA, fname_HPVincid_PSA_base,
#   fname_HPVincid_bymain_PSA_base, fname_HPVincid_bymain_PSA_sprintf
#   HPVincid_bymain_PSA  (list: [[sex]][[iVaccSet]][[iset]] -> matrix)
#   irow_ref_noVacc


# match the vaccination setting
Vcover_vec = 85;
Vcover_M_vec = c(0, 85, 50, 25);
Vdur_vec = c(100, 20);

VaccSetting_grid = expand.grid(Vcover_vec, Vcover_M_vec, Vdur_vec);

VaccSetting_base_string_vec = apply(VaccSetting_grid, 1, function(x) sprintf('VC%d%d_age1212_dur%d', x[1], x[2], x[3]));


# load parameters for non-cervical cancer modeling
fname_paramest <- paste0(output_folder_fit_CxInc, sprintf("%s_nonCeCx_PSApara.xlsx", outputdate_string))

CxInc_type_sex_full <- list(
	Female = c("vagina_vulva", "OPC", "anus"),
	Male   = c("OPC", "penis", "anus")
)
CxInc_type_sex = with(CxInc_type_sex_full, list(Female = Female, # change if needed
	Male = Male
	));
sex_vec       <- c("Female", "Male")
sex_vec_short <- substr(sex_vec, 1, 1)
sex_vec_index <- seq_along(sex_vec)

# list of the param for each cancer
param_CxIncfit_CxName <- c(sapply(sex_vec_index, function(i)
	paste(sex_vec_short[i], CxInc_type_sex_full[[i]], sep="_")))

param_CxIncfit_list <- sapply(param_CxIncfit_CxName, simplify=FALSE, function(s)
	as.matrix(read.xlsx(fname_paramest, sheet=paste0(s, "_param"), rowNames=TRUE)))

num_CxIncfit = length(param_CxIncfit_CxName);



# Derived settings (not stored in RData)

# time horizon in the model projection
time_horizon_vec = 1:169;
num_time_horizon = max(time_horizon_vec);

# vector of HPV incidence input to read (PSA)
vec_HPV_set = 1:100; # PSA

# parameter estimation
n_delay = 1;
n_scaling = 2;
n_scaling_age = 0;
scaling_age = 7;
n_RRprog = 0;
n_set = 3;
paraest_setting = c(n_delay=n_delay, n_scaling=n_scaling, n_scaling_age=n_scaling_age, scaling_age=scaling_age, n_RRprog=n_RRprog, n_set=n_set);

n_paramfitted = 2*n_delay + n_scaling + n_scaling_age;

# population data or lifetable data
fitwPopdata_YN = FALSE;

## HPV incidence setting
HPVincid_agegp_org = c(seq(10, 85, by=5), 85); # age endpt
HPVincid_agegp_new = c(seq(10, 85, by=5), 85);
HPVincid_min = 0.001;

# outcome (CxInc) agegp
outcome_agegp_start = seq(10, 85, by=5);
outcome_agegp_num = length(outcome_agegp_start);
outcome_agegp_max = 100;
outcome_agegp = c(outcome_agegp_start, outcome_agegp_max); 
outcome_agegp_diff = tail(outcome_agegp,-1) - head(outcome_agegp,-1);
outcome_agegp_diff_map = outcome_agegp_diff/5;

pop_agegp_Female = c(sapply(1:(outcome_agegp_num-1), function(i) sum(pop_age_1yr_Female[outcome_agegp[i]:(outcome_agegp[i+1]-1)])), tail(pop_age_1yr_Female,1));
pop_agegp_Male = c(sapply(1:(outcome_agegp_num-1), function(i) sum(pop_age_1yr_Male[outcome_agegp[i]:(outcome_agegp[i+1]-1)])), tail(pop_age_1yr_Male,1));

idx_CxInc_outcomeAgegp_data = 4:16;
idx_CxInc_outcomeAgegp_calibrate_basic = 4:16;
idx_CxInc_outcomeAgegp_calibrate_M_Tonsil = idx_CxInc_outcomeAgegp_calibrate_basic;

HPVincid_input = list(header = FALSE, 
	numAgegp = length(HPVincid_agegp_org)-1
	);

idx_HPVincid_outcomeAgegp = 1:16;


# ASR weight
ASRweight_Seig1960_age1085 = c(9, 9, 8, 8, 6, 6, 6, 6, 5, 4, 4, 3, 2, 1, 0.5, 0.5)


## project the cancer incidence
for (isex in sex_vec_index){
	CxInc_type_sex_temp = CxInc_type_sex[[sex_vec[isex]]];
	num_CxInc_type_sex_temp = length(CxInc_type_sex_temp);

	if (sex_vec_short[isex]=="F"){
		pop_age_1yr_temp = pop_age_1yr_Female;
		pop_agegp_temp = pop_agegp_lifetable_Female;
		mat_ratioAlive_temp = mat_ratioAlive_Female;
	} else if (sex_vec_short[isex]=="M"){
		pop_age_1yr_temp = pop_age_1yr_Male;
		pop_agegp_temp = pop_agegp_lifetable_Male;
		mat_ratioAlive_temp = mat_ratioAlive_Male;
	}
	
	
	for (itype in 1:num_CxInc_type_sex_temp){

		CxInc_name = CxInc_type_sex_temp[itype];
		CxInc_name_out = paste(sex_vec[isex], CxInc_name, sep='_');
		CxInc_name_short = paste(sex_vec_short[isex], CxInc_name, sep='_');
		
		paramest_itype_mat = param_CxIncfit_list[[CxInc_name_short]];
		
		HPVattrib_itype = data_HPVattrib[, CxInc_name_short]
		
		print( CxInc_name_short );
		
		# ASR of projCxInc, one file for each CxInc
		fname_out_ASRprojCxInc = sprintf('ASRprojCxInc_%s_%s.xlsx', CxInc_name_short, output_suffix);
		ASRprojCxInc_out = setNames(rep(list(NULL), length(VaccSetting_base_string_vec)), VaccSetting_base_string_vec)
		
		
		for (iVaccSet in seq(VaccSetting_base_string_vec)){
		
		xVaccSet_base = VaccSetting_base_string_vec[iVaccSet];
		xVaccSet = VaccSetting_suffix_string_list[[iVaccSet]];
		fname_out_projCxInc = sprintf('projCxInc_%s_%s_%s.xlsx', xVaccSet_base, CxInc_name_short, output_suffix);
		
		
		CxInc_proj_list_iVaccSet = rep(list(NULL), max(vec_HPV_set));
		names(CxInc_proj_list_iVaccSet)[vec_HPV_set] = paste0("PSA_", vec_HPV_set);
		
		CxInc_pRC_proj_list_iVaccSet = rep(list(NULL), max(vec_HPV_set));
		names(CxInc_pRC_proj_list_iVaccSet)[vec_HPV_set] = paste0("PSA_", vec_HPV_set);
		
		
		time1 = Sys.time();
		for (iset in vec_HPV_set){ # (PSA)
			if (iset%%10==0){ print( sprintf("%s, %s, iset %d", CxInc_name_short, gsub("_Catchup(.*)", "", xVaccSet[iset]), iset) )}
		
			paramest_itype_iset = paramest_itype_mat[iset, ]
			HPVattrib_itype_iset = HPVattrib_itype[iset]
			
			# use pre-loaded HPVincid matrix instead of reading CSV
			HPVincid_read = HPVincid_bymain_PSA[[isex]][[iVaccSet]][[iset]]
			
			CxInc_out_all_itimehor = rep(list(NULL), num_time_horizon);
			CxInc_pRC_out_all_itimehor = rep(list(NULL), num_time_horizon);
			for (itimehor in time_horizon_vec){
			
				HPVincid_itimehor = HPVincid_read[itimehor, ]

				if (fitwPopdata_YN){
					pop_agegp_xdata = pop_age_1yr_temp;
				} else{
					pop_agegp_xdata = pop_agegp_temp;
				}

				xdata = list('HPVincid_data'=HPVincid_itimehor/100, 
					'outcome_agegp_diff_map'=outcome_agegp_diff_map, 
					'pop_agegp'=pop_agegp_xdata,
					'paraest_setting'=paraest_setting,
					'mat_ratioAlive'=mat_ratioAlive_temp
					);
					
				extInput_9vHR = list('set9vHR_list'=set9vHR_list)
					
				extInput = list('xdata' = xdata,
					'extInput_9vHR' = extInput_9vHR)

				CxInc_itime = fun_gen_CxInc_100k_byType(paramest_itype_iset, extInput, gen_pRC_YN=TRUE);

				CxInc_out_all_itimehor[[itimehor]] = Reduce("+", CxInc_itime$CxInc);
				CxInc_pRC_out_all_itimehor[[itimehor]] = c(CxInc_itime$CxInc_pRC, HPV_all9vHR=sum(CxInc_itime$CxInc_pRC[c("HPV_16","HPV_18_o9vHR")]));
			
			} # itimehor

			CxInc_out_all_itimehor = do.call(rbind, CxInc_out_all_itimehor);
			# account for HPV attribution
			CxInc_out_all_itimehor_HPVattrib = HPVattrib_itype_iset*CxInc_out_all_itimehor + (1-HPVattrib_itype_iset)*pracma::repmat(CxInc_out_all_itimehor[irow_ref_noVacc,], n=nrow(CxInc_out_all_itimehor), m=1)
			

			CxInc_proj_list_iVaccSet[[iset]] = CxInc_out_all_itimehor_HPVattrib;
			CxInc_pRC_proj_list_iVaccSet[[iset]] = do.call(rbind, CxInc_pRC_out_all_itimehor);
			
			
		} # iset (PSA)
		
		# ASRprojCxInc
		ASR_CxInc_proj_iVaccSet = sapply(CxInc_proj_list_iVaccSet, function(xCxInc) apply(xCxInc, 1, function(x) sum(x*ASRweight_Seig1960_age1085)/100))
		ASRprojCxInc_out[[iVaccSet]] = ASR_CxInc_proj_iVaccSet
		
		time2 = Sys.time();
		print(time2 - time1);


		write.xlsx(CxInc_proj_list_iVaccSet, file=paste0(output_folder_proj_CxInc_new, fname_out_projCxInc), colNames=FALSE);
		
		fname_out_CxInc_pRC = gsub("projCxInc", "CxInc_pRC", fname_out_projCxInc)
		write.xlsx(CxInc_pRC_proj_list_iVaccSet, file=paste0(output_folder_proj_CxInc_new, fname_out_CxInc_pRC), colNames=FALSE);
		
		
		} # iVaccSet
		
		# write ASRprojCxInc
		write.xlsx(ASRprojCxInc_out, file=paste0(output_folder_proj_CxInc_new, fname_out_ASRprojCxInc), colNames=TRUE);

		
	} # itype

} # for-isex
save.image( paste0(output_folder_proj_CxInc_new, sprintf("export_projCxInc_%s.RData", outputdate_string))  )
