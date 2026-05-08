# Estimate the CE of including HPV-related non-cervical cancers

# update the projected changes in nonCeCx incidence 
# folder_proj_nonCeCx, projCxInc_output_date, fname_proj_nonCeCx_base

# update vaccination scenarios

xVdur_boys_1M = 20; # 20 or 30   # <-- check this
run_2F1M_YN = TRUE;     # <-- check this
# herdProt_GWart_YN;  # controlled outside the script


### load data
# Load PSA summary data for cervical cancer outputs which accounts for cervical screening (from C++)
# structure - data_PSAsumm_list, data_PSAsumm_noS_list - [[iVdur]][[iVC]][[iPSA]]
suffix_base_loaddata = "2F1M"    ## <- change, e.g., GNV, 2F1M, Vyr30
load( gsub("(/){2,}", "/", paste(folder_main, sprintf("data/CEA_CeCxScr_PSAsumm_data_%s.RData", suffix_base_loaddata), sep="/")) )

load( gsub("(/){2,}", "/", paste(folder_main, sprintf("data/CEA_data_HPVincid_gwarts_%s.RData", suffix_base_loaddata), sep="/")) )


# Load projected incidence of non-cervical cancers following HPV vaccination
folder_proj_nonCeCx_base = suffix_base_loaddata
# folder_proj_nonCeCx_base = "XXX"     # <-- change if needed
projCxInc_output_date = outputdate_string
folder_proj_nonCeCx = gsub("(/){2,}","/", paste(output_folder_proj_CxInc, folder_proj_nonCeCx_base, "/", sep="/"));
fname_proj_nonCeCx_base = sprintf("projCxInc_VC85%%d_age1212_dur%%d_%%s_%s.xlsx", projCxInc_output_date); # VC_boy (0, 25, 50, 85), Vdur (20, 100), nonCeCx (6 types); each xlsx file contains sheets for each PSA


# folder for noMaleCx_FOV, no vaccine uptake for female - account for MSM
folder_proj_noMaleCxFOV_base = "noMaleCx"    ## <- change if needed
projCxInc_noMaleCxFOV_output_date = outputdate_string
folder_proj_noMaleCxFOV = gsub("(/){2,}","/", paste(output_folder_proj_CxInc, folder_proj_noMaleCxFOV_base, "/", sep="/"));
fname_proj_noMaleCxFOV_base = sprintf("projCxInc_VC0%%d_age1212_dur%%d_%%s_%s.xlsx", projCxInc_noMaleCxFOV_output_date); # VC_boy (0, 25, 50, 85), Vdur (20, 100), nonCeCx (6 types); each xlsx file contains sheets for each PSA


## output settings
outputdate_use = outputdate_string # may change to another date for not to overwrite previous output
output_folder_suffix = suffix_base_loaddata
# output_folder_suffix = "XXX"     # <-- change if needed
output_folder_CEA_new = gsub("(/){2,}","/", paste(output_folder_CEA, output_folder_suffix, "/", sep="/"));
dir.create(output_folder_CEA_new, showWarnings = TRUE)


## settings
# whether to consider herd protection for genital warts (sensitivity analysis)
# default is TRUE, set at FALSE for sensitivity analysis
if (herdProt_GWart_YN==FALSE){
	print( "No herd effect on genital warts" )
	print( "confirm before running" )
	# browser()

	# if there is no herd effect, vaccinated people may also be protected after vaccine efficacy wanes for a short period, but will have no protection after that period
	# assume partial protection 10 more years = two 5-year age group >> lower protection during partial protection, but no protection after partial protection
	herdProt_set = list(Vdur_20 = c(jcol_endVProtect = 5, jcol_partProtect = 7, relative_partProtect = 0.82), 
		Vdur_30 = c(jcol_endVProtect = 7, jcol_partProtect = 9, relative_partProtect = 0.92)
	)
}

# whether to run the 2F1M strategy
if (run_2F1M_YN==TRUE){
	print( sprintf( "Vdur for boys = %d", xVdur_boys_1M) )
	print(  "run 2F1M. confirm before running" )
	# browser()
}


## vaccination scenario settings
Vage_girls = 12
VC_girls = 85
VC_boys_vec = c(0, 25, 50, 85)# include VC=0 for comparison
Vdur_vec = 100 # set at Vdur=100 for 2F1M

VC_boys_num = length(VC_boys_vec);
Vdur_num = length(Vdur_vec);


library(openxlsx) # write.xlsx
repmat = pracma::repmat


## settings
# MSM subgroup
prop_Mpenis_MSM = 0.046
prop_MOPC_MSM = 0.18
prop_Manus_MSM = 0.4
prop_Mpenis_FOV = 1 - prop_Mpenis_MSM
prop_MOPC_FOV = 1 - prop_MOPC_MSM
prop_Manus_FOV = 1 - prop_Manus_MSM

adj_otherCx_MSM_YN = TRUE # adjustment if excluding the MSM subgroup
Mcancer_adj_vec = c("M_penis", "M_OPC", "M_anus")

# function to adjust for attribution from MSM subgroup
fun_CxInc_noAttrib = function(CxInc_wAttrib, CxInc_noVacc, prop_Attrib){
	# get back the CxInc before accounting for prop_Attrib
	# CxInc_wAttrib, the CxInc after accounting for prop_Attrib
	return( (CxInc_wAttrib - (1-prop_Attrib)*CxInc_noVacc)/prop_Attrib)
} # fun_CxInc_noAttrib

fun_CxInc_adjpropMSM_Vacc_MSM = function(CxInc_noVacc, CxInc_Vacc, CxInc_MSM, prop_Attrib, prop_MSM){
	# estimate the CxInc with Vacc, accounting for prop_MSM; input CxInc_MSM
	return(  (1-prop_Attrib)*CxInc_noVacc + prop_Attrib*((1-prop_MSM)*CxInc_Vacc + prop_MSM*CxInc_MSM)  )
} # fun_CxInc_adjpropMSM_Vacc_MSM



## setting for calculation
Yr_GNVstart = 2026;
Yr_eval_num = 100; # year for time horizon

MinAge_cohort = 10
MaxAge_cohort = 79
Yr_eval_num_bycohort =  (MaxAge_cohort - MinAge_cohort + 1) + Yr_eval_num - 1; # +1 for the number of cohort in Yr1, -1 for excluding the cohort of age 10 in Yr1


# screening uptake
pScrRate_base_beforeCSP = 0.50;
pScrRate_base_afterCSP = 0.50;
pScrRate_base_progChange = pScrRate_base_afterCSP;
pScrRate_age_beforeCSP = 70; # based on Yr_current, Yr_CSP, Yr_ScreenHK
pScrRate_age_afterCSP = 45;
pScrRate_age_progChange = 25;



# get population size for the results on cervical cancer with cervical screening, from the stochastic part (Cpp)
# Use the first element of data_PSAsumm_list to get structure
iVdur_init = 1
iVC_init = 1
iPSA_init = 1
data_PSAsumm_init = data_PSAsumm_list[[iVdur_init]][[iVC_init]][[iPSA_init]]
irow_pick = which(data_PSAsumm_init$Vacc_EvaYear<=Yr_eval_num)

pop_proj_CeCx = matrix(0, length(irow_pick), 1)
pop_proj_VaccBoy = matrix(0, length(irow_pick), 1)
for (irow in irow_pick){
	yr_temp = Yr_GNVstart + data_PSAsumm_init$Vacc_EvaYear[irow] - 1
	age_temp = data_PSAsumm_init$CohortEvaAge[irow]
	if (yr_temp<=max(data_pop_proj_F$Year)){
		irow_datapop = which(data_pop_proj_F$Year==yr_temp)
	} else {
		irow_datapop = nrow(data_pop_proj_F)
	}
	# cohort size among females
	pop_proj_CeCx[irow] = data_pop_proj_F[irow_datapop, age_temp]
	# cohort size among males, for cost of vaccinating boys
	pop_proj_VaccBoy[irow] = data_pop_proj_M[irow_datapop, age_temp]
}

# SubjAtRisk for CeCx output generated from Cpp
SubjAtRisk_vec_CeCx = rep(1, nrow(SubjAtRisk_cohort_list[["Female"]]))


## output of cancer incidence for non-cervical cancers
# 16 columns, each for 5-year age group, from ages 10-14, 15-19, to 80-84, 85+

# potential PSA: attributable fraction, treatment cost, prop_vagina_vulva, 5yrQUtil
# the PSA for 5yrQUtil refer to those used for cervical cancer 

shname_base = "PSA_%d"; # per iPSA
nonCeCx_vec = c("F_vagina_vulva", "F_OPC", "F_anus", "M_OPC", "M_penis", "M_anus");
# nonCeCx_vec = nonCeCx_vec[1]; # test
nonCeCx_num = length(nonCeCx_vec)

# index of  _Mcancer_adj_vec_  in  _nonCeCx_cex_
inonCeCx_adjMSM_vec = match(Mcancer_adj_vec, nonCeCx_vec)


# cost; account for cancer treatment cost, costs for specialist consultation, diagnosis and biopsy
# to be used as fixed value if PSA_costQALY_YN==FALSE
nonCeCx_cost_raw = c("vagina_vulva"=152846, "OPC"=278088, "anus"=246594, "penis"=193780)
prop_vagina_vulva = 33/(33+51); # proportion of vagina among (vagina+vulva);
# relative survival, refer to "SEER - relative survival.xlsx"; account for stage distribution, QALY weight for during (6 months) and post-treatment (4.5 years)
nonCeCx_5yrRS = c("F_vagina_vulva"=0, "F_vagina"=50.99, "F_vulva"=70.61, "F_OPC"=56.77, "F_anus"=70, "M_OPC"=56.90, "M_penis"=66.2, "M_anus"=70.98)/100;
nonCeCx_QUtil_5yr = c("F_vagina_vulva"=0, "F_vagina"=2.31, "F_vulva"=3.28, "F_OPC"=2.55, "F_anus"=3.2, "M_OPC"=2.55, "M_penis"=3.08, "M_anus"=3.26);
nonCeCx_5yrRS["F_vagina_vulva"] = prop_vagina_vulva*nonCeCx_5yrRS["F_vagina"] + (1-prop_vagina_vulva)*nonCeCx_5yrRS["F_vulva"]
nonCeCx_QUtil_5yr["F_vagina_vulva"] = prop_vagina_vulva*nonCeCx_QUtil_5yr["F_vagina"] + (1-prop_vagina_vulva)*nonCeCx_QUtil_5yr["F_vulva"]
nonCeCx_QLost_5yr = 5-nonCeCx_QUtil_5yr; 

# function to calculate LY lost due to nonCeCx death
fun_LYlost_nonCeCx = function(disc_factor, deathAge, maxAge=85){
	if (disc_factor>1){
		stop( sprintf('disc_factor = %.3f >1; disc_rate = %.3f. check disc_factor.', disc_factor, 1/disc_factor-1) );
	} else if (disc_factor==1){
		return( maxAge - deathAge )
	} else{
		return( (1-disc_factor^(maxAge-deathAge))/(1-disc_factor) )
	}
} # fun_LYlost_nonCeCx
# function to calculate LY lost due to nonCeCx death, account for relative survival, i.e., not all people will reach (survive) till  maxAge
fun_LYlost_nonCeCx_relSurv = function(disc_factor, deathAge, relSurv_mat, maxAge=85){
	if (deathAge==maxAge){
		LYlost_temp = 0
	} else{	
		# same as fun_LYlost_nonCeCx() if relSurv_temp = rep(1, (maxAge-deathAge))
		relSurv_temp = relSurv_mat[deathAge, (deathAge+1):maxAge]
		LYlost_temp = rep(1, (maxAge-deathAge))
		LYlost_temp = sum(LYlost_temp * relSurv_temp * disc_factor^(0:(length(LYlost_temp)-1)));
	}
	return( LYlost_temp )
} # fun_LYlost_nonCeCx_relSurv


nonCeCx_agegp = cbind(seq(10, 85, by=5), c((seq(15, 85, by=5)-1), 85))
nonCeCx_agegp_num = nrow(nonCeCx_agegp)
nonCeCx_agegp_midpt = apply(nonCeCx_agegp,1, mean)


# population data for F/M by age group, do once only; recall: Female = 1, Male = 2
# give number of pop in each age group
# by-cohort pop data
numAgeYr_AgeGp = 5; # number of age-year / interval of each age group
popAgeGp_cohort_nonCeCx_FM_list = lapply(list(pop_proj_CeCx, pop_proj_VaccBoy), function(x) numAgeYr_AgeGp*x)


discRate_cost = 0.03;
discRate_health = 0.03;
discRate_cost_factor = 1/(1+discRate_cost);
discRate_health_factor = 1/(1+discRate_health);
discRate_cost_factor_vec = (1/(1+discRate_cost))^((1:Yr_eval_num)-1)
discRate_health_factor_vec = (1/(1+discRate_health))^((1:Yr_eval_num)-1)

# LY lost for a death case at deathAge, _accounted_ for discRate_health
# LYlost_nonCeCxDeath_byAgeGp_disc = fun_LYlost_nonCeCx(disc_factor=discRate_health_factor, deathAge=nonCeCx_agegp_midpt)
# LY lost for a death case at deathAge, accounted for relative survival
LYlost_nonCeCxDeath_byAgeGp_disc_relSurv_FM_list = sapply(1:2, simplify=FALSE, function(sex) sapply(nonCeCx_agegp_midpt, function(xmidpt) fun_LYlost_nonCeCx_relSurv(deathAge=xmidpt, maxAge=85, disc_factor=discRate_health_factor, relSurv_mat=data_lifetab_FM_list[[sex]])))
discRate_health_factor_vec_cohort_LYlost_nonCeCxDeath = (1/(1+discRate_health))^(discRate_import_cohortInfo$Eval_VaccYr-1) # for calculation by-cohort

irow_ref_noVacc_cohort = which.max(discRate_import_cohortInfo$CohortAge) # the oldest cohort

# discount rate for costVacc for boys
yrSince_Yr_GNVstart = with(data_PSAsumm_init, Vage_girls-CohortEvaAge + Vacc_EvaYear - 1)
yrSince_Yr_GNVstart[yrSince_Yr_GNVstart<0] = 0
costVacc_discRate_cost_factor_vec = discRate_cost_factor^yrSince_Yr_GNVstart;


iPSA_vec = 1:100
iPSA_num = max(iPSA_vec)


## whether to include PSA for individual cost and QALY items for CeCx and nonCeCx
PSA_costQALY_YN = TRUE
PSA_costQALY_num = iPSA_num;
PSA_costQALY_num = ifelse(PSA_costQALY_YN, PSA_costQALY_num, 1)
nonCeCx_vary_cxStageProp_YN = TRUE;

matchname_PSA_nonCeCx = c("F_vagina_vulva"="F_VV", "F_OPC"="F_OPC", "F_anus"="F_anus", "M_OPC"="M_OPC", "M_penis"="M_penis", "M_anus"="M_anus")
if (!identical(nonCeCx_vec, names(matchname_PSA_nonCeCx))){ print('check the order of matchname_PSA_nonCeCx') }
input_cost_PSA_nonCeCx = sapply(matchname_PSA_nonCeCx, simplify=FALSE, function(xnonCeCx) input_PSA_nonCeCx_param[, paste0("cost_",xnonCeCx)])
input_util5yr_PSA_nonCeCx = sapply(matchname_PSA_nonCeCx, simplify=FALSE, function(xnonCeCx) input_PSA_nonCeCx_param[, paste0("util5yr_",xnonCeCx)])


	# jcol in PSA_Summ
	jcol_PSAsumm_Cost = c("sumCostTx_Cx1", "sumCostTx_Cx2", "sumCostTx_Cx3", "sumCostTx_Cx4", "sumCostCIN2Tx", "sumCostCIN3Tx", "sumCostvisitScreen", "sumCostHPVTest", "sumCostCytoTest", "sumCostColpoTest", "sumCostPalliCare_Cx1", "sumCostPalliCare_Cx2", "sumCostPalliCare_Cx3", "sumCostPalliCare_Cx4")
	jcol_PSAsumm_Cost_screening = 7:10;
	jcol_PSAsumm_Cost_treatCancer = 1:4;
	jcol_PSAsumm_Cost_treatCancerPalliCare = 11:14;
	jcol_PSAsumm_Cost_treatCIN = 5:6;
	
	jcol_PSAsumm_QALY = c("QoL_NORMAL", "QoL_Init_NormalCyto", "QoL_ASCUS", "QoL_Colpo_Normal", "QoL_CIN1", "QoL_CIN2", "QoL_CIN3", "QoL_Cx1", "QoL_Cx2", "QoL_Cx3", "QoL_Cx4", "QoL_CAsurv_Cx1", "QoL_CAsurv_Cx2", "QoL_CAsurv_Cx3", "QoL_CAsurv_Cx4", "QoL_HPVneg", "QoL_HPVpos", "QoL_ASCUS_HPVneg", "QoL_RepeatNormalCyto", "QoL_LSIL")
	jcol_PSAsumm_QALY_normal = 1;
	jcol_PSAsumm_QALY_screening = c(2:4, 16:20);
	jcol_PSAsumm_QALY_treatCancer = 8:15;
	jcol_PSAsumm_QALY_treatCIN = 5:7;


# get cohort age in PSAsumm
jcol_CohortAge_Summ = 56;
CohortAge_T = data_PSAsumm_init[,jcol_CohortAge_Summ];
nrow_CohortAge = length(CohortAge_T);

# screen uptake per cohort age
pScrRate_vec = rep(0, nrow(data_PSAsumm_init)); # set all 
pScrRate_vec[which(CohortAge_T>=pScrRate_age_beforeCSP)] = pScrRate_base_beforeCSP
pScrRate_vec[which(pScrRate_age_afterCSP<=CohortAge_T & CohortAge_T<pScrRate_age_beforeCSP)] = approx(x=c(pScrRate_age_afterCSP, pScrRate_age_beforeCSP), y=c(pScrRate_base_afterCSP, pScrRate_base_beforeCSP), xout=(pScrRate_age_afterCSP):(pScrRate_age_beforeCSP-1))$y
pScrRate_vec[which(CohortAge_T<=pScrRate_age_afterCSP)] = pScrRate_base_afterCSP
pScrRate_vec[which(CohortAge_T<=pScrRate_age_progChange)] = pScrRate_base_progChange
irow_pick = which(data_PSAsumm_init$Vacc_EvaYear<=Yr_eval_num)
pScrRate_vec = pScrRate_vec[irow_pick];

# costs and QALYs relating to screening
pScrRate_mat_Cost = do.call(cbind, rep(list(pScrRate_vec), length(jcol_PSAsumm_Cost)))
pScrRate_mat_QALY = do.call(cbind, rep(list(pScrRate_vec), length(jcol_PSAsumm_QALY)))


# separate cost and QALY for CeCx and nonCeCx
# a single overall value for each VC/Vdur/iPSA
output_template_overall_count = array(0, dim=c(VC_boys_num, Vdur_num, iPSA_num)); 
output_template_overall_PSAcostQALY = array(0, dim=c(VC_boys_num, Vdur_num, iPSA_num*PSA_costQALY_num));

Cost_array = output_template_overall_PSAcostQALY
QALY_array = Cost_array; # QALY for cervical cancer related, from Cpp with screening and the transmission dynamic model
Cost_vaccGirls_array = Cost_array;

# breakdown variables for CeCx
Cost_CeCx_screening_array = Cost_array;
Cost_CeCx_treatCancer_array = Cost_array;
Cost_CeCx_treatCancerPalliCare_array = Cost_array;
Cost_CeCx_treatCIN_array = Cost_array;
QALY_CeCx_normal_array = QALY_array;
QALY_CeCx_screening_array = QALY_array;
QALY_CeCx_treatCancer_array = QALY_array;
QALY_CeCx_treatCIN_array = QALY_array;


## estimate the changes of non-cervical cancer incidence
# comparing by time; first get the count of rate x pop_size, then compare count at each time
count_nonCeCx_array = output_template_overall_count;

# first compare by rate, then multiply by pop_size
Caseprev_nonCeCx_array = count_nonCeCx_array # case prevented, compared to time1
Costsave_nonCeCx_array = Cost_array # cost saved, compared to time1
QALYgain_nonCeCx_array = QALY_array; # QALY gain, compared to time1
LYgain_nonCeCxDeath_array = QALY_array # LY gain, compared to time1


# by nonCeCxType
# cancer cases
count_nonCeCx_array_byType = rep(list(count_nonCeCx_array), nonCeCx_num);
Caseprev_nonCeCx_array_byType = count_nonCeCx_array_byType;
# cost/QALY
Costsave_nonCeCx_array_byType = rep(list(Cost_array), nonCeCx_num);
QALYgain_nonCeCx_array_byType = Costsave_nonCeCx_array_byType
LYgain_nonCeCxDeath_array_byType = Costsave_nonCeCx_array_byType


# age-standardized rate; byTime, one more dimension
ASR_nonCeCx_array_byType_byTime = rep(list(array(0, dim=c(VC_boys_num, Vdur_num, iPSA_num, Yr_eval_num_bycohort))), nonCeCx_num) # by-cohort
Rate_nonCeCx_byTime_byAge_byType = rep(list(array(list(), dim=c(VC_boys_num, Vdur_num, iPSA_num))), nonCeCx_num) 
Ratediff_nonCeCx_byTime_byAge_byType = Rate_nonCeCx_byTime_byAge_byType
count_nonCeCx_byTime_byAge_byType = Rate_nonCeCx_byTime_byAge_byType
Caseprev_nonCeCx_byTime_byAge_byType = Rate_nonCeCx_byTime_byAge_byType


	# outcomes when no changes in MaleCx from MSM subgroup for FOV
	count_nonCeCx_array_noMaleCxFOV = count_nonCeCx_array;
	Caseprev_nonCeCx_array_noMaleCxFOV = Caseprev_nonCeCx_array;

	Costsave_nonCeCx_array_noMaleCxFOV = Costsave_nonCeCx_array;
	QALYgain_nonCeCx_array_noMaleCxFOV = QALYgain_nonCeCx_array;
	LYgain_nonCeCxDeath_array_noMaleCxFOV = LYgain_nonCeCxDeath_array
	
	Costsave_nonCeCx_array_byType_noMaleCxFOV = Costsave_nonCeCx_array_byType;
	QALYgain_nonCeCx_array_byType_noMaleCxFOV = QALYgain_nonCeCx_array_byType
	LYgain_nonCeCxDeath_array_byType_noMaleCxFOV = LYgain_nonCeCxDeath_array_byType

	# cancer cases
	count_nonCeCx_array_byType_noMaleCxFOV = count_nonCeCx_array_byType;
	Caseprev_nonCeCx_array_byType_noMaleCxFOV = Caseprev_nonCeCx_array_byType;

	# age-standardized rate; byTime, one more dimension
	ASR_nonCeCx_array_byType_byTime_noMaleCxFOV = ASR_nonCeCx_array_byType_byTime



# cost to vaccinate boys; because assuming that vaccine uptake is unchanged for girls, the cost will not be needed
Cost_vaccBoys_array = Cost_array;

# PSACost_visitVacc = PSACost_GPconsult / NumDose_Vacc + CountTime_Vacc * PSACost_PatientTime / 2 + CountTransport_Vacc * PSACost_Transport,
Cost_Vacc_perdose = 1100; # not to change, cost in Cpp (PSAsumm)
Cost_GPconsult = 250;
Cost_Transport =  50;
Cost_Vacc_use_base = Cost_Vacc_perdose + Cost_GPconsult + Cost_Transport; # include Cost_Transport per Cpp setting but this may be changed later

numDose_Vacc_Cpp = 1; # number of dose for each individual
numDose_Vacc = 2; # number of dose for each individual
# assume the same number of vaccination dose for each individual for both boys and girls


	# _genital warts_
	# variables related to gwarts starting from here
	# incidence rate in 100,000 person-years, from Lin BMC Infect Dis 2010
	gwarts_incid_100kpersonyr = c(Female=124.86, Male=292.19)
	gwarts_cost_org = 1100 # ref Cheung 2023
	gwarts_utilityloss_org = 0.018 # per episode, ref. Woodhall 2011, also used in Jit BMJ 2011;
	gwarts_propHPV611 = 0.9; # attributable to HPV-6/11
	gwarts_agerange_1yr = 20:69;
	gwarts_agerange = range(gwarts_agerange_1yr);
	gwarts_agerange_startage_5yr = floor(gwarts_agerange/5)*5;
	gwarts_age_CostQALY = 30; # assuming that age cohorts older than this age at the current year would have minimum effect due to vaccination (i.e., up to age 34 for age grouping 30-34); not related to the cohorts in younger ages in subsequent years
	gwarts_idx_5yrAgegp_fromAge10 = match(gwarts_agerange_startage_5yr, seq(10, 85, by=5));
	gwarts_idx_5yrAgegp_fromAge10 = gwarts_idx_5yrAgegp_fromAge10[1]:gwarts_idx_5yrAgegp_fromAge10[2]
	gwarts_CostQALY_idx_5yrAgegp_fromAge10 = match(gwarts_age_CostQALY, seq(10, 85, by=5)); 
	gwarts_Caseprev_age_vec = c(20, 15); # when considering age cohorts that have a higher vaccine uptake (pop-based and opportunistic vaccination); i.e., age cohorts of 20-24 for age=20, and 15-19 for age=15
	gwarts_Caseprev_age_vec_num = length(gwarts_Caseprev_age_vec);
	gwarts_Caseprev_idx_5yrAgegp_fromAge10_age_list = sapply(gwarts_Caseprev_age_vec, simplify=FALSE, function(xage) match(xage, seq(10, 85, by=5))); 
	
	gwarts_stdpop_pick = stdpop_pick[gwarts_idx_5yrAgegp_fromAge10]
	gwarts_stdpop_pick_pct = gwarts_stdpop_pick/100; # stdpop in pct
	gwarts_popAgeGp_cohort = lapply(popAgeGp_cohort_nonCeCx_FM_list, function(x) x*length(gwarts_idx_5yrAgegp_fromAge10))

	gwarts_disc_CostQALY = discRate_import[, gwarts_CostQALY_idx_5yrAgegp_fromAge10]
	gwarts_disc_Caseprev_age_list = sapply(gwarts_Caseprev_idx_5yrAgegp_fromAge10_age_list, simplify=FALSE, function(xage) discRate_import[, xage])
	
	gwarts_cost_PSAvec = rep(gwarts_cost_org, PSA_costQALY_num)
	gwarts_utilityloss_PSAvec = rep(gwarts_utilityloss_org, PSA_costQALY_num)

	# compared to noVacc in girls/boys, similar to nonCeCx items
	Caseprev_gwarts_array_bysex = rep(list(count_nonCeCx_array), sex_num); 
	Costsave_gwarts_array_bysex = rep(list(Cost_array), sex_num);
	QALYgain_gwarts_array_bysex = Costsave_gwarts_array_bysex;
	# initial number of cases of gwarts, should not be affected by VC and Vdur; to calculate the relative reduction of gwarts by vaccines
	Caseinit_allAges_gwarts_array_bysex = rep(rep(list(NULL), iPSA_num), sex_num); # cases among all ages
	Caseinit_gwarts_array_bysex = rep(list(array(0, dim=c(1, 1, iPSA_num))), sex_num); # among age cohorts that are affected by vaccines only
	relative_Caseprev_gwarts_array_bysex = Caseprev_gwarts_array_bysex;  # relative changes, calculate outside the for-loops
	# referring to gwarts_age_caseprev
	Caseprev_gwarts_array_bysex_age_list = rep(list(Caseprev_gwarts_array_bysex), gwarts_Caseprev_age_vec_num);
	Caseinit_gwarts_array_bysex_age_list = rep(list(Caseinit_gwarts_array_bysex), gwarts_Caseprev_age_vec_num);
	relative_Caseprev_gwarts_array_bysex_age_list = rep(list(relative_Caseprev_gwarts_array_bysex), gwarts_Caseprev_age_vec_num);



## start retrieving results
for (iVdur in 1:Vdur_num){
for (iVC in 1:VC_boys_num){
time1 = Sys.time()

## cost/QALY per PSA
# CeCx
Cost_array_iPSA = rep(list(rep(0, PSA_costQALY_num)), iPSA_num)
QALY_array_iPSA = Cost_array_iPSA
# vaccCost
Cost_vaccGirls_array_iPSA = Cost_array_iPSA
Cost_vaccBoys_array_iPSA = Cost_array_iPSA

# nonCeCx, byType
Costsave_nonCeCx_array_byType_iPSA = rep(list(Cost_array_iPSA), nonCeCx_num)
QALYgain_nonCeCx_array_byType_iPSA = Costsave_nonCeCx_array_byType_iPSA;
LYgain_nonCeCxDeath_array_byType_iPSA = Costsave_nonCeCx_array_byType_iPSA;

# nonCeCx, noMaleCx_FOV, byType
Costsave_nonCeCx_array_byType_noMaleCxFOV_iPSA = Costsave_nonCeCx_array_byType_iPSA;
QALYgain_nonCeCx_array_byType_noMaleCxFOV_iPSA = Costsave_nonCeCx_array_byType_iPSA;
LYgain_nonCeCxDeath_array_byType_noMaleCxFOV_iPSA = Costsave_nonCeCx_array_byType_iPSA;


# breakdown variables for CeCx
Cost_CeCx_screening_array_iPSA = Cost_array_iPSA;
Cost_CeCx_treatCancer_array_iPSA = Cost_array_iPSA;
Cost_CeCx_treatCancerPalliCare_array_iPSA = Cost_array_iPSA;
Cost_CeCx_treatCIN_array_iPSA = Cost_array_iPSA;
QALY_CeCx_normal_array_iPSA = QALY_array_iPSA;
QALY_CeCx_screening_array_iPSA = QALY_array_iPSA;
QALY_CeCx_treatCancer_array_iPSA = QALY_array_iPSA;
QALY_CeCx_treatCIN_array_iPSA = QALY_array_iPSA;


# genital warts
Costsave_gwarts_array_bysex_iPSA = rep(list(Cost_array_iPSA), sex_num)
QALYgain_gwarts_array_bysex_iPSA = Costsave_gwarts_array_bysex_iPSA

for (iPSA in iPSA_vec){

xVC = VC_boys_vec[iVC]
xVdur = Vdur_vec[iVdur]

if (iPSA%%10 == 1){
	print( sprintf("xVdur = %d, xVC = %d, iPSA = %d", xVdur, xVC, iPSA ))
}


# cervical cancer
	data_PSAsumm = data_PSAsumm_list[[iVdur]][[iVC]][[iPSA]]
	data_PSAsumm_noS = data_PSAsumm_noS_list[[iVdur]][[iVC]][[iPSA]]

	irow_pick = which(data_PSAsumm$Vacc_EvaYear<=Yr_eval_num)

	# assume fixed vaccine cost within each iPSA on natural history, i.e., no PSA_costQALY for vaccine cost

	# vaccination cost for girls is independent of screen uptake
	Cost_vaccGirls_temp_wS = with(data_PSAsumm, (numDose_Vacc/numDose_Vacc_Cpp*sumCostVacc)[irow_pick]) # sumCostVacc is discounted
	Cost_vaccGirls_temp_wS = Cost_vaccGirls_temp_wS / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;
	sum_Cost_vaccGirls_temp = sum(Cost_vaccGirls_temp_wS);
	Cost_vaccGirls_array_iPSA[[iPSA]] = rep(sum_Cost_vaccGirls_temp, PSA_costQALY_num)
	
	# cost to vaccinate boys
	if (xVC==0){
		sum_Cost_vaccBoys_temp = 0

	} else{
		TotDose_girls_cohort = data_PSAsumm$sumTotVacc[irow_pick]
		vaccProg_girls_YN = 1*(data_PSAsumm$sumCostVacc[irow_pick]!=0) # there was background (opportunistic) vaccination for older girls, so consider those vaccinated under the program only
		TotDose_boys = TotDose_girls_cohort/VC_girls * xVC * numDose_Vacc/numDose_Vacc_Cpp; # VC_girls and xVC are both in the same magnitude (percentage), so no need to "* or / 100"

		sum_Cost_vaccBoys_temp = sum( Cost_Vacc_use_base * TotDose_boys * vaccProg_girls_YN * costVacc_discRate_cost_factor_vec[irow_pick] / Cpp_cohortsize * pop_proj_VaccBoy * SubjAtRisk_vec_CeCx);
	} 
	Cost_vaccBoys_array_iPSA[[iPSA]] = rep(sum_Cost_vaccBoys_temp, PSA_costQALY_num)


if (PSA_costQALY_YN==FALSE){
	# with screening
	Cost_temp_wS = with(data_PSAsumm, (sumCostTest + sumCostTx)[irow_pick])
	
	Cost_temp_wS = Cost_temp_wS / Cpp_cohortsize * pop_proj_CeCx
	QALY_temp_wS = data_PSAsumm$sumQALY[irow_pick] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;

	# no screening
	Cost_temp_noS = with(data_PSAsumm_noS, (sumCost - sumCostVacc)[irow_pick])
	
	Cost_temp_noS = Cost_temp_noS / Cpp_cohortsize * pop_proj_CeCx
	QALY_temp_noS = data_PSAsumm_noS$sumQALY[irow_pick] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;

	Cost_temp = screenUptake*Cost_temp_wS + (1-screenUptake)*Cost_temp_noS;
	QALY_temp = screenUptake*QALY_temp_wS + (1-screenUptake)*QALY_temp_noS;

	Cost_array_iPSA[[iPSA]] = sum(Cost_temp) # rep(sum(Cost_temp), PSA_costQALY_num);
	QALY_array_iPSA[[iPSA]] = sum(QALY_temp) # rep(sum(QALY_temp), PSA_costQALY_num);
	
	# PSA_costQALY_YN == FALSE
} else if (PSA_costQALY_YN == TRUE){
	 # update PSA Cost & QALY 
	# Cost = costTest + costTx = sumCostTx_Cx1, sumCostTx_Cx2, sumCostTx_Cx3, sumCostTx_Cx4, sumCostCIN2Tx, sumCostCIN3Tx, sumCostvisitScreen, sumCostHPVTest, sumCostCytoTest, sumCostColpoTest, sumCostPalliCare_Cx1, sumCostPalliCare_Cx2, sumCostPalliCare_Cx3, sumCostPalliCare_Cx4
	# QALY = QoL_NORMAL, QoL_Init_NormalCyto, QoL_ASCUS, QoL_Colpo_Normal, QoL_CIN1, QoL_CIN2, QoL_CIN3, QoL_Cx1, QoL_Cx2, QoL_Cx3, QoL_Cx4, QoL_CAsurv_Cx1, QoL_CAsurv_Cx2, QoL_CAsurv_Cx3, QoL_CAsurv_Cx4, QoL_HPVneg, QoL_HPVpos, QoL_ASCUS_HPVneg, QoL_RepeatNormalCyto, QoL_LSIL
	
	# change sumCostTest, sumCostTx, and sumQALY in data_PSAsumm, data_PSAsumm_noS
	
	
	Cost_PSAsumm_temp_wS = data_PSAsumm[irow_pick, jcol_PSAsumm_Cost] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;
	Cost_PSAsumm_temp_noS = data_PSAsumm_noS[irow_pick, jcol_PSAsumm_Cost] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;
	Cost_PSAsumm_temp = pScrRate_mat_Cost*Cost_PSAsumm_temp_wS + (1-pScrRate_mat_Cost)*Cost_PSAsumm_temp_noS; # pScrRate_mat_Cost
	Cost_PSAsumm_temp = colSums(Cost_PSAsumm_temp)
	
	QALY_PSAsumm_temp_wS = data_PSAsumm[irow_pick, jcol_PSAsumm_QALY] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;
	QALY_PSAsumm_temp_noS = data_PSAsumm_noS[irow_pick, jcol_PSAsumm_QALY] / Cpp_cohortsize * pop_proj_CeCx * SubjAtRisk_vec_CeCx;
	QALY_PSAsumm_temp = pScrRate_mat_QALY*QALY_PSAsumm_temp_wS + (1-pScrRate_mat_QALY)*QALY_PSAsumm_temp_noS; # pScrRate_mat_QALY
	QALY_PSAsumm_temp = colSums(QALY_PSAsumm_temp)
	
	Cost_temp_iPSA = repmat(Cost_PSAsumm_temp[jcol_PSAsumm_Cost], m=1, n=PSA_costQALY_num) * input_PSA_CeCx_param_ratio[, jcol_PSAsumm_Cost]
	Cost_array_iPSA[[iPSA]] = rowSums(Cost_temp_iPSA)
	
	QALY_temp_iPSA = repmat(QALY_PSAsumm_temp[jcol_PSAsumm_QALY], m=1, n=PSA_costQALY_num) * input_PSA_CeCx_param_ratio[, jcol_PSAsumm_QALY]
	QALY_array_iPSA[[iPSA]] = rowSums(QALY_temp_iPSA)
	
	# breakdown variables for CeCx, add "drop=FALSE" as jcol_PSAsumm may be one element only
	Cost_CeCx_screening_array_iPSA[[iPSA]] = rowSums(Cost_temp_iPSA[, jcol_PSAsumm_Cost_screening, drop=FALSE])
	Cost_CeCx_treatCancer_array_iPSA[[iPSA]] = rowSums(Cost_temp_iPSA[, jcol_PSAsumm_Cost_treatCancer, drop=FALSE])
	Cost_CeCx_treatCancerPalliCare_array_iPSA[[iPSA]] = rowSums(Cost_temp_iPSA[, jcol_PSAsumm_Cost_treatCancerPalliCare, drop=FALSE])
	Cost_CeCx_treatCIN_array_iPSA[[iPSA]] = rowSums(Cost_temp_iPSA[, jcol_PSAsumm_Cost_treatCIN, drop=FALSE])
	QALY_CeCx_normal_array_iPSA[[iPSA]] = rowSums(QALY_temp_iPSA[, jcol_PSAsumm_QALY_normal, drop=FALSE])
	QALY_CeCx_screening_array_iPSA[[iPSA]] = rowSums(QALY_temp_iPSA[, jcol_PSAsumm_QALY_screening, drop=FALSE])
	QALY_CeCx_treatCancer_array_iPSA[[iPSA]] = rowSums(QALY_temp_iPSA[, jcol_PSAsumm_QALY_treatCancer, drop=FALSE])
	QALY_CeCx_treatCIN_array_iPSA[[iPSA]] = rowSums(QALY_temp_iPSA[, jcol_PSAsumm_QALY_treatCIN, drop=FALSE])	
	
	# PSA_costQALY_YN == TRUE
}
## end CeCx


	# prop of vaccinated, to estimate the case of no herd effect 
	F_enterVacc_cohort = data_PSAsumm$enterVaccGp
	F_enterVacc_cohort = F_enterVacc_cohort[irow_pick]
	M_enterVacc_cohort = rep(0, length(F_enterVacc_cohort))
	M_enterVacc_cohort[CohortAge_T[irow_pick]<=12] = xVC/100
	FM_enterVacc_cohort = list("F"=F_enterVacc_cohort, "M"=M_enterVacc_cohort)


## nonCeCx, non-cervical cancers
nonCeCx_temp_allZero = rep(0, nonCeCx_num);
# not related to cost/QALY
count_nonCeCx_temp = nonCeCx_temp_allZero
Caseprev_nonCeCx_temp = count_nonCeCx_temp;
ASR_nonCeCx_temp = rep(list(rep(0, Yr_eval_num)), nonCeCx_num)

for (inonCeCx in 1:nonCeCx_num){
	xnonCeCx = nonCeCx_vec[inonCeCx];
	str_nonCeCx = gsub("(F|M)_","",xnonCeCx)

	sex_pick = sex_vec[substr(xnonCeCx,1,1)]
	data_pop_cohort_nonCeCx = popAgeGp_cohort_nonCeCx_FM_list[[sex_pick]] # by-cohort

	# data_proj_nonCeCx => projected incidence rate; data_pop_nonCeCx => projected population size; count_nonCeCx => incidence count of nonCeCx
	# data_proj_nonCeCx, projected incidence rate per 100,000
	fname_proj_nonCeCx_pick = sprintf(fname_proj_nonCeCx_base, xVC, xVdur, xnonCeCx)
	shname = paste0("PSA_",iPSA)
	data_proj_nonCeCx = as.matrix( read.xlsx(paste0(folder_proj_nonCeCx, fname_proj_nonCeCx_pick), sheet=shname, colNames=FALSE) ) # convert to matrix

	data_pop_cohort_nonCeCx_inAgeGp = pracma::repmat(data_pop_cohort_nonCeCx, 1, ncol(data_proj_nonCeCx)) # pop data by cohort


	# projected nonCeCx should already be based on the cohort involved
	count_nonCeCx = data_proj_nonCeCx/100000 * data_pop_cohort_nonCeCx_inAgeGp; # by-cohort
	count_nonCeCx_temp[inonCeCx] = sum(count_nonCeCx);
	count_nonCeCx_byTime_byAge_byType[[inonCeCx]][[iVC, iVdur, iPSA]] = count_nonCeCx;

	Rate_nonCeCx_byTime_byAge_byType[[inonCeCx]][[iVC, iVdur, iPSA]] = data_proj_nonCeCx

	# ASR by time, unrelated to population size
	ASR_nonCeCx_temp[[inonCeCx]] = apply(data_proj_nonCeCx,1, function(x) sum(x*stdpop_pick/100))

	# Casediff, Costsave, QALYgain, LYgain, based on difference in rate then multiplied by popsize
	# prevented case of nonCeCx when compared to time1; "relative" count
	# preventedCount_nonCeCx_byTime_byAge_temp = pracma::repmat(as.vector(count_nonCeCx[1,]), nrow(count_nonCeCx), 1) - as.matrix(count_nonCeCx) # compare count_nonCeCx at time1; population size data_pop_nonCeCx may affect
	# calculate the change in incidRate_nonCeCx first, then apply for data_pop_nonCeCx
	
	# Casediff, difference in number of cases, (time_1- time_i) x pop_size_time_i; each age group per 100,000 individuals
	
	irow_ref_noVacc = irow_ref_noVacc_cohort;
	Ratediff_nonCeCx_byTime_byAge_temp = pracma::repmat(data_proj_nonCeCx[irow_ref_noVacc,], nrow(data_proj_nonCeCx), 1) - as.matrix(data_proj_nonCeCx); # difference in incidence rate
	Ratediff_nonCeCx_byTime_byAge_temp = pmax(Ratediff_nonCeCx_byTime_byAge_temp, 0) # order of x/y in pmax(x,y) would make the output different

	
	# diff in rate x popsize
	# Casediff_nonCeCx_byTime_byAge_temp = Ratediff_nonCeCx_byTime_byAge_temp/100000 * data_pop_nonCeCx; # cross-sectional
	Casediff_nonCeCx_byTime_byAge_temp = Ratediff_nonCeCx_byTime_byAge_temp/100000 * data_pop_cohort_nonCeCx_inAgeGp; # by-cohort

	Caseprev_nonCeCx_temp[inonCeCx] = sum(Casediff_nonCeCx_byTime_byAge_temp)
	Caseprev_nonCeCx_byTime_byAge_byType[[inonCeCx]][[iVC, iVdur, iPSA]] = Casediff_nonCeCx_byTime_byAge_temp;
	
	Ratediff_nonCeCx_byTime_byAge_byType[[inonCeCx]][[iVC, iVdur, iPSA]] = Ratediff_nonCeCx_byTime_byAge_temp


	# clinical outcomes, no discounting;
	# case count; case prevented, ASR
	count_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, iPSA] = count_nonCeCx_temp[inonCeCx]
	Caseprev_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, iPSA] = Caseprev_nonCeCx_temp[inonCeCx]
	ASR_nonCeCx_array_byType_byTime[[inonCeCx]][iVC, iVdur, iPSA, ] = ASR_nonCeCx_temp[[inonCeCx]]

	# cost / QALY
	Costsave_nonCeCx_temp_type = sum(Casediff_nonCeCx_byTime_byAge_temp * discRate_import_cost)
	QALYgain_nonCeCx_temp_type = sum(Casediff_nonCeCx_byTime_byAge_temp * discRate_import_health)

	LYlost_nonCeCxDeath_byAgeGp_temp = LYlost_nonCeCxDeath_byAgeGp_disc_relSurv_FM_list[[sex_pick]]
	LYlost_nonCeCxDeath_byAgeGp_timehor = pracma::repmat(LYlost_nonCeCxDeath_byAgeGp_temp, nrow(data_proj_nonCeCx), 1)
	LYgain_nonCeCxDeath_temp_type = sum(rowSums( as.matrix(Casediff_nonCeCx_byTime_byAge_temp * LYlost_nonCeCxDeath_byAgeGp_timehor) ) * discRate_health_factor_vec_cohort_LYlost_nonCeCxDeath)
	LYgain_nonCeCxDeath_temp_type = (1-nonCeCx_5yrRS[xnonCeCx]) * LYgain_nonCeCxDeath_temp_type
	

	# save the outputs per iVC, iVdur, iPSA
	# by Casediff_nonCeCx
	if (PSA_costQALY_YN==FALSE){
		
		Costsave_nonCeCx_temp_type = nonCeCx_cost_raw[str_nonCeCx] * Costsave_nonCeCx_temp_type;
		QALYgain_nonCeCx_temp_type = nonCeCx_QLost_5yr[xnonCeCx] * QALYgain_nonCeCx_temp_type;

		# PSA_costQALY_YN==FALSE
	} else if (PSA_costQALY_YN==TRUE){

		Costsave_nonCeCx_temp_type = input_cost_PSA_nonCeCx[[xnonCeCx]] * Costsave_nonCeCx_temp_type;
		QALYgain_nonCeCx_temp_type = input_util5yr_PSA_nonCeCx[[xnonCeCx]] * QALYgain_nonCeCx_temp_type;
		LYgain_nonCeCxDeath_temp_type = rep( LYgain_nonCeCxDeath_temp_type, PSA_costQALY_num) # no PSA on relative survival (RS), so rep per PSA_costQALY_num

		# PSA_costQALY_YN==TRUE
	}
	# save to iPSA
	Costsave_nonCeCx_array_byType_iPSA[[inonCeCx]][[iPSA]] = Costsave_nonCeCx_temp_type;
	QALYgain_nonCeCx_array_byType_iPSA[[inonCeCx]][[iPSA]] = QALYgain_nonCeCx_temp_type;
	LYgain_nonCeCxDeath_array_byType_iPSA[[inonCeCx]][[iPSA]] = LYgain_nonCeCxDeath_temp_type
	

	# sensitivity analysis to account for MaleCx among MSM subgroup
	if (nonCeCx_vec[inonCeCx] %in% Mcancer_adj_vec){ # nonCeCx remains unchanged (& xVC%in%xVC_MaleCx_vec)
	
		prop_MaleCx_MSM = c("F_vagina_vulva" = NA, "F_OPC" = NA, "F_anus" = NA, "M_OPC" = prop_MOPC_MSM, "M_penis" = prop_Mpenis_MSM, "M_anus" = prop_Manus_MSM);

		p_HPVattrib = data_HPVattrib[iPSA, xnonCeCx]

		# CxInc with vaccination
		Cx_Vacc_wgt = data_proj_nonCeCx
		Cx_noVacc = pracma::repmat(Cx_Vacc_wgt[irow_ref_noVacc,], n=nrow(Cx_Vacc_wgt), m=1)
		Cx_Vacc_noAttrib = fun_CxInc_noAttrib(Cx_Vacc_wgt, Cx_noVacc, p_HPVattrib)
		
		# CxInc for male, no female vaccination
		fname_proj_noMaleCxFOV_pick = sprintf(fname_proj_noMaleCxFOV_base, xVC, xVdur, xnonCeCx)
		shname = paste0("PSA_",iPSA)
		Cx_MSM_wgt = as.matrix( read.xlsx(paste0(folder_proj_noMaleCxFOV, fname_proj_noMaleCxFOV_pick), sheet=shname, colNames=FALSE) ) # convert to matrix
		Cx_MSM_noAttrib = fun_CxInc_noAttrib(Cx_MSM_wgt, Cx_noVacc, p_HPVattrib)
		
		data_proj_nonCeCx_noMaleCxFOV = fun_CxInc_adjpropMSM_Vacc_MSM(Cx_noVacc, Cx_Vacc_noAttrib, Cx_MSM_noAttrib, p_HPVattrib, prop_MaleCx_MSM[xnonCeCx])


		# update the calculation using data_proj_nonCeCx_noMaleCxFOV
		count_nonCeCx_noMaleCxFOV = data_proj_nonCeCx_noMaleCxFOV/100000 * data_pop_cohort_nonCeCx_inAgeGp;
		count_nonCeCx_temp_noMaleCxFOV = sum(count_nonCeCx_noMaleCxFOV);
		
		# Casediff, difference in number of cases, (time_1- time_i) x pop_size_time_i; each age group per 100,000 individuals
		Ratediff_nonCeCx_byTime_byAge_temp = pracma::repmat(data_proj_nonCeCx_noMaleCxFOV[irow_ref_noVacc,], nrow(data_proj_nonCeCx_noMaleCxFOV), 1) - as.matrix(data_proj_nonCeCx_noMaleCxFOV) # difference in incidence rate
		Ratediff_nonCeCx_byTime_byAge_temp = pmax(Ratediff_nonCeCx_byTime_byAge_temp, 0) # order of x/y in pmax(x,y) would make the output different
		
		# diff in rate x popsize
		Casediff_nonCeCx_byTime_byAge_temp = Ratediff_nonCeCx_byTime_byAge_temp/100000 * data_pop_cohort_nonCeCx_inAgeGp; # by-cohort 
		
		Caseprev_nonCeCx_temp_noMaleCxFOV = sum(Casediff_nonCeCx_byTime_byAge_temp)

		ASR_nonCeCx_temp_noMaleCxFOV = apply(data_proj_nonCeCx_noMaleCxFOV,1, function(x) sum(x*stdpop_pick/100))
		
		count_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA] = count_nonCeCx_temp_noMaleCxFOV
		Caseprev_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA] = Caseprev_nonCeCx_temp_noMaleCxFOV
		ASR_nonCeCx_array_byType_byTime_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA, ] = ASR_nonCeCx_temp_noMaleCxFOV


		# save to main array
		Costsave_nonCeCx_temp_type_noMaleCxFOV = sum(Casediff_nonCeCx_byTime_byAge_temp * discRate_import_cost)
		QALYgain_nonCeCx_temp_type_noMaleCxFOV = sum(Casediff_nonCeCx_byTime_byAge_temp * discRate_import_health)

		LYlost_nonCeCxDeath_byAgeGp_temp = LYlost_nonCeCxDeath_byAgeGp_disc_relSurv_FM_list[[sex_pick]]
		LYlost_nonCeCxDeath_byAgeGp_timehor = pracma::repmat(LYlost_nonCeCxDeath_byAgeGp_temp, nrow(data_proj_nonCeCx_noMaleCxFOV), 1)
		LYgain_nonCeCxDeath_temp_type_noMaleCxFOV = sum(rowSums(as.matrix(Casediff_nonCeCx_byTime_byAge_temp * LYlost_nonCeCxDeath_byAgeGp_timehor)) * discRate_health_factor_vec_cohort_LYlost_nonCeCxDeath)
		LYgain_nonCeCxDeath_temp_type_noMaleCxFOV = (1-nonCeCx_5yrRS[xnonCeCx]) * LYgain_nonCeCxDeath_temp_type_noMaleCxFOV
		

		if (PSA_costQALY_YN==FALSE){
			Costsave_nonCeCx_temp_type_noMaleCxFOV = nonCeCx_cost_raw[str_nonCeCx] * Costsave_nonCeCx_temp_type_noMaleCxFOV
			QALYgain_nonCeCx_temp_type_noMaleCxFOV = nonCeCx_QLost_5yr[xnonCeCx] * QALYgain_nonCeCx_temp_type_noMaleCxFOV

		} else if (PSA_costQALY_YN==TRUE){
			Costsave_nonCeCx_temp_type_noMaleCxFOV = input_cost_PSA_nonCeCx[[xnonCeCx]] * Costsave_nonCeCx_temp_type_noMaleCxFOV
			QALYgain_nonCeCx_temp_type_noMaleCxFOV = input_util5yr_PSA_nonCeCx[[xnonCeCx]] * QALYgain_nonCeCx_temp_type_noMaleCxFOV

			LYgain_nonCeCxDeath_temp_type_noMaleCxFOV = rep(LYgain_nonCeCxDeath_temp_type_noMaleCxFOV, PSA_costQALY_num)
		}
		# save to iPSA
		Costsave_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = Costsave_nonCeCx_temp_type_noMaleCxFOV;
		QALYgain_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = QALYgain_nonCeCx_temp_type_noMaleCxFOV
		LYgain_nonCeCxDeath_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = LYgain_nonCeCxDeath_temp_type_noMaleCxFOV
		
	# nonCeCx_vec[inonCeCx] %in% Mcancer_adj_vec
	} else{ 
		# for other nonCeCx, the values for noMaleCx_FOV would be the same
		count_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA] = count_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, iPSA]
		Caseprev_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA] = Caseprev_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, iPSA]
		ASR_nonCeCx_array_byType_byTime_noMaleCxFOV[[inonCeCx]][iVC, iVdur, iPSA, ] = ASR_nonCeCx_array_byType_byTime[[inonCeCx]][iVC, iVdur, iPSA, ]

		# save to main array
		Costsave_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = Costsave_nonCeCx_array_byType_iPSA[[inonCeCx]][[iPSA]]
		QALYgain_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = QALYgain_nonCeCx_array_byType_iPSA[[inonCeCx]][[iPSA]]
		LYgain_nonCeCxDeath_array_byType_noMaleCxFOV_iPSA[[inonCeCx]][[iPSA]] = LYgain_nonCeCxDeath_array_byType_iPSA[[inonCeCx]][[iPSA]]
		
	} # if-else MaleCx

	
} # for-inonCeCx

count_nonCeCx_array[iVC, iVdur, iPSA] = sum(sapply(count_nonCeCx_array_byType, function(x) x[iVC, iVdur, iPSA]))
Caseprev_nonCeCx_array[iVC, iVdur, iPSA] = sum(sapply(Caseprev_nonCeCx_array_byType, function(x) x[iVC, iVdur, iPSA]))

# for noMaleCx_FOV
count_nonCeCx_array_noMaleCxFOV[iVC, iVdur, iPSA] = sum(sapply(count_nonCeCx_array_byType_noMaleCxFOV, function(x) x[iVC, iVdur, iPSA]))
Caseprev_nonCeCx_array_noMaleCxFOV[iVC, iVdur, iPSA] = sum(sapply(Caseprev_nonCeCx_array_byType_noMaleCxFOV, function(x) x[iVC, iVdur, iPSA]))


	# genital warts
	dataread_HPVincid_list = dataread_HPVincid_list_Vdur_VC_PSA[[iVdur]][[iVC]][[iPSA]]
	
	# irow_ref_noVacc, use the same defined for nonCeCx
	jcol_dataread_HPVincid = 1:16; # indices within each HPV type
	ncol_dataread_HPVincid_eachHPV = 19; # number of columns for each HPV type; +3 columns for VaccYr, CohortAge, and CalendarYr
	iHPV_incidChange_gwarts = 1:2 # use average of HPV-16/18
	for (isex in 1:sex_num){
		x_gwarts_incid = gwarts_incid_100kpersonyr[isex]/100000 # gwarts of isex
		x_popAgeGp_cohort = gwarts_popAgeGp_cohort[[isex]]

		# HPVincid by age group
		data_HPVincid_temp_list = sapply(1:2, simplify=FALSE, function(iHPVtype) 
			dataread_HPVincid_list[[isex]][, (iHPVtype-1)*ncol_dataread_HPVincid_eachHPV+jcol_dataread_HPVincid, drop=FALSE]);
			
		if (herdProt_GWart_YN==FALSE){
			xx_isex = names(sex_vec)[isex]
			xVdur_FM_TT = c("F"=xVdur, "M"=xVdur)
			if (run_2F1M_YN == TRUE){
				xVdur_FM_TT = c("F"=100, "M"=xVdur_boys_1M)
			}
			xVdur_TT = xVdur_FM_TT[xx_isex]
		
			for (ilist_TT in 1:length(data_HPVincid_temp_list)){
				data_TT = data_HPVincid_temp_list[[ilist_TT]];
				data_TT_ref_noVacc = data_TT[irow_ref_noVacc,]
				jcol_eq0 = which(data_TT[1, ]==0)
				
				VEff_TT = data_VE_input[iPSA, ilist_TT]; # vaccine efficacy
				for (jcol_TT in jcol_dataread_HPVincid){
					if (jcol_TT %in% jcol_eq0){ next} # 
					adjfactor_waning = 1;
					if (xVdur_TT!=100){ # account for waning efficacy
						if (xVdur_TT==20){
							herdProt_set_TT = herdProt_set$Vdur_20
						} else if (xVdur_TT==30){
							herdProt_set_TT = herdProt_set$Vdur_30
						}
						if (jcol_TT > herdProt_set_TT["jcol_endVProtect"]){
							adjfactor_waning = herdProt_set_TT["relative_partProtect"]
						}
						if (jcol_TT > herdProt_set_TT["jcol_partProtect"]){
							adjfactor_waning = 0
						}
					} # end adjfactor_waning
					
					data_TT_jcol = data_TT[, jcol_TT]					
					data_TT_ref_noVacc_jcol = data_TT_ref_noVacc[jcol_TT]
					
					data_TT_jcol_new = ((1-FM_enterVacc_cohort[[xx_isex]]) + FM_enterVacc_cohort[[xx_isex]]*(1-VEff_TT*adjfactor_waning)) * data_TT_ref_noVacc_jcol
					data_TT[, jcol_TT] = data_TT_jcol_new
				} # for- jcol_TT
				
				data_HPVincid_temp_list[[ilist_TT]] = data_TT;
			} # for- ilist_TT
		}
		
		data_HPVincid_temp = Reduce("+", data_HPVincid_temp_list)
		std_HPVincid_temp = apply(data_HPVincid_temp, 1, function(x) sum(x[gwarts_idx_5yrAgegp_fromAge10]*gwarts_stdpop_pick_pct)) # standardized HPVincid
		relChange_HPVincid = std_HPVincid_temp/std_HPVincid_temp[irow_ref_noVacc]
		relReduce_HPVincid = pmax(1 - relChange_HPVincid, 0) # relative reduction, should be 0 at irow_ref_noVacc
		# Caseprev
		Caseprev_gwarts_temp_0 = (x_gwarts_incid * x_popAgeGp_cohort) * (gwarts_propHPV611*relReduce_HPVincid)
		Caseprev_gwarts_temp = Caseprev_gwarts_temp_0 * (gwarts_disc_CostQALY>0); # assuming that some age cohorts are more affected by vaccination
		Caseprev_gwarts_array_bysex[[isex]][iVC, iVdur, iPSA] = sum(Caseprev_gwarts_temp)
		# vary by age cohort
		for (iage_Caseprev_gw in 1:gwarts_Caseprev_age_vec_num){
			Caseprev_gwarts_temp = Caseprev_gwarts_temp_0 * (gwarts_disc_Caseprev_age_list[[iage_Caseprev_gw]]>0);
			Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][iVC, iVdur, iPSA] = sum(Caseprev_gwarts_temp)
		} # for- iage_Caseprev_gw
		# Costsave / QALYgain
		Costsave_gwarts_array_bysex_iPSA[[isex]][[iPSA]] = gwarts_cost_PSAvec * sum(Caseprev_gwarts_temp * gwarts_disc_CostQALY)
		QALYgain_gwarts_array_bysex_iPSA[[isex]][[iPSA]] = gwarts_utilityloss_PSAvec * sum(Caseprev_gwarts_temp * gwarts_disc_CostQALY)
		
		if (iVC==1 & iVdur==1){ # calculate the initial cases of gwarts
			Caseinit_gwarts_temp_0 = (x_gwarts_incid * x_popAgeGp_cohort); # no need to *std_HPVincid_temp because std_HPVincid_temp is used to estimate the relative changes only
			Caseinit_allAges_gwarts_array_bysex[[isex]][[iPSA]] = Caseinit_gwarts_temp_0; # all ages
			# at gwarts_age_CostQALY
			Caseinit_gwarts_temp = Caseinit_gwarts_temp_0 * (gwarts_disc_CostQALY>0); # same 
			Caseinit_gwarts_array_bysex[[isex]][1,1, iPSA] = sum(Caseinit_gwarts_temp); # no need to consider gwarts_propHPV611, because the initial cases should include all regardless of HPV types
			# at gwarts_age_Caseprev
			for (iage_Caseprev_gw in 1:gwarts_Caseprev_age_vec_num){
				Caseinit_gwarts_temp = Caseinit_gwarts_temp_0 * (gwarts_disc_Caseprev_age_list[[iage_Caseprev_gw]]>0); # same 
			Caseinit_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][1,1, iPSA] = sum(Caseinit_gwarts_temp);
			} # for- iage_Caseprev_gw
		}
	}

} # for-iPSA


# update, change from list() to vector for PSA_costQALY_YN
Cost_vaccGirls_array[iVC, iVdur, ] = do.call(c, Cost_vaccGirls_array_iPSA)
Cost_vaccBoys_array[iVC, iVdur, ] = do.call(c, Cost_vaccBoys_array_iPSA)

Cost_array[iVC, iVdur, ] = do.call(c, Cost_array_iPSA)
QALY_array[iVC, iVdur, ] = do.call(c, QALY_array_iPSA)
# breakdown variables for CeCx
Cost_CeCx_screening_array[iVC, iVdur, ] = do.call(c, Cost_CeCx_screening_array_iPSA)
Cost_CeCx_treatCancer_array[iVC, iVdur, ] = do.call(c, Cost_CeCx_treatCancer_array_iPSA)
Cost_CeCx_treatCancerPalliCare_array[iVC, iVdur, ] = do.call(c, Cost_CeCx_treatCancerPalliCare_array_iPSA)
Cost_CeCx_treatCIN_array[iVC, iVdur, ] = do.call(c, Cost_CeCx_treatCIN_array_iPSA)
QALY_CeCx_normal_array[iVC, iVdur, ] = do.call(c, QALY_CeCx_normal_array_iPSA)
QALY_CeCx_screening_array[iVC, iVdur, ] = do.call(c, QALY_CeCx_screening_array_iPSA)
QALY_CeCx_treatCancer_array[iVC, iVdur, ] = do.call(c, QALY_CeCx_treatCancer_array_iPSA)
QALY_CeCx_treatCIN_array[iVC, iVdur, ] = do.call(c, QALY_CeCx_treatCIN_array_iPSA)


# nonCeCx outcomes, by nonCeCx_Type
for (inonCeCx in 1:nonCeCx_num){
	Costsave_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, ] = do.call(c, Costsave_nonCeCx_array_byType_iPSA[[inonCeCx]]);
	QALYgain_nonCeCx_array_byType[[inonCeCx]][iVC, iVdur, ] = do.call(c, QALYgain_nonCeCx_array_byType_iPSA[[inonCeCx]]);
	LYgain_nonCeCxDeath_array_byType[[inonCeCx]][iVC, iVdur, ] = do.call(c, LYgain_nonCeCxDeath_array_byType_iPSA[[inonCeCx]]);
	
	Costsave_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ] = do.call(c, Costsave_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]]);
	QALYgain_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ] = do.call(c, QALYgain_nonCeCx_array_byType_noMaleCxFOV_iPSA[[inonCeCx]]);
	LYgain_nonCeCxDeath_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ] = do.call(c, LYgain_nonCeCxDeath_array_byType_noMaleCxFOV_iPSA[[inonCeCx]]);
}

# nonCeCx outcomes, overall
Costsave_nonCeCx_array[iVC, iVdur, ] = Reduce("+", lapply(Costsave_nonCeCx_array_byType, function(xlist) xlist[iVC, iVdur, ]))
QALYgain_nonCeCx_array[iVC, iVdur, ] = Reduce("+", lapply(QALYgain_nonCeCx_array_byType, function(xlist) xlist[iVC, iVdur, ]))
LYgain_nonCeCxDeath_array[iVC, iVdur, ] = Reduce("+", lapply(LYgain_nonCeCxDeath_array_byType, function(xlist) xlist[iVC, iVdur, ]))

Costsave_nonCeCx_array_noMaleCxFOV[iVC, iVdur, ] = Reduce("+", lapply(Costsave_nonCeCx_array_byType_noMaleCxFOV, function(xlist) xlist[iVC, iVdur, ]))
QALYgain_nonCeCx_array_noMaleCxFOV[iVC, iVdur, ] = Reduce("+", lapply(QALYgain_nonCeCx_array_byType_noMaleCxFOV, function(xlist) xlist[iVC, iVdur, ]))
LYgain_nonCeCxDeath_array_noMaleCxFOV[iVC, iVdur, ] = Reduce("+", lapply(LYgain_nonCeCxDeath_array_byType_noMaleCxFOV, function(xlist) xlist[iVC, iVdur, ]))


# genital warts
for (isex in 1:sex_num){
	Costsave_gwarts_array_bysex[[isex]][iVC, iVdur, ] = do.call(c, Costsave_gwarts_array_bysex_iPSA[[isex]])
	QALYgain_gwarts_array_bysex[[isex]][iVC, iVdur, ] = do.call(c, QALYgain_gwarts_array_bysex_iPSA[[isex]])
	
	relative_Caseprev_gwarts_array_bysex[[isex]][iVC, iVdur, ] = Caseprev_gwarts_array_bysex[[isex]][iVC, iVdur, ]/Caseinit_gwarts_array_bysex[[isex]][1,1, ];
	for (iage_Caseprev_gw in 1:gwarts_Caseprev_age_vec_num){
		relative_Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][iVC, iVdur, ] = Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][iVC, iVdur, ]/Caseinit_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][1,1, ];
	} # for- iage_Caseprev_gw
} # for- isex


time2 = Sys.time()
print(time2 - time1)
} # for-iVC
} # for-iVdur


outputdate = outputdate_use;

fname_RData_out = sprintf("%scalc_CEA_nonCeCx_%s.RData", output_folder_CEA_new, outputdate)
save( list=setdiff(ls(), c('data_PSAsumm_list', 'data_PSAsumm_noS_list', 'dataread_HPVincid_list_Vdur_VC_PSA')), file=fname_RData_out )


# different calculations for Cost and QALY for nonCeCx
# comparison of one-dose should be calculated differently
# 1-dose, when compared to 2-dose (a) lifelong protection but lower VE or (b) shorter protection but same VE


# print the summary stat
print_summstat = function(x, form="median", lv=90, big.mark=",", digit=1){
	if (form=="mean"){
		out = unlist(t.test(x, conf.level=ifelse(lv<1, lv, lv/100))[c("estimate", "conf.int")])
	} else {
		out = quantile(x, probs=c(0.5, 0.5+c(-1,1)*ifelse(lv<1, lv,lv/100)/2))
	}
	out = round(out, digit=digit)
	out = format(out, big.mark=",", scientific=FALSE, trim=TRUE) # digits in format() is not the decimal digit # or use trimws() to remove white space
	out = do.call(sprintf, c(fmt = "%s (%s, %s)", as.list(out)));
	return(out)
} # print_summstat


# reference vaccine cost, to calculate the relative change in TVC
Cost_Vacc_perdose_tender_plusAdmin_USD = 177; # new vaccination cost, vaccine cost based on the tender price plus admin expenses
Cost_Vacc_perdose_inclAdmin_ref = Cost_Vacc_perdose_tender_plusAdmin_USD*7.8; # new vaccination cost, including admin expenses
Cost_Vacc_perdose_AdminCost_assumed = 250;


# threshold vaccine cost (TVC) ratio 
# ** varying cost on vaccination
Cost_Vacc_perdose_new = Cost_Vacc_perdose_inclAdmin_ref - Cost_Vacc_perdose_AdminCost_assumed; 

consider_Cost_Transport_YN = FALSE; # no Cost_Transport for school-based program
num_subplot = VC_boys_num - 1;
subplot_iVC_YN = TRUE; # TRUE = 1 subplot for each iVC
if (!subplot_iVC_YN){ num_subplot = 1;}


# **recall: 2F1M vs 2F, numDose for girls and boys are fixed

numDose_Vacc_new = 1; # note. numDose_Vacc=2; originally vacc schedule = 2
numDose_Vacc_new_vec = 2; # girls, fixed
numDose_Vacc_new_Boys = 1; # boys, fixed # for 2F1M


# opt_width_default = 140 # options()$width
options(width = 200)
fname_output_txt = sprintf("%scalcCEA_output_%s.txt", output_folder_CEA_new, outputdate)
fname_output_pdf = gsub("txt", "pdf", fname_output_txt)
sink(file=fname_output_txt)
plot_PDF_YN = TRUE;
if (plot_PDF_YN) {pdf(fname_output_pdf, width=num_subplot*5, height=5)}



calc_gwarts_YN_vec = TRUE # c(FALSE, TRUE) # to include scenarios not accounting for genital warts at all (i.e., considering cervical and non-cervical cancers only)
gwarts_cost_use = gwarts_cost_org;
gwarts_utilityloss_use = gwarts_utilityloss_org;
print( sprintf("gwarts_cost_use = %.0f", gwarts_cost_use) )
print( sprintf("gwarts_utilityloss_use = %.3f", gwarts_utilityloss_use) )
gwarts_cost_ratio = gwarts_cost_use/gwarts_cost_org;
gwarts_utilityloss_ratio = gwarts_utilityloss_use/gwarts_utilityloss_org;


## adjust the proportion of MaleCx that can be prevented by FOV.
# prop_MaleCx_FOV = MSW / (MSW + MSM), MSW = men who have sex with women (heterosexual)
# prop_MaleCx_FOV=1 >> all MaleCx cancer can be prevented by FOV; prop_MaleCx_FOV=0 >> no MaleCx cancer can be prevented by FOV
run_SensAnaly_propMaleCx_YN_vec = list(FALSE, c(FALSE, TRUE))[[2]]; # set to TRUE to run the sensitivity analysis by prop_MaleCx_FOV or prop_MaleCx_MSM

inonCeCx_Manus = which(nonCeCx_vec=="M_anus");
inonCeCx_Mpenis = which(nonCeCx_vec=="M_penis");
inonCeCx_MOPC = which(nonCeCx_vec=="M_OPC");
inonCeCx_adjMSM_vec = NA;
if (adj_otherCx_MSM_YN){
	inonCeCx_adjMSM_vec = c(inonCeCx_Manus, inonCeCx_Mpenis, inonCeCx_MOPC);
}

	prop_MaleCx_MSM_num = 1 # one set of prop_MaleCx_MSM

	output_template_overall = array(0, dim=dim(Cost_array)); # a single overall value for each VC/Vdur/iPSA
	output_template_overall_bypMaleCx = rep(list(output_template_overall), prop_MaleCx_MSM_num);

	Costsave_nonCeCx_array_adjMaleCxFOV_bypMaleCx = output_template_overall_bypMaleCx
	QALYgain_nonCeCx_array_adjMaleCxFOV_bypMaleCx = output_template_overall_bypMaleCx
	LYgain_nonCeCxDeath_array_adjMaleCxFOV_bypMaleCx = output_template_overall_bypMaleCx

	count_nonCeCx_array_adjMaleCxFOV_bypMaleCx = output_template_overall_bypMaleCx
	Caseprev_nonCeCx_array_adjMaleCxFOV_bypMaleCx = output_template_overall_bypMaleCx

	# byType, for Costsave_nonCeCx_array_byType_forSummary
	Costsave_nonCeCx_array_byType_bypMaleCx = rep(list(output_template_overall_bypMaleCx), nonCeCx_num)
	QALYgain_nonCeCx_array_byType_bypMaleCx = Costsave_nonCeCx_array_byType_bypMaleCx
	LYgain_nonCeCxDeath_array_byType_bypMaleCx = Costsave_nonCeCx_array_byType_bypMaleCx
	
	Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx = setNames(rep(list(rep(list(NULL), prop_MaleCx_MSM_num)), nonCeCx_num), nonCeCx_vec);
	QALYgain_MaleCx_array_adjMaleCxFOV_bypMaleCx = Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx
	LYgain_MaleCxDeath_array_adjMaleCxFOV_bypMaleCx = Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx



	for (iprop_MaleCx in 1:prop_MaleCx_MSM_num){
		
		
		# MaleCx / Cx for MSM 
		for (inonCeCx in 1:nonCeCx_num){
		
			if (! (inonCeCx %in% inonCeCx_adjMSM_vec)){
				# no need to change it's not related to MaleCx
				next
			} 
			# use the findings adjusted for prop_MaleCx_MSM

			Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]] = Costsave_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]]
			QALYgain_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]] = QALYgain_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]]
			LYgain_MaleCxDeath_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]] = LYgain_nonCeCxDeath_array_byType_noMaleCxFOV[[inonCeCx]]

		} # for- inonCeCx
		
		# all nonCeCx
		for (iVdur in 1:Vdur_num){
		for (iVC in 1:VC_boys_num){
			
			Costsave_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]][iVC, iVdur, ] = rowSums( do.call(cbind, sapply(1:nonCeCx_num, simplify=FALSE, function(inonCeCx)
				Costsave_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
				) ))
			QALYgain_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]][iVC, iVdur, ] = rowSums( do.call(cbind, sapply(1:nonCeCx_num, simplify=FALSE, function(inonCeCx)
				QALYgain_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
				) ))
			LYgain_nonCeCxDeath_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]][iVC, iVdur, ] = rowSums( do.call(cbind, sapply(1:nonCeCx_num, simplify=FALSE, function(inonCeCx)
				LYgain_nonCeCxDeath_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
				) ))

				
			count_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]][iVC, iVdur, ] = rowSums( do.call(cbind, sapply(1:nonCeCx_num, simplify=FALSE, function(inonCeCx)
				count_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
				) ))
			Caseprev_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]][iVC, iVdur, ] = rowSums( do.call(cbind, sapply(1:nonCeCx_num, simplify=FALSE, function(inonCeCx)
				Caseprev_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
				) ))


		} # for-iVC
		} # for-iVdur
		# finish calculating the adjustment for each cost/QALY related item
		
		for (inonCeCx in 1:nonCeCx_num){
			if (!(inonCeCx %in% inonCeCx_adjMSM_vec)){ next}
		
			for (iVdur in 1:Vdur_num){
			for (iVC in 1:VC_boys_num){
			
			Costsave_nonCeCx_array_byType_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][iVC, iVdur, ] = Costsave_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
			QALYgain_nonCeCx_array_byType_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][iVC, iVdur, ] = 
				QALYgain_nonCeCx_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]
			LYgain_nonCeCxDeath_array_byType_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][iVC, iVdur, ] = 
				LYgain_nonCeCxDeath_array_byType_noMaleCxFOV[[inonCeCx]][iVC, iVdur, ]

			} # for-iVC
			} # for-iVdur
			
		} # for- inonCeCx

	} # for-iprop_MaleCx


# save all ICERs and TVCs by calc_gwarts_YN_vec(2), numDose_Vacc_new_vec (2), pMaleCx_FOV, Vdur
ICERall_out = array(list(NULL), dim=c(2, 2, 1+prop_MaleCx_MSM_num, Vdur_num));
dimnames(ICERall_out) = list(paste0("gwarts_",0:1), paste0("numDose_Vacc_",1:2), c("pMaleCx_FOV_noChange", "pMaleCx_FOV_adjust"), paste0("Vdur_",Vdur_vec) )
TVCall_out = ICERall_out;


for (calc_gwarts_YN in calc_gwarts_YN_vec){
ii_gwarts_ICER_all_out = calc_gwarts_YN+1; # ICERall_out/TVCall_out

for (numDose_Vacc_new in numDose_Vacc_new_vec){

for (run_SensAnaly_propMaleCx_YN in run_SensAnaly_propMaleCx_YN_vec){

# recall
prop_MaleCx_MSM = c("F_vagina_vulva" = NA, "F_OPC" = NA, "F_anus" = NA, "M_OPC" = prop_MOPC_MSM, "M_penis" = prop_Mpenis_MSM, "M_anus" = prop_Manus_MSM);



for (iprop_MaleCx in 1:prop_MaleCx_MSM_num){

if (run_SensAnaly_propMaleCx_YN){
	writeLines( sprintf("\nprop_MaleCx_MSM: %s", paste(sapply(which(!is.na(prop_MaleCx_MSM)), simplify=TRUE, function(inonCeCx) 
sprintf("%s %.2f", names(prop_MaleCx_MSM)[inonCeCx], prop_MaleCx_MSM[inonCeCx])), collapse="; ")

) )
}

if (run_SensAnaly_propMaleCx_YN==FALSE){ # ICERall_out/TVCall_out
	ii_pMaleCx_ICERall_out = 1;
} else if (run_SensAnaly_propMaleCx_YN==TRUE){
	ii_pMaleCx_ICERall_out = 1 + iprop_MaleCx; 
}


Vdur_vec_comp = 1:Vdur_num # vector of Vdur to compare; numDose_Vacc_new may change Vdur_vec_comp

cat("\n\n"); print( sprintf("numDose_Vacc_new = %d", numDose_Vacc_new) )

if (calc_gwarts_YN==1){
	print( "_sensitivity analysis_ include genital warts" )
}

if (numDose_Vacc_new==2){
	diff_calc_1dose_YN = FALSE;
} else if (numDose_Vacc_new==1){
	diff_calc_1dose_YN = TRUE;
}


## compare cost/QALY by different outcomes (CeCx, nonCeCx, gwarts)
# for CeCx, higher VC, lower HPV infection, QALY_normal increases, QALY_cecx decreases, QALY_screening/treatCIN uncertain
diff_Cost_VaccGirls_byVdur = rep(list(NULL), Vdur_num);
diff_Cost_VaccBoys_byVdur = rep(list(NULL), Vdur_num);

diff_Cost_CeCx_vsNoVacc_byVdur = rep(list(NULL), Vdur_num);
diff_QALYgain_CeCx_vsNoVacc_byVdur = diff_Cost_CeCx_vsNoVacc_byVdur;

diff_Cost_CeCx_treatCxOnly_vsNoVacc_byVdur = rep(list(NULL), Vdur_num);
diff_QALYgain_CeCx_treatCxOnly_vsNoVacc_byVdur = diff_Cost_CeCx_vsNoVacc_byVdur;

diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur = rep(list(rep(list(NULL), nonCeCx_num)), Vdur_num);
diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur= diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur;
diff_LYgain_nonCeCxDeath_vsNoVacc_byType_byVdur= diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur;
diff_QALYgain_plusLYgainDeath_nonCeCx_vsNoVacc_byType_byVdur= diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur;

diff_Costsave_gwarts_vsNoVacc_bySex_byVdur = rep(list(rep(list(NULL), sex_num)), Vdur_num);
diff_QALYgain_gwarts_vsNoVacc_bySex_byVdur = diff_Costsave_gwarts_vsNoVacc_bySex_byVdur


# WTPthreshold = Cost_diff / QALY_diff = (Cost_diff_woVacc + Cost_diff_onlyVacc) / QALY_diff
# Cost_diff_onlyVacc = (WTPthrehold * QALY_diff - Cost_diff_woVacc)
# TVCratio = (WTPthrehold * QALY_diff - Cost_diff_woVacc) / Cost_diff_onlyVacc
TVC_array = array(0, dim=dim(Cost_array)) # list(NULL)



	# direct comparison between GNV of different uptakes for boys vs FOV of 0% uptake of boys

	Cost_Vacc_use_new = Cost_Vacc_perdose_new + Cost_GPconsult + consider_Cost_Transport_YN*Cost_Transport;
	Cost_Vacc_adjfactor = (numDose_Vacc_new/numDose_Vacc)*(Cost_Vacc_use_new/Cost_Vacc_use_base);

	print( sprintf("Cost_Vacc_perdose_new = %.1f", Cost_Vacc_perdose_new) )
	print( sprintf("Cost_Vacc_use_new = %.1f", Cost_Vacc_use_new) )
	print( sprintf("numDose_Vacc_new = %d", numDose_Vacc_new) )
	print( sprintf("Cost_Vacc_adjfactor = %.2f", Cost_Vacc_adjfactor) )

	# cost factor for 1dose for boys
	Cost_Vacc_adjfactor_1dose_boys = (1/numDose_Vacc) * (Cost_Vacc_use_new/Cost_Vacc_use_base); # should be same as Cost_Vacc_adjfactor
																
	print( sprintf("Cost_Vacc_adjfactor_1dose_boys = %.2f", Cost_Vacc_adjfactor_1dose_boys) )


	# need to apply the Cost_Vacc_adjfactor to Cost_vaccGirls too because Cost_vaccGirls was calculated based on (numDose_Vacc/numDose_Vacc_Vacc_Cpp)

	Cost_vaccGirls_array_TEMP = Cost_Vacc_adjfactor*Cost_vaccGirls_array
	Cost_vaccBoys_array_TEMP = Cost_Vacc_adjfactor_1dose_boys*Cost_vaccBoys_array																						 
	if (run_SensAnaly_propMaleCx_YN){
		Cost_overall_array = Cost_array + Cost_vaccGirls_array_TEMP + Cost_vaccBoys_array_TEMP - Costsave_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]]
		QALY_overall_array = QALY_array + QALYgain_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]] + LYgain_nonCeCxDeath_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]]
		
		Cost_overall_array_woVacc = Cost_array - Costsave_nonCeCx_array_adjMaleCxFOV_bypMaleCx[[iprop_MaleCx]] # overall cost without vaccination cost
	} else{
		Cost_overall_array = Cost_array + Cost_vaccGirls_array_TEMP + Cost_vaccBoys_array_TEMP - Costsave_nonCeCx_array
		QALY_overall_array = QALY_array + QALYgain_nonCeCx_array + LYgain_nonCeCxDeath_array
		
		Cost_overall_array_woVacc = Cost_array - Costsave_nonCeCx_array # overall cost without vaccination cost
	}
	
	Cost_overall_array_onlyVacc_girlsboys_noadjust = Cost_vaccGirls_array + Cost_vaccBoys_array # overall cost for vaccination cost only, no adjustment

	# genital warts
	if (calc_gwarts_YN==1){
		Cost_overall_array = Cost_overall_array - Reduce("+", lapply(Costsave_gwarts_array_bysex, function(xlist) gwarts_cost_ratio*xlist))
		QALY_overall_array = QALY_overall_array + Reduce("+", lapply(QALYgain_gwarts_array_bysex, function(xlist) gwarts_utilityloss_ratio*xlist))

		Cost_overall_array_woVacc = Cost_overall_array_woVacc - Reduce("+", lapply(Costsave_gwarts_array_bysex, function(xlist) xlist))
	}
	
	# diff_Costsave/QALYgain_vsNoVacc; benefit of VCboys
	VC_boys_comp_T = 1:VC_boys_num; # keep all VC to maintain the indexing for iVC
	for (iVdur in 1:Vdur_num){
		# CeCx, _Cost_ is the absolute cost, so the difference (benefit of VCboys) is VCboys_0 - VCboys_others; _QALY_ is the absolute QALY, the diff (benefit of VCboys) is VCboys_others - VCboys_0
		diff_Cost_CeCx_vsNoVacc_byVdur[[iVdur]] = pracma::repmat(Cost_array[1, iVdur, ], n=VC_boys_num,m=1) - Cost_array[VC_boys_comp_T, iVdur, ]
		diff_QALYgain_CeCx_vsNoVacc_byVdur[[iVdur]] = QALY_array[VC_boys_comp_T, iVdur, ] - pracma::repmat(QALY_array[1, iVdur, ], n=VC_boys_num,m=1)
		# CeCx, related to cancer treatment only; higher VCboys has better benefit with fewer CeCx cases, so QALYgain = VCboys_0 - VCboys_others
		diff_Cost_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]] = pracma::repmat(Cost_CeCx_treatCancer_array[1, iVdur, ], n=VC_boys_num,m=1) - Cost_CeCx_treatCancer_array[VC_boys_comp_T, iVdur, ]
		diff_QALYgain_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]] = pracma::repmat(QALY_CeCx_treatCancer_array[1, iVdur, ], n=VC_boys_num,m=1) - QALY_CeCx_treatCancer_array[VC_boys_comp_T, iVdur, ]
		
		# nonCeCx, gwarts, _Costsave_ and _QALYgain_ are recorded, so the benefit of VCboys is VCboys_others - VCboys_0
		for (inonCeCx in 1:nonCeCx_num){
			if (inonCeCx %in% inonCeCx_adjMSM_vec & run_SensAnaly_propMaleCx_YN==TRUE){
				diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(Costsave_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][1, iVdur, ], n=VC_boys_num,m=1)
				diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		QALYgain_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(QALYgain_MaleCx_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][1, iVdur, ], n=VC_boys_num,m=1)
				diff_LYgain_nonCeCxDeath_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		LYgain_MaleCxDeath_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(LYgain_MaleCxDeath_array_adjMaleCxFOV_bypMaleCx[[inonCeCx]][[iprop_MaleCx]][1, iVdur, ], n=VC_boys_num,m=1)
			} else{
				diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		Costsave_nonCeCx_array_byType[[inonCeCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(Costsave_nonCeCx_array_byType[[inonCeCx]][1, iVdur, ], n=VC_boys_num,m=1)
				diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		QALYgain_nonCeCx_array_byType[[inonCeCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(QALYgain_nonCeCx_array_byType[[inonCeCx]][1, iVdur, ], n=VC_boys_num,m=1)
				diff_LYgain_nonCeCxDeath_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = 		LYgain_nonCeCxDeath_array_byType[[inonCeCx]][VC_boys_comp_T, iVdur, ] - pracma::repmat(LYgain_nonCeCxDeath_array_byType[[inonCeCx]][1, iVdur, ], n=VC_boys_num,m=1)
		}
		diff_QALYgain_plusLYgainDeath_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] = diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]] + diff_LYgain_nonCeCxDeath_vsNoVacc_byType_byVdur[[iVdur]][[inonCeCx]];
		} # end inonCeCx
		if (calc_gwarts_YN==1){
			for (isex in 1:sex_num){
				diff_Costsave_gwarts_vsNoVacc_bySex_byVdur[[iVdur]][[isex]] = 		Costsave_gwarts_array_bysex[[isex]][VC_boys_comp_T, iVdur, ] - pracma::repmat(Costsave_gwarts_array_bysex[[isex]][1, iVdur, ], n=VC_boys_num,m=1)
				diff_QALYgain_gwarts_vsNoVacc_bySex_byVdur[[iVdur]][[isex]] = 		QALYgain_gwarts_array_bysex[[isex]][VC_boys_comp_T, iVdur, ] - pracma::repmat(QALYgain_gwarts_array_bysex[[isex]][1, iVdur, ], n=VC_boys_num,m=1)
			} # end isex
		} # if-calc_gwarts_YN==1
	} # for-iVdur
	

# summary() work on matrix/data.frame. so apply summary_new() to vector only
quant_probs = c(0.5, 0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 1)
summary_new = function(x) { if (TRUE) as.table(c("Mean"=mean(x), quantile(x, quant_probs))) else summary(x) }


# clinical outcomes of genital warts, gwarts prevented; no need to repeat for SensAnalys_pMaleCx
if (run_SensAnaly_propMaleCx_YN==FALSE){
if (calc_gwarts_YN==1){
	for (iVdur in Vdur_vec_comp){
		writeLines("")
		print( "gwarts prevented, compared to irow_ref_noVacc" )
		print( sprintf("iVdur = %d; xVdur = %d", iVdur, Vdur_vec[iVdur]) )
		if (diff_calc_1dose_YN==TRUE & numDose_Vacc_new==1){
			print( "diff_calc_1dose_YN==TRUE & numDose_Vacc_new==1" )
			print( "1d-GNV is compared to 2d-FOV; so the case for VC_boys==0 refers to 2d-FOV" )
		}
		for (isex in 1:sex_num){
			print( names(sex_full_vec)[isex] )
			
			print( sprintf("at age gwarts_age_CostQALY: %d", gwarts_age_CostQALY) )
			print( "gwarts prevented, compared to irow_ref_noVacc" )
			Caseprev_gwarts_summ_Temp = t(apply(Caseprev_gwarts_array_bysex[[isex]][, iVdur, ],1, summary_new))
			rownames(Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec);
			print( Caseprev_gwarts_summ_Temp );
			print( "relative gwarts prevented, compared to irow_ref_noVacc" )
			relative_Caseprev_gwarts_summ_Temp = t(apply(relative_Caseprev_gwarts_array_bysex[[isex]][, iVdur, ],1, summary_new))
			rownames(relative_Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec);
			print( relative_Caseprev_gwarts_summ_Temp );
			print( "additional relative gwarts prevented, compared to VCboys=0" )
			relative_Caseprev_gwarts_summ_Temp = do.call(rbind, sapply(2:VC_boys_num, simplify=FALSE, function(iiVCboys)
				summary_new(relative_Caseprev_gwarts_array_bysex[[isex]][iiVCboys, iVdur, ] - relative_Caseprev_gwarts_array_bysex[[isex]][1, iVdur, ])
			))
			rownames(relative_Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec[2:VC_boys_num]);
			print( relative_Caseprev_gwarts_summ_Temp );
			print( t(apply(relative_Caseprev_gwarts_summ_Temp, 1, function(x) round(x*100, digit=1))) ) # not to use print_summstat because the input ("_summ_Temp") is already the summ_stat
			
			for (iage_Caseprev_gw in 1:gwarts_Caseprev_age_vec_num){
				print( sprintf("at age gwarts_age_Caseprev: %d", gwarts_Caseprev_age_vec[iage_Caseprev_gw]) )
				print( "gwarts prevented, compared to irow_ref_noVacc" )
				Caseprev_gwarts_summ_Temp = t(apply(Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][, iVdur, ],1, summary_new))
				rownames(Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec);
				print( Caseprev_gwarts_summ_Temp );
				print( "relative gwarts prevented, compared to irow_ref_noVacc" )
				relative_Caseprev_gwarts_summ_Temp = t(apply(relative_Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][, iVdur, ],1, summary_new))
				rownames(relative_Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec);
				print( relative_Caseprev_gwarts_summ_Temp );
				print( "additional relative gwarts prevented, compared to VCboys=0" )

				relative_Caseprev_gwarts_summ_Temp = do.call(rbind, sapply(2:VC_boys_num, simplify=FALSE, function(iiVCboys)
					summary_new(relative_Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][iiVCboys, iVdur, ] - relative_Caseprev_gwarts_array_bysex_age_list[[iage_Caseprev_gw]][[isex]][1, iVdur, ])
				))
				rownames(relative_Caseprev_gwarts_summ_Temp) = paste0("VCboys_", VC_boys_vec[2:VC_boys_num]);
				print( relative_Caseprev_gwarts_summ_Temp );
				print( t(apply(relative_Caseprev_gwarts_summ_Temp, 1, function(x) round(x*100, digit=1))) )
			} # for- iage_Caseprev_gw
			
		} # for- isex
	} # for- iVdur, gwarts prevented
} # if- calc_gwarts_YN 
} # if- run_SensAnaly_propMaleCx_YN==FALSE


# calculate ICER
WTPthreshold = 406538 # be aware of unit/currency
plot_Cost_USD_YN = TRUE;


ICER_byVdur = rep(list(NULL), Vdur_num)
names(ICER_byVdur) = Vdur_vec
ICERsumm_byVdur = ICER_byVdur
for (iVdur in Vdur_vec_comp){ # Vdur_vec_comp would change for 1-dose GNV vs 2-dose FOV
	writeLines( sprintf('\niVdur %d', iVdur) )
	Cost_overall_array_iVdur = Cost_overall_array[,iVdur,]
	QALY_overall_array_iVdur = QALY_overall_array[,iVdur,]
	
	names_stg = sapply(VC_boys_vec[-1], function(x) sprintf("%dvs%s", x, VC_boys_vec[1]))
	
	
	# use repmat to retain the dimension; expand rows only, so change n
	Cost_diff_vsNoVacc = Cost_overall_array_iVdur - pracma::repmat(Cost_overall_array_iVdur[1,], m=1, n=nrow(Cost_overall_array_iVdur))
	QALY_diff_vsNoVacc = QALY_overall_array_iVdur - pracma::repmat(QALY_overall_array_iVdur[1,], m=1, n=nrow(Cost_overall_array_iVdur))
	if (numDose_Vacc_new==2){ # non-negative QALY changes for 2dose GNV (incl 2F1M) vs 2dose FOV
		QALY_diff_vsNoVacc = pmax(QALY_diff_vsNoVacc, 0)
	}	
	rownames(Cost_diff_vsNoVacc) = sprintf("VCboys_%dvs0", VC_boys_vec)
	rownames(QALY_diff_vsNoVacc) = rownames(Cost_diff_vsNoVacc)
	NMB_vsNoVacc = WTPthreshold*QALY_diff_vsNoVacc - Cost_diff_vsNoVacc	# net monetary benefit
	

	print( "summary stat of Cost_overall_array_iVdur" )
	print( apply(Cost_overall_array_iVdur,1, summary) )
	print( "summary stat of QALY_overall_array_iVdur" )
	print( apply(QALY_overall_array_iVdur,1, summary) )

	print( "summary stat of Cost_diff_vsNoVacc" )
	print( apply(Cost_diff_vsNoVacc,1, summary) )
	print( "summary stat of QALY_diff_vsNoVacc" )
	print( apply(QALY_diff_vsNoVacc,1, summary) )	
	print( "summary stat of NMB_vsNoVacc" )
	print( apply(NMB_vsNoVacc,1, summary) )	
	
	ICER_vsNoVacc = t(sapply(2:VC_boys_num, function(x) Cost_diff_vsNoVacc[x,]/QALY_diff_vsNoVacc[x,]))
	rownames(ICER_vsNoVacc) = names_stg
	ICER_byVdur[[iVdur]] = ICER_vsNoVacc

	ICERsumm_temp = t( apply(ICER_vsNoVacc, 1, summary_new) )
	rownames(ICERsumm_temp) = names_stg
	ICERsumm_byVdur[[iVdur]] = ICERsumm_temp
	
	# export ICERs of all simulations
	# browser()
	ICERall_out[[ii_gwarts_ICER_all_out, numDose_Vacc_new, ii_pMaleCx_ICERall_out, iVdur]] = ICER_vsNoVacc;
	
	
	# _YN for ICER < WTP;  Cost_diff>0 (use more money); QALY_diff>0 (better health)
	ICER_vsNoVacc_lessWTP_YN = (ICER_vsNoVacc < WTPthreshold)
	Cost_diff_vsNoVacc_positive_YN = (Cost_diff_vsNoVacc[2:VC_boys_num,, drop=FALSE] > 0)
	QALY_diff_vsNoVacc_positive_YN = (QALY_diff_vsNoVacc[2:VC_boys_num,, drop=FALSE] > 0)
	NMB_vsNoVacc_positive_YN = (NMB_vsNoVacc[2:VC_boys_num,, drop=FALSE] > 0)

	Quadrant1234 = 1*(QALY_diff_vsNoVacc_positive_YN & Cost_diff_vsNoVacc_positive_YN) + 2*(!QALY_diff_vsNoVacc_positive_YN & Cost_diff_vsNoVacc_positive_YN) + 3*(!QALY_diff_vsNoVacc_positive_YN & !Cost_diff_vsNoVacc_positive_YN) + 4*(QALY_diff_vsNoVacc_positive_YN & !Cost_diff_vsNoVacc_positive_YN)
	print( "table(Quadrant1234)" )
	apply(Quadrant1234, 1, table, simplify=FALSE)


	print( "prop%(QALY_diff_vsNoVacc > 0)" )
	print( rowMeans(QALY_diff_vsNoVacc_positive_YN) )
	print( "prop%(Cost_diff_vsNoVacc > 0)" )
	print( rowMeans(Cost_diff_vsNoVacc_positive_YN) )
	print( "prop%(ICER_vsNoVacc < WTPthreshold & QALY_diff_vsNoVacc >0)" )
	print( rowMeans(ICER_vsNoVacc_lessWTP_YN & QALY_diff_vsNoVacc_positive_YN) )
	print( "prop%(ICER_vsNoVacc < WTPthreshold); suitable for comparisons if all points are in Q1" )
	print( rowMeans(ICER_vsNoVacc_lessWTP_YN) )
	print( "prop%(NMB_vsNoVacc >0); applicable to with or without points in Q3" )
	print( rowMeans(NMB_vsNoVacc_positive_YN) )


	print( "prop%(QALY_diff_vsNoVacc > 0 & Cost_diff_vsNoVacc > 0); Q1, the normal quadrant" )
	print( rowMeans(Quadrant1234==1) )
	print( "prop%(QALY_diff_vsNoVacc < 0 & Cost_diff_vsNoVacc > 0); Q2, the strongly dominated quadrant (unfavored)" )
	print( rowMeans(Quadrant1234==2) )
	print( "prop%(QALY_diff_vsNoVacc < 0 & Cost_diff_vsNoVacc < 0); Q3, the uncertain quadrant" )
	print( rowMeans(Quadrant1234==3) )
	print( "prop%(QALY_diff_vsNoVacc > 0 & Cost_diff_vsNoVacc < 0); Q4, the cost-saving quadrant" )
	print( rowMeans(Quadrant1234==4) )

	print( "n and prop%(ICER_vsNoVacc < WTPthreshold & Q1); cost-effective" )
		print( rowSums(ICER_vsNoVacc_lessWTP_YN & Quadrant1234==1) )
		print( rowMeans(ICER_vsNoVacc_lessWTP_YN & Quadrant1234==1) )
		cat( "prop% among points in Q1 only\n" )
		print( rowSums(ICER_vsNoVacc_lessWTP_YN & Quadrant1234==1)/rowSums(Quadrant1234==1) )
	print( "n and prop%(NMB_vsNoVacc >0 & Q3); preferred for points in Q3" )
		print( rowSums(NMB_vsNoVacc_positive_YN & Quadrant1234==3) )
		print( rowMeans(NMB_vsNoVacc_positive_YN & Quadrant1234==3) )
		cat( "prop% among points in Q3 only\n" )
		print( rowSums(NMB_vsNoVacc_positive_YN & Quadrant1234==3)/rowSums(Quadrant1234==3) )	
	

	cat("\n")
	print( "ICER summary" )
	print( ICERsumm_temp )
	print( "ICER summary (USD)" )
	print( ICERsumm_temp/7.8 )

	
	## scatter plot or hist2d
	plot_hist2d_YN = 1; # 0 for a scatter plot
	# setting for hist2d
	hist2d_set = within(list(), {
	  breaks_color = c(5, seq(50, 500, length.out=10-1), 10000)
	  color = rev(hcl.colors(length(breaks_color)-1, "Blues 3")); #rev(c(heat.colors(9))); # transparent
	  nbins = 20;
	  })

	# setting for scatter
	pch_plot = 16;

	# point for the mean point of cost and QALY
	col_meanpt = "orange"; # ifelse(plot_hist2d_YN==1, "blue", "orange");
	pch_meanpt = 17; # pch=4 for "X", =17 for a "triangle"
	cex_meanpt = 2.25; # =3 for a "X", =2.25 for a triangle 
	lwd_meanpt = 1; # =3 for a "X", =1 for a triangle 

	
	Cost_plotunit = 1/1000000
	lab_Cost_plotunit = " (per 1M)"
	if (plot_Cost_USD_YN){ # USD
		Cost_plotunit = Cost_plotunit/7.8
		lab_Cost_plotunit = " (per 1M USD)"	
	}	
	Cost_diff_vsNoVacc_plot = Cost_diff_vsNoVacc * Cost_plotunit
	ylab_plot = paste0("Difference in discounted cost", lab_Cost_plotunit)
	xlab_plot = "Difference in discounted QALY" #for GNV vs FOV"

	xlim_plot = range(QALY_diff_vsNoVacc[2:VC_boys_num,])
	xlim_plot = range(pretty(xlim_plot, n=2))
	ylim_plot = range(Cost_diff_vsNoVacc_plot[2:VC_boys_num,])
	ylim_plot = range(pretty(ylim_plot, n=2))
	
	
	if (numDose_Vacc_new_Boys==1){
		# ICER among simulations with positive QALY, exclude negative QALY
		ICER_vsNoVacc_posQALY = sapply(1:nrow(ICER_vsNoVacc), simplify=FALSE, function(x) ICER_vsNoVacc[x, QALY_diff_vsNoVacc_positive_YN[x,]])
		names(ICER_vsNoVacc_posQALY) = rownames(ICER_vsNoVacc)
		print( "stat for ICER_vsNoVacc _among_ QALY_diff_vsNoVacc>0" )
		print( "prop%(ICER_vsNoVacc < WTPthreshold _among_ QALY_diff_vsNoVacc >0)" )
		print( sapply(ICER_vsNoVacc_posQALY, function(x) mean(x<WTPthreshold)) )
		print( "summary for ICER_vsNoVacc _among_ QALY_diff_vsNoVacc>0" )
		print( do.call(rbind, sapply(ICER_vsNoVacc_posQALY, simplify=FALSE, function(x) c(n=length(x), summary_new(x)))) )
		print( "summary for ICER_vsNoVacc _among_ QALY_diff_vsNoVacc>0 (USD)" )
		print( do.call(rbind, sapply(ICER_vsNoVacc_posQALY, simplify=FALSE, function(x) c(n=length(x), summary_new(x)/7.8))) )


		if (subplot_iVC_YN){ # scatter plot
			# subplots
			if (plot_PDF_YN==FALSE){windows(width=num_subplot*5, height=1*5)}
			par(mar=c(4,5.75,3,0.75))
			layout(matrix(1:num_subplot, nrow=1))
			pch_col_alpha = adjustcolor('blue', alpha=0.25)

			for (iVC in 2:VC_boys_num){
				xplot = QALY_diff_vsNoVacc[iVC,]
				yplot = Cost_diff_vsNoVacc_plot[iVC,]
				if (plot_hist2d_YN==0){
					plot(x=xplot, y=yplot, xlim=xlim_plot, ylim=ylim_plot, xlab="", ylab="", bty="L", las=1, cex=0.9, col=pch_col_alpha, pch=pch_plot, axes=FALSE)
					abline(v=0, h=0, col='grey', lty=2) # reference lines
					abline(a=0, b=WTPthreshold*Cost_plotunit, col='green')
				} else if (plot_hist2d_YN==1){
					plot(NA, xlim=xlim_plot, ylim=ylim_plot, xlab="", ylab="", axes=FALSE)
					abline(v=0, h=0, col='grey', lty=2) # reference lines
					nullout <- with(hist2d_set, gplots::hist2d(x=xplot, y=yplot, add=TRUE,
				                    breaks=breaks_color, col=color, nbins = nbins));
					abline(a=0, b=WTPthreshold*Cost_plotunit, col='green')
				}
				axis(1, at=par()$usr[1:2], label=FALSE, col.ticks=NA)
				axis(1, col=NA, col.ticks='black', padj=-0.15, cex.axis=1.5)
				axis(2, at=par()$usr[3:4], label=FALSE, col.ticks=NA)
				axis(2, col=NA, col.ticks='black', las=1, hadj=0.95, cex.axis=1.5)
				mtext(side=1, xlab_plot, cex=1, line=2.5)
				mtext(side=2, ylab_plot, cex=1, line=3.75)

				mtext(text=sprintf("GNV %d%% uptake for boys", VC_boys_vec[iVC]), side=3, line=0.15) # , adj=0.05 
				mtext(text="2F1M-GNV vs 2-dose FOV", side=3, line=1.5) # , adj=0.05 
				points(x=mean(xplot), y=mean(yplot), pch=pch_meanpt, col=col_meanpt, cex=cex_meanpt, lwd=lwd_meanpt)
				# confidence ellipse
				# car::dataEllipse(xplot, yplot, levels=0.9, lwd=1, col='blue', lty=2, plot.points=FALSE)

				# % prop of CE
				leg_text = sprintf("CE: %.1f%%", 100*mean(ICER_vsNoVacc_posQALY[[iVC-1]]<WTPthreshold))
				with(par(), text(x=xaxp[2], y=0+diff(yaxp[1:2])/10, labels=leg_text, adj=c(1,0), cex=1.45))

				
				mtext_at = par()$usr[1]
				if (iVC==2){ # add text when plotting the first plot
					mtext(text=sprintf('iVdur=%d',iVdur), side=3, line=1.75, cex=0.7, adj=1, at=mtext_at)
					if (TRUE){ # calc_gwarts_YN){
						mtext(text=sprintf('gwarts=%s',ifelse(calc_gwarts_YN,"Y","N")), side=3, line=1, cex=0.7, adj=1, at=mtext_at)
					}
					if (run_SensAnaly_propMaleCx_YN){
						mtext(text="pMaleCx_adj", side=3, line=0.25, adj=-0.2, cex=0.7)
					}
				}
			} # for-iVC
		} else {

			plot(x=NA, y=NA, xlim=xlim_plot, ylim=ylim_plot, xlab="", ylab="", bty="L", las=1, axes=FALSE)
				abline(v=0, h=0, col='grey', lty=2)
				abline(a=0, b=WTPthreshold*Cost_plotunit, col='green')
			axis(1, at=par()$usr[1:2], label=FALSE, col.ticks=NA)
			axis(1, col=NA, col.ticks='black', padj=-0.15, cex.axis=1)
			axis(2, at=par()$usr[3:4], label=FALSE, col.ticks=NA)
			axis(2, col=NA, col.ticks='black', las=1, hadj=0.95, cex.axis=1)
			mtext(side=1, xlab_plot, cex=1, line=2.5)
			mtext(side=2, ylab_plot, cex=1, line=3.5)
			mtext(text="2F1M-GNV vs 2-dose FOV", side=3, line=1.5) # , adj=0.05 
			mtext(text=sprintf('Basecase\niVdur=%d',iVdur), side=3, line=1.5, adj=-0.10, cex=0.8)

			pch_col = c('blue', 'darkgreen', 'purple')
			pch_col_alpha = adjustcolor(pch_col, alpha=0.25)

			for (iVC in 2:VC_boys_num){
				xplot = QALY_diff_vsNoVacc[iVC,]
				yplot = Cost_diff_vsNoVacc_plot[iVC,]
				points(xplot, yplot, col=pch_col_alpha[iVC-1], pch=pch_plot, cex=0.9)
				points(x=mean(xplot), y=mean(yplot), pch=pch_meanpt, col=col_meanpt, cex=cex_meanpt, lwd=lwd_meanpt)
			} # for-iVC
			
			legend("topleft", legend=paste0(VC_boys_vec[2:VC_boys_num], "%"), fill=adjustcolor(pch_col, alpha=0.75), title="Vaccine uptake\nin boys", bty="n", title.cex=0.85)
			
		} # scatter plot
		
	} # if-numDose_Vacc_new==1
	
	
	
	# Threshold vaccination cost
	Cost_overall_array_woVacc_iVdur = Cost_overall_array_woVacc[,iVdur,]
	Cost_overall_array_onlyVacc_iVdur = Cost_overall_array_onlyVacc_girlsboys_noadjust[,iVdur,]
	TVCratio_T = WTPthreshold * QALY_diff_vsNoVacc[-1,, drop=FALSE] # exclude VCboys=0
	
	Cost_diff_vsNoVacc_woVacc_iVdur = Cost_overall_array_woVacc_iVdur - pracma::repmat(Cost_overall_array_woVacc_iVdur[1,], m=1, n=nrow(Cost_overall_array_woVacc_iVdur))
	Cost_diff_vsNoVacc_onlyVacc_iVdur = Cost_overall_array_onlyVacc_iVdur - pracma::repmat(Cost_overall_array_onlyVacc_iVdur[1,], m=1, n=nrow(Cost_overall_array_onlyVacc_iVdur))
	Cost_diff_vsNoVacc_woVacc_iVdur = Cost_diff_vsNoVacc_woVacc_iVdur[-1,, drop=FALSE] # exclude VCboys=0
	Cost_diff_vsNoVacc_onlyVacc_iVdur = Cost_diff_vsNoVacc_onlyVacc_iVdur[-1,, drop=FALSE] # exclude VCboys=0
	
	TVCratio_Temp = (TVCratio_T - Cost_diff_vsNoVacc_woVacc_iVdur)/Cost_diff_vsNoVacc_onlyVacc_iVdur
	TVC_inclAdminfee_Temp = TVCratio_Temp * Cost_Vacc_use_base # ** no need to be divided by "/ (numDose_Vacc_new/numDose_Vacc)" because numDose_Vacc has been accounted in Cost_overall_array_onlyVacc_girlsboys_noadjust 
	# TVC including (Cost_GPconsult + consider_Cost_Transport_YN*Cost_Transport)
	TVC_perdose_Temp = TVC_inclAdminfee_Temp - (Cost_GPconsult + consider_Cost_Transport_YN*Cost_Transport) # TVC per dose
	TVC_array[-1,iVdur,] = TVC_inclAdminfee_Temp
	
	print( sprintf("summary stat of TVC_inclAdminfee_array at iVdur=%s", iVdur) )
	print( t(apply(TVC_inclAdminfee_Temp,1, summary_new)) )
	print( sprintf("summary stat of TVC_perdose_array at iVdur=%s", iVdur) )
	print( t(apply(TVC_perdose_Temp,1, summary_new)) )
	print( sprintf("summary stat of TVC_inclAdminfee_array at iVdur=%s (USD)", iVdur) )
	print( t(apply(TVC_inclAdminfee_Temp,1, summary_new))/7.8 )
	print( sprintf("summary stat of TVC_perdose_array at iVdur=%s (USD)", iVdur) )
	print( t(apply(TVC_perdose_Temp,1, summary_new))/7.8 )
	print( "rowMeans((TVC_inclAdminfee_Temp > Cost_Vacc_use_new))" )
	print( rowMeans((TVC_inclAdminfee_Temp > Cost_Vacc_use_new)) ) # should match print( rowMeans(ICER_vsNoVacc_lessWTP_YN) ) [1-x] but does not match for cases of cost-saving

	# export TVCall_out of all simulations
	# browser()
	TVCall_out[[ii_gwarts_ICER_all_out, numDose_Vacc_new, ii_pMaleCx_ICERall_out, iVdur]] = TVC_inclAdminfee_Temp;

	
	
	# print summary stat
	summstat_Temp = cbind("ICER_USD"=apply(ICER_vsNoVacc/7.8, 1, print_summstat, digit=0),
		"TVC_inclAdminfee_USD"=apply(TVC_inclAdminfee_Temp/7.8, 1, print_summstat),
		"relChange_TVC"=apply(TVC_inclAdminfee_Temp/Cost_Vacc_perdose_inclAdmin_ref*100, 1, print_summstat),
		"relReduce_TVC"=apply((1-TVC_inclAdminfee_Temp/Cost_Vacc_perdose_inclAdmin_ref)*100, 1, print_summstat)
		)
	print( "summary stat" )
	print(summstat_Temp)	
	

	# CEAC
	par(mar=c(4,4.75,3,2))
	layout(matrix(1:num_subplot, nrow=1))
	WTP_vec_max = 100000
	WTP_vec = seq(0, WTP_vec_max, by=100);
	col_ceac_vec = c(NA, "red", "blue", "darkgreen")
	ylab_plot = "Probability of being adoptable" # cost-effective
	xlab_plot = paste0("WTP threshold", ifelse(plot_Cost_USD_YN, " (USD)", ""))

	# one plot for each iVC
	for (iVC in 2:VC_boys_num){
		ceac_cost = Cost_diff_vsNoVacc_plot[iVC,];
		ceac_health = QALY_diff_vsNoVacc[iVC,];
		ceac_nmb = do.call(cbind, sapply(WTP_vec, simplify=FALSE, function(xWTP) (xWTP*Cost_plotunit*ifelse(plot_Cost_USD_YN,7.8,1))*ceac_health - ceac_cost));
		ceac_probCE = colMeans(1*(ceac_nmb>0))
		
		plot(NA, ylim=c(0,1), xlim=c(0, WTP_vec_max), xlab="",ylab="", axes=FALSE, xaxs="i", yaxs="i" )
		segments(x0=WTPthreshold/ifelse(plot_Cost_USD_YN,7.8,1), y0=0,y1=1, col="grey80", lwd=1.5, lty=2)
		xtick_at = with(par(), seq(xaxp[1], xaxp[2], length.out=xaxp[3]+1))
		xtick_at_lab = format(xtick_at, big.mark=",", scientific=FALSE, trim=TRUE)
		axis(1, at=xtick_at, labels=xtick_at_lab, cex.axis=1.35)
		mtext(side=1, text=xlab_plot, line=2.5)
		ytick_at = seq(0,1, by=0.2)
		axis(2, at=ytick_at, cex.axis=1.35, las=1)
		mtext(side=2, text=ylab_plot, line=3)
		
		mtext(text=sprintf("GNV %d%% uptake for boys", VC_boys_vec[iVC]), side=3, line=0.5) # , adj=0.05 
		lines(x=WTP_vec, y=ceac_probCE, lwd=2, col=col_ceac_vec[iVC], xpd=TRUE)
	} # for-iVC


	# one plot for all iVC
	par(mar=c(4,4.75,1.5,1.5))
	
	scale_oneplotCEAC = 2;
	
	layout(matrix(1:scale_oneplotCEAC, nrow=1))
	
	plot(NA, ylim=c(0,1), xlim=c(0, WTP_vec_max), xlab="",ylab="", axes=FALSE, xaxs="i", yaxs="i" )
	segments(x0=WTPthreshold/ifelse(plot_Cost_USD_YN,7.8,1), y0=0,y1=1, col="grey80", lwd=1.5, lty=2)
	xtick_at = with(par(), seq(xaxp[1], xaxp[2], length.out=xaxp[3]+1))
	xtick_at_lab = format(xtick_at, big.mark=",", scientific=FALSE, trim=TRUE)
	axis(1, at=xtick_at, labels=xtick_at_lab, cex.axis=1.35)
	mtext(side=1, text=xlab_plot, line=2.5, cex=1.45)
	ytick_at = seq(0,1, by=0.2)
	axis(2, at=ytick_at, cex.axis=1.35, las=1)
									 
	mtext(side=2, text=ylab_plot, line=3, cex=1.45)
	
	
	for (iVC in 2:VC_boys_num){
		ceac_cost = Cost_diff_vsNoVacc_plot[iVC,];
		ceac_health = QALY_diff_vsNoVacc[iVC,];
		ceac_nmb = do.call(cbind, sapply(WTP_vec, simplify=FALSE, function(xWTP) (xWTP*Cost_plotunit*ifelse(plot_Cost_USD_YN,7.8,1))*ceac_health - ceac_cost));
		ceac_probCE = colMeans(1*(ceac_nmb>0))
		lines(x=WTP_vec, y=ceac_probCE, lwd=2, col=col_ceac_vec[iVC], xpd=TRUE)
	} # for-iVC

	legend(ifelse(numDose_Vacc_new==1, "topleft","bottomright"), legend=sprintf("%d%%", VC_boys_vec[2:VC_boys_num]), title="Uptake among boys", inset=c(ifelse(numDose_Vacc_new==1, 0.025,0),0.025), col=col_ceac_vec[2:VC_boys_num], lwd=2, cex=1.25)


	nullout = replicate(scale_oneplotCEAC-1, plot.new()) # create blank plot so that subsequent plots can align with the plot settings for each plot types
	
	
	# plot diff_Costsave/QALYgain_vsNoVacc by health outcome at iVdur

	# get stat_only, for pie plots, pick col_1 for the median/mean (the point estimate)
	# overall, includes cost_screening, QALY_normal, QALY_screening, LYgain_nonCeCxDeath
	stat_diff_Costsave_Temp_byVC = rep(list(NULL), VC_boys_num);
	stat_diff_QALYgain_Temp_byVC = stat_diff_Costsave_Temp_byVC;

	# related to treatment (cancer / gwarts) only; no LYgain_nonCeCxDeath
	stat_diff_Costsave_treatment_Temp_byVC = stat_diff_Costsave_Temp_byVC;
	stat_diff_QALYgain_treatment_Temp_byVC = stat_diff_Costsave_Temp_byVC;

	quantile_probs =c(0.5, 0.05,0.95)

	for (iVC in 2:VC_boys_num){
		
		# nonCeCx
		stat_diff_Costsave_Temp_byVC_nonCeCx = do.call(rbind, lapply(diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]], function(xmat) quantile(xmat[iVC,], probs=quantile_probs)))

		
		# QALYgain_plusLYgainDeath
		stat_diff_QALYgain_plusLYgainDeath_Temp_byVC_nonCeCx = do.call(rbind, lapply(diff_QALYgain_plusLYgainDeath_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]], function(xmat) quantile(xmat[iVC,], probs=quantile_probs)))
	
		# QALYgain related to treating Cx only
		stat_diff_QALYgain_treatmentOnly_Temp_byVC_nonCeCx = do.call(rbind, lapply(diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]], function(xmat) quantile(xmat[iVC,], probs=quantile_probs)))

		
		# CeCx related, appended by nonCeCx
		# overall
		stat_diff_Costsave_Temp_byVC[[iVC]] = rbind(quantile(diff_Cost_CeCx_vsNoVacc_byVdur[[iVdur]][iVC,], probs=quantile_probs), 
			stat_diff_Costsave_Temp_byVC_nonCeCx)
		stat_diff_QALYgain_Temp_byVC[[iVC]] = rbind(quantile(diff_QALYgain_CeCx_vsNoVacc_byVdur[[iVdur]][iVC,], probs=quantile_probs), 
			stat_diff_QALYgain_plusLYgainDeath_Temp_byVC_nonCeCx)


		# treatment only
		stat_diff_Costsave_treatment_Temp_byVC[[iVC]] = rbind(quantile(diff_Cost_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]][iVC,], probs=quantile_probs), 
			stat_diff_Costsave_Temp_byVC_nonCeCx)
		stat_diff_QALYgain_treatment_Temp_byVC[[iVC]] = rbind(quantile(diff_QALYgain_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]][iVC,], probs=quantile_probs), 
			stat_diff_QALYgain_treatmentOnly_Temp_byVC_nonCeCx)
			
		
		# gwarts
		if (calc_gwarts_YN==1){
			stat_diff_Costsave_Temp_byVC_gwarts = do.call(rbind, lapply(diff_Costsave_gwarts_vsNoVacc_bySex_byVdur[[iVdur]], function(xmat) quantile(xmat[iVC,], probs=quantile_probs)))
			stat_diff_QALYgain_Temp_byVC_gwarts = do.call(rbind, lapply(diff_QALYgain_gwarts_vsNoVacc_bySex_byVdur[[iVdur]], function(xmat) quantile(xmat[iVC,], probs=quantile_probs)))
		
		
			stat_diff_Costsave_Temp_byVC[[iVC]] = rbind(stat_diff_Costsave_Temp_byVC[[iVC]],
				stat_diff_Costsave_Temp_byVC_gwarts)
			stat_diff_QALYgain_Temp_byVC[[iVC]] = rbind(stat_diff_QALYgain_Temp_byVC[[iVC]],
				stat_diff_QALYgain_Temp_byVC_gwarts)

			stat_diff_Costsave_treatment_Temp_byVC[[iVC]] = rbind(stat_diff_Costsave_treatment_Temp_byVC[[iVC]],
				stat_diff_Costsave_Temp_byVC_gwarts)
			stat_diff_QALYgain_treatment_Temp_byVC[[iVC]] = rbind(stat_diff_QALYgain_treatment_Temp_byVC[[iVC]],
				stat_diff_QALYgain_Temp_byVC_gwarts)
		} # calc_gwarts_YN==1
	} # iVC


	# violin plot of the diff
	diff_Costsave_all_list = c(list(diff_Cost_CeCx_vsNoVacc_byVdur[[iVdur]]), diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]])
	diff_QALYgain_all_list = c(list(diff_QALYgain_CeCx_vsNoVacc_byVdur[[iVdur]]), diff_QALYgain_plusLYgainDeath_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]])

	diff_Costsave_Treatment_list = c(list(diff_Cost_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]]), diff_Costsave_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]])
	diff_QALYgain_Treatment_list = c(list(diff_QALYgain_CeCx_treatCxOnly_vsNoVacc_byVdur[[iVdur]]), diff_QALYgain_nonCeCx_vsNoVacc_byType_byVdur[[iVdur]])
	if (calc_gwarts_YN==1){
		diff_Costsave_all_list = c(diff_Costsave_all_list, diff_Costsave_gwarts_vsNoVacc_bySex_byVdur[[iVdur]])
		diff_QALYgain_all_list = c(diff_QALYgain_all_list, diff_QALYgain_gwarts_vsNoVacc_bySex_byVdur[[iVdur]])
		diff_Costsave_Treatment_list = c(diff_Costsave_Treatment_list, diff_Costsave_gwarts_vsNoVacc_bySex_byVdur[[iVdur]])
		diff_QALYgain_Treatment_list = c(diff_QALYgain_Treatment_list, diff_QALYgain_gwarts_vsNoVacc_bySex_byVdur[[iVdur]])
	}	


	# pie / bar chart
	col_pie_vec = RColorBrewer::brewer.pal(n=9, name = 'Set3')
	outlab_pie = c("F_cervix", nonCeCx_vec, "F_gwarts", "M_gwarts"); # outcome labels
	outlab_pie = gsub("F_vagina_vulva", "F_vv", outlab_pie)

	iout_plot = 1:7;
	if (calc_gwarts_YN==1){iout_plot = 1:9};

	xplot_lab_vec = c("Cost saved", "QALY gained", "Cost saved (treatment-related only)", "QALY gained (treatment-related only)")
	xplot_list = list(stat_diff_Costsave_Temp_byVC, stat_diff_QALYgain_Temp_byVC, stat_diff_Costsave_treatment_Temp_byVC, stat_diff_QALYgain_treatment_Temp_byVC)
	xviolin_list = list(diff_Costsave_all_list, diff_QALYgain_all_list, diff_Costsave_Treatment_list, diff_QALYgain_Treatment_list)


	# pie-chart
	for (iplot in 1:length(xplot_list)){
		layout(matrix(1:num_subplot, nrow=1))
		par(mar=c(4,1.5,3,1.5))
		xpie_plot = xplot_list[[iplot]] # pick col_1 only
		xpie_lab = xplot_lab_vec[iplot]
		for (iVC in 2:VC_boys_num){
			if (all(xpie_plot[[iVC]][,1]>0)){
				pie(xpie_plot[[iVC]][,1], col=col_pie_vec, labels=outlab_pie, cex=1.35, init.angle=90, clockwise=FALSE, radius=0.9)
			} else {
				# no pie chart for negative values
				plot(NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=FALSE)
				text(x=0,y=0, labels="At least one of the point estimates is negative.\nNo pie chart is available.", cex=1.35)
			}
			mtext(text=sprintf("GNV %d%% uptake for boys", VC_boys_vec[iVC]), side=3, line=-0.5) # , adj=0.05 
			if (iVC==2){
				mtext(text=xpie_lab, side=3, line=1.25, adj=0, at=par()$usr[1], cex=0.8)
			}
			
		} # iVC
		# one common legend at the bottom
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
		legend(x=0,y=-1.04, legend=outlab_pie[iout_plot], fill=col_pie_vec[iout_plot], horiz=TRUE, xjust=0.5, yjust=0, xpd=TRUE, cex=1.8)		
	} # iplot

	# bar chart
	par(mar=c(4,5.75,3,0.75))
	layout(matrix(1:num_subplot, nrow=1))
	# bar chart
	xbar_xlab_vec = c("Difference in discounted cost (1M USD)", "Difference in discounted QALY", "Difference in discounted cost (1M USD)", "Difference in discounted QALY")
	col_xbar_vec = c("F"="red", "M"="blue")
	pch_xbar_vec = c("F"=1, "M"=2)
	for (iplot in 1:length(xplot_list)){
		xbar_plot = xplot_list[[iplot]]
		xviolin_plot = xviolin_list[[iplot]]
		xbar_lab = xplot_lab_vec[iplot]
		xbar_xlab = xbar_xlab_vec[iplot]
		
		col_xbar = rep(col_xbar_vec["F"], length(iout_plot))
		col_xbar[grep("^M_", outlab_pie[iout_plot])] = col_xbar_vec["M"]
		pch_xbar = rep(pch_xbar_vec["F"], length(iout_plot))
		pch_xbar[grep("^M_", outlab_pie[iout_plot])] = pch_xbar_vec["M"]
		
		for (iVC in 2:VC_boys_num){
			xbar_plot_iVC = xbar_plot[[iVC]]
			xviolin_plot_iVC = do.call(cbind, lapply(xviolin_plot, function(xmat) xmat[iVC,]))																	
			lowbound0_YN = !(diff_calc_1dose_YN==TRUE & numDose_Vacc_new==1);
			if (lowbound0_YN){
				# non-zero 
				xbar_plot_iVC = pmax(xbar_plot_iVC, 0); 
			}
			if (grepl("cost", xbar_lab, ignore.case=TRUE)){
				xbar_plot_iVC = xbar_plot_iVC/7.8/1e6;
				xviolin_plot_iVC = xviolin_plot_iVC/7.8/1e6
			}
			
			if (lowbound0_YN){
				xlim_plot = range(pretty(c(0, max(xbar_plot_iVC)), n=3))
			} else {
				xlim_plot = range(pretty(range(xbar_plot_iVC), n=3))
			}
			ylim_plot = rev(range(iout_plot) + 0.15*c(-1,1))
			plot(NA, xlim=xlim_plot, ylim=ylim_plot, axes=FALSE, xlab="",ylab="")
			axis(1, cex.axis=1.45)
			if (FALSE){
				# bar plot
				points(x=xbar_plot_iVC[,1], y=iout_plot, col=col_xbar, pch=pch_xbar, cex=2)
				segments(x0=xbar_plot_iVC[,2], x1=xbar_plot_iVC[,3], y0=iout_plot, col=col_xbar, lwd=2.5)
			} else{ 
				# violin plot; may not be good to use vioplot which involves kernel density, because the guassian kernel may give negative value i.e., lower than the min (0)
				y_yplot = xviolin_plot_iVC # use "y_" because vioplot plot horizontal=FALSE by default
				probs_boxplot = c(0.05, 0.25, 0.5, 0.75, 0.95)
				summ_pct_out = apply(y_yplot, 2, quantile, probs=probs_boxplot)
				if (lowbound0_YN){ # non-zero 
					y_yplot = pmax(y_yplot, 0); 
					summ_pct_out = pmax(summ_pct_out, 0); 
				}
				col_xbar_alpha = adjustcolor(col_xbar, alpha.f=0.5)
				sepdraw_box_YN = TRUE
				vioplot::vioplot(y_yplot, at=iout_plot, add=TRUE, horizontal=TRUE, col=col_xbar_alpha, wex=0.3, axes=FALSE, drawRect=(!sepdraw_box_YN))
				summ_pct_out = list(stats=summ_pct_out, n=rep(100, ncol(summ_pct_out))) # $n is needed fo bxp but have no effect to the plot
				bxp(summ_pct_out, at=iout_plot, horizontal=TRUE, add=TRUE, axes=FALSE, boxfill="grey80", medpch=NA, medcol="black", medcex=1.5, medlwd=2, boxwex=0.125, whisklty=1, whisklwd=2, whiskcol="green3", staplelty="blank") # medlty="blank" if not to draw the line for median; medpch=16 for using a circle to indicate the median
			}
			
			mtext(side=1, text=xbar_xlab, line=2.75)
			text(x=par()$usr[1], y=iout_plot, labels=outlab_pie[iout_plot], adj=1, xpd=TRUE, cex=1.5)

			mtext(text=sprintf("GNV %d%% uptake for boys", VC_boys_vec[iVC]), side=3, line=-0.5) # , adj=0.05 
			if (iVC==2){
				mtext(text=xbar_lab, side=3, line=1.25, adj=0, at=with(par(), usr[1]-diff(usr[1:2])/8), cex=0.8)
			}
		} # iVC
	} # iplot
	# end plot diff_Costsave/QALYgain_vsNoVacc by health outcome at iVdur


	
} # for-iVdur


# export Cost_overall_array and QALY_overall_array
fname_RData_out2 = gsub('.RData', sprintf('_CostQALY_overallOnly_gw%s.RData',ifelse(calc_gwarts_YN,"Y","N")), fname_RData_out, fixed=TRUE)
save(list=c("Cost_overall_array", "QALY_overall_array"), file=fname_RData_out2) # use save(), not save.iamge for saving some variables only


} # for-iprop_MaleCx

} # for-run_SensAnaly_propMaleCx_YN


} # for-numDose_Vacc_new
} # for-calc_gwarts

sink()
if (plot_PDF_YN) {dev.off()}


# save Cost QALY of 2F1M, for 2F1M-GNV vs 1-dose GNV
Cost_2F1M_list = list(
	Cost_vaccGirls_array_TEMP = Cost_Vacc_adjfactor*Cost_vaccGirls_array,
	Cost_vaccBoys_array_TEMP = Cost_Vacc_adjfactor_1dose_boys*Cost_vaccBoys_array,

	Cost_array = Cost_array,
	Costsave_nonCeCx_array = Costsave_nonCeCx_array,
	
	Cost_overall_array = Cost_array + Cost_vaccGirls_array_TEMP + Cost_vaccBoys_array_TEMP - Costsave_nonCeCx_array,
	Cost_overall_array_woVacc = Cost_array - Costsave_nonCeCx_array, 
															
	Cost_overall_array_onlyVacc_girlsboys_noadjust = (numDose_Vacc_new/numDose_Vacc)*Cost_vaccGirls_array + (numDose_Vacc_new_Boys/numDose_Vacc)*Cost_vaccBoys_array,
	
	Costsave_nonCeCx_array_byType = Costsave_nonCeCx_array_byType,
	Costsave_nonCeCx_array_byType_noMaleCxFOV = Costsave_nonCeCx_array_byType_noMaleCxFOV # noMaleCxFOV
)	
QALY_2F1M_list = list(
	QALY_array = QALY_array,
	QALYgain_nonCeCx_array = QALYgain_nonCeCx_array,
	LYgain_nonCeCxDeath_array = LYgain_nonCeCxDeath_array,

	QALY_overall_array = QALY_array + QALYgain_nonCeCx_array + LYgain_nonCeCxDeath_array,
	
	QALYgain_nonCeCx_array_byType = QALYgain_nonCeCx_array_byType,
	QALYgain_nonCeCx_array_byType_noMaleCxFOV = QALYgain_nonCeCx_array_byType_noMaleCxFOV,

	LYgain_nonCeCxDeath_array_byType = LYgain_nonCeCxDeath_array_byType,
	LYgain_nonCeCxDeath_array_byType_noMaleCxFOV = LYgain_nonCeCxDeath_array_byType_noMaleCxFOV
)

# gwarts
Cost_2F1M_list = modifyList(Cost_2F1M_list, 
	list(Cost_overall_array_gwarts = Cost_2F1M_list$Cost_overall_array - Reduce("+", lapply(Costsave_gwarts_array_bysex, function(xlist) xlist)), 
		Cost_overall_array_woVacc_gwarts = Cost_2F1M_list$Cost_overall_array_woVacc - Reduce("+", lapply(Costsave_gwarts_array_bysex, function(xlist) xlist)),
		Costsave_gwarts_array_sumsex = Reduce("+", lapply(Costsave_gwarts_array_bysex, function(xlist) xlist)) 
	)
 
)
QALY_2F1M_list = modifyList(QALY_2F1M_list, 
	list(QALY_overall_array_gwarts = QALY_2F1M_list$QALY_overall_array + Reduce("+", lapply(QALYgain_gwarts_array_bysex, function(xlist) xlist)),
		QALYgain_gwarts_array_sumsex = Reduce("+", lapply(QALYgain_gwarts_array_bysex, function(xlist) xlist))
	)
)

fname_RData_out2 = gsub('.RData', '_CostQALY.RData', fname_RData_out, fixed=TRUE)
save(list=c("Cost_2F1M_list", "QALY_2F1M_list"), file=fname_RData_out2) # use save(), not save.iamge for saving some variables only


# export ICERs and TVCs for all simulations
fname_RData_out3 = gsub('.RData', '_ICER_TVC.RData', fname_RData_out, fixed=TRUE)
save(list=c("ICERall_out", "TVCall_out", "prop_MaleCx_MSM"), file=fname_RData_out3)
