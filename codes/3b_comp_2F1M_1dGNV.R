# compare 2F1M-GNV vs 1-dose GNV

# variable names in RData_2F1M and RData_1dGNV should be different

# 2F1M, saved in iVdur==1 (xVdur 100)
iVdur_2F1M = 1;
outputdate_2F1M = outputdate_string
folder_2F1M = gsub("(/){2,}", "/", sprintf("%s/2F1M/", output_folder_CEA))
fnameRData_2F1M = sprintf("calc_CEA_nonCeCx_%s_CostQALY.RData", outputdate_2F1M);


# 1dGNV, saved in iVdur==2
iVdur_1dGNV = 2;
outputdate_1dGNV = outputdate_string
folder_1dGNV = gsub("(/){2,}", "/", sprintf("%s/GNV/", output_folder_CEA))
fnameRData_1dGNV = sprintf("calc_CEA_nonCeCx_%s_CostQALY.RData", outputdate_1dGNV);

# output setting
outputdate_comp = outputdate_string
folder_compare = gsub("(/){2,}", "/", sprintf("%s/comp_2F1M_1dGNV/", output_folder_CEA))
fname_out = sprintf("%scomp_2F1M_1dGNV_%s.txt", folder_compare, outputdate_comp)
fname_plot_pdf = gsub('.txt', '.pdf', fname_out, fixed=TRUE)
if (!dir.exists(folder_compare)){
	dir.create(folder_compare)
}



# load variables
load(paste0(folder_2F1M,fnameRData_2F1M))
load(paste0(folder_1dGNV,fnameRData_1dGNV))


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


# settings
VC_boys_vec = c(0, 25, 50, 85)# include VC=0 for comparison
VC_boys_num = length(VC_boys_vec)
iVC_comp = 2:VC_boys_num
xVC_comp = VC_boys_vec[iVC_comp]
xVC_comp_num = length(xVC_comp);

nonCeCx_num = 6;
Vdur_num = 2;
iPSA_num = 100;
PSA_costQALY_num = 100;


# calculate ICER
WTPthreshold = 406538; # be aware of unit/currency
plot_Cost_USD_YN = TRUE;


# reference vaccine cost, to calculate the relative change in TVC
Cost_Vacc_perdose_inclAdmin_ref = 1700; # Family Planning Association $1700 for vaccine itself + admin costs (assume $250) >> Cost_Vacc_perdose_new = 1450

	Cost_Vacc_perdose_tender_plusAdmin_USD = 177; # <-- change this; new vaccination cost, vaccine cost based on the tender price plus admin
	Cost_Vacc_perdose_inclAdmin_ref = Cost_Vacc_perdose_tender_plusAdmin_USD*7.8; # new vaccination cost

# ratio of vaccination cost used in the RData
Cost_vacc_use_new = 1700; # cost in RData, 
Cost_vacc_use_base = 1400; # not to change, cost in Cpp

	Cost_vacc_use_new = Cost_Vacc_perdose_inclAdmin_ref;

ratioVaccCost = Cost_vacc_use_new/Cost_vacc_use_base;

# summary() work on matrix/data.frame too. so apply summary_new() to vector only
quant_probs = c(0.5, 0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 1)
summary_new = function(x) { if (TRUE) as.table(c("Mean"=mean(x), quantile(x, quant_probs))) else summary(x) }

	# checking
	if (TRUE){ # check ratioVaccCost

		ratioVaccCost_2F1M = ratioVaccCost;
		Cost_T = Cost_2F1M_list$Cost_overall_array[, iVdur_2F1M, ];
		Cost_woVacc_T = Cost_2F1M_list$Cost_overall_array_woVacc[, iVdur_2F1M, ]
		Cost_onlyVacc_T = Cost_2F1M_list$Cost_overall_array_onlyVacc_girlsboys_noadjust[, iVdur_2F1M, ]
		print( "check cost, 2F1M" )
		print( t(apply((ratioVaccCost_2F1M*Cost_onlyVacc_T + Cost_woVacc_T) - Cost_T,1, summary)) )
		

		# check the version of the code "calc_CEA_xx.R", newly version set Cost_onlyVacc_noadjust_1dGNV = Cost with 1dGNV already
		ratioVaccCost_1dGNV = 1; # 1/2*ratioVaccCost
		Cost_T = Cost_1dGNV_list$Cost_overall_array[, iVdur_1dGNV, ];
		Cost_woVacc_T = Cost_1dGNV_list$Cost_overall_array_woVacc[, iVdur_1dGNV, ]
		Cost_onlyVacc_T = Cost_1dGNV_list$Cost_overall_array_onlyVacc_girlsboys_noadjust[, iVdur_1dGNV, ]

		print( "check cost, 1dGNV" )
		print( t(apply((ratioVaccCost_2F1M*Cost_onlyVacc_T + Cost_woVacc_T) - Cost_T,1, summary)) )
		print( t(apply((1/2*ratioVaccCost_2F1M*Cost_onlyVacc_T + Cost_woVacc_T) - Cost_T,1, summary)) )
		print( t(apply((ratioVaccCost_1dGNV*Cost_onlyVacc_T + Cost_woVacc_T) - Cost_T,1, summary)) )


		
		# the cost for 2dFOV should be the same in 2F1M and 1dGNV
		print( 'show the cost for 2dFOV in 2F1M and 1dGNV' )
		print( 'cost, 2d-FOV, 2F1M' ) # the cost of 2dFOV in the 2F1M output
		print( summary(Cost_2F1M_list$Cost_overall_array[1, iVdur_2F1M, ]) )
		print( 'cost, 2d-FOV, 1dGNV' ) # the cost of 2dFOV in the 1dGNV output
		print( summary(Cost_1dGNV_list$Cost_overall_array[1, iVdur_1dGNV, ]) )

		print( 'the cost for 2dFOV should be the same in 2F1M and 1dGNV' )
		print( "check the outcomes before further proceeding (ratioVaccCost)" )
		browser()
	}


	if (TRUE){ # check cost, QALY, gwarts
		print( 'show the cost and QALY regarding gwarts' )
		Cost_diff_gwarts_wogwarts_2F1M = Cost_2F1M_list$Cost_overall_array_gwarts[iVC_comp, iVdur_2F1M, ] - Cost_2F1M_list$Cost_overall_array[iVC_comp, iVdur_2F1M, ]
		Cost_diff_gwarts_wogwarts_1dGNV = Cost_1dGNV_list$Cost_overall_array_gwarts[iVC_comp, iVdur_1dGNV, ] - Cost_1dGNV_list$Cost_overall_array[iVC_comp, iVdur_1dGNV, ]
		print( "Cost_diff_gwarts_wogwarts_2F1M" )
		print( t(apply(Cost_diff_gwarts_wogwarts_2F1M,1, summary_new)) )
		print( "Cost_diff_gwarts_wogwarts_1dGNV" )
		print( t(apply(Cost_diff_gwarts_wogwarts_1dGNV,1, summary_new)) )

		QALY_diff_gwarts_wogwarts_2F1M = QALY_2F1M_list$QALY_overall_array_gwarts[iVC_comp, iVdur_2F1M, ] - QALY_2F1M_list$QALY_overall_array[iVC_comp, iVdur_2F1M, ]
		QALY_diff_gwarts_wogwarts_1dGNV = QALY_1dGNV_list$QALY_overall_array_gwarts[iVC_comp, iVdur_1dGNV, ] - QALY_1dGNV_list$QALY_overall_array[iVC_comp, iVdur_1dGNV, ]
		print( "QALY_diff_gwarts_wogwarts_2F1M" )
		print( t(apply(QALY_diff_gwarts_wogwarts_2F1M,1, summary_new)) )
		print( "QALY_diff_gwarts_wogwarts_1dGNV" )
		print( t(apply(QALY_diff_gwarts_wogwarts_1dGNV,1, summary_new)) )
		
		print( '1dGNV should have fewer QALYgain_diff than 2F1M' )
		print( "check the outcomes before further proceeding (cost, QALY, gwarts)" )
		browser()
	}	


output_template_overall_PSAcostQALY = array(0, dim=c(VC_boys_num, Vdur_num, iPSA_num*PSA_costQALY_num)); 


# opt_width_default = 140 # options()$width
options(width = 200)
sink(fname_out)
print( 'compare 2F1M vs 1dGNV' )

plot_PDF_YN = TRUE;
num_subplot = VC_boys_num - 1;
if (plot_PDF_YN) {pdf(fname_plot_pdf, width=num_subplot*5, height=5)}


calc_gwarts_YN_vec = c(0, 1);


# save all ICERs and TVCs by calc_gwarts_YN_vec(2), numDose_Vacc_new_vec (2), Vdur (2)
Vdur_vec = c(100, 20)
ICERall_out = array(list(NULL), dim=c(2, 2, 2));
dimnames(ICERall_out) = list(paste0("gwarts_",0:1), paste0("numDose_Vacc_",1:2), paste0("Vdur_",Vdur_vec) )

TVCall_out = ICERall_out;

# indexing for output
ii_numDose_Vacc = 2;
ii_Vdur = 1;


for (ixlim_plot in 1:2){ # whether to use the same xlim_plot
# ixlim_plot==2 for using the same xlim_plot


for (calc_gwarts_YN in calc_gwarts_YN_vec){
	ii_gwarts_ICER_all_out = calc_gwarts_YN+1; # ICERall_out/TVCall_out

	if (ixlim_plot==1){
		cat("\n"); print( sprintf('gwarts=%s', ifelse(calc_gwarts_YN,'Y','N')) )
	}
	
	# pick VC = iVC_comp, exclude VCboys=0
	if (calc_gwarts_YN==0){
		Cost_2F1M = Cost_2F1M_list$Cost_overall_array[iVC_comp, iVdur_2F1M, ]
		QALY_2F1M = QALY_2F1M_list$QALY_overall_array[iVC_comp, iVdur_2F1M, ]
		Cost_1dGNV = Cost_1dGNV_list$Cost_overall_array[iVC_comp, iVdur_1dGNV, ]
		QALY_1dGNV = QALY_1dGNV_list$QALY_overall_array[iVC_comp, iVdur_1dGNV, ]
	} else if (calc_gwarts_YN==1){
		Cost_2F1M = Cost_2F1M_list$Cost_overall_array_gwarts[iVC_comp, iVdur_2F1M, ]
		QALY_2F1M = QALY_2F1M_list$QALY_overall_array_gwarts[iVC_comp, iVdur_2F1M, ]
		Cost_1dGNV = Cost_1dGNV_list$Cost_overall_array_gwarts[iVC_comp, iVdur_1dGNV, ]
		QALY_1dGNV = QALY_1dGNV_list$QALY_overall_array_gwarts[iVC_comp, iVdur_1dGNV, ]
	}
	
	
	# for TVC calculation, Cost woVacc and onlyVacc
	Cost_onlyVacc_2F1M = ratioVaccCost_2F1M*Cost_2F1M_list$Cost_overall_array_onlyVacc_girlsboys_noadjust[iVC_comp, iVdur_2F1M, ]
	Cost_woVacc_2F1M = Cost_2F1M - Cost_onlyVacc_2F1M;
	Cost_onlyVacc_noadj_2F1M = Cost_onlyVacc_2F1M/ratioVaccCost;
	
	
	Cost_onlyVacc_1dGNV = ratioVaccCost_1dGNV*Cost_1dGNV_list$Cost_overall_array_onlyVacc_girlsboys_noadjust[iVC_comp, iVdur_1dGNV, ]
	Cost_woVacc_1dGNV = Cost_1dGNV - Cost_onlyVacc_1dGNV;
	Cost_onlyVacc_noadj_1dGNV = Cost_onlyVacc_1dGNV/ratioVaccCost;

	
	Cost_diff = Cost_2F1M - Cost_1dGNV;
	QALY_diff = QALY_2F1M - QALY_1dGNV;
	ICER = Cost_diff / QALY_diff;
	rownames(Cost_diff) = xVC_comp; # VC_boys_vec
	rownames(QALY_diff) = rownames(Cost_diff);
	rownames(ICER) = rownames(Cost_diff)
	
	
	if (ixlim_plot==1){ # print once only; ixlim_plot varies the xlim_plot for plotting
		print( "summary, Cost_2F1M" )
		print( t(apply(Cost_2F1M,1, summary_new)) )
		print( "summary, QALY_2F1M" )
		print( t(apply(QALY_2F1M,1, summary_new)) )
		print( "summary, Cost_1dGNV" )
		print( t(apply(Cost_1dGNV,1, summary_new)) )
		print( "summary, QALY_1dGNV" )
		print( t(apply(QALY_1dGNV,1, summary_new)) )	
		
		print( "summary, Cost_diff" )
		print( t(apply(Cost_diff,1, summary_new)) )
		print( "summary, QALY_diff" )
		print( t(apply(QALY_diff,1, summary_new)) )
		print( "summary ICER ")
		print( t(apply(ICER,1, summary_new)) )
		print( "summary ICER (USD) ")
		print( t(apply(ICER,1, summary_new))/7.8 )
		print( "prop%(ICER<WTPthreshold)" )
		print( rowMeans(ICER<WTPthreshold) )

		# TVC calculation
		QALYmoney = WTPthreshold * QALY_diff;
		Cost_diff_woVacc = Cost_woVacc_2F1M - Cost_woVacc_1dGNV;
		Cost_diff_onlyVacc = Cost_onlyVacc_2F1M - Cost_onlyVacc_1dGNV;
		Cost_diff_onlyVacc_noadj = Cost_onlyVacc_noadj_2F1M - Cost_onlyVacc_noadj_1dGNV
		# identical(Cost_diff, (Cost_diff_woVacc + Cost_diff_onlyVacc))
		
		TVCratio_T = (QALYmoney - Cost_diff_woVacc)/(Cost_diff_onlyVacc_noadj);
		TVC = TVCratio_T * Cost_vacc_use_base;
		print( "TVC" )
		print( t(apply(TVC, 1, summary_new)) )
		print( "TVC (USD)" )
		print( t(apply(TVC, 1, summary_new))/7.8 )
		# print( t(apply(TVCratio_T, 1, summary_new)) )
		
		
		# export ICERs and TVCs for all simulations
		ICERall_out[[ii_gwarts_ICER_all_out, ii_numDose_Vacc, ii_Vdur]] = ICER
		TVCall_out[[ii_gwarts_ICER_all_out, ii_numDose_Vacc, ii_Vdur]] = TVC
		
		
		if (FALSE){ # should equal WTPthreshold at TVC 
			newICER_check = (Cost_diff_woVacc + TVC2/Cost_vacc_use_base * Cost_diff_onlyVacc_noadj)/QALY_diff
			print( t(apply(newICER_check, 1, summary_new)) )
		}

		# print summary stat
		summstat_Temp = cbind("ICER_USD"=apply(ICER/7.8, 1, print_summstat, digit=0),
			"TVC_USD"=apply(TVC/7.8, 1, print_summstat),
			"relChange_TVC"=apply(TVC/Cost_Vacc_perdose_inclAdmin_ref*100, 1, print_summstat),
			"relReduce_TVC"=apply((1-TVC/Cost_Vacc_perdose_inclAdmin_ref)*100, 1, print_summstat)
			)
		print( "summary stat" )
		print(summstat_Temp)	
	
	} # if- (ixlim_plot==1)


	# plot
	# scatter plot
	pch_plot = 16; 1; 

	Cost_plotunit = 1/1000000
	lab_Cost_plotunit = " (per 1M)"
	if (plot_Cost_USD_YN){ # USD
		Cost_plotunit = Cost_plotunit/7.8
		lab_Cost_plotunit = " (per 1M USD)"	
	}	

	Cost_diff_plot = Cost_diff * Cost_plotunit
	ylab_plot = paste0("Difference in discounted cost", lab_Cost_plotunit)
	xlab_plot = "Difference in discounted QALY" #  for GNV vs FOV

	xlim_plot = range(QALY_diff)
	xlim_plot_prettyn = (2:3)[which.min(sapply(2:3, function(n) max(pretty(xlim_plot, n=n))))]
	xlim_plot = range(pretty(xlim_plot, n=xlim_plot_prettyn))
	ylim_plot = range(Cost_diff_plot)
	ylim_plot_prettyn = (2:3)[which.min(sapply(2:3, function(n) max(pretty(ylim_plot, n=n))))]
	ylim_plot = range(pretty(ylim_plot, n=ylim_plot_prettyn))
	if (ixlim_plot==2){ # use the same range for 2F1M-GNV vs 2d-FOV
		xlim_plot = c(0, 10000); 
		ylim_plot = c(0, 150);
	}

	if (plot_PDF_YN==FALSE){windows(width=num_subplot*5, height=1*5)}
	par(mar=c(4,5.75,3,0.75))
	layout(matrix(1:num_subplot, nrow=1))
	pch_col_alpha = adjustcolor('blue', alpha=0.25)


	for (iVC in 1:xVC_comp_num){
		xplot = QALY_diff[iVC,]
		yplot = Cost_diff_plot[iVC,]
		plot(x=xplot, y=yplot, xlim=xlim_plot, ylim=ylim_plot, xlab="", ylab="", bty="L", las=1, cex=0.9, col=pch_col_alpha, pch=pch_plot, axes=FALSE)
			abline(v=0, h=0, col='grey', lty=2) # reference lines
			abline(a=0, b=WTPthreshold*Cost_plotunit, col='green')
		axis(1, at=par()$usr[1:2], label=FALSE, col.ticks=NA)
		axis(1, col=NA, col.ticks='black', padj=-0.15, cex.axis=1.5)
		axis(2, at=par()$usr[3:4], label=FALSE, col.ticks=NA)
		axis(2, col=NA, col.ticks='black', las=1, hadj=0.95, cex.axis=1.5)
		mtext(side=1, xlab_plot, cex=1, line=2.5)
		mtext(side=2, ylab_plot, cex=1, line=4)

		mtext(text=sprintf("GNV %d%% uptake for boys", xVC_comp[iVC]), side=3, line=0.15) # , adj=0.05 
		mtext(text="2F1M-GNV vs 1-dose GNV", side=3, line=1.5) # , adj=0.05 
		points(x=mean(xplot), y=mean(yplot), pch=17, col="orange", cex=3)
		
		# cost-effective, Q1
		# % prop of CE
		leg_text = sprintf("CE: %.1f%%", 100*mean((yplot/Cost_plotunit)/xplot<WTPthreshold))
		with(par(), text(x=xaxp[2], y=0+diff(yaxp[1:2])/10, labels=leg_text, adj=c(1,0), cex=1.45))
		
		
		mtext_at = par()$usr[1]
		if (iVC==1){ # add text when plotting the first plot
			if (TRUE){ # calc_gwarts_YN){
				mtext(text=sprintf('gwarts=%s',ifelse(calc_gwarts_YN,"Y","N")), side=3, line=1, cex=0.7, adj=1, at=mtext_at)
			}
			
		}
	} # for-iVC	
	
} # for-calc_gwarts_YN



} # for- xlim_plot 


sink()
if (plot_PDF_YN) { dev.off()}


# export ICERs and TVCs
fname_out = gsub("(/){2,}", "/", sprintf("%s/comp_2F1M_1GNV_ICERs_TVCs_%s.RData", folder_compare, outputdate_comp))
save(list=c("ICERall_out", "TVCall_out"), file=fname_out)
