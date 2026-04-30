# do the plots for MCMC outcomes


# 2D histogram; https://www.r-bloggers.com/2014/09/5-ways-to-do-2d-histograms-in-r/

##### Addendum: 2D Histogram + 1D on sides (from Computational ActSci w R) #######
#http://books.google.ca/books?id=YWcLBAAAQBAJ&pg=PA60&lpg=PA60&dq=kde2d+log&source=bl&ots=7AB-RAoMqY&sig=gFaHSoQCoGMXrR9BTaLOdCs198U&hl=en&sa=X&ei=8mQDVPqtMsi4ggSRnILQDw&redir_esc=y#v=onepage&q=kde2d%20log&f=false

hist2d_1dsides = function(dataf){
	h1 <- hist(dataf$x, breaks=25, plot=F)
	h2 <- hist(dataf$y, breaks=25, plot=F)
	top <- max(h1$counts, h2$counts)
	k <- kde2d(dataf$x, dataf$y, n=25)

	# margins
	oldpar <- par()
	par(mar=c(3,3,1,1))
	layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
	image(k, col=r) #plot the image
	par(mar=c(0,2,1,0))
	barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
	par(mar=c(2,0,0.5,1))
	barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
} # end hist2d_1dsides


##### OPTION 3: stat_bin2d from package 'ggplot' #######
hist2d_gg = function(dataf, varnames){

	p <- ggplot(dataf, aes(.data[[varnames[1]]], .data[[varnames[2]]]))
	h3 <- p + stat_bin2d(bins=50) + scale_fill_gradientn(colours=rev(hcl.colors(10, palette="YlOrRd")))
	# scale_fill_gradient(low='yellow', high='red')
	h3 <- h3 + theme_classic()

	return(h3)
} # hist2d_gg


library(ggplot2) # for 2d histogram with legend


RDatafolder_input = output_folder_fit_CxInc 

RData_input_mat = as.data.frame( cbind(outputdate_string, rbind(c("F", "VV"), c("F", "OPC"), c("F", "anus"),   c("M", "penis"), c("M", "OPC"), c("M", "anus"))) )
Rversion = do.call(paste, c(sessionInfo()$R.version[c('major','minor')], sep='.'))
if (Rversion=="4.0.3"){
	# simplify=TRUE/FALSE for apply() available after 4.1.0
	RData_input_list = sapply(1:nrow(RData_input_mat), simplify=FALSE, function(i) as.vector(RData_input_mat[i,]))
} else {
	RData_input_list = apply(RData_input_mat,1, simplify=FALSE, function(x) as.vector(x)) # simplify=T/F for apply() available after 4.4.0
}

RData_vec = sapply(RData_input_list, function(x) do.call(sprintf, c(fmt="estoutput_%s_fit_CxInc_%s_%s.RData", as.list(x))))

# check if RData exists
check_RData = sapply(RData_vec, function(xRData) file.exists(paste0(RDatafolder_input, xRData)));
if (any(!check_RData)){
	stop( sprintf("The following RData file does not exist: %s\nCheck before continue.\n", paste(  RData_vec[which(!check_RData)], collapse=",")) )
}



# export extracted parameter sets
export_Extpara = rep(list(NULL), length(RData_input_list));
names(export_Extpara) = sapply(RData_input_list, function(x) paste(x[2], x[3], sep="_"))
names(export_Extpara) = gsub("VV", "vagina_vulva", names(export_Extpara))
fname_export_Extpara = sprintf("%s%s_nonCeCx_PSApara.xlsx", RDatafolder_input, outputdate_string)


# start iRData
for (iRData in 1:length(RData_vec)){


RDataname_input = RData_vec[iRData]
time1_load = Sys.time() # the RData may also contain the variable time1 which will then overwrite the original time1 if they are of the same name
load( gsub("//","/", paste0(RDatafolder_input, RDataname_input)) )
time2_load = Sys.time()
timediff = time2_load - time1_load
cat( sprintf("Time to load RData: %f %s\n", as.numeric(timediff), attributes(timediff)$units) )

# base fname for plots
fname_pdfplot_base = sprintf("%s%%s_%s.pdf", RDatafolder_input, output_suffix)


# note: for gamma distribution, alpha=shape, beta=rate=1/scale; mean = shape x scale or alpha/beta
# varname_vec = as.vector( do.call(c, sapply(paste0('HPV_',1:3), simplify=FALSE, function(x) paste(x, c('gamma_shape', 'gamma_scale', 'scaling1', 'scaling2'), sep="_"))) )

varname_vec = as.vector( do.call(c, sapply(paste0('HPV_',1:3), simplify=FALSE, function(x) paste(x, c('gamma_shape', 'gamma_scale', paste0('scaling',1:n_scaling)), sep="_"))) )

iparamest_list_vec = 1; #RData may contain fitting of more than 1 CxInc


# call back mcmc setting
# vec_HPV_set = 1:100 # retrieve the one used in the .RData, unless another set is to be used
MCMCout_estparam_temp = MCMCout_estparam_list[[1]]
nrow_allMCMC = nrow(MCMCout_estparam_temp)
num_vec_HPV_set = length(vec_HPV_set)
nrow_eachMCMC = nrow_allMCMC/num_vec_HPV_set
irow_eachMCMC = sapply(1:num_vec_HPV_set, simplify=FALSE, function(iPSA) nrow_eachMCMC*(iPSA-1) + (1:nrow_eachMCMC)) # irow of each MCMC in the original output (before thinning))

mcmc_iterations = nrow_allMCMC/num_vec_HPV_set
mcmc_thinning = min(50, mcmc_iterations);
mcmc_thin_vec = seq(from=mcmc_thinning, to=nrow_allMCMC, by=mcmc_thinning);
mcmc_thin_YN = TRUE;

# irow index for picking from each MCMC output
ipick_temp = floor(length(mcmc_thin_vec)/num_vec_HPV_set)
irow_pick_from_MCMC = seq(ipick_temp, by=ipick_temp, length.out=num_vec_HPV_set)


pctstat_vec = c(0.5, 0.025,0.975, 0.05,0.95, 0.1,0.9, 0.25,0.75);
fun_pctstat = function(x, pct=pctstat_vec){ quantile(x, pctstat_vec) };


nparam_iset =  2*n_delay + n_scaling + n_RRprog

MCMCout_estparam_temp_list_all = rep(list(NULL), length(iparamest_list_vec))

adjrelative_param_YN = FALSE; # some model settings consider relative parameterts


## export the parameters
for (ilist in iparamest_list_vec){

	MCMCout_estparam_temp = MCMCout_estparam_list[[ilist]]
	MCMCout_estparam_temp = MCMCout_estparam_temp[mcmc_thin_vec,]

	stat_MCMCout_estparam_temp = apply(MCMCout_estparam_temp, 2, fun_pctstat)
	rownames(stat_MCMCout_estparam_temp) = paste0("pct_", pctstat_vec)
	MCMCout_estparam_temp_list = list(estparam=MCMCout_estparam_temp, stat_estparam=stat_MCMCout_estparam_temp)

	if (adjrelative_param_YN){ # be careful
		print( "update relative parameters" )
		adjust_MCMCout_estparam_temp = MCMCout_estparam_temp;
		adjust_MCMCout_estparam_temp[, 4:6] = adjust_MCMCout_estparam_temp[, 1:3] * adjust_MCMCout_estparam_temp[, 4:6]
		adjust_MCMCout_estparam_temp[, 7:9] = adjust_MCMCout_estparam_temp[, 1:3] * adjust_MCMCout_estparam_temp[, 7:9]
		# calculate the mean duration
		adjust_MCMCout_estparam_temp = cbind(adjust_MCMCout_estparam_temp 
			, adjust_MCMCout_estparam_temp[,1]*adjust_MCMCout_estparam_temp[,2]
			, adjust_MCMCout_estparam_temp[,4]*adjust_MCMCout_estparam_temp[,5]
			, adjust_MCMCout_estparam_temp[,7]*adjust_MCMCout_estparam_temp[,8])
	
		stat_adjust_MCMCout_estparam_temp = apply(adjust_MCMCout_estparam_temp, 2, fun_pctstat)
		rownames(stat_adjust_MCMCout_estparam_temp) = paste0("pct_", pctstat_vec)
		MCMCout_estparam_temp_list = append(MCMCout_estparam_temp_list, list(adjestparam=adjust_MCMCout_estparam_temp, stat_adjestparam=stat_adjust_MCMCout_estparam_temp))
	} else{
		# if not fitting relative param, calculate the mean duration and relative proportion of progression directly from the estparam
		meandur_MCMCout_estparam_temp = cbind(
			MCMCout_estparam_temp[,1]*MCMCout_estparam_temp[,2]
			,MCMCout_estparam_temp[,5]*MCMCout_estparam_temp[,6]
			,MCMCout_estparam_temp[,9]*MCMCout_estparam_temp[,10] 
			
			,MCMCout_estparam_temp[,3]/MCMCout_estparam_temp[,7]
			,MCMCout_estparam_temp[,4]/MCMCout_estparam_temp[,8]
			,MCMCout_estparam_temp[,3]/MCMCout_estparam_temp[,11]
			,MCMCout_estparam_temp[,4]/MCMCout_estparam_temp[,12]
			)
		
		stat_meandur_MCMCout_estparam_temp = apply(meandur_MCMCout_estparam_temp, 2, fun_pctstat)
		rownames(stat_meandur_MCMCout_estparam_temp) = paste0("pct_", pctstat_vec)
		MCMCout_estparam_temp_list = append(MCMCout_estparam_temp_list, list(meandur=meandur_MCMCout_estparam_temp, stat_meandur=stat_meandur_MCMCout_estparam_temp))
	}

	if (TRUE){ # generate the proportion of HPV types among cancers 
		HPVincid_ilist = HPVincid_tofit_list[[ilist]]
		CxInc_allPSA_list = rep(list(NULL), length=nrow(HPVincid_ilist))
		for (iPSA in vec_HPV_set){
			HPVincid_iPSA = HPVincid_ilist[iPSA,]
			x_iPSA = MCMCout_estparam_temp[ irow_pick_from_MCMC[iPSA], ]

			xdata_iPSA = list('HPVincid_data'=HPVincid_iPSA/100, 
				'outcome_agegp_diff_map'=outcome_agegp_diff_map, 
				'pop_agegp'=pop_agegp_temp,
				'paraest_setting'=paraest_setting,
				'mat_ratioAlive'=mat_ratioAlive_temp
				);

			extInput_iPSA = list('xdata'=xdata_iPSA, 'CxInc_data'=CxInc_datafit,
				'idx_CxInc_agegp_data'=idx_CxInc_outcomeAgegp_data, 'idx_CxInc_agegp_calib'=idx_CxInc_outcomeAgegp_calibrate,
				'logL_const'=logpoisspdf_C(CxInc_datafit), # log of the constant in poisson distribution
				'extInput_9vHR'=extInput_9vHR
				);
		
			CxInc_iPSA = fun_gen_CxInc_100k_byType(x_iPSA, extInput_iPSA, gen_pRC_YN=TRUE);
			
			
			CxInc_iPSA_temp = Reduce("+", CxInc_iPSA$CxInc) 
			# CxInc_iPSA_pRC_temp = with(as.list(CxInc_iPSA$CxInc_pRC), c(HPV_16=HPV_16, HPV_all9vHR=HPV_16+HPV_18_o9vHR, HPV_nonV=HPV_nonV))
			CxInc_iPSA_pRC_temp = c(CxInc_iPSA$CxInc_pRC, "HPV_all9vHR"=Reduce("sum", CxInc_iPSA$CxInc[1:3])/Reduce("sum", CxInc_iPSA$CxInc))
			
			CxInc_allPSA_list[[iPSA]] = c(CxInc_iPSA_temp, CxInc_iPSA_pRC_temp)

		} # iPSA
		CxInc_allPSA_list = do.call(rbind, CxInc_allPSA_list)

		stat_CxInc_allPSA_list = apply(CxInc_allPSA_list, 2, fun_pctstat)
		rownames(stat_CxInc_allPSA_list) = paste0("pct_", pctstat_vec)

	}
		
		MCMCout_estparam_temp_list = append(MCMCout_estparam_temp_list, list(CxInc=CxInc_allPSA_list, stat_CxInc=stat_CxInc_allPSA_list))
		MCMCout_estparam_temp_list_all[[ilist]] = MCMCout_estparam_temp_list;

		fnameout_estparm_xlsx = sprintf("%s%s_thin%d.xlsx", RDatafolder_input, output_suffix, mcmc_thinning)
		openxlsx::write.xlsx(MCMCout_estparam_temp_list, fnameout_estparm_xlsx, rowNames=TRUE)
	
	# extract parameter sets for each vec_HPV_set
	export_Extpara_temp = MCMCout_estparam_temp[ irow_pick_from_MCMC, ]
	
	stat_export_Extpara_temp = apply(export_Extpara_temp, 2, fun_pctstat)
	rownames(stat_export_Extpara_temp) = paste0("pct_", pctstat_vec)
	
	export_Extpara_temp = list(export_Extpara_temp, stat_export_Extpara_temp)
	names(export_Extpara_temp) = sprintf("%s_%s", names(export_Extpara)[iRData], c("param", "stat_param"))
	
	if (adjrelative_param_YN){
		adjust_export_Extpara_temp = adjust_MCMCout_estparam_temp[ irow_pick_from_MCMC, ]
		stat_adjust_export_Extpara_temp = apply(adjust_export_Extpara_temp, 2, fun_pctstat)
		rownames(stat_adjust_export_Extpara_temp) = paste0("pct_", pctstat_vec)

		adjust_export_Extpara_temp = list(adjust_export_Extpara_temp, stat_adjust_export_Extpara_temp)
		names(adjust_export_Extpara_temp) = sprintf("%s_%s", names(export_Extpara)[iRData], c("adjparam", "stat_adjparam"))
	
		export_Extpara_temp = append(export_Extpara_temp, adjust_export_Extpara_temp)
	}
	
	export_Extpara[[iRData]] = export_Extpara_temp;
	
} # for-ilist


# for some versions of the fit_CxInc_xxx.R.Version, no outcome stored in MCMCout_CxInc_list, so rerun CxInc
if (ncol(MCMCout_CxInc_list[[1]])==0){
	orgnames_MCMCout_CxInc_list = names(MCMCout_CxInc_list)
	MCMCout_CxInc_list = rep(list(NULL), length(iparamest_list_vec))
	MCMCout_CxInc_pRC_list = rep(list(NULL), length(iparamest_list_vec))
	names(MCMCout_CxInc_list) = orgnames_MCMCout_CxInc_list # retrieve name
	names(MCMCout_CxInc_pRC_list) = orgnames_MCMCout_CxInc_list
	
	for (ilist in iparamest_list_vec){
		timeXX1 = Sys.time()
		
		MCMCout_estparam_temp = MCMCout_estparam_list[[ilist]]
		HPVincid_ilist = HPVincid_tofit_list[[ilist]]
		CxInc_allPSA_list = rep(list(NULL), length=nrow(MCMCout_estparam_temp))
		
		for (iPSA in vec_HPV_set){
			if (iPSA%%10==1){ print(sprintf("iPSA = %d", iPSA))}
			HPVincid_iPSA = HPVincid_ilist[iPSA,]
			irow_iPSA = irow_eachMCMC[[iPSA]]
			CxInc_iPSA_list = rep(list(NULL), nrow_eachMCMC)
			
			for (irow in irow_iPSA){
				if (irow%%mcmc_thinning!=1){ next} # whether to re-calculate some to save time
				x_irow = MCMCout_estparam_temp[ irow, ]
				xdata_irow = list('HPVincid_data'=HPVincid_iPSA/100, 
					'outcome_agegp_diff_map'=outcome_agegp_diff_map, 
					'pop_agegp'=pop_agegp_temp,
					'paraest_setting'=paraest_setting,
					'mat_ratioAlive'=mat_ratioAlive_temp
					);

				extInput_irow = list('xdata'=xdata_irow, 'CxInc_data'=CxInc_datafit,
					'idx_CxInc_agegp_data'=idx_CxInc_outcomeAgegp_data, 'idx_CxInc_agegp_calib'=idx_CxInc_outcomeAgegp_calibrate,
					'logL_const'=logpoisspdf_C(CxInc_datafit), # log of the constant in poisson distribution
					'extInput_9vHR'=extInput_9vHR
					);
			
				CxInc_irow = fun_gen_CxInc_100k_byType(x_irow, extInput_irow, gen_pRC_YN=TRUE);

				CxInc_irow_temp = Reduce("+", CxInc_irow$CxInc) 
				CxInc_irow_pRC_temp = c(CxInc_irow$CxInc_pRC, "HPV_all9vHR"=Reduce("sum", CxInc_irow$CxInc[1:3])/Reduce("sum", CxInc_irow$CxInc))
				
				CxInc_iPSA_list[[irow]] = c(CxInc_irow_temp, CxInc_irow_pRC_temp)
			} # irow
			
			CxInc_allPSA_list[[iPSA]] = do.call(rbind, CxInc_iPSA_list)
			
		} # iPSA 
		
		# update MCMCout_CxInc_list
		CxInc_allPSA_temp = do.call(rbind, CxInc_allPSA_list)
		MCMCout_CxInc_list[[ilist]] = CxInc_allPSA_temp[, 1:outcome_agegp_num]
		MCMCout_CxInc_pRC_list[[ilist]] = CxInc_allPSA_temp[, -(1:outcome_agegp_num)]
		
		timeXX2 = Sys.time()
		print( timeXX2 - timeXX1 )
		
	} # ilist
} # check for recalculating MCMCout_CxInc


## plot
time1 = Sys.time()
pdf_fname = sprintf(fname_pdfplot_base, "plot_mcmcout")
pdf(pdf_fname, width=5, height=5)

# CxInc
par(mar=c(3,4,2,1))
plotx = seq(10, 85, by=5)
xlim_plot = range(plotx) # c(10, 85)
range_plot = 90;

if (any(1:3 %in% idx_CxInc_outcomeAgegp_data)){
	idx_plot_CxInc = 4:16
} else{
	idx_plot_CxInc = 1:length(idx_CxInc_outcomeAgegp_data)
}


for (ilist in iparamest_list_vec){
	CxInc_estparam_temp = MCMCout_CxInc_list[[ilist]][, idx_CxInc_outcomeAgegp_data] # check if needed
	
	if (mcmc_thin_YN){ # thinning of MCMC outputs
		if (nrow(CxInc_estparam_temp) == nrow(MCMCout_estparam_list[[ilist]])){ # if not equal, then MCMCout_CxInc_list may be re-calcuated and no need to thin
			CxInc_estparam_temp = CxInc_estparam_temp[mcmc_thin_vec, ]
		}
	}

	plotx_CxInc = HPVincid_agegp_new[idx_CxInc_outcomeAgegp_data];
	ploty_CxInc = t( apply(CxInc_estparam_temp,2, quantile, probs=c(0.5, 0.5+range_plot/100/2*c(-1,1))) )

	CxInc_datafit_temp = CxInc_datafit_list[[ilist]];
	ploty_CxInc_datafit = t( CxInc_datafit_temp ); 
	
	plotx_CxInc = plotx_CxInc[idx_plot_CxInc]
	ploty_CxInc = ploty_CxInc[idx_plot_CxInc,]
	ploty_CxInc_datafit = ploty_CxInc_datafit[idx_plot_CxInc,]

	ylim_plot = c(0, ceiling(max(max(ploty_CxInc), max(CxInc_datafit_temp))))

	plot(NA, xlim=xlim_plot, ylim=ylim_plot, axes=FALSE, xlab="", ylab="Cancer incidence (per 100,000)")
	mtext(side=1, text='Age (years)', line=2)
	axis(1, at=plotx)
	axis(2, las=1)
	
	matlines(x=plotx_CxInc, y=ploty_CxInc_datafit, col='blue', lwd=c(2,1.5,1.5), lty=c(1,2,2))
	matlines(x=plotx_CxInc, y=ploty_CxInc, col='orange', lwd=2, lty=c(1,2,2))

	mtext(text=names(MCMCout_CxInc_list)[ilist], side=3, line=0.5)

	legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot), 'Empirical (10-year average)', 'Empirical (10-year lowest/highest)'), col=c(rep('orange',2),rep('blue',2)), lty=rep(c(1,2),2), lwd=2, cex=0.8)

	rm(CxInc_estparam_temp)
	gc(reset=TRUE)


} # for-ilist, CxInc


# scaling_factor
if (n_scaling!=1){
par(mar=c(3,5,2,1))
for (ilist in iparamest_list_vec){

	MCMCout_estparam_temp = MCMCout_estparam_list[[ilist]]
	if (mcmc_thin_YN){ # thinning of MCMC outputs
		MCMCout_estparam_temp = MCMCout_estparam_temp[mcmc_thin_vec, ]
	}
	
	for (iHPV in 1:n_set){
	ploty = apply(MCMCout_estparam_temp[, nparam_iset*(iHPV-1) + (2*n_delay+(1:n_scaling))], 1, fun_beta_scaling_map);


	plotx = 1:16;
	plotx_at = seq(1, 16, by=2);
	plotx_lab = seq(10, 80, by=10);
	
	ploty_median = t( apply(ploty,1, quantile, c(0.5,0.05,0.95)) )
	ylim_plot = range(pretty(c(0, max(ploty_median))))
	
	# matplot(x=plotx, y=ploty, type='l', xlab='', ylab='', col='grey80', las=1, axes=FALSE, ylim=ylim_plot)
	# matlines(plotx, ploty_median, col='blue', lwd=2, lty=c(1,2,2))
	matplot(plotx, ploty_median, col='blue', lwd=2, lty=c(1,2,2), type='l', xlab='', ylab='', axes=FALSE, ylim=ylim_plot)
	axis(1, at=plotx_at, labels=plotx_lab)
	axis(2, las=1)
	mtext(side=1, text='Age (years)', line=2)
	mtext(side=2, text='scaling factor', line=4);
	mtext(side=3, text=sprintf("HPV%d_%s", iHPV, gsub("MCMCout_param_","scaling_", names(MCMCout_estparam_list)[ilist])), line=0.5)
	legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot)), col=rep('blue',2), lty=c(1,2), lwd=2, cex=0.8)
	
	} # for-iset

	rm(MCMCout_estparam_temp, ploty)
	gc(reset=TRUE)


} # end for-ilist, scaling_map
} # if n_scaling!=1


# traceplot
par(mar=c(3,3,1,1))
for (ilist in iparamest_list_vec){
	MCMCout_estparam_temp = coda::as.mcmc(MCMCout_estparam_list[[ilist]])
	if (mcmc_thin_YN){
		MCMCout_estparam_temp = coda::as.mcmc(MCMCout_estparam_temp[mcmc_thin_vec,])
	}
	
	layout(matrix(1:2, nrow=2, byrow=TRUE))
	nullout = sapply(1:ncol(MCMCout_estparam_temp), function(x) plot(MCMCout_estparam_temp[,x], auto.layout=FALSE, main=varname_vec[x]))
	
	rm(MCMCout_estparam_temp)
	gc(reset=TRUE)
} # end for-ilist, traceplot



# 2D histogram
for (ilist in iparamest_list_vec){

	MCMCout_estparam_temp = MCMCout_estparam_list[[ilist]]
	if (mcmc_thin_YN){
		MCMCout_estparam_temp = MCMCout_estparam_temp[mcmc_thin_vec,]
	}
	MCMCout_estparam_temp = as.data.frame(MCMCout_estparam_temp)
	colnames(MCMCout_estparam_temp) = varname_vec

	for (iHPV in 1:3){
		hout = hist2d_gg(MCMCout_estparam_temp, sprintf("HPV_%d_%s", iHPV, c("gamma_shape","gamma_scale")))
		print( hout )
		
		hout = hist2d_gg(MCMCout_estparam_temp, sprintf("HPV_%d_%s", iHPV, c("scaling1","scaling2")))
		print( hout )
	}

	rm(MCMCout_estparam_temp)
	gc(reset=TRUE)
} # end for-ilist, 2D histogram



if (FALSE){
	# correlation of parameters
	# http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs
	par(mar=c(3,3,3,3)) # does not affect the pattern, but will affect the line for mtext 
	for (ilist in iparamest_list_vec){

		MCMCout_estparam_temp = MCMCout_estparam_list[[ilist]]
		if (TRUE){ # thinning of MCMC outputs
			# need thinning otherwise will use up many space
			MCMCout_estparam_temp = MCMCout_estparam_temp[mcmc_thin_vec, ]
		}

		psych::pairs.panels(MCMCout_estparam_temp, method = "pearson", hist.col = "#00AFBB", density = TRUE, ellipses = TRUE)
		mtext(side=3, text=gsub("MCMCout_param_","_correlation_", names(MCMCout_estparam_list)[ilist]), line=1.5)

		rm(MCMCout_estparam_temp)
		gc(reset=TRUE)
		
	} # end-ilist, correlation
}


# when all plots completed
dev.off()
time2 = Sys.time()
timediff = time2 - time1
cat( sprintf("Time to plot: %f %s\n", as.numeric(timediff), attributes(timediff)$units) )



} # end for-iRData


# one file to plot
if (iparamest_list_vec==1){
	print( "plot one CxInc plot" )

	pdf_fname_CxInc_byCx = sprintf("%s%s_%s.pdf", RDatafolder_input, "plot_CxInc_allCx", outputdate_string)
	pdf(pdf_fname_CxInc_byCx, width=6*2, height=4*ceiling(length(RData_vec)/2))

	layout(matrix(1:length(RData_vec), ncol=2, byrow=FALSE))
	par(mar=c(3.5,4,2.5,1))

	iplottype_vec = 1:2
	for (iplottype in iplottype_vec){

		for (iRData in 1:length(RData_vec)){

		RDataname_input = RData_vec[iRData]
		load( gsub("//","/", paste0(RDatafolder_input, RDataname_input)) )
		
			# plotx may be included in the loading .RData
			plotx = seq(25, 85, by=5)
			xlim_plot = range(plotx) # c(10, 85)
			range_plot = 90;
			col_empirical = 'blue'
			col_model = 'orange' 
			pch_empirical = 0
			pch_model = 2

			
			if (any(1:3 %in% idx_CxInc_outcomeAgegp_data)){
				idx_plot_CxInc = 4:16
			} else{
				idx_plot_CxInc = 1:length(idx_CxInc_outcomeAgegp_data)
			}			

		for (ilist in iparamest_list_vec){
			CxInc_estparam_temp = MCMCout_CxInc_list[[ilist]][, idx_CxInc_outcomeAgegp_data] # check if needed
			
			if (mcmc_thin_YN){ # thinning of MCMC outputs
				if (nrow(CxInc_estparam_temp) == nrow(MCMCout_estparam_list[[ilist]])){ # if not equal, then MCMCout_CxInc_list may be re-calcuated and no need to thin
					CxInc_estparam_temp = CxInc_estparam_temp[mcmc_thin_vec, ]
				}
			}

			plotx_CxInc = HPVincid_agegp_new[idx_CxInc_outcomeAgegp_data];
			ploty_CxInc = t( apply(CxInc_estparam_temp,2, quantile, probs=c(0.5, 0.5+range_plot/100/2*c(-1,1))) )
			
			CxInc_datafit_temp = CxInc_datafit_list[[ilist]];
			ploty_CxInc_datafit = t( CxInc_datafit_temp ); 
			
			# to plot ages 25+ only
			plotx_CxInc = plotx_CxInc[idx_plot_CxInc]
			ploty_CxInc = ploty_CxInc[idx_plot_CxInc,]
			ploty_CxInc_datafit = ploty_CxInc_datafit[idx_plot_CxInc,]

			ylim_plot = c(0, ceiling(max(max(ploty_CxInc), max(CxInc_datafit_temp))))

			plot(NA, xlim=xlim_plot, ylim=ylim_plot, axes=FALSE, xlab="", ylab="")
			mtext(side=1, text='Age (years)', line=2.25, cex=1.2)
			mtext(side=2, text="Cancer incidence (per 100,000)", line=2.75)
			mtext(text=names(MCMCout_CxInc_list)[ilist], side=3, line=0.5)
			axis(1, at=plotx, cex.axis=1.2)
			axis(2, las=1, cex.axis=1.2)
			
			if (iplottype==1){			
				# line plot
			
				matlines(x=plotx_CxInc, y=ploty_CxInc_datafit, col=col_empirical, lwd=c(2,1.5,1.5), lty=c(1,2,2))
				matlines(x=plotx_CxInc, y=ploty_CxInc, col=col_model, lwd=2, lty=c(1,2,2))

				legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot), 'Empirical (10-year average)', 'Empirical (10-year lowest/highest)'), col=c(rep(col_model,2),rep(col_empirical,2)), lty=rep(c(1,2),2), lwd=1.5, cex=1, inset=c(0.025,0))
			
			} else if (iplottype==2){
				# marker, error bar plot
			
				xshift = 0.75
				matpoints(x=plotx_CxInc-xshift, y=ploty_CxInc_datafit[,1], col=col_empirical, lwd=1.5, pch=pch_empirical)
				segments(x0=plotx_CxInc-xshift, y0=ploty_CxInc_datafit[,2], y1=ploty_CxInc_datafit[,3], col=col_empirical, lwd=1.5)
				
				matpoints(x=plotx_CxInc+xshift, y=ploty_CxInc[,1], col=col_model, lwd=1.5, pch=pch_model)
				segments(x0=plotx_CxInc+xshift, y0=ploty_CxInc[,2], y1=ploty_CxInc[,3], col=col_model, lwd=1.5)


				legend('topleft', legend=c('Median of estimates', sprintf('%d%% range of estimates', range_plot), 'Empirical (10-year average)', 'Empirical (10-year lowest/highest)'), col=c(rep(col_model,2),rep(col_empirical,2)), lty=c(NA,1,NA,1), pch=c(pch_model, NA, pch_empirical, NA), lwd=1, cex=1, inset=c(0.025,0))
			}

			rm(CxInc_estparam_temp)
			gc(reset=TRUE)

		} # for-ilist, CxInc

		} # end for-iRData

	} # for-iplottype


	dev.off()

} # if iparamest_list_vec==1

export_Extpara = do.call(c, export_Extpara)
names(export_Extpara) = gsub("^(.*)\\.","", names(export_Extpara))

# export the parameter sets for PSA
openxlsx::write.xlsx(export_Extpara, file=fname_export_Extpara, rowNames=TRUE)

