library("Rcpp")
library("Matrix")
library("dplyr")
library("magrittr")
library("stringr")
library("broom")
load("splitdata.RData")
load("transitions_types.RData")
sourceCpp("real_data.cpp")
# ==========================================================================
# Description of variables
# ==========================================================================
# individual_presenting_complaint - Patients presenting complaint
# individual_cameback - whether the patient returned in 7 days or note
# ESILevel - Severity index. 1 is most severe, and 5 is least severe
# patient_id - scrambled IBEX numbers
# Disposition - Final status of the patient
# Max_time - Last time point recorded
# alldata - data frame with (Event, previous_state, tran, time) info
# ==========================================================================
# save(alldata, individual_presenting_complaint, individual_cameback, patient_id, ESILevel, Disposition, Max_time, file = "splitdata.RData")
# selec_tran = c(
# "BedDate-FirstProv",  
# "arr-BedDate",        
# "Radiology-Departure",
# "FirstProv-Lab",      
# "Lab-Radiology",      
# "Lab-Nurse",          
# "Lab-Departure",      
# "EKG-Lab",            
# "Nurse-Radiology",    
# "FirstProv-EEG")
# # selec_tran = "arr-BedDate"

selection_criteria = (Max_time <= 480.00) &
	(individual_presenting_complaint == "Chest Pain") &  
	(Disposition != "Expired")

idx = which(selection_criteria & (individual_cameback == 0)) %>% unname

test_idx = which(selection_criteria & (individual_cameback == 1)) %>% unname


idx %>% length
dat = alldata[idx]
N = dat%>% length
Nbasis = 1000

# transitions_types %<>% filter(tran %in% selec_tran)
M = nrow(transitions_types)

PenMatdiag = bandSparse(Nbasis,k=c(-1,0,1), diagonals = list(rep(-1,Nbasis-1), c(1,rep(2,Nbasis-2),1),rep(-1,Nbasis-1)))
PenMat = Diagonal(Nbasis)

coefficients = list()
for(I in seq(M)) {
	X = get_n_y(dat, transitions_types$previous_state[I], transitions_types$tran[I], Nbasis)
	y = X$g_ij
	nij = X$n_ij
	G = Diagonal(x = y)
	beta_ij = solve(G + N*(0.01*PenMat+10*PenMatdiag), nij)
	coefficients[[I]] = beta_ij[,1]
	# filename = paste0("AllESI_1_2/",transitions_types$tran[I],".RData")
	# save(beta_ij, file = filename)
}
names(coefficients) = transitions_types$tran

test_dat = alldata[test_idx]

stat0 = numeric()
for(i in seq(length(dat))) {
	temp = rep(0,M)
	for(I in seq(M)) {
		temp[I] = generalized_test(dat[[i]], coefficients[[I]],transitions_types$previous_state[I], transitions_types$tran[I], 0.01, Nbasis)
	}
	stat0[i] = max(temp)
}
stat1 = numeric()
for(i in seq(length(test_dat))) {
	temp = rep(0,M)
	for(I in seq(M)) {
		temp[I] = generalized_test(test_dat[[i]], coefficients[[I]],transitions_types$previous_state[I], transitions_types$tran[I], 0.01, Nbasis)
	}
	stat1[i] = max(temp)
}


print("---- Beta error FR ----")
print(sum(stat1 <= quantile(stat0,.9))/(length(stat1)))


stat0_los = numeric()
for(i in seq(length(dat))) {
	x = dat[[i]] %>% filter(Event == "BedDate")
	stat0_los[i] = x$Time
}

stat1_los = numeric()
for(i in seq(length(test_dat))) {
	x = dat[[i]] %>% filter(Event == "BedDate")
	stat1_los[i] = x$Time
}

print("---- Beta error simple case----")
print(sum(stat1_los <= quantile(stat0_los[stat0_los<quantile(stat0_los,.9)],.9))/(length(stat1_los)))











