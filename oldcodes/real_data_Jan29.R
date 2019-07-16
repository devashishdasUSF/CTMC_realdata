library("Rcpp")
library("Matrix")
library("dplyr")
library("magrittr")
library("stringr")
library("broom")

folder = "/shares/ddasrg/CTMC_control/real_data/CTMC_real_data/"
sourceCpp(paste0(folder, "real_data.cpp"))


# ==========================================================================
# DATA FOLER "/shares/ddasrg/CTMC_control/real_data/dataJan31
# ==========================================================================
dat_f = "/shares/ddasrg/CTMC_control/real_data/dataJan31/"
load(paste0(dat_f, "split.RData"))
load(paste0(dat_f, "info.RData"))
load(paste0(dat_f, "transition_types.RData"))

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


# ==========================================================================
# old selection criteria
# ==========================================================================
# selection_criteria = #(info$LOS <= 480.00) & 
# 	# (info$ESILevel == "ESI Level 3") &
# 	(info$PresentingComplaint == "Chest Pain") &  
# 	(info$Disposition != "Expired")

# train = which(selection_criteria & (info$Camebackin14days == FALSE)) %>% unname

# test = which(selection_criteria & (info$Camebackin14days == TRUE)) %>% unname

##
# 
#
# transition_types = c("arr-BedDate", "BedDate-FirstProv", "BedDate-EKG")
# dat = splitdata
train = (info$PresentingComplaint == "Chest Pain" & (info$ESILevel == "ESILevel 1" | info$ESILevel == "ESILevel 2") & month <= 6)  & (info$LOS <= 1500.00) 
test = (info$PresentingComplaint == "Chest Pain") & (info$LOS <= 1500.00)

which(train) %>% length
which(test )%>% length

dat = splitdata[train]
testdat = splitdata[test]
N1 = testdat %>% length
N = dat%>% length
Nbasis = 1000

# transitions_types %<>% filter(tran %in% selec_tran)
PreviousState = transition_types %>% sapply(function(x)(strsplit(x, split="-"))) %>% lapply(first) %>% unlist %>% unname
M = length(transition_types)

PenMatdiag = bandSparse(Nbasis,k=c(-1,0,1), diagonals = list(rep(-1,Nbasis-1), c(1,rep(2,Nbasis-2),1),rep(-1,Nbasis-1)))
PenMat = Diagonal(Nbasis)

coefficients = list()
for(I in seq(M)) {
	X = get_n_y(dat, PreviousState[I], transition_types[I], Nbasis)
	y = X$g_ij
	nij = X$n_ij
	G = Diagonal(x = y)
	beta_ij = solve(G + N*(0.0001*PenMat+0.001*PenMatdiag), nij)
	coefficients[[I]] = beta_ij[,1]
	# filename = paste0("AllESI_1_2/",transition_types[I],".RData")
	# save(beta_ij, file = filename)
}
names(coefficients) = transition_types






stat0 = numeric()
for(i in seq(length(dat))) {
	temp = rep(0,M)
	for(I in seq(M)) {
		temp[I] = generalized_test(dat[[i]], coefficients[[I]],PreviousState[I], transition_types[I], 0.0001, Nbasis)
	}
	stat0[i] = max(temp)
}
reject = stat0>=quantile(stat0,.95)
stat0_new = stat0[!reject]



stat1 = numeric()
for(i in seq(length(testdat))) {
	temp = rep(0,M)
	for(I in seq(M)) {
		temp[I] = generalized_test(testdat[[i]], coefficients[[I]],PreviousState[I], transition_types[I], 0.001, Nbasis)
	}
	stat1[i] = max(temp)
}


print("---- Beta error FR ----")
print(sum(stat1 <= quantile(stat0_new,.9))/(length(stat1)))



stat0_los = numeric()
for(i in seq(length(dat))) {
	x = dat[[i]] %>% filter(Event == "Departure")
	stat0_los[i] = x$Time
}

stat1_los = numeric()
for(i in seq(length(testdat))) {
	x = testdat[[i]] %>% filter(Event == "Departure")
	stat1_los[i] = x$Time
}

print("---- Beta error simple case----")
print(sum(stat1_los <= quantile(stat0_los[stat0_los<quantile(stat0_los,.9)],.9))/(length(stat1_los)))






#  1 Abdominal Pain                  4900
#  2 Chest Pain                      3723
#  3 Shortness of Breath             1789
#  4 Fall                            1766
#  5 Fever                           1536
#  6 Dyspnea                         1143
#  7 Headache                         915
#  8 Dizziness                        909
#  9 (BH) Suicidal Ideation/Attempt   876
# 10 Cough (masked)                   824




