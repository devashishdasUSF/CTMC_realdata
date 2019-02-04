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
train = (info$PresentingComplaint == "Abdominal Pain")  & (info$LOS <= 1500.00) & (info$Camebackin14days == FALSE)
test = (info$PresentingComplaint == "Abdominal Pain") & (info$LOS <= 1500.00) & (info$Camebackin14days == TRUE)

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


