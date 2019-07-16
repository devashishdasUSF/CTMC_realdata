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






train1 = (

	info$PresentingComplaint == "Chest Pain" & 
	(info$ESILevel == "ESI Level 1" | info$ESILevel == "ESI Level 2") 
	& info$Month <= 6 
	& info$Camebackin14days == TRUE
	
	)  & (info$LOS <= 1500.00)

train2 = (

	info$PresentingComplaint == "Chest Pain" & 
	(info$ESILevel == "ESI Level 1" | info$ESILevel == "ESI Level 2") 
	& info$Month <= 6 
	& info$Camebackin14days == FALSE
	
	)  & (info$LOS <= 1500.00)




transition_types = c("arr-BedDate")#, "BedDate-FirstProv", "BedDate-EKG")
PreviousState = transition_types %>% sapply(function(x)(strsplit(x, split="-"))) %>% lapply(first) %>% unlist %>% unname

PenMatdiag = bandSparse(Nbasis,k=c(-1,0,1), diagonals = list(rep(-1,Nbasis-1), c(1,rep(2,Nbasis-2),1),rep(-1,Nbasis-1)))
PenMat = Diagonal(Nbasis)

M = length(transition_types)

### average cumulative intensity functions

dat = splitdata[train2]
### jitter data
jitter_times = function(df) {
	noise = sort(abs(rnorm(nrow(df))/1000))
	df$Time = df$Time + noise
	df
}
dat = lapply(dat,jitter_times)
N = dat%>% length
Nbasis = 1000

I = 1
X = get_n_y(dat, PreviousState[I], transition_types[I], Nbasis, 100)
y = X$g_ij
nij = X$n_ij
G = Diagonal(x = y)
beta_ij = solve(G + N*(0.0001*PenMat), nij)


plot(beta_ij, type = "l", col = "green")
lines((nij/y))

plot((nij/y), type = "l", xlim = c(0,100), ylim = c(0,45))
lines(beta_ij, , col = "green")

coefficients = list()
for(I in seq(M)) {
	X = get_n_y(dat, PreviousState[I], transition_types[I], Nbasis)
	y = X$g_ij
	nij = X$n_ij
	G = Diagonal(x = y)
	beta_ij = solve(G + N*(0.001*PenMatdiag), nij)
	coefficients[[I]] = beta_ij[,1]
	# filename = paste0("AllESI_1_2/",transition_types[I],".RData")
	# save(beta_ij, file = filename)
}
names(coefficients) = transition_types

plot(coefficients[[1]], type = "l")







