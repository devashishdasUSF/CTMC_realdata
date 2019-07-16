library("magrittr")

folder = "/shares/ddasrg/CTMC_control/real_data/CTMC_real_data/"
rawdata = read.table(file = paste0("deidentified.txt"), header = T)

PresentingComplaint = (rawdata %>% group_by(ind) %>% summarize(complaint = first(PresentingComplaint)))$complaint
ESILevel = (rawdata %>% group_by(ind) %>% summarize(esi = first(ESILevel)))$esi
Month = (rawdata %>% group_by(ind) %>% summarize(month = first(month)))$month
Camebackin14days = (rawdata %>% group_by(ind) %>% summarize(testtrain = first(Camebackin7days)))$testtrain
patient_id = (rawdata %>% group_by(ind) %>% summarize(id = first(ind)))$id
Disposition = (rawdata %>% group_by(ind) %>% summarize(Disposition = first(DispoType)))$Disposition
LOS = (rawdata %>% group_by(ind) %>% summarize(LOS= last(Time)))$LOS

info = data.frame(patient_id, ESILevel, Month, LOS, Disposition, PresentingComplaint, Camebackin14days)

splitdata = split(rawdata, rawdata$ind) %>%
	lapply(function(x){x %>% select(Event, Time, transition, previous) %>%
		mutate(tran = transition, previous_state = previous) %>% 
		select(-transition, -previous)
})


transition_types = unique(rawdata$transition)

save(splitdata, file="split.RData")
save(info, file = "info.RData")
save(transition_types, file="transition_types.RData")