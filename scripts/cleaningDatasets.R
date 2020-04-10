# SET UP ------------------------------------------------------------------

library(tidyverse)
library(stringr)
library(data.table)
library(traitdataform)



# DATA CLEANING -----------------------------------------------------------


# * Read in raw files -------------------------------------------------------


iucnSp <- read_csv('data/finalIUCNsp.csv') %>% 
  mutate(state = tolower(state)) %>% 
  distinct(state, species, .keep_all = TRUE) %>% 
  filter(!is.na(state)) %>%
  mutate(
    state = case_when(
      state == 'new mex' ~ 'new mexico',
      state == 'tennesee' ~ 'tennessee',
      state == 'louisianna' ~ 'louisiana',
      state == 'massachusettes' ~ 'massachusetts',
      state == 'virgina' ~ 'virginia',
      TRUE ~ state
    )
  ) %>% 
  filter(!(state %in% c('all oceanic states', 'navajo nation', 'district of columbia'))) %>% 
  mutate(species = word(species, 1, 2, sep = " "))



stateAll <- read_csv('data/stateListedSpecies.csv') %>% 
  select(state:Notes) %>% 
  mutate(state = tolower(state)) %>% 
  distinct(state,species, .keep_all = TRUE) %>% 
  filter(state != 'NA')


stateImp <- stateAll %>% 
  filter(!(Ranking %in% c('PRNG', 'X')), !(state == 'new york' & Ranking == 'SC'),
         !(state %in% c('utah', 'north dakota')),
         !(state == 'washington' & Ranking %in% c('T', 'S'))) %>% 
  mutate(species = word(species, 1, 2, sep = " ")) %>% #subset to just sp. name (no subsp)
  distinct(state, species, .keep_all = TRUE) 


fedSp <- read_csv('data/fedListedSp.csv') %>% 
  select(state:recorder) %>% 
  mutate(state = tolower(state)) %>% 
  filter(state != 'NA', state != 'district of columbia') %>% 
  mutate(species = word(species, 1, 2, sep = " ")) #just keep genus and species



# * get accepted names so datasets match -------------------------------------------------

## use taxize and traitdataform to get consistent list of species names



## reduce lists to just unique species and remove those just to genus level

stateList <- stateImp %>% distinct(species) %>% 
  filter(!(str_detect(species, '\\s(sp).$')), !(str_detect(species, '\\s(sp)$')))

fedList <- fedSp %>% distinct(species) %>% 
  filter(!(str_detect(species, '\\s(sp).$')), !(str_detect(species, '\\s(sp)$')))

iucnList <- iucnSp %>% distinct(species) %>% 
  filter(!(str_detect(species, '\\s(sp).$')), !(str_detect(species, '\\s(sp)$')))

# do for full U.S. Red list too
redList <- read_csv('data/IUCN_US_Species.csv') %>% rename(species = result.scientific_name) %>% 
  mutate(species = word(species, 1, 2, sep = " "))%>% 
  distinct(species) %>% 
  filter(!(str_detect(species, '\\s(sp).$')), !(str_detect(species, '\\s(sp)$')))

redListFull <- read_csv('data/IUCN_US_Species.csv') %>% rename(species = result.scientific_name) %>% 
  mutate(species = word(species, 1, 2, sep = " "))%>% 
  distinct(species, .keep_all = TRUE) %>% select(country:result.category)

# get list of accepted names for each group

## state list; some species names not found in database, so need to run it
## this way to continue search function through the errors
sp2 <- vector("list", length = nrow(stateList))

for (i in 1:nrow(stateList)){
  tryCatch({
    sp2[[i]] <- get_gbif_taxonomy(
      paste0(stateList[i, ]),
      subspecies = FALSE,
      higherrank = FALSE,
      fuzzy = TRUE,
      verbose = TRUE,
      conf_threshold = 10
    )
    print(i)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}


stateListClean <- rbindlist(sp2, fill = TRUE)
#write.csv(stateListClean, 'data/stateSpeciesNames.csv')


## fed list; all species found in database
fedListClean <- fedList %>% pull(species) %>% 
  get_gbif_taxonomy(., verbose = TRUE)
# write.csv(fedListClean,'data/federalSpeciesNames.csv')


## iucn list
sp <- vector("list", length = nrow(iucnList))

for (i in 1:nrow(iucnList)){
  tryCatch({
    sp[[i]] <- get_gbif_taxonomy(
      paste0(iucnList[i, ]),
      subspecies = FALSE,
      higherrank = FALSE,
      fuzzy = TRUE,
      verbose = TRUE,
      conf_threshold = 10
    )
    print(i)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}


iucnListClean <- rbindlist(sp, fill = TRUE)
#write.csv(iucnListClean, 'data/IUCNSpeciesNames.csv')


## full U.S. IUCN Red List (not just imperiled list)
sp <- vector("list", length = nrow(redList))

for (i in 1:nrow(redList)){
  tryCatch({
    sp[[i]] <- get_gbif_taxonomy(
      paste0(redList[i, ]),
      subspecies = FALSE,
      higherrank = FALSE,
      fuzzy = TRUE,
      verbose = TRUE,
      conf_threshold = 10
    )
    print(i)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}


redListClean <- rbindlist(sp, fill = TRUE) %>% as_tibble()
redListClean$scientificName <- as.character(redListClean$scientificName)
redListClean$scientificNameStd <- as.character(redListClean$scientificNameStd)


redListFull$acceptedName <- vector(mode = "character", length = nrow(redListFull))

for (i in 1:nrow(redListFull)){
  redListFull$acceptedName[i] <-  
    ifelse(redListFull[i,3] %in% redListClean$scientificName, redListClean$scientificNameStd[redListClean$scientificName == paste(redListFull[i,3])],
           paste(redListFull[i,3]))
}

#write.csv(redListFull, 'data/US_IUCN_acceptedName.csv')



# * clean up final lists of species -----------------------------------------

stateNames <- read_csv('data/stateSpeciesNames.csv')
fedNames <- read_csv('data/federalSpeciesNames.csv')
iucnNames <- read_csv('data/iucnSpeciesNames.csv')

# need to add new species names to original species lists, for each sp name,
# match with sp list, if syn = true then make name accepted name, if NA
# then keep as current name, if false then keep accepted name

# cleaning before accidentally removed all species names that also started with 'sp...',
# Clearly my regex skills need some work

# need to filter those out, redo search, and then rbind back to state Names

missingState <- stateImp %>% filter(!(species %in% stateNames$scientificName)) %>%
  filter(!(str_detect(species, '\\s(sp).$')), !(str_detect(species, '\\s(sp)$')),
         !(str_detect(species, '\\s(ssp)')), !(str_detect(species, '\\s(var).$')),
         !(str_detect(species, '\\s(sp.s)')), !(str_detect(species, '\\s(spp.)')),
         !(str_detect(species, '\\s(species)'))) %>% distinct(species)

## re-do state species list for missing state species
sp2 <- vector("list", length = nrow(missingState))

for (i in 1:nrow(missingState)){
  tryCatch({
    sp2[[i]] <- get_gbif_taxonomy(
      paste0(missingState[i, ]),
      subspecies = FALSE,
      higherrank = FALSE,
      fuzzy = TRUE,
      verbose = TRUE,
      conf_threshold = 10
    )
    print(i)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}


missingNames <- rbindlist(sp2, fill = TRUE) %>% as_tibble() %>% 
  mutate(species = scientificName)

stateNamesNew <- stateNames %>% select(scientificName:warnings) %>% 
  mutate(species = scientificName) %>% rbind(missingNames)
#write.csv(stateNamesNew, 'data/stateNamesUpdated.csv')

# now add new column in stateImp with the accepted name

stateImp$acceptedName <- vector(mode = "character", length = nrow(stateImp))

for (i in 1:nrow(stateImp)){
  stateImp$acceptedName[i] <-  
    ifelse(stateImp[i,2] %in% stateNamesNew$scientificName, stateNamesNew$scientificNameStd[stateNamesNew$scientificName == paste(stateImp[i,2])],
           paste(stateImp[i,2]))
}


#now do same as above for federal and iucn

##federal fixes
missingFed <- fedSp %>% filter(!(species %in% fedNames$scientificName)) %>%
  filter(!(str_detect(species, '\\s(sp.)$')),
         !(str_detect(species, '\\s(spp.)'))) %>% 
  pull(species) %>% 
  get_gbif_taxonomy(., verbose = TRUE)

fedNamesNew <- fedNames %>% select(scientificName:warnings) %>% 
  rbind(missingFed) 
#write.csv(fedNamesNew, 'data/fedNamesUpdated.csv')

fedSp$acceptedName <- vector(mode = "character", length = nrow(fedSp))

for (i in 1:nrow(fedSp)){
  fedSp$acceptedName[i] <-  
    ifelse(fedSp[i,2] %in% fedNamesNew$scientificName, fedNamesNew$scientificNameStd[fedNamesNew$scientificName == paste(fedSp[i,2])],
           paste(fedSp[i,2]))
}


## iucn
missingIUCN <- iucnSp %>% filter(!(species %in% iucnNames$scientificName)) %>% 
  filter(!(str_detect(species, '\\s(sp.)$'))) %>% distinct(species)

sp <- vector("list", length = nrow(missingIUCN))

for (i in 1:nrow(missingIUCN)){
  tryCatch({
    sp[[i]] <- get_gbif_taxonomy(
      paste0(missingIUCN[i, ]),
      subspecies = FALSE,
      higherrank = FALSE,
      fuzzy = TRUE,
      verbose = TRUE,
      conf_threshold = 10
    )
    print(i)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}



missingIUCN <- rbindlist(sp, fill = TRUE) %>% as_tibble()


iucnNamesNew <- iucnNames %>% select(scientificName:warnings) %>% 
  rbind(missingIUCN)
#write.csv(iucnNamesNew, 'data/iucnNamesUpdated.csv')

iucnSp$acceptedName <- vector(mode = "character", length = nrow(iucnSp))

for (i in 1:nrow(iucnSp)){
  iucnSp$acceptedName[i] <-  
    ifelse(iucnSp[i,1] %in% iucnNamesNew$scientificName, iucnNamesNew$scientificNameStd[iucnNamesNew$scientificName == paste(iucnSp[i,1])],
           paste(iucnSp[i,1]))
}

###There are NA values in accepted name where there is NA for accepted name, some look funky
## so filter out NA, rename as 'acceptedName', then rbind back with full dataset

#state NA species

stateNA <- stateImp %>% filter(is.na(acceptedName)) %>% select(state:Notes) %>% 
  mutate(acceptedName = species) %>% mutate(acceptedName = str_remove(acceptedName, "\\?")) %>% 
  mutate(acceptedName = str_remove(acceptedName, "\\.")) %>% 
  filter(!(str_detect(species, '\\s(new)$')), !(str_detect(species, 'Turtle')),
         !(str_detect(species, '\\s(d.)$')), !(str_detect(species, '\\s(m.)$')),
         !(str_detect(species, '\\[')))

stateImp <- stateImp %>% filter(!is.na(acceptedName)) %>% rbind(stateNA)

# fed NA

fedNA <- fedSp %>% filter(is.na(acceptedName)) %>% select(state:recorder) %>% 
  mutate(acceptedName = species)

fedSp <- fedSp %>% filter(!is.na(acceptedName)) %>% rbind(fedNA)

#iucn NA

iucnNA <- iucnSp %>% filter(is.na(acceptedName)) %>% select(species:state) %>% 
  mutate(acceptedName = species)

iucnSp <- iucnSp %>% filter(!is.na(acceptedName)) %>% rbind(iucnNA)

##SAVED FINAL FILES AS '_ACCEPTEDNAME'

### final fixes

stateImp <- read_csv('data/stateSpeciesAcceptedName.csv') %>% select(-X1)

fedSp <- read_csv('data/fedSpeciesAcceptedName.csv') %>% select(-X1)

iucnSp <- read_csv('data/iucnSpeciesAcceptedName.csv') %>% select(-X1)

fedSp <- fedSp %>% mutate(state = if_else(state == 'conneticut', 'connecticut', state)) %>% 
  distinct(state, acceptedName, .keep_all = TRUE)

# add TN federally-listed species to stateImp file; TN lists them but they were not on the original
# list for some reason
tnFed <- fedSp %>% filter(state == 'tennessee') %>% mutate(Ranking = status, Notes = 'NA') %>% 
  select(state:taxon_specific, Ranking, acceptedName, Notes)

stateImp <- stateImp %>% rbind(tnFed) %>% distinct(state, acceptedName, .keep_all = TRUE)

#all files saved as iucnFINAL, fedFINAL, and stateFINAL for analyses

