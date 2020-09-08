
# SET UP ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(stringr)
library(data.table)
library(ggpubr)
library(Hmisc)
library(MASS)
library(DescTools)
library(fiftystater)


## download state shapefile for maps
data("fifty_states")


## read in most recent files
iucnSp <- read_csv('data/iucnFINAL_revisions.csv')
fedSp <- read_csv('data/fedFINAL.csv')
stateImp <- read_csv('data/stateFINAL.csv')

# make summary files (# sp per state)
stateAllSp <- stateImp %>% add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% rename(totalProtectedSpecies = n) %>% ungroup() %>% 
  mutate(totalProtectedSpecies = if_else(totalProtectedSpecies == 1, NA_integer_, totalProtectedSpecies))

fedAllSp <- fedSp %>% 
  group_by(state) %>% count() %>% rename(totalListedSpecies = n)

iucnAllSp <- iucnSp %>% 
  group_by(state) %>% count() %>% rename(totalListedSpecies = n)



# ANALYSES ----------------------------------------------------------------


# * IUCN Coverage ---------------------------------------------------------


# * * state coverage of IUCN ----------------------------------------------------------


## of imperiled species that are not protected by any state
iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% 
  write.csv('data/iucnStateMissing.csv')

#total not listed = 849 / 1624 imperiled species
# 52% imperiled species not protected by any state in the US
## of these 849, 633 are animals and 194 are plants

# separate analyses by phylum
iucnPhyla <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(phylum = tolower(phylum)) %>% 
  group_by(phylum) %>% count() %>% rename(totalIUCNsp = n)

iucnStatePhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp) %>% 
  ungroup() %>% 
  # combine phyla of kindom fungi and plants
  mutate(phylum = if_else(phylum %in% c('ascomycota', 'basidiomycota'), 'fungi', phylum)) %>% 
  mutate(phylum = if_else(phylum %in% c('bryophyta', 'marchantiophyta', 'tracheophyta'), 'plants', phylum)) %>% 
  group_by(phylum) %>% 
  mutate(notListed = sum(notListed), totalIUCNsp = sum(totalIUCNsp), proportionNotListed = notListed/totalIUCNsp) %>% 
  distinct(phylum, .keep_all = TRUE) %>% 
  filter(totalIUCNsp >= 10) %>% 
  mutate(listed = totalIUCNsp-notListed)


# separate by class
iucnClass <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(totalIUCNsp = n)

iucnStateClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)


## ** gtest for state list

### phyla g-test

iucnMat <- iucnStatePhylaSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'phylum') %>% 
  data.matrix()

GTest(iucnMat)
#p-value < 2.2e-16

# write function to run pairwise tests
FUN = function(i,j){
  GTest(matrix(c(iucnMat[i,1], iucnMat[i,2],
                 iucnMat[j,1], iucnMat[j,2]),
               nrow = 2,
               byrow = TRUE),
        correct = 'none')$ p.value
}

phyla.state.gtable <- pairwise.table(FUN,
                                    rownames(iucnMat),
                                    p.adjust.method = 'none')
#           arthropoda        fungi       plants     chordata  cnidaria
# fungi    2.482906e-01           NA           NA           NA        NA
# plants   0.000000e+00 2.676305e-08           NA           NA        NA
# chordata 0.000000e+00 1.308959e-06 1.279604e-02           NA        NA
# cnidaria 5.796136e-03 3.696696e-03 3.136989e-01 6.710513e-01        NA
# mollusca 3.526877e-08 1.435768e-03 2.623457e-13 4.046280e-07 0.3496248
 


### vertebrate class g-test

iucnStateClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% mutate(phylum = tolower(phylum),
                                                                class = tolower(class)) %>% 
  filter(phylum == 'chordata') %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp, listed = totalIUCNsp - notListed) %>% 
  filter(totalIUCNsp >= 10)

iucnMatChor <- iucnStateClassSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'class') %>% 
  data.matrix()

GTest(iucnMatChor)
#p-value = 7.457e-12

# pairwise tests
FUN2 = function(i,j){
  GTest(matrix(c(iucnMatChor[i,1], iucnMatChor[i,2],
                 iucnMatChor[j,1], iucnMatChor[j,2]),
               nrow = 2,
               byrow = TRUE),
        correct = 'none')$ p.value
}

class.gtest <- pairwise.table(FUN2,
                              rownames(iucnMatChor),
                              p.adjust.method = 'none')
#                 actinopterygii     amphibia         aves chondrichthyes  mammalia
# amphibia         9.528538e-01           NA           NA             NA        NA
# aves             1.589119e-01 2.782451e-01           NA             NA        NA
# chondrichthyes   2.004952e-12 6.906491e-10 7.732035e-09             NA        NA
# mammalia         8.565065e-01 9.183808e-01 2.555050e-01   1.238099e-09        NA
# reptilia         3.268430e-02 7.910533e-02 4.263264e-03   6.586953e-13 0.1095269


#Now for iucn species that are on at least one state list, assess the number of states they
# are found in and how many of those states list them


iucnByState <- iucnSp %>% filter(acceptedName %in% stateImp$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

imperiledPerState <- stateImp %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesListedIn = n) %>% left_join(iucnByState, by = 'acceptedName') %>% 
  mutate(proportionListed = statesListedIn/statesOccurIn)

#1/3 (248/775: 32%) are not protected by all the states they occur in, making total number of sp
#not protected in their entire U.S. range = 68%.
#of these 248, 215 are animals and 32 are plants

imperiledPerState %>% filter(proportionListed >= 1) %>% # = 527
  filter(statesOccurIn == 1) # = 483
# Of the sp. that are protected in their entire range, 483/527: 92% are state endemics 
## (i.e. are only found in 1 state) of the endemics in one state, 303 are plants (63%)

## redo phyla and class analyses with just the species that are state listed in their entire range

stateEntireRange <- imperiledPerState %>% filter(proportionListed >= 1)

iucnEntireStatePhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateEntireRange$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)

iucnEntireStateClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateEntireRange$acceptedName)) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)


# ** Partition by state 
## find the proportion of iucn imperiled species that are listed in each state

iucnStateListed <- map_dfr(unique(iucnSp$state),
                function(x) {
                  st <- stateImp %>%
                    filter(state == x)
                  iucnSp %>%
                    filter(state == x) %>%
                    filter(acceptedName %in% st$acceptedName)
                }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(iucnStateListed = n) 

iucnStateListedSummary <-  iucnSp %>% group_by(state) %>% count() %>% rename(totalImperiled = n) %>% 
  left_join(iucnStateListed, by = 'state') %>% mutate(proportionListed = iucnStateListed/totalImperiled,
                                           proportionNotListed = (totalImperiled - iucnStateListed)/totalImperiled) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))

#RESULTS
## Hawaii protects highest proportion of species, which is only 60%
## only 3 states protect more than 50% of iucn imperiled species
## So 94% of states protect 50% or less
## 37/50 states protect less than 1/3 of imperiled species
## The average proportion of iucn species not listed by states is 76%


#now reverse, calculate how many state listed species are not on IUCN list

stateNotIUCN <- stateImp %>% group_by(state) %>% filter(!(acceptedName %in% iucnSp$acceptedName))


# check if these species have been assessed by IUCN

usIUCN <- read_csv('data/US_IUCN_acceptedName.csv')

notIUCNlisted <- stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% group_by(acceptedName) %>% 
  count() %>% rename(notIUCNlisted = n)

stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  group_by(taxon_general) %>%  count() %>% rename(notIUCNlisted = n) %>% ungroup()
#4839 species total, 3445 are plants:
# 1 Amphibian                 2
# 2 Amphibians               74
# 3 Birds                   183
# 4 Fish                    256
# 5 Fungi                    20
# 6 Invertebrates           631
# 7 Mammals                 113
# 8 Plants                 3445
# 9 Reptiles                115

notIUCNlisted %>% filter(!(acceptedName %in% usIUCN$acceptedName))
##only 32% of the missing species have been assessed by IUCN,
##3314/4839 not assessed

notIUCNphyla <- stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% 
  filter(!(acceptedName %in% usIUCN$acceptedName)) %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  group_by(taxon_general) %>% count()
# 2823 of these species (85%) that are listed by states but not IUCN assessed are plants



# * * federal (ESA) coverage of IUCN --------------------------------------------------

#overall # of imperiled species that are not protected by the ESA
iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName))

#total not listed =  989 / 1624 total imperiled species
# 61% imperiled species not protected by the ESA

# separate analyses by phylum

iucnFedPhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp) %>% 
  ungroup() %>% 
  # combine phyla of kindom fungi and plants
  mutate(phylum = if_else(phylum %in% c('ascomycota', 'basidiomycota'), 'fungi', phylum)) %>% 
  mutate(phylum = if_else(phylum %in% c('bryophyta', 'marchantiophyta', 'tracheophyta'), 'plants', phylum)) %>% 
  group_by(phylum) %>% 
  mutate(notListed = sum(notListed), totalIUCNsp = sum(totalIUCNsp), proportionNotListed = notListed/totalIUCNsp) %>% 
  distinct(phylum, .keep_all = TRUE) %>% 
  filter(totalIUCNsp >= 10) %>% 
  mutate(listed = totalIUCNsp-notListed)

# separate by class

iucnFedClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName)) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)

# ** g-tests

##Phyla g-test

iucnMat <- iucnFedPhylaSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'phylum') %>% 
  data.matrix()

GTest(iucnMat)
#p-value < 2.2e-16


phyla.fed.gtable <- pairwise.table(FUN,
                                     rownames(iucnMat),
                                     p.adjust.method = 'none')
#           arthropoda        fungi       plants   chordata  cnidaria
# fungi    3.867101e-01           NA           NA         NA        NA
# plants   0.000000e+00 3.758564e-07           NA         NA        NA
# chordata 1.741829e-12 9.014512e-04 2.476950e-10         NA        NA
# cnidaria 9.952794e-01 5.702326e-01 7.417483e-04 0.04794384        NA
# mollusca 1.177300e-06 7.618115e-03 4.829470e-14 0.04635084 0.1417872


### vertebrate class g-test

iucnFedClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName)) %>% mutate(phylum = tolower(phylum),
                                                                class = tolower(class)) %>% 
  filter(phylum == 'chordata') %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp, listed = totalIUCNsp - notListed) %>% 
  filter(totalIUCNsp >= 10)

iucnMatChor <- iucnFedClassSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'class') %>% 
  data.matrix()

GTest(iucnMatChor)
#p-value = 8.377e-06

class.gtest <- pairwise.table(FUN2,
                              rownames(iucnMatChor),
                              p.adjust.method = 'none')
#               actinopterygii     amphibia         aves chondrichthyes  mammalia
# amphibia         6.163506e-01           NA           NA             NA        NA
# aves             6.004746e-01 9.518674e-01           NA             NA        NA
# chondrichthyes   9.768812e-08 1.123128e-05 2.401790e-06             NA        NA
# mammalia         2.553170e-01 5.823544e-01 5.066868e-01   1.233327e-04        NA
# reptilia         4.137846e-01 2.970152e-01 2.762208e-01   5.720850e-07 0.1316592


#assess the number of states IUCN ESA listed species are found in compared to 
#the number of states that ESA lists them in

iucnByStateFed <- iucnSp %>% filter(acceptedName %in% fedSp$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

imperiledPerStateFed <- fedSp %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesESAListedIn = n) %>% left_join(iucnByStateFed, by = 'acceptedName') %>% 
  mutate(proportionListed = statesESAListedIn/statesOccurIn)

#79 / 635 (12%) of ESA listed IUCN imperiled species are not listed in all the states they are
#believed to occur in. May be due to differences in distribution info or delisting of specific subpop/subspec, 
# making total number of sp not protected in their entire U.S. range = 66%.

imperiledPerStateFed %>% filter(proportionListed >= 1) %>% # = 556
  filter(statesOccurIn == 1) # = 464
# Of the sp. that are protected in their entire range, 464/556: 83% are state endemics (i.e. are only found in 1 state)



## parition by state
iucnFedListed <- map_dfr(unique(iucnSp$state),
                           function(x) {
                             st <- fedSp %>%
                               filter(state == x)
                             iucnSp %>%
                               filter(state == x) %>%
                               filter(acceptedName %in% st$acceptedName)
                           }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(iucnFedListed = n) 

iucnFedListedSummary <-  iucnSp %>% group_by(state) %>% count() %>% rename(totalImperiled = n) %>% 
  left_join(iucnFedListed, by = 'state') %>% mutate(proportionListed = iucnFedListed/totalImperiled,
                                                      proportionNotListed = (totalImperiled - iucnFedListed)/totalImperiled) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))
# ave proporiton of IUCN sp NOT listed per state by ESA is 82%



# * * State and Fed Combined Coverage ------------------------------------------------------------------------

## get total number of species not protected
iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName) & !(acceptedName %in% stateImp$acceptedName))
# 795 / 1624 unprotected (49%)

## now combine state and ESA lists to see total coverage by both fed and state laws

iucnStateFedPhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName) & !(acceptedName %in% stateImp$acceptedName)) %>% #795 / 1624
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp) %>% 
  ungroup() %>% 
  # combine phyla of kindom fungi and plants
  mutate(phylum = if_else(phylum %in% c('ascomycota', 'basidiomycota'), 'fungi', phylum)) %>% 
  mutate(phylum = if_else(phylum %in% c('bryophyta', 'marchantiophyta', 'tracheophyta'), 'plants', phylum)) %>% 
  group_by(phylum) %>% 
  mutate(notListed = sum(notListed), totalIUCNsp = sum(totalIUCNsp), proportionNotListed = notListed/totalIUCNsp) %>% 
  distinct(phylum, .keep_all = TRUE) %>% 
  filter(totalIUCNsp >= 10) %>% 
  mutate(listed = totalIUCNsp-notListed)

iucnStateFedClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName) & !(acceptedName %in% stateImp$acceptedName)) %>% mutate(phylum = tolower(phylum),
                                                                                                          class = tolower(class)) %>% 
  filter(phylum == 'chordata') %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp, listed = totalIUCNsp - notListed) %>% 
  filter(totalIUCNsp >= 10)

### ** g-test

## phlya g-test

iucnMat <- iucnStateFedPhylaSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'phylum') %>% 
  data.matrix()

GTest(iucnMat)
#p-value < 2.2e-16


phyla.fed.state.gtable <- pairwise.table(FUN,
                                   rownames(iucnMat),
                                   p.adjust.method = 'none')
#          arthropoda        fungi       plants     chordata  cnidaria
# fungi    1.106385e-01           NA           NA           NA        NA
# plants   0.000000e+00 1.753708e-08           NA           NA        NA
# chordata 0.000000e+00 1.095135e-07 2.785398e-01           NA        NA
# cnidaria 1.899932e-02 3.696696e-03 2.861741e-01 4.175090e-01        NA
# mollusca 4.117733e-08 3.075887e-04 6.774744e-10 2.701492e-07 0.5686835


### vertebrate class g-test


iucnMatChor <- iucnStateFedClassSummary %>% dplyr::select(notListed, listed) %>% 
  column_to_rownames(var = 'class') %>% 
  data.matrix()

GTest(iucnMatChor)
#p-value = 1.488e-14

class.fed.state.gtest <- pairwise.table(FUN2,
                              rownames(iucnMatChor),
                              p.adjust.method = 'none')

#                 actinopterygii     amphibia         aves chondrichthyes  mammalia
# amphibia         1.863791e-01           NA           NA             NA        NA
# aves             4.537139e-02 8.718873e-03           NA             NA        NA
# chondrichthyes   2.475797e-14 1.189049e-13 2.420620e-09             NA        NA
# mammalia         9.927275e-01 3.149303e-01 1.649913e-01   1.338559e-10        NA
# reptilia         6.870915e-02 5.635347e-01 3.262532e-03   1.288969e-13 0.1406868


#combined, how many species are not protected by all the states they occur in
iucnByStateComb <- iucnSp %>% filter(acceptedName %in% fedSp$acceptedName | acceptedName %in% stateImp$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

FedStateComb <- fedSp %>% mutate(Notes = NA) %>% rename(Ranking = status) %>% 
  rbind(stateImp) %>% distinct(acceptedName, state, .keep_all = TRUE)

imperiledPerStateComb <- FedStateComb %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesCombListedIn = n) %>% left_join(iucnByStateComb, by = 'acceptedName') %>% 
  mutate(proportionListed = statesCombListedIn/statesOccurIn)

#183 / 829 (22%) of the species are not ESA listed in all the states they are
#believed to occur in, making total number of sp
#not protected in their entire U.S. range = 60%.

imperiledPerStateComb %>% filter(proportionListed >= 1) %>% # = 646
  filter(statesOccurIn == 1) # = 530
# Of the sp. that are protected in their entire range, 530/646: 82% are state endemics (i.e. are only found in 1 state)



## parition by state
iucnCombListed <- map_dfr(unique(iucnSp$state),
                         function(x) {
                           st <- FedStateComb %>%
                             filter(state == x)
                           iucnSp %>%
                             filter(state == x) %>%
                             filter(acceptedName %in% st$acceptedName)
                         }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(iucnCombListed = n) 

iucnCombListedSummary <-  iucnSp %>% group_by(state) %>% count() %>% rename(totalImperiled = n) %>% 
  left_join(iucnCombListed, by = 'state') %>% mutate(proportionListed = iucnCombListed/totalImperiled,
                                                    proportionNotListed = (totalImperiled - iucnCombListed)/totalImperiled) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))
# ave proporiton of IUCN sp unlisted per state by esa and state combined is 71%



# * State coverage of ESA ------------------------------------------------

fedSp %>% filter(!(acceptedName %in% stateImp$acceptedName)) %>% distinct(acceptedName, .keep_all = TRUE) 
  # %>%  write.csv('data/fedMissingSpecies.csv')
# 219 / 1386 species not on any state list (16%)

## partition by state
FedStateListed <- map_dfr(unique(fedSp$state),
                           function(x) {
                             st <- stateImp %>%
                               filter(state == x)
                             fedSp %>%
                               filter(state == x & acceptedName %in% st$acceptedName)
                           }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(fedStateListed = n) 

FedStateListedSummary <-  fedSp %>% group_by(state) %>% count() %>% rename(totalFedListed = n) %>% 
  left_join(FedStateListed, by = 'state') %>% mutate(proportionListed = fedStateListed/totalFedListed,
                                           proportionNotListed = (totalFedListed - fedStateListed)/totalFedListed) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))

## currently 1/3 of states (16/50) protect 50% or less of fed species in their state

## on average, 37% of fed species are not listed by states (over 1/3) 



# * Consequences of ESA removal ---------------------------

## read in legislation file

stateLeg <- read_csv('data/stateLegislationSummary.csv') %>% dplyr::select(state:last_updated) %>% 
  mutate(state = tolower(state))

## remove esa species that are automatically included by states

st <- stateLeg %>% filter(ESA_included == 'Y') %>% pull(state)
  # 32 states that include some 'ESA - listed species' on state lists automatically

#filter states by ESA categories they list (index # if states in alphabetical order)

allESA <- st[c(1:9, 11:13, 15, 16, 18, 20:27,29, 31, 32)] 

fe <- st[c(10, 14, 19, 28)]

animals <- st[17]

feAnimals <- st[30]

# now filter out those species from those state lists

removeAllESA <- map_dfr(allESA,
                        function(x) {
                          fsp <- fedSp %>%
                            filter(state == x)
                          stateImp %>%
                            filter(state == x) %>% 
                            filter(!(acceptedName %in% fsp$acceptedName))
                        })

removeFEOnly <- map_dfr(fe,
                        function(x) {
                          fsp <- fedSp %>%
                            filter(state == x & status == 'E')
                          stateImp %>%
                            filter(state == x) %>% 
                            filter(!(acceptedName %in% fsp$acceptedName))
                        })

removeESAAnimals <- map_dfr(animals,
                            function(x) {
                              fsp <- fedSp %>%
                                filter(state == x & taxon_general != 'Plants')
                              stateImp %>%
                                filter(state == x) %>% 
                                filter(!(acceptedName %in% fsp$acceptedName))
                            })

removeFEAnimals <-  map_dfr(feAnimals,
                            function(x) {
                              fsp <- fedSp %>%
                                filter(state == x & taxon_general != 'Plants' & status == 'E')
                              stateImp %>%
                                filter(state == x) %>% 
                                filter(!(acceptedName %in% fsp$acceptedName))
                            })

indepLists <- stateImp %>% filter(!(state %in% st))

# now rbind all lists

stateImp_noESA <- rbind(removeAllESA, removeFEOnly, removeESAAnimals, removeFEAnimals, indepLists)


# * * IUCN coverage NEW --------------------------------------------------------------------------

#overall # of imperiled species that are not protected by any state
iucnSpNoESA <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp_noESA$acceptedName))

#total not listed = 850 -> 1309 / 1624 total imperiled species = 30% increase in unprotected sp
# 52% -> 81%  imperiled species not protected by any state in the US

#how many were protected with state and esa
iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% FedStateComb$acceptedName))
# 795 were not protected before, so 1309 - 795 = 514 specieslose protection

#separate new unprotected by phylum

iucnStatePhylaSummary_noESA <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp_noESA$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed_noESA = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  ungroup() %>% 
  # combine phyla of kindom fungi and plants
  mutate(phylum = if_else(phylum %in% c('ascomycota', 'basidiomycota'), 'fungi', phylum)) %>% 
  mutate(phylum = if_else(phylum %in% c('bryophyta', 'marchantiophyta', 'tracheophyta'), 'plants', phylum)) %>% 
  group_by(phylum) %>% 
  mutate(notListed_noESA = sum(notListed_noESA), totalIUCNsp = sum(totalIUCNsp)) %>% 
  ungroup() %>% 
  distinct(phylum, .keep_all = TRUE) %>% 
  filter(totalIUCNsp >= 10) %>% 
  mutate(proportionNotListed = notListed_noESA/totalIUCNsp, listed = totalIUCNsp - notListed_noESA)
  

# separate by class
iucnStateClassSummary_noESA <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp_noESA$acceptedName)) %>% mutate(phylum = tolower(phylum),
                                                                      class = tolower(class)) %>%
  filter(phylum == 'chordata') %>% 
  group_by(class) %>% count() %>% rename(notListed_noESA = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed_noESA = notListed_noESA/totalIUCNsp) %>% 
  ungroup() %>% 
  filter(totalIUCNsp >= 10) %>% 
  mutate(listed = totalIUCNsp - notListed_noESA)



## g-tests

### phyla g-test

iucnMat <- iucnStatePhylaSummary_noESA %>% dplyr::select(phylum, notListed_noESA, listed) %>% 
  column_to_rownames(var = 'phylum') %>% 
  data.matrix()

GTest(iucnMat)
#p-value < 2.2e-16

phyla.state.gtable <- pairwise.table(FUN,
                                     rownames(iucnMat),
                                     p.adjust.method = 'none')
#           arthropoda        fungi       plants     chordata  cnidaria
# fungi    3.203429e-02           NA           NA           NA        NA
# plants   4.467456e-01 4.909020e-02           NA           NA        NA
# chordata 1.110223e-16 4.293312e-06 0.000000e+00           NA        NA
# cnidaria 1.095643e-02 7.571306e-04 4.148389e-03 9.534122e-01        NA
# mollusca 1.208740e-04 9.030777e-04 4.337624e-08 4.076825e-06 0.2277775


### class g-test

iucnMatChor <- iucnStateClassSummary_noESA %>% dplyr::select(class,notListed_noESA, listed) %>% 
  column_to_rownames(var = 'class') %>% 
  data.matrix()

GTest(iucnMatChor)
# p-value = 1.714e-11

#pairwise tests

class.gtest <- pairwise.table(FUN2,
                              rownames(iucnMatChor),
                              p.adjust.method = 'none')
#               actinopterygii     amphibia         aves chondrichthyes  mammalia
# amphibia         4.236644e-02           NA           NA             NA        NA
# aves             2.746841e-02 6.060031e-04           NA             NA        NA
# chondrichthyes   5.266098e-06 1.277840e-07 2.232630e-03             NA        NA
# mammalia         4.852769e-03 4.501210e-01 4.732285e-05   9.858050e-09        NA
# reptilia         4.360624e-04 1.423739e-01 3.320414e-06   7.209324e-10 0.4710849



#Now compare, for iucn state listed species, the number of states they live in compared to 
#the number of states that list them
iucnByStateNoESA <- iucnSp %>% filter(acceptedName %in% stateImp_noESA$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

imperiledPerStateNoESA <- stateImp_noESA %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesListedIn = n) %>% left_join(iucnByStateNoESA, by = 'acceptedName') %>% 
  mutate(proportionListed = statesListedIn/statesOccurIn)
##210 sp not protected in all the states they occur in


## find the proportion of iucn imperiled species that are listed in each state

iucnStateListedNoESA <- map_dfr(unique(iucnSp$state),
                           function(x) {
                             sta <- stateImp_noESA %>%
                               filter(state == x)
                             iucnSp %>%
                               filter(state == x) %>%
                               filter(acceptedName %in% sta$acceptedName)
                           }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(iucnStateListed_noESA = n) 

iucnStateListedSummaryNoESA <-
  iucnSp %>% group_by(state) %>% count() %>% rename(totalImperiled = n) %>%
  left_join(iucnStateListedNoESA, by = 'state') %>%
  mutate(
    iucnStateListed_noESA = if_else(
      is.na(iucnStateListed_noESA),
      as.integer(0),
      iucnStateListed_noESA
    )
    ,
    proportionListed = iucnStateListed_noESA / totalImperiled,
    proportionNotListed = (totalImperiled - iucnStateListed_noESA) /
      totalImperiled
  )
#RESULTS
## now highest is VA with 44% protected
## now no states protect more than half of the imperiled species in their borders
## The average proportion of iucn species not listed by states is 75% -> now 85%

## * * Fed coverage new --------------------------------------------------------------------


fedSp %>% filter(!(acceptedName %in% stateImp_noESA$acceptedName)) %>%
  distinct(acceptedName, .keep_all = TRUE)

#1017 / 1386 = 73% 
FedStateListedNoESA <- map_dfr(unique(fedSp$state),
                          function(x) {
                            st <- stateImp_noESA %>%
                              filter(state == x)
                            fedSp %>%
                              filter(state == x & acceptedName %in% st$acceptedName)
                          }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(fedStateListed_noESA = n) 

FedStateListedSummaryNoESA <-
  fedSp %>% group_by(state) %>% count() %>% rename(totalFedListed = n) %>%
  left_join(FedStateListedNoESA, by = 'state') %>%
  mutate(
    fedStateListed_noESA = if_else(is.na(fedStateListed_noESA), as.integer(0), fedStateListed_noESA),
    proportionListed = fedStateListed_noESA / totalFedListed,
    proportionNotListed = (totalFedListed - fedStateListed_noESA) / totalFedListed
  )
# now average proportion of ESA species not listed in each state is 84% (b/c lots of states are 0 now)


## new senate bill, remove all ESA species that only occur in 1 state, what are those consequences?

### ESA species that are only found in 1 state

fedSingleState <- fedSp  %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n) %>% filter(statesOccurIn == 1)
# 1064

# left bind with original data
fedSp %>% filter(acceptedName %in% fedSingleState$acceptedName) %>% left_join(fedSingleState, by = 'acceptedName')

# how many of these would remain protected by states

fedSingleState %>% filter(acceptedName %in% stateImp$acceptedName)
# 862 still protected by state, 19% would lose total protection

## but, what if those species were only on state lists because they were ESA listed...
fedSingleState %>% filter(acceptedName %in% stateImp_noESA$acceptedName)
## now only 202 independently protected by states 

#analyze the ones that would lose total protection under current state laws
fedSp %>% filter(acceptedName %in% fedSingleState$acceptedName) %>% left_join(fedSingleState, by = 'acceptedName') %>%
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% group_by(taxon_general) %>% count()
# taxon_general     n
#  1 Amphibians       10
# 2 Birds             5
# 3 Fish             17
# 4 Invertebrates    61
# 5 Mammals           8
# 6 Plants          100
# 7 Reptiles          1



## * Model Selection --------------------------------------------------------------------------------------------

#read in state info file

stateSummary <- read_csv('data/state_info_modSelection.csv')

### create refined set of vars by removing correlated ones

stateSumRefine <- stateSummary %>% dplyr::select(-state)

## transform non-normal variables
stateSumRefine <- stateSumRefine %>% mutate(logTotalProtected = log10(totalProtectedSpecies + 1),
                                            sqrtTotalProtected = sqrt(totalProtectedSpecies), #sqrt has better shapiro results
                                            log_atRisk = log10(percent_atrisk), 
                                            logEndemic = log10(endemic_sp + 1), logPlants = log10(percent_plants),
                                            percentBirds = log10(percent_birds), 
                                            logArea = log10(area_sqkm), logTotal_sp = log10(total_sp),
                                            logPop_density = log10(popDensity), logHunt_fish = log10(percent_HuntFish),
                                            sqrtFarmland = sqrt(percent_Farmland), logFedland = log10(percent_FedLand),
                                            logIncome = log10(Income_2019), logTotalImperiled = log10(totalImperiled))

## check correlations
mat <- stateSumRefine %>% dplyr::select(-c(totalProtectedSpecies, logTotalProtected, sqrtTotalProtected))

corr <- rcorr(as.matrix(mat), type = 'pearson')


### create full model (without correlated vars)
full.mod <- lm(logTotalProtected ~ party_trend_rank + logTotal_sp +
                 percent_mammals + percent_reptiles + percent_amphib + logPlants + 
                  logPop_density + logArea + logHunt_fish + sqrtFarmland + logFedland +
                 logIncome + logTotalImperiled, 
               data = stateSumRefine)

  
#now run model selection based on AIC using stepAIC function from MASS package

## start with backwards

stepAIC(full.mod, direction = 'backward', criteria = "AICc")
#best = AIC -47.42 ; logTotalProtected ~ party_trend_rank + logArea + logTotal_sp + 
#percent_reptiles + logPop_density

##now do forwards
modZero <- lm(logTotalProtected ~ 1, data = stateSumRefine)
stepAIC(modZero, direction = 'forward', scope = list(lower = modZero, upper = full.mod))
# best = AIC -46.21; logTotalProtected ~ logPop_density + party_trend_rank + logArea

##now do both
stepAIC(full.mod, direction = 'both', k = length(nrow(stateSumRefine)), 
                  scope = list(lower = modZero, upper = full.mod))
# best is same as backwards
## but total sp and %reptiles are pretending variables, corrected AIC gives model
## same as forward stepwise model

modAIC <- MASS::stepAIC(modZero, direction = 'forward', scope = list(lower = modZero, upper = full.mod))
summary(modAIC)
#best model

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -3.141356   1.548633  -2.028   0.0483 *  
#   logPop_density    0.949014   0.212463   4.467 5.13e-05 ***
#   party_trend_rank -0.017363   0.006859  -2.532   0.0148 *  
#   logArea           0.617591   0.247580   2.495   0.0163 * 
#                   2.5 %      97.5 %
#   (Intercept)     -6.2585929 -0.02412000
# logPop_density    0.5213478  1.37668061
# party_trend_rank -0.0311691 -0.00355717
# logArea           0.1192378  1.11594462

 

# # test if corr still significant without Hawaii (outlier)
stateSumNoHI <- stateSummary %>% filter(state != 'hawaii') %>% 
  dplyr::select(-state) %>% 
  mutate(logTotalProtected = log10(totalProtectedSpecies + 1),
                                            sqrtTotalProtected = sqrt(totalProtectedSpecies), #sqrt has better shapiro results
                                            log_atRisk = log10(percent_atrisk), 
                                            logEndemic = log10(endemic_sp + 1), logPlants = log10(percent_plants),
                                            percentBirds = log10(percent_birds), 
                                            logArea = log10(area_sqkm), logTotal_sp = log10(total_sp),
                                            logPop_density = log10(popDensity), logHunt_fish = log10(percent_HuntFish),
                                            sqrtFarmland = sqrt(percent_Farmland), logFedland = log10(percent_FedLand),
                                            logIncome = log10(Income_2019), logTotalImperiled = log10(totalImperiled))

full.mod <- lm(logTotalProtected ~ party_trend_rank + logTotal_sp + coast + percentBirds +
                 percent_mammals + percent_reptiles + percent_amphib + logPlants + 
                 logPop_density + logArea + logHunt_fish + sqrtFarmland + logFedland +
                 logIncome + logTotalImperiled, 
               data = stateSumNoHI)

modAIC <- MASS::stepAIC(full.mod, direction = "both")
#results the same

## now separate by wildlife and plants and re-run model selection

### wildlife
stateWildlife <- stateImp %>% filter(taxon_general != 'Plants') %>%
  add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% dplyr::select(state, totalWildlife = n) %>% ungroup() %>% 
  mutate(totalWildlife = if_else(totalWildlife == 1, as.integer(0) ,totalWildlife)) %>% left_join(stateSummary, by = 'state') %>% 
  dplyr::select(-state) %>% 
  mutate(sqrtTotalWildlife = sqrt(totalWildlife),
         logTotalWildlife = log10(totalWildlife + 1), #sqrt has better shapiro results
         log_atRisk = log10(percent_atrisk), 
         logEndemic = log10(endemic_sp + 1), logPlants = log10(percent_plants),
         percentBirds = log10(percent_birds), 
         logArea = log10(area_sqkm), logTotal_sp = log10(total_sp),
         logPop_density = log10(popDensity), logHunt_fish = log10(percent_HuntFish),
         sqrtFarmland = sqrt(percent_Farmland), logFedland = log10(percent_FedLand),
         logIncome = log10(Income_2019), logTotalImperiled = log10(totalImperiled))

wildMod <- lm(totalWildlife ~ party_trend_rank + logTotal_sp + coast + percentBirds +
                percent_mammals + percent_reptiles + percent_amphib + logPlants + 
                logPop_density + logArea + logHunt_fish + sqrtFarmland + logFedland +
                logIncome + logTotalImperiled, data = stateWildlife)

wildModZero <- lm(logTotalWildlife ~ 1, data = stateWildlife)

stepModWild <- MASS::stepAIC(wildMod, direction = 'both', scope = list(lower = wildModZero, 
                                                                       upper = wildMod))
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)       1.551662   2.237409   0.694  0.49164   
# party_trend_rank -0.017484   0.006644  -2.631  0.01168 * 
#   logTotal_sp      -1.728155   0.909232  -1.901  0.06391 . 
# percent_reptiles  0.477219   0.174556   2.734  0.00898 **
#   logPop_density    0.744644   0.235692   3.159  0.00286 **
#   logArea           0.781392   0.289325   2.701  0.00979 **


# plants

stateWithPlants <- stateImp %>% filter(taxon_general == 'Plants') %>% pull(state) %>% unique()
stateNoPlant <- stateImp %>% filter(!(state %in% stateWithPlants)) %>% pull(state) %>% unique()

statePlants <- stateImp %>% 
  filter(taxon_general == 'Plants') %>%
  add_row(state = c(stateNoPlant,'west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>%
  count() %>% ungroup() %>% 
  dplyr::select(state, totalPlants = n) %>% 
  mutate(totalPlants = if_else(totalPlants == 1, as.integer(0),totalPlants)) %>% left_join(stateSummary, by = 'state') %>% 
  dplyr::select(-state) %>% 
  mutate(logTotalPlants = log10(totalPlants + 1),
         sqrtPlants = sqrt(totalPlants),
         log_atRisk = log10(percent_atrisk), 
         logEndemic = log10(endemic_sp + 1), logPlants = log10(percent_plants),
         percentBirds = log10(percent_birds), 
         logArea = log10(area_sqkm), logTotal_sp = log10(total_sp),
         logPop_density = log10(popDensity), logHunt_fish = log10(percent_HuntFish),
         sqrtFarmland = sqrt(percent_Farmland), logFedland = log10(percent_FedLand),
         logIncome = log10(Income_2019), logTotalImperiled = log10(totalImperiled))

plantMod <- lm(logTotalPlants ~ party_trend_rank + logTotal_sp + percentBirds + coast +
                 percent_mammals + percent_reptiles + percent_amphib + logPlants + 
                 logPop_density + logArea + logHunt_fish + sqrtFarmland + logFedland +
                 logIncome + logTotalImperiled, data = statePlants)
plantModZero <- lm(logTotalPlants ~ 1, data = statePlants)


stepModPlant <- MASS::stepAIC(plantMod, direction = 'both', scope = list(lower = plantModZero,
                                                                         upper = plantMod))
summary(stepModPlant)
# top model same but mammals (-0.57) instead of pop density
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)      -17.06479   12.38432  -1.378  0.17519   
# party_trend_rank  -0.03014    0.01076  -2.802  0.00753 **
#   percent_mammals   -0.57285    0.20280  -2.825  0.00709 **
#   logArea            0.70153    0.32258   2.175  0.03507 * 
#   logHunt_fish      -1.23112    0.89449  -1.376  0.17568   
# logIncome          3.76155    2.50573   1.501  0.14045 

#NOTE: Model averaging was performed using JMP software for final model results reported in manuscript

# FIGURES -------------------------------------------------------------

## * correlation figures -------------

### protected sp vs. imperiled sp

corplot1 <- ggscatter(stateSummary, x = "totalImperiled", y = "totalProtectedSpecies",
                      size = 7, alpha = 0.4, color = 'black',
                      add = "reg.line", conf.int = TRUE,
                      xlim = c(20, 600),
                      xlab = "IUCN Imperiled Species", ylab = "State Listed Species",
                      font.x = c(18,'bold', 'black'), ylim = c(0, 801), font.y = c(18, 'bold', 'black'),
                      add.params = list(color = "black",
                                        fill = "lightgray"))+
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600), labels = c(0, 100, 200, 300, 400, 500, 600))+
  theme(plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 450, label.y = 810, size = 5
  )

# see difference without Hawaii (outlier)

stateSummary %>% filter(state != 'hawaii') %>% 
  ggscatter(x = "totalImperiled", y = "totalProtectedSpecies",
            size = 7, alpha = 0.4, color = 'black',
            add = "reg.line", conf.int = TRUE,
            # cor.coef = TRUE, cor.method = "spearman",
            # cor.coeff.args = list(label.x = 35, label.y = 840, size = 7),
            xlab = "IUCN Imperiled Species", ylab = "State Listed Species",
            font.x = c(21,'bold', 'black'), font.y = c(21, 'bold', 'black'),
            add.params = list(color = "black",
                              fill = "lightgray"))+
  theme(plot.background = element_rect(fill = "#F5F3EA", color = NA), 
        panel.background = element_rect(fill = "#F5F3EA", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 10, label.y = 810, size = 5
  )

# total protected vs. political party ranking
corplot2 <- ggscatter(stateSummary, x = "party_trend_rank", y = "totalProtectedSpecies",
                      color = 'party_trend_rank', size = 7,
                      add = "reg.line", conf.int = TRUE,
                      
                      xlab = "", ylab = "State Listed Species",
                      add.params = list(color = "black",
                                        fill = "lightgray"),
                      ylim = c(0, 801),
                      font.y = c(18, 'bold', 'black'),
                      xlim = c(-25, 25))+
  gradient_color(c('blue', 'purple', 'red'))+
  scale_x_continuous(breaks = c(-25, 0, 25), labels = c("-25" = "Most Democratic State", 
                                                        "0" = "Neutral", "25"="Most Republican State"))+
  theme(legend.position = 'none', axis.text.x = element_text(face = 'bold', 
                                                             size = '14', hjust = c(0.2, 0.5, 0.9)),
        axis.ticks.x = element_line(size = 2),
        axis.ticks.length.x = unit(0.2, "cm"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 10, label.y = 810, size = 5
  )

## political party ranking with logged response (total protected)
ggscatter(stateSumRefine, x = "party_trend_rank", y = "logTotalProtected",
          color = 'party_trend_rank', size = 7,
          add = "reg.line", conf.int = TRUE,
          
          xlab = "", ylab = "log(Total State Listed Species)",
          add.params = list(color = "black",
                            fill = "lightgray"),
          ylim = c(0, 3),
          font.y = c(18, 'bold', 'black'),
          xlim = c(-25, 25))+
  gradient_color(c('blue', 'purple', 'red'))+
  scale_x_continuous(breaks = c(-25, 0, 25), labels = c("-25" = "Most Democratic State", 
                                                        "0" = "Neutral", "25"="Most Republican State"))+
  theme(legend.position = 'none', axis.text.x = element_text(face = 'bold', 
                                                             size = '14', hjust = c(0.2, 0.5, 0.9)),
        axis.ticks.x = element_line(size = 2),
        axis.ticks.length.x = unit(0.2, "cm"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 10, label.y = 3, size = 5
  )

# total protected vs. state area
ggscatter(stateSumRefine, x = "logArea", y = "logTotalProtected",
          size = 7, alpha = 0.4, color = 'black',
          add = "reg.line", conf.int = TRUE,
          xlim = c(3.5, 6.5),
          xlab = "Staet Area (Log)", ylab = "State Listed Species",
          font.x = c(18,'bold', 'black'), ylim = c(0, 3), font.y = c(18, 'bold', 'black'),
          add.params = list(color = "black",
                            fill = "lightgray"))+
  scale_x_continuous(breaks = c(3, 4, 5, 6), labels = c(3, 4, 5, 6))+
  theme(plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 6, label.y = 3, size = 5
  )

#total protected vs. human population density
ggscatter(stateSumRefine, x = "logPop_density", y = "totalProtectedSpecies",
          size = 10, alpha = 0.4, color = 'black',
          add = "reg.line", conf.int = TRUE,
          xlim = c(0, 3.05),
          xlab = "log(Population Density)", ylab = "State Listed Species",
          font.x = c(18,'bold', 'black'), ylim = c(0, 805), font.y = c(18, 'bold', 'black'),
          add.params = list(color = "black",
                            fill = "lightgray"))+
  scale_x_continuous(breaks = c(0, 1, 2, 3), labels = c(0, 1, 2, 3))+
  theme(plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 2.2, label.y = 805, size = 5
  )




# * State-listed maps ----------------------------------------


# * * All state-listed species -------------------------------------------------------------


s1 <- ggplot(stateAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalProtectedSpecies), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#E1ECF9' ,
                       mid = '#236AB9',
                       high = '#091D34',
                       midpoint = 400, na.value = 'white', limits = c(0,801),
                      breaks = c(0,200, 400, 600, 800),
                      name = 'State Listed Species',
                      guide = guide_colorbar(title.position = 'top',
                                             direction = 'horizontal', barwidth = unit(10, 'cm'),
                                             barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# * * Split state listed animals and plants -----------------------------------------------


## state-listed animal species
stateWildlife <- stateImp %>% filter(taxon_general != 'Plants') %>%
  add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% dplyr::select(totalWildlife = n) %>% ungroup() %>% 
  mutate(totalWildlife = if_else(totalWildlife == 1, NA_integer_,totalWildlife))


sw <- ggplot(stateWildlife, aes(map_id = state))+
  geom_map(aes(fill = totalWildlife), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#fe9929', high = '#662506',
                       midpoint = 118, na.value = 'white', limits = c(0,239),
                       breaks = c(0,75, 150, 225),
                       name = 'State Listed Animal Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


## state-listed plant species
stateWithPlants <- stateImp %>% filter(taxon_general == 'Plants') %>% pull(state) %>% unique()
stateNoPlant <- stateImp %>% filter(!(state %in% stateWithPlants)) %>% pull(state) %>% unique()

statePlants <- stateImp %>% 
  filter(taxon_general == 'Plants') %>%
  add_row(state = c(stateNoPlant,'west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>%
  count() %>% ungroup() %>% 
  dplyr::select(state, totalPlants = n) %>% 
  mutate(totalPlants = if_else(totalPlants == 1, NA_integer_,totalPlants))


sp <- ggplot(statePlants, aes(map_id = state))+
  geom_map(aes(fill = totalPlants), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#78c679', high = '#004529',
                       midpoint = 356, na.value = 'white', limits = c(0,750),
                       breaks = c(0, 250, 500, 750),
                       name = 'State Listed Plant Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))



### * Fed maps -------------------------------------

# total ESA listed species in each state
f <- ggplot(fedAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalListedSpecies), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#DFE4E7' ,
                       mid = '#6088A1',
                       high = '#033757', 
                       midpoint = 65, na.value = '#011d2e', limits = c(0,130),
                       breaks = c(0,25,50,75,100, 125),
                       name = 'ESA Listed Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))+
  annotate("text", x = c(-99,-125), y = c(24,35),
           label = c(489, 270),
           color = 'black', size = 3, fontface = 'bold')+
  annotate('segment', x = -101, y = 24, xend = -103, yend = 24, 
           arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'), size = 0.2)+
  annotate('segment', x = -123, y = 35, xend = -121, yend = 35,
           arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'), size = 0.2)


#Total ESA listed animals in each state

fedWildlife <- fedSp %>% filter(taxon_general != 'Plants') %>%
  group_by(state) %>% count() %>% dplyr::select(state, totalWildlife = n)

fw <-  ggplot(fedWildlife, aes(map_id = state))+
  geom_map(aes(fill = totalWildlife), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#fe9929', high = '#662506',
                       midpoint = 54, na.value = '#F5F3EA', limits = c(0,107),
                       breaks = c(0, 25, 50, 75, 100),
                       name = 'ESA Listed Animal Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


### Total ESA listed plants in each state

fedPlants <- fedSp %>% filter(taxon_general == 'Plants') %>% 
  group_by(state) %>% count() %>% dplyr::select(state, totalPlants = n)


fp <- ggplot(fedPlants, aes(map_id = state))+
  geom_map(aes(fill = totalPlants), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#78c679', high = '#004529',
                       midpoint = 31, na.value = '#0d3317', limits = c(0,62),
                       breaks = c(0, 15, 30, 45, 60),
                       name = 'ESA Listed Plant Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))+
  annotate("text", x = c(-99,-125), y = c(24,35),
           label = c(410, 171),
           color = 'black', size = 3, fontface = 'bold')+
  annotate('segment', x = -101, y = 24, xend = -103, yend = 24, 
           arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'), size = 0.2)+
  annotate('segment', x = -123, y = 35, xend = -121, yend = 35,
           arrow = arrow(length = unit(0.2, 'cm'), type = 'closed'), size = 0.2)



## proportion of ESA species that are also state-listed

fedStateListedSummary <- FedStateListedSummary %>% mutate(percentProtected = proportionListed*100) %>% 
  mutate(percentProtected = if_else(percentProtected == 0, NA_integer_, as.integer(percentProtected)))

fedProp <- ggplot(fedStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = percentProtected), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#DFE4E7' , mid = '#6088A1', high = '#02273E',
                       midpoint = 50, na.value = 'white', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = '% of ESA Species State Listed',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(7, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.09),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,0,2), "lines"))


## proportion of ESA animal species also state-listed
FedAnimalStateListed <- map_dfr(unique(fedSp$state),
                          function(x) {
                            st <- stateImp %>%
                              filter(state == x)
                            fedSp %>%
                              filter(taxon_general != 'Plants') %>% 
                              filter(state == x & acceptedName %in% st$acceptedName)
                          }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(fedStateListed = n) 

FedAnimalStateListedSummary <-
  fedSp %>% filter(taxon_general != 'Plants') %>% group_by(state) %>% count() %>% rename(totalFedListed = n) %>%
  left_join(FedAnimalStateListed, by = 'state') %>%
  mutate(
    proportionListed = fedStateListed / totalFedListed,
    proportionNotListed = (totalFedListed - fedStateListed) /
      totalFedListed
  ) %>%
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>%
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed)) %>%
  mutate(percentProtected = proportionListed * 100) %>%
  mutate(percentProtected = if_else(
    percentProtected == 0,
    NA_integer_,
    as.integer(percentProtected)
  ))


fedAnimalProp <- ggplot(FedAnimalStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = percentProtected), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#fe9929', high = '#662506',
                       midpoint = 50, na.value = 'white', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = '% of ESA Animal Species State Listed',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(7, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.09),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,0,2), "lines"))

## proportion of ESA plant species also state-listed

FedPlantStateListed <- map_dfr(unique(fedSp$state),
                                function(x) {
                                  st <- stateImp %>%
                                    filter(state == x & taxon_general == 'Plants')
                                  fedSp %>%
                                    filter(taxon_general == 'Plants') %>% 
                                    filter(state == x & acceptedName %in% st$acceptedName)
                                }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(fedStateListed = n) 

FedPlantStateListedSummary <-
  fedSp %>% filter(taxon_general == 'Plants') %>% group_by(state) %>% count() %>% rename(totalFedListed = n) %>%
  left_join(FedPlantStateListed, by = 'state') %>%
  mutate(
    proportionListed = fedStateListed / totalFedListed,
    proportionNotListed = (totalFedListed - fedStateListed) /
      totalFedListed
  ) %>%
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>%
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed)) %>%
  mutate(percentProtected = proportionListed * 100) %>%
  mutate(percentProtected = if_else(
    percentProtected == 0,
    NA_integer_,
    as.integer(percentProtected)
  ))


fedPLantProp <- ggplot(FedPlantStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = percentProtected), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#78c679', high = '#004529',
                       midpoint = 50, na.value = 'white', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = '% of ESA Animal Species State Listed',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(7, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.09),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,0,2), "lines"))



### * IUCN maps ----------------------

### plot total IUCN imperiled species that occur in each state

iucnAllSp <- iucnSp %>% 
  group_by(state) %>% count() %>% dplyr::select(totalListedSpecies = n)


i <- ggplot(iucnAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalListedSpecies), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
   # scale_fill_gradientn2(low = '#fff2f2',
   #                     mid = '#c74242',
   #                     high = '#4f0702',
   #                     midpoint = 100, na.value = 'gray', limits = c(0,200),
   #                     breaks = c(0, 50, 100,150, 200),
   #                     name = 'Total IUCN Species',
   #                     guide = guide_colorbar(title.position = 'top',
   #                                            direction = 'horizontal', barwidth = unit(10, 'cm'),
   #                                            barheight = unit(0.3, 'cm'))) +
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       limits = c(0,234), breaks = c(0,25, 50,75, 100,125,150,175,200, 225), name = 'Total IUCN Species',
                       guide = guide_colorbar(title.position = 'top', direction = 'horizontal',
                                              barwidth = unit(10, 'cm'), barheight = unit(0.3, 'cm')))+
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

##separate by iucn imperiled wildlife and plants


### iucn wildlife

iucnWildlife <- iucnSp %>% mutate(kingdon = tolower(kingdom)) %>% 
  filter(kingdon != 'plantae') %>% group_by(state) %>% 
  count() %>% select(state, totalWildlife = n)

iw <- ggplot(iucnWildlife, aes(map_id = state))+
  geom_map(aes(fill = totalWildlife), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#fe9929', high = '#662506',
                       midpoint = 108, na.value = '#F5F3EA', limits = c(0,216),
                       breaks = c(0,50, 100, 150, 200),
                       name = 'IUCN Imperiled Animal Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(4.5, 'cm'),
                                              barheight = unit(1, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.88, 0.15),
        plot.background = element_rect(fill = "#F5F3EA", color = NA), 
        panel.background = element_rect(fill = "#F5F3EA", color = NA), 
        legend.background = element_rect(fill = "#F5F3EA", color = NA),
        legend.title = element_text(family = "sans", size = 9, face = 'bold', color = '#22211d'),
        legend.text = element_text(family = 'sans', size = 12, face = 'bold',
                                   color = '#22211d'))

## iucn plants

iucnPlants <- iucnSp %>% mutate(kingdom = tolower(kingdom)) %>% 
  filter(kingdom == 'plantae') %>% group_by(state) %>% 
  count() %>% select(state, totalPlants = n)


ip <- ggplot(iucnPlants, aes(map_id = state))+
  geom_map(aes(fill = totalPlants), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#78c679', high = '#004529',
                       midpoint = 20, na.value = '#0d3317', limits = c(0,40),
                       breaks = c(0, 10, 20, 30, 40),
                       name = 'IUCN Imperiled Plant Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(4.5, 'cm'),
                                              barheight = unit(1, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.88, 0.15),
        plot.background = element_rect(fill = "#F5F3EA", color = NA), 
        panel.background = element_rect(fill = "#F5F3EA", color = NA), 
        legend.background = element_rect(fill = "#F5F3EA", color = NA),
        legend.title = element_text(family = "sans", size = 9, face = 'bold', color = '#22211d'),
        legend.text = element_text(family = 'sans', size = 12, face = 'bold',
                                   color = '#22211d'))+
  annotate("text", x = -100, y = 24,
           label = 407,
           color = 'black', size = 4.5, fontface = 'bold')+
  annotate('segment', x = -101, y = 24, xend = -103, yend = 24, 
           arrow = arrow(length = unit(0.2, 'cm')))




## plot iucn imperiled species that are state-listed in each state

iucnState <- ggplot(iucnStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                      na.value = 'gray', limits = c(0,.6008),
                       breaks = c(0,.1, .2, .3, .4,.5, .6),
                       name = 'State Listed IUCN Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# plot iucn imperiled species that are ESA listed in each state

iucnFed <- ggplot(iucnFedListedSummary, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,.6),
                       breaks = c(0,.1, .2, .3, .4,.5, .6),
                       name = 'ESA Listed IUCN Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


#plot iucn imperiled species listed by either state or ESA (combined state and fed lists) in each state

iucnComb <- ggplot(iucnCombListedSummary, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,.6008),
                       breaks = c(0,.1, .2, .3, .4,.5, .6),
                       name = 'Combined State and ESA Listed IUCN Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))



# * Losing ESA protections -------------------------------------------

# fed current

fedProp <- ggplot(FedStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white', '#d2fbd4','#a5dbc2','#7bbcb0','#559c9e',
                                  '#3a7c89','#235d72','#123f5a', '#071822'),
                       na.value = 'gray', limits = c(0,1),
                       breaks = c(0,.2, .4, .6, .8, 1),
                       name = 'Proportion of ESA species State listed',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# fed proportion protected without ESA

fedNoESA <- ggplot(FedStateListedSummaryNoESA, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#d9d9d9", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white', '#d2fbd4','#a5dbc2','#7bbcb0','#559c9e',
                                  '#3a7c89','#235d72','#123f5a', '#071822'), limits = c(0,1),
                       breaks = c(0,.2, .4, .6, .8, 1),
                       name = 'Proportion of ESA species State listed',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

## total species losing protection

speciesLost <- FedStateListedSummaryNoESA %>% 
  mutate(speciesLost = totalFedListed - fedStateListed_noESA) %>% 
  mutate(proportionLost = speciesLost/totalFedListed*100)

fedLost <- ggplot(speciesLost, aes(map_id = state))+
  geom_map(aes(fill = speciesLost), map = fifty_states, color = "#ebeae8", size = 0.003)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white', '#d2fbd4','#a5dbc2','#7bbcb0','#559c9e',
                                  '#3a7c89','#235d72','#123f5a', '#071822'), na.value = 'gray', limits = c(0,130),
                       breaks = c(0, 25, 50, 75, 100, 125),
                       name = 'Number of Species Losing Protection',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# proportion of fed sp. losing protection

ggplot(speciesLost, aes(map_id = state))+
  geom_map(aes(fill = proportionLost), map = fifty_states, color = "#ebeae8", size = 0.003)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white', '#d2fbd4','#a5dbc2','#7bbcb0','#559c9e',
                                  '#3a7c89','#235d72','#123f5a', '#0d2e42'), na.value = 'gray', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = 'Proportion of Species Losing Protection',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))



#iucn proportion currently

iucnStateFig3 <- ggplot(iucnStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,.6008),
                       breaks = c(0,.1, .2, .3, .4,.5, .6),
                       name = 'State Listed IUCN Species',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

#iucn proportion without ESA

iucnNoESA <- ggplot(iucnStateListedSummaryNoESA, aes(map_id = state))+
  geom_map(aes(fill = proportionListed), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,.6008),
                       breaks = c(0,.1, .2, .3, .4,.5, .6),
                       name = 'Proportion of IUCN Species Listed Without the ESA',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


## total iucn species losing protection

speciesLostIUCN <- iucnStateListedSummaryNoESA %>% 
  dplyr::select(state, iucnStateListed_noESA) %>% 
  left_join(iucnCombListedSummary, by = 'state') %>% 
  mutate(speciesLost = iucnCombListed - iucnStateListed_noESA) %>% 
  mutate(proportionLost = speciesLost/iucnCombListed*100)


iucnLost <- ggplot(speciesLostIUCN, aes(map_id = state))+
  geom_map(aes(fill = speciesLost), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,90),
                       breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),
                       name = 'Number of IUCN Species Losing Protection',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# proportion IUCN species losing protection

ggplot(speciesLostIUCN, aes(map_id = state))+
  geom_map(aes(fill = proportionLost), map = fifty_states, color = "#A8A592", size = 0.3)+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradientn(colors = c('white','#FFE68E','#FFB45C','#ff0000',
                                  '#B80000','#850000','#390000'),
                       na.value = 'gray', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = 'Number of IUCN Species Losing Protection',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(10, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_blank(),
        legend.title = element_text(family = "sans", size = 11.5, face = 'bold', color = '#22211d',
                                    hjust = 1, debug = FALSE, margin = unit(c(0,0,0,0), "lines")),
        legend.text = element_text(family = 'sans', size = 10.5, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


