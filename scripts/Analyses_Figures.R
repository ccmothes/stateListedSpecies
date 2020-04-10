
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
iucnSp <- read_csv('data/iucnFINAL.csv')
fedSp <- read_csv('data/fedFINAL.csv')
stateImp <- read_csv('data/stateFINAL.csv')

# make summary files (# sp per state)
stateAllSp <- stateImp %>% add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% select(state, totalProtectedSpecies = n) %>% ungroup() %>% 
  mutate(totalProtectedSpecies = if_else(totalProtectedSpecies == 1, NA_integer_,totalProtectedSpecies))

fedAllSp <- fedSp %>% 
  group_by(state) %>% count() %>% select(state, totalListedSpecies = n)

iucnAllSp <- iucnSp %>% 
  group_by(state) %>% count() %>% select(state, totalListedSpecies = n)




# ANALYSES ----------------------------------------------------------------


# * IUCN Coverage ---------------------------------------------------------

##get category info for iucn species

###NOTE there were about 48 species where the species as a whole was not imperiled,
###BUT it had either a sub species or sub population that was imperiled, so we counted that


#overall # of imperiled species that are not protected by any state
iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName))

#total not listed = 850 / 1625 total imperiled species
# 52% imperiled species not protected by any state in the US

# separate analyses by phylum
iucnPhyla <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(phylum = tolower(phylum)) %>% 
  group_by(phylum) %>% count() %>% rename(totalIUCNsp = n)

iucnStatePhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)

# separate by class
iucnClass <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(totalIUCNsp = n)

iucnStateClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp)


#Now compare, for iucn state listed species, the number of states they live in compared to 
#the number of states that list them

iucnByState <- iucnSp %>% filter(acceptedName %in% stateImp$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

imperiledPerState <- stateImp %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesListedIn = n) %>% left_join(iucnByState, by = 'acceptedName') %>% 
  mutate(proportionListed = statesListedIn/statesOccurIn)

#1/3 are not protected by all the states they occur in, making total number of sp
#not protected in their entire U.S. range = 68%.
# Of the sp. that are protected in their entire range, 93% are state endemics


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
## The average proportion of iucn species not listed by states is 75% !!


#now reverse, calculate how many state listed species are not on IUCN list

stateNotIUCN <- stateImp %>% group_by(state) %>% filter(!(acceptedName %in% iucnSp$acceptedName))


# check if they are on the iucn list at all

usIUCN <- read_csv('data/US_IUCN_acceptedName.csv')

notIUCNlisted <- stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% group_by(acceptedName) %>% 
  count() %>% rename(notIUCNlisted = n)

stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  group_by(taxon_general) %>%  count() %>% rename(notIUCNlisted = n)
#4841 species total, 3447 are plants:
# 1 Amphibian                 2
# 2 Amphibians               74
# 3 Birds                   183
# 4 Fish                    256
# 5 Fungi                    20
# 6 Invertebrates           631
# 7 Mammals                 113
# 8 Plants                 3447
# 9 Reptiles                115

notIUCNlisted %>% filter(!(acceptedName %in% usIUCN$acceptedName))
##only 32% of the missing species are even assessed by IUCN,
##3316/4841 not assessed

notIUCNphyla <- stateImp %>% filter(!(acceptedName %in% iucnSp$acceptedName)) %>% 
  filter(!(acceptedName %in% usIUCN$acceptedName)) %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  group_by(taxon_general) %>% count()
# 2825 of these species (85%) that are listed by states but not IUCN assessed are PLANTS


# * Federal Coverage ------------------------------------------------

fedSp %>% filter(!(acceptedName %in% stateImp$acceptedName)) %>% distinct(acceptedName)
# 217 species not on any state list (16%)

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



# New Results if ESA was removed ---------------------------

## read in legislation file

stateLeg <- read_csv('data/stateLegislationSummary.csv') %>% dplyr::select(state:last_updated) %>% 
  mutate(state = tolower(state))

## remove esa species from states that automatically include them
st <- stateLeg %>% filter(ESA_includes == 'Y') %>% pull(state)
  # 34 states that include 'ESA - listed species' to define state listed species,
 # but 3 of these states don't have any imperiled sp laws (i.e. they define the species, but
# don't have any legislative provisions protecting them)

ESAstates <- stateImp %>% filter(state %in% st) %>% group_by(state) %>% 
  filter(!(acceptedName %in% fedSp$acceptedName)) %>% ungroup()

NonESAstates <- stateImp %>% filter(!(state %in% st))

# add back FT species for states that don't automatically include FT (i.e. they listed them independently)

statesNoT <- c('south carolina', 'texas', 'mississippi', 'kentucky')

esaThreatenedMissing <- map_dfr(statesNoT,
                                function(x) {
                                  fsp <- fedSp %>%
                                    filter(state == x & status == 'T')
                                  stateImp %>%
                                    filter(state == x & acceptedName %in% fsp$acceptedName)
                                })


stateImp_noESA <- rbind(ESAstates, NonESAstates, esaThreatenedMissing)


# * * IUCN coverage NEW --------------------------------------------------------------------------

#overall # of imperiled species that are not protected by any state
iucnSpNoESA <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp_noESA$acceptedName))

#total not listed = 850 -> 1326 / 1625 total imperiled species = 30% increase in unprotected sp
# 52% -> 82%  imperiled species not protected by any state in the US


#Now compare, for iucn state listed species, the number of states they live in compared to 
#the number of states that list them
iucnByStateNoESA <- iucnSp %>% filter(acceptedName %in% stateImp_noESA$acceptedName) %>% 
  group_by(acceptedName) %>% count() %>% rename(statesOccurIn = n)

imperiledPerState <- stateImp_noESA %>% filter(acceptedName %in% iucnSp$acceptedName) %>% group_by(acceptedName) %>% 
  count() %>% rename(statesListedIn = n) %>% left_join(iucnByStateNoESA, by = 'acceptedName') %>% 
  mutate(proportionListed = statesListedIn/statesOccurIn)
##main result: of all imperiled species listed by at least one state, 1/3 -> 67% of them are not protect by
##all the states they occur in


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
  rename(iucnStateListed = n) 

iucnStateListedSummaryNoESA <-  iucnSp %>% group_by(state) %>% count() %>% rename(totalImperiled = n) %>% 
  left_join(iucnStateListedNoESA, by = 'state') %>% mutate(proportionListed = iucnStateListed/totalImperiled,
                                                      proportionNotListed = (totalImperiled - iucnStateListed)/totalImperiled) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))

#RESULTS
## now highest is VA with 44% protected
## now no states protect more than half of the imperiled species in their borders
## The average proportion of iucn species not listed by states is 75% -> now 87%!

## * * Fed coverage new --------------------------------------------------------------------


fedSp %>% filter(!(acceptedName %in% stateImp_noESA$acceptedName)) %>%
  distinct(acceptedName, .keep_all = TRUE)

#1079 / 1386 = 78% 
FedStateListedNoESA <- map_dfr(unique(fedSp$state),
                          function(x) {
                            st <- stateImp_noESA %>%
                              filter(state == x)
                            fedSp %>%
                              filter(state == x & acceptedName %in% st$acceptedName)
                          }) %>% 
  group_by(state) %>% 
  count() %>% 
  rename(fedStateListed = n) 

FedStateListedSummaryNoESA <-  fedSp %>% group_by(state) %>% count() %>% rename(totalFedListed = n) %>% 
  left_join(FedStateListedNoESA, by = 'state') %>% mutate(proportionListed = fedStateListed/totalFedListed,
                                                     proportionNotListed = (totalFedListed - fedStateListed)/totalFedListed) %>% 
  mutate(proportionNotListed = if_else(is.na(proportionNotListed), 1, proportionNotListed)) %>% 
  mutate(proportionListed = if_else(is.na(proportionListed), 0, proportionListed))

# now average proportion of ESA species not listed in each state is 86% (b/c lots of states are 0 now)

# * Correlations -------------------------------------------------------

## load file with state predictor variables
stateInfo <- read_csv('data/state_info.csv') %>% mutate(state = tolower(state))

#now join this and iucn summary (to get total # iucn imperiled species) with state summary

stateSummary <- stateAllSp %>% left_join(stateInfo, by = 'state') %>% left_join(iucnAllSp, by = 'state') %>% 
  dplyr::select(state:endemic_sp, totalImperiled = totalListedSpecies) %>% 
  mutate(totalProtectedSpecies = if_else(is.na(totalProtectedSpecies), as.integer(0), totalProtectedSpecies))

## create correlation matrix to remove highly correlated variables
mat <- stateSummary %>% select(totalProtectedSpecies, party_trend_rank, area_sqkm, 
                               diversity_rank, total_sp, risk_rank, endemism_rank, endemic_sp, totalImperiled)

rcorr(as.matrix(mat), type = 'spearman')
# remove endemic sp. and proportion at risk ('risk_rank') due to high correlation (r > 0.7)
# with total species (i.e. biodiversity)



## ** model selection --------------------------------------------------------------------------------------------

### create refined set of vars by removing correlated ones

stateSumRefine <- stateSummary %>% dplyr::select(totalProtectedSpecies, party_trend_rank, area_sqkm, total_sp,
                                          totalImperiled)

### create full model

full.mod <- lm(totalProtectedSpecies ~ party_trend_rank + area_sqkm + total_sp + totalImperiled, 
               data = stateSumRefine)

  
#now run model selection based on AIC
stepMod <- stepAIC(full.mod, direction = 'both')
stepMod$anova
## final model includes party trend and total imperiled sp


finalMod <- lm(totalProtectedSpecies ~ party_trend_rank + totalImperiled, 
               data = stateSumRefine)
anova(finalMod)

# test if corr still significant without Hawaii (outlier)
stateSumNoHI <- stateSummary %>% filter(state != 'hawaii') 
anova(lm(totalProtectedSpecies ~ party_trend_rank + totalImperiled, data = stateSumNoHI))

## now separate by wildlife and plants and re-run model selection

### wildlife
stateWildlife <- stateImp %>% filter(taxon_general != 'Plants') %>%
  add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% dplyr::select(state, totalWildlife = n) %>% ungroup() %>% 
  mutate(totalWildlife = if_else(totalWildlife == 1, NA_integer_,totalWildlife)) %>% left_join(stateInfo, by = 'state') %>% left_join(iucnAllSp, by = 'state') %>% 
  dplyr::select(state:endemic_sp, totalImperiled = totalListedSpecies) %>% 
  mutate(totalWildlife = if_else(is.na(totalWildlife), as.integer(0), totalWildlife))

wildMod <- lm(totalWildlife ~ party_trend_rank + area_sqkm + total_sp 
              + totalImperiled, data = stateWildlife)

stepModWild <- stepAIC(wildMod, direction = 'both')
stepModWild$anova
## wildlife model the same as full model with the addition of total species as
## an important predictor
wildlFinalMod <- lm(totalWildlife ~ party_trend_rank + total_sp + totalImperiled,
                    data = stateWildlife)
anova(wildlFinalMod)

# plants

stateWithPlants <- stateImp %>% filter(taxon_general == 'Plants') %>% pull(state) %>% unique()
stateNoPlant <- stateImp %>% filter(!(state %in% stateWithPlants)) %>% pull(state) %>% unique()

statePlants <- stateImp %>% 
  filter(taxon_general == 'Plants') %>%
  add_row(state = c(stateNoPlant,'west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>%
  count() %>% ungroup() %>% 
  dplyr::select(state, totalPlants = n) %>% 
  mutate(totalPlants = if_else(totalPlants == 1, NA_integer_,totalPlants)) %>% 
  left_join(stateInfo, by = 'state') %>% left_join(iucnAllSp, by = 'state') %>% 
  dplyr::select(state:endemic_sp, totalImperiled = totalListedSpecies) %>% 
  mutate(totalPlants = if_else(is.na(totalPlants), as.integer(0), totalPlants))

plantMod <- lm(totalPlants ~ party_trend_rank + area_sqkm + total_sp +
                   totalImperiled, data = statePlants)

stepModPlant <- stepAIC(plantMod, direction = 'both')
stepModPlant$anova
# plant model results the same as full model


## ** gtest -------------------------------------------------------------------------------


## iucn coverage gtest
iucnPhyla <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(phylum = tolower(phylum)) %>% 
  group_by(phylum) %>% count() %>% rename(totalIUCNsp = n)

iucnStatePhylaSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% 
  mutate(phylum = tolower(phylum)) %>% group_by(phylum) %>% count() %>% 
  rename(notListed = n) %>% left_join(iucnPhyla, by = 'phylum') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp, listed = totalIUCNsp - notListed) %>% 
  filter(totalIUCNsp >= 10) #removed those with less than 10 sp listed

iucnMat <- iucnStatePhylaSummary %>% select(notListed, listed) %>% 
  column_to_rownames(var = 'phylum') %>% 
  data.matrix()

GTest(iucnMat)
#p-value < 2.2e-16

#run pairwise tests
FUN = function(i,j){
  GTest(matrix(c(iucnMat[i,1], iucnMat[i,2],
                 iucnMat[j,1], iucnMat[j,2]),
               nrow = 2,
               byrow = TRUE),
          correct = 'none')$ p.value
}

phyla.gtest.table <- pairwise.table(FUN,
               rownames(iucnMat),
               p.adjust.method = 'none')
#                 arthropoda basidiomycota     chordata  cnidaria     mollusca
# basidiomycota 1.426100e-02            NA           NA        NA           NA
# chordata      0.000000e+00  2.042432e-07           NA        NA           NA
# cnidaria      5.796136e-03  2.133429e-04 6.710513e-01        NA           NA
# mollusca      3.526877e-08  8.680616e-05 4.046280e-07 0.3496248           NA
# tracheophyta  0.000000e+00  5.788582e-09 6.484118e-03 0.2872709 6.117329e-14

## Break down chordates by class
iucnClass <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% mutate(class = tolower(class)) %>% 
  group_by(class) %>% count() %>% rename(totalIUCNsp = n)
iucnStateClassSummary <- iucnSp %>% distinct(acceptedName, .keep_all = TRUE) %>% 
  filter(!(acceptedName %in% stateImp$acceptedName)) %>% mutate(phylum = tolower(phylum),
                                                                class = tolower(class)) %>% 
  filter(phylum == 'chordata') %>% 
  group_by(class) %>% count() %>% rename(notListed = n) %>% left_join(iucnClass, by = 'class') %>% 
  mutate(proportionNotListed = notListed/totalIUCNsp, listed = totalIUCNsp - notListed) %>% 
  filter(totalIUCNsp >= 10)

iucnMatChor <- iucnStateClassSummary %>% select(notListed, listed) %>% 
  column_to_rownames(var = 'class') %>% 
  data.matrix()

GTest(iucnMatChor)
#p-value = 7.457e-12

#pairwise tests
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



# VIZUALIZATION -------------------------------------------------------------

## ** correlation figures -------------

### protected sp vs. imperiled sp

corplot1 <- ggscatter(stateSummary, x = "totalImperiled", y = "totalProtectedSpecies",
                      size = 7, alpha = 0.4, color = 'black',
                      add = "reg.line", conf.int = TRUE,
                      
                      xlab = "IUCN Imperiled Species", ylab = "State Listed Species",
                      font.x = c(18,'bold', 'black'), ylim = c(0, 801), font.y = c(18, 'bold', 'black'),
                      add.params = list(color = "black",
                                        fill = "lightgray"))+
  
  theme(plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 400, label.y = 810, size = 5
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

# protected sp vs. political party ranking
corplot2 <- ggscatter(stateSummary, x = "party_trend_rank", y = "totalProtectedSpecies",
                      color = 'party_trend_rank', size = 7,
                      add = "reg.line", conf.int = TRUE,
                      
                      xlab = "", ylab = "State Listed Species",
                      add.params = list(color = "black",
                                        fill = "lightgray"),
                      ylim = c(-40, 801),
                      font.y = c(18, 'bold', 'black'),
                      xlim = c(-24, 25))+
  gradient_color(c('blue', 'purple', 'red'))+
  theme(legend.position = 'none', axis.text.x = element_blank(),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA),
        axis.ticks.x = element_blank())+
  annotate("text", x = c(-15, 15), y = c(-40, -40),
           label = c('Most Democratic State', 'Most Republican State'),
           color = 'black', size = 6, fontface = 'bold')+
  annotate('segment', x = -5, y = -35, xend = 5, yend = -35,
           colour = 'black', size = 1, arrow = arrow(length = unit(0.3,"cm"),
                                                     type = 'closed'))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 10, label.y = 810, size = 5
  )

# * * Figure 2 -------------------------------------------------------------------------
ggarrange(corplot2, corplot1, ncol = 1, labels = c('a)', 'b)'), 
          label.y = 1.02,label.x = -0.01,
          font.label = list(size = 22))



# * state maps ----------------------------------------


# * * All species -------------------------------------------------------------


### get total num species listed per state

stateAllSp <- stateImp %>% add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% select(state, totalProtectedSpecies = n) %>% ungroup() %>% 
  mutate(totalProtectedSpecies = if_else(totalProtectedSpecies == 1, NA_integer_,totalProtectedSpecies))


# make plot
s1 <- ggplot(stateAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalProtectedSpecies), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ECEBF7' , mid = '#748DC8', high = '#04013F',
                       midpoint = 400, na.value = 'white', limits = c(0,801),
                      breaks = c(0,200, 400, 600, 800),
                      name = 'Total State Listed Species',
                      guide = guide_colorbar(title.position = 'top',
                                             direction = 'horizontal', barwidth = unit(10, 'cm'),
                                             barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.05),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(family = "sans", size = 11, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# * * Split state listed wildlife and plants -----------------------------------------------



## state wildlife

stateWildlife <- stateImp %>% filter(taxon_general != 'Plants') %>%
  add_row(state = c('west virginia', 'wyoming', 'utah', 'north dakota', 'alabama')) %>% 
  group_by(state) %>% count() %>% dplyr::select(totalWildlife = n) %>% ungroup() %>% 
  mutate(totalWildlife = if_else(totalWildlife == 1, NA_integer_,totalWildlife))


sw <- ggplot(stateWildlife, aes(map_id = state))+
  geom_map(aes(fill = totalWildlife), map = fifty_states, color = "#A8A592")+
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
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(family = "sans", size = 11, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))


## state plant species
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
  geom_map(aes(fill = totalPlants), map = fifty_states, color = "#A8A592")+
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
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(family = "sans", size = 11, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 10, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,2.5,0), "lines"))

# * * Figure 1 ----------------

ggarrange(s1, sw, sp, ncol = 1, labels = c('a)', 'b)', 'c)'), hjust = -8, vjust = 3)



### fed map

### get total num species listed per state

fedAllSp <- fedSp %>% 
  group_by(state) %>% count() %>% select(state, totalListedSpecies = n)


f <- ggplot(fedAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalListedSpecies), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ECEBF7' , mid = '#748DC8', high = '#04013F', 
                       midpoint = 400, na.value = 'white', limits = c(0,801),
                      breaks = c(0,200,400,600,800),
                      name = 'Federally Listed Species',
                      guide = guide_colorbar(title.position = 'top',
                                             direction = 'horizontal', barwidth = unit(4.5, 'cm'),
                                             barheight = unit(1, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.88, 0.15),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(family = "sans", size = 12, face = 'bold', color = '#22211d'),
        legend.text = element_text(family = 'sans', size = 13, face = 'bold',
                                   color = '#22211d'))

### iucn map

### get total imperiled species per state

iucnAllSp <- iucnSp %>% 
  group_by(state) %>% count() %>% select(totalListedSpecies = n)


i <- ggplot(iucnAllSp, aes(map_id = state))+
  geom_map(aes(fill = totalListedSpecies), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ECEBF7' , mid = '#748DC8', high = '#04013F', 
                       midpoint = 400, na.value = '#F5F3EA', limits = c(0,801),
                      breaks = c(0,200,400,600,800),
                      name = 'IUCN Imperiled Species',
                      guide = guide_colorbar(title.position = 'top',
                                             direction = 'horizontal', barwidth = unit(4.5, 'cm'),
                                             barheight = unit(1, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.88, 0.15),
        plot.background = element_rect(fill = "#F5F3EA", color = NA), 
        panel.background = element_rect(fill = "#F5F3EA", color = NA), 
        legend.background = element_rect(fill = "#F5F3EA", color = NA),
        legend.title = element_text(family = "sans", size = 12, face = 'bold', color = '#22211d'),
        legend.text = element_text(family = 'sans', size = 13, face = 'bold',
                                   color = '#22211d'))


## fed wildilfe and plants


### wildlife

fedWildlife <- fedSp %>% filter(taxon_general != 'Plants') %>%
  group_by(state) %>% count() %>% select(state, totalWildlife = n)

fw <-  ggplot(fedWildlife, aes(map_id = state))+
  geom_map(aes(fill = totalWildlife), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#fe9929', high = '#662506',
                       midpoint = 54, na.value = '#F5F3EA', limits = c(0,107),
                       breaks = c(0, 25, 50, 75, 100),
                       name = 'Federally Listed Animal Species',
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
### fed plants

fedPlants <- fedSp %>% filter(taxon_general == 'Plants') %>% 
  group_by(state) %>% count() %>% select(state, totalPlants = n)


fp <- ggplot(fedPlants, aes(map_id = state))+
  geom_map(aes(fill = totalPlants), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#ffffe5' , mid = '#78c679', high = '#004529',
                       midpoint = 31, na.value = '#0d3317', limits = c(0,62),
                       breaks = c(0, 15, 30, 45, 60),
                       name = 'Federally Listed Plant Species',
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
  annotate("text", x = c(-100,-124), y = c(24,35),
           label = c(410, 171),
           color = 'black', size = 4.5, fontface = 'bold')+
  annotate('segment', x = -101, y = 24, xend = -103, yend = 24, 
           arrow = arrow(length = unit(0.2, 'cm')))+
  annotate('segment', x = -123, y = 35, xend = -121, yend = 35,
           arrow = arrow(length = unit(0.2, 'cm')))


##iucn wildlife and plants


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






# * * proportion maps ---------------------------------------------------------

# iucn Map
iucnStateListedSummary <- iucnStateListedSummary %>% mutate(percentProtected = proportionListed*100) %>% 
 mutate(percentProtected = if_else(percentProtected == 0, NA_integer_, as.integer(percentProtected)))
iucnProp <- ggplot(iucnStateListedSummary, aes(map_id = state)) +
  geom_map(aes(fill = percentProtected), map = fifty_states, color = "#A8A592") +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map() +
  theme_void() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_gradient2(
    low = '#DFE4E7' ,
    mid = '#6088A1',
    high = '#02273E',
    midpoint = 30,
    na.value = 'white',
    limits = c(0, 60),
    breaks = c(0, 15, 30, 45, 60),
    name = '% of Imperiled Species Protected',
    guide = guide_colorbar(
      title.position = 'top',
      direction = 'horizontal',
      barwidth = unit(7, 'cm'),
      barheight = unit(0.3, 'cm')
    )
  ) +
  labs(x = "", y = "") +
  theme(
    legend.position = c(0.5,-0.09),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(
      family = "sans",
      size = 11,
      face = 'bold',
      color = '#22211d',
      hjust = 1
    ),
    legend.text = element_text(
      family = 'sans',
      size = 10,
      face = 'bold',
      color = '#22211d'
    ),
    plot.margin = unit(c(0,0,0,2), "lines")
  )



## fed porportion map

fedStateListedSummary <- FedStateListedSummary %>% mutate(percentProtected = proportionListed*100) %>% 
  mutate(percentProtected = if_else(percentProtected == 0, NA_integer_, as.integer(percentProtected)))
fedProp <- ggplot(fedStateListedSummary, aes(map_id = state))+
  geom_map(aes(fill = percentProtected), map = fifty_states, color = "#A8A592")+
  expand_limits(x = fifty_states$long, y = fifty_states$lat)+
  coord_map()+
  theme_void()+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  scale_fill_gradient2(low = '#DFE4E7' , mid = '#6088A1', high = '#02273E',
                       midpoint = 50, na.value = 'white', limits = c(0,100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = '% of ESA Species Protected',
                       guide = guide_colorbar(title.position = 'top',
                                              direction = 'horizontal', barwidth = unit(7, 'cm'),
                                              barheight = unit(0.3, 'cm'))) +
  labs(x = "", y = "")+
  theme(legend.position = c(0.5, -0.09),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA), 
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(family = "sans", size = 11, face = 'bold', color = '#22211d',
                                    hjust = 1),
        legend.text = element_text(family = 'sans', size = 9, face = 'bold',
                                   color = '#22211d'),
        plot.margin = unit(c(0,0,0,2), "lines"))


# * * Figure 3: -------------------------------------------------------------------------

ggarrange(iucnProp, fedProp, nrow= 1, labels = c('a)', 'b)'), vjust = 33, hjust = -1.5)

