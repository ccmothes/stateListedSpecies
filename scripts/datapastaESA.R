# using datapasta to read in esa species lists

library(tidyverse)
library(datapasta)
library(stringr)

# alabama

# alaska

ak <- tibble::tribble(
  ~ Status,                                                       ~Species,
    
  "E", "Albatross, short-tailed Wherever found (Phoebastria (=Diomedea) albatrus)",
  "T",                              "Bear, polar Wherever found (Ursus maritimus)",
  "T",                       "Bison, wood Wherever found (Bison bison athabascae)",
  "E",                         "Curlew, Eskimo Wherever found (Numenius borealis)",
  "T",                     "Eider, spectacled Wherever found (Somateria fischeri)",
  "T",                   "Eider, Steller's AK breeding pop. (Polysticta stelleri)",
  "T",         "Otter, Northern Sea Southwest Alaska DPS (Enhydra lutris kenyoni)",
  "E",           "Whale, sperm Wherever found (Physeter catodon (=macrocephalus))"
  
) %>%
  rbind(
    tibble::tribble(
      ~ Status,
      ~ Species,
      "E",
      "Fern, Aleutian shield (Polystichum aleuticum)"
    )
  ) %>% 
  mutate(Species = str_extract_all(Species, "\\(.*")) %>% 
  mutate(Species = str_replace_all(Species, "[()]", "")) %>% 
  mutate(Species = str_replace_all(ak$Species, "=\\w+ *", "")
)


CT <- tibble::tribble(
        ~Status,                                                                                                                                         ~Species/Listing.Name,
            "E",                                                                                                                "Bat, Indiana Wherever found (Myotis sodalis)",
            "T",                                                                                            "Bat, Northern long-eared Wherever found (Myotis septentrionalis)",
            "T",                                                                                                            "Knot, red Wherever found (Calidris canutus rufa)",
            "T", "Plover, piping [Atlantic Coast and Northern Great Plains populations] - Wherever found, except those areas where listed as endangered. (Charadrius melodus)",
            "E",                                                                                               "Sea turtle, hawksbill Wherever found (Eretmochelys imbricata)",
            "E",                                                                                               "Sea turtle, leatherback Wherever found (Dermochelys coriacea)",
            "E",                                                                                "Tern, roseate Northeast U.S. nesting population (Sterna dougallii dougallii)",
            "T",                                                                                                   "Tiger beetle, Puritan Wherever found (Cicindela puritana)",
            "T",                                                                                "turtle, bog Wherever found, except GA, NC, SC, TN, VA (Clemmys muhlenbergii)",
            "E",                                                                                                   "Wedgemussel, dwarf Wherever found (Alasmidonta heterodon)"
        )


