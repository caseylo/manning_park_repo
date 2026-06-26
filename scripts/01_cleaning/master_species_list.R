#### Creating a master species list

library(tidyverse)
library(dplyr)

BEC <- read_csv("data/raw/BEC_species_code_and_latin_names.csv")

data <- read_csv("data/raw/species_code_and_latin_names_2025.csv")


####

## Checking for species codes not in BEC that are in data so I can add them into a master list
not_in_BEC <- data %>%
  anti_join(BEC, by = "BC_species_code")

not_in_BEC

species_lookup <- bind_rows(
  BEC %>% select(BC_species_code, Latin_name, Functional_group),
  data %>% select(BC_species_code, Latin_name, Functional_group)
) %>%
  distinct(BC_species_code, .keep_all = TRUE)

species_lookup

species_lookup_final <- species_lookup %>%
  mutate(
    BC_species_code = recode(
      BC_species_code,
      "ABIELASD"   = "ABIELAS",
      "ACERGLA1"    = "ACERGLA",
      "ARNICA"      = "ARNICA SP.",
      "ACHIMIL2"    = "ACHIMIL",
      "ALNUALN2"    = "ALNUALN",
      "ANTEALPI"   = "ANTEALP",
      "ARENARIA"   = "ARENARIA SP.",
      "ANTENARIA"  = "ANTENARIA SP.",
      "ASTER"      = "ASTER SP.",
      "BRACHYT"    = "BRACHYTHECIUM SP.",
      "CALAMAGR"   = "CALAMAGROTIS SP.",
      "CAREX"      = "CAREX SP.",
      "CEPHALOZ"   = "CEPHALOZIA SP.",
      "CIRSIUM"    = "CIRSIUM SP.",
      "CLADONIA"   = "CLADONIA SP.",
      "CORNUNA"    = "CORNCAN",
      "DICRANUM"   = "DICRANUM SP.",
      "DRABA"      = "DRABA SP.",
      "EQUISETU"   = "EQUISETUM SP.",
      "ERIGERON"   = "ERIGERON SP.",
      "FESTUCA"    = "FESTUCA SP.",
      "HIERACIU"   = "HIERACIUM SP.",
      "LETHARIA"   = "LETHARIA SP.",
      "LICHEN"     = "LICHEN SP.",
      "LUPILEP2"   = "LUPILEP",
      "LOPHOZIA"   = "LOPHOZIA SP.",
      "MITELLA"    = "MITELLA SP.",
      "PICEENGD"   = "PICEENG",
      "PICEENED"   = "PICEENG",
      "PINUCON1"   = "PINUCON",
      "PINUCON2"   = "PINUCON",
      "PLAGIOMN"   = "PLAGIOMNIUM SP.",
      "POA_CUS2"   = "POA_CUS",
      "POACEAE_"    = "POACEAE SP.",
      "POLYGONU"   = "POLYGONUM SP.",
      "PSEUMEN1"   = "PSEUMEN",
      "PSEUMEN2"   = "PSEUMEN",
      "PSEUMEND"   = "PSEUMEN",
      "POLYTRIC"    = "POLYTRICHUM SP.",
      "POTENTIL"   = "POTENTILLA SP.",
      "POTEUNIF"   = "POTEUNI",
      "RACOMITR"   = "RACOMITRIUM SP.",
      "PETAFRIG"   = "PETAFRI",
      "PELTIGER"   = "PELTIGERA SP.",
      "ROSA"       = "ROSA SP.",
      "RHYTIDIO"   = "RHYTIDIOPSIS SP.",
      "SALIX"      = "SALIX SP.",
      "SENECIO"    = "SENECIO SP.",
      "LINOINV"    = "LONIINV",
      "SABURUB"    = "SAMBRAC",
      "POHLIA"     = "POHLIA SP.",
      "TIARTRI2"   = "TIARTRI",
      "VIOLA"      = "VIOLA SP."  
    )
  )

species_lookup_final

## SAVE
write_csv(species_lookup_final,"data/cleaned/master_species_list.csv")
