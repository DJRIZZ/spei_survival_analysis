# This script:
# (1) matches plastic band codes read from bands on breeding female Spectacled Eiders at Kigigak Island between 2019 and 2025
# with the number on their metal band (not readable in the field without capturing the bird) based on recent and past banding records.
# (2) uses that table to create an encounter history with 1 record for each band and 1 column for each year

# the output file is: "ec_na_1992-2015_ad_ducklings_2019-2025_adults.csv" and includes data from 2 periods, pre- and post-2015
# the pre-2015 data (1992-2015) includes individuals banded with plastic bands as both adults and ducklings
# the post-2015 data (2019, 2021-2025) includes only adults (only 25 female ducklins were banded with plastic bands in 2023, thus the resight period was short, only 2025)
# all birds marked with PTTs or TDRs are removed (TDR birds are included after TDR is retrived and removed)
# all birds first marked with yellow plastic bands in 2025 are excluded from the encounter history because they provide no data for the analysis given 2025 is the last encounter occasion
# We have limited data on banding before 2019, as the data shared by the Yukon Delta NWR is not complete (folders of data from each year 1994-2015;
# e.g., the 2014 folder has no MARK (capture and banding) data.

# We rely on those annual files in addition summary tables created by Yukon Delta NWR staff that summarize banding data,
# and our own (Endangered Species Recovery Program) data for bands deployed between 2019 and 2025.

# Got packages?
library(RODBC) # connect to MS Access data base
library(tidyverse)

# First, read-in pre- and post-2015 data from the Access database ####
# set path and file name for Kigigak SPEI data base that contains data on bands deployed and bands resighted 2019-2025
db <- "../spei_kigigak/data/database/db_kig_nesting_waterfowl_2019-2025_compiled_20250724.accdb"

# connect to the Access database via RODBC package
con_db <- odbcConnectAccess2007(db)

# get list of data tables available in the Access database
sqlTables(con_db,tableType="TABLE")$TABLE_NAME

# Load tables with resighted band codes 2019-2025 that need to be matched to their metal band numbers
nest <- sqlFetch(con_db, "tbl_Nest") # band codes associated with birds resighted at their nests
resight <- sqlFetch(con_db, "tbl_Resight") # band codes from birds seen away from nests, e.g., roosting flocks

# Load tables that contain banding data (plastic band codes and associated metal band numbers)
band <- sqlFetch(con_db, "tbl_Band") # banding data from birds captured 2019-2025, including recaptures AND new bands
old_bands <- sqlFetch(con_db, "tbl_spei_kig_bands_1992-2015") # bands deployed and resighted before 2019 when ES took over the project from file "ALL BANDS_Revised_November2015.xls"
odbcClose # close the connection the data base

# Compile band resights 2019-2025 ####
nest.p <- subset(nest, id_plasticBand != "na") # subset to nest records with plastic band codes, i.e., those not entered as NA
nest.p$id_year <- format(nest.p$dt_found, "%Y") # create capture year variable from the date the nest was found
#keep.nest.p <- c("id_yrNest", "dt_found", "id_plasticBand", "cat_resightMethod") # vector of relevant columns to keep
keep.nest.p <- c("id_plasticBand", "id_year") # vector of relevant columns to keep
nest.p <-nest.p[ ,keep.nest.p] # reduce to relevant columns
rm(keep.nest.p)
resight$id_year <- format(resight$dt_resight, "%Y") # for resights not associated with nests, create year variable from date of resight
keep.resight <- c("id_plasticBand", "id_year") # vector of relevant columns to keep
resight <- resight[ ,keep.resight] # reduce resight table to relevant columns
rm(keep.resight)
resight.all <- rbind(nest.p, resight) # combine nest resights and non-nest resights

# Compile banding data for the 2 periods: pre- and post-2015
# post-2015
bands.recent <- subset(band, id_plasticBand != "na") # subset to banding records that have plastic band codes (males were not banded with plastic bands)
bands.recent <- subset(bands.recent, cat_age != "local") # remove individs banded as ducklings; at Kig most were PTT in 2023; non-PTT locals from 2023 had 0 resights 2024-2025
#bands.recent <- subset(bands.recent, is_transmitter != "yes") # remove birds that were implanted with satellite transmitters
bands.recent$id_year <- format(bands.recent$dt_capt, "%Y") # create capture year variable
#plastic <- subset(plastic, cat_tdrStatus != "deployed") # keep TDR birds for now, can removed later - all TDR color bands began with "A"
#keep.band.p <- c("id_metalBand", "id_yrNest", "dt_capt", "is_recap", "cat_age", "cat_sex", "id_plasticBand", "is_bandReplaced", "id_oldBand") # relevant columns to keep
keep.bands.recent <- c("id_metalBand","id_plasticBand", "id_year") # relevant columns to keep
bands.recent <- bands.recent[ ,keep.bands.recent] # reduce to relevant columns
rm(keep.bands.recent)
# pre-2015
bands.past <- old_bands[!is.na(old_bands$id_plasticBand), ] # select records with plastic band codes
colnames(bands.past)[colnames(bands.past) == "id_yearRecap"] <- "id_year" # rename year recap to id_year
keep.old_bands <- c("id_metalBand", "id_year", "id_plasticBand") # relevant columns to keep
bands.past <- bands.past[ ,keep.old_bands] # reduce to relevant columns
rm(keep.old_bands)

# combine all resights pre- and post-2015
bands.all <- rbind(bands.recent, bands.past) # combine recent and past band data
ct_bands_all <- length(unique(bands.all$id_metalBand))
ct_bands_all
ct_resights <- table(bands.all$id_metalBand)
ct_resights

# Merge metal bands from bands.all into resight.all by id_plasticBand
bands_all_merged <- merge(resight.all, bands.all, by = "id_plasticBand", all.x = TRUE)

# these are the plastic bands resighted at nests 2019-2025 that don't have a known corresponding metal band numbers
nst_unkn_plastic_bands <- bands_all_merged[is.na(bands_all_merged$id_metalBand), ]
nst_unkn_plastic_bands # should be empty

# merge metal band number in to table of resights 2019-2025
bands_reftable <- cbind(bands_all_merged[1], bands_all_merged[3]) # reduce to just plastic and metal band info
bands_reftable <- unique(bands_reftable) # removed duplicate records related to resights

test.resight.all <- merge(resight.all, bands_reftable, by = "id_plasticBand", all.x = TRUE)

# compare resight.all with metal band
quick.compare <-transform(
                test.resight.all,
                multiple_products = +(ave(match(id_metalBand, unique(id_metalBand)), id_plasticBand, FUN = var) > 0
                  )
                    )
# add id_metalBand to resight.all
resight_mp <- merge(test.resight.all, bands_reftable, by = "id_plasticBand", all.x = TRUE)
resight_mp <- resight_mp[ ,1:3]
colnames(resight_mp) <- c("id_plasticBand", "id_year", "id_metalBand")

# rbind the updated resight.all with bands.recent
resight_tab <- rbind(resight_mp, bands.recent)
resight_tab <- unique(resight_tab) # removed duplicate records related to multiple resights within a single year

# vector of id_plasticBand for birds with implanted PTTs in 2018
ptt <- c("00V", "44G", "V37", "TKK", "V26", "V29", "V38", "V57", "V60", "V62", "V69", "V79", "V82", "V83", "V88", "V91", "V94", "V97", "V98" )

# add ducklings from Kig marked with PTTs
ptt <- c(ptt, "26C", "22-", "00-", "08-", "03-", "51-", "46-", "11-", "43-", "47-", "59-")

# remove PTT birds due to potential effects of implanted transmitter and the very low resight effort in 2018 (only arrival mist net captures)
resight_tab_noPTT <- resight_tab[!(resight_tab$id_plasticBand %in% ptt), ]

# how many banded birds resighted 2019-2015?
length(unique(resight_tab_noPTT$id_metalBand))

resight_tab <- resight_tab_noPTT

# vector of plastic band codes for birds with band-mounted TDRs (note: when TDRs are removed, these birds are banded with yellow plastic to enter the survival sample)
tdr <- c("A01", "A09", "A12", "A14", "A17", "A20", "A29", "A32", "A33", "A37", "A40", "A41", "A42", "A43", "A45", "A47", "A48", "A49", "A50", "A55", "A56", "A57", "A70", "A71", "A73", "A76", "A79", "A81", "A82", "A83", "A88", "A95")

# remove birds with TDRs on bands
resight_tab_noTDR <- resight_tab[!(resight_tab$id_plasticBand %in% tdr), ]

resight_tab <- resight_tab_noTDR

table(resight_tab$id_year)

# vector of plastic bands 2019-2025 that were put on ducklings, excluding ducklings implanted with PTTs, consider removing these or creating an indicator variable for age class
ducklings <- c("00C", "08C", "14-", "15C", "21C", "22-", "26C", "31C", "52C", "57C", "60C", "61C", "64-", "65C", "73C", "74C", "79C", "81C", "84C", "85-", "90-", "90C", "91C", "92C", "99C")
# NB only 25 female ducklings were marked with plastic bands in 2023 and none were resighted in 2025; not enough data to include ducklings

resight_tab_noDucklings <- resight_tab[!(resight_tab$id_plasticBand %in% ducklings), ]

resight_tab <- resight_tab_noDucklings

table(resight_tab$id_year)

# vector if individuals first marked in the last occasion (these provide no data for the analysis, can remove)
first_last <- c(250707630, 250707633, 200762260, 200762347) # this list excludes individs first marked in 2025 with TDR bands, which were removed above (A29, A82, A76, A57); note 200762260 & 200762347 were a TDR birds that had their TDRs removed in 2025 and were re-banded yellow at that time

# remove individs first marked during the last occasion
resight_tab_noFirstLast <- resight_tab[!(resight_tab$id_metalBand %in% first_last), ]

resight_tab <- resight_tab_noFirstLast

# make encounter history with counts per year per band and then transform to 1's and 0's
# define sessions (years)
sessions <- sort(unique(resight_tab$id_year))

# create id × session data frame
wide <- resight_tab %>%
  mutate(present = 1L) %>%
  distinct(id_metalBand, id_year, .keep_all = TRUE) %>%
  arrange(id_year) %>% 
  pivot_wider(id_cols = id_metalBand,
              names_from = id_year,
              names_prefix = "s",
              values_from = present,
              values_fill = 0L)

# Collapse session columns into a single capture-history string
#wide <- wide %>%
#  mutate(ch = pmap_chr(select(., starts_with("s")), ~paste0(c(...), collapse = ""))) %>%
#  select(id, ch)

# combine the 2019, 2-21-2025 encounter history with the 1992-2015 encounter history with year columns with no resight effort (2016, 2017, 2018, and 2020) filled with NAs

# load encounter history for 1992-2015
eh1_all <- read.csv("data/data_christie_paper/Encounter_hist_SPEI.csv")

# change column names
past.years <- seq(1992, 2015, 1)
past.years <- paste0("s", past.years)
eh1.names <- c("id_metalBand", past.years, "is_duckling", "is_adult")
colnames(eh1_all) <- eh1.names

# for now, drop the age columns in eh2
eh1 <- eh1_all[ ,1:25]

# encounter history for 2019, 2021-2025
eh2 <- wide

# define years and effort
all_years <- 1992:2025
no_effort_years <- c(2015:2018, 2020)
effort_years <- setdiff(all_years, no_effort_years)
all_cols <- paste0("s", all_years)

# Fill in missing years and columns in each encounter history
fill_missing_years <- function(df, all_cols) {
  missing_cols <- setdiff(all_cols, names(df))
  df[missing_cols] <- NA_integer_
  df <- df[, c("id_metalBand", all_cols)]
  return(df)
}

eh1_full <- fill_missing_years(eh1, all_cols)
eh2_full <- fill_missing_years(eh2, all_cols)

# merge the two data frames by id
combined <- full_join(eh1_full, eh2_full, by = "id_metalBand", suffix = c(".x", ".y"))

# coalesce values for each year and fill zeros where appropriate
for (year in all_years) {
  col_x <- paste0("s", year, ".x")
  col_y <- paste0("s", year, ".y")
  combined[[paste0("s", year)]] <- dplyr::coalesce(combined[[col_x]], combined[[col_y]])
}

# remove intermediate .x and .y columns
combined <- combined[, c("id_metalBand", paste0("s", all_years))]

# replace NAs with 0 for years with effort, leave NA for no-effort years
for (year in effort_years) {
  col <- paste0("s", year)
  combined[[col]][is.na(combined[[col]])] <- 0L
}
# no-effort years remain NA

# save as csv
write.csv(combined, "output/ec_na_1992-2015_ad_ducklings_2019-2025_adults.csv", row.names = FALSE)
write.csv(combined, "data/ec_na_1992-2015_ad_ducklings_2019-2025_adults.csv", row.names = FALSE)
