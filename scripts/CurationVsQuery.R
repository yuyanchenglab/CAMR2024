#!/bin/env Rscript3

# JM 7/31/2024
# How do my statistics methods compare to a literature search when looking for proper cell markers?

# Setup ----

library(dplyr) # pipe
library(data.table) # fast table manip
library(ggplot2)
library(stringr)

cheng_theme = new.env() # Note: Use dollar sign operator to access values
source("/project/hipaa_ycheng11lab/software/apps/github/figures/theme.cheng.module.R", local = cheng_theme)

analysisName = "CurationVsQuery"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)

curatedPath = "spreadsheets/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "spreadsheets/normalized_majorclass_mean_expression_sensitive.csv"
queriedMinorPath = "spreadsheets/normalized_sub_mean_expression_sensitive.csv"

plotPath = paste0("figures/", analysisName, '/')
dir.create(plotPath, showWarnings = FALSE)

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath)
queriedMinor <- fread(queriedMinorPath)

# Munge ----

# Columns: Name, Parent_Name, Marker, Marker_Function, Marker_Source, Marker_Notes

## Curated ----
### Replace apostrophe with ASCII
### Search for non-ASCII characters
### Search for blank-like cells -- none
### Reformat General Cell Markers into their own column and ensure that all subclasses start with subtype name or n/a
### Data -> from Table/Range
#### Remove blank rows in table (340 -> 333)
#### Fill down in name column
#### Replace n/a with R's NA
#### Replace empty with NA
#### Fixed rows with more than 1 gene (338)
#### Renamed generic subtypes
### Close and Load

# Move markers as left as possible
curated[is.na(`Type within Subclass`) & # The minor type is undefined
        !is.na(`Type-Specific Markers`) & # There is a value in the field
        (is.na(Markers) | # There is nothing in this field so it can move up
         Markers == `Type-Specific Markers`), ] # This field is a copy and can be re-written

curated[is.na(Subclass) & !is.na(Markers), # The minor type is undefined but There is a value in its field
        `Cell Marker` := Markers]
curated[is.na(Subclass) & !is.na(Markers), # The minor type is undefined but There is a value in its field
        Markers := NA]

# Make sure genes are in mouse format
gene_cols = c("Cell Marker", "Markers", "Type-Specific Markers")
curated[, cbind(gene_cols) := lapply(.SD, str_to_title), .SDcols = gene_cols ]

# Descend the tree, this can probably be looped
major_curated = curated[is.na(`Type-Specific Markers`) & is.na(Markers), ] %>%
  reframe(Name = `Cell Type`,
          Parent_Name = "Root",
          Level = 1,
          Marker = `Cell Marker`,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

middl_curated = curated[is.na(`Type-Specific Markers`) & is.na(`Cell Marker`), ] %>%
  reframe(Name = Subclass,
          Parent_Name = `Cell Type`,
          Level = 2,
          Marker = Markers,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

minor_curated = curated[is.na(Markers) & is.na(`Cell Marker`), ] %>%
  reframe(Name = `Type within Subclass`,
          Parent_Name = Subclass,
          Level = 3,
          Marker = `Type-Specific Markers`,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

# Are all rows accounted for? YES!
nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated) == nrow(curated) # TRUE

clean_curated = rbind(major_curated, middl_curated, minor_curated)
clean_curated$Set = "Clean"

# Next time, I think that only a few changes should be made to the original file
# to make it human understandable and machine understandable
# 1) columns: MajorClass, Class, Minor, Gene, Marker_Source, Marker_Annotation, Gene_Annotation (Optional), Notes
# 2) Be careful with copy and paste with regards to non-ASCII characters
# 3) No need for blank rows (Optional)
# 4) No headers with spaces, underscores as the only non-alphanumeric
# 5) NA instead of n/a
# extra) Columns should have 1 type of thing. It doesn't matter what it is, as long as it's the same.

## Queried ----

majorIndex = queriedMajor[, lapply(.SD, which.max), .SDcols = 2:ncol(queriedMajor)] %>% unlist()
minorIndex = queriedMinor[, lapply(.SD, which.max), .SDcols = 2:ncol(queriedMinor)] %>% unlist()

clean_queried =
  rbind(cbind(Name = queriedMajor[majorIndex, majorclass],
              Marker = colnames(queriedMajor)[-1] %>% str_to_title()),
        cbind(queriedMinor[minorIndex, author_cell_type],
              colnames(queriedMinor)[-1] %>% str_to_title())) %>%
  as.data.frame()
clean_queried$Set = "Query"

## Harmonize the datasets ----
## It would take longer to do this programmatically than manually, so here's my manual effort.
## I will use my query as the base and conform the curated set to it.
##
## At some free time, check whether there was a column with matching names
##
## 1) Remove plurals
## 2) Group names to title case
## 3) Written out to shorthand
##
##

unique(clean_curated$Name) %>% length() # 173 Groups # reduced to 87???
unique(clean_queried$Name) %>% length() # 109 Groups

# No additional string sanitation needed with regards to capitalization...
# intersect(unique(clean_curated$Name), unique(clean_queried$Name)) %>% length() ==
#   intersect(toupper(unique(clean_curated$Name)), toupper(unique(clean_queried$Name))) %>% length()

num_matching_groups <- function() {
  intersect(unique(clean_curated$Name), unique(clean_queried$Name)) %>% length()
}

num_matching_groups() # 16, 13 after clean_curated is changed

unique(clean_curated$Marker) %>% length() # 248 Markers
unique(clean_queried$Marker) %>% length() # 493 Markers

# Plurals to Single

endsWith(clean_curated$Name, 's') %>% sum() # 114
endsWith(clean_queried$Name, 's') %>% sum() # 0

clean_curated$Name = sub("s$", "", clean_curated$Name)
clean_curated$Parent_Name = sub("s$", "", clean_curated$Parent_Name)
num_matching_groups() # 20
endsWith(clean_curated$Name, 's') %>% sum() # 0
endsWith(clean_curated$Parent_Name, 's') %>% sum() # 0
unique(clean_curated$Name) %>% length() # 173 Groups

# Group names to title case

clean_queried$Name = toupper(clean_queried$Name)
clean_curated$Name = toupper(clean_curated$Name)
clean_curated$Parent_Name = toupper(clean_curated$Parent_Name)
num_matching_groups() # 20 # 19 after curated was cleaned
unique(clean_curated$Name) %>% length() # 173 Groups

# Long to Shorthand Curate, lots of manual here
# I want to leave out some query stuff because it is kinda cringe to have so many numbers
# AC and alphaon stuff and number tricks
source("CurationVsQuery.util.R")

length(CURATE2QUERY) # 6 to 7

# Test
cbind(from = clean_curated$Name[clean_curated$Name %in% names(CURATE2QUERY)],
      to   = CURATE2QUERY[match(clean_curated$Name, names(CURATE2QUERY)) %>% na.omit()])

# Do
clean_curated$Name[clean_curated$Name %in% names(CURATE2QUERY)] =
  CURATE2QUERY[match(clean_curated$Name, names(CURATE2QUERY)) %>% na.omit()]

# Test
cbind(from = clean_curated$Parent_Name[clean_curated$Parent_Name %in% names(CURATE2QUERY)],
      to   = CURATE2QUERY[match(clean_curated$Parent_Name, names(CURATE2QUERY)) %>% na.omit()])

# Do
clean_curated$Parent_Name[clean_curated$Parent_Name %in% names(CURATE2QUERY)] =
  CURATE2QUERY[match(clean_curated$Parent_Name, names(CURATE2QUERY)) %>% na.omit()]

num_matching_groups() # 26 = 20 + length(CURATE2QUERY) - Ganglion Cell

# AC Query

sub("AC_", "AC", clean_queried$Name) %>% unique() %>% length() # 109 groups, still fine
clean_queried$Name = sub("AC_", "AC", clean_queried$Name)
num_matching_groups() # 78 / 109

# Let's move numbers to after the name in Query

grep("[0-9]*_", clean_queried$Name %>% unique(), value = TRUE) %>% # Manually checked, all these fit the pattern
  { sub("([0-9]*)_(.*)", "\\2_\\1", .) }

clean_queried$Name = sub("([0-9]*)_(.*)", "\\2_\\1", clean_queried$Name)
num_matching_groups() # 85 / 109 # Not unexpected, this was mostly so it was easier for my human eyes to look at all of the classes

length(QUERY2CURATE) # 17

# Test
cbind(from = clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)],
      to   = QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()])

# Do
clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)] =
  QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()]
unique(clean_queried$Name) %>% length() # 108 Groups, we lost 1! W3D1 used to be two
num_matching_groups() # 101 / 108

clean_queried$Name %>% unique() %>% setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique())
# "RGC"        "NOVEL_10"   "OODS_NT_12" "OODS_DV_16" "NOVEL_30"   "M1_33"      "M1DUP_40"
# 7 Classes

# Add Original Names?

# Let's add a layer to the tree
# 1) Add the new node with any of its markers or even none to curated
# 2) Point its children at it in curated
# 3) Update the name in Query
# 4) Done

# RGC
## children: grep("RGC")

# OODS_NT_12
## OODS_N OODST

# OODS_DV_16
## OODSV OODSD

# M1_33 M1DUP40
## Combine into just M1

# "GC"->"Other"-> Novel_10 & Novel_30

# Add Original Names?

# Stop and Review ----

# Explore and Compare ----
# clean_curated, clean_queried

# Query had fewer groups, but more markers...

Marker_Intersection = merge(clean_curated, clean_queried, by = c("Name", "Marker"), all = F)
dplyr::inner_join(clean_curated, clean_queried, by = c("Name", "Marker"))
# Scratch ----

if (FALSE) {
  # Are all rows accounted for?
  nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated) == nrow(curated) # FALSE
  nrow(curated) - (nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated)) # 12

  # OK... where are the rest of the rows?

  curated[, sapply(.SD, is.na) %>% rowSums() %>% {. == 2} %>% all(), .SDcols = gene_cols] # FALSE
  missing = curated[, sapply(.SD, is.na) %>% rowSums() %>% {. == 2}, .SDcols = gene_cols]

  curated[!missing, .SD, .SDcols = gene_cols]
}
