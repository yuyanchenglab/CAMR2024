#!/bin/env Rscript3

# JM 08/12/2024
# How do my statistics methods compare to a literature search when looking for proper cell markers?

# Setup ----

library(dplyr) # pipe
library(data.table) # fast table manip
library(ggplot2)
library(stringr)
library(ggvenn) # venn
library(ggVennDiagram) # upset plot
library(eulerr)

cheng_theme = new.env() # Note: Use dollar sign operator to access values
source("/project/hipaa_ycheng11lab/software/apps/github/figures/theme.cheng.module.R", local = cheng_theme)

analysisName = "CurationVsQueryByCell"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)

curatedPath = "spreadsheets/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "spreadsheets/ovr_top_filtered_genes_majorclass_coefficients_sensitive.csv"

plotPath = paste0("figures/", analysisName, '/')
dir.create(plotPath, showWarnings = FALSE)

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

majorclass = c("AC", "ASTROCYTE", "BC", "CONE", "ENDOTHELIAL", "HC", "MG", "MICROGLIA", "PERICYTE", "RGC", "ROD", "RPE")

needed_major = c("AC", "BC", "Microglia", "RGC")

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath,
                      select = c("Gene", "Cell Type"),
                      col.names = c("Marker", "Name"))
queriedMinor <- sapply(needed_major,
                       \(major) { fread(paste0("spreadsheets/ovr_top_", "complete", "_filtered_markers_", major, "_coefficients_sensitive.csv"),
                                               select = c("Gene", "Cell", "Subclass")) },
                       USE.NAMES = TRUE, simplify = FALSE) %>%
  rbindlist()

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
### Close and Load Here

# Move markers as left as possible
curated[is.na(`Type within Subclass`) & !is.na(`Type-Specific Markers`), ] # None

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
          Parent_Name = `Cell Type`,
          Marker = `Cell Marker`,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

middl_curated = curated[is.na(`Type-Specific Markers`) & is.na(`Cell Marker`), ] %>%
  reframe(Name = Subclass,
          Parent_Name = `Cell Type`,
          Marker = Markers,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

minor_curated = curated[is.na(Markers) & is.na(`Cell Marker`), ] %>%
  reframe(Name = `Type within Subclass`,
          Parent_Name = `Cell Type`,
          Marker = `Type-Specific Markers`,
          Marker_Function = `Function of Marker`,
          Marker_Source = `Image Link`,
          Marker_Notes = `Notes`)

# Are all rows accounted for? YES!
nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated) == nrow(curated) # TRUE

clean_curated = rbind(major_curated, middl_curated, minor_curated) %>%
  mutate(Set = "Curated",
         Original_Marker = Marker,
         Original_Name = Name,
         Original_Parent_Name = Parent_Name)
setDT(clean_curated)

## Queried ----

queriedMinor = queriedMinor %>%
  reframe(Name = Subclass,
          Marker = str_to_title(Gene),
          Parent_Name = Cell)

queriedMajor$Parent_Name = queriedMajor$Name

clean_queried =
  rbindlist(list(Major = queriedMajor, Minor = queriedMinor),
            idcol = "MajorMinor", fill = TRUE) %>%
  summarise(.by = c(Name, Marker), Parent_Name = unique(Parent_Name),
            MajorMinor = if (n() > 1) "Major&Minor" else MajorMinor)  %>%
  mutate(Set = "Query",
         Original_Marker = Marker,
         Original_Name = Name,
         Original_Parent_Name = Parent_Name)

## Harmonize the datasets ----
## I will try to use my query as the base and conform the curated set to it.

unique(clean_curated$Name) %>% length() # 173 Groups
unique(clean_queried$Name) %>% length() # 135 Groups

num_matching_groups <- function() {
  intersect(unique(clean_curated$Name), unique(clean_queried$Name)) %>% length()
}
num_matching_groups() # 16

## 1) Markers to Mouse Gene formatting

### Get rid of spaces in markers, added later 5:14
clean_curated$Marker = gsub(" ", "", clean_curated$Marker)
unique(clean_curated$Marker) %>% length() # 248->247 Markers

clean_queried$Marker = str_to_title(clean_queried$Marker)
clean_curated$Marker = str_to_title(clean_curated$Marker)
num_matching_groups() # 16
unique(clean_curated$Marker) %>% length() # 248 Markers
unique(clean_queried$Marker) %>% length() # 624 Markers

clean_curated$Marker[clean_curated$Marker == "Bhlhe22/Bhlhb5"] = "Bhlhe22"
clean_curated$Marker[clean_curated$Marker == "Prkca/Pkca"] = "Prkca"
clean_curated$Marker[clean_curated$Marker == "Lect1"] = "Cnmd"
clean_curated$Marker[clean_curated$Marker == "Seripini1"] = "Serpini1"
clean_curated$Marker[clean_curated$Marker == "Lgsf3"] = "Igsf3"
clean_curated$Marker[clean_curated$Marker == "Lgfbp5"] = "Igfbp5"
clean_curated$Marker[clean_curated$Marker == "Anza3"] = "Anxa3" # Fig2D of Tran et al
clean_curated$Marker[clean_curated$Marker == "A230065h16rik"] = "Lbhd2" # https://www.phosphosite.org/proteinAction?id=5141911
clean_curated$Marker[clean_curated$Marker == "Fam19a4"] = "Tafa4" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
clean_curated$Marker[clean_curated$Marker == "4833423e24rik"] = "Fads2b" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
unique(clean_curated$Marker) %>% length() # 246 Markers

## 2) Remove plurals, curate

endsWith(clean_curated$Name, 's') %>% sum() # 114
endsWith(clean_queried$Name, 's') %>% sum() # 0
clean_queried$Name = sub("s$", "", clean_queried$Name)
clean_curated$Name = sub("s$", "", clean_curated$Name)
clean_curated$Parent_Name = sub("s$", "", clean_curated$Parent_Name)
clean_queried$Parent_Name = sub("s$", "", clean_queried$Parent_Name)
num_matching_groups() # 20
endsWith(clean_curated$Name, 's') %>% sum() # 0
endsWith(clean_curated$Parent_Name, 's') %>% sum() # 0
unique(clean_curated$Name) %>% length() # 173 Groups, Clobber check
unique(clean_queried$Name) %>% length() # 173 Groups

# 3) Group names to title case

clean_queried$Name = toupper(clean_queried$Name)
clean_curated$Name = toupper(clean_curated$Name)
clean_curated$Parent_Name = toupper(clean_curated$Parent_Name)
clean_queried$Parent_Name = toupper(clean_queried$Parent_Name)
num_matching_groups() # 20
unique(clean_curated$Name) %>% length() # 173 Groups, Clobber check
unique(clean_queried$Name) %>% length() # 135 Groups

# 4) Manual here for curate
source("scripts/CurationVsQuery.util.R")

length(CURATE2QUERY) # 6 to 7

# Test
cbind(from = clean_curated$Name[clean_curated$Name %in% names(CURATE2QUERY)],
      to   = CURATE2QUERY[match(clean_curated$Name, names(CURATE2QUERY)) %>% na.omit()]) %>% data.frame() %>% distinct()
# Do
clean_curated$Name[clean_curated$Name %in% names(CURATE2QUERY)] =
  CURATE2QUERY[match(clean_curated$Name, names(CURATE2QUERY)) %>% na.omit()]

# Test
cbind(from = clean_curated$Parent_Name[clean_curated$Parent_Name %in% names(CURATE2QUERY)],
      to   = CURATE2QUERY[match(clean_curated$Parent_Name, names(CURATE2QUERY)) %>% na.omit()]) %>% data.frame() %>% distinct()
# Do
clean_curated$Parent_Name[clean_curated$Parent_Name %in% names(CURATE2QUERY)] =
  CURATE2QUERY[match(clean_curated$Parent_Name, names(CURATE2QUERY)) %>% na.omit()]
num_matching_groups() # 27 = 20 + length(CURATE2QUERY)

# 5) AC for Query

sub("AC_", "AC", clean_queried$Name) %>% unique() %>% length() # 135 groups, still fine
clean_queried$Name = sub("AC_", "AC", clean_queried$Name)
num_matching_groups() # 89 / 135 # 87 / 135

# 6) numbers to after name in Query
grep("[0-9]*_", clean_queried$Name %>% unique(), value = TRUE) %>% # Manually checked, all these fit the pattern
  { sub("([0-9]*)_(.*)", "\\2_\\1", .) }
clean_queried$Name = sub("([0-9]*)_(.*)", "\\2_\\1", clean_queried$Name)
num_matching_groups() # 106 / 135 # Not unexpected, this was mostly so it was easier for my human eyes to look at all of the classes

# 7) Manual for Query
length(QUERY2CURATE) # 25

# Test
cbind(from = clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)],
      to   = QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()]) %>% as.data.frame() %>% distinct()
# Do
clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)] =
  QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()]
unique(clean_queried$Name) %>% length() # 134 Groups, we lost 1! W3D1 used to be two
num_matching_groups() # 127 / 135

clean_queried$Name %>% unique() %>% setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique())
# "NOVEL_10"    "NT-OODSGC"   "NOVEL_13"    "OODS_CCK_14" "DV-OODSGC"   "NOVEL_24"    "NOVEL_30"

# Remove any odd duplicates
nrow(clean_queried) == nrow(clean_queried %>% summarize(.by = c(Name, Marker)) %>% distinct()) # 2321
nrow(clean_curated) == nrow(clean_curated %>% summarize(.by = c(Name, Marker)) %>% distinct()) # 336, mismatch
nrow(clean_curated) == nrow(clean_curated %>% distinct()) # 337
clean_curated = clean_curated %>%
  filter(!duplicated(paste0(Name, Marker))) #] = "JM 08/06/2024 Duplicate removed."

# Stop and Review ----

fwrite(clean_curated, "spreadsheets/curated_table.txt")
fwrite(clean_queried, "spreadsheets/queried_by_cell_table.txt")

rm(curated, major_curated, middl_curated, minor_curated, queriedMajor, queriedMinor,
   QUERY2CURATE, CURATE2QUERY); gc()

# Explore and Compare ----
# clean_curated, clean_queried

clean_curated$Name %>% unique %>% length
clean_queried$Name %>% unique %>% length

clean_curated$Marker %>% unique %>% length # 248
clean_queried$Marker %>% unique %>% length # 624, no clobbering of Markers

gene_list = list(Curated = unique(clean_curated$Marker), Queried = unique(clean_queried$Marker))
cell_list = list(Curated = unique(clean_curated$Name), Queried = unique(clean_queried$Name))

euler_gene = euler(gene_list)
vd = plot(euler_gene, quantities = list(fontsize = 20), main = list(fontsize = 20, label = "Gene List Overlap", vjust = 1), labels = list(fontsize = 20))
ggsave(paste0(plotPath, "VD.Genes.png"), vd)

upGene = ggVennDiagram(gene_list, force_upset = TRUE)
ggsave(paste0(plotPath, "Upset.Genes.png"), upGene, width = 6, height = 3)

euler_cell = euler(cell_list)
vd = plot(euler_cell, quantities = list(fontsize = 20), main = list(fontsize = 20, label = "Cell List Overlap", vjust = 1), labels = list(fontsize = 20))
ggsave(paste0(plotPath, "VD.Cells.png"), vd)

upName = ggVennDiagram(cell_list, force_upset = TRUE)
ggsave(paste0(plotPath, "Upset.Cells.png"), upName, width = 6, height = 3)

## Merge ----
# Take a Union of the tables to get a table with cellxgene and
# bool column fromQuery, fromcurated
# JM 8/12/2024

clean_merged =
  merge(clean_curated, clean_queried, by = c("Marker", "Name"), all = TRUE) %>%
  summarize(.by = c(Name, Marker),
            Parent_Name = ifelse(is.na(Parent_Name.x), Parent_Name.y, Parent_Name.x),
            # Marker_Methods_Unofficial =
            #   if (!is.na(Set.x) & !is.na(Set.y)) { paste0(Set.x, "&", MajorMinor) }
            #   else if (is.na(MajorMinor)) { Set.x }
            #   else { MajorMinor },
            # Marker_Methods_Official = sub("Major&Minor|Major|Minor", "Queried", Marker_Methods_Unofficial),
            Curated = ifelse(!is.na(Set.x), "Curated", NA),
            Queried = ifelse(!is.na(Set.y), "Queried", NA),
            Parent_Name,
            Queried_Name = Original_Name.x,
            Queried_Parent_Name = Original_Parent_Name.x)

all(unique(clean_merged$Parent_Name) %in% majorclass) # FALSE
all(majorclass %in% unique(clean_merged$Parent_Name)) # TRUE

m2a = fread("spreadsheets/major_author.csv", header = TRUE) # query
m2a_curated = clean_curated %>% summarise(.by = c(Name,Parent_Name)) # curated
library(magrittr)

clean_merged$Major_Name = NA

parent_is_majorclass = clean_merged$Parent_Name %in% majorclass
clean_merged$Major_Name[parent_is_majorclass] = clean_merged$Parent_Name[parent_is_majorclass] %>% toupper()
sum(is.na(clean_merged$Major_Name)) # 89

is_majorclass = is.na(clean_merged$Major_Name) &
  clean_merged$Name %in% majorclass
clean_merged$Major_Name[is_majorclass] = clean_merged$Name[is_majorclass] %>% toupper()
sum(is.na(clean_merged$Major_Name)) # 58

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("EPITHELIUM", clean_merged$Name)] = "EPITHELIAL"
sum(is.na(clean_merged$Major_Name)) # 54

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("DENDRITIC CELL|MONOCYTE|T CELL|B CELL|NK CELL|MACROPHAGE|BAM|FIBROBLAST|NEUTROPHIL", clean_merged$Name)] = "IMMUNE"
sum(is.na(clean_merged$Major_Name)) # 25

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("ASTROCYTE", clean_merged$Name)] = "ASTROCYTE" #= "MELANOCYTE"
sum(is.na(clean_merged$Major_Name)) # 12

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl(" MG", clean_merged$Name)] = "MG"
sum(is.na(clean_merged$Major_Name)) # 5

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("MICROGLIA", clean_merged$Name)] = "MG"
sum(is.na(clean_merged$Major_Name)) # 2

clean_merged$Major_Name[is.na(clean_merged$Major_Name)] = "MONOCYTE"
sum(is.na(clean_merged$Major_Name)) # 0!

clean_merged = clean_merged %>%
  arrange(Major_Name, Name, Marker) %>%
  reframe(Name, Marker, Curated, Queried, Major_Name) %T>%
  fwrite("spreadsheets/merged_curated-queried_ByCell_markers.txt", sep = '\t')

# break full table into files by major class
for (cell_type in unique(clean_merged$Major_Name)) {
  clean_merged %>%
    filter(Major_Name == cell_type) %>%
    reframe(Name, Marker, Curated, Queried) %>%
    fwrite(paste0("spreadsheets/merged_curated-queried_", cell_type, "_markers.txt"), sep = '\t')
}

# How did the search space overlap with the curated set to begin with?

variable_genes = fread("data/CAMR_all_genes.csv")
clean_curated$Marker %>% unique() %>% length() # 247 genes
curated_markers = unique(clean_curated$Marker)
variable_genes_formatted = str_to_title(variable_genes$feature_name)
curated_markers[curated_markers %in% variable_genes_formatted] # 231
curated_markers[!curated_markers %in% variable_genes_formatted] # 16 # This one
data.frame(Curated = curated_markers,
           Highly_Variable = curated_markers %in% variable_genes_formatted) %>%
  fwrite("spreadsheets/genes_in_curated_markers.txt", sep = '\t')

sum(major_curated$Marker %>% unique() %>% str_to_title() %in% variable_genes_formatted)
sum(middl_curated$Marker %>% unique() %>% str_to_title() %in% variable_genes_formatted)

# Scratch ----

if (FALSE) {

}
