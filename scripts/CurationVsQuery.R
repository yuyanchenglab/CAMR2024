#!/bin/env Rscript3

# JM 7/31/2024
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

analysisName = "CurationVsQuery"
analysisPath = "/project/hipaa_ycheng11lab/atlas/CAMR2024/"
setwd(analysisPath)

curatedPath = "spreadsheets/CuratedMouseRetinaMarkers.txt"
queriedMajorPath = "spreadsheets/ovr_top_filtered_genes_majorclass_coefficients_sensitive.csv"
queriedMinorPath = "spreadsheets/ovr_top_filtered_genes_author_cell_type_coefficients_sensitive.csv"

plotPath = paste0("figures/", analysisName, '/')
dir.create(plotPath, showWarnings = FALSE)

## Reproducible Results ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

curated <- fread(curatedPath)
queriedMajor <- fread(queriedMajorPath,
                      select = c("Gene", "Cell Type"),
                      col.names = c("Marker", "Name"))
queriedMinor <- fread(queriedMinorPath,
                      select = c("Gene", "Cell Type"),
                      col.names = c("Marker", "Name"))

## Check queried genes ----

# queried_source_major = fread("spreadsheets/ovr_top_20_genes_by_cell_type_reproduction.csv")
# queried_original_major = fread("spreadsheets/ovr_top_20_genes_by_cell_type.csv")
# queried_source_minor = fread("spreadsheets/ovr_top_20_genes_by_sub_cell_type_reproduction.csv")
#
# all(colnames(queriedMajor)[-1] %in% unique(queried_source_major$Gene))
# unique(queried_source_major$Gene)[!colnames(queriedMajor)[-1] %in% unique(queried_source_major$Gene)]
#
# all(colnames(queriedMajor)[-1] %in% unique(queried_original_major$Gene))
# unique(queried_original_major$Gene)[!colnames(queriedMajor)[-1] %in% unique(queried_original_major$Gene)]
#
# all(colnames(queriedMinor)[-1] %in% unique(queried_source_minor$Gene))
#
# ### Pos only?
#
# queried_source_major_pos = queried_source_major[Coefficient > 0, ]
# queried_source_minor_pos =queried_source_minor[Coefficient > 0, ]
#
# qsjp = unique(queried_source_major_pos$Gene)
# all(colnames(queriedMajor)[-1] %in% qsjp)
# colnames(queriedMajor)[-1][!colnames(queriedMajor)[-1] %in% qsjp]
#
# # Just Cd74
#
# all(colnames(queriedMinor)[-1] %in% unique(queried_source_minor_pos$Gene))


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
curated[is.na(`Type within Subclass`) & # The minor type is undefined
        !is.na(`Type-Specific Markers`),# & # There is a value in the field
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
          Parent_Name = Subclass,
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

# Next time, I think that only a few changes should be made to the original file
# to make it human understandable and machine understandable
# 1) columns: MajorClass, Class, Minor, Gene, Marker_Source, Marker_Annotation, Gene_Annotation (Optional), Notes
# 2) Be careful with copy and paste with regards to non-ASCII characters
# 3) No need for blank rows (Optional)
# 4) No headers with spaces, underscores as the only non-alphanumeric
# 5) NA instead of n/a
# extra) Columns should have 1 type of thing. It doesn't matter what it is, as long as it's the same.

## Queried ----

clean_queried =
  rbindlist(list(Major = queriedMajor, Minor = queriedMinor), idcol = "MajorMinor") %>%
  summarise(.by = c(Name, Marker),
            MajorMinor = if (n() > 1) "Major&Minor" else MajorMinor)  %>%
  mutate(Set = "Query",
         Original_Marker = Marker,
         Original_Name = Name) #%>% {table(.$MajorMinor)}

## Harmonize the datasets ----
## It would take longer to do this programmatically than manually, so here's my manual effort.
## I will use my query as the base and conform the curated set to it.
##
## At some free time, check whether there was a column with matching names
##
## 1) Remove plurals, curate
## 2) Group names to uppercase, curate
## 3) Individual manual adjustment -- Written out to shorthand, curate
## 4) Remove AC underscores, query
## 5) Put numbers at end, query
## 6) Add intermediates to tree, curate
##

unique(clean_curated$Name) %>% length() # 173 Groups # reduced to 87???
unique(clean_queried$Name) %>% length() # 135 Groups

# No additional string sanitation needed with regards to capitalization...
# intersect(unique(clean_curated$Name), unique(clean_queried$Name)) %>% length() ==
#   intersect(toupper(unique(clean_curated$Name)), toupper(unique(clean_queried$Name))) %>% length()

num_matching_groups <- function() {
  intersect(unique(clean_curated$Name), unique(clean_queried$Name)) %>% length()
}

num_matching_groups() # 16

clean_queried$Name = str_to_title(clean_queried$Name)
clean_curated$Name = str_to_title(clean_curated$Name)
clean_curated$Parent_Name = str_to_title(clean_curated$Parent_Name)

clean_queried$Marker = str_to_title(clean_queried$Marker)
clean_curated$Marker = str_to_title(clean_curated$Marker)

num_matching_groups() # 16

unique(clean_curated$Marker) %>% length() # 248 Markers
unique(clean_queried$Marker) %>% length() # 493 Markers

# Plurals to Single

endsWith(clean_curated$Name, 's') %>% sum() # 114
endsWith(clean_queried$Name, 's') %>% sum() # 38

clean_queried$Name = sub("s$", "", clean_queried$Name)
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
num_matching_groups() # 20
# Clobber check
unique(clean_curated$Name) %>% length() # 173 Groups
unique(clean_queried$Name) %>% length() # 135 Groups

# Long to Shorthand Curate, lots of manual here
# I want to leave out some query stuff because it is kinda cringe to have so many numbers
# AC and alphaon stuff and number tricks
source("CurationVsQuery.util.R")

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

num_matching_groups() # 26 = 20 + length(CURATE2QUERY) - Ganglion Cell

# AC Query

sub("AC_", "AC", clean_queried$Name) %>% unique() %>% length() # 135 groups, still fine
clean_queried$Name = sub("AC_", "AC", clean_queried$Name)
num_matching_groups() # 87 / 135

# Let's move numbers to after the name in Query

grep("[0-9]*_", clean_queried$Name %>% unique(), value = TRUE) %>% # Manually checked, all these fit the pattern
  { sub("([0-9]*)_(.*)", "\\2_\\1", .) }

clean_queried$Name = sub("([0-9]*)_(.*)", "\\2_\\1", clean_queried$Name)
num_matching_groups() # 104 / 135 # Not unexpected, this was mostly so it was easier for my human eyes to look at all of the classes

length(QUERY2CURATE) # 21

# Test
cbind(from = clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)],
      to   = QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()]) %>% as.data.frame() %>% distinct()

# Do
clean_queried$Name[clean_queried$Name %in% names(QUERY2CURATE)] =
  QUERY2CURATE[match(clean_queried$Name, names(QUERY2CURATE)) %>% na.omit()]
unique(clean_queried$Name) %>% length() # 134 Groups, we lost 1! W3D1 used to be two
num_matching_groups() # 126 / 135

clean_queried$Name %>% unique() %>% setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique())
# "RGC"        "NOVEL_10"   "OODS_NT_12" "OODS_DV_16" "NOVEL_30"  "NOVEL_13" ""NOVEL_30" "DV-OODSGC"| also before taking plunge:  "M1_33"      "M1DUP_40"
# 8 Classes

clean_curated$Name[clean_curated$Name == "GC"] = "RGC"
clean_curated$Parent_Name[clean_curated$Parent_Name == "GC"] = "RGC"

clean_queried$Name %>% unique() %>% setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique())
# 7 classes now not include RGC

# Let's add a layer to the tree
# 1) Add the new node with any of its markers or even none to curated
# 2) Point its children at it in curated
# 3) Update the name in Query
# 4) Done

if (FALSE) {
  new_cells = c("RGC", "NT-OODSGC", "DV-OODSGC")
  new_parents = c("GC","OODSGC","OODSGC")
  clean_curated =
    lapply(1:length(new_cells), \(i) {
      clean_curated[clean_curated$Name == new_parents[i], ] %>%
        # Inherit markers from above, not below though
        mutate(Parent_Name = Name,
               Name = new_cells[i],
               Marker_Source = "JM 08/02/2024 Query", # Override source
               Marker_Notes = "JM 08/02/2024 Added Manually",
               Original_Marker = NA) # Override notes
    }) %>% t() %>%
    rbindlist() %>%
    rbind(clean_curated)

  # RGC: except J-RGC that are a subsubtype
  clean_curated$Parent_Name[grep("^[^J]*RGC", clean_curated$Name)] = "RGC" # JM Ran @ 3:17

  # NT-OODSGC: OODS_N but not OODST # They don't express the OODSGC marker
  clean_curated$Parent_Name[clean_curated$Name == "N-OODSGC"] = "NT-OODSGC"
  clean_curated$Parent_Name[clean_curated$Name == "MAYBE T-OODSGC"] = "GC"

  # DV-OODSGC: OODSV OODSD
clean_curated$Parent_Name[clean_curated$Name %in% c("V-OODSGC", "D-OODSGC")] = "DV-OODSGC"
}
# Stop and Review ----

# 3:35 -- It all looks good to me now.

# Whoopsies, there are duplicates in query...

nrow(clean_queried) == nrow(clean_queried %>% summarize(.by = c(Name, Marker)) %>% distinct()) # 2321
nrow(clean_curated) == nrow(clean_curated %>% summarize(.by = c(Name, Marker)) %>% distinct()) # 336, mismatch
nrow(clean_curated) == nrow(clean_curated %>% distinct()) # 337
clean_curated = clean_curated %>%
  filter(!duplicated(paste0(Name, Marker))) #] = "JM 08/06/2024 Duplicate removed."

# Let's also add some missing edges

fwrite(clean_curated, "spreadsheets/curated_table.txt")
fwrite(clean_queried, "spreadsheets/queried_table.txt")

rm(curated, major_curated, middl_curated, minor_curated, queried_original_major,
   queried_source_major, queried_source_major_pos, queried_source_minor,
   queried_source_minor_pos, queriedMajor, queriedMinor,  majorIndex, minorIndex,
   new_cells, new_parents, qsjp, QUERY2CURATE, CURATE2QUERY); gc()

# Explore and Compare ----
# clean_curated, clean_queried

clean_curated$Name %>% unique %>% length
clean_queried$Name %>% unique %>% length

clean_curated$Marker %>% unique %>% length
clean_queried$Marker %>% unique %>% length

gene_list = list(Curated = unique(clean_curated$Marker), Queried = unique(clean_queried$Marker))
cell_list = list(Curated = unique(clean_curated$Name), Queried = unique(clean_queried$Name))

# ggvenn(gene_list, auto_scale = T, text_size = 5)
# ggsave(paste0(plotPath, "VD.Genes.png"))

euler_gene = euler(gene_list)
vd = plot(euler_gene, quantities = list(fontsize = 20), main = list(fontsize = 20, label = "Gene List Overlap", vjust = 1), labels = list(fontsize = 20))
ggsave(paste0(plotPath, "VD.Genes.png"), vd)

upGene = ggVennDiagram(gene_list, force_upset = TRUE)
ggsave(paste0(plotPath, "Upset.Genes.png"), upGene, width = 6, height = 3)

# ggvenn(cell_list, auto_scale = T, text_size = 5)
# ggsave(paste0(plotPath, "VD.Cells.png"))

upName = ggVennDiagram(cell_list, force_upset = TRUE)
ggsave(paste0(plotPath, "Upset.Cells.png"), upName, width = 6, height = 3)

euler_cell = euler(cell_list)
vd = plot(euler_cell, quantities = list(fontsize = 20), main = list(fontsize = 20, label = "Cell List Overlap", vjust = 1), labels = list(fontsize = 20))
ggsave(paste0(plotPath, "VD.Cells.png"), vd)

## Merge ----
# Try a Union approach? Take a Union of the tables to get a table with cellxgene and
# bool column fromQuery, fromcurated
# JM 8/6/2024

clean_merged =
  merge(clean_curated, clean_queried, by = c("Marker", "Name"), all = TRUE) %>%
  summarize(.by = c(Name, Marker),
            Parent_Name,
            Marker_Methods_Unofficial =
              if (!is.na(Set.x) & !is.na(Set.y)) { paste0(Set.x, "&", MajorMinor) }
              else if (is.na(MajorMinor)) { Set.x }
              else { MajorMinor },
            Marker_Methods_Official = sub("Major&Minor|Major|Minor", "Queried", Marker_Methods_Unofficial),
            Queried_Name = Original_Name.y,
            Marker_Notes) %>%
  mutate(Curated = c(NA, "Curated")[1+grepl("Curated", Marker_Methods_Official)],
         Queried = c(NA, "Queried")[1+grepl("Queried", Marker_Methods_Official)]) %>%
  reframe(Name, Marker, Curated, Queried, Parent_Name, Queried_Name)

majorclass = c("AC", "ASTROCYTE", "BC", "CONE", "ENDOTHELIAL", "HC", "MG", "MICROGLIA", "PERICYTE", "RGC", "ROD", "RPE")
# all(majorclass %in% c(clean_merged$GrandParent_Name, clean_merged$Parent_Name, clean_merged$Name)) # FALSE

m2a = fread("spreadsheets/major_author.csv", header = TRUE) # query

m2a_curated = clean_curated %>% summarise(.by = c(Name,Parent_Name)) # curated

clean_merged$Major_Name = NA
is_subtype = clean_merged$Queried_Name %in% m2a$author_cell_type
clean_merged$Major_Name[is_subtype] = m2a$majorclass[match(clean_merged$Queried_Name[is_subtype], m2a$author_cell_type)] %>% toupper()

parent_is_majorclass = is.na(clean_merged$Major_Name) & toupper(clean_merged$Parent_Name) %in% majorclass
clean_merged$Major_Name[parent_is_majorclass] = clean_merged$Parent_Name[parent_is_majorclass] %>% toupper()

name_is_majorclass = is.na(clean_merged$Major_Name) & toupper(clean_merged$Name) %in% majorclass
clean_merged$Major_Name[name_is_majorclass] = clean_merged$Name[name_is_majorclass] %>% toupper()

parent_2_parent = is.na(clean_merged$Major_Name)
any(m2a_curated$Parent_Name[match(clean_merged$Name[parent_2_parent], m2a_curated$Name)] %in% majorclass) # FALSE
P1 = m2a_curated$Parent_Name[match(clean_merged$Name[parent_2_parent], m2a_curated$Name)]
P2 = m2a_curated$Parent_Name[match(P1, m2a_curated$Name)]
clean_merged$Major_Name[parent_2_parent] = P2

clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("NOVEL", clean_merged$Name)] = "RGC"
clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("^BC", clean_merged$Name)] = "BC"
clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("EPITHELIUM", clean_merged$Name)] = "EPITHELIAL"
clean_merged$Major_Name[is.na(clean_merged$Major_Name) & grepl("DENDRITIC CELL|MONOCYTE|T CELL|B CELL|NK CELL|MACROPHAGE|BAM|FIBROBLAST|NEUTROPHIL", clean_merged$Name)] = "IMMUNE"
clean_merged$Major_Name[is.na(clean_merged$Major_Name)] = "MELANOCYTE"

# clean_merged$Major_Name %>% is.na() %>% sum() # TRUE
# all(clean_merged$Parent_Name.x == clean_merged$Parent_Name.y, na.rm = TRUE)

# clean_merged$Parent_Name = clean_curated$Parent_Name[match(clean_merged$Name, clean_curated$Name)]
# clean_merged$GrandParent_Name = clean_curated$Parent_Name[match(clean_merged$Parent_Name, clean_curated$Name)]

clean_merged %>%
  reframe(Name, Marker, Curated, Queried, Major_Name) %>%
  fwrite("spreadsheets/merged_curated-queried_cell_markers.txt", sep = '\t')

# clean_merged = clean_merged %>%
#   rowwise() %>%
#   mutate(Major_Name = if (any(majorclass %in% Parent_Name)) {Parent_Name}
#          else if (any(majorclass %in% Name)) {Name}
#          else if (any(majorclass %in% GrandParent_Name)) {GrandParent_Name}
#          else {NA})


# break into major class
for (cell_type in unique(clean_merged$Major_Name)) {
  clean_merged %>%
    filter(Major_Name == cell_type) %>%
    reframe(Name, Marker, Curated, Queried) %>%
    fwrite(paste0("spreadsheets/merged_curated-queried_", cell_type, "_markers.txt"), sep = '\t')
}

# Scratch ----

if (FALSE) {
  # Are all rows accounted for?
  nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated) == nrow(curated) # FALSE
  nrow(curated) - (nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated)) # 12

  # OK... where are the rest of the rows?

  curated[, sapply(.SD, is.na) %>% rowSums() %>% {. == 2} %>% all(), .SDcols = gene_cols] # FALSE
  missing = curated[, sapply(.SD, is.na) %>% rowSums() %>% {. == 2}, .SDcols = gene_cols]

  curated[!missing, .SD, .SDcols = gene_cols]

  # 08/02/2024

  # "GC"->"Other"-> Novel_10 & Novel_30 from query to match curated

  clean_curated$Parent_Name[clean_curated$Name == "NOVEL_30"]  # They don't express the OODSGC marker
  # library(VennDiagram)
  #
  # venn.diagram(gene_list, filename = paste0(plotPath, "VD.Genes.png"), main = "Gene List Overlap", force.unique = TRUE)

  ## Inheritance ----

  # Try to cascade the gene

  # DFS
  while (nrow(sentinel)) {
    clean_curated = rbindlist(list(clean_curated, clean_curated$Marker[clean_curated$Parent_Name[sentinel]]))
  }

}
