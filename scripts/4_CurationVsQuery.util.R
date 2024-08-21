# Munge Curate ----

## Excel ----

## Pre-upload Excel Work
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

## R ----

munge_curate <- function(curated, verbose = FALSE) {
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
  if (verbose) {
    print("Do the number of rows match?")
    print(nrow(minor_curated) + nrow(major_curated) + nrow(middl_curated) == nrow(curated)) # TRUE
  }

  clean_curated = rbind(major_curated, middl_curated, minor_curated) %>%
    mutate(Curated = "Curated",
           Curated_Marker = Marker,
           Curated_Name = Name,
           Curated_Parent_Name = Parent_Name)
  setDT(clean_curated)

  return(clean_curated)

}

# Gene Naming ----

markers_to_mouse <- function(clean, verbose = FALSE) {

  start = unique(clean$Marker) %>% length()
  clean$Marker = gsub(" ", "", clean$Marker)
  clean$Marker = str_to_title(clean$Marker)

  end = unique(clean$Marker) %>% length()
  if (verbose) {
    print(paste("Starting Markers:", start))
    print(paste("Ending Markers:  ", end))
  }

  return(clean)
}

clean_genes <- function(clean, verbose = FALSE) {
  start_genes = unique(clean$Marker) %>% length()
  clean$Marker[clean$Marker == "Bhlhe22/Bhlhb5"] = "Bhlhe22"
  clean$Marker[clean$Marker == "Prkca/Pkca"] = "Prkca"
  clean$Marker[clean$Marker == "Lect1"] = "Cnmd"
  clean$Marker[clean$Marker == "Seripini1"] = "Serpini1"
  clean$Marker[clean$Marker == "Lgsf3"] = "Igsf3"
  clean$Marker[clean$Marker == "Lgfbp5"] = "Igfbp5"
  clean$Marker[clean$Marker == "Anza3"] = "Anxa3" # Fig2D of Tran et al
  clean$Marker[clean$Marker == "A230065h16rik"] = "Lbhd2" # https://www.phosphosite.org/proteinAction?id=5141911
  clean$Marker[clean$Marker == "Fam19a4"] = "Tafa4" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
  clean$Marker[clean$Marker == "4833423e24rik"] = "Fads2b" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
  end_genes = unique(clean$Marker) %>% length()
  if (verbose) {
    print(paste("Starting Genes:", start_genes))
    print(paste("Ending Genes:  ", end_genes))
  }

  return(clean)
}

# Cell Naming ----

cell_clobber_check <- function() {
  curated_cells = unique(clean_curated$Name)
  queried_cells = unique(clean_queried$Name)
  overlap_cells = intersect(curated_cells, queried_cells)
  print(paste("Curated Cells:", curated_cells %>% length()))
  print(paste("Queried Cells:", queried_cells %>% length()))
  print(paste("Overlap Cells:", overlap_cells %>% length()))
}

make_CURATE2QUERY <- function() {
  # View(clean_queried$Name %>% unique() %>% as.data.frame())
  clean_curated$Name %>% unique() %>% dput()
}

depluralize <- function(clean, verbose = FALSE) {
  name_plural = endsWith(clean$Name, 's') %>% sum() # 114
  clean$Name = sub("s$", "", clean$Name)
  parent_plural = endsWith(clean$Parent_Name, 's') %>% sum() # 317
  clean$Parent_Name = sub("s$", "", clean$Parent_Name)
  if (verbose) {
    print(paste("Plural Names:  ", name_plural))
    print(paste("Plural Parents:", parent_plural))
  }

  return(clean)
}

set_name_case <- function(clean) {
  clean$Name = toupper(clean$Name)
  clean$Parent_Name = toupper(clean$Parent_Name)
  return(clean)
}

CURATE2QUERY =
  c("BIPOLAR CELL", "BC",
    "HORIZONTAL CELL", "HC",
    "AMACRINE CELL", "AC",
    "ROD BCS (RBCS)", "RBC",
    "MULLER GLIA", "MG",
    "ENDOTHELIAL CELL", "ENDOTHELIAL",
    "GANGLION CELL", "RGC",
    "CYCLING MICROGLIA", "CYCLING_MICROGLIA",
    "PAN-BC", "UNASSIGNED_BC"
    # "M3, M5, OR M6", # This needed to be fixed from the %s:/, /,\n/g
) %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }

check_CURATE2QUERY <- function() {
  all(names(CURATE2QUERY) %in% clean_curated$Name)
  all(CURATE2QUERY %in% clean_queried$Name)
}
set_cellname_manually <- function(clean, dictionary, verbose = FALSE) {

  if (verbose) {
    print(length(dictionary)) # 7
    # Test
    print(cbind(from = clean$Name[clean$Name %in% names(dictionary)],
          to   = dictionary[match(clean$Name, names(dictionary)) %>% na.omit()]) %>% data.frame() %>% distinct())
    print(cbind(from = clean$Parent_Name[clean$Parent_Name %in% names(dictionary)],
          to   = dictionary[match(clean$Parent_Name, names(dictionary)) %>% na.omit()]) %>% data.frame() %>% distinct())
  }

    # Do
  clean$Name[clean$Name %in% names(dictionary)] =
    dictionary[match(clean$Name, names(dictionary)) %>% na.omit()]
  clean$Parent_Name[clean$Parent_Name %in% names(dictionary)] =
    dictionary[match(clean$Parent_Name, names(dictionary)) %>% na.omit()]

  return(clean)
}

# Manual Query ----
# Written after ACs were fixed
make_QUERY2CURATE <- function() {
  # View(clean_curated %>% reframe(Name, Parent_Name) %>% distinct())
  clean_queried$Name %>% unique() %>%
    setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()) %>%
    dput()
}
QUERY2CURATE =
c(#"RGC", # There is no plain "RGC" weird in curated
  # "NOVEL_10",
  "OODS_NT_12", "NT-OODSGC", # Many to One, not in curated
  # "NOVEL_13",
  # OODS_CCK_14,
  "OODS_DV_16", "DV-OODSGC", # Many to One, not in curated
  "TBR1_S1_17", "TBR1-S1",
  "W3D1.1_1", "W3D1", # Executive Decision to combine
  "TBR1_S2_21", "TBR1-S2",
  "M5_22", "M3, M5, OR M6", # This one is painful
  "W3D2_23", "W3D2",
  "FMIDIOFF_28", "F-MIDI-OFF",
  "W3D1.2_2", "W3D1", # Executive Decision
  # "NOVEL_30",
  "M2_31", "M2",
  "F_NOVEL_32", "F-NOVEL",
  "M1_33", "M1B", # based solely on expression in the query of Opn4 in the other one
  "FMIDION_38", "F-MIDI-ON",
  "M1DUP_40", "M1A", # based solely on expression in the query of Opn4
  "FMINION_3", "F-MINI-ON",
  "ALPHAONT_41", "ALPHAON-T",
  "ALPHAOFFS_42", "ALPHAOFF-S", # Process of elimination
  "ALPHAONS_43", "ALPHAON-S", # Process of elimination
  "ALPHAOFFT_45", "ALPHAOFF-T",
  "FMINIOFF_4", "F-MINI-OFF",
  "J-RGC_5", "J-RGC",
  "W3B_6", "W3B",
  "TBR1_NOVEL_9", "TBR1-NOVEL",
  "AC12", "MOUSE AC", # Could be wrong... be weary of this one. Double check it later
  "MICROGLIA(CYCLING)", "CYCLING_MICROGLIA")  %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }

check_QUERY2CURATE <- function() {
  all(names(QUERY2CURATE) %in% clean_queried$Name)
  all(QUERY2CURATE %in% c(clean_curated$Name, clean_curated$Parent_Name))
  QUERY2CURATE[!QUERY2CURATE %in% c(clean_curated$Name, clean_curated$Parent_Name)]
}

query_sub <- function(clean) {
  clean$Name = sub("AC_", "AC", clean$Name)
  clean$Name = sub("([0-9]+)_(.*)", "\\2_\\1", clean$Name)

  return(clean)
}

plot_gene_cell_venn <- function(clean_curated, clean_queried, plotPath) {
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
}

# TODO: Adjust for previous pipeline changes
get_major_name <- function(clean, verbose = FALSE) {

  clean$Major_Name = NA

  na_count = c(total = nrow(clean))

  parent_is_majorclass = clean$Parent_Name %in% majorclass
  clean$Major_Name[parent_is_majorclass] = clean$Parent_Name[parent_is_majorclass] %>% toupper()
  na_count = c(na_count, Parent_Name = sum(is.na(clean$Major_Name))) # 89

  is_majorclass = is.na(clean$Major_Name) & clean$Name %in% majorclass
  clean$Major_Name[is_majorclass] = clean$Name[is_majorclass] %>% toupper()
  na_count = c(na_count, Name = sum(is.na(clean$Major_Name))) # 58

  clean$Major_Name[is.na(clean$Major_Name) & grepl("EPITHELIUM", clean$Name)] = "EPITHELIAL"
  na_count = c(na_count, Epithelial = sum(is.na(clean$Major_Name))) # 54

  clean$Major_Name[is.na(clean$Major_Name) & grepl("DENDRITIC CELL|MONOCYTE|T CELL|B CELL|NK CELL|MACROPHAGE|BAM|FIBROBLAST|NEUTROPHIL", clean$Name)] = "IMMUNE"
  na_count = c(na_count, Immune = sum(is.na(clean$Major_Name))) # 25

  clean$Major_Name[is.na(clean$Major_Name) & grepl("ASTROCYTE", clean$Name)] = "ASTROCYTE" #= "MELANOCYTE"
  na_count = c(na_count, Astrocyte = sum(is.na(clean$Major_Name))) # 12

  clean$Major_Name[is.na(clean$Major_Name) & grepl(" MG", clean$Name)] = "MG"
  na_count = c(na_count, MG = sum(is.na(clean$Major_Name))) # 5

  clean$Major_Name[is.na(clean$Major_Name) & grepl("MICROGLIA", clean$Name)] = "MICROGLIA"
  na_count = c(na_count, Microglia = sum(is.na(clean$Major_Name))) # 2

  clean$Major_Name[is.na(clean$Major_Name)] = "MONOCYTE"
  na_count = c(na_count, Monocyte = sum(is.na(clean$Major_Name))) # 0!

  if (verbose) {
    print("Removal by Step:")
    print(na_count)
  }

  return(clean)
}
