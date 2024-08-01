
# Let's see everything at once
# View(clean_queried$Name %>% unique() %>% as.data.frame())
# clean_curated$Name %>% unique() %>% dput()

# One Side ----

CURATE2QUERY =
  c("BIPOLAR CELL", "BC",
    "HORIZONTAL CELL", "HC",
    "AMACRINE CELL", "AC",
    "ROD BCS (RBCS)", "RBC",
    "MULLER GLIA", "MG",
    "ENDOTHELIAL CELL", "ENDOTHELIAL",
    "GANGLION CELL", "GC"
    # "M3, M5, OR M6", # This needed to be fixed from the %s:/, /,\n/g
) %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }

# Did I mispell something?
# all(names(CURATE2QUERY) %in% clean_curated$Name)
# all(CURATE2QUERY %in% clean_queried$Name)

# The Other Side ----
# Written after ACs were fixed
View(clean_curated %>% reframe(Name, Parent_Name) %>% distinct())
clean_queried$Name %>% unique() %>%
  setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()) %>%
  dput()

QUERY2CURATE =
c(#"RGC", "GC", # There is no plain "RGC" weird
  #"OODS_NT_12", "OODSGC", Many to One
  #"OODS_DV_16", "OODSGC", Many to One
  "W3D1.1_1", "W3D1",
  "TBR1_S2_21", "TBR1-S1",
  "M5_22", "M3, M5, OR M6",
  "W3D2_23", "W3D2",
  "FMIDIOFF_28", "F-MIDI-OFF",
  "W3D1.2_2", "W3D1",
  "M2_31", "M2",
  "F_NOVEL_32", "F-NOVEL",
  # "M1_33", Many to One
  "FMINION_3", "F-MIDI-ON",
  "ALPHAONT_41", "ALPHAON-T",
  "ALPHAOFFS_42", "ALPHAOFF-S",
  "ALPHAONS_43", "ALPHAON-S",
  "ALPHAOFFT_45", "ALPHAOFF-T",
  "W3B_6", "W3B",
  "TBR1_NOVEL_9", "TBR1-NOVEL",
  "AC12", "MOUSE AC", # Could be wrong... be weary of this one. Double check it later
  "MICROGLIA(CYCLING)", "CYCLING MICROGLIA")  %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }
