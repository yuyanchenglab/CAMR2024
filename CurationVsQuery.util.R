
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
# View(clean_curated %>% reframe(Name, Parent_Name) %>% distinct())
# clean_queried$Name %>% unique() %>%
#   setdiff(c(clean_curated$Parent_Name, clean_curated$Name) %>% unique()) %>%
#   dput()

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
  "ALPHAOFF_42", "ALPHAOFF-S", # Process of elimination
  "ALPHAON_43", "ALPHAON-S", # Process of elimination
  "ALPHAOFFT_45", "ALPHAOFF-T",
  "FMINIOFF_4", "F-MINI-OFF",
  "J-RGC_5", "J-RGC",
  "W3B_6", "W3B",
  "TBR1_NOVEL_9", "TBR1-NOVEL",
  "AC12", "MOUSE AC", # Could be wrong... be weary of this one. Double check it later
  "MICROGLIA(CYCLING)", "CYCLING MICROGLIA")  %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }

# all(names(QUERY2CURATE) %in% clean_queried$Name)
# all(QUERY2CURATE %in% c(clean_curated$Name, clean_curated$Parent_Name))
# QUERY2CURATE[!QUERY2CURATE %in% c(clean_curated$Name, clean_curated$Parent_Name)]
