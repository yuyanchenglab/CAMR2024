# JM 08/27/2024
# Let's analyze this cell chat object.

# Setup ----

library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(dplyr) # pipe
library(data.table) # fast table manip
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(magrittr)
library(ChengLabThemes)

analysisName = "06_Curated_Immune"
analysisPath = "/project/ycheng11lab/jfmaurer/mouse_retina_atlas_chen_2024/"
setwd(analysisPath)

outPath = paste0(analysisPath, analysisName, "/V3_Checks/")
dir.create(outPath, showWarnings = FALSE)

## Reproducibility ----

now = \() Sys.time() %>% format('%F_%H-%M-%S-%Z')
sessionInfoPath = paste0("logs/", analysisName, ".", now(), ".sessionInfo")
sessionInfo() %>% capture.output(file = sessionInfoPath)

set.seed(1)

# Data ----

time_points = setNames(c("Ctrl", "0.5", "1", "2", "4", "7", "14"),
                       c("Ctrl", "0.5dpc",  "01dpc",  "02dpc",  "04dpc",  "07dpc",  "14dpc"))

all_cells = c("Plasma B cell", "B cell", "T cell", # Lymphocyte
                  "Monocyte", "Neutrophil", # Myeloid
                  "Endothelial", "Pericyte", # Vascular and retina below
                  "Astrocyte", "RGC", "Amacrine cell", "Microglia", "Muller glia", "Bipolar cell", "Horizontal cell", "Photoreceptor")

immune_cell_names = c("B CELL", "T CELL", "NK CELL", "MONOCYTE", "NEUTROPHIL", "MACROPHAGE", "DENDRITIC CELL")

ccData = readRDS("/project/ycheng11lab/jfmaurer/pub-rgcs_cellChat/data/checkpoint.RDS")
clean_curated = fread("09_Designer_Analysis/PanelDesignerV3.txt")

# Gene Harmony ----

setdiff(clean_curated$Marker, str_to_title(rownames(ccData))) %>% length() # 8
#[1] "Cd90"     "Vmn2r122" "Lbhd2"    "Tafa4"    "Fads2b"   "Igfbp"    "Ceacam"   "Cd39"
# formatted_genes = str_to_title(rownames(ccData))

clean_curated$Marker[clean_curated$Marker == "Lbhd2"] = "A230065H16Rik" # https://www.phosphosite.org/proteinAction?id=5141911
clean_curated$Marker[clean_curated$Marker == "Tafa4"] = "Fam19a4" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
clean_curated$Marker[clean_curated$Marker == "Fads2b"] = "4833423E24Rik" # https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAFA4
clean_curated$Marker[clean_curated$Marker == "Cd90"] = "Thy1" # "Cd90", https://www.genecards.org/cgi-bin/carddisp.pl?gene=THY1
clean_curated$Marker[clean_curated$Marker == "Cd39"] = "Entpd1" # "Cd39", https://en.wikipedia.org/wiki/ENTPD1
# "Igfbp" which one, Igfbp1 to Igfbp7?
# "Ceacam", again, which one?
clean_curated %<>% filter(!Marker %in% c("Igfbp", "Ceacam", "Vmn2r122"))
# Capitalization
clean_curated$Marker[clean_curated$Marker == "2010007h06rik"] = "2010007H06Rik"

# Major ---

# Harmonize ----

# clean_curated$Major_Name %>% unique() %>% dput()

ccData$Major_Name = toupper(ccData$cell_type_granular) # changed from majorclass on 9/6
# ccData$Major_Name %>% unique() %>% dput()
{
data_harmony =
  c(#"MULLER GLIA",
# "MULLER GLIA (IFN)",
# "ASTROCYTE (IGF2)",
# "MULLER GLIA (REACTIVE)",
#     "ASTROCYTE (IFN)",
# "ASTROCYTE",
# "ENDOTHELIAL",
# "MULLER GLIA/ASTROCYTE (DOUBLET)",
#     "MULLER GLIA/NEURON (DOUBLET)",
"NK CELL", "NK CELL",
"NEUTROPHIL", "NEUTROPHIL",
"LY6C+ MONOCYTE (CHIL3)", "MONOCYTE",
    # "MICROGLIA",
"T CELL", "T CELL",
"B CELL", "B CELL",
"MACROPHAGE (MS4A7 MHC-LO)", "MACROPHAGE",
    # "AMACRINE CELL",
"MACROPHAGE (MS4A7 MHC-HI)", "MACROPHAGE",
# "PHOTORECEPTOR",
#     "ASTROCYTE/MICROGLIA (DOUBLET)",
# "MICROGLIA (IFN)",
# "RGC",
"MON/MAC", "MON/MAC",
    # "MICROGLIA (CYCLING)",
"LY6C- MONOCYTE (EAR2)", "MONOCYTE",
# "MAST CELL",
    "MACROPHAGE (GPNMB)", "MACROPHAGE",
"PLASMA B CELL", "B CELL",
# "PERICYTE",
# "SMOOTH MUSCLE CELL",
    # "BIPOLAR CELL",
# "ASTROCYTE (REACTIVE)",
# "MICROGLIA (CYCLING II)",
    "PDC", "DENDRITIC CELL",
"MO-DC", "DENDRITIC CELL",
# "HORIZONTAL CELL",
"MON/DC", "MON/DC",
"DC (CCR7)", "DENDRITIC CELL"
# "ASTROCYTE (CYCLING)",
#     "EPITHELIAL CELL",
# "3_FMINION",
# "18_NOVEL",
# "23_W3D2",
# "6_W3B",
#     "12_OODS_NT",
# "10_NOVEL",
# "11_NOVEL",
# "31_M2",
# "8_NOVEL",
# "24_NOVEL",
#     "4_FMINIOFF",
# "16_OODS_DV",
# "13_NOVEL",
# "2_W3D1.2",
# "39_NOVEL",
#     "27_NOVEL",
# "33_M1",
# "15_NOVEL",
# "14_OODS_CCK",
# "1_W3D1.1",
# "7_NOVEL",
#     "9_TBR1_NOVEL",
# "5_J-RGC",
# "17_TBR1_S1",
# "29_NOVEL",
# "25_NOVEL",
#     "19_NOVEL",
# "30_NOVEL",
# "26_NOVEL",
# "35_NOVEL",
# "20_NOVEL",
# "21_TBR1_S2",
#     "32_F_NOVEL",
# "44_NOVEL",
# "34_NOVEL",
# "42_ALPHAOFFS",
# "45_ALPHAOFFT",
#     "UNASSIGNED",
# "36_NOVEL",
# "28_FMIDIOFF",
# "38_FMIDION",
# "37_NOVEL",
#     "22_M5",
# "43_ALPHAONS",
# "41_ALPHAONT",
# "40_M1DUP")
) %>%
    matrix(nrow = 2) %>% t() %>%
    { setNames(.[,2], .[,1]) }
}
# select(clean_curated, Name, Major_Name) %>% distinct() %>% filter(Major_Name == "IMMUNE") %>% select(Name) %>% dput()

curated_Major_Name =
  c("B CELL", "B CELL",
    "BAM", "MACROPHAGE",
    "DENDRITIC CELL", "DENDRITIC CELL", # Combine with Macrophage for Monocyte?
    "EAR2+ MONOCYTE", "MONOCYTE",
    "MACROPHAGE", "MACROPHAGE", # Combine with Dendritic for Monocyte?
    "MONOCYTE", "MONOCYTE",
    "NEUTROPHIL", "NEUTROPHIL",
    "PLASMA B CELL", "B CELL",
    "SCLERAL FIBROBLAST", "IMMUNE",
    "T CELL", "T CELL",
    "T/NK CELL", "T CELL",
    "CONE", "PHOTORECEPTOR", # Legit?
    "ROD", "PHOTORECEPTOR") %>%
  matrix(nrow = 2) %>% t() %>%
  { setNames(.[,2], .[,1]) }

clean_curated$Major_Name[clean_curated$Name %in% names(curated_Major_Name)] =
  curated_Major_Name[clean_curated$Name[clean_curated$Name %in% names(curated_Major_Name)]]

ccData$Major_Name[ccData$Major_Name %in% names(data_harmony)] =
  data_harmony[ccData$Major_Name[ccData$Major_Name %in% names(data_harmony)]]

# mnd = unique(ccData$Major_Name)
# mnc = unique(clean_curated$Major_Name)
# setdiff(mnd, mnc)
# setdiff(mnc, mnd)
# intersect(mnc, mnd) # majorclass_for_plotting

# Filter ----

which_correct = ccData$majorclass %in% c("RGC", "Astrocyte", "Muller glia")
ccData$Major_Name[which_correct] = toupper(ccData$majorclass[which_correct])
ccData$Major_Name[grep("MICROGLIA", ccData$Major_Name)] = "MICROGLIA"
ccData = subset(ccData, subset = Major_Name %notin% c("HORIZONTAL CELL", "MAST CELL", "SMOOTH MUSCLE CELL", "EPITHELIAL CELL"))

# Plot ----

p = DotPlot(ccData, features = unique(clean_curated$Marker), group.by = "Major_Name") +
  SparseBubbleTheme() +
  labs(x = "", y = "") +
  theme(axis.ticks = element_blank(),
        panel.grid = element_line(color = "gray90", linetype = "dashed"))
ggsave(paste0(outPath, "6_dotplot_majorclass_unfiltered_normExpression.pdf"), p, height = 6, width = 50, limitsize = FALSE)

p = DotPlot(ccData, features = unique(clean_curated$Marker), group.by = "Major_Name", scale = FALSE) +
  SparseBubbleTheme() +
  labs(x = "", y = "") +
  theme(axis.ticks = element_blank(),
        panel.grid = element_line(color = "gray90", linetype = "dashed")) +
  scale_color_continuous(limits = c(0,4), oob = scales::oob_squish) +
  guides(color = guide_colorbar(title = "Average Raw\nExpression"))
ggsave(paste0(outPath, "6_dotplot_majorclass_unfiltered_rawExpression.pdf"), p, height = 6, width = 50, limitsize = FALSE)

# Scratch ----

q()

ccData$Major_Name = toupper(ccData$majorclass)
original_harmonizing_table =
c("MULLER GLIA", "MG",
  # "ASTROCYTE",
  # "ENDOTHELIAL",
  # "T CELL",
  # "NEUTROPHIL",
  # "MONOCYTE",
  # "MICROGLIA",
  # "B CELL",
  "AMACRINE CELL", "AC",
  # "PHOTORECEPTOR",
  # "RGC",
  "PLASMA B CELL", "B CELL",
  # "PERICYTE",
  "BIPOLAR CELL", "BC",
  "HORIZONTAL CELL", "HC") %>%
