options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_ROOT <- Sys.getenv(
  "THESIS_ROOT",
  "/Users/mac/Desktop/TCGA_Thesis"
)

RESULT_DIR <- file.path(PROJECT_ROOT, "results")
TABLE_DIR <- file.path(RESULT_DIR, "tables")
FIGURE_DIR <- file.path(RESULT_DIR, "figures")
RDATA_DIR <- file.path(RESULT_DIR, "rdata")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDATA_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ReactomePA",
  "DOSE",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "pheatmap",
  "RColorBrewer",
  "minfi",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：",
    paste(missing_packages, collapse = ", ")
  )
}

library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(DOSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ==================================================================
# 一、CpG位点与基因匹配
# ==================================================================

sig_probes <- readRDS(
  file.path(
    RDATA_DIR,
    "limma_deltabeta_sig_probes.rds"
  )
)

epic_anno <- getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19
)

hm450_anno <- getAnnotation(
  IlluminaHumanMethylation450kanno.ilmn12.hg19
)

extract_annotation <- function(annotation) {
  annotation <- as.data.frame(annotation) %>%
    tibble::rownames_to_column("ProbeID")

  annotation %>%
    transmute(
      ProbeID = ProbeID,
      chr = chr,
      start = pos,
      end = pos,
      Gene = UCSC_RefGene_Name,
      Gene_Region = UCSC_RefGene_Group,
      CpG_Island = Relation_to_Island
    )
}

combined_anno <- bind_rows(
  extract_annotation(epic_anno),
  extract_annotation(hm450_anno)
) %>%
  distinct(ProbeID, .keep_all = TRUE)

cpg_gene_map <- combined_anno %>%
  filter(ProbeID %in% sig_probes) %>%
  mutate(
    Gene = ifelse(
      Gene == "",
      NA_character_,
      Gene
    )
  ) %>%
  separate_rows(
    Gene,
    sep = ";"
  ) %>%
  filter(!is.na(Gene)) %>%
  transmute(
    CpG_Site = ProbeID,
    Chromosome = chr,
    Position = start,
    Target_Gene = Gene,
    Gene_Region = Gene_Region
  ) %>%
  arrange(match(CpG_Site, sig_probes)) %>%
  distinct(
    Target_Gene,
    .keep_all = TRUE
  )

write.csv(
  cpg_gene_map,
  file.path(
    TABLE_DIR,
    "Top_CpG_Core_Gene_Map_hg19.csv"
  ),
  row.names = FALSE
)

saveRDS(
  cpg_gene_map,
  file.path(
    RDATA_DIR,
    "Top_CpG_Core_Gene_Map_hg19.rds"
  )
)

# ==================================================================
# 二、基因ID转换
# ==================================================================

core_genes <- unique(
  cpg_gene_map$Target_Gene
)

gene_ids <- bitr(
  core_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

cpg_gene_map <- cpg_gene_map %>%
  left_join(
    gene_ids,
    by = c("Target_Gene" = "SYMBOL")
  ) %>%
  rename(
    Entrez_ID = ENTREZID
  )

write.csv(
  cpg_gene_map,
  file.path(
    TABLE_DIR,
    "Top_CpG_Core_Gene_Map_hg19.csv"
  ),
  row.names = FALSE
)

entrez_ids <- unique(
  na.omit(cpg_gene_map$Entrez_ID)
)

# ==================================================================
# 三、GO和KEGG富集分析
# ==================================================================

go_enrich <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr",
  minGSSize = 1,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = FALSE
)

kegg_enrich <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pAdjustMethod = "fdr",
  minGSSize = 1,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  use_internal_data = FALSE
)

write.csv(
  as.data.frame(go_enrich),
  file.path(
    TABLE_DIR,
    "GO_Enrichment_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(kegg_enrich),
  file.path(
    TABLE_DIR,
    "KEGG_Enrichment_Results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  go_enrich,
  file.path(
    RDATA_DIR,
    "GO_Enrichment_Results.rds"
  )
)

saveRDS(
  kegg_enrich,
  file.path(
    RDATA_DIR,
    "KEGG_Enrichment_Results.rds"
  )
)

# ==================================================================
# 四、Reactome和DO富集
# ==================================================================

reactome_enrich <- enrichPathway(
  gene = entrez_ids,
  organism = "human",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  minGSSize = 1,
  readable = TRUE
)

do_enrich <- enrichDO(
  gene = entrez_ids,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(
  as.data.frame(reactome_enrich),
  file.path(
    TABLE_DIR,
    "Reactome_Enrichment_Results.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(do_enrich),
  file.path(
    TABLE_DIR,
    "DO_Enrichment_Results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  reactome_enrich,
  file.path(
    RDATA_DIR,
    "Reactome_Enrichment_Results.rds"
  )
)

saveRDS(
  do_enrich,
  file.path(
    RDATA_DIR,
    "DO_Enrichment_Results.rds"
  )
)

# ==================================================================
# 五、GSEA
# ==================================================================

load(
  file.path(
    RDATA_DIR,
    "limma_analysis_objects.RData"
  )
)

results_limma <- readRDS(
  file.path(
    RDATA_DIR,
    "limma_deltabeta_full_results.rds"
  )
)

results_limma <- results_limma %>%
  filter(
    Probe_ID %in% sig_probes
  )

gene_level_stats <- results_limma %>%
  left_join(
    combined_anno,
    by = c("Probe_ID" = "ProbeID")
  ) %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  separate_rows(
    Gene,
    sep = ";"
  ) %>%
  group_by(Gene) %>%
  slice_min(
    order_by = adj.P.Val,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

gsea_gene_ids <- bitr(
  gene_level_stats$Gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

gene_list <- gene_level_stats %>%
  inner_join(
    gsea_gene_ids,
    by = c("Gene" = "SYMBOL")
  ) %>%
  group_by(ENTREZID) %>%
  slice_max(
    order_by = abs(logFC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  pull(
    logFC,
    ENTREZID
  ) %>%
  sort(
    decreasing = TRUE
  )

gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 1,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  verbose = FALSE
)

write.csv(
  as.data.frame(gsea_kegg),
  file.path(
    TABLE_DIR,
    "KEGG_GSEA_Results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  gsea_kegg,
  file.path(
    RDATA_DIR,
    "KEGG_GSEA_Results.rds"
  )
)

# ==================================================================
# 六、免疫相关GO/KEGG热图
# ==================================================================

load(
  file.path(
    RDATA_DIR,
    "train_data.RData"
  )
)

train_data <- train_data_final

val_data <- readRDS(
  file.path(
    RDATA_DIR,
    "val_limma_filtered.rds"
  )
)

load(
  file.path(
    RDATA_DIR,
    "GSE119144_test_data.RData"
  )
)

test_data <- final_df

process_heatmap_data <- function(data, probes) {
  common_probes <- intersect(
    probes,
    colnames(data)
  )

  data %>%
    select(
      Sample_Name,
      Sample_Group,
      all_of(common_probes)
    ) %>%
    mutate(
      Sample_Group = factor(
        ifelse(
          Sample_Group %in% c("R", 1, "1"),
          "Responder",
          "Non-responder"
        )
      )
    )
}

all_samples_df <- bind_rows(
  process_heatmap_data(train_data, sig_probes),
  process_heatmap_data(val_data, sig_probes),
  process_heatmap_data(test_data, sig_probes)
)

immune_keywords <- c(
  "immun", "T cell", "B cell", "leukocyte",
  "lymphocyte", "cytokine", "chemokine",
  "interleukin", "interferon", "antigen",
  "MHC", "HLA", "NK cell", "natural killer",
  "macrophage", "dendritic", "monocyte",
  "complement", "toll-like", "NOD-like",
  "RIG-I", "inflammasome", "Th1", "Th2",
  "Th17", "Treg", "cytotox", "autoimmun",
  "inflammatory", "hematopoietic",
  "antibody", "immunoglobulin", "PD-1",
  "PD-L1", "CTLA", "checkpoint",
  "CD4", "CD8", "TCR", "BCR",
  "phagocyt", "adaptive immune", "innate immune",
  "immune response", "cell killing", "NF-kappaB"
)

immune_pattern <- paste(
  immune_keywords,
  collapse = "|"
)

select_unique_probe <- function(pathway_row,
                                core_gene_map,
                                used_probes = character()) {
  genes_in_pathway <- unlist(
    strsplit(
      as.character(pathway_row$geneID),
      "/",
      fixed = TRUE
    )
  )

  candidates <- core_gene_map %>%
    filter(
      Entrez_ID %in% genes_in_pathway
    ) %>%
    pull(CpG_Site) %>%
    unique()

  candidates <- intersect(
    candidates,
    colnames(all_samples_df)
  )

  candidates <- setdiff(
    candidates,
    used_probes
  )

  if (length(candidates) == 0) {
    return(list(
      probe = NA_character_,
      used_probes = used_probes
    ))
  }

  list(
    probe = candidates[1],
    used_probes = c(
      used_probes,
      candidates[1]
    )
  )
}

make_heatmap <- function(probe_map,
                         data,
                         filename,
                         main_title,
                         row_annotation_columns) {
  if (nrow(probe_map) == 0) {
    return(invisible(NULL))
  }

  probes <- unique(probe_map$Probe)

  heatmap_matrix <- data %>%
    select(
      Sample_Name,
      all_of(probes)
    ) %>%
    tibble::column_to_rownames("Sample_Name") %>%
    as.matrix() %>%
    t()

  annotation_col <- data %>%
    select(
      Sample_Name,
      Sample_Group
    ) %>%
    distinct(Sample_Name, .keep_all = TRUE) %>%
    tibble::column_to_rownames("Sample_Name")

  annotation_row <- probe_map %>%
    select(
      Probe,
      all_of(row_annotation_columns)
    ) %>%
    tibble::column_to_rownames("Probe")

  annotation_colors <- list(
    Sample_Group = c(
      Responder = "#E1341E",
      `Non-responder` = "#3B8AB8"
    )
  )

  for (annotation_name in row_annotation_columns) {
    values <- unique(
      annotation_row[[annotation_name]]
    )

    colors <- colorRampPalette(
      brewer.pal(8, "Set2")
    )(length(values))

    annotation_colors[[annotation_name]] <-
      setNames(colors, values)
  }

  pdf(
    filename,
    width = 14,
    height = 8
  )

  pheatmap(
    heatmap_matrix,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    main = main_title,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    fontsize_row = 9,
    color = colorRampPalette(
      c("navy", "white", "firebrick3")
    )(50)
  )

  dev.off()
}

core_gene_map <- cpg_gene_map

# -------------------- KEGG免疫通路 --------------------

kegg_df <- as.data.frame(kegg_enrich)

kegg_immune <- kegg_df %>%
  filter(
    grepl(
      immune_pattern,
      Description,
      ignore.case = TRUE
    )
  )

kegg_probe_map <- data.frame()

used_probes <- character(0)

for (i in seq_len(nrow(kegg_immune))) {
  selected <- select_unique_probe(
    kegg_immune[i, ],
    core_gene_map,
    used_probes
  )

  if (!is.na(selected$probe)) {
    kegg_probe_map <- bind_rows(
      kegg_probe_map,
      data.frame(
        Probe = selected$probe,
        Pathway = kegg_immune$Description[i],
        stringsAsFactors = FALSE
      )
    )

    used_probes <- selected$used_probes
  }
}

write.csv(
  kegg_probe_map,
  file.path(
    TABLE_DIR,
    "KEGG_Immune_Probe_Map.csv"
  ),
  row.names = FALSE
)

make_heatmap(
  kegg_probe_map,
  all_samples_df,
  file.path(
    FIGURE_DIR,
    "KEGG_Immune_Probes_Heatmap.pdf"
  ),
  "Methylation Heatmap of Immune-Related KEGG Pathway Probes",
  "Pathway"
)

# -------------------- GO免疫条目 --------------------

go_df <- as.data.frame(go_enrich)

go_immune <- go_df %>%
  filter(
    grepl(
      immune_pattern,
      Description,
      ignore.case = TRUE
    )
  ) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup()

go_probe_map <- data.frame()
used_probes <- character(0)

for (i in seq_len(nrow(go_immune))) {
  selected <- select_unique_probe(
    go_immune[i, ],
    core_gene_map,
    used_probes
  )

  if (!is.na(selected$probe)) {
    go_probe_map <- bind_rows(
      go_probe_map,
      data.frame(
        Probe = selected$probe,
        Ontology = go_immune$ONTOLOGY[i],
        Term = go_immune$Description[i],
        stringsAsFactors = FALSE
      )
    )

    used_probes <- selected$used_probes
  }
}

write.csv(
  go_probe_map,
  file.path(
    TABLE_DIR,
    "GO_Immune_Probe_Map.csv"
  ),
  row.names = FALSE
)

make_heatmap(
  go_probe_map,
  all_samples_df,
  file.path(
    FIGURE_DIR,
    "GO_Immune_Probes_Heatmap.pdf"
  ),
  "Methylation Heatmap of Immune-Related GO Term Probes",
  c("Ontology", "Term")
)

# ==================================================================
# 七、CpG基因组区域富集
# ==================================================================

load_annotation_for_enrichment <- function(annotation) {
  as.data.frame(annotation) %>%
    tibble::rownames_to_column("ProbeID") %>%
    transmute(
      ProbeID = ProbeID,
      chr = chr,
      start = pos,
      end = pos,
      Gene = UCSC_RefGene_Name,
      Gene_Region = UCSC_RefGene_Group,
      CpG_Island = Relation_to_Island
    )
}

annotation_for_enrichment <- bind_rows(
  load_annotation_for_enrichment(epic_anno),
  load_annotation_for_enrichment(hm450_anno)
) %>%
  distinct(ProbeID, .keep_all = TRUE)

train_filtered <- readRDS(
  file.path(
    RDATA_DIR,
    "train_limma_filtered.rds"
  )
)

background_probe_matrix <- train_filtered %>%
  select(
    -Sample_Name,
    -Sample_ID,
    -Sample_Group
  )

background_cv <- apply(
  background_probe_matrix,
  2,
  function(x) {
    sd(x, na.rm = TRUE) /
      mean(x, na.rm = TRUE)
  }
)

background_probes <- names(
  background_cv[background_cv > 0.1]
)

sig_annotation <- annotation_for_enrichment %>%
  filter(ProbeID %in% sig_probes)

background_annotation <- annotation_for_enrichment %>%
  filter(ProbeID %in% background_probes)

sig_gene_region <- sig_annotation %>%
  select(ProbeID, Gene_Region) %>%
  separate_rows(Gene_Region, sep = ";") %>%
  count(Gene_Region, name = "count_sig") %>%
  mutate(
    freq_sig = count_sig / length(sig_probes)
  )

background_gene_region <- background_annotation %>%
  select(ProbeID, Gene_Region) %>%
  separate_rows(Gene_Region, sep = ";") %>%
  count(Gene_Region, name = "count_bg") %>%
  mutate(
    freq_bg = count_bg / length(background_probes)
  )

gene_region_enrichment <- full_join(
  sig_gene_region,
  background_gene_region,
  by = "Gene_Region"
) %>%
  mutate(
    across(
      everything(),
      ~ replace_na(.x, 0)
    )
  ) %>%
  mutate(
    EnrichmentFold = freq_sig / freq_bg
  ) %>%
  filter(
    !is.na(Gene_Region),
    Gene_Region != ""
  )

sig_island <- sig_annotation %>%
  count(
    CpG_Island,
    name = "count_sig"
  ) %>%
  mutate(
    freq_sig = count_sig / length(sig_probes)
  )

background_island <- background_annotation %>%
  count(
    CpG_Island,
    name = "count_bg"
  ) %>%
  mutate(
    freq_bg = count_bg / length(background_probes)
  )

cpg_island_enrichment <- full_join(
  sig_island,
  background_island,
  by = "CpG_Island"
) %>%
  mutate(
    across(
      everything(),
      ~ replace_na(.x, 0)
    )
  ) %>%
  mutate(
    Group = case_when(
      grepl("Shore", CpG_Island) ~ "Shore",
      grepl("Shelf", CpG_Island) ~ "Shelf",
      TRUE ~ CpG_Island
    )
  ) %>%
  group_by(Group) %>%
  summarise(
    freq_sig = sum(freq_sig),
    freq_bg = sum(freq_bg),
    .groups = "drop"
  ) %>%
  mutate(
    EnrichmentFold = freq_sig / freq_bg
  ) %>%
  filter(
    !is.na(Group),
    Group != ""
  )

write.csv(
  gene_region_enrichment,
  file.path(
    TABLE_DIR,
    "CpG_Genomic_Feature_Enrichment.csv"
  ),
  row.names = FALSE
)

write.csv(
  cpg_island_enrichment,
  file.path(
    TABLE_DIR,
    "CpG_Island_Relation_Enrichment.csv"
  ),
  row.names = FALSE
)

# ==================================================================
# 八、GSEA、Reactome和DO结果保存
# ==================================================================

message("07_enrichment.R运行完成。")