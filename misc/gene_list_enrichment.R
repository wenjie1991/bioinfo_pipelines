library(stringr)
library(magrittr)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)

#' Transform gene symbol vector to entrez id
#'
#' @param x A character scalar. Contains gene symbol list.g
#' @return A character scalar, contains gene entrez id. Be careful, the order may change.
symbol2entrez = function(x) {
    ensembl = str_split_fixed(x, "\\|", 2) %>% extract(, 1) %>% str_split_fixed("\\.", 2) %>% extract(, 1)
    select(org.Mm.eg.db, ensembl, "ENTREZID", "SYMBOL")[["ENTREZID"]] %>% na.omit
}

#' Perform function enrichment analysis
#'
#' @param query A string scalr, contains gene entrezid.
#' @param species A string, defines the species database to use.
#' @param background A string scalar, optional, defines the background gene list.
#' @return A data.table.
#'   DataBase: where the gene set comes from
#'   ID: gene set id
#'   Description Detail description of gene set
#'   GeneRatio Query Gene in pathway / query gene 
#'   BgRatio Pathway gene / background
#'   pvalue
#'   p.adjust
#'   qvalue
#'   geneID
#'   Count Qury gene in pathway
gene_set_enrichment = function(query, species = "mouse", background = NULL) {

  # TODO: Human datasets
  species_dict = c("mouse" = "Mus musculus", "human" = "Homo ...")
  species_official_name = species_dict[species]

  m_df = msigdbr(species = species_official_name, category = "H")
  m_HGENE = m_df %>% dplyr::select(gs_id, entrez_gene) %>% as.data.frame()
  m_HNAME = m_df %>% dplyr::select(gs_id, gs_name) %>% as.data.frame()
  #   m_df = msigdbr(species = species_official_name, category = "C2")
  #   m_c2GENE = m_df %>% dplyr::select(gs_id, entrez_gene) %>% as.data.frame()
  #   m_c2NAME = m_df %>% dplyr::select(gs_id, gs_name) %>% as.data.frame()
  #   m_df = msigdbr(species = species_official_name, category = "C4")
  #   m_c4GENE = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
  #   m_c4NAME = m_df %>% dplyr::select(gs_id, gs_name) %>% as.data.frame()

  enrichemnt = list()
  # TODO: auto transform the gene name type
  entrezid = query

  if (is.null(background)) {
    enrichemnt$KEGG = enrichKEGG(entrezid, organism = 'mmu', pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, universe = background)
    enrichemnt$HALLMARK = enricher(entrezid, universe = background, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = m_HGENE, TERM2NAME = m_HNAME)
    #  enrichemnt$C4 = enricher(entrezid, universe = background, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = m_c4GENE, TERM2NAME = m_c4NAME)
    #  enrichemnt$C2 = enricher(entrezid, universe = background, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = m_c2GENE, TERM2NAME = m_c2NAME)
    enrichemnt = lapply(enrichemnt, function(x) {setReadable(x, OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>% data.frame}) %>% ldply(.id = "DataBase") %>% data.table
  }
  enrichemnt
}

plot_enrichment_result = function(x) {
  g = ggplot(x) + aes(x = -log10(pvalue), y = Description, color = Count) + geom_point()
  g + theme_bw() + scale_color_continuous()
}


## EXAMPLE:
# d = fread("../data/As_result/WATER_VS_DSS_AS_diffpsi0.1.tsv")
# up_symbol_v = d[change_type == "up", symbol] %>% symbol2entrez
# down_symbol_v = d[change_type == "down", symbol] %>% symbol2entrez
# both_v = d$symbol %>% symbol2entrez

# dir.create("../data/GeneSet_result")
# x = gene_set_enrichment(up_symbol_v)[pvalue <= 0.1]
# g1 =plot_enrichment_result(x) + labs(title = "UP genes")
# write_tsv(x, "../data/GeneSet_result/Water_VS_DSS_UP_pathway.tsv")
