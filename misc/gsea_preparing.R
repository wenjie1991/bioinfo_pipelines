library(stringr)

######
#  README  #
######


## read expression tables
d = fread("../data/table.exp/tablecounts_tpm.csv")

## get sample names
sample_names = paste0(d[2, -1] %>% unlist, "_"  ,d[1, -1] %>% unlist)
names(d)[-1] = sample_names
names(d)[1] = "symbol"

## reorder samples & gene annotation
m_unpair = d[-c(1:2), c(1, 1:10)]
names(m_unpair)[1:2] = c("NAME", "DESCRIPTION")

## removing low expressed genes
m_unpair = m_unpair[M1_cMT > 0 | M1_cMWT > 0 | M2_cMT > 0 | M2_cMWT > 0 | M3_cMT > 0 | M3_cMWT > 0 | M4_T > 0 | M5_T > 0 | M6_T > 0]
write_tsv(m_unpair, "../data/GSEA/un_paired_expression_matrix_9_samples.txt")

## prepare phenotype
phenotype_unpair = "6 2 1\n# cMT cMWT\ncMT cMWT cMT cMWT cMT cMWT\n"
write(phenotype_unpair, "../data/GSEA/un_paired_expression_phenotype.cls")

## prepare custom gene signature
mRepair = fread("../data/mRepair_all.tsv.txt")
mFetal = fread("../data/mFetal_all.tsv.txt")

x_Repair_up = mRepair[DSS.Sca.Plus.Non.DSS.fdr < 0.05][DSS.Sca.Plus.Non.DSS.logFC > 2]
x_Repair_down = mRepair[DSS.Sca.Plus.Non.DSS.fdr < 0.05][DSS.Sca.Plus.Non.DSS.logFC < -2]

r_up = x_Repair_up$GeneName %>% na.omit %>% paste(collapse = "\t") %>% paste0("REPAIR_UP\t\t", .)
r_down = x_Repair_down$GeneName %>% na.omit %>% paste(collapse = "\t") %>% paste0("REPAIR_DOWN\t\t", .)
write(c(r_up, r_down), "../data/mRepair.gmt", sep = "\n")
