BiocManager::install("TCGAbiolinks") # fonctionne pas
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

BiocManager::install('GenomicDataCommons')
library(GenomicDataCommons)


q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-GBM' &
           data_type=='Gene Expression Quantification')
q %>% facet('analysis.workflow_type') %>% aggregations()


