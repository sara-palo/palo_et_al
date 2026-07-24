library(easyshiny)
library(Seurat)
library(showtext)

setwd('~/work/FL_adult_ATAC/2024_FL_ATAC/scripts/shiny_app')
rds_path <- '../../rds_files/'

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
FL

DimPlot(FL)

FL[['peaks']] <- NULL

FL[['chromvar']] <- FL[['cisbpchromvar']]
DefaultAssay(FL) <- 'chromvar'
FL[['cisbpchromvar']] <- NULL
FL[['cvejic_umap']] <- NULL
FL[['CITEumap']] <- NULL

FL

md <- FL@meta.data
colnames(md)

md.keep <- md[,c(2:3, 10, 18, 20, 22, 26, 16, 34:35, 38:43, 50:54)]
FL@meta.data <- md.keep

FL_conf <- create_config(FL)

annot_colors <- list(
  cluster_cols = c(
    'HSC' = '#E9CB53',
    'MPP-I' = '#52BE80',
    'MPP-II' = '#E7A235',
    'MPP-III' = '#16A085',
    'Ly-I' = '#A569BD',
    'Ly-II' = '#793885',
    'DC' = '#4e3666',
    'GMP-I' = '#56ace5',
    'GMP-II' = '#2874A6',
    'MEP-I' = '#DC7633',
    'MEP-II' = '#E46457',
    'Ery-I' = '#B03A2E',
    'Ery-II' = '#78281F',
    'Ba/M' = '#F06292'
  ),
  sample_numberFL = c('FL1' = '#93bedd', 'FL2' = '#1c6494'),
  seurat_cluster_cols = c(
    '0' = '#E9CB53',
    '1' = '#52BE80',
    '2' = '#16A085',
    '3' = '#A569BD',
    '4' = '#793885',
    '5' = '#E7A235',
    '6' = '#3498DB',
    '7' = '#58558f',
    '8' = '#2874A6',
    '9' = '#F06292',
    '10' = '#DC7633',
    '11' = '#E46457',
    '12' = '#78281F',
    '13' = '#4e342e'
  )
)

FL_conf <- mod_colour(FL_conf, m = 'seurat_clusters', nc = unname(annot_colors$seurat_cluster_cols))
FL_conf <- mod_colour(FL_conf, m = 'pred.id.BM', nc = unname(annot_colors$cluster_cols))
FL_conf <- mod_colour(FL_conf, m = 'sample_number', nc = unname(annot_colors$sample_numberFL))

# manually edited the ui.R and server.R scripts afterwards
# make_app(FL, FL_conf, gex.assay = 'chromvar', gex.slot = 'data', shiny.title = 'scATAC-seq of 1st trimester FL', 
#          shiny.dir = '/Users/sa5488pa/work/FL_adult_ATAC/2024_FL_ATAC/scripts/shiny_app')

make_file(FL, FL_conf, gex.assay = 'chromvar', gex.slot = 'data', default.multigene = 
            c('HOXB7', 'EBF1', 'TCF3', 'CEBPA', 'CEBPB', 'GATA1', 'NFE2', 'IRF8', 'HOXB8', 'KLF5'), 
         default.gene1 = 'EBF1', default.gene2 = 'CEBPA', 
         shiny.dir = '/Users/sa5488pa/work/FL_adult_ATAC/2024_FL_ATAC/scripts/shiny_app')

list.files('.')


