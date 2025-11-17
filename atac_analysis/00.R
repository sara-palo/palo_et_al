setwd("~/work/FL_adult_ATAC/2024_FL_ATAC")

set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'dplyr', 'ComplexHeatmap', 'GenomicRanges', 'EnsDb.Hsapiens.v86') 
invisible(lapply(packages, library, character.only = TRUE))

rawdata_path <- '../merged_FL12_BM12/outs/'
rds_path <- 'rds_files/'
fig_path <- 'figs/'
csv_path <- 'csv_files/'
ext_data_path <- 'external_data/'

csv_path_da <- paste0(csv_path, 'da_peaks')
csv_path_ch <- paste0(csv_path, 'da_chromvar')
sup_table_path <- paste0(csv_path, 'sup')

# create folders for output
dirs <- c(rds_path, fig_path, csv_path, csv_path_ch, csv_path_da, sup_table_path)

for(dir in dirs){
  if(! dir.exists(dir)){
    dir.create(dir)
  }
}

main_figs <- paste0('fig', 1:4)
sup_figs <- paste0('sup', 1:4)

for(i in c(main_figs, sup_figs, 'additional_plots')){
  dir <- paste0(fig_path, i)
  if(! dir.exists(dir)){
    dir.create(dir)
  }
}

for(i in c('vs_all', 'vs_prim', 'FLvsBM')){
  dir <- paste0(csv_path_da, '/', i)
  if(! dir.exists(dir)){
    dir.create(dir)
  }
}

for(i in c('vs_all', 'vs_prim')){
  dir <- paste0(csv_path_ch, '/', i)
  if(! dir.exists(dir)){
    dir.create(dir)
  }
}

# plot colors
annot_colors <- list(
  cluster_cols = c(
    'HSC' = '#E9CB53',
    'MPP-I' = '#52BE80',
    'MPP-III' = '#16A085',
    'Ly-I' = '#A569BD',
    'Ly-II' = '#793885',
    'DC' = '#4e3666',
    'GMP-II' = '#2874A6',
    'GMP-I' = '#56ace5',
    'MPP-II' = '#E7A235',
    'MEP-I' = '#DC7633',
    'MEP-II' = '#E46457',
    'Ery-I' = '#B03A2E',
    'Ery-II' = '#78281F',
    'Ba/M' = '#F06292'
  ),
  tissue_cols = c('FL' = '#1c6494', 'BM' = '#1b4b1f'),
  sample_numberFL = c('FL1' = '#93bedd', 'FL2' = '#1c6494'),
  sample_numberBM = c('BM1' = '#78a177', 'BM2' = '#1b4b1f'),
  sample_number2nd = c('h5a' = '#b3ccf3', 'h5b' = '#3fa2ba', 'h6a' = '#0f3f5e'),
  da_cols =
    c(
      'down' = '#081b80',
      'neither' = 'black',
      'up' = '#840000'
    ),
  da_cols2 =
    c(
      'down in BM' = '#081b80',
      'up in BM' = '#840000',
      'down in FL' = '#081b80', 
      'up in FL' = '#840000' 
    ),
  heatmap_cols = c('#9E1000', '#FFFEF9'),
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
  ),
  expr_cols_1 = c('lightgrey', 'navy'),
  expr_cols_2 = c('white', 'navy'),
  cite_cols = c(
    'HSC' = '#E9CB53',
    'MPP-II' = '#52BE80',
    'Ly-I' = '#A569BD',
    'Ly-II' = '#793885',
    'DC-I' = '#58558f',
    'Ly-III' = '#4e3666',
    'DC-Mono' = '#2874A6',
    'GMP' = '#56ace5',
    'MPP-I' = '#E7A235',
    'MEP' = '#78281F',
    'Cyc' = '#533511'
  )
)

clust_order <- c('HSC', 'MPP-I', 'MPP-II', 'MPP-III', 'Ly-I', 'Ly-II', 'DC', 'GMP-I', 'GMP-II', 'MEP-I', 'MEP-II', 'Ery-I', 'Ery-II', 'Ba/M')
clust_order_CITE <- c('HSC', 'MPP-I', 'MPP-II', 'Ly-I', 'Ly-II', 'Ly-III', 'DC-I', 'Cyc', 'DC-Mono', 'GMP', 'MEP')

source('scripts/helpers.R')
