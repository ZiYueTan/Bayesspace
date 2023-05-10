library(usethis)
library(assertthat)
library(reticulate)
library(Giotto)
library(GiottoData)
library(BayesSpace)

# Discover isolated spots in sce before creating giotto object,
# or "key error" occurs when using "do_HMRF" function due to some spots lacking neighbors
find.neighbors.iso <- function(sce, platform){
  if (platform == "ST"){
    #ST platform: squares
    #neighbors: below, above, right, left
    origin <- data.frame(origin.row = c(-1, 1, 0, 0),
                         origin.col = c(0, 0, 1, -1))
  }
  if (platform == "Visium")
    #Visium platform: hexagon
    #neighbors: left, right, upper left, upper right, lower left, lower right
    origin <- data.frame(origin.row = c(0, 0, 1, 1, -1, -1), 
                         origin.col =  c(-2, 2, -1, 1, -1, 1))
  
  
  #get coordinates
  positions <- as.data.frame(colData(sce)[, c("row", "col")])
  positions$spot.index <- c(1 : nrow(positions))
  n <- nrow(positions)
  
  #get each possible spot neighbor
  neighbor.cor <- merge(positions, origin)
  neighbor.cor$neighbor.row <- neighbor.cor$row + neighbor.cor$origin.row
  neighbor.cor$neighbor.col <- neighbor.cor$col + neighbor.cor$origin.col
  
  #get each possible spot neighbor
  neighbor.cor <- merge(positions, origin)
  neighbor.cor$neighbor.row <- neighbor.cor$row + neighbor.cor$origin.row
  neighbor.cor$neighbor.col <- neighbor.cor$col + neighbor.cor$origin.col
  
  #match each possible spot neighbor with existing spots
  spots.neighbors <- merge(neighbor.cor, positions, by.x = c("neighbor.row", "neighbor.col"), by.y = c("row", "col"), suffixes = c(".spots", ".neighbors"))
  spots.neighbors <- lapply(seq_len(n), function(x) as.vector(spots.neighbors[which(spots.neighbors$spot.index.spots == x), ]$spot.index.neighbors))
  
  
  isolate <- which(purrr::reduce((purrr::map(spots.neighbors, length)), c) < 1)
  
  return(isolate)
}

GiottoVisium <- function(sce, platform, q, d, betas, betas_to_add,  result_path, add_file = FALSE, save_name){
  
  # Ensure the Python environment for Giotto has been installed.
  genv_exists = checkGiottoEnvironment()
  if(!genv_exists){
    # The following command need only be run once to install the Giotto environment.
    installGiottoEnvironment()
  }
  
  # Set up Giotto Environment
  # 1.set working directory
  results_folder = result_path
  
  # 2. set giotto python path
  # set python path to your preferred python version path
  # set python path to conda env/bin/ directory if manually installed 
  # Giotto python dependencies by conda
  # python_path = '/path_to_conda/.conda/envs/giotto/bin/python'
  # set python path to NULL if you want to automatically install 
  # (only the 1st time) and use the giotto miniconda environment
  python_path = NULL
  if(is.null(python_path)) {
    installGiottoEnvironment()
  }
  
  # 3. Create Giotto Instructions
  instrs = createGiottoInstructions(save_dir = results_folder,
                                    save_plot = FALSE,
                                    show_plot = FALSE,
                                    python_path = python_path)
  
  
  # Discover isolated spots in sce
  isolate <- find.neighbors.iso(sce, platform)
  
  # Remove isolated spots before creating giotto object
  if(length(isolate) > 0){
    sce_exprs <- assay(sce, "counts")[, -isolate]
    sce_locs <- data.frame(colData(sce)[-isolate,c("col", "row", "barcode")])
  }else{
    sce_exprs <- assay(sce, "counts")
    sce_locs <- data.frame(colData(sce)[ ,c("col", "row")])
    sce_locs$barcode <- colnames(sce)
  }
  
  #Create giotto object
  sce_gobject <- createGiottoObject(expression = sce_exprs, 
                                    spatial_locs = sce_locs,
                                    instructions = instrs)
  
  # Filtered genes expressed in fewer than 10 spots
  sce_gobject <- filterGiotto(gobject = sce_gobject,
                              expression_threshold = 1,
                              feat_det_in_min_cells = 10,
                              min_det_feats_per_cell = 0,
                              expression_values = c('raw'),
                              verbose = T)
  # Normalize
  sce_gobject <- normalizeGiotto(gobject = sce_gobject, verbose = T)
  
  # Add gene & cell statistics
  sce_gobject <- addStatistics(gobject = sce_gobject)
  
  # Highly variable features
  sce_gobject <- calculateHVF(gobject = sce_gobject)
  
  # Run PCA on expression value
  gene_metadata = fDataDT(sce_gobject)
  featgenes = gene_metadata[hvf == 'yes']$feat_ID
  sce_gobject<- Giotto::runPCA(gobject = sce_gobject,
                               feats_to_use = featgenes, ncp = d)
  
  # Create Delanuay network
  sce_gobject <- createSpatialNetwork(gobject = sce_gobject,
                                      method = 'Delaunay', 
                                      maximum_distance_delaunay = "auto",
                                      name = 'spatial_network')
  
  
  hmrf_folder = paste0(results_folder,'/HMRF')
  if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
  
  # do HMRF
  HMRF_spatial_genes = Giotto::doHMRF(gobject = sce_gobject,
                                      expression_values = 'scaled',
                                      spatial_genes = featgenes, k = q, 
                                      spatial_network_name="spatial_network",
                                      betas = betas,   spatial_dimensions = c("sdimx", "sdimy"),
                                      output_folder = gsub("x", q, paste0(hmrf_folder, '/', 'Spatial_genes_kx_scaled')))
  
  # add HMPF results to giotto object
  sce_gobject = addHMRF(gobject = sce_gobject, HMRFoutput = HMRF_spatial_genes, 
                        k = q, betas_to_add = betas_to_add,
                        hmrf_name = 'HMRF')
  
  # get HMRF results dataframe
  results <- pDataDT(sce_gobject)
  
  #add isolated spots
  if(length(isolate) > 0){
    len <- length(isolate)
    add_isolate <- data.frame("cell_ID" = colData(sce)$barcode[isolate],  
                              "nr_feats"=numeric(len), "perc_feats" = numeric(len), 
                              "total_expr" = numeric(len),"cluster" = sample(c(1 : q),  len))
    
    colnames(add_isolate) <- colnames(results)
    
    results <- rbind(add_isolate, results)
  }
  
  
  if(add_file){
    write.csv(results, file = paste0(results_folder, "/", save_name))
  }else{
    # return results
    return(results)
  }
  
}
