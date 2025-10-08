Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
library(data.table)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type="character", help="Required String input file prefix (stored in RCTD_input)"),
    make_option(c("-o", "--output"), type="character", help="Required String output file prefix"),
    make_option(c("-t", "--threads"), type="numeric", default=8, help="Number of threads"),
    make_option(c("-r", "--ref"), type="character", help="One of the following: human_tonsil, human_thymus, mouse")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$ref)) {
    print_help(opt_parser)
    stop("All of -i, -o, and -r are required")
}

input <- opt$input
output <- opt$output
numCores <- opt$threads
ref <- opt$ref

cat("input prefix:", opt$input, "\n")
cat("output prefix:", opt$output, "\n")
cat("number of threads:", opt$threads, "\n")
cat("ref:", opt$ref, "\n")

if (ref == "human_tonsil") {
  dge <- "/home/kevin/RCTD_ref/tonsil_atlas_ref/raw_CD4_NBC_MBC_subtype_human_tonsil_dge.csv"
  metadata <- "/home/kevin/RCTD_ref/tonsil_atlas_ref/CD4_NBC_MBC_subtype_human_tonsil_metadata.csv"
} else if (ref == "human_thymus") {
  dge <- "/home/kevin/RCTD_ref/thymus_atlas_ref/thymus_2024_dge.csv"
  metadata <- "/home/kevin/RCTD_ref/thymus_atlas_ref/thymus_2024_metadata.csv"
} else if (ref == "mouse") {
  dge <- "/home/kevin/RCTD_ref/tonsil_atlas_ref/mouse_orthologs_human_tonsil_dge.csv"
  metadata <- "/home/kevin/RCTD_ref/tonsil_atlas_ref/mouse_orthologs_human_tonsil_metadata.csv"
} else {
  stop("Invalid input: must be one of 'human_tonsil', 'human_thymus', or 'mouse'")
}

output_dir <- paste0("/home/kevin/RCTD_ref/RCTD_output/",output,"_RCTD/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


datadir <- "." # directory for slideseq and reference data
savedir <- '.'
if(!dir.exists(savedir))
  dir.create(savedir)

cl <- makeCluster(numCores)
registerDoParallel(cl)
### Load in/preprocess your data, this might vary based on your file type

counts <- read.csv(dge,check.names=FALSE) # load in counts matrix
rownames(counts) <- counts[,1]; 
counts[,1] <- NULL # Move first column to rownames
meta_data <- read.csv(metadata,check.names=FALSE) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$cell_type; 
names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; 
names(nUMI) <- meta_data$barcode # create nUMI named list

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
#> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
#> colSums of counts. If this is unintended, please correct this discrepancy. If
#> this is intended, there is no problem.

## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix
#> [1] 384 475
table(reference@cell_types) #number of occurences for each cell type
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 
#> 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25

## Save RDS object (optional)
# saveRDS(reference, file.path(datadir,'SCRef.rds'))


counts <- read.csv(file.path(datadir,"RCTD_input",paste0(input,"_MappedDGEForR.csv")),check.names=FALSE) # load in counts matrix
coords <- read.csv(file.path(datadir,"RCTD_input",paste0(input,"_BeadLocationsForR.csv")), row.names=1)
rownames(counts)<-counts[,1]; 
counts[,1] <- NULL # Move first column to rownames
`rownames<-`(coords, coords$barcodes) 
#coords$barcodes
#rownames(coords)
coords$barcodes <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

## Examine SpatialRNA object (optional)
print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI


print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

myRCTD <- create.RCTD(puck, reference, max_cores = numCores, CELL_MIN_INSTANCE=9)

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

saveRDS(myRCTD, file.path(paste0(output,'_myRCTD.rds')))

intersection_output <- as.data.frame(as.matrix(myRCTD@results[["results_df"]]))
fwrite(x = intersection_output, file = paste0(output_dir,output,'_output.csv'), row.names = TRUE)
weights <- as.data.frame(as.matrix(myRCTD@results$weights))
norm_weights <- normalize_weights(weights)
fwrite(x = norm_weights, file = paste0(output_dir,output,'_norm_weights.csv'), row.names = TRUE)
