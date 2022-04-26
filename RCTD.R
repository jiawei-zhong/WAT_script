### Load in/preprocess your data, this might vary based on your file type
refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
# counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
# rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames

counts <- scsncomb@assays$RNA@counts



# meta_data <- read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
# cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
# cell_types <- as.factor(cell_types) # convert to factor data type
# nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list


meta_data <- data.frame(barcode=scsncomb_meta$barcode, 
                        cluster=scsncomb_meta$cluster_anno,
                        nUMI=scsncomb$nCount_RNA)
colnames(meta_data) <- c("barcode", "cluster", "nUMI")
meta_data$nUMI <- as.integer(meta_data$nUMI)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)




## Examine reference object (optional)
print(dim(reference@counts)) #observe Digital Gene Expression matrix


table(reference@cell_types) #number of occurences for each cell type


###########################----s49 start----###########################

## Spatial s49
## requires 1. coords 2. counts 3. nUMI

# datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
# counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
# coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
# rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
# rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
# nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI


st_s49_counts <- st_s49@assays$Spatial@counts
st_s49_coords <- data.frame(xcoord=st_s49@meta.data$imagecol,ycoord= (-st_s49@meta.data$imagerow))
rownames(st_s49_coords) <- rownames(st_s49@meta.data)
st_s49_nUMI <- colSums(st_s49_counts)
### Create SpatialRNA object
st_s49_puck <- SpatialRNA(st_s49_coords, st_s49_counts, st_s49_nUMI)

## Examine SpatialRNA object (optional)
print(dim(st_s49_puck@counts)) # observe Digital Gene Expression matrix
hist(log(st_s49_puck@nUMI,2)) # histogram of log_2 nUMI




print(head(st_s49_puck@coords)) # start of coordinate data.frame
barcodes <- colnames(st_s49_puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(st_s49_puck, barcodes, st_s49_puck@nUMI, ylimit = c(0,round(quantile(st_s49_puck@nUMI,0.9))), 
                     title ='plot of nUMI') 



st_s49_myRCTD <- create.RCTD(st_s49_puck, reference, max_cores = 20, MAX_MULTI_TYPES = 8)
st_s49_myRCTD <- run.RCTD(RCTD = st_s49_myRCTD, doublet_mode = 'multi')





results <- st_s49_myRCTD@results
# normalize the cell type proportions to sum to 1.
st_s49_norm_weights = normalize_weights(results$weights)
st_s49_norm_weights <- as.data.frame(st_s49_norm_weights)


###########################----s49 end----###########################








###########################----s50 start----###########################

## Spatial s50
## requires 1. coords 2. counts 3. nUMI

# datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
# counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
# coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
# rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
# rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
# nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI


st_s50_counts <- st_s50@assays$Spatial@counts
st_s50_coords <- data.frame(xcoord=st_s50@meta.data$imagecol,ycoord= (-st_s50@meta.data$imagerow))
rownames(st_s50_coords) <- rownames(st_s50@meta.data)
st_s50_nUMI <- colSums(st_s50_counts)
### Create SpatialRNA object
st_s50_puck <- SpatialRNA(st_s50_coords, st_s50_counts, st_s50_nUMI)

## Examine SpatialRNA object (optional)
print(dim(st_s50_puck@counts)) # observe Digital Gene Expression matrix
hist(log(st_s50_puck@nUMI,2)) # histogram of log_2 nUMI




print(head(st_s50_puck@coords)) # start of coordinate data.frame
barcodes <- colnames(st_s50_puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(st_s50_puck, barcodes, st_s50_puck@nUMI, ylimit = c(0,round(quantile(st_s50_puck@nUMI,0.9))), 
                     title ='plot of nUMI') 



st_s50_myRCTD <- create.RCTD(st_s50_puck, reference, max_cores = 20, MAX_MULTI_TYPES = 8)
st_s50_myRCTD <- run.RCTD(RCTD = st_s50_myRCTD, doublet_mode = 'multi')





results <- st_s50_myRCTD@results
# normalize the cell type proportions to sum to 1.
st_s50_norm_weights = normalize_weights(results$weights)
st_s50_norm_weights <- as.data.frame(st_s50_norm_weights)


###########################----s50 end----###########################










###########################----s51 start----###########################

## Spatial s51
## requires 1. coords 2. counts 3. nUMI

# datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
# counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
# coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
# rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
# rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
# nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI


st_s51_counts <- st_s51@assays$Spatial@counts
st_s51_coords <- data.frame(xcoord=st_s51@meta.data$imagecol,ycoord= (-st_s51@meta.data$imagerow))
rownames(st_s51_coords) <- rownames(st_s51@meta.data)
st_s51_nUMI <- colSums(st_s51_counts)
### Create SpatialRNA object
st_s51_puck <- SpatialRNA(st_s51_coords, st_s51_counts, st_s51_nUMI)

## Examine SpatialRNA object (optional)
print(dim(st_s51_puck@counts)) # observe Digital Gene Expression matrix
hist(log(st_s51_puck@nUMI,2)) # histogram of log_2 nUMI




print(head(st_s51_puck@coords)) # start of coordinate data.frame
barcodes <- colnames(st_s51_puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(st_s51_puck, barcodes, st_s51_puck@nUMI, ylimit = c(0,round(quantile(st_s51_puck@nUMI,0.9))), 
                     title ='plot of nUMI') 



st_s51_myRCTD <- create.RCTD(st_s51_puck, reference, max_cores = 20, MAX_MULTI_TYPES = 8)
st_s51_myRCTD <- run.RCTD(RCTD = st_s51_myRCTD, doublet_mode = "multi")





results <- st_s51_myRCTD@results
# normalize the cell type proportions to sum to 1.
st_s51_norm_weights = normalize_weights(results$weights)
st_s51_norm_weights <- as.data.frame(st_s51_norm_weights)


###########################----s51 end----###########################
