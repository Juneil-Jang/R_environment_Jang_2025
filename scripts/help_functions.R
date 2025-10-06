#' Inverted versions of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' Utility function for NULL coalescing
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata") {
  if (!is.null(col_name)) {
    if (col_name %!in% df_colnames) {
      stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir.path) {
  if (!dir.exists(dir.path)) {dir.create(dir.path)}
}

#' Wrapper for missing packages
#'
#' @noRd
check_package <- function(package, repo = "CRAN", git_repo = "") {
  
  if (repo == "CRAN") {
    install_function <- "install.packages('"
  } else if (repo == "github") {
    install_function <- paste0("devtools::install_github('", git_repo, "/")
  } else if (repo == "Bioc") {
    install_function <- "BiocManager::install('"
  }
  
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      paste0("Package ", package, " is not installed.\n",
             "Please run: ", install_function, package, "')"))
  }
  requireNamespace(package, quietly = TRUE)
}

missing_package <- function(...) {
  check_package(...)
}

#'rename the column name with metadata
#'
#' @noRd
rename_columns <- function(df, metadata) {
  # Loop through each column name in the dataframe
  for (i in 1:dim(metadata)[1]) {
    col = metadata$Channel.name[i]
    idx = grep(col, colnames(df))
    col
    idx
    matching_channel <- metadata$Channel.name[grepl( col, metadata$Channel.name)]
    matching_channel
    if (length(idx) > 0) {
      # Get the corresponding marker
      marker <- metadata$markers[metadata$Channel.name == matching_channel]
      # Rename the column with the marker
      colnames(df)[idx] <- marker
    } else {stop(paste(sep = " ", i,idx,col,"\n", "colnames of data/ and ORIGINAL MARKERS.csv are not matched"))}
  }
  return(df)
}

read.cytofFiles <- function (file.loc = getwd(), file.type = ".csv", nrows = NULL, 
                             do.embed.file.names = TRUE, header = TRUE) 
{
  if (!is.element("Spectre", installed.packages()[, 1])) 
    stop("Spectre is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) 
    stop("data.table is required but not installed")
  require(Spectre)
  require(data.table)
  orig_wd <- getwd()
  if (!dir.exists(paste(orig_wd, file.loc, sep = "/")) & !dir.exists(file.loc)) {
    warning("We were not able to find the directory specified by file.loc. Are you sure that location exists?")
  }
  setwd(file.loc)
  wd <- getwd()
  if (length(list.files(path = wd, pattern = file.type)) == 
      0) {
    warning("We did not find any files in that directory, are you sure this is the right place?")
  }
  data.list = list()
  ncol.check = list()
  colName.check = list()
  nrow.check = list()
  if (file.type == ".csv") {
    file.names <- list.files(path = wd, pattern = file.type)
    for (file in file.names) {
      if (is.null(nrows)) {
        tempdata <- data.table::fread(file, check.names = FALSE, 
                                      header = header)
      }
      if (!is.null(nrows)) {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        tempdata <- data.table::fread(file, check.names = FALSE, 
                                      header = header, nrows = nrows)
      }
      file <- gsub(".csv", "", file)
      data.list[[file]] <- tempdata
    }
    rm(tempdata)
    msg <- "CSV files have been imported into a list"
  }
  if (file.type == ".fcs") {
    if (!is.element("flowCore", installed.packages()[, 1])) 
      stop("flowCore is required but not installed")
    require(flowCore)
    file.names <- list.files(path = wd, pattern = file.type)
    for (file in file.names) {
      if (is.null(nrows)) {
        x <- flowCore::read.FCS(file, transformation = FALSE,truncate_max_range = FALSE)
      }
      if (!is.null(nrows)) {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        x <- flowCore::read.FCS(file, transformation = FALSE, 
                                which.lines = nrows,truncate_max_range = FALSE)
      }
      nms <- vector()
      for (o in c(1:nrow(x@parameters@data))) {
        pr <- x@parameters@data$name[[o]]
        st <- x@parameters@data$desc[[o]]
        if (!is.na(st)) {
          nms <- c(nms, paste0(pr, "_", st))
        }
        else {
          nms <- c(nms, pr)
        }
      }
      tempdata <- exprs(x)
      tempdata <- tempdata[1:nrow(tempdata), 1:ncol(tempdata)]
      tempdata <- as.data.table(tempdata)
      names(tempdata) <- nms
      file <- gsub(".fcs", "", file)
      data.list[[file]] <- tempdata
    }
    rm(tempdata)
    msg <- "FCS files have been imported into a list"
  }
  if (do.embed.file.names == TRUE) {
    all.file.names <- c(names(data.list))
    all.file.names
    all.file.nums <- c(1:(length(data.list)))
    all.file.nums
    for (a in all.file.names) {
      data.list[[a]]$FileName <- a
    }
    for (i in all.file.nums) {
      data.list[[i]]$FileNo <- i
    }
  }
  setwd(orig_wd)
  return(data.list)
  message(msg)
}

make.cytofheatmap <- function (dat, sample.col, plot.cols, annot.cols = NULL, feature.annots = NULL, 
                               annotation_colors = NULL, file.name = paste0("Pheatmap by ", 
                                                                            sample.col, ".png"), plot.title = paste0(sample.col, 
                                                                                                                     " heatmap"), transpose = FALSE, normalise = TRUE, is.fold = FALSE, 
                               fold.range = NULL, dendrograms = "both", cutree_rows = 1, 
                               cutree_cols = 1, row.sep = c(), col.sep = c(), cell.size = 15, 
                               standard.colours = "BuPu", fold.colours = "Spectre", path = NULL) 
{
  if (!is.element("pheatmap", installed.packages()[, 1])) 
    stop("pheatmap is required but not installed")
  if (!is.element("RColorBrewer", installed.packages()[, 1])) 
    stop("RColorBrewer is required but not installed")
  if (!is.element("scales", installed.packages()[, 1])) 
    stop("scales is required but not installed")
  require(pheatmap)
  require(RColorBrewer)
  require(scales)
  
  #install.packages("circlize")
  if (standard.colours == "rev(RdBu)") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "RdBu"))(31))
    colour.palette <- rev(colour.palette)
  }
  if (standard.colours == "BuPu") {
    colour.palette <- (colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                 "BuPu"))(31))
  }
  dat <- as.data.frame(dat)
  heatmap.data <- dat
  rownames(heatmap.data) <- t(dat[sample.col])
  heatmap.data
  if (is.null(annot.cols) == FALSE) {
    annot <- heatmap.data[annot.cols]
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (is.null(annot.cols) == TRUE) {
    annot <- NULL
    heatmap.data <- heatmap.data[plot.cols]
    heatmap.data
  }
  if (transpose == TRUE) {
    heatmap.data.t <- as.data.frame(t(heatmap.data))
    heatmap.data <- heatmap.data.t
  }
  if (normalise == TRUE) {
    if (is.fold == FALSE) {
      row.nam <- row.names(heatmap.data)
      col.nam <- names(heatmap.data)
      norm.fun <- function(x) {
        (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
                                      min(x, na.rm = TRUE))
      }
      heatmap.data.norm <- as.data.frame(lapply(heatmap.data, 
                                                norm.fun))
      names(heatmap.data.norm) <- col.nam
      max(heatmap.data.norm)
      heatmap.data.norm <- as.matrix(heatmap.data.norm)
      heatmap.data <- heatmap.data.norm
      rownames(heatmap.data) <- row.nam
    }
  }
  heatmap.data <- as.matrix(heatmap.data)
  if (dendrograms == "none") {
    row.clustering <- FALSE
    col.clustering <- FALSE
  }
  if (dendrograms != "none") {
    hclustfunc <- function(x) hclust(x, method = "complete")
    distfunc <- function(x) dist(x, method = "euclidean")
    if (dendrograms == "both") {
      row.clustering <- TRUE
      col.clustering <- TRUE
    }
    if (dendrograms == "column") {
      row.clustering <- FALSE
      col.clustering <- TRUE
    }
    if (dendrograms == "row") {
      row.clustering <- TRUE
      col.clustering <- FALSE
    }
  }
  
  if (is.fold == FALSE) {
    map.colour <- colour.palette
    sym.key <- FALSE
    sym.breaks <- FALSE
    heatmap.data
    my.max <- function(x) ifelse(!all(is.na(x)), max(x, 
                                                     na.rm = T), NA)
    my.min <- function(x) ifelse(!all(is.na(x)), min(x, 
                                                     na.rm = T), NA)
    # my.breaks <- c(seq(my.min(heatmap.data), median(heatmap.data),length.out=15),seq(median(heatmap.data),my.max(heatmap.data),length.out=17))
    
    my.breaks <- c(seq(my.min(heatmap.data),quantile(heatmap.data,c(0.5)),length.out=15),
                   seq(quantile(heatmap.data,c(0.55)),my.max(heatmap.data),length=17))
  }
  scale.set <- "none"
  title.text <- plot.title
  if (is.null(path)) {
    flnm <- file.name
  }
  if (!is.null(path)) {
    flnm <- paste0(path, "/", file.name)
  }
  
  
  pheatmap::pheatmap(mat = as.matrix(heatmap.data), main = title.text, 
                     cellwidth = cell.size, cellheight = cell.size, cluster_rows = row.clustering, 
                     cluster_cols = col.clustering, breaks = my.breaks, cutree_rows = cutree_rows, 
                     cutree_cols = cutree_cols, gaps_row = row.sep, gaps_col = col.sep, 
                     annotation_row = annot, annotation_col = feature.annots, 
                     annotation_colors = annotation_colors, color = map.colour, 
                     filename = flnm)
  message(paste0("A pheatmap has been saved to your working directory", 
                 paste0(path, file.name)))
}
