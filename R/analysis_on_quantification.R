#' Format data for analysis
#'
#' Create an useful object for analysis with omics data and informations
#' about the omics (metabolites/buckets...) and study design.
#'
#' @param data a data frame contaning omics data with samples in columns and
#' features of interest in rows (metabolites/buckets...).
#' @param design a data frame describing the colums of \code{data} with at
#' least two columns and the first one corresponding to colnames of \code{data}.
#' Default to NULL (colnames of \code{data} is used for study design).
#' @param feature_info a data frame describing the rows of \code{data} with
#' at least two columns and the first one corresponding to rownames of
#' \code{data}. Default to NULL (rownames of \code{data} is used for feature
#' informations).
#' @param zero.threshold remove features with more than \code{zero.threshold}
#' percent of zeros. Default to \code{100}.
#' @param zero.group variable name of design data frame specifying the group
#' variable used to remove features with more than \code{zero.threshold}
#' percent of zeros. Default to \code{NULL}, no group.
#' @param outliers names of outliers to remove.
#'
#' @return An object of type \code{\link{SummarizedExperiment}} useful for
#' further analysis.
#'
#' @examples
#' # Import quantification results
#' library(ASICSdata)
#' quantif_path <- system.file("extdata", "results_ASICS.txt",
#'                             package = "ASICSdata")
#' quantification <- read.table(quantif_path, header = TRUE, row.names = 1)
#'
#' # Import design
#' design <- read.table(system.file("extdata", "design_diabete_example.txt",
#'                                  package = "ASICSdata"), header = TRUE)
#'
#' # Create object for analysis and remove features with more than 25% of
#' #zeros
#' analysis_obj <- formatForAnalysis(quantification,
#'                                   design = design,
#'                                   zero.threshold = 25,
#'                                   zero.group = "condition")
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
formatForAnalysis <- function(data, design = NULL, feature_info = NULL,
                              zero.threshold = 100, zero.group = NULL,
                              outliers = NULL){

  if ((is.null(design) & !is.null(zero.group)) |
      (!is.null(design) && !is.null(zero.group) &&
       !(zero.group %in% colnames(design)))) {
    stop("zero.group need to be a variable name of design data frame")
  }

  # if design and feature_info are NULL use rownames and colnames
  #else same order than for data
  if (is.null(design)) {
    design <- data.frame(sample.name = colnames(data))
  } else {
    if (!all(colnames(data) %in% design[, 1])){
      stop(paste("First column of 'design' matrix must contain all",
                 "colnames of 'data' matrix."))
    }
    design <- design[design[, 1] %in% colnames(data), ]
    design <- design[match(colnames(data), design[, 1]), ]
  }
  if (is.null(feature_info)) {
    feature_info <- data.frame(feature.name = rownames(data))
  } else {
    if (!all(rownames(data) %in% feature_info[, 1])){
      stop(paste("First column of 'feature_info' matrix must contain all",
                 "rownames of 'data' matrix."))
    }
    feature_info <- feature_info[feature_info[, 1] %in% rownames(data), ]
    feature_info <- feature_info[match(rownames(data), feature_info[, 1]), ]
  }

  raw_data <- SummarizedExperiment(assays = as.matrix(data),
                                   rowData = feature_info,
                                   colData = design)

  ## clean data
  clean_data <- raw_data
  if (zero.threshold != 100) {
    if (is.null(zero.group)) {
      zero_percent <- apply(assay(clean_data, 1), 1,
                            function(x) sum(x == 0)/length(x))
      to_remove <- zero_percent > zero.threshold / 100

      clean_data <- clean_data[!to_remove, ]
    } else {
      zero_percent <- apply(assay(clean_data, 1), 1,
                            function(x) by(x, colData(clean_data)[, zero.group],
                                           function(z) sum(z == 0)/length(z)))
      to_remove <- apply(zero_percent, 2,
                         function(x) all(x > zero.threshold / 100))

      clean_data <- clean_data[!to_remove, ]
    }
  }
  if (!is.null(outliers)) {
    clean_data <- clean_data[, !(colData(clean_data)[, 1] %in% outliers)]
  }

  return(clean_data)
}


#' Principal Component Analysis (PCA) on \code{\link{SummarizedExperiment}}
#' object
#'
#' Perform a PCA with the function of \code{\link{ropls}} package on a
#' \code{\link{SummarizedExperiment}} object obtains with the
#' \code{\link{formatForAnalysis}} function
#'
#' @param analysis_data a \code{\link{SummarizedExperiment}} object obtains with
#' the \code{\link{formatForAnalysis}} function.
#' @param scale.unit if \code{TRUE}, data are centered and scaled to unit
#' variance
#' @param type.data type of data used for the analyses (\emph{e. g.}
#' \code{"quantifications"}, \code{"buckets"}...). Default to
#' \code{"quantifications"}.
#' @param condition the name of the design variable (two level factor)
#' specifying the groups.

#'
#' @return A S4 object of class \linkS4class{AnalysisResults} containing PCA
#' results.
#' @seealso \linkS4class{AnalysisResults}
#'
#' @examples
#' # Import quantification results
#' library(ASICSdata)
#' quantif_path <- system.file("extdata", "results_ASICS.txt",
#'                             package = "ASICSdata")
#' quantification <- read.table(quantif_path, header = TRUE, row.names = 1)
#'
#' # Create object for analysis and remove features with more than 25% of
#' #zeros
#' analysis_obj <- formatForAnalysis(quantification,
#'                                     zero.threshold = 25)
#' res_pca <- pca(analysis_obj)
#'
#' @importFrom ropls opls
#' @importFrom SummarizedExperiment assay colData
#' @importFrom gridExtra grid.arrange
#' @export
pca <- function(analysis_data, scale.unit = TRUE,
                type.data = "quantifications", condition = NULL){

  if (!(is.null(condition)) &&
      !(condition %in% names(colData(analysis_data)))) {
    stop("'condition' need to be a variable name of design data frame")
  }

  if (scale.unit) {
    scale_unit <- "standard"
  } else {
    scale_unit <- "none"
  }

  # PCA
  resPCA <- opls(x = t(assay(analysis_data, 1)),
                 predI = min(10, dim(analysis_data)[2]), scaleC = scale_unit,
                 printL = FALSE, plotL = FALSE)

  resPCA@suppLs$yModelMN <- colData(analysis_data)

  mean_by_group <- data.frame()
  if (!is.null(condition)) {
    mean_by_group <- aggregate(t(assay(analysis_data, 1)),
                               list(colData(analysis_data)[, condition]), mean)
    rownames(mean_by_group) <- mean_by_group[, 1]
    mean_by_group[, 1] <- NULL
    mean_by_group <- as.data.frame(t(mean_by_group))
  }

  resPCA_obj <- new(Class = "AnalysisResults",
                    type.analysis = "PCA",
                    type.data = type.data,
                    dataset = analysis_data,
                    results = resPCA,
                    cv.error = numeric(0),
                    mean.by.group = mean_by_group)

  return(resPCA_obj)
}



#' @importFrom gridExtra grid.arrange
.plotPCA <- function(res.pca, graph, add.label = TRUE,
                     axes = c(1, 2), col.ind = NULL){

  if (!(is.null(col.ind)) &&
      !(col.ind %in% colnames(res.pca@suppLs$yModelMN))) {
    stop("'condition' need to be a variable name of design data frame")
  }

  eigen_value <- data.frame(dim = factor(seq_len(10)),
                            eigen_var = res.pca@pcaVarVn,
                            eigen_perc = res.pca@modelDF$R2X * 100)

  to_plot <- list()
  pos <- 1

  # graph
  if ("eig" %in% graph) {
    to_plot[[pos]] <- .plotEigen(eigen_value)
    pos <- pos + 1
  }

  if ("ind" %in% graph) {
    indiv_coord <- data.frame(res.pca@scoreMN[, axes])
    colnames(indiv_coord) <- c("x", "y")

    conditions <- NULL
    if (!is.null(col.ind)) {
      conditions <- res.pca@suppLs$yModelMN[, col.ind]
    }

    to_plot[[pos]] <- .plotIndiv(indiv_coord, eigen_value, add.label, axes,
                                 conditions)
    pos <- pos + 1
  }

  if ("var" %in% graph) {
    var_coord <- data.frame(res.pca@loadingMN[, axes])
    colnames(var_coord) <- c("x", "y")

    # variable contributions on each axes
    var_contrib <- (var_coord ^ 2 * 100) / colSums(var_coord ^ 2)

    # contributions of variables on selected axes
    var_color <- (var_contrib[, 1] * eigen_value$eigen_var[axes[1]] +
                    var_contrib[, 2] * eigen_value$eigen_var[axes[2]]) /
      (eigen_value$eigen_var[axes[1]] + eigen_value$eigen_var[axes[2]])

    to_plot[[pos]] <- .plotVar(var_coord, var_color, eigen_value, axes)
  }

  if (length(graph) == 1) {
    to_return <- to_plot[[1]]
    return(to_return)
  } else {
    grid.arrange(grobs = to_plot, ncol = length(graph))
  }
}


.plotEigen <- function(eigen_value){
  # label for percentage of explained variances
  text_labels <- paste0(round(eigen_value$eigen_perc, 1), "%")

  # graph
  p_eigen <- ggplot(eigen_value, aes_string("dim", "eigen_perc", group = 1)) +
    geom_bar(stat = "identity", fill = "#377DB8", color = "#0B65B2") +
    geom_text(label = text_labels, vjust = -0.4) +
    labs(title = "PCA - Scree plot", x = "Dimensions",
         y = "Percentage of explained variances") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

  return(p_eigen)
}


.plotVar <- function(var_coord, var_color, eigen_value, axes){

  # add label for the 10 most contributive variables
  coord_lab <- var_coord[var_color >= sort(var_color, decreasing = TRUE)[10], ]
  coord_lab$label <- rownames(coord_lab)

  # graph limits
  limit <- max(abs(var_coord$x), abs(var_coord$y)) + 0.05


  # graph
  p_var <-
    ggplot(var_coord) +
    geom_segment(aes_string(x = 0, y = 0, xend = "x", yend = "y",
                            color = "var_color"),
                 arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
                 alpha = ifelse(rownames(var_coord) %in% rownames(coord_lab),
                                1, 0.3)) +
    geom_text(data = coord_lab, aes_string(x = "x", y = "y",
                                           label = "label"),
              hjust = 0.5, vjust = ifelse(coord_lab$y > 0, -0.5, 1.2),
              size = 3, col = "gray30") +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2)  +
    labs(title = "PCA - Plot of variables",
         x = paste0("Dimension ", axes[1], " (",
                    round(eigen_value$eigen_perc[axes[1]], 1), "%)"),
         y = paste0("Dimension ", axes[2], " (",
                    round(eigen_value$eigen_perc[axes[2]], 1), "%)"),
         color = "Contributions") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(colours = c("#1A9641", "#F7F76C", "#D7191C")) +
    coord_fixed() + scale_x_continuous(limits = c(- limit, limit)) +
    scale_y_continuous(limits = c(- limit, limit))

  return(p_var)
}


.plotIndiv <- function(indiv_coord, eigen_value, add.label = TRUE,
                       axes = c(1, 2), condition = NULL){

  # labels to add
  if (add.label) {
    label <- rownames(indiv_coord)
    # label colors
    col.label <- geom_text(aes(label = label), hjust = 0.5, vjust = -0.5,
                           size = 3, col = "gray30", show.legend = FALSE)
  } else {
    label <- ""
  }


  # if color on individuals
  ellipse <- NULL
  plot_centroids <- NULL

  if (!is.null(condition)) {
    # label colors (same color as points)
    col.label <- geom_text(aes(label = label, color = condition),
                           hjust = 0.5, vjust = -0.5, size = 3,
                           show.legend = FALSE)
    # add ellipses and centroids
    ellipse <-
      stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.1)
    centroids <- aggregate(cbind(x, y)~condition, data = indiv_coord, mean)
    plot_centroids <-
      geom_point(data = centroids, aes_string(x = "x", y = "y",
                                              color = "condition"), size = 4)
  }

  # graph
  p_indiv <-
    ggplot(indiv_coord, aes_string(x = "x", y = "y", color = condition)) +
    geom_point() + col.label + ellipse + plot_centroids +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(title = "PCA - Plot of individuals",
         x = paste0("Dimension ", axes[1], " (",
                    round(eigen_value$eigen_perc[axes[1]], 1), "%)"),
         y = paste0("Dimension ", axes[2], " (",
                    round(eigen_value$eigen_perc[axes[2]], 1), "%)"),
         color = "Conditions", fill = "Conditions") +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") + theme_minimal()  +
    theme(plot.title = element_text(hjust = 0.5))

  return(p_indiv)
}

#' Orthogonal projections to latent structures discriminant analysis (OPLS-DA)
#' on \code{\link{SummarizedExperiment}} object
#'
#' Perform an OPLS-DA with the function of \code{\link{ropls}} package on a
#' \code{\link{SummarizedExperiment}} object obtains with the
#' \code{\link{formatForAnalysis}} function
#'
#' @param analysis_data a \code{\link{SummarizedExperiment}} object obtains with
#' the \code{\link{formatForAnalysis}} function.
#' @param condition the name of the design variable (two level factor)
#' specifying the response to be modelled.
#' @param scale.unit if \code{TRUE}, data are centered and scaled to unit
#' variance
#' @param cross.val number of cross validation folds.
#' @param thres.VIP a number specifying the VIP threshold used to identify
#' influential variables.
#' @param orthoI parameter of the \code{\link{opls}} function (package
#' \code{\link{ropls}}) specifying the number of orthogonal components ; when
#' set to NA, OPLS is performed and the number of orthogonal components is
#' automatically computed by using the cross-validation (with a maximum of 9
#' orthogonal components). Default to \code{NA}.
#' @param type.data type of data used for the analyses (\emph{e. g.}
#' \code{"quantifications"}, \code{"buckets"}...). Default to
#' \code{"quantifications"}.
#'
#' @return A S4 object of class \linkS4class{AnalysisResults} containing OPLS-DA
#' results.
#' @seealso \linkS4class{AnalysisResults}
#'
#' @references Thevenot, E.A., Roux, A., Xu, Y., Ezan, E., Junot, C. 2015.
#' Analysis of the human adult urinary metabolome variations with age, body mass
#' index and gender by implementing a comprehensive workflow for univariate and
#' OPLS statistical analyses. \emph{Journal of Proteome Research}.
#' \strong{14}:3322-3335.
#'
#' @examples
#' # Import quantification results
#' library(ASICSdata)
#' quantif_path <- system.file("extdata", "results_ASICS.txt",
#'                             package = "ASICSdata")
#' quantification <- read.table(quantif_path, header = TRUE, row.names = 1)
#'
#' # Import design
#' design <- read.table(system.file("extdata", "design_diabete_example.txt",
#'                                  package = "ASICSdata"), header = TRUE)
#'
#' # Create object for analysis and remove features with more than 25% of
#' #zeros
#' analysis_obj <- formatForAnalysis(quantification,
#'                                   zero.threshold = 25, design = design)
#' res_oplsda <- oplsda(analysis_obj, "condition", orthoI = 1)
#'
#' @export
#' @importFrom ropls opls getSubsetVi predict getVipVn
#' @importFrom stats aggregate
oplsda <- function(analysis_data, condition, scale.unit = TRUE,
                   cross.val = 1, thres.VIP = 1, orthoI = NA,
                   type.data = "quantifications"){

  if (!(condition %in% names(colData(analysis_data)))) {
    stop("'condition' need to be a variable name of design data frame")
  }

  ## scale.unit
  if (scale.unit) {
    scale_unit <- "standard"
  } else {
    scale_unit <- "none"
  }

  ## opls-da and cross-validation
  if (cross.val > 1) {
    folds <- sample(cut(seq_len(ncol(analysis_data)), breaks = cross.val,
                        labels = FALSE), ncol(analysis_data))
    cv_oplsda <-
      sapply(seq_len(cross.val),
             function(i) opls(t(assay(analysis_data, 1)),
                              colData(analysis_data)[, condition], permI = 0,
                              predI = 1, orthoI = orthoI, scaleC = scale_unit,
                              subset = which(folds != i),
                              plotL = FALSE, printL = FALSE))
  } else {
    cv_oplsda <-
      sapply(1:1,
             function(i) opls(t(assay(analysis_data, 1)),
                              colData(analysis_data)[, condition], permI = 0,
                              predI = 1, orthoI = orthoI, scaleC = scale_unit,
                              subset = "odd",
                              plotL = FALSE, printL = FALSE))
  }

  ## prediction error
  cv_error <- round(mean(sapply(cv_oplsda, .errorPred,
                                analysis_data, condition)), 2)

  #VIP
  VIP_all <- sapply(cv_oplsda, getVipVn)
  mean_VIP <- apply(VIP_all, 1, mean)

  ## Influencial feature
  mean_by_group <- aggregate(t(assay(analysis_data, 1)),
                             list(colData(analysis_data)[, condition]), mean)
  rownames(mean_by_group) <- mean_by_group[, 1]
  mean_by_group[, 1] <- NULL
  mean_by_group <- as.data.frame(t(mean_by_group))
  mean_by_group$VIP <- mean_VIP
  mean_by_group$influential <- mean_VIP >= thres.VIP
  mean_by_group <- mean_by_group[order(mean_by_group$VIP, decreasing = TRUE),]

  resOPLSDA_obj <- new(Class = "AnalysisResults",
                     type.analysis = "OPLS-DA",
                     type.data = type.data,
                     dataset = analysis_data,
                     results = cv_oplsda,
                     cv.error = cv_error,
                     mean.by.group = mean_by_group)


  return(resOPLSDA_obj)

}

.errorPred <- function(x, analysis_data, condition){
  trainVi <- getSubsetVi(x)
  comparison <- table(colData(analysis_data)[, condition][-trainVi],
                      predict(x, t(assay(analysis_data, 1))[-trainVi, ]))
  return(sum(comparison[row(comparison) != col(comparison)]) / sum(comparison))
}

#' @importFrom SummarizedExperiment assay
#' @importFrom zoo na.approx
.plotOPLSDA <- function(res.oplsda, graph, add.label = TRUE, xlim, ylim){

  to_plot <- list()
  pos <- 1

  # eigen values
  eigen_value <-
    data.frame(dim = factor(seq_len(2)), eigen_var = c(1, 1),
               eigen_perc = res.oplsda@results[[1]]@modelDF$R2X[1:2] * 100)


  if ("ind" %in% graph) {
    indiv_coord <-
      data.frame(x = res.oplsda@results[[1]]@scoreMN,
                 y = res.oplsda@results[[1]]@orthoScoreMN,
                 row.names = rownames(res.oplsda@results[[1]]@scoreMN))
    colnames(indiv_coord) <- c("x", "y")



    to_plot[[pos]] <-
      .plotIndiv(indiv_coord, eigen_value, add.label = add.label,
                 axes = c(1, 2),
                 condition = res.oplsda@results[[1]]@suppLs$y[getSubsetVi(
                   res.oplsda@results[[1]])]) +
      coord_cartesian() + ylab("Orthogonal component 1") +
      xlab("Dimension 1") +
      ggtitle("OPLS-DA - Plot of individuals")

    pos <- pos + 1
  }

  if ("var" %in% graph) {
    var_coord <-
      data.frame(x = res.oplsda@results[[1]]@loadingMN,
                 y = res.oplsda@results[[1]]@orthoLoadingMN,
                 row.names = rownames(res.oplsda@results[[1]]@loadingMN))
    colnames(var_coord) <- c("x", "y")

    # contributions of variables on selected axes
    var_color <- getVipVn(res.oplsda@results[[1]])

    to_plot[[pos]] <- .plotVar(var_coord, var_color, eigen_value,
                               axes = c(1, 2)) +
      labs(title = "OPLS-DA - Plot of variables",
           y = "Orthogonal component 1",
           x = "Dimension 1",
           color = "VIP")
  }

  if ("buckets" %in% graph) {
    data_wide <- as.data.frame(assay(res.oplsda@dataset, 1))
    buckets_num <- as.numeric(rownames(data_wide))
    # Compute median spectrum
    spect_med <- apply(t(data_wide), 2, median)

    # Compute mean of VIP
    vip <- res.oplsda@mean.by.group[
      order(as.numeric(rownames(res.oplsda@mean.by.group))),]
    vip$VIP[!vip$influential] <- NA

    # Value (in p.p.m) of the bucket correspond to the center so middle points
    #between all buckets are needed to visualize influential buckets with
    #another color
    data_res <- data.frame(buckets = c(buckets_num,
                                       buckets_num - diff(buckets_num)[1]/2),
                           spectrum = c(spect_med, rep(NA, length(buckets_num))),
                           important = rep(vip$influential, 2),
                           VIP = rep(vip$VIP, 2))


    data_res <- data_res[order(data_res$buckets), ]
    data_res[1, 2] <- 0
    data_res$spectrum <- na.approx(data_res$spectrum, data_res$buckets)

    # Plots
    if (is.null(ylim)) ylim <- c(0, max(data_res$spectrum))
    p_buckets <-
      ggplot(data_res, aes_string("buckets", "spectrum", colour = "VIP")) +
      geom_line(aes(group=1), na.rm = TRUE) + theme_bw() +
      scale_x_reverse(limits = rev(c(xlim[1], xlim[2]))) +
      ylim(c(ylim[1], ylim[2])) +
      labs(x = "Chemical shift (in p.p.m.)", y = "Intensity") +
      scale_colour_gradientn(colours = c("#F7F76C", "#D7191C"))

    to_plot[[pos]] <- p_buckets
  }


  if (length(graph) == 1) {
    to_return <- to_plot[[1]]
    return(to_return)
  } else {
    grid.arrange(grobs = to_plot, ncol = length(graph))
  }
}



#' Kruskal-Wallis rank sum test on \code{\link{SummarizedExperiment}} object
#'
#' Perform a Kruskal-Wallis test on a \code{\link{SummarizedExperiment}} object
#' obtains with the \code{\link{formatForAnalysis}} function
#'
#' @param analysis_data a \code{\link{SummarizedExperiment}} object obtains with
#' the \code{\link{formatForAnalysis}} function.
#' @param condition the name of the design variable (two level factor)
#' specifying the group of each sample.
#' @param correction p-value correction method, see \code{\link{p.adjust}}.
#' Default to \code{"BH"}.
#' @param alpha cutoff value for adjusted p-values. Default to 0.05.
#' @param type.data type of data used for the analyses (\emph{e. g.}
#' \code{"quantifications"}, \code{"buckets"}...). Default to
#' \code{"quantifications"}.
#'
#' @return A S4 object of class \linkS4class{AnalysisResults} containing tests
#' results.
#' @seealso \linkS4class{AnalysisResults}
#'
#' @examples
#' # Import quantification results
#' library(ASICSdata)
#' quantif_path <- system.file("extdata", "results_ASICS.txt",
#'                             package = "ASICSdata")
#' quantification <- read.table(quantif_path, header = TRUE, row.names = 1)
#'
#' # Import design
#' design <- read.table(system.file("extdata", "design_diabete_example.txt",
#'                                  package = "ASICSdata"), header = TRUE)
#'
#' # Create object for analysis and remove features with more than 25% of
#' #zeros
#' analysis_obj <- formatForAnalysis(quantification,
#'                                     zero.threshold = 25, design = design)
#' res_tests <- kruskalWallis(analysis_obj, "condition")
#'
#' @importFrom stats p.adjust kruskal.test
#' @importFrom SummarizedExperiment colData colData<-
#' @export
kruskalWallis <- function(analysis_data, condition,
                          correction = c("BH", "bonferroni", "BY", "fdr",
                                         "holm", "hochberg", "hommel"),
                          alpha = 0.05, type.data = "quantifications"){

  if (!(condition %in% names(colData(analysis_data)))) {
    stop("'condition' need to be a variable name of design data frame")
  }

  correction <- try(match.arg(correction), silent = TRUE)
  if(is(correction, "try-error")) {
    stop(paste("'correction' should be one of 'BH', 'bonferroni', 'BY', 'fdr',",
               "'holm', 'hochberg', 'hommel'"))
  }

  kruskal_pval <-
    apply(t(assay(analysis_data, 1)), 2,
          function(x)
            kruskal.test(x~colData(analysis_data)[, condition])$p.value)

  p_adj <- p.adjust(kruskal_pval, method = correction)

  mean_by_group <- aggregate(t(assay(analysis_data, 1)),
                             list(colData(analysis_data)[, condition]), mean)
  rownames(mean_by_group) <- mean_by_group[, 1]
  mean_by_group[, 1] <- NULL
  mean_by_group <- as.data.frame(t(mean_by_group))
  mean_by_group$significant <- p_adj < alpha

  res_tests <- data.frame(Feature = rownames(assay(analysis_data, 1)),
                          "Adjusted p-value" = p_adj)
  res_tests <- res_tests[order(res_tests$Adjusted.p.value,
                               decreasing = FALSE), ]
  rownames(res_tests) <- NULL

  # Add condition to analysis_data with the variable name: conditionTest
  colData(analysis_data)$conditionTest <- colData(analysis_data)[, condition]

  resTest_obj <- new(Class = "AnalysisResults",
                    type.analysis = "Kruskal-Wallis tests",
                    type.data = type.data,
                    dataset = analysis_data,
                    results = res_tests,
                    cv.error = numeric(0),
                    mean.by.group = mean_by_group)

  return(resTest_obj)
}

#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assays
#' @importFrom grDevices boxplot.stats
.plotTests <- function(res.tests, graph, xlim, ylim){
  data_wide <- as.data.frame(assay(res.tests@dataset, 1))

  to_plot <- list()
  pos <- 1

  # graph
  if ("boxplot" %in% graph) {
    #Keep only significant features
    sig <- res.tests@mean.by.group$significant
    data_wide <- data_wide[rownames(data_wide) %in%
                             rownames(res.tests@mean.by.group)[sig], ]

    #Reshape the data frame
    data_long <- reshape(data_wide,
                         idvar = "feature", ids = row.names(data_wide),
                         times = names(data_wide), timevar = "sample",
                         varying = list(names(data_wide)), direction = "long",
                         v.names = "quantif")
    data_long$group <- rep(colData(res.tests@dataset)[, "conditionTest"],
                           each = nrow(data_wide))

    #Boxplot stats
    stat <- tapply(data_long$quantif, list(data_long$feature, data_long$group),
                   function(x) boxplot.stats(x))
    stats <- unlist(tapply(data_long$quantif, list(data_long$feature,
                                                   data_long$group),
                           function(x) boxplot.stats(x)$stats))

    #New labels with p-values
    padj <- res.tests@results[order(res.tests@results$Feature), ]
    padj <- padj[padj$Feature %in% unlist(dimnames(stat)[1]), ]
    padj <- ifelse(padj$Adjusted.p.value < 0.0001, "< 0.0001",
                   round(padj$Adjusted.p.value, 4))
    feature_with_padj <- paste0(unlist(dimnames(stat)[1]), " (adj. p-value: ",
                              padj, ")")

    df <- data.frame(
      feature = rep(rep(feature_with_padj, each = 5),
                    length(unlist(dimnames(stat)[2]))),
      group = rep(unlist(dimnames(stat)[2]),
                  each = length(unlist(dimnames(stat)[1])) * 5),
      quantif = unlist(stats))

    p_boxplot <-
      ggplot(df, aes_string(x = "group", y = "quantif", fill = "group")) +
      geom_boxplot(position = "dodge", show.legend = FALSE) +
      facet_wrap(~feature, scales = "free") +
      labs(x = "Groups", y = "Relative quantification") + theme_bw()

    to_plot[[pos]] <- p_boxplot
    pos <- pos + 1
  }

  if ("buckets" %in% graph) {
    buckets_num <- as.numeric(rownames(data_wide))
    # Compute median spectrum
    spect_med <- apply(t(data_wide), 2, median)

    # Compute direction of change of significant buckets
    difference <- res.tests@mean.by.group[, 2] - res.tests@mean.by.group[, 1]
    dir_change <- ifelse(res.tests@mean.by.group$significant,
                         ifelse(difference > 0,
                                colnames(res.tests@mean.by.group)[2],
                                colnames(res.tests@mean.by.group)[1]),
                         NA)

    # Value (in p.p.m) of the bucket correspond to the center so middle points
    #between all buckets are needed to visualize influential buckets with
    #another color
    data_res <- data.frame(buckets = c(buckets_num,
                                       buckets_num - diff(buckets_num)[1]/2),
                           spectrum = c(spect_med, rep(NA, length(buckets_num))),
                           dir_change = as.factor(rep(dir_change, 2)))


    data_res <- data_res[order(data_res$buckets), ]
    data_res[1, 2] <- 0
    data_res$spectrum <- na.approx(data_res$spectrum, data_res$buckets)

    # Plots
    if (is.null(ylim)) ylim <- c(0, max(data_res$spectrum))
    p_buckets <-
      ggplot(data_res, aes_string("buckets", "spectrum",
                                  colour = "dir_change")) +
      geom_line(aes(group=1), na.rm = TRUE, size = 0.7) + theme_bw() +
      scale_x_reverse(limits = rev(c(xlim[1], xlim[2]))) +
      ylim(c(ylim[1], ylim[2])) +
      labs(x = "Chemical shift (in p.p.m.)", y = "Intensity") +
      scale_colour_manual(name = "Significatively higher in:",
                          values = c("blue", "red"), na.translate = TRUE,
                          na.value = "grey70")

    to_plot[[pos]] <- p_buckets
  }


  if (length(graph) == 1) {
    to_return <- to_plot[[1]]
    return(to_return)
  } else {
    grid.arrange(grobs = to_plot, ncol = length(graph))
  }
}

