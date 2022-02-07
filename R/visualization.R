#' Create a principal component plot or spectral map for HDX-MS deuteration levels
#'
#' @param dataset The HDXMS dataset in long format with at least following columns: peptide, label_time, state and deuterium_level
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#' @param type The type of plot: principal component (pca) or spectral map (sma). Default value is sma.
#'
#' @return A plot.mpm object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom mpm mpm
#' @importFrom stats reshape
#' @export
#'
exploratoryDataAnalysis <- function(dataset, nbRep=NULL, type="sma"){

  if(is.null(nbRep)){
    stop("No number of replicates specified")
  }
  label_time <- NULL
  tmp <- subset(dataset, label_time!=0)
  tmp2 <- tmp[, c("peptide", "label_time", "state", "deuterium_level")]
  df <- convertLongToWide(tmp2, nbRep)

  if(type=="pca"){
    mdsOut <- mpm(df, center = "column", normal = "column")

  }else if(type=="sma"){
    mdsOut <- mpm(df, logtrans = FALSE, row.weight = "mean", col.weight = "mean")

  }else{
    stop("Unknown type specified")
  }

  return(mdsOut)
}

#' create a protein coverage plot
#'
#' @param dataset The HDXMS dataset from an external software tool with at least following columns: Sequence, Start, and End
#'
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot theme labs geom_rect scale_y_discrete
#' @export
#'
coveragePlot <- function(dataset){

  Start <- End <- NULL
  redDataset <- dataset[!duplicated(dataset[,'Sequence']),]
  nbObs <- dim(redDataset)[1]
  pCoverage <- ggplot2::ggplot() +
    geom_rect(data = redDataset, aes(xmin = Start, xmax = End, ymin = rep(-0.45,nbObs), ymax=rep(0.45,nbObs), fill=TRUE), alpha=0.5)+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      legend.position = "none",
      axis.text.y=element_blank())+
    scale_y_discrete(breaks=NULL) +
    labs(x = "Residue Number", y = "")

  return(pCoverage)
}

#' create a butterfly-plot
#'
#' @param dataset output from calculateDeuteriumLevels
#' @param states  vector indicating which protein states should be compared
#'
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot theme labs geom_line geom_point
#' @export
#'
createButterflyPlot <- function(dataset, states = NULL){

  diff2 <- Deuteration_Difference <- adj.P.Val <- state <- label_time <-
    protein_state_A <- protein_state_B <- NULL
  if (is.null(states))
    stop("No protein states specified")
  if (length(states)!=2)
    stop("Incorrect number of states specified")

  position <- factor(paste(dataset$start, dataset$end,sep="-"), levels=unique(paste(dataset$start, dataset$end, sep="-")))
  newData <- cbind(dataset, position)

  selData <- subset(newData, state==states[1] & label_time!=0)
  selData2 <- subset(newData, state==states[2] & label_time!=0)

  mergedData <- data.frame("position"=selData$position, "label_time"=as.factor(selData$label_time), "protein_state_A"= selData$deuterium_level, "protein_state_B" = selData2$deuterium_level)

  pButterfly <- ggplot2::ggplot(data = mergedData) +
    geom_point(aes(x = position, y = protein_state_A, colour = label_time, group = label_time))+
    geom_point(aes(x = position, y = -protein_state_B, colour = label_time, group = label_time))+
    geom_line(aes(x = position, y = protein_state_A, colour = label_time, group = label_time))+
    geom_line(aes(x = position, y = -protein_state_B, colour = label_time, group = label_time))+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey70", size = 0.3),
      axis.text.x = element_text(angle=90, size=6))+
    labs(y = "Deuterium Level", x = "Peptide") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"))

  return(pButterfly)
}

#' create a deuteration-plot
#'
#' @param dataset output from calculateDeuteriumLevels
#' @param peptides  vector indicating which peptides should be shown
#'
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot theme labs geom_smooth geom_point facet_wrap
#' @export
#'
createDeuterationPlot <- function(dataset, peptides = NULL){

  peptide <- label_time <- deuterium_level <- state <- NULL

  if (is.null(peptides))
    stop("No peptides specified")


  selData <- subset(dataset, peptide%in%peptides)
  selData$peptide <- factor(selData$peptide, levels=peptides)
  pDeuterationPlot <- ggplot2::ggplot(selData,
    aes(x = label_time, y = deuterium_level, group = state, colour = state)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ log(x+0.01), se=FALSE) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey70", size = 0.3),
      axis.text.x = element_text(size=8))+
    labs(y = "Deuterium Level", x = "Time") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"))+
    facet_wrap(~ peptide, ncol = 3)
}


#' create a volcano plot from the output of the moderated test-statistics
#'
#' @param diffDeutResult data.frame from one of the moderated test-statistics
#' @param alpha significance-level. Default is 0.05
#'
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot geom_point theme labs aes element_text element_rect element_line
#' @export
#'
createVolcanoPlot <- function(diffDeutResult, alpha = 0.05) {

  diff2 <- Deuteration_Difference <- adj.P.Val <- NULL
  diffDeutResult$diff2 <- "no"
  diffDeutResult$diff2[diffDeutResult$adj.P.Val <= alpha] <- "yes"

  pVolcano <- ggplot2::ggplot(data = diffDeutResult,
    aes(x = Deuteration_Difference, y = -log10(adj.P.Val), group = diff2)) +
    geom_point(aes(colour = diff2, shape = diff2), size=4) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey70", size = 0.3)
    ) +
    labs(x = "Differences in Deuteration", y = "log10(adjusted p-values)") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"),
      legend.title = element_blank())

  return(pVolcano)

}

#' create a barplot showing the number of differential deuterated peptides
#'
#' @param diffDeutResult output from a moderated test-statistics
#' @param alpha significance-level. Default is 0.05
#' @param col a vector of colors for the bars. Default is NULL
#' @param xlabels a vector of names for the labels of the x-axis. Default is NULL
#' @param groupedByTime logical variable indicating if the bars should be organized by labeling time. Default is TRUE
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot geom_bar geom_text theme labs scale_fill_manual scale_x_discrete element_blank
#' @export
#'
createBarPlot <- function(diffDeutResult, alpha=0.05, col=NULL, xlabels=NULL, groupedByTime=TRUE){

  comparison <- dep <- NULL
  df <- data.frame(comparison = names(diffDeutResult), dep = unlist(lapply(diffDeutResult, function(df, a){length(which(df$adj.P.Val<= a))},alpha)))
  timepoints <- paste0("time_",unique(grep("([0-9]$)",unlist(strsplit(df$comparison, split="_time_")),value=TRUE)))
  if(groupedByTime==TRUE)
    df$comparison <- factor(df$comparison, levels=df$comparison[grep(paste(timepoints,collapse="|"), df$comparison)])

  pBarplot <- ggplot2::ggplot(data=df, aes(x=comparison, y=dep, fill=comparison)) +
        geom_bar(stat="identity", width=0.75) +
        geom_text(aes(label=dep), vjust=1.6, color="black", size=3.5)+
        labs(x = "Comparison", y = "Number of differentially deuterated peptides") +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
          panel.border = element_rect(colour = "black", fill = "NA"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey70", size = 0.3),
          axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90),
          legend.position="none")

  if(!is.null(col))
    pBarplot <- pBarplot + scale_fill_manual(values=col)
  if(!is.null(xlabels)){
    pBarplot <- pBarplot + scale_x_discrete(labels=xlabels)
  }
  return(pBarplot)
}

#' create a Woods plot
#'
#' @param diffDeutResult data.frame from one of the moderated test-statistics
#' @param seqInfo data.frame containing peptide sequence, start position and end position
#' @param alpha significance level
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot theme labs geom_segment scale_y_discrete
#' @export
#'
createWoodsPlot <- function(diffDeutResult, seqInfo, alpha = 0.05){

  start <- end <- NULL

  diff2 <- Deuteration_Difference <- adj.P.Val <- NULL
  mergedData <- merge(diffDeutResult,seqInfo, by.x=0, by.y="peptide")
  mergedData$diff2 <- "not significant"
  mergedData$diff2[mergedData$adj.P.Val <= alpha & mergedData$Deuteration_Difference > 0] <- "pos"
  mergedData$diff2[mergedData$adj.P.Val <= alpha & mergedData$Deuteration_Difference < 0] <- "neg"

  pWoods <- ggplot2::ggplot(data = mergedData,
    aes(x=start, xend=end, y=Deuteration_Difference, yend=Deuteration_Difference)) +
    geom_segment(aes(colour = diff2), size=2)+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank())+
    labs(y = "Differences in Deuteration", x = "Residue Number") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"),
      legend.title = element_blank())

  return(pWoods)
}

#' create a MA-plot showing the difference in deuteration versus the average deuteration-value
#'
#' @param diffDeutResult  data.frame from one of the moderated test-statistics
#' @param alpha significance-level. Default is 0.05
#' @return A ggplot object
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom ggplot2 ggplot geom_point geom_text theme labs
#' @export
#'
createMAPlot <- function(diffDeutResult, alpha=0.05){

  diff2 <- Average_Deuteration <- Deuteration_Difference <- adj.P.Val <- NULL
  diffDeutResult$diff2 <- "not significant"
  diffDeutResult$diff2[diffDeutResult$adj.P.Val <= alpha & diffDeutResult$Deuteration_Difference > 0] <- "pos"
  diffDeutResult$diff2[diffDeutResult$adj.P.Val <= alpha & diffDeutResult$Deuteration_Difference < 0] <- "neg"

  pMAplot <- ggplot2::ggplot(data = diffDeutResult,
    aes(x = Average_Deuteration, y = Deuteration_Difference, col = diff2, shape = diff2)) +
    geom_point(size=4) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey70", size = 0.3)
    ) +
    labs(x = "Average Deuteration-levels", y = "Differences in deuteration") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"),
      legend.title = element_blank())

  return(pMAplot)
}

#' create clustered heatmaps
#'
#' @param dataset output from convertLongToWide
#' @param colors  vectors of colors used in heatmap
#' @return A pheatmap
#'
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom pheatmap pheatmap
#' @export
#'
createHeatmap <- function(dataset, colors = NULL){

  if(is.null(colors)){
    pHeatmap <- pheatmap::pheatmap(dataset)
  }else{
    pHeatmap <- pheatmap::pheatmap(dataset, color = colors)
  }
  return(pHeatmap)
}
