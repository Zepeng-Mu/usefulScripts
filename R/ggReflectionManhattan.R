#' Creates a reflection manhattan plot with ggplot2
#' 
#' Creates a reflection manhattan plot from association studies (or any data frame with 
#' chromosome, position, and p-value).
#' 
#' @param x A data.frame with columns "BP," "CHR," "P," and optionally, "SNP."
#' @param chr A string denoting the column name for the chromosome. Defaults to 
#'   PLINK's "CHR." Said column must be numeric. If you have X, Y, or MT 
#'   chromosomes, be sure to renumber these 23, 24, 25, etc.
#' @param bp A string denoting the column name for the chromosomal position. 
#'   Defaults to PLINK's "BP." Said column must be numeric.
#' @param p A string denoting the column name for the p-value. Defaults to 
#'   PLINK's "P." Said column must be numeric.
#' @param snp A string denoting the column name for the SNP name (rs number). 
#'   Defaults to PLINK's "SNP." Said column should be a character.
#' @param col A character vector indicating which colors to alternate.
#' @param chrlabs A character vector equal to the number of chromosomes
#'   specifying the chromosome labels (e.g., \code{c(1:22, "X", "Y", "MT")}).
#' @param suggestiveline Where to draw a "suggestive" line. Default 
#'   -log10(1e-5). Set to NULL to disable.
#' @param genomewideline Where to draw a "genome-wide sigificant" line. Default 
#'   -log10(5e-8). Set to NULL to disable.
#' @param highlight A character vector of SNPs in your dataset to highlight. 
#'   These SNPs should all be in your dataset.
#' @param logp If TRUE, the -log10 of the p-value is plotted. It isn't very 
#'   useful to plot raw p-values, but plotting the raw value could be useful for
#'   other genome-wide plots, for example, peak heights, bayes factors, test 
#'   statistics, other "scores," etc.
#' @param ... Arguments passed on to other plot/points functions
#'   
#' @return A ggplot object reflection manhattan plot.
#'   
#' @keywords visualization manhattan
#'   
#' @import utils
#' @import graphics
#' @import stats
#' @import tidyverse
#' @import ggrepel
#' 
#' @examples
#' ggreflectionManhattan(gwasResults)
#'   
#' @export

ggReflectionManhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                                  col = c("gray80", "gray90"), chrlabs = NULL,
                                  suggestiveline = -log10(1e-7), genomewideline = NULL, 
                                  highlight1 = NULL, highlight2 = NULL,
                                  highlight.snp = "lociName", highlight.text = "gene",
                                  text1 = F, text2 = F, logp = TRUE, annotatePval = NULL, ...) {
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) {
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], pos  =  NA,
                    index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  } else {
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],
                    stringsAsFactors = FALSE, pos = NA, index = NA)
  }
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR, length))
  
  nchr <- length(unique(d$CHR))
  if (nchr == 1) { ## For a single chromosome
    d$pos <- d$BP
    xlabel <- paste('Chromosome', unique(d$CHR), 'position')
  } else { ## For multiple chromosomes
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d$pos[d$index == i] <- d$BP[d$index == i]
      } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
        lastbase <- lastbase + max(d$BP[d$index == (i - 1)])
        d$BP[d$index == i] <- d$BP[d$index == i] - min(d$BP[d$index == i]) + 1
        d$pos[d$index == i] <- d$BP[d$index == i] + lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel <- 'Chromosome'
    labs <- unique(d$CHR)
    #Only show odd chromosome numbers to save space
    labs[labs %% 2 == 0] <- ""
  }
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Create a vector of alternatiting colors
  col <- rep_len(col, max(d$index))
  
  g <- ggplot() +
    theme_classic(base_size = 8) +
    xlab("Position") +
    ylab(expression(-log[10](P))) +
    geom_text(aes(x = max(d$pos), y = max(d$logp), label = "eQTL"),
              vjust = "inward", hjust = "outward", size = 3, col = "tomato") +
    geom_text(aes(x = max(d$pos), y = -max(d$logp), label = "sQTL"),
              vjust = "inward", hjust = "outward", size = 3, col = "purple") +
    theme(legend.position = "none")
  
  # Add points to the plot
  if (nchr == 1) {
    g <- g + geom_point(data = d, mapping = aes(x = pos, y = logp),
                        color = col[1], size = 0.3) +
      geom_point(data = d, mapping = aes(x = pos, y = logp),
                 color = col[1], size = 0.3)
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    g <- g + geom_point(data = d, mapping = aes(x = pos, y = logp, color = as.factor(index)),
                        size = 0.3) +
      geom_point(data = d, mapping = aes(x = pos, y = -logp, color = as.factor(index)),
                 size = 0.3) +
      scale_color_manual(values = col) +
      scale_x_continuous(breaks = ticks, labels = labs)
  }
  
  # We need to plot absolute values for y labels
  # Therefore we need to extract y labels first
  # And then get the absolute value
  ylabs <- ggplot_build(g)$layout$panel_params[[1]]$y.major_source
  
  g <- g + scale_y_continuous(breaks = ylabs, labels = abs(ylabs)) +
    geom_hline(yintercept = 0, color = "black", linetype = 1, lwd = 0.3)
  
  # Add suggestive and genomewide lines
  if (!is.null(suggestiveline)) {
    g <- g + geom_hline(yintercept = suggestiveline,
                        col = "blue", linetype = 2, lwd = 0.5) +
      geom_hline(yintercept = -suggestiveline,
                 col = "blue", linetype = 2, lwd = 0.5)
  }
  
  if (!is.null(genomewideline)){
    g <- g + geom_hline(yintercept = genomewideline,
                        col = "red", linetype = 2, lwd = 0.5) +
      geom_hline(yintercept = -genomewideline,
                 col = "red", linetype = 2, lwd = 0.5)
  }
  
  ## Plot gene text first then highlight SNPs.
  ## So that SNPs are not masked by link in ggrepel
  # Highlight snps from a character vector
  if (!is.null(highlight1)) {
    highlight1 <- data.frame(SNP = highlight1[[highlight.snp]],
                             text = highlight1[[highlight.text]],
                             stringsAsFactors = F)
    
    if (any(!(highlight1$SNP %in% d$SNP))) {
      warning("You're trying to highlight SNPs that don't exist in your results.")
    }
    
    d.highlight1 <- d %>% filter(SNP %in% highlight1$SNP)
  }
  
  # Highlight snps in the lower part 
  if (!is.null(highlight2)) {
    highlight2 <- data.frame(SNP = highlight2[[highlight.snp]],
                             text = highlight2[[highlight.text]],
                             stringsAsFactors = F)
    
    if (any(!(highlight2[[highlight.snp]] %in% d$SNP))) {
      warning("You're trying to highlight SNPs that don't exist in your results.")
    }
    d.highlight2 <- d %>% filter(SNP %in% highlight2$SNP)
  }
  
  # Add gene name labels to the plot
  if (text1) {
    d.text1 <- inner_join(highlight1, d.highlight1, by = "SNP")
    g <- g + geom_text_repel(data = d.text1, mapping = aes(x = pos, y = logp, label = text),
                             color = "grey60", size = 4, segment.size = 0.5, nudge_y = 2,
                             nudge_x = 5, min.segment.length = 0.2)
  }
  
  # Add gene name labels to the lower part of the plot
  if (text2) {
    d.text2 <- inner_join(highlight2, d.highlight2, by = "SNP")
    g <- g + geom_text_repel(data = d.text2, mapping = aes(x = pos, y = -logp, label = text),
                             color = "grey60", size = 4, segment.size = 0.5, nudge_y = -2,
                             nudge_x = -5, min.segment.length = 0.2)
  }
  
  # Add tomato and purple to eQTLs and sQTLs, respectively
  g <- g + geom_point(data = d.highlight1, mapping = aes(x = pos, y = logp),
                      color = "tomato", size = 0.5) +
    geom_point(data = d.highlight2, mapping = aes(x = pos, y = -logp),
               color = "purple", size = 0.5)
  
  # Add black dots to those shared between eQTLs and sQTLs
  common.snp <- intersect(d.highlight1$pos, d.highlight2$pos)
  d.highlight1 <- d.highlight1 %>% filter(pos %in% common.snp)
  d.highlight2 <- d.highlight2 %>% filter(pos %in% common.snp)
  
  g <- g + geom_point(data = d.highlight1, mapping = aes(x = pos, y = logp),
                      color = "black", size = 0.5) +
    geom_point(data = d.highlight2, mapping = aes(x = pos, y = -logp),
               color = "black", size = 0.5)
  
  return(g)
}
