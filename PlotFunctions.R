library(ggplot2)
library(ggrepel)
library(grid)


PlotGOTermsPhase <- function(dat, jtitle, cbPalette, to.append=FALSE, phasestr = "tmid", ampstr = "minuslogpval", amp.max=NULL){
  if (is.null(amp.max)){
    amp.max <- ceiling(max(dat[[ampstr]]))
  } else {
    dat[[ampstr]] <- sapply(dat[[ampstr]], function(a) ifelse(a > amp.max, amp.max, a))
  }
  
  if (is.logical(to.append)){
    m <- ggplot(dat,
                aes_string(x = phasestr, y = ampstr, fill = "Term")) + 
      geom_polygon(alpha = 0.3) + 
      coord_polar(theta = "x") + 
      scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
      theme_bw() +
      ggtitle(jtitle) + 
      geom_hline(yintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
      geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
      theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
            panel.border = element_blank(),
            legend.key = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank())
  } else {
    m <- to.append + geom_polygon(mapping = aes_string(y = phasestr, x = ampstr, fill = "Term", label = NULL), data = dat, alpha = 0.3)
      # coord_polar(theta = "y") + 
      # theme_bw() +
      # theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
      #       panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
      #       panel.border = element_blank(),
      #       legend.key = element_blank(),
      #       axis.ticks = element_blank(),
      #       panel.grid  = element_blank())
  }

  if (!missing(cbPalette)){
    m <- m + scale_fill_manual(values = cbPalette)
  }
  return(m)
}

gg_color_hue <- function(n) {
  # http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

is_top_N <- function(x, N){
  top.n <- sort(x, decreasing = TRUE)[N]
  return(x >= top.n)
}

PlotGeneByRhythmicParameters <- function(fits.best, dat.long, jgene, amp.filt = 0.15, jtitle="", facet.rows = 1, jcex=24, pointsize=1, center=TRUE){  
  source('~/projects/tissue-specificity/scripts/functions/DataHandlingFunctions.R')
  tissue <- as.character(unique(dat.long$tissue))
  # plot expression, but collapse into rhythmic parameters
  dat.sub <- subset(dat.long, gene == jgene)
  if (center){
    dat.sub <- dat.sub %>%
      group_by(tissue) %>%
      mutate(exprs = scale(exprs, center = TRUE, scale = FALSE))
  }
  
  fits.sub <- subset(fits.best, gene == jgene)
  params.lst <- fits.sub$param.list[[1]]
  params <- names(params.lst)
  params.amp <- params.lst[params[grepl("amp", params)]]
  params.amp <- params.amp[which(params.amp > amp.filt)]
  # rhyths <- sapply(names(sort(params.lst[params[grepl("amp", params)]], decreasing = TRUE)), function(n) strsplit(n, "\\.")[[1]][[1]], USE.NAMES = FALSE)
  rhyths <- sapply(names(sort(params.amp, decreasing = TRUE)), function(n) strsplit(n, "\\.")[[1]][[1]], USE.NAMES = FALSE)
  # annotate dat.sub with rhythm param
  # rhyths <- sapply(params.amp, function(p) strsplit(p, "\\.")[[1]][[1]], USE.NAMES = FALSE)
  rhyth.hash <- hash(tissue, sapply(tissue, function(tiss){
    # indx <- rhyths[which(rhyths == tiss)]
    # label as rhythmic group with same parameters
    rhyth.group <- MatchGroup(tiss, rhyths, no.match.str="Flat")
    return(rhyth.group)
  }))
  # annotate dat.sub
  # order by amplitude
  dat.sub$grp <- factor(sapply(as.character(dat.sub$tissue), function(tiss) rhyth.hash[[tiss]], USE.NAMES = FALSE), levels = c(rhyths, "Flat"))
  m <- ggplot(dat.sub, aes(x = time, y = exprs, group = tissue, colour = grp)) + facet_wrap(~grp, nrow = facet.rows) + geom_point(size = pointsize) + geom_line() + 
    theme_bw(jcex) + theme(aspect.ratio = 1, legend.position = "none") + xlab("Time (CT)") + ylab("Log2 mRNA Abundance") + ggtitle(jtitle)
  return(m)
}

FilterGenesByTissue <- function(dat.long, genes.lst, tissues){
  # filter genes by tissue.
  # genes.lst, list of genes, with names matching to tissues
  dat.sub <- lapply(tissues, function(tiss){
    gene.lst <- genes[[tiss]]
    return(subset(dat.long, gene %in% gene.lst & tissue == tiss))
  })
  dat.sub <- do.call(rbind, dat.sub)
  return(dat.sub)
}

PlotOverlayTimeSeries <- function(dat.long, genes, tissues, jscale = T, jalpha = 0.05, jtitle = "", jnrow = 1){
  # jnrow: only for more than 1 tissue
  # if more than 1 tissue, genes is a list of vectors with names matching to tissues
  if (length(tissues) > 1){
    # expect list of genes with tissues as names
    dat.sub <- lapply(tissues, function(tiss){
      gene.lst <- genes[[tiss]]
      return(subset(dat.long, gene %in% gene.lst & tissue == tiss))
    })
    dat.sub <- do.call(rbind, dat.sub)
  } else {
    # genes are just a vector
    dat.sub <- subset(dat.long, gene %in% genes & tissue %in% tissues)  
  }
  # scale and center
  dat.sub <- dat.sub %>%
    group_by(gene, experiment, tissue) %>%
    mutate(exprs.scaled = scale(exprs, center = T, scale = jscale))
  
  if (length(tissues) == 1){
    m <- ggplot(subset(dat.sub), aes(x = time, y = exprs.scaled, group = gene)) + geom_line(alpha = jalpha) + facet_wrap(~experiment) + ggtitle(jtitle)
    m <- m + theme_bw(12) + 
      theme(aspect.ratio=1,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
  } else {
    m <- ggplot(subset(dat.sub), aes(x = time, y = exprs.scaled, group = gene)) + geom_line(alpha = jalpha) + facet_wrap(~tissue, nrow = jnrow) + ggtitle(jtitle)
    m <- m + theme_bw(12) + 
      theme(aspect.ratio=1,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
  }
  if (jscale){
    .ylab <- "Exprs (scaled)"
  } else {
    .ylab <- "Exprs (centered)"
  }
  m <- m + ylab(.ylab) + xlab("Time (CT)")
}

PlotMeanExprsOfModel <- function(dat.mean, genes, jmodel, sorted = TRUE, avg.method = "mean", desc=T){
  dat.mean.sub <- subset(dat.mean, gene %in% genes)
  if (sorted){
    # sort by highesst to lowest mean expression
    if (avg.method == "mean"){
      dat.mean.sum <- dat.mean.sub %>%
        group_by(tissue) %>%
        summarise(mean.exprs = mean(exprs.mean))  
    } else if (avg.method == "median"){
      dat.mean.sum <- dat.mean.sub %>%
        group_by(tissue) %>%
        summarise(mean.exprs = median(exprs.mean))  
    }
    dat.mean.sum <- dat.mean.sum[order(dat.mean.sum$mean.exprs, decreasing = desc), ]
    dat.mean.sub$tissue <- factor(as.character(dat.mean.sub$tissue), levels = as.character(dat.mean.sum$tissue))
  }
  m <- ggplot(dat.mean.sub, aes(x = tissue, y = exprs.mean)) + geom_boxplot() + theme_bw(18) + 
    theme(aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("") + ylab("Log2 mRNA Abundance") +
    ggtitle(paste0("Mean exprs of ", length(genes), " genes in ", jmodel, " module"))
  return(m)
}

GetHistCounts <- function(dat, br = 0:24){
  # get histogram counts to plot histogram circular ggplot2
  h <- hist(dat$phase.avg, br=br,plot=FALSE)
  br[which(br == 24)] <- 0  # loop it
  co <- make_circ_coord(br[-1],h$counts)
  return(data.frame(angles = br[-1], heights = h$counts))
}

make_circ_coord = function(t,x,ttot=24)
{
  dt=(t[2]-t[1])*.45
  a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
  h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
  list(angles=a,heights=h)
}
circular_phase24H_histogram = function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5, jtitle = "", cex.main = 1)
{
  # from Jingkui
  library(plotrix)
  library(circular)
  
  #color_DHS = rgb(0.6,0,0.2)
  par(lwd=lwd,cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
  #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
  br=0:24
  h=hist(x, br=br,plot=FALSE)
  co=make_circ_coord(br[-1],h$counts)
  radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2, rp.type='p',poly.col=color_hist, main = jtitle)
}

OrderDecreasing <- function(dat, jfactor, jval){
  # Reorder factors by decreasing value of jval
  # used in conjunction with barplots makes life easier
  dat[[jfactor]] <- factor(dat[[jfactor]], levels = dat[[jfactor]][order(dat[[jval]], decreasing = TRUE)])
  return(dat)
}

PlotFirstNComponents <- function(svd.complex, comps = c(1), period = 24){
  # Plot first N.comps
  jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)
  for (comp in comps){
    # GetEigens from SvdFunctions.R
    eigens <- GetEigens(svd.complex, period=24, comp=comp, adj.mag = TRUE)
    multiplot(eigens$v.plot, eigens$u.plot, layout = jlayout)
  }
}

PlotAmpPhase <- function(dat, amp_str = "amp", phase_str = "phase", lab_str = "gene", textsize = "auto"){
  # Expect amp and phase in data frame column names.
  # label as gene names
  m <- ggplot(data = dat, aes_string(x = amp_str, y = phase_str, label = lab_str)) + 
    geom_point() + 
    coord_polar(theta = "y") +
    ylab("Phase") +
    xlab("Amp") +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  if (textsize == "auto"){
    m <- m + geom_text(aes_string(x = amp_str, y = phase_str, size = amp_str))
  } else {
    m <- m + geom_text(aes_string(x = amp_str, y = phase_str), size = textsize)
  }
  return(m)
}

PlotAmpPhaseTissue <- function(dat, ampsize = 5){
  # Expect amp and phase in data frame column names.
  # label as tissues
  ggplot(data = dat, aes(x = amp, y = phase, label = tissue)) + 
    geom_point() + 
    geom_text(aes(x = amp, y = phase), size = ampsize) + 
    coord_polar(theta = "y") +
    ylab("Phase") +
    xlab("Amp") +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
}

PlotAmpPhaseAllTissues <- function(dat, jtissue){
  # Expect amp and phase in data frame column names.
  # label as gene names
  if (!missing(jtissue)){
    dat <- subset(dat, tissue == jtissue)
  }
    ggplot(data = dat, aes(x = amp, y = phase, label = gene)) + 
      geom_point(size=0.5) + 
      geom_text(aes(x = amp, y = phase, size = amp)) + 
      coord_polar(theta = "y") +
      xlab("Phase") +
      ylab("Amp") +
      facet_wrap(~tissue) + 
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2)) 
}

PlotComplex <- function(complex.matrix, gene.list, labels,
                        axis.min, axis.max, col="HSV",
                        main='Plot title', 
                        rotate=0,
                        add.text.plot=TRUE,
                        jpch=20,
                        threshold=0,
                        verbose=FALSE,
                        jcex=1){
  # Plot genes on complex plane
  # 
  # ARGS:
  # complex matrix Gene expression. Genes are rows, samples are columns.
  #   expect rownames in complex matrix whcih are gene names (and match genelist)
  # gene.list Optionally filter by gene list
  # colors Colors, if "HSV" compute HSV from angles using PhaseToHsv
  # axis.min: x and y axis min
  # axis.max: x and y axis max
  # main: title
  # rotate: rotate counter clockwise by how many radians? This gives start.angle
  # add.text.plot: add an optional text plot: useful if your gene list is huge.
  # jpch: pch plot symbols '.' for dot, 1 for circles, 20 for filled circles.
  # threshold: if many datapoints, show text label only if magnitude is greater than threshold
  # jcex: adjusting size of main, label and axis of plot. Useful for presentations.
  # 
  # requires wordcloud package for textplot function
  
  library(PhaseHSV)  # for colors around the circle
  
  if (add.text.plot){
    library(wordcloud)  # install.packages("wordcloud") 
  }
  
  # check if data is a matrix. If data frame, coerce to matrix.
  # this prevents errors: "Error non-numeric argument to function
  if (is.data.frame(complex.matrix)){
    complex.matrix <- as.matrix(complex.matrix)
  }
  
  if (missing(gene.list)){
    dat <- complex.matrix  
  } else {
    dat <- complex.matrix[gene.list, ]
  }
  if (missing(labels)){
    text.labels <- rownames(dat)  
  } else {
    text.labels <- labels
  }
  if (missing(axis.min) & missing(axis.max)){
    jmax <- max(Mod(complex.matrix))
    axis.min = -jmax
    axis.max = jmax
  }
  
  # rotate 
  rotation <- complex(modulus = 1, argument = rotate)
  dat <- dat * rotation
  
  plot.colors <- hsv(h=PhaseToHsv(Arg(dat), -pi, pi), s=1, v=1)
  
  plot(dat, col=plot.colors, 
       xlim=c(axis.min, axis.max), 
       ylim=c(axis.min, axis.max), 
       pch=jpch,
       main=main,
       xlab="Real",
       ylab="Complex",
       cex.main=jcex,
       cex.axis=jcex,
       cex.lab=jcex)
  abline(v=0)
  abline(h=0)
  
  # too many data points, only show labels for "large" datapoints
  filter.i <- which(Mod(dat) > threshold)
  text.labels[filter.i]
  text(dat[filter.i],
       labels=text.labels[filter.i], 
       pos=3)
  
  if (verbose){
    cat(paste0(text.labels[filter.i], collapse='", "'))
    # cat(paste0(text.labels[filter.i], collapse="\n"))
    cat("\n")
  }
  
  if (add.text.plot){
    if (is.null(dim(dat)) || nrow(dat) == 1) {
      warning("Data is not a matrix with >1 row. Not adding text plot")
    } else {
      textplot(Re(dat), Im(dat), text.labels, main=main,
               xlim=c(axis.min, axis.max),
               ylim=c(axis.min, axis.max),
               xlab="Real",
               ylab="Complex",
               cex=0.6,
               cex.main=jcex,
               cex.axis=jcex,
               cex.lab=jcex)
      abline(v=0)
      abline(h=0)  
    }
  }
}



PlotRnaMicroarrayFit <- function(tissue, gene, coeff.mat, array.exprs, rna.seq.exprs, rna.tissue){
  # Plotting function when fitting array with expression RNA Seq with microarray.
  # 
  # Args:
  # tissue: string representing tissue name in ARRAY e.g. "Liver"
  # gene: gene to plot
  # coeff.mat: contains intercept and slope for each tissue for every gene. Generated from
  # fit_array_with_exprs.rna.seq.vs.array.R
  # 
  # array.exprs: expression of array
  # 
  # rna.seq.exprs: expression in rnaseq
  # 
  # rna.tissue: string representing tissue name in rna.seq .eg. "Liv"
  # 
  # 
  
  # plot for one tissue only
  t.grep <- tissue
  t.grep.rnaseq <- rna.tissue
  
  intercept <- coeff.mat[gene, paste0(tissue, '_intercept')]
  slope <- coeff.mat[gene, paste0(tissue, '_slope')]
  
  x <- as.matrix(array.exprs[gene, which(grepl(t.grep, colnames(array.exprs)))])
  y <- as.matrix(rna.seq.exprs[gene, which(grepl(t.grep.rnaseq, colnames(rna.seq.exprs)))])
  plot(x, y, xlab='Array exprs',
       ylab='RNA exprs',
       main=paste(gene, 'Intercept=', signif(intercept, digits=3), 'Slope=', signif(slope, digits=3)))
  abline(intercept, slope, lty='dotted')
}

PlotComplexCircle <- function(complex.vector,
                              jlabels,
                              filter = 0,
                              period = 24,
                              rotate = 0,  # by hours
                              xlabel = "Magnitude",
                              ylabel = "Phase (hr)",
                              size = 24,
                              textsize = 6){
  library(ggplot2)
  
  theme_set(theme_grey(base_size = size))
  
  if (missing(jlabels)){
    no.label <- TRUE
  } else {
    no.label <- FALSE
  }
  
  jnames <- names(complex.vector)
  if (is.null(jnames)){
    jnames <- rep(NA, length(complex.vector))
  }
  complex.vector <- as.matrix(complex.vector)
  magnitude <- Mod(complex.vector)
  phase <- Arg(complex.vector)
  phase <- phase * (period / (2 * pi)) + rotate
  phase <- phase %% 24
  
  if (missing(jlabels)){
    jlabels <- jnames
  } else if (filter != 0){
    jlabels <- ifelse(magnitude < filter, "", jnames)
  } 
  
  dat <- data.frame(Magnitude = as.numeric(magnitude), Phase = as.numeric(phase), Labels = jlabels)
  jbreaks <- seq(length.out = period / 2, by = 2, to = period)
  if (!no.label){
    ggplot(data = dat, aes(x = Magnitude, y = Phase, label = Labels)) + 
      geom_point() + 
      geom_text(size = textsize) + 
      coord_polar(theta = "y") +
      xlab(xlabel) +
      ylab(ylabel) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
#   } else if (!no.label & filter > 0) {
#     ggplot(data = dat, aes(x = Magnitude, y = Phase, label = Labels)) + 
#       geom_point() + 
#       geom_text(size = textsize) + 
#       coord_polar(theta = "y") +
#       xlab(xlabel) +
#       ylab(ylabel) +
#       scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  }  else if (no.label) {
    ggplot(data = dat, aes(x = Magnitude, y = Phase)) + 
      geom_point() + 
      coord_polar(theta = "y") +
      xlab(xlabel) +
      ylab(ylabel) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  }
}

PlotFitDiagnostics <- function(array.exprs.full,
                               array.exprs.subset,
                               rna.seq.exprs, 
                               clockgenes, 
                               coeff.mat,
                               outpath, 
                               tissue, 
                               tissue.rna.seq,
                               allow.negs=FALSE){
  if (missing(tissue.rna.seq) == TRUE){
    tissue.rna.seq <- tissue
  }
    
  pdf(outpath)
  par(mfrow=c(2, 2))
  for (gene in clockgenes){
    slope <- coeff.mat[gene, paste0(tissue, "_slope")]
    intercept <- coeff.mat[gene, paste0(tissue, "_intercept")]
    
    # get microarray and rnaseq
    x <- array.exprs.full[gene, grepl(tissue, colnames(array.exprs))]
    y.rna.seq <- rna.seq.exprs[gene, grepl(tissue.rna.seq, colnames(rna.seq.exprs))]
    # get indices of samples with both microarray and rnaseq
    x.rna.seq.i <- colnames(array.exprs.subset[gene, grepl(tissue, colnames(array.exprs.subset))])  # subset only ones that have mRNA matched
    # convert to indices
    x.rna.seq.i <- names(x) %in% x.rna.seq.i  # logical True/False
    # get subset x, containing microarray and rnaseq
    x.subset <- x[x.rna.seq.i]
    
    plot.symbols <- sapply(x.rna.seq.i, function(true.false){
      if (true.false == TRUE){
        symbol <- 1
      } else {
        symbol <- 8
      }
      })
    y <- slope * x + intercept
    if (allow.negs == FALSE){
      y[which(y < 1)] <- 1
    }
    plot(2^x, 2^y, main=paste("o=samp w/ RNASeq+array", gene), pch=c(plot.symbols), xlab="Observed microarray (normal)", ylab="Predicted expression (normal)")
    abline(h=0)
    
    plot(x, y, main=paste("log2", gene, tissue), pch=c(plot.symbols), xlab="Observed microarray (log2)", ylab="Predicted expression (log2)")
    abline(h=0)
    
    # Plot raw on normal scale
    
    plot(2^x.subset, 2^y.rna.seq, main=paste("RNASeq vs Array: normal scale"), xlab="Microarray normal scale", ylab="DESeq normalized counts")
    
    # Plot raw on log2 scale
    plot(x.subset, y.rna.seq, 
         main=paste(gene, 'Intercept=', signif(intercept, digits=3), 
                    'Slope=', signif(slope, digits=3)), 
         xlab="Microarray log2 scale", 
         ylab="DESeq normalized counts (log2)") 
    abline(intercept, slope, lty='dotted')
    }
  dev.off() 
}

PlotBeforeAfter <- function(gene, array.before, array.after, rna.seq, y.max=14, convert.log2=FALSE, N.TISSUES=12){
  par(mfrow=c(3,1))
  if (convert.log2 == FALSE){
    array.before.gene <- as.numeric(array.before[gene, ])
    array.after.gene <- as.numeric(array.after[gene, ])
    rna.seq.gene <- as.numeric(rna.seq[gene, ])    
  } else if (convert.log2 == TRUE){
    array.before.gene <- log2(as.numeric(array.before[gene, ]))
    array.after.gene <- log2(as.numeric(array.after[gene, ]) + 1)
    rna.seq.gene <- log2(as.numeric(rna.seq[gene, ]) + 1)
  }

  # Array before adjustment
  plot(array.before.gene, main=paste(gene, 'log2 expression: array before adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # Array after adjustment
  plot(array.after.gene, main=paste(gene, 'log2 exprs: array after adjustment'),
       col=rep(1:N.TISSUES, each=24), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  # RNA Seq
  plot(rna.seq.gene, main=paste(gene, 'log2 exprs: rnaseq'),
       col=rep(1:N.TISSUES, each=8), type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  par(mfrow=c(1,1))
}

PlotAgainstRnaSeq <- function(gene, rna.seq, array.exprs.adjusted, 
                              common.samples, y.max=14){
  rna.seq.full <- matrix(NA, nrow=1, ncol=ncol(array.exprs.adjusted), 
                    dimnames=list(gene, colnames(array.exprs.adjusted)))
  rna.seq.full[gene, common.samples] <- as.matrix(rna.seq[gene, ])
  
  plot(seq(1:length(rna.seq.full)), rna.seq.full, main=paste(gene, 'black=rnaseq, red=array after adjust'),
       col=1, type='b', ylim=c(0, y.max), ylab="log2 exprs", 
       xlab=paste(tissue.names, collapse=" "))
  lines(as.matrix(array.exprs.adjusted[gene, ]), col=2, pch=22, type='o', cex=0.25)
}

GetFullR <- function(gene, rna.seq.exprs, common.samples){
  R.full <- matrix(0, nrow=1, ncol=ncol(array.exprs), 
                   dimnames=list(gene, colnames(array.exprs)))
  R.full[gene, common.samples] <- as.matrix(rna.seq.exprs[gene, common.samples])
  return(R.full)
}

GetUnobsObsSymbol <- function(all.samples, common.samples, unobs=8, obs=1){
  # create plot symbols. 8 = * = unobserved. 1 = o = observed.
  unobs.symbol <- unobs
  obs.symbol <- obs
  plot.symbols <- matrix(unobs.symbol, nrow=1, ncol=ncol(array.exprs),
                         dimnames=list(gene, colnames(array.exprs)))
  plot.symbols[gene, common.samples] <- obs.symbol
  return(plot.symbols)
}

PlotDiagnostics <- function(gene, array.exprs, rna.seq.exprs, 
                            common.samples, slope, int){
  # use full array.exprs, fill missing rna.seq.exprs with 0s and label with *
  # slope and int comes from coeff.mat
  # create R vs M full 288, R = 0 for "missing" values... for plotting
  
  # plot 2 by 1
  par(mfrow=c(2,1))
  
  R.full <- matrix(0, nrow=1, ncol=ncol(array.exprs), 
                   dimnames=list(gene, colnames(array.exprs)))
  R.full[gene, common.samples] <- as.matrix(rna.seq.exprs[gene, common.samples])
  M.full <- as.matrix(array.exprs[gene, ])
  # create M and R for fitting...
  R <- as.matrix(rna.seq.exprs[gene, common.samples])
  M <- as.matrix(array.exprs.subset.common.g[gene, ])
  # create plot symbols. 8 = * = unobserved. 1 = o = observed.
  unobs.symbol <- 8
  obs.symbol <- 1
  unobs.size <- 0.25
  obs.size <- 1
  plot.symbols <- matrix(unobs.symbol, nrow=1, ncol=ncol(array.exprs),
                         dimnames=list(gene, colnames(array.exprs)))
  plot.symbols[gene, common.samples] <- obs.symbol
  plot.cex <- sapply(plot.symbols, function(x){
    if(x == unobs.symbol){
      # 8 is unobserved 
      return(unobs.size)  # make half size
    } else {
      # 1 is observed
      return(obs.size)  # don't change size
    }
  })
  # int <- coeff.mat[gene, "intercept"]
  # slope <- coeff.mat[gene, "slope"]
  plot(M.full, R.full, main=paste(gene, "slope=", int, "int=", slope), 
       xlab="Microarray (log2)", 
       ylab="RNA-Seq (log2)",
       pch=plot.symbols,
       cex=plot.cex)
  abline(int, slope)
  # plot data on normal scale
  m.norm <- 2^M.full
  r.norm = 2^R.full - 1
  # draw its line in normal scale
  # convert log(y) = a * log(x) + b to normal scale
  f.r.norm <- function(slope, int, m) 2^(int) * m ^ slope - 1
  m.norm.predict <- seq(min(m.norm), max(m.norm), 10)
  r.norm.predict <- f.r.norm(slope, int, m.norm.predict)
  y.max <- max(c(r.norm, r.norm.predict))
  
  plot(m.norm, r.norm, main="norm. scale data. o=observed, *=unobserved",
       xlab="Microarray (normal scale)",
       ylab="RNA-seq (DESeq-normalized count",
       pch=plot.symbols,
       cex=plot.cex,
       ylim=c(0, y.max))
  lines(m.norm.predict, r.norm.predict)
  
  par(mfrow=c(1,1))
}

PlotArgsMatrix <- function(complex.mat, colors, main = "Title", jcex = 1){
  if (missing(colors)){
    hues <- seq(from=0, to=1, length.out=100)
    colors <- hsv(h=hues, s=1, v=1)
  }
  # --------- BEGIN: PLOT MATRIX OF PHASE ANGLES ----------------- # 
  # 
  # order by phase angle
  y <- 1:ncol(complex.mat)  # length of 12, genes are mix of these.
  x <- 1:nrow(complex.mat)  # length of top.N. samples are mix of these.
  par(mar=c(6.1, 5.1, 4.1, 2.1))
  image(x, y, Arg(complex.mat),
        col=colors,
        main=main, 
        axes=FALSE, xlab="", ylab="",
        cex.lab = jcex,
        cex.main = jcex,
        cex.axis = jcex)
  axis(1, at=x, labels=FALSE, tick=FALSE)
  axis(2, at=y, labels=FALSE, tick=FALSE)
  # Now draw the textual axis labels
  # gene labels
  if (jcex != 1){
    gene.cex <- 1
  } else {
    gene.cex <- 0.5
  }
  text(x, par("usr")[3] - 1.5,
       labels=rownames(complex.mat),
       srt=90, 
       pos=4,
       offset=0,
       xpd=TRUE, 
       cex=gene.cex) 
  # sample labels
  text(par("usr")[1] - .3, y, 
       labels = colnames(complex.mat), 
       srt=0, 
       pos=2, 
       offset=0,
       xpd=TRUE, 
       cex=1.5) 
  # ---------- END: PLOT MATRIX OF PHASE ANGLES ------------------- # 
}

PlotLoadings <- function(Loadings, title="Plot title", plot.colors, cex = 1) {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  if (missing(plot.colors)){
    plot.colors <- rep(1:12, each=24)
  }
  plot(Loadings, main=title, col=plot.colors, type='o',
       cex.axis = cex,
       cex.main = cex,
       cex.lab = cex)
}

PCbiplot <- function(PC, jtitle="Plot title", x="PC1", y="PC2") {
  # http://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  data$labsize <- data[[x]] ^ 2 + data[[y]]^2
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.3, aes(label=obsnames, size = labsize))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), color="red")
  plot <- plot + ggtitle(jtitle)
  return(plot)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Add plot activities functions here --------------------------------------

PlotMeanActivitiesWithSE.singleintercept <- function(df){
  jgene <- unique(df$gene)
  ggplot(df,
         aes(x = tissue, y = exprs)) + 
    geom_line() + 
    geom_errorbar(aes(ymax = exprs + se, ymin = exprs - se)) +
    ylab("Activity") +
    ggtitle(jgene)
}

ConvertArgToPhase <- function(phase.rads, omega){
  # convert phase in radians to phase in time, using omega.
  # expect phase to be between -pi to pi, convert that
  # to 0 to 24 hours.
  
  # convert range from -pi to pi to 0 to 2pi
  phase.rads[which(phase.rads < 0)] <- phase.rads[which(phase.rads < 0)] + 2 * pi
  phase <- phase.rads / omega
  return(phase)
}

PlotComplex2 <- function(vec.complex, labels, omega = 2 * pi / 24, 
                         title = "My title", xlab = "Amplitude of activity", ylab = "Phase of activity (CT)", 
                         ampscale = 2, constant.amp = FALSE, dot.col = "gray85", jsize = 22, dotsize = 1.5, dotshape = 18, 
                         disable.text=FALSE,
                         add.arrow=FALSE, disable.repel=FALSE, amp.max="auto", axis.line.col = "black", text.col = "black"){
  # Convert complex to amplitude (2 * fourier amplitude) and phase, given omega.
  # then plot in polar coordinates
  # fourier amplitudes are half-amplitudes of the sine-wave
  # http://www.prosig.com/signal-processing/FourierAnalysis.pdf
  # Default: ampscale = 2 because fourier amplitude is more like half-amplitude by default
  # 
  # Add colours
  # http://stackoverflow.com/questions/20808009/remove-extra-space-and-ring-at-the-edge-of-a-polar-plot
  df <- data.frame(amp = Mod(vec.complex) * ampscale,
                   phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                   label = labels)

  if (amp.max == "auto"){
    amp.max <- ceiling(max(df$amp) * 2) / 2
  } else {
    amp.max <- as.numeric(amp.max)
    if (is.na(amp.max)) warning("Amp max must be numeric or 'auto'")
  }
  if (!is.character(dot.col)){
    # create vector of colors by matching label to colour name
    print("Assuming dot.col is a hash table... with keys as labels")
    dot.col <- sapply(labels, function(l) ifelse(!is.null(dot.col[[l]]), dot.col[[l]], "gray"))
  }
  if (!is.numeric(dotshape)){
    print("Assuming dotshape is a hash table... with keys as labels")
    dotshape <- sapply(labels, function(l) dotshape[[l]])
  }
  if (!is.numeric(dotsize)){
    print("Assuming dotsize is a hash table... with keys as labels")
    dotsize <- sapply(labels, function(l){
		      if (l == ""){
		        return(0.1)
		      } else {
			if (!is.null(dotsize[[l]])){
			  return(dotsize[[l]])
			} else {
			  return(0.1)
			}
		      }
    })
  }
  m <- ggplot(data = df, aes(x = amp, y = phase, label = label)) + 
    geom_point(size = dotsize, colour = dot.col, shape = dotshape) +
    coord_polar(theta = "y") + 
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
    scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
    theme_bw(jsize) + 
    geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = axis.line.col, size = 0.2, linetype = "dashed") +
    geom_hline(yintercept = seq(6, 24, by = 6), colour = axis.line.col, size = 0.2, linetype = "solid") +
    theme(panel.grid.minor = element_blank(), 
          # panel.grid.major = element_line(size = 0.5, colour = "grey"), 
          panel.background = element_blank(), 
          # axis.line = element_line(colour = "black"),
          legend.position="bottom",
          panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          line = element_blank())
  
  if (!disable.text){
    # add text
    df.txt <- subset(df, label != "")
    if (constant.amp != FALSE){
      if (!disable.repel){
        m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp, color = text.col)
      } else {
        m <- m + geom_text(data = df.txt, aes(x = max(amp), y = phase, label = label), size = constant.amp, color = text.col)
      }
    } else {
      m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label), color = text.col)
    }
  }
  if (add.arrow){
    m <- m + geom_segment(aes(x=0, xend=amp, yend=phase), color = "gray85")
  }
  return(m)
}

PlotComplexLong <- function(dat.complex, omega, jtitle = "Title", jxlab = "xlab", jylab = "ylab", ampscale = 2){
  # Plot circle with long form dat input
  if (missing(omega)){
    omega <- 2 * pi / 24
  }
  dat.complex.ampphase <- data.frame(amp = Mod(dat.complex$exprs.transformed) * ampscale,
                   phase = ConvertArgToPhase(Arg(dat.complex$exprs.transformed), omega = omega),
                   jlabel = as.character(dat.complex$tissue))
  m <- ggplot(data = dat.complex.ampphase, aes(x = amp, y = phase, label = jlabel)) + 
    geom_point(size = 0.5) +
    coord_polar(theta = "y") +
    xlab(jxlab) +
    ylab(jylab) +  
    geom_text(aes(x = amp, y = phase, size = amp), vjust = 0) +
    ggtitle(jtitle) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(2, 24, 2))
  return(m)
}

PlotEigengene <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  eigengene <- svd.obj$v[, comp]
  if (rotate){
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
    # rotate eigengene by -phase ref
    eigengene <- eigengene * Conj(rotate.factor)
  }
  v.plot <- PlotComplex2(eigengene, labels = rownames(svd.obj$v), omega = omega, title = paste("Right singular value", comp))
  print(v.plot)
}

PlotEigensamp <- function(svd.obj, comp, omega = 2 * pi / 24, rotate=TRUE){
  if (rotate){
    eigengene <- svd.obj$v[, comp]
    # rotate to phase of largest magnitude in sample of eigengene
    phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    rotate.factor <- complex(modulus = 1, argument = phase.reference)
  }
  eigensamp <- svd.obj$u[, comp]
  eigensamp <- eigensamp * rotate.factor
  u.plot <- PlotComplex2(eigensamp, labels = rownames(svd.obj$u), omega = omega, title = paste("Left singular value", comp))   
  print(u.plot)
}


# Heatmaps ----------------------------------------------------------------

PlotHeatmapGeneList <- function(gene.lst.ordered, dat.long, blackend, minval, maxval, jtiss="", jtitle=""){
  jsub <- subset(dat.long, tissue == jtiss & gene %in% gene.lst.ordered & experiment == "array")
  jsub <- jsub %>%
    group_by(gene) %>%
    mutate(exprs.scaled = scale(exprs, center = TRUE, scale = TRUE))
  M <- dcast(jsub, gene ~ time, value.var = "exprs.scaled"); rownames(M) <- M$gene; M$gene <- NULL
  # reorder
  M <- M[gene.lst.ordered, ]
  PlotExprsHeatmap(M, jtitle = jtitle, blackend=blackend, minval = minval, maxval = maxval, jdendro = "none")  
}

PlotRelampHeatmap <- function(M, jtitle = "Plot Title", blackend = 0.15, yellowstart = 0.151, minval = 0, maxval = 1, dist.method = "manhattan", jdendro = "both"){
  library(gplots)
  my.palette <- colorRampPalette(c("black", "yellow"))(n = 300)
  # # (optional) defines the color breaks manually for a "skewed" color transition
  col.breaks = c(seq(minval, blackend, length=150),  # black
                 seq(yellowstart, maxval, length=151))  # yellow
  # par(mar=c(17,8,4,8)+0.1) 
  par(mar=c(5,4,4,2)+0.1)
  heatmap.2(as.matrix(M), 
            col=my.palette, 
            breaks = col.breaks, 
            scale="none", 
            key=T, 
            keysize=1.5, 
            density.info = "density", 
            trace="none", 
            cexCol=1.7, 
            labRow=NA, 
            dendrogram = jdendro,
            main = jtitle,
            margins=c(8,30),
            distfun = function(x) dist(x, method = dist.method))
}

PlotExprsHeatmap <- function(M, jtitle = "Plot Title", blackend = 0.3, minval = 0, maxval = 1, jdendro = "none"){
  library(gplots)
  my.palette <- colorRampPalette(c("black", "yellow"))(n = 300)
  # # (optional) defines the color breaks manually for a "skewed" color transition
  col.breaks = c(seq(minval, blackend, length=150),  # black
                 seq(blackend + 0.01, maxval, length=151))  # yellow
  heatmap.2(as.matrix(M), 
            col=my.palette, 
            breaks = col.breaks, 
            scale="none", 
            key=T, 
            keysize=1.5, 
            density.info = "density", 
            trace="none", 
            cexCol=1.4, 
            labRow=NA, 
            Rowv=FALSE,
            Colv=FALSE,
            dendrogram = jdendro,
            main = jtitle,
            margins=c(5, 10))
}


PlotHeatmapNconds <- function(fits.best.sub, dat.long, filt.tiss, jexperiment = "array", 
                              blueend = -0.5, blackend = 0.5, min.n = -3, max.n = 3, remove.gene.labels=FALSE, jlabRow = NULL, jlabCol = NULL,
                              jmin.col = "blue", jmax.col = "yellow", jscale=FALSE){
  
  # fits.best.sub <- subset(fits.best, model == jmodel)
  # fits.best.sub <- subset(fits.best, gene %in% genes.tw )
  
  genes <- as.character(fits.best.sub$gene)
  
  jmodel <- as.character(fits.best.sub$model[[1]])
  
  dat.sub <- subset(dat.long, gene %in% genes & experiment == jexperiment & !tissue %in% filt.tiss)
  # dat.sub <- dat.long
  
  # center and scale
  dat.sub <- dat.sub %>%
    group_by(gene, tissue) %>%
    mutate(exprs.scaled = scale(exprs, center = TRUE, scale = jscale))
  

  
  
  mat <- dcast(dat.sub, gene ~ tissue + time, value.var = "exprs.scaled")
  rownames(mat) <- mat$gene; mat$gene <- NULL
  head(mat)
  
  # sort by phases
  phases.dic.keys <- as.character(fits.best.sub$gene)
  phases.dic.vals <- sapply(fits.best.sub$param.list, function(p) p[grep("phase", names(p))][[1]])  # take first phase we see: useful if tissue-wide
  #   p[[paste0(ref, ".phase")]])  # specify by ref
  phases.dat <- data.frame(gene = phases.dic.keys, phase = phases.dic.vals)
  phases.dat <- phases.dat[order(phases.dat$phase), ]
  
  # add phases of individual tissues as a check
  
  # load("Robjs/dat.fit.Robj")
  # source("scripts/functions/GrepRikGenes.R")
  # dat.fit$gene <- FixRikGenes(dat.fit$gene)
  
#   dat.fit.sub <- subset(dat.fit, gene %in% genes)
#   
#   t1.dic.keys <- paste(dat.fit.sub$gene, dat.fit.sub$tissue, sep = ",")
#   t1.dic.vals <- dat.fit.sub$phase
#   t1.dic <- hash(t1.dic.keys, t1.dic.vals)
  
  tiss <- strsplit(gsub(pattern = ";", replacement = ",", x = jmodel), ",")[[1]]
#   for (tis in tiss){
#     phases.dat[[tis]] <- sapply(as.character(phases.dat$gene), function(jgene) t1.dic[[paste(jgene, tis, sep = ",")]])
#   }
  
  # sort by first tissue
  # phases.dat <- phases.dat[order(phases.dat[[tiss[[1]]]]), ]
  phases.dat <- phases.dat[order(phases.dat$phase), ]
  
  # sort by phases
  mat <- as.matrix(mat[as.character(phases.dat$gene), ])
  head(mat)
  
  # heatmap
  n.co <- length(unique(dat.sub$tissue))
  time <- rep(unique(subset(dat.sub, tissue == unique(dat.sub$tissue)[[1]])$time), n.co)
  # condi_name <- as.character(unique(dat.sub$tissue))
  condi_name <- unique(sapply(colnames(mat), function(m) strsplit(m, "_")[[1]][[1]]))  # match from mat
  cond_begins <- seq(1, length(time), length(time) / n.co)# represent index start of next condition. 
  cond_mid <- mean(c(1, length(time) / n.co))  # mid distance between first sample in cond1 and last sample in cond1.
  VLINE_X_LOC <<- cond_begins - 0.5  # Subtract 0.5 to get it between two samples.
  TEXT_X_LOC <<- cond_begins + cond_mid  # a vector in middle of sample, can put text labels conveniently.
  TEXT_Y_LOC <<- 0.99 * length(genes)  # number of genes
  C_NAME <<- unique(condi_name)
  # HLINE_Y_LOC <- 
    
    # mat[1, grep(tiss[[1]], colnames(mat))] <- 1
    # mat <- matrix(1, nrow = length(genes), ncol = length(unique(dat.sub$tissue)) * length(unique(dat.sub$time)))
    
    #   min.n <- -3
    #   max.n <- 3
  
  my.palette <- colorRampPalette(c(jmin.col, "black", jmax.col))(n = 300)
  # # (optional) defines the color breaks manually for a "skewed" color transition
  blackstart <- blueend + 0.01;
  redstart <- blackend + 0.01;
  col.breaks <- c(seq(min.n, blueend, length=100),
                  seq(blackstart, blackend, length=101),
                  seq(redstart, max.n, length = 100))
  # col.breaks = c(seq(0, blackend, length=150),  # black
  #                seq(yellowstart, maxval, length=151))  # yellow
  
  # mat[1, grep(tiss[[1]], colnames(mat))] <- 100
#   jgenes <- c("Nr1d1", "Npas2", "Dbp", "Arntl")
#   for (jgene in jgenes){
#     rowi <- which(rownames(mat) == jgene)
#     rowi.end <- rowi + 10
#     mat[rowi:rowi.end, grep(tiss[[1]], colnames(mat))] <- 10 
#   }

  heatmap.2 (mat,
             # dendrogram control
             Rowv = FALSE,
             Colv=FALSE,
             dendrogram = "none",
             symm = FALSE,
             
             # data scaling
             scale = "none",
             na.rm=TRUE,
             
             # image plot
             add.expr = c(abline(v = VLINE_X_LOC, 
                                 col = 'white'),
                          text(x=TEXT_X_LOC, 
                               y=TEXT_Y_LOC, 
                               labels = C_NAME,
                               col = 'white')),
             
             # mapping data to colors
             # breaks,
             # symbreaks=min(x < 0, na.rm=TRUE) || scale!="none",
             
             # colors
             col = my.palette,
             breaks = col.breaks,
             
             # level trace
             trace="none",
             
             # Row/Column Labeling
             margins = c(5, 5),
             # margins = jmargins,
             labRow = jlabRow,
             labCol = jlabCol,
             srtRow = NULL,
             srtCol = NULL,
             adjRow = c(0,NA),
             adjCol = c(NA,0),
             offsetRow = 0.5,
             offsetCol = 0.5,
             
             # color key + density info
             key = TRUE,
             keysize = 1.5,
             density.info="none",
             
             # plot labels
             main = NULL,
             xlab = NULL,
             ylab = NULL,
             
             # plot layout
             lmat = NULL,
             lhei = NULL,
             lwid = NULL,
             
             # row size col size
             # cexRow = nrow(mat),
             # cexCol = 2 * ncol(mat)
  )
  rownames(phases.dat) <- phases.dat$gene; phases.dat$gene <- NULL
  #   phases.dat$range <- apply(phases.dat, 1, function(row){
  #     # adjust so min = 0
  #     diff(range(row))) 
  #   }
  return(list(mat=mat, phases.dat=phases.dat))
}

PlotHeatmapAmpPhasePval <- function(gene.list, fits.relamp, use.alpha, title="", single.tissue=FALSE){
  fits.sub <- subset(fits.relamp, gene %in% gene.list) 
  # sort by phases of jtiss
  fits.sub <- fits.sub %>%
    group_by(tissue) %>%
    arrange(phase)
  # order by sum of amplitude
  fits.sumamp <- fits.sub %>%
    group_by(tissue) %>%
    summarise(amp.sum = sum(amp)) %>%
    arrange(desc(amp.sum))
  
  tissue.order <- as.character(as.character(fits.sumamp$tissue))
  gene.order <- as.character(subset(fits.sub, tissue == jtiss)$gene)
  
  fits.sub$cols <- as.factor(mapply(PhaseAmpPvalToColor, fits.sub$phase, fits.sub$amp, fits.sub$pval, MoreArgs = list(rotate.hr = -8)))
  # add pvalue for alpha
  fits.sub$cols.alpha <- factor(mapply(function(hex, pval) AddAlphaToHexColor(hex, alpha = SaturationCurve(-log10(pval), Vmax = 1, k = 0.25, x0 = 0)),
                                       as.character(fits.sub$cols), as.numeric(fits.sub$pval)))
  
  if (use.alpha){
    fits.sub$cols.i <- seq(nrow(fits.sub))
    colhash <- hash(as.character(fits.sub$cols.i), as.character(fits.sub$cols.alpha))
  } else {
    # fits.sub$cols.i <- as.numeric(fits.sub$cols)
    fits.sub$cols.i <- seq(nrow(fits.sub))
    colhash <- hash(as.character(fits.sub$cols.i), as.character(fits.sub$cols))
  }
  
  fits.mat <- dcast(fits.sub, formula = gene ~ tissue, value.var = "cols.i")
  rownames(fits.mat) <- fits.mat$gene
  fits.mat$gene <- NULL
  # rreoder
  fits.mat <- fits.mat[gene.order, ]
  fits.mat <- fits.mat[, tissue.order]
  fits.mat <- as.matrix(fits.mat)
  if (!single.tissue){
    par(mar=c(5.1,4.1,6.1,2.1))
    image(t(fits.mat), col = sapply(sort(unlist(fits.mat)), FactorToHex, colhash), yaxt = "n", xaxt = "n", axes=FALSE,xlab="",ylab="",srt=0, main=title)
    axis(3, at = seq(0, 1, length.out = ncol(fits.mat)), labels=colnames(fits.mat), srt=0, tick=FALSE)  
    par(mar=c(5.1,4.1,4.1,2.1))
  } else {
    tissue <- colnames(fits.mat)[[1]]
    fits.mat <- as.matrix(fits.mat)[, 1]  # take single tissue
    # par(mar=c(5.1,20.1,4.1,20.1))
    par(mar=c(5.1,10.1,4.1,10.1))
    image(t(fits.mat), col = sapply(sort(unlist(fits.mat)), FactorToHex, colhash), yaxt = "n", xaxt = "n", axes=FALSE,xlab="",ylab="",srt=0, main=title)
    axis(3, at = seq(0, 1, length.out = 1), labels=tissue, srt=0, tick=FALSE)  
    par(mar=c(5.1,4.1,4.1,2.1))
  }
}
