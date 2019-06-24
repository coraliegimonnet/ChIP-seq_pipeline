#!/usr/bin/env/Rscript

# packages loading

packages <- c("optparse",
				"GenomicRanges",
				"magrittr",
				"dplyr",
				"tidyr",
				"stringr",
				"ggplot2")

for (package in packages)
{
	if (!(package %in% rownames(installed.packages())))
	{
		print(paste("Need to install the package:", package, sep=""))
		stop("Package missing\n", call. = FALSE)
	}
	do.call('library', list(package))
}

# parameters list

option_list = list (
	make_option(c("-i","--input"), type = "character", default = NULL, metavar = "character",
		help = "Bed file for input"),
	make_option(c("--sharp"), type = "logical", default = FALSE, metavar = "character",
		help = "Union between sharp peaks"),
	make_option(c("--broad"), type = "logical", default = FALSE, metavar = "character",
		help = "Union between broad peaks"),
	make_option(c("--experiment"), type = "character", default = "experiment",
		help = "Name of experiment"),
	make_option(c("--mark"), type = "character", default = "H3K27me3",
		help = "Histone marks of interest"))

opt_parser <- OptionParser(usage = "USAGE : %prog [options]",
						option_list = option_list)

opt <- parse_args(opt_parser)

# parameters check

if (is.null(opt$input))
{
	print_help(opt_parser)
}

# Functions
if (isTRUE(opt$sharp)) {
	qw <- function(words, simplify=T){
	  l.words <- strsplit(words, "\\s+")
	  if (length(l.words) & simplify) return(l.words[[1]])
	  return(l.words)
	}

	write_saf <- function(g, file, geneid.prefix="peak_loc_", append=F, write=T){
	  if (!write) return()
	  col.names <- T
	  if (append) col.names <- F
	  g %>% as.data.frame %>% mutate(GeneID=sprintf("%s%i", geneid.prefix, 1:n()), Strand=".") %>%
	    dplyr::rename(Chr=seqnames, Start=start, End=end) %>%
	    dplyr::select(c(GeneID, Chr, Start, End, Strand)) %>%
	    write.table(., fsaf, sep="\t", quote=F, row.names=F, append=append, col.names=col.names)
	}

	study_hits <- function(hits, g1, g2, label1='q', label2='s', label.equal='='){
  #hits <- findOverlaps(g1, g2)
	  hits.ranges <- ranges(hits, ranges(g1), ranges(g2))
	  hits.df <- hits %>% as.data.frame %>%
	    mutate(seqnames=as.vector(seqnames(g1[queryHits(hits)])),
           start=start(hits.ranges),
           end=end(hits.ranges),
           width=width(hits.ranges),                
           qstart=start(g1[queryHits(hits)]), qend=end(g1[queryHits(hits)]),
           q.p.over=width/((qend-qstart)+1),
           sstart=start(g2[subjectHits(hits)]), send=end(g2[subjectHits(hits)]),
           s.p.over=width/((send-sstart)+1)) %>%
    #l pour less, m pour more
	    mutate(e5=ifelse(qstart==sstart, label.equal, ifelse(qstart==start, label1, label2)),
           e3=ifelse(qend==send, label.equal, ifelse(qend==end, label1, label2)),
           pair.type=sprintf("%s%s", e5, e3))
	  return(hits.df)
	}

	#data
	fpeaks_by_cond_peaks <- opt$input
	#read data
	pepr.colnames <- qw("seqnames start0 end name")
	col_names <- c(pepr.colnames, "exp")
	peaks <- read.table(fpeaks_by_cond_peaks, sep="\t", col.names=col_names)
	peak_type <- "s"
	gpeaks <- with(peaks, GRanges(seqnames=seqnames, IRanges(start=start0+1, end=end), strand="*", 
                              dplyr::select(peaks, -c(seqnames, start0, end)) %>% 
                                mutate(meth=opt$mark,
                                      # jour='J35',  # str_extract(exp, "J1|J35"),
                                       phe=opt$experiment,
                                       peak.type=peak_type, #str_extract(exp, "(b|s)$"),
                                       ipeak=1:n())                            
	))

	gpeaks.rd <- reduce(gpeaks)
	write.really <- T
	# peaks per condition
	gpeaks %>% as.data.frame %>% 
	  group_by(exp) %>% 
	  summarize(n=n(), swidth=sprintf("%.3f", sum(width)/1e6)) %>% 
	  as.data.frame %>% rename("n pics"=n, "taille en Mb"=swidth)

	# peaks fusion
	df.by.p <- data.frame()
	for(m in levels(as.factor(gpeaks$meth))){
	    for (b in levels(as.factor(gpeaks$peak.type))){
	      g <- with(gpeaks, gpeaks[meth==m & peak.type==b])
	      rd <-  reduce(g)
	      if (length(rd)){
	        fbase <- sprintf(opt$experiment,"_%s-%s-union", m, b)
	        fsaf <- sprintf("%s.saf", fbase)
	        fcount <- sprintf("%s.count", fbase)
	        #cat(sprintf("fsaf=%s\n", fsaf))
	        df.by.p <- rbind(df.by.p, 
	                         data.frame(meth=m, peak.type=b, 
	                                    n=length(g),
	                                    width=sprintf("%.2f", sum(width(g))/1e6),
	                                    n.rd=length(rd),
	                                    width.rd=sprintf("%.2f", sum(width(rd))/1e6),
	                                    ft=sprintf("link:%s[]", fsaf),
	                                    count=sprintf("link:%s[]", fcount)
	                         ))
	        write_saf(rd, fsaf, geneid.prefix=sprintf("%s_peak-%s", m, b), write=write.really)
	      }
	    }
	}
	df.by.p  <- df.by.p %>% rename("n pics"=n, "taille en Mb"=width, "n pics rd"=n.rd, "taille rd en Mb"=width.rd)
	#ascii(df.by.p, include.rownames = F, width="90", header=T, grid='none', frame='topbot', format=qw("s s s d f f s s"))#, align='center')
}

if (isTRUE(opt$broad)) {
	qw <- function(words, simplify=T){
	  l.words <- strsplit(words, "\\s+")
	  if (length(l.words) & simplify) return(l.words[[1]])
	  return(l.words)
	}

	write_saf <- function(g, file, geneid.prefix="peak_loc_", append=F, write=T){
	  if (!write) return()
	  col.names <- T
	  if (append) col.names <- F
	  g %>% as.data.frame %>% mutate(GeneID=sprintf("%s%i", geneid.prefix, 1:n()), Strand=".") %>%
	    dplyr::rename(Chr=seqnames, Start=start, End=end) %>%
	    dplyr::select(c(GeneID, Chr, Start, End, Strand)) %>%
	    write.table(., fsaf, sep="\t", quote=F, row.names=F, append=append, col.names=col.names)
	}

	study_hits <- function(hits, g1, g2, label1='q', label2='s', label.equal='='){
	  #hits <- findOverlaps(g1, g2)
	  hits.ranges <- ranges(hits, ranges(g1), ranges(g2))
	  hits.df <- hits %>% as.data.frame %>%
	    mutate(seqnames=as.vector(seqnames(g1[queryHits(hits)])),
	           start=start(hits.ranges),
	           end=end(hits.ranges),
	           width=width(hits.ranges),                
	           qstart=start(g1[queryHits(hits)]), qend=end(g1[queryHits(hits)]),
	           q.p.over=width/((qend-qstart)+1),
	           sstart=start(g2[subjectHits(hits)]), send=end(g2[subjectHits(hits)]),
	           s.p.over=width/((send-sstart)+1)) %>%
	    #l pour less, m pour more
	    mutate(e5=ifelse(qstart==sstart, label.equal, ifelse(qstart==start, label1, label2)),
	           e3=ifelse(qend==send, label.equal, ifelse(qend==end, label1, label2)),
	           pair.type=sprintf("%s%s", e5, e3))
	  return(hits.df)
	}
	  #data 
	fpeaks_by_cond_peaks <- opt$input

	  # read data
	pepr.colnames <- qw("seqnames start0 end name")
	col_names <- c(pepr.colnames, "exp")
	peaks <- read.table(fpeaks_by_cond_peaks, sep="\t", col.names=col_names)
	peak_type <- "b"
	gpeaks <- with(peaks, GRanges(seqnames=seqnames, IRanges(start=start0+1, end=end), strand="*", 
	                              dplyr::select(peaks, -c(seqnames, start0, end)) %>% 
	                                mutate(meth=opt$mark,
	                                       #jour='J35',  # str_extract(exp, "J1|J35"),
	                                       phe=opt$experiment,
	                                       peak.type=peak_type, #str_extract(exp, "(b|s)$"),
	                                       ipeak=1:n())                            
	))

	gpeaks.rd <- reduce(gpeaks)
	write.really <- T

	#peak per condition
	gpeaks %>% as.data.frame %>% 
		group_by(exp) %>% 
		summarize(n=n(), swidth=sprintf("%.3f", sum(width)/1e6)) %>% 
		as.data.frame %>% rename("n pics"=n, "taille en Mb"=swidth)

	# peak union
	df.by.p <- data.frame()
	for(m in levels(as.factor(gpeaks$meth))){
	    for (b in levels(as.factor(gpeaks$peak.type))){
	      g <- with(gpeaks, gpeaks[meth==m & peak.type==b])
	      rd <-  reduce(g)
	      if (length(rd)){
	        fbase <- sprintf(opt$experiment,"%s_%s-union", m, b)
	        fsaf <- sprintf("%s.saf", fbase)
	        fcount <- sprintf("%s.count", fbase)
	        #cat(sprintf("fsaf=%s\n", fsaf))
	        df.by.p <- rbind(df.by.p, 
	                         data.frame(meth=m, jour=j, peak.type=b, 
	                                    n=length(g),
	                                    width=sprintf("%.2f", sum(width(g))/1e6),
	                                    n.rd=length(rd),
	                                    width.rd=sprintf("%.2f", sum(width(rd))/1e6),
	                                    ft=sprintf("link:%s[]", fsaf),
	                                    count=sprintf("link:%s[]", fcount)
	                         ))
	        write_saf(rd, fsaf, geneid.prefix=sprintf("%s_%s_peak", m, b), write=write.really)
	      }
	    }
	}
	df.by.p  <- df.by.p %>% rename("n pics"=n, "taille en Mb"=width, "n pics rd"=n.rd, "taille rd en Mb"=width.rd)
}