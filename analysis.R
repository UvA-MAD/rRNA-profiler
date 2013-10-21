library(ShortRead)
library(ggbio)
library(Rsamtools)
library(rtracklayer)
library(ggplot2)
library(yaml)
library(reshape)

# parse configuration in yaml
conf <- yaml.load_file('config.yaml')
# data source files
## raw reads
fq.path <- conf$fq.path
## rRNA annotations
gff.path <- conf$gff.path
## rRNA reference
rRNA.ref <- conf$rRNA.ref
rRNA5S.ref <- conf$rRNA5S.ref
# data analysis path
data.home <- conf$data.home

## trmimmed split by len 
fq.short.path <- paste(data.home, 'short.fastq', sep='')
fq.long.path <- paste(data.home, 'long.fastq', sep='')
## bam files with alignments
bam.short <- paste(data.home, 'short.bam', sep='')
bam.long <- paste(data.home, 'long.bam', sep='')


# create dir for plots
plot.dir <- paste(data.home, 'plots', sep='')
dir.create(plot.dir)

# trim and split into long and short reads file
if (! file.exists(fq.short.path)) {
    # -- trim reads -- #
    # 3 nucleotides in a window need to be of specified qual
    window.len <- 3 
    # minimal qual
    min.qual <- 20
    min.qual.char <- "5"

    # use trimTailw to trim the reads
    ## read fastq
    fq <- readFastq(fq.path)
    fq.trimmed <- trimTails(fq, k=window.len, a=min.qual.char, successive=T)

    # split reads into 3 len categories:
    # shortest (not used): [1:14] short: [15:25] long: [26:]
    fq.short <- fq.trimmed[width(fq.trimmed) %in% 15:25]
    fq.long <- fq.trimmed[width(fq.trimmed) > 25]
    # write them to the file:
    writeFastq(fq.short, fq.short.path)
    writeFastq(fq.long, fq.long.path)
}

# run alignment
if (! file.exists(bam.short)){
    tmap.command.short <- paste('tmap mapall -a1 -g 3 -n 30',
                          '-f', rRNA.ref,
                          '-r', fq.short.path,
                          '-o 1 -v stage1 map1 map2 map3',
                          '>', bam.short)
    system(tmap.command.short, intern=T)
    tmap.command.long <- paste('tmap mapall -a1 -g 3 -n 30',
                          '-f', rRNA.ref,
                          '-r', fq.long.path,
                          '-o 1 -v stage1 map1 map2 map3',
                          '>', bam.long)
    system(tmap.command.long, intern=T)

    # sort and index bam files
    sortBam(bam.short, strsplit(bam.short, split="\\.")[[1]][1])
    sortBam(bam.long, strsplit(bam.long, split="\\.")[[1]][1])
    indexBam(bam.short)
    indexBam(bam.long)
}

len.categories <- c('short', 'long')

short.bam <- paste(data.home, 'short_low_mm.bam', sep='')
long.bam <- paste(data.home, 'long_low_mm.bam', sep='')

# filter alignments according to criteria
# 0.1 * read width MM
if (! file.exists(short.bam)) {
    ## run python script to filter the read
    system(paste('python ', data.home, 'filter_by_n_mm.py', sep=''))

    ## index resulting file
    indexBam(short.bam)
    indexBam(long.bam)
}

# read in the reads
short.reads <- readGappedAlignments(short.bam, use.name=T)
long.reads <- readGappedAlignments(long.bam, use.name=T)

# combind both categories into one Granges object with
# extra meta column 'length category' 
mcols(short.reads)$length.category <- rep('short',
                                            length(short.reads))
mcols(long.reads)$length.category <- rep('long',
                                           length(long.reads))
all.reads <- c(short.reads, long.reads)

# plotting (for short long and combined)
## consensus with rRNA annotation track 

### plotting reads on consensus with annotation track

### annotation
annots <- import(gff.path, asRangedData=F)
annot.tr <- autoplot(annots[seqnames(annots)=='EMBOSS_001|consensus'],
                        aes(color=group, fill=group))

### short reads
reads.tr <- autoplot(short.reads[seqnames(short.reads) == 'EMBOSS_001|consensus'])

### plot and save 
png(filename=paste(data.home,
                   'plots/short_reads_on_consensus',
                   '.png', sep=''),
    width=1024, height=512)
print(tracks(coverage=reads.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()

### long reads
reads.tr <- autoplot(long.reads[seqnames(long.reads) == 'EMBOSS_001|consensus'])

### plot and save 
png(filename=paste(data.home,
                   'plots/long_reads_on_consensus',
                   '.png', sep=''),
    width=1024, height=512)
print(tracks(coverage=reads.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()

### long and short reads together on consensus
reads.tr <- autoplot(all.reads[seqnames(all.reads) == 'EMBOSS_001|consensus'],
                     aes(color=length.category, fill=length.category))

### plot and save 
png(filename=paste(data.home,
                   'plots/combined_reads_on_consensus',
                   '.png', sep=''),
    width=1024, height=512)
print(tracks(coverage=reads.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()

# plot coverage on consensus
cov.short.reads <- coverage(short.reads)

plot.fname <- paste(data.home, 'plots/short_read_coverage_on_consensus.png', sep='')
png(filename=plot.fname, width=1024, height=512)
cov.tr <- autoplot(cov.short.reads$'EMBOSS_001|consensus', binwidth=1,
               main='short read coverage on consensus')
print(tracks(coverage=cov.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()


cov.long.reads <- coverage(long.reads)

plot.fname <- paste(data.home, 'plots/long_read_coverage_on_consensus.png', sep='')
png(filename=plot.fname, width=1024, height=512)
cov.tr <- autoplot(cov.long.reads$'EMBOSS_001|consensus', binwidth=1,
               main='long read coverage on consensus')
print(tracks(coverage=cov.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()

cov.all.reads <- coverage(all.reads)

plot.fname <- paste(data.home, 'plots/combined_read_coverage_on_consensus.png', sep='')
png(filename=plot.fname, width=1024, height=512)
cov.tr <- autoplot(cov.all.reads$'EMBOSS_001|consensus', binwidth=1,
               main='combined read coverage on consensus')
print(tracks(coverage=cov.tr, annotation=annot.tr, heights=c(5,1)))
dev.off()

### all rRNA species plotted separately
### (short reads, long reads, combined, as reads and as coverage)

nameFromGroup <- function(n) {
    strsplit(n, split='=')[[1]][2]
}

library(plyr)
rRNA.names <- aaply(as(annots$group, 'character'), 1,
                    nameFromGroup) 

slevels <- c(rep('EMBOSS_001|consensus', 6),
                 'gi|15079186:1020-1971',
                 'gi|15079186:2043-3725')

short.read.overlaps <- findOverlaps(annots, short.reads)
long.read.overlaps <- findOverlaps(annots, long.reads)
all.read.overlaps <- findOverlaps(annots, all.reads)

## iterate over all rRNA species in rRNA names
#for testing this should go in a loop 
for (i in c(1:8)){
    rrna <- rRNA.names[[i]]
    slevel <- slevels[i]

    # fetch reads in overlapping with current annotation
    loc.short.reads <- short.reads[subjectHits(short.read.overlaps)[queryHits(short.read.overlaps) == i]]
    loc.long.reads <- long.reads[subjectHits(long.read.overlaps)[queryHits(long.read.overlaps) == i]]
    loc.all.reads <- all.reads[subjectHits(all.read.overlaps)[queryHits(all.read.overlaps) == i]]

    seqlevels(loc.short.reads) <- slevel
    seqlevels(loc.long.reads)  <- slevel
    seqlevels(loc.all.reads) <- slevel

    # plot reads
    plot.fname <- paste(data.home, 'plots/short_reads_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(loc.short.reads,
                   main=paste('short reads on', rrna)))
    dev.off()

    plot.fname <- paste(data.home, 'plots/long_reads_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(loc.long.reads,
                   main=paste('long reads on', rrna)))
    dev.off()

    plot.fname <- paste(data.home, 'plots/combined_reads_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(loc.all.reads,
                   aes(fill=length.category, color=length.category),
                   main=paste('short and long reads on', rrna)))
    dev.off()

    # plot coverage
    cov.short.reads <- coverage(loc.short.reads)
    cov.short.reads <- cov.short.reads[[1]][start(annots[i]):end(annots[i])]
    cov.long.reads <- coverage(loc.long.reads)
    cov.long.reads <- cov.long.reads[[1]][start(annots[i]):end(annots[i])]
    cov.all.reads <- coverage(loc.all.reads)
    cov.all.reads <- cov.all.reads[[1]][start(annots[i]):end(annots[i])]

    plot.fname <- paste(data.home, 'plots/short_read_coverage_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(cov.short.reads, binwidth=1,
                   main=paste('short read coverage on', rrna)))
    dev.off()

    plot.fname <- paste(data.home, 'plots/long_read_coverage_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(cov.long.reads, binwidth=1,
                   main=paste('long read coverage on', rrna)))
    dev.off()

    plot.fname <- paste(data.home, 'plots/combined_read_coverage_on_', rrna, '.png', sep='')
    png(filename=plot.fname, width=1024, height=512)
    print(autoplot(cov.all.reads, binwidth=1,
                   main=paste('all reads coverage on', rrna)))
    dev.off()
} 
# end of the loop

# barplots
short.read.counts <- countOverlaps(annots, short.reads)
long.read.counts <- countOverlaps(annots, long.reads)
all.read.counts <- countOverlaps(annots, all.reads)

# 5S read counts need some more work
rRNA5S.fa <- readDNAStringSet(rRNA5S.ref)
rRNA5S.ids <- names(rRNA5S.fa)
rRNA5S.short.read.counts  <- length(short.reads[seqnames(short.reads) %in% rRNA5S.ids])
rRNA5S.long.read.counts  <- length(long.reads[seqnames(long.reads) %in% rRNA5S.ids])
rRNA5S.all.read.counts  <- length(all.reads[seqnames(all.reads) %in% rRNA5S.ids])
rRNA5S.len <- median(width(rRNA5S.fa))
rRNA5S.read.counts.df <- data.frame(rRNA= rep('5S_rRNA', 3),
                                    count=c(rRNA5S.short.read.counts,
                                            rRNA5S.long.read.counts,
                                            rRNA5S.all.read.counts),
                                    length.category=c('short', 'long', 'all'),
                                    rRNA.length=rRNA5S.len)

# get lengths of rRNA species
rRNA.lengths <- width(annots)

# put the numbers into dataframes
raw.short.read.counts <- data.frame(rRNA=rRNA.names,
                                    count=short.read.counts,
                                    length.category='short',
                                    rRNA.length=rRNA.lengths)
raw.long.read.counts <- data.frame(rRNA=rRNA.names,
                                   count=long.read.counts,
                                   length.category='long',
                                   rRNA.length=rRNA.lengths)
raw.all.read.counts <- data.frame(rRNA=rRNA.names,
                                  count=all.read.counts,
                                  length.category='all',
                                  rRNA.length=rRNA.lengths)


read.counts <- rbind(raw.short.read.counts,
                     raw.long.read.counts,
                     raw.all.read.counts,
                     rRNA5S.read.counts.df)

plot.fname <- paste(data.home, 'plots/raw_counts_per_rRNA.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <-  ggplot(read.counts, aes(x=rRNA, y=count, fill=length.category))
p <- p + ggplot2::geom_bar(position='dodge', stat='identity') 
p <- p + ggtitle('raw read counts')
print(p)
dev.off()

# normalize counts per 1k of rRNA length
read.counts$counts.per.1k <- read.counts$count * 1000 / read.counts$rRNA.length

# plot normalized
plot.fname <- paste(data.home, 'plots/normalized_counts_per_rRNA.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <-  ggplot(read.counts, aes(x=rRNA, y=counts.per.1k, fill=length.category)) 
p <- p + ggplot2::geom_bar(position='dodge', stat='identity') 
p <- p + ggtitle('read counts normalized for rRNA length')
print(p)
dev.off()

# count reads in each 5S rRNA 
rRNA5S.short.reads <- short.reads[seqnames(short.reads) %in% rRNA5S.ids]
rRNA5S.long.reads <- long.reads[seqnames(long.reads) %in% rRNA5S.ids]
rRNA5S.all.reads <- all.reads[seqnames(all.reads) %in% rRNA5S.ids]
rRNA5S.short.split.counts <- data.frame(table(seqnames(rRNA5S.short.reads)))
colnames(rRNA5S.short.split.counts) <- c('seqid','short.read.count')
rRNA5S.long.split.counts <- data.frame(table(seqnames(rRNA5S.long.reads)))
colnames(rRNA5S.long.split.counts) <- c('seqid','long.read.count')
rRNA5S.all.split.counts <- data.frame(table(seqnames(rRNA5S.all.reads)))
colnames(rRNA5S.all.split.counts) <- c('seqid','all.read.count')
rRNA5S.counts <- merge(rRNA5S.short.split.counts,
                          rRNA5S.long.split.counts)
rRNA5S.counts  <- merge(rRNA5S.counts, rRNA5S.all.split.counts)

# make barplots
plot.fname <- paste(data.home, 'plots/short_reads_5S_rRNA_counts.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <- ggplot(rRNA5S.counts, aes(x=seqid, y=short.read.count))
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('short reads count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank())
print(p)
dev.off() 
plot.fname <- paste(data.home, 'plots/long_reads_5S_rRNA_counts.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <- ggplot(rRNA5S.counts, aes(x=seqid, y=long.read.count))
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('long reads count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank())
print(p)
dev.off() 
plot.fname <- paste(data.home, 'plots/all_reads_5S_rRNA_counts.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <- ggplot(rRNA5S.counts, aes(x=seqid, y=all.read.count)) 
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('all reads count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank())
print(p) 
dev.off() 

# make barplots for top 10 5s RNA sequences
short.top.ten.5S <- rRNA5S.counts[order(rRNA5S.counts$short.read.count,
                                     decreasing=T),][1:10,1:2]
long.top.ten.5S <- rRNA5S.counts[order(rRNA5S.counts$long.read.count,
                                     decreasing=T),][1:10,c(1,3)]
all.top.ten.5S <- rRNA5S.counts[order(rRNA5S.counts$all.read.count,
                                     decreasing=T),][1:10,c(1,4)]

# calculate the fraction of reads mapping to paricular seq
total.short.5S <- sum(rRNA5S.counts$short.read.count)
total.long.5S <- sum(rRNA5S.counts$long.read.count)
total.all.5S <- sum(rRNA5S.counts$all.read.count)

short.top.ten.5S$read.fraction <- round(short.top.ten.5S$short.read.count / total.short.5S, digits=2)

plot.fname <- paste(data.home, 'plots/short_top_5S.png', sep='')
png(filename=plot.fname, width=1024, height=1024)
p <- ggplot(short.top.ten.5S, aes(x=seqid, y=short.read.count))
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('short reads top ten count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_text(angle=90))
p <- p + ggplot2::stat_bin(geom='text',
                           color='white',
                           aes(label=read.fraction, vjust=2))
print(p)
dev.off()

long.top.ten.5S$read.fraction <- round(long.top.ten.5S$long.read.count / total.long.5S, digits=2)

plot.fname <- paste(data.home, 'plots/long_top_5S.png', sep='')
png(filename=plot.fname, width=1024, height=1024)
p <- ggplot(long.top.ten.5S, aes(x=seqid, y=long.read.count))
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('long reads top ten count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_text(angle=90))
p <- p + ggplot2::stat_bin(geom='text',
                           color='white',
                           aes(label=read.fraction, vjust=2))
print(p)
dev.off()

all.top.ten.5S$read.fraction <- round(all.top.ten.5S$all.read.count / total.all.5S, digits=2)

plot.fname <- paste(data.home, 'plots/all_top_5S.png', sep='')
png(filename=plot.fname, width=1024, height=1024)
p <- ggplot(all.top.ten.5S, aes(x=seqid, y=all.read.count))
p <- p + ggplot2::geom_bar(stat='identity')
p <- p + ggtitle('all reads top ten count for 5S rRNA sequences')
p <- p + theme(axis.text.x=element_text(angle=90))
p <- p + ggplot2::stat_bin(geom='text',
                           color='white',
                           aes(label=read.fraction, vjust=2))
print(p)
dev.off()

# make barplots for with stacked long/short bars of union of top ten 
short.top.ten.seqids <- as(short.top.ten.5S$seqid, 'character')
long.top.ten.seqids <- as(long.top.ten.5S$seqid, 'character')
union.top.ten.seqids <- union(short.top.ten.seqids, long.top.ten.seqids)
union.top.ten <- rRNA5S.counts[rRNA5S.counts$seqid %in% union.top.ten.seqids,1:3]
colnames(union.top.ten)[2:3] <- c('short', 'long')
union.top.ten <- melt(union.top.ten)
colnames(union.top.ten)[2:3] <- c('length.category', 'count')

plot.fname <- paste(data.home, 'plots/top_covered_5S.png', sep='')
png(filename=plot.fname, width=1024, height=512)
p <- ggplot(union.top.ten, aes(x=seqid, y=count, fill=length.category)) +
    ggplot2::geom_bar(stat='identity') +
    ggtitle('top covered 5S rRNA sequences') +
    coord_flip() +
    ggplot2::annotate('text',
             x=length(union.top.ten.seqids) - 1,
             y=0.5*max(union.top.ten$count),
             label="5S rRNA",
             size=13, color='#3399FF')
print(p)
dev.off()

# find and count border spanning reads
borders <- c('ETS|28S', '28S|ITS2', 'ITS2|5.8S', '5.8S|ITS1', 'ITS1|18S')
neighbours <- list(c(1,2), c(2,3), c(3,4), c(4,5), c(5,6))
border.idx <- 1:5
border.counts <- data.frame(border=factor(),
                           count=numeric(),
                           length.category=factor())

# short reads
for (i in border.idx) {
    b <- borders[i]             
    n <- neighbours[[i]]
    left.reads <- short.reads[subjectHits(short.read.overlaps)[queryHits(short.read.overlaps) == n[1]]]
    right.reads <- short.reads[subjectHits(short.read.overlaps)[queryHits(short.read.overlaps) == n[2]]]
    count <- sum(names(left.reads) %in% names(right.reads))
    count.row = data.frame(border=b, count=count, length.category='short')
    border.counts <- rbind(border.counts, count.row)
}

# long reads
for (i in border.idx) {
    b <- borders[i]             
    n <- neighbours[[i]]
    left.reads <- long.reads[subjectHits(long.read.overlaps)[queryHits(long.read.overlaps) == n[1]]]
    right.reads <- long.reads[subjectHits(long.read.overlaps)[queryHits(long.read.overlaps) == n[2]]]
    count <- sum(names(left.reads) %in% names(right.reads))
    count.row = data.frame(border=b, count=count, length.category='long')
    border.counts <- rbind(border.counts, count.row)
}

# all reads
for (i in border.idx) {
    b <- borders[i]             
    n <- neighbours[[i]]
    left.reads <- all.reads[subjectHits(all.read.overlaps)[queryHits(all.read.overlaps) == n[1]]]
    right.reads <- all.reads[subjectHits(all.read.overlaps)[queryHits(all.read.overlaps) == n[2]]]
    count <- sum(names(left.reads) %in% names(right.reads))
    count.row = data.frame(border=b, count=count, length.category='all')
    border.counts <- rbind(border.counts, count.row)
}

plot.fname <- paste(data.home, 'plots/border_spanning_read_counts.png', sep='')
png(filename=plot.fname, width=1024, height=1024)
p <-  ggplot(border.counts, aes(x=border, y=count, fill=length.category))
p <- p + ggplot2::geom_bar(position='dodge', stat='identity') 
p <- p + ggtitle('border spanning read counts')
print(p)
dev.off()
