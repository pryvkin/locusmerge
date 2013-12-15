#!/usr/bin/env Rscript

args <- commandArgs(T)

if ( length(args) < 3 ) {
    write("USAGE: prioritize_annot.R raw_annot coverage outfile", stderr())
    q(status=1)
}

raw.annot.filename <- args[1]
cov.filename <- args[2]
annot.out.filename <- args[3]

out.priorities.filename <- sprintf("%s.priorities",
                                   raw.annot.filename)
out.cooccurrence.filename <- sprintf("%s.cooccurrence",
                                     raw.annot.filename)

raw.annot <- read.table(raw.annot.filename, sep="\t", as.is=T,
                    col.names=c("locus.id", "annot.cls", "annot.id"))
coverage <- read.table(cov.filename, sep="\t", as.is=T, row.names=4)
coverage <- coverage[,6,drop=F]
colnames(coverage) <- "coverage"

## take frequency of each class as its priority
## (higher priorities take precedence)

# NOTE: ties could be an issue
## class.freq <- table(raw.annot$annot.cls)
## class.freq["intergenic"] <- 0
## class.freq <- sort(class.freq, decreasing=T)
## tied.freq <- duplicated(class.freq)

## if (any(tied.freq)) {
##     tied.freq.values <- class.freq[tied.freq]
##     tied.cls <- names(class.freq)[class.freq %in% tied.freq.values]
##     tied.freq.values.bycls <- class.freq[class.freq %in% tied.freq.values]
##     write(sprintf("WARNING: Some classes tied for frequency: %s",
##                   paste(sprintf("%s (%d)", tied.cls, tied.freq.values.bycls),
##                         collapse=",")), stderr())
## }
## cls.pri <- class.freq

# instead, use coverage - more reads = higher pri
# NOTE: ties could be an issue
raw.annot$coverage <- coverage[raw.annot$locus.id,"coverage"]
cls.cov <- aggregate(coverage ~ annot.cls, raw.annot, median, na.rm=T)
row.names(cls.cov) <- cls.cov$annot.cls
tied.cov <- duplicated(cls.cov$coverage)
if (any(tied.cov)) {
    tied.cov.values <- cls.cov$coverage[tied.cov]
    tied.cls <- rownames(cls.cov)[cls.cov$coverage %in% tied.cov.values]
    tied.cov.values.bycls <- cls.cov$coverage[cls.cov$coverage %in% tied.cov.values]
    write(sprintf("WARNING: Some classes tied for median counts: %s",
                  paste(sprintf("%s (%d)", tied.cls, tied.cov.values.bycls),
                        collapse=",")), stderr())
}

cls.pri <- cls.cov[,2]
names(cls.pri) <- rownames(cls.cov)
cls.pri["intergenic"] <- 0

write.table(file=out.priorities.filename, data.frame(cls.pri),
            col.names=F, row.names=T, quote=F,
            sep="\t")

## output, for information purposes, the co-occurrence
## of annotation classes at loci
annot.pairs <- do.call(rbind, by(raw.annot, raw.annot$locus.id, function(a) {
    uniq.cls <- unique(a$annot.cls)
    if (length(uniq.cls) == 1)
        return(NULL)
    else
        t(combn(uniq.cls, 2))
} ) )

every.cls <- unique(raw.annot$annot.cls)
cooccurrence <- table(factor(annot.pairs[,1], levels=every.cls),
                      factor(annot.pairs[,2], levels=every.cls) )
write.table(file=out.cooccurrence.filename, cooccurrence,
            col.names=NA, row.names=T, quote=F,
            sep="\t")

## prioritize annotations: for each locus, take the highest pri cls
annot.pri <- do.call(rbind, by(raw.annot, raw.annot$locus.id, function(a) {
    all.cls <- paste(a$annot.cls, collapse=",")
    all.id <- paste(a$annot.id, collapse=",")
    if (length(unique(a$annot.cls)) > 1)
        raw.cls <- "multi"
    else
        raw.cls <- a$annot.cls[1]
    which.maxpri <- which.max(cls.pri[a$annot.cls])
    pri.cls <- a$annot.cls[which.maxpri]
    pri.id <- a$annot.id[which.maxpri]
    data.frame(a$locus.id[1], raw.cls, pri.cls, all.cls, all.id, pri.id)
} ) )

write.table(file=annot.out.filename, annot.pri,
            col.names=F, row.names=F, quote=F,
            sep="\t")

