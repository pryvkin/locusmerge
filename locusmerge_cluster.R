#!/usr/bin/env Rscript

args <- commandArgs(T)

if (length(args) < 6) {
 cat("USAGE: ./locusmerge_cluster.R in.similarity in.annot in.cov in.ucov allowed_overlaps.txt out.clusters\n")
 q(status=1)
}

suppressPackageStartupMessages(library('igraph'))
suppressPackageStartupMessages(library('fastcluster'))

in.sim.filename <- args[1]
in.annot.filename <- args[2]
in.cov.filename <- args[3]
in.ucov.filename <- args[4]
in.allowed.filename <- args[5]
out.filename <- args[6]

write("Loading locus similarity scores...", stderr());
locus.similarities <- read.table(in.sim.filename, sep="\t", as.is=T)
write("Loading locus annotations...", stderr());
annot <- read.table(in.annot.filename, sep="\t", row.names=1, as.is=T,
                    col.names=c("locus.id", "raw.cls", "cls", "all.cls", "all.id", "id"))
write("Loading coverage files...", stderr());
load.coverage <- function(filename) {
    read.table(filename, sep="\t", row.names=4, as.is=T)[,6,drop=F]
}
coverage <- load.coverage(in.cov.filename)
uniqcov <- load.coverage(in.ucov.filename)
write("Loading allowed_class_pairs...", stderr());
allowed.cls.pairs <- read.table(in.allowed.filename, sep="\t", as.is=T,
                                comment.char="#", col.names=c("cls1", "cls2"))

# build class penalty matrix; if x_ij = 1 then avoid mixing
# those classes within a cluster
write("Building class penalty matrix...", stderr())
all.classes = unique(unlist(strsplit(annot$all.cls, ",")))
avoid.cls.pairs = outer(rep(1, length(all.classes)), rep(1, length(all.classes)))
rownames(avoid.cls.pairs) = all.classes
colnames(avoid.cls.pairs) = all.classes

specified.cls <- c(allowed.cls.pairs$cls1, allowed.cls.pairs$cls2)
specified.cls <- setdiff(specified.cls, "*")
missing.cls <- !(specified.cls %in% all.classes)
if (any(missing.cls)) {
    write(sprintf("WARNING: Non-existent classes were specified in %s: %s",
                  in.allowed.filename,
                  paste(specified.cls[missing.cls], collapse=",")), stderr())
}

allowed.cls.pairs <- subset(allowed.cls.pairs,
                            cls1 %in% all.classes & cls2 %in% all.classes)

cls.edge.list = as.matrix(allowed.cls.pairs[,c(1,2)])
allowed.cls.graph = graph.edgelist(cls.edge.list, FALSE)
components = decompose.graph(allowed.cls.graph)
dummy <- lapply(components, function(g) {
    x <- unique(V(g)$name)
    if (any(x == "*")) {
        # handle wildcards
        non.wild <- x[x != "*"]
        all.others <- setdiff(all.classes, non.wild)
        pairs <- expand.grid(non.wild, all.others)
    } else {
        pairs <- t(combn(x, 2))
    }
    apply(pairs, 1, function(p) {
        avoid.cls.pairs[p[1], p[2]] <<- 0
        avoid.cls.pairs[p[2], p[1]] <<- 0
    } )
})

relevant.cls <- setdiff(all.classes, "intergenic")

# compute symmetric similarity by taking the max of the nonsymmetric similarities
write("Symmetrizing locus similarity scores...", stderr())
rownames(locus.similarities) = apply(locus.similarities[,c(1,2)], 1, paste, collapse="@@@")

symm.sim <- c()
i <- 0
prev.pct.done <- 0
for(pair.name in rownames(locus.similarities)) {
  rev.pair.name <- paste(rev(unlist(strsplit(pair.name, "@@@"))), collapse="@@@")
  symm.sim[pair.name] <-
    max( c(locus.similarities[pair.name,4], locus.similarities[rev.pair.name,4]))
  pct.done <- round(100*i/nrow(locus.similarities))
  if (pct.done != prev.pct.done && (pct.done %% 2) == 0)
      cat(".")
  i <- i + 1
  prev.pct.done <- pct.done
}
cat("\n")

# save(symm.sim, file=file.path(dirname(out.filename), 'symmsim.Rdata'))
# load(file=file.path(dirname(out.filename), 'symmsim.Rdata'))

# remove redundant edges (since we'll be using a symmetric relation)
edge.list = as.matrix(locus.similarities[,c(1,2)])

edge.list2 = unlist(strsplit(unique(apply(edge.list,1, function(r) {
  paste(sort(r), collapse="@@@")
} ) ),"@@@"))
edge.list2 = matrix(edge.list2, nrow=length(edge.list2)/2, byrow=T)
# note: keep self-edges so we can keep singleton components in the list
rm(edge.list)


# construct the graph

g = graph.edgelist(edge.list2, FALSE)
components = decompose.graph(g)
rm(g)

# set vertex attributes
maxcov = rep(NA, length(components))
for(i in 1:length(components)) {
  g = components[[i]]
  V(g)$coverage = coverage[V(g)$name,1]
  V(g)$uniq.coverage = uniqcov[V(g)$name,1]
  V(g)$label = annot[V(g)$name,]$cls
  components[[i]] = g
  maxcov[i] = max(V(g)$coverage)
  rm(g)
}

# expressed components
min.cov = 20
xmap.err.rate = 0.01
expr.components = components[which(maxcov >= min.cov)]

is.singleton = unlist(lapply(expr.components, vcount)) == 1
singletons = expr.components[is.singleton]

# remove singletons from list of components before proceeding
expr.components = expr.components[!is.singleton]

# build tables
expr.component.tables = list()
for(i in 1:length(expr.components)) {
  expr.component.tables[[i]] = data.frame(id=V(expr.components[[i]])$name,
                         coverage=V(expr.components[[i]])$coverage,
                         uniq.cov=V(expr.components[[i]])$uniq.coverage,
                         annot.cls=V(expr.components[[i]])$label,
                         annot.id=annot[V(expr.components[[i]])$name,]$id)
  expr.component.tables[[i]] =
    expr.component.tables[[i]][order(-expr.component.tables[[i]]$coverage),]
}

write("Computing hierarchical clusterings...", stderr())
# hierarchical clustering into trees
htrees = list()
for(i in 1:length(expr.component.tables)) {
  ect = expr.component.tables[[i]]
  ids = as.character(ect$id)
  g = expr.components[[i]]

  dm = outer(1:vcount(g), 1:vcount(g))
  rownames(dm) = ids
  colnames(dm) = ids

  # populate distance matrix
  for(j in 1:nrow(ect)) {
    neibs = V(g)[nei(ids[j])]$name
    # distance = 1 - symmetric_similarity
    dm[ids[j], neibs] = 1 - symm.sim[paste(ids[j], neibs, sep="@@@")]
    rm(neibs)
  }

  dm = as.dist(dm)
  htrees[[i]] = hclust(dm, method="single")
  rm(dm)
}

# now cut the each tree i into K_i clusters
# where K_i satisfies some criterion:
# namely, we want the lowest K_i with the
# fewest "rogue" clusters

count.rogue.elements = function(ht, k, annot, avoid.cls.pairs, relevant.cls) {
  clu = cutree(ht, k=k)
  nr = 0
  for(j in 1:k) {
    annots = unique(annot[names(clu)[clu == j],]$cls)
    relevant.annots = annots[annots %in% relevant.cls]
    if (length(relevant.annots) < 2)
      next
    relevant.pairs = unique(combn(relevant.annots,2, FUN=paste, collapse="@@@"))
    offending.pairs = relevant.pairs[
      sapply(relevant.pairs, function(r) {
        s = unlist(strsplit(r, "@@@"))
        avoid.cls.pairs[s[1], s[2]] == 1
      } )]
    offending.cls = unlist(strsplit(offending.pairs, "@@@"))
    nr = nr + sum(annots %in% offending.cls)
    ## if (length(relevant.annots) > 1) {
    ##   for (u in 1:length(relevant.annots)) {
    ##     if (any(avoid.cls.pairs[relevant.annots[u], relevant.annots] == 1)) {
    ##       nr = nr + 1
    ##       break
    ##     }
    ##   }
    ## }
  }
  return( nr )
}
count.rogue.elements.in.cluster = function(ids, annot, avoid.cls.pairs, relevant.cls) {
  annots = unique(annot[ids,]$cls)
  relevant.annots = annots[annots %in% relevant.cls]
  if (length(relevant.annots) < 2)
    return (0)
  relevant.pairs = unique(combn(relevant.annots,2, FUN=paste, collapse="@@@"))
  offending.pairs = relevant.pairs[
    sapply(relevant.pairs, function(r) {
                        s = unlist(strsplit(r, "@@@"))
                        avoid.cls.pairs[s[1], s[2]] == 1
                      } )]
  offending.cls = unlist(strsplit(offending.pairs, "@@@"))
  return ( sum(annots %in% offending.cls) )
}

write("Optimizing k...", stderr())

pdf(file.path(dirname(out.filename), "k_optimization.pdf"), points=12)
par(mfrow=c(2,2))
clusterings = list()
for(i in 1:length(expr.component.tables)) {
  ect = expr.component.tables[[i]]
  ht = htrees[[i]]
  # let max_k be (1+) the # of loci with rel. classes or 0.5 * no. of loci
  max_k = min( 1 + sum(ect$annot.cls %in% relevant.cls),
    round(nrow(ect)/2) )

  # do a binary search for the last point where n decreases (nr[k] < nr[k-1])
  min_nr = count.rogue.elements(ht, max_k, annot, avoid.cls.pairs, relevant.cls)
  rogue_by_k = rep(NA, max_k)
  rogue_by_k[max_k] = min_nr

  optimal.k = 1
  
  lower_k = 1
  upper_k = max_k
  k = ceiling(max_k/2)
  while(lower_k < upper_k) {
    nr = count.rogue.elements(ht, k, annot, avoid.cls.pairs, relevant.cls)
    #cat(sprintf("%d    k=%d, nr=%d; (%d, %d)\n", i, k, nr, lower_k, upper_k))
    rogue_by_k[k] = nr
    k_dir = 0
    if (nr <= min_nr) {
      # we've shrunk the upper bound on k
      upper_k = k
      k_dir = -1
    } else if (nr > min_nr) {
      # we've increased the lower bound on k
      lower_k = k
      k_dir = 1
    }
    if ((lower_k == upper_k))
      break
    if ( (k > 1) && !any(is.na(rogue_by_k[c(k, k-1)])) && (rogue_by_k[k] < rogue_by_k[k-1])) {
      optimal.k = k
      break
    }
    if ( (k < max_k) && !any(is.na(rogue_by_k[c(k, k+1)])) && (rogue_by_k[k] > rogue_by_k[k+1])) {
      optimal.k = k+1
      break
    }
    k = ceiling((lower_k + upper_k) / 2)
    # ensure we move to a new k if this one was already computed
    if (!is.na(rogue_by_k[k]))
      k = k + k_dir
  }

  clu <- cutree(ht, k=optimal.k)
  clusterings[[i]] <- clu
  
  # compute inter-cluster scores (proxy for crossmapping between clusters)
  # we want this to be low.
  ## if (length(unique(clu)) > 1) {
  ##     clu.pairs <- combn(unique(clu), 2) 
  ##     clu.max.pairwise.sim <- apply(clu.pairs, 2, function(cluster.pair) {
  ##         clu.u <- names(clu)[clu == cluster.pair[1]]
  ##         clu.v <- names(clu)[clu == cluster.pair[2]]
  ##         locus.id.pairs <- as.vector(outer(clu.u, clu.v, function(locus.a, locus.b) {
  ##             z <- t(apply(cbind(locus.a, locus.b), 1, sort))
  ##             z <- paste(z[,1], z[,2], sep="@@@")
  ##             } ))
  ##         locus.sim.scores <- locus.similarities[locus.id.pairs , 4]
  ##         locus.sim.scores[is.na(locus.sim.scores)] <- 0
  ##         max(locus.sim.scores)
  ##     } )
  ##     hist(clu.max.pairwise.sim)
  ## } else {
  ##     plot(0, 0, pch='')
  ## }
  
  
  plot.x <- (1:length(rogue_by_k))[!is.na(rogue_by_k)]
  plot.y <- rogue_by_k[!is.na(rogue_by_k)]

  plot(plot.x, plot.y,
       main=sprintf("Expr component %d, n=%d", i, nrow(ect)),
       xlab="k", ylab="N_rogue", ylim=c(0,max(plot.y)))
  lines(plot.x, plot.y)
  abline(v = optimal.k, lty=2, col=2)
  
}

dev.off()

cluster.annot = data.frame(clu.id=NULL, nmemb=NULL, nrogue=NULL,
  main.id=NULL, main.cov=NULL, main.uniq.cov=NULL,
  main.annot.cls = NULL, main.annot.id=NULL,
  other.id=NULL, low.cov=NULL)

comma.sep = function(x) { paste(x, collapse=",") }

nclu = 1
for(i in 1:length(expr.component.tables)) {
  ect = expr.component.tables[[i]]
  clu = clusterings[[i]]
  for(j in 1:max(clu)) {
    which.clu = clu == j
    clu.memb = ect[which.clu,]
    which.main = which.max(clu.memb$coverage)
    other.id = NA
    if (sum(which.clu) > 1) {
      other.id =  comma.sep(clu.memb[-which.main, 'id'])
    }
    cluster.annot = rbind(cluster.annot, data.frame(
      clu.id=sprintf("CL%0.6d.%0.6d", i, j),
      nmemb = sum(clu==j),
      nrogue = count.rogue.elements.in.cluster(as.character(clu.memb$id), annot,
        avoid.cls.pairs, relevant.cls),
      main.id = as.character(clu.memb[which.main, 'id']),
      main.cov = clu.memb[which.main, 'coverage'],
      main.uniq.cov = clu.memb[which.main, 'uniq.cov'],
      main.annot.cls = as.character(clu.memb[which.main, 'annot.cls']),
      main.annot.id = as.character(clu.memb[which.main, 'annot.id']),
      other.id = other.id,
      low.cov = clu.memb[which.main, 'coverage'] < min.cov,
      stringsAsFactors=F))
    nclu = 1 + nclu
  }
}

# add singletons
for(i in 1:length(singletons)) {
  id = V(singletons[[i]])[1]$name
  cluster.annot = rbind(cluster.annot, data.frame(
    clu.id=sprintf("CL%0.6d.%0.6d", length(expr.component.tables) + i, j),
    nmemb = 1, nrogue = 0, main.id = id,
    main.cov = coverage[id,1], main.uniq.cov = uniqcov[id,1],
    main.annot.cls = annot[id,]$cls, main.annot.id = annot[id,]$id,
    other.id = NA, low.cov = coverage[id,1] < min.cov, stringsAsFactors=F))
  nclu = 1 + nclu
}

write.table(cluster.annot, file=out.filename, sep="\t", col.names=F, row.names=F, quote=F)


