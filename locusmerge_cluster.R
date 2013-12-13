#!/usr/bin/env Rscript

if (length(args) < 5) {
 cat("USAGE: ./locusmerge_cluster.R in.similarity in.annot in.cov in.ucov out.clusters\n")
 q(status=1)
}

args <- commandArgs(T)

library('igraph')
library('fastcluster')

in.sim.filename <- args[1]
in.annot <- args[2]
in.cov.filename <- args[3]
in.ucov.filename <- args[4]
out.filename <- args[5]

annot=read.table('test.annot',sep="\t",as.is=T,row.names=1)
coverage=read.table('test2.cov',sep="\t",row.names=1)
uniqcov=read.table('test2.uniqcov',sep="\t",row.names=1)

locus.similarities = read.table("/home/pry/data/tmp/locusmerge/test_all.sim", as.is=T, sep="\t");

# avoid clusters containing 2+ of the following locus classes
all.classes = unique(unlist(strsplit(annot[,3], ",")))
avoid.cls.pairs = outer(rep(0, length(all.classes)), rep(0, length(all.classes)))
rownames(avoid.cls.pairs) = all.classes

colnames(avoid.cls.pairs) = all.classes
avoid.cls.pairs["miRNA",   c("tRNA", "mt-tRNA", "snoRNA", "snRNA", "rRNA")] = 1
avoid.cls.pairs["tRNA",    c("miRNA", "snoRNA", "snRNA", "rRNA")] = 1
avoid.cls.pairs["mt-tRNA", c("miRNA", "snoRNA", "snRNA", "rRNA")] = 1
avoid.cls.pairs["snoRNA",  c("miRNA", "tRNA", "mt-tRNA", "snRNA","rRNA")] = 1
avoid.cls.pairs["snRNA",   c("miRNA", "tRNA", "mt-tRNA", "snoRNA", "rRNA")] = 1
avoid.cls.pairs["rRNA",    c("miRNA", "tRNA", "mt-tRNA", "snoRNA", "snRNA")] = 1
relevant.cls = union(rownames(avoid.cls.pairs)[rowSums(avoid.cls.pairs) > 0],
  colnames(avoid.cls.pairs)[colSums(avoid.cls.pairs) > 0] )

# compute symmetric similarity by taking the max of the nonsymmetric similarities
rownames(locus.similarities) = apply(locus.similarities[,c(1,2)], 1, paste, collapse="@@@")

symm.sim = c()
i = 0
for(pair.name in rownames(locus.similarities)) {
  rev.pair.name = paste(rev(unlist(strsplit(pair.name, "@@@"))), collapse="@@@")
  symm.sim[pair.name] =
    max( c(locus.similarities[pair.name,4], locus.similarities[rev.pair.name,4]))
  if (i %% 10 == 0)
    print(100*i/nrow(locus.similarities))
  i = i +1
}

#save(symm.sim, file='/home/pry/data/tmp/locusmerge/symmsim.Rdata')
#load(file='/home/pry/data/tmp/locusmerge/symmsim.Rdata')

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
  V(g)$label = annot[V(g)$name,1]
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
                         annot.id=annot[V(expr.components[[i]])$name,2])
  expr.component.tables[[i]] =
    expr.component.tables[[i]][order(-expr.component.tables[[i]]$coverage),]
}

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
    annots = unique(annot[names(clu)[clu == j],]$V2)
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
  annots = unique(annot[ids,]$V2)
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
    cat(sprintf("%d    k=%d, nr=%d; (%d, %d)\n", i, k, nr, lower_k, upper_k))
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
  clusterings[[i]] = cutree(ht, k=optimal.k)
}

####
# TODO: ADD SINGLETONS
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
    main.annot.cls = annot[id, 1], main.annot.id = annot[id, 2],
    other.id = NA, low.cov = coverage[id,1] < min.cov, stringsAsFactors=F))
  nclu = 1 + nclu
}

write.table(cluster.annot, file=out.filename, sep="\t", col.names=F, row.names=F, quote=F)


