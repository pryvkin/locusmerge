CXX = g++
CXXFLAGS = -Wall -O2

SAMTOOLS_CXXFLAGS = -Isamtools
SAMTOOLS_LDFLAGS = -Lsamtools -lbam -lz -lpthread

all: locusmerge_pairs locusmerge_similarity call_loci bam_filter_by_nh

locusmerge_pairs: locusmerge_pairs.cpp
	$(CXX) $(CXXFLAGS)  $< -o $@

locusmerge_similarity: locusmerge_similarity.cpp
	$(CXX) $(CXXFLAGS)  $< -o $@

call_loci: call_loci.cpp
	$(CXX) $(CXXFLAGS)  $< -o $@

bam_filter_by_nh: bam_filter_by_nh.cpp samtools/libbam.a
	$(CXX) $(CXXFLAGS) $(SAMTOOLS_CXXFLAGS)  $< -o $@ $(SAMTOOLS_LDFLAGS)

samtools/libbam.a: samtools/libbam.a
	cd samtools && make lib

clean: 
	rm -f locusmerge_pairs locusmerge_similarity call_loci bam_filter_by_nh
	cd samtools && make clean

