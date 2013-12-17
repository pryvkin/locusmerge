#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sam.h"

using namespace std;

// counts the total number of reads in a BAM file

int main(int argc, char **argv) {
  if (argc < 4) {
    cerr << "USAGE: " << argv[0] << " in.bam  max_nh  out.bam\n";
    return 1;
  }

  // open BAM file
  bamFile fp;
  if ((fp = bam_open(argv[1], "r")) == 0) {
    cerr << "Failed to open BAM file " << argv[1] << "\n";
    return 1;
  }

  int max_nh(atoi(argv[2]));

  bamFile out_bam;
  if ((out_bam = bam_open(argv[3], "wb")) == 0) {
    cerr << "Failed to open BAM file " << argv[3] << "\n";
    return 1;
  }

  bam_header_t *hdr = bam_header_read(fp);

  bam_header_write(out_bam, hdr);

  bam1_t *b = bam_init1();

  int n_alns_filtered(0);
  int n_reads_filtered(0);

  while(bam_read1(fp, b) > 0) {
    uint8_t *aux_data = bam_aux_get(b, "NH");
    int nh = bam_aux2i( aux_data );

    if (nh <= max_nh)
      bam_write1(out_bam, b);
    else {
      ++n_alns_filtered;
      int hi = bam_aux2i(bam_aux_get(b, "HI"));
      if (hi == 0)
	++n_reads_filtered;
    }
  }

  bam_close(out_bam);

  cout << n_alns_filtered << " alignments of " << n_reads_filtered << " reads with >" << max_nh << " hits removed.\n";

  bam_destroy1(b);
  bam_header_destroy(hdr);
  bam_close(fp);

  return 0;
}
