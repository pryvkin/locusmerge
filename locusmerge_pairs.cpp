#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

void process_locus_ids(vector<string> &locus_ids) {
  sort(locus_ids.begin(), locus_ids.end());
  const unsigned int n( locus_ids.size() );
  for(unsigned int i=0; i < n; ++i) {
    for(unsigned int j=0; j < n; ++j) {
      //if (i != j)
      cout << locus_ids[i] << "\t" << locus_ids[j] << "\n";
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " bam_vs_loci.intersect\n";
    return 1;
  }

  istream *inp;
  ifstream infile;
  if (string(argv[1]) == "-") {
    inp = &cin;
  } else {
    infile.open(argv[1]);
    inp = &infile;
  }
  istream &in(*inp);

  // track every pair of cross-mapping loci
  string line;
  string prev_read_id;
  vector<string> locus_ids; // all loci that current read maps to
  while(getline(in, line)) {
    istringstream line_stream(line);
    string read_id, locus_id;
    getline(line_stream, read_id, '\t'); // read id

    if (!prev_read_id.empty() && (read_id != prev_read_id)) {
      // we're at the next read id
      process_locus_ids(locus_ids);
      locus_ids.clear();
    }

    getline(line_stream, locus_id, '\t'); // locus id
    locus_ids.push_back(locus_id);

    prev_read_id = read_id;
  }
  process_locus_ids(locus_ids);

  return 0;
}
