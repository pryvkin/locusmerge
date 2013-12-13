#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

struct pair_entry {
  string locus_i;
  string locus_j;
  unsigned int count;
};

void process_pairs(vector<pair_entry> &pairs, unsigned int max_count) {
  const unsigned int n = pairs.size();
  for(unsigned int i=0; i < n; ++i)
    cout << pairs[i].locus_i << "\t"
	 << pairs[i].locus_j << "\t"
	 << max_count << "\t"
	 << double(pairs[i].count) / double(max_count) << "\n";
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " pairs\n";
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
  string prev_locus_i;
  unsigned int max_count(0);
  vector<pair_entry> pairs;

  while(getline(in, line)) {
    istringstream line_stream(line);
    string dummy;
    pair_entry pe;
    getline(line_stream, pe.locus_i, '\t'); // locus id

    if (!prev_locus_i.empty() && (pe.locus_i != prev_locus_i)) {
      // we're at the next locus_i
      process_pairs(pairs, max_count);
      pairs.clear();
      max_count = 0;
    }

    getline(line_stream, pe.locus_j, '\t'); // locus id
    getline(line_stream, dummy, '\t'); // count
    pe.count = atoi(dummy.c_str());
    pairs.push_back(pe);

    if (pe.count > max_count)
      max_count = pe.count;

    prev_locus_i = pe.locus_i;
  }
  process_pairs(pairs, max_count);

  return 0;
}
