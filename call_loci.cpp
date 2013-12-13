#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;

// gets a string from a tab delimited stream
void getstr(istream &is, string &s) {
  getline(is, s, '\t');
}

// gets an int from a tab delimited stream
void getint(istream &is, int &i) {
  string s;
  getline(is, s, '\t');
  i = atoi(s.c_str());
}

bool parse_argument(const string &arg, string &key, string &value) {
  key.clear();
  value.clear();
  if (arg.size() >= 3 && arg[0] == '-' && arg[1] == '-') {
    int i;
    for (i=2; i < arg.size() && arg[i] != '='; ++i)
	key += arg[i];
    for(++i; i < arg.size(); ++i)
      value += arg[i];
    return(true);
  } else
    return(false);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "USAGE: [OPTIONS] " << argv[0] << " in.bedgraph\n";
    cerr << "  --threshold=N  minimum read count [10]\n";
    cerr << "  --minrun=N     minimum consecutive coverage [14]\n";
    cerr << "  --maxgap=N     max gap between segments [44]\n";

    return(1);
  }

  string bamfn( argv[1] );
  int threshold = 10;
  int minrun = 14;
  int maxgap = 44;

  // parse command line arguments
  for (int i=1; i < argc; ++i) {
    string key, value;
    if (parse_argument(argv[i], key, value)) {
      if (key == "threshold") {
	threshold = atoi(value.c_str());
	if (threshold < 1) {
	  cerr << "Invalid threshold (" << value << "): must be >= 1\n";
	  return(1);
	}
      } else if (key == "minrun") {
	minrun = atoi(value.c_str());
	if (minrun < 1) {
	  cerr << "Invalid minrun (" << value << "): must be >= 1\n";
	  return(1);
	}
      } else if (key == "maxgap") {
	maxgap = atoi(value.c_str());
	if (maxgap < 1) {
	  cerr << "Invalid maxgap (" << value << "): must be >= 1\n";
	  return(1);
	}
      } else {
	cerr << "Invalid option: " << key << "\n";
	return(1);
      }
    }
  }

  cerr << "Using threshold=" << threshold << "\n"
       << "      minrun=" << minrun << "\n"
       << "      maxgap=" << maxgap << "\n";


  ifstream bedgraph(bamfn.c_str());
  if (!bedgraph.is_open()) {
    cerr << "Could not open file " << bamfn << "\n";
    return(1);
  }

  string line;
  string prev_chr;
  int first_start = -1;
  int prev_end = -1;

  while(getline(bedgraph, line)) {
    string chr;
    int start, end, count;
    istringstream line_str(line);
    getstr(line_str, chr);
    getint(line_str, start);
    getint(line_str, end);
    getint(line_str, count);

    if (!prev_chr.empty() && chr != prev_chr) {
      // end of chromosome: reset everything and process last run
      if (first_start >= 0 && ((prev_end - first_start) >= minrun) )
	cout << prev_chr << "\t" << first_start << "\t" << prev_end << "\n";
      first_start = -1;
    }

    if (count >= threshold) {
      if (first_start == -1) {
	// not currently in a run, so start one here
	first_start = start;
	//cout << "Starting (" << start << "," << end << "," << count << ")\n";
      } else {
	// currently in a run
	if ((start - prev_end) < maxgap) {
	  //cout << "Including (" << start << "," << end << "," << count <<  ")\n";
	  // include this chunk in the run
	} else {
	  // we've hit a gap too big, end the run
	  // output it if the run is long enough
	  if ( (prev_end - first_start) >= minrun )
	    cout << prev_chr << "\t" << first_start << "\t" << prev_end << "\n";
	  first_start = start;
	  //cout << "Ending run at (" << start << "," << end << "," << count <<  ")\n";
	}
      }
      prev_end = end;
    } else {
      //cout << "Skipping (" << start << "," << end << "," << count <<  ")\n";
    }

    prev_chr = chr;
  }
  return(0);
}
