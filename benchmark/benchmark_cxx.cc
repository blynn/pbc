#include <fstream>
#include <sys/time.h>
#include <time.h>
#include "pbcxx.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using namespace pbc;

double pbc_get_time(void) {
  static struct timeval last_tv, tv;
  static int first = 1;
  static double res = 0;

  if (first) {
    gettimeofday(&last_tv, NULL);
    first = 0;
    return 0;
  } else {
    gettimeofday(&tv, NULL);
    res += tv.tv_sec - last_tv.tv_sec;
    res += (tv.tv_usec - last_tv.tv_usec) / 1000000.0;
    last_tv = tv;

    return res;
  }
}

// Helper function to read files, mainly for pairing params
string read_file(const char* filename) {
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) {
    cerr << "Could not read file " << filename << ", error: " << errno << endl;
    return "";
  }

  // Allocate just the size of the file
  string result;
  in.seekg(0, std::ios::end);
  result.resize(in.tellg());
  in.seekg(0, std::ios::beg);

  in.read(&result[0], result.size());
  in.close();
  return result;
}

int main(int argc, char **argv) {
  if (argc <= 1) {
    cerr << "Usage: " << argv[0] << " [params_file]" << endl;
    return 1;
  }

  string pairing_params = read_file(argv[1]);
  if (pairing_params.length() <= 0) {
    cerr << "No pairing params" << endl;
    return 2;
  }

  Pairing pairing(pairing_params);
  G1Element x(pairing);
  G2Element y(pairing);
  GTElement r(pairing);

  double t0, t1, ttotal;


  ttotal = 0.0;
  const size_t n = 10;
  for (size_t i = 0; i < n; ++i) {
    x.set_random();
    y.set_random();

    t0 = pbc_get_time();
    r = pairing(x, y);
    t1 = pbc_get_time();
    ttotal += t1 - t0;

    x.debug_print("x =");
    y.debug_print("y =");
    r.debug_print("r =");
  }
  printf("average pairing time = %f\n", ttotal / n);
  return 0;
}
