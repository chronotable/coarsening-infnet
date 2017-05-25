#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>

#include <omp.h>


using namespace std;

typedef unsigned int Node;
// wall clock time
inline double gettimeofday_sec(){
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
}
// user time
inline double getutime(){
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return (double)r.ru_utime.tv_sec + (double)r.ru_utime.tv_usec * 1e-6;
}
// system time
inline double getstime(){
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return (double)r.ru_stime.tv_sec + (double)r.ru_stime.tv_usec * 1e-6;
}


class CoarseningInfluenceGraph{
private:
  int scc_inmemory(const Node n, vector< pair<Node, Node> >& es, vector<long long int>& at_e, vector<Node>& comp);
public:
  bool run_linear(const int num_threads, const vector<pair<pair<Node, Node>, double> > &es, const int R, const string output_prefix);
};


// Random Number Generator
class Xorshift{
public:
  Xorshift(int seed) {
    x = _(seed, 0);
    y = _(x, 1);
    z = _(y, 2);
    w = _(z, 3);
  }

  int _(int s, int i) {
    return 1812433253 * (s ^ (s >> 30)) + i + 1;
  }

  // 32bit signed
  inline int nextInt() {
    unsigned int t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    return w = w ^ (w >> 19) ^ t ^ (t >> 8);
  }

  // error = O(n*2^-32)
  inline int nextInt(int n) {
    return (int) (n * nextDouble());
  }

  // [0, 1) (53bit)
  inline double nextDouble() {
    unsigned int a = ((unsigned int) nextInt()) >> 5, b =
	 		((unsigned int) nextInt()) >> 6;
    return (a * 67108864.0 + b) * (1.0 / (1LL << 53));
  }

  // 64bit signed
  inline int nextLong(){
    return (((long long int)nextInt()) << 32ll) | ((long long int)nextInt());
  }

private:
  unsigned int x, y, z, w;
};


