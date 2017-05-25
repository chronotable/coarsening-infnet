/**
 * the sub-linear-space implementation of our paper "Coarsening Massive Influence Networks for Scalable Diffusion Analysis"
 *
 * usage: ./main graph_file R output_prefix num_threads
 *
 * @param graph_file: name of graph file
 ** format (in text): each line contains "src dst edge-probability"
 ** src and dst are indexed 0 to N-1 (N is the number of vertices)
 * @param R: r of r-robust scc
 * @param output_prefix: resulting output file name prefix
 ** "[output_prefix]_coarsened-graph.txt" is the coarsened graph
 ** "[output_prefix]_mapping.txt" is the correspondence mapping π : V -> W (see algorithm 2 in our paper)
 *** header contains the number of vertices in the original graph, then next lines contain the mapping of each vertex
 * @param num_threads: the number of running threads
 *
 * attention: this program generates some intermediate files (random_graph_[number].bin and swsscc_tid[number]_[number].bin). remove them after the execution if necessary.
 **/

#include "coarsening.hpp"

using namespace std;


int main(int argc, char* argv[]){
	if (argc < 5) {
    printf("usage: %s graph R output_prefix num_threads\n", argv[0]);
		exit(1);
	}

  double start_time = gettimeofday_sec();

  char* file_name = argv[1];
	int R = atoi(argv[2]);
  string output_prefix = argv[3];
  const int num_threads = atoi(argv[4]);

  FILE* file_fp = fopen(file_name, "r"); //open(file_name, O_RDONLY);
  if(file_fp == NULL){
    fprintf(stderr, "cannot open the file: %s\n", file_name);
    return 1;
  }
	Node u, v;
	double p;
  Node n = 0;
  long long int m = 0;
  for(; ~fscanf(file_fp, "%u%u%lf", &u, &v, &p); ){
    if(u == v)
      continue;
    n = max(n, max(u, v));
    m++;
  }
  n++;
  fclose(file_fp);

  CoarseningInfluenceGraph cig;
  cig.run_sublinear(num_threads, file_name, R, output_prefix, n, m);

  // finalyze

  double wc_time = gettimeofday_sec() - start_time;
  printf("wall_clock_time=%.10f\n", wc_time);
  double stime = getstime();
  printf("system_time=%.10f\n", stime);
  double utime = getutime();
  printf("user_time=%.10f\n", utime);

  return 0;
}
