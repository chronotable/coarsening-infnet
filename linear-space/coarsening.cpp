#include "coarsening.hpp"

using namespace std;

inline long long int make_pair_value(const Node n, const int p0, const int p1){
  return ((long long int)n) * p0 + p1;
}

bool CoarseningInfluenceGraph::run_linear(const int num_threads, const vector<pair<pair<Node, Node>, double> > &es, const int R, const string output_prefix){
  const double start_time = gettimeofday_sec();
  Node n = 0;
  long long int m = es.size();
  for (long long int i = 0; i < m; i++)
    n = max(n, (Node)max(es[i].first.first, es[i].first.second) + 1);

  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  //sort(es.begin(), es.end());

  ///// initialize
  vector<vector<long long int> > at_e(num_threads, vector<long long int>(n + 1, 0));
  vector<vector< pair<Node, Node> > > r_es(num_threads);

  pair<Node, Node>* id_index = new pair<Node, Node>[n]; // id and index
  for(Node i = 0; i < n; i++){
    id_index[i].first = 0; // initially, all the nodes belong to a same group
    id_index[i].second = i;
  }
  vector< long long int > id_pair(n);
  vector< vector<Node> > comp(num_threads, vector<Node>(n, -1)); // component id for SCC decomposition

  // parallel
  const int num_loop = (R + num_threads - 1) / num_threads;
  cout << "start creation of random graphs" << endl;
  for(int nl = 0; nl < num_loop; nl++){
    const int threads_need = min(num_threads, R - nl * num_threads);
    omp_set_num_threads(threads_need);
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      const int t = tid + nl * num_threads;
      if(t < R){
        Xorshift xs = Xorshift(t);
        r_es[tid].clear();

        at_e[tid].assign(n + 1, 0);
        for(long long int i = 0; i < m; i++){
          if (xs.nextDouble() < es[i].second){
            r_es[tid].push_back(es[i].first);
            at_e[tid][es[i].first.first + 1]++;
          }
        }
        at_e[tid][0] = 0;
        for(Node i = 1; i <= n; i++)
          at_e[tid][i] += at_e[tid][i - 1];

        comp[tid].assign(n, -1);
        scc_inmemory(n, r_es[tid], at_e[tid], comp[tid]);
        r_es[tid].clear();
      }
#pragma omp barrier
    }

    // tournament merge
    int merge_loop = 0;
    for(int i = 1; i < num_threads; i<<=1, merge_loop++);
    int step = 1;
    for(int ml = 0; ml < merge_loop; ml++){
      const int next = step << 1;
      const int nt = (threads_need + next - 1) / next;
      omp_set_num_threads(nt);
#pragma omp parallel
      {
        const int tid = omp_get_thread_num();
        const int index = tid * next;
        const int neighbor = index + step;
        if(neighbor < threads_need){
          vector<long long int> cp(n); // comp pair id
          for(Node u = 0; u < n; u++)
            cp[u] = make_pair_value(n, comp[index][u], comp[neighbor][u]);
          sort(cp.begin(), cp.end());
          cp.erase(unique(cp.begin(), cp.end()));
          unordered_map< long long int, int> pair2nid; // pair to new id
          int id_count = 0;
          for(vector<long long int >::iterator it = cp.begin(); it != cp.end(); it++){
            pair2nid[*it] = id_count++;
          }
          // new index
          for(Node u = 0; u < n; u++)
            comp[index][u] = pair2nid[make_pair_value(n, comp[index][u], comp[neighbor][u])];
        }
#pragma omp barrier
      }
      step = next;
    }

    // sum-up merge
    id_pair.resize(n);
    for(Node u = 0; u < n; u++)
      id_pair[u] = make_pair_value(n, id_index[u].first, comp[0][u]);
    sort(id_pair.begin(), id_pair.end());
    id_pair.erase(unique(id_pair.begin(), id_pair.end()), id_pair.end());
    unordered_map< long long int, int> pair2nid; // pair to new id
    int id_count = 0;
    for(vector<long long int >::iterator it = id_pair.begin(); it != id_pair.end(); it++){
      pair2nid[*it] = id_count++;
    }
    // new index
    for(Node u = 0; u < n; u++)
      id_index[u].first = pair2nid[make_pair_value(n, id_index[u].first, comp[0][u])];
    fflush(stdout);
	}

  sort(id_index, id_index + n);
  vector<Node> mapping(n, -1);
  vector<Node> rscc_size;
  Node num_cnodes = 0; // the number of nodes after coarsening
  for(Node u = 0; u < n;){
    const Node id = id_index[u].first;
    Node size = 0;
    for(; u < n && id_index[u].first == id; u++){
      mapping[id_index[u].second] = num_cnodes;
      size++;
    }
    num_cnodes++;
    rscc_size.push_back(size);
  }
  printf("num_cnodes=%d (%.10f)\n", num_cnodes, 1.0 * num_cnodes / n);

  // reduced outer edges (SCC to SCC edges)
  unordered_map<long long int, pair<int, double> > sccedges; // key(edge), <frequency, (1-p)^freq>
  for(long long int i = 0; i < m; i++){
    const Node u = mapping[es[i].first.first], v = mapping[es[i].first.second];
    if(u != v && (rscc_size[u] > 1 || rscc_size[v] > 1)){
      long long int key = 1ll * u * n + v;
      if(sccedges.find(key) == sccedges.end()){
        sccedges[key].first = 0;
        sccedges[key].second = 1.0;
      }
      sccedges[key].first++;
      sccedges[key].second *= 1.0 - es[i].second;
    }
  }

  /////
  // write new graph
  const string graph_name = output_prefix + "_coarsened-graph.txt";
  const size_t vbuf_size = 1<<25;
  static char vbuf[vbuf_size];
  FILE* fp = fopen(graph_name.c_str(), "w");
  if(fp == NULL){
    cerr << "cannot open the file: " << graph_name << endl;
    return false;
  }
  setvbuf(fp, vbuf, _IOFBF, vbuf_size);
  vector< pair<pair<Node, Node>, double> > edgelist;
  // edges between nodes with rscc size 1
  for(long long int i = 0; i < m; i++){
    const Node u = mapping[es[i].first.first], v = mapping[es[i].first.second];
    const double p = es[i].second;
    if(u != v && rscc_size[u] == 1 && rscc_size[v] == 1){
      edgelist.push_back(make_pair(make_pair(u, v), p));
    }
  }
  // otherwise
  for(unordered_map< long long int, pair<int, double> >::iterator it = sccedges.begin(); it != sccedges.end(); it++){
    const Node u = (Node)(it->first / n);
    const Node v = (Node)(it->first % n);
    const double p = 1.0 - it->second.second;
    edgelist.push_back(make_pair(make_pair(u, v), p));
  }
  sort(edgelist.begin(), edgelist.end());
  for(size_t i = 0; i < edgelist.size(); i++)
    fprintf(fp, "%u %u %.10e\n", edgelist[i].first.first, edgelist[i].first.second, edgelist[i].second);
  fclose(fp);
  printf("num_cedges=%lld\n", (long long int)edgelist.size());

  // write node mapping (old node ID to new node ID)
  const string mapping_name = output_prefix + "_mapping.txt";
  fp = fopen(mapping_name.c_str(), "w");
  if(fp == NULL){
    cerr << "cannot open the file: " << mapping_name << endl;
    return false;
  }
  setvbuf(fp, vbuf, _IOFBF, vbuf_size);
  fprintf(fp, "%u\n", n);
  for(Node i = 0; i < n; i++){
    fprintf(fp, "%u\n", mapping[i]);
  }
  fclose(fp);

  // finalize
  delete [] id_index;

  printf("all finished: %.5f\n", gettimeofday_sec() - start_time);
  fflush(stdout);

  return true;
}

bool less_pair_second(pair<Node, Node> left, pair<Node, Node> right){
  if(left.second != right.second)
    return left.second < right.second;
  return left.first < right.first;
}

// side effect: es and at_e will be modified
int CoarseningInfluenceGraph::scc_inmemory(const Node n, vector< pair<Node, Node> >& es, vector<long long int>& at_e, vector<Node>& comp) {
  vector<bool> vis(n, false);
  vector< pair<Node, int> > S;
  vector<int> lis(n, 0);
  Node lis_size = 0;
  int k = 0;
  for(Node i = 0; i < n; i++)
    S.push_back(make_pair(i, 0));
  for(; !S.empty();) {
    const Node v = S.back().first;
    int state = S.back().second;
    S.pop_back();
    if(state == 0){
      if(vis[v])
        continue;
      vis[v] = true;
      S.push_back(make_pair(v, 1));
      for(long long int i = at_e[v]; i < at_e[v + 1]; i++){
        const Node u = es[i].second;
        S.push_back(make_pair(u, 0));
      }
    }else{
      lis[lis_size++] = v;
    }
  }
  for(Node i = 0; i < n; i++)
    S.push_back(make_pair(lis[i], -1));

  // reverse the edges
  sort(es.begin(), es.end(), less_pair_second);
  for(Node i = 0; i <= n; i++)
    at_e[i] = 0;
  for(size_t i = 0; i < es.size(); i++)
    at_e[es[i].second + 1]++;
  at_e[0] = 0;
  for(Node i = 1; i <= n; i++)
    at_e[i] += at_e[i - 1];
  vis.assign(n, false);
  for(; !S.empty();){
    Node v = S.back().first;
    int arg = S.back().second;
    S.pop_back();
    if(vis[v])
      continue;
		vis[v] = true;
		comp[v] = arg == -1 ? k++ : arg;
    for(long long int i = at_e[v]; i < at_e[v + 1]; i++){
      Node u = es[i].first;
      S.push_back(make_pair(u, comp[v]));
    }
  }
  return k;
}
