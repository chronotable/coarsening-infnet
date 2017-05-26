#include "coarsening.hpp"

using namespace std;


// Semi-w-streaming SCC Union-Find
// see the paper: L. Laura and F. Santaroni. Computing strongly connected components in the streaming model. In TAPAS, pages 193–205, 2011.
class SSUF{
public:
  SSUF(Node _n){
    n = _n;
    a = new Node[n+1];
    t = new Node[n+1];
    init();
  }
  ~SSUF(){
    delete [] a;
    delete [] t;
  }
  void init(){
    for(Node i = 0; i <= n; i++){
      a[i] = n; // dummy root
      t[i] = i;
    }
  }
  inline Node translate(const Node v){
    if(t[v] != v)
      t[v] = translate(t[v]);
    return t[v];
  }
  inline void setTranslate(const Node v, const Node nt){
    t[translate(v)] = translate(nt);
  }
  inline Node parent(const Node v){
    assert(v <= n);
    //return a[v];
    return translate(a[translate(v)]);
  }
  inline void setParent(const Node v, const Node par){
    a[translate(v)] = translate(par);
  }
  bool existParent(const Node v, const Node target){
    assert(v < n && target < n);
    Node cur = translate(v);
    while(cur != n){
      if(cur == translate(target))
        return true;
      cur = parent(translate(cur));
    }
    return false;
  }
  int distance2root(const Node v){
    assert(v < n);
    int ret = 0;
    Node cur = translate(v);
    while(cur != n){
      cur = parent(translate(cur));
      ret++;
    }
    return ret;
  }
  bool loop(const Node v){
    map<Node, int> used;
    Node cur = translate(v); //v; //translate(v);
    while(cur != n){
      if(used[cur])
        return true;
      used[cur] = 1;
      cur = parent(translate(cur));
    }
    return false;
  }
private:
  Node n;
  Node* a; // parent
  Node* t; // translate
};
int runSWSSCC(const char* file_name, const Node num_vertices, int* cid, BinaryGraphReader<Node>* _bgr, BinaryWriter<Node>* _bw, const int tid){
  BinaryGraphReader<Node>& bgr = *_bgr;
  BinaryWriter<Node>& bw = *_bw;
  int infd = open(file_name, O_RDONLY);
  if(infd < 0){
    fprintf(stderr, "cannot open the file: %s\n", file_name);
    exit(EXIT_FAILURE);
  }
  SSUF uf(num_vertices);

  // start
  bgr.setfd(infd);
  Node src, dst;
  //const int max_pass = 100;
  int pass = 1;
  //for(pass = 1; pass < max_pass; pass++){
  for(pass = 1; ; pass++){
    bool changed = false;
    const string next_file = "swsscc_tid" + to_string((long long int)tid) + "_" + to_string((long long int)(pass % 2)) + ".bin";
    int nfd = open(next_file.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE);
    if(nfd < 0){
      fprintf(stderr, "cannot open the file: %s\n", next_file.c_str());
      exit(EXIT_FAILURE);
    }
    bw.setfd(nfd);
    while(bgr.nextEdge()){
      src = bgr.getSrc();
      dst = bgr.getDst();
      assert(src < num_vertices && dst < num_vertices);
      const Node u = uf.translate(src);
      const Node v = uf.translate(dst);

      if(u == v){ // case 5: self loop
        continue;
      }
      if(uf.existParent(v, u)){ // case 4: forward edge
        continue;
      }
      if(uf.existParent(u, v)){ // case 1: backward edge
        Node cur = u;
        while(cur != v){
          Node next = uf.parent(cur);
          uf.setTranslate(cur, v);
          cur = next;
        }
        changed = true;
        continue;
      }
      if(uf.parent(v) == num_vertices){ // case 0
        uf.setParent(v, u);
        changed = true;
        continue;
      }
      const int dis_u = uf.distance2root(u);
      const int dis_v = uf.distance2root(v);
      if(dis_u < /*<=*/ dis_v){ // case 2: cross forward edge
        bw.bwwrite(&u, 1);
        bw.bwwrite(&v, 1);
      }else{ // case 3: cross NON forward edge
        const Node pv = uf.parent(v); // parent of v
        assert(pv != num_vertices);
        bw.bwwrite(&pv, 1);
        bw.bwwrite(&v, 1);
        uf.setParent(v, u);
        changed = true;
      }
    }
    bw.bwflush();
    close(infd);
    close(nfd);
    if(!changed)
      break;
    infd = open(next_file.c_str(), O_RDONLY);
    bgr.setfd(infd);
  }
  //if(pass == max_pass){fprintf(stderr, "could not finish withing %d passes\n", max_pass); exit(EXIT_FAILURE);}

  // output
  map<Node, Node> node2cid;
  Node counter = 0;
  for(Node i = 0; i < num_vertices; i++){
    Node t = uf.translate(i);
    if(node2cid.find(t) == node2cid.end())
      node2cid[t] = counter++;
    cid[i] = node2cid[t];
  }
  return (int)counter;
}
/////


//inline long long int make_pair_value(const Node n, const int p0, const int p1){
inline long long int make_pair_value(const Node n, const Node p0, const Node p1){
  return ((long long int)n) * p0 + p1;
}

bool CoarseningInfluenceGraph::run_sublinear(const int num_threads, const char* file_name, const int R, const string output_prefix, const Node n, const long long int m) {
  assert(num_threads > 0);
  const double start_time = gettimeofday_sec();
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;

  pair<Node, Node>* id_index = new pair<Node, Node>[n]; // id and index
  for(Node i = 0; i < n; i++){
    id_index[i].first = 0; // initially, all the nodes belong to a same group
    id_index[i].second = i;
  }
  vector< long long int > id_pair(n);
  int** comp = new int*[num_threads];
  for(int i = 0; i < num_threads; i++)
    comp[i] = new int[n];

  // parallel
  BinaryWriter<Node>** random_graph_bws = new BinaryWriter<Node>*[num_threads];
  for(int i = 0; i < num_threads; i++)random_graph_bws[i] = new BinaryWriter<Node>(1ll<<24);
  vector<string> random_graph_files(R);
  vector<int> nfds(R);
  Xorshift** xss = new Xorshift*[R];
  for(int i = 0; i < R; i++)
    xss[i] = new Xorshift(i);
  for(int i = 0; i < R; i++){
    random_graph_files[i] = "random_graph_" + to_string((long long int)i) + ".bin";    
    nfds[i] = open(random_graph_files[i].c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE);
    if(nfds[i] < 0){
      fprintf(stderr, "cannot open the random_graph_file: %s\n", random_graph_files[i].c_str());
      exit(EXIT_FAILURE);
    }
  }
  vector<FILE*> filefps(R);
  for(int i = 0; i < R; i++){
    filefps[i] = fopen(file_name, "r");
    if(filefps[i] == NULL){
      fprintf(stderr, "cannot open the file: %s\n", file_name);
      return 1;
    }
  }

  BinaryGraphReader<Node>** scc_bgrs = new BinaryGraphReader<Node>*[num_threads];
  for(int i = 0; i < num_threads; i++)scc_bgrs[i] = new BinaryGraphReader<Node>(1ll<<24);
  BinaryWriter<Node>** scc_bws = new BinaryWriter<Node>*[num_threads];
  for(int i = 0; i < num_threads; i++)scc_bws[i] = new BinaryWriter<Node>(1ll<<24);
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
        random_graph_bws[tid]->setfd(nfds[t]);
        fseek(filefps[tid], 0, SEEK_SET);
        Node src, dst;
        double p;
        for(; ~fscanf(filefps[tid], "%u%u%lf", &src, &dst, &p); ){
          if(src == dst)
            continue;
          if(xss[t]->nextDouble() < p){
            random_graph_bws[tid]->bwwrite(&src, 1);
            random_graph_bws[tid]->bwwrite(&dst, 1);
          }
        }
        random_graph_bws[tid]->bwflush();
        close(nfds[t]);
        runSWSSCC(random_graph_files[t].c_str(), n, comp[tid], scc_bgrs[tid], scc_bws[tid], tid);
      }
#pragma omp barrier
    }

    // tournament merge
    //double merge_start_time = gettimeofday_sec();
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
          vector< long long int > cp(n); // comp pair id
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
	}
  for(int i = 0; i < num_threads; i++)fclose(filefps[i]);
  for(int i = 0; i < R; i++)delete xss[i];
  delete [] xss;
  for(int i = 0; i < num_threads; i++){
    delete scc_bgrs[i];
    delete scc_bws[i];
  }
  delete [] scc_bgrs;
  delete [] scc_bws;
  for(int i = 0; i < num_threads; i++)delete [] comp[i];
  delete [] comp;

  // calc the number of reduced nodes
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
  FILE* fp = fopen(file_name, "r");
  if(fp == NULL){
    fprintf(stderr, "cannot open the file: %s\n", file_name);
    return false;
  }
  Node src, dst;
  double p;
  fseek(fp, 0, SEEK_SET);
  for(; ~fscanf(fp, "%u%u%lf", &src, &dst, &p); ){
    if(src == dst)
      continue;
    const Node u = mapping[src], v = mapping[dst];
    if(u != v && (rscc_size[u] > 1 || rscc_size[v] > 1)){
      long long int key = 1ll * u * n + v;
      if(sccedges.find(key) == sccedges.end()){
        sccedges[key].first = 0;
        sccedges[key].second = 1.0;
      }
      sccedges[key].first++;
      sccedges[key].second *= 1.0 - p;
    }
	}
  unordered_map< long long int, pair<int, double> >::iterator it;

  /////
  // write new graph
  const string graph_name = output_prefix + "_coarsened-graph.txt";
  const size_t vbuf_size = 1<<25;
  static char vbuf[vbuf_size];
  FILE* graph_fp = fopen(graph_name.c_str(), "w");
  if(graph_fp == NULL){
    cerr << "cannot open the file: " << graph_name << endl;
    return false;
  }
  setvbuf(graph_fp, NULL, _IOFBF, 1<<22);

  // only scc size 1 node to node
  long long int num_cedges = 0;
  fseek(fp, 0, SEEK_SET);
  for(; ~fscanf(fp, "%u%u%lf", &src, &dst, &p); ){
    if(src == dst)
      continue;
    const Node u = mapping[src], v = mapping[dst];
    if(u != v && rscc_size[u] == 1 && rscc_size[v] == 1){
      fprintf(graph_fp, "%u %u %.10e\n", u, v, p);
      num_cedges++;
    }
  }
  // otherwise
  for(it = sccedges.begin(); it != sccedges.end(); it++){
    Node u = (Node)(it->first / n);
    Node v = (Node)(it->first % n);
    double np = 1.0 - it->second.second;
    fprintf(graph_fp, "%u %u %.10e\n", u, v, np);
    num_cedges++;
  }
  fclose(graph_fp);
  printf("num_cedges=%lld\n", num_cedges);

  // write node mapping (old node ID to new node ID)
  const string mapping_name = output_prefix + "_mapping.txt";
  FILE* mapping_fp = fopen(mapping_name.c_str(), "w");
  if(mapping_fp == NULL){
    cerr << "cannot open the file: " << mapping_name << endl;
    return false;
  }
  setvbuf(mapping_fp, vbuf, _IOFBF, vbuf_size);
  fprintf(mapping_fp, "%u\n", n);
  for(Node i = 0; i < n; i++){
    fprintf(mapping_fp, "%u\n", mapping[i]);
  }
  fclose(mapping_fp);

  // finalize
  delete [] id_index;

  printf("all finished: %.5f\n", gettimeofday_sec() - start_time);
  fflush(stdout);

  return true;
}
