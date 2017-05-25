# Coarsening (Massive Influence Graph)

Coarsening (reduction) algorithms with linear-space and sub-linear-space for influence graphs under the independent cascade model


## Usage
* Given a graph with edge probabilities, it generates a coarsened graph and vertex mapping (vetex ID in the original graph to vertex ID in the coarsened graph).
* OpenMP is necessary.

    $ cd linear-space/sublinear-space
    $ make
    $ ./main graph R output_prefix num_threads

* graph: input file (see below)
* R: the number of random graphs (larger = more accurate)
* output_prefix: resulting output file name prefix
    * "[output_prefix]_coarsened-graph.txt" is the coarsened graph
    * "[output_prefix]_mapping.txt" is the correspondence mapping π : V -> W (see algorithm 1 and 2 in our paper)
* num_threads: the number of running threads


### Example

    $ make
    $ ./main ../sample_graphs/sample_graph0.tsv 16 ./output 2


### Format of the input graph and output coarsened graph

    u_1	v_1	p_1
    ...
    u_i	v_i	p_i
    ...
    u_m	v_m	p_m

* The i-th line means an directed edge (u_i, v_i) with propagation probability p_i (see `sample_graph.tsv`).
* Vertices should be described by integers starting from zero.


### Format of the output mapping file
    n(the number of vertices in the original input graph)
    w_1
    ...
    w_n
* The first line stands for the number of vertices in the input graph.
* Next n lines stand for the correspondence mapping π.


## References

* Naoto Ohsaka, Tomohiro Sonobe, Sumio Fujita, and Ken-ichi Kawarabayashi. **[Coarsening Massive Influence Networks for Scalable Diffusion Analysis](http://dl.acm.org/citation.cfm?id=3064045)**.
In SIGMOD, pages 635--650, 2017.

