# IITM Global Router

A GPU-accelerated, congestion-aware 3D global router. Instead of starting from 2D FLUTE Steiner trees, it directly constructs a 3D rectilinear minimum spanning tree (3D-RMST) for every multi-pin net using a parallel GPU implementation of Borůvka's algorithm with congestion-dependent edge costs. All nets are routed in a single pass, without rip-up and reroute.

## Overview

For each batch of nets the router performs:

1. **Batch generation** (`part_net/`) – nets are sorted by bounding-box area and grouped into batches of nets with non-overlapping bounding boxes (Boost.Geometry), so every net in a batch can be routed in parallel. Batch selection is itself parallelized by splitting nets into x-interval buckets.
2. **3D L-pattern generation** (`l_pattern_files/`) – for every pin pair of every net in the batch, all candidate 3D L-shaped routes are costed on the GPU and the cheapest is kept. The congestion cost of a GCELL is `exp(d + 1 - c)` (`d` = demand, `c` = capacity).
3. **Parallel 3D-RMST construction** (`brvsk_files/`, `mst_files/`) – Borůvka's algorithm on the GPU, with a disjoint-set-union (DSU) structure, picks `n-1` of the `n(n-1)/2` candidate edges of each net.
4. **Parallel edge commit** (`mst_files/mst_bit_commit.hpp`, `bit_util.hpp`) – redundant GCELLS from overlapping L-routes are removed, and GCELL capacities and congestion costs are updated with binary indexed trees (Fenwick trees) in a single kernel using atomic updates.
5. **Route writing** (`writer_files/`) – global routes are written out on CPU threads while the next batches are routed (pipelined execution).

## Requirements

Tested configuration:

| Component | Version |
|---|---|
| OS | Ubuntu 22.04 (64-bit Linux) |
| GPU | NVIDIA RTX 4080 (compute capability 8.9) |
| CUDA toolkit | 13.1 (includes Thrust) |
| Host compiler | g++ 11.4, C++20, OpenMP |
| Boost | Boost.Geometry (header-only) |
| spdlog | header-only |

Boost and spdlog headers are expected under `/usr/local/include`. Adjust the `-I` path in `build.sh` if they are installed elsewhere.

## Build

```bash
./build.sh
```

which runs

```bash
nvcc main.cpp -I/usr/local/include -std=c++20 -x cu -arch=sm_80 -O3 -Xcompiler='-fopenmp' -o route
```

and produces the `route` binary. For a different GPU you may set `-arch` to its compute capability (e.g. `-arch=sm_89` for RTX 40-series, `-arch=sm_90` for H100).

## Benchmarks

The router reads the [ISPD 2025 contest](https://github.com/liangrj2014/ISPD25_contest) benchmarks (`.cap` and `.net` files):

```
ISPD2025_benchmarks/
├── visible/<design>/<design>.cap, <design>.net
└── blind/<design>/<design>.cap, <design>.net
```

Designs: `ariane`, `bsg_chip`, `mempool_tile_wrap`, `mempool_group`, `mempool_cluster`, `NV_NVDLA_partition_c`.

## Usage

Configuration used for the results in the paper:

```bash
BENCH=~/ISPD2025_benchmarks
./route -cap $BENCH/visible/ariane/ariane.cap \
        -net $BENCH/visible/ariane/ariane.net \
        -new_fmt -chip_name ariane_visible -gr_out ariane_visible.route \
        -mst_cpu_threads 10 -bit_pmodel -batch_part -num_xdim_parts 4 \
        -wrt_mrg -rseg_fmt -oprd_ppatch
```

For the other designs, replace `visible/ariane/ariane` with the corresponding `visible/<design>/<design>` or `blind/<design>/<design>` paths and set `-chip_name` accordingly.

The number of x-interval buckets (`-num_xdim_parts`) sets the parallelism of batch generation: more buckets give more parallelism and a lower runtime, especially on large designs. Its value was chosen per design in the paper. For example, for the large `mempool_cluster` design:

```bash
./route -cap $BENCH/visible/mempool_cluster/mempool_cluster.cap \
        -net $BENCH/visible/mempool_cluster/mempool_cluster.net \
        -new_fmt -chip_name mempool_cluster_visible -gr_out mempool_cluster_visible.route \
        -mst_cpu_threads 10 -bit_pmodel -batch_part -num_xdim_parts 15 \
        -wrt_mrg -rseg_fmt -oprd_ppatch
```

### Options

| Option | Description |
|---|---|
| `-cap <file>` | Routing resource (`.cap`) file |
| `-net <file>` | Net (`.net`) file |
| `-gr_out <file>` | Output global-route file |
| `-chip_name <name>` | Design name (used for logging) |
| `-new_fmt` | Parse the ISPD25 `.cap` and `.net` format |
| `-bit_pmodel` | Use binary indexed trees for congestion-cost updates and queries |
| `-batch_part` | Parallel batch generation using x-interval buckets |
| `-num_xdim_parts <k>` | Number of x-interval buckets for parallel batch generation (default 4) |
| `-mst_cpu_threads <t>` | Number of CPU threads for MST-related host work |
| `-writer_threads <t>` | Number of route-writer threads |
| `-wrt_mrg` | Merge collinear route segments when writing |
| `-rseg_fmt` | Write routes in the ISPD25 guide format (DBU coordinates and metal layer names) |

## Output

The output file contains the global route of each net in the ISPD25 global-routing format. It can be evaluated with the ISPD25 contest evaluator and the OpenROAD-based timing and power scripts, or converted into routing corridors for a detailed router.

## Evaluation in the paper

- Detailed routing, DRC and open-net checks: Synopsys IC Compiler II (ICC2), with the global routes given as routing corridors.
- TNS, WNS and total power: the modified OpenROAD and evaluation scripts from the ISPD25 contest.
- Baseline: InstantGR (Xiao et al., IEEE TCAD 2026) at commit `cfce09f`.

## License

Copyright 2026 Sammith K A and S Ramprasath.

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for details.

## Contact

- Sammith K A – ee24s008@smail.iitm.ac.in
- S Ramprasath – ramprasath@ee.iitm.ac.in

Department of Electrical Engineering, Indian Institute of Technology Madras, Chennai, India.
