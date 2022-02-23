# Sciara-fv2 CUDA
[![License: GPL](https://img.shields.io/badge/License-GPL-blue.svg)](/LICENSE)  
*SCIARA is a deterministic Macroscopic Cellular Automata model for simulating lava flows.*
  
In this paper we presents a CUDA parallelization of the deterministic Macroscopic Cellular Automata model SCIARA for simulating lavaflows. Furthermore, we present the results of an experimental analysiswe conducted in order to evaluate performance-wise impact. Specifically, given benchmarks of our parallel implementation with respect toserial algorithm, we observe how parallel implementation is very com-petitive and even faster.

## :zap: How to run?

To run Sciara-fv2 CUDA backend:
```shell script
$ make run-$VERSION
# Where $VERSION can be 'cuda-straightforward', 'cuda-tiled-halo', 'cuda-tiled-no-halo', 'cuda-multi-gpu'
```

**NOTE:** Building and running requires [CUDA](https://developer.nvidia.com/cuda-zone) compatible devices with 5.2 or higher compute capabilities.


## :fire: Performance

<table style="border: 0 !important">
  <tr>
    <td><img src="/data/speedups/bar.png" alt="" /></td>
    <td><img src="/data/speedups/line.png" alt="" /></td>
  </tr>
</table>

## :camera: Screenshots

<table style="border: 0 !important">
  <tr>
    <td><img src="/data/plots/flow.png" alt="" /></td>
    <td><img src="/data/plots/temp.png" alt="" /></td>
    <td><img src="/data/plots/morph.png" alt="" /></td>
  </tr>
</table>

## :globe_with_meridians: Third-Party Software
Sciara-fv2 CUDA uses and depends on third-party open-source tools and libraries which are outside of this repository.

## :page_with_curl: License

Copyright (c) Ransome CA. All rights reserved.

Licensed under the [GPL-3.0](/LICENSE) license.
