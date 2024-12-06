===
FP3
===

-----------------------------------
Faster Parallel Processing Pipeline
-----------------------------------

Introduction
============

The usual DIALS pipeline is to:

- import data
- find spots on all images
- index
- (optionally) determine Bravais lattice
- refine
- integrate
- (optionally) pre-scale the data
- determine symmetry
- scale

This uses all of the data for every stage in the process, to give a good global model of the experiment, and some parallelism is employed throughout the process. While it is effective it is not necessarily _quick_, and processing can take a large amount of time for more substantial data sets. The aim of the ``fp3`` pipeline is to streamline some of the analysis without significantly compromising the quality of the results, in the same spirit as the ``fast_dp`` pipeline, which uses ``XDS``.

Workflow
========

The workflow of ``fp3`` is very similar to that of ``fast_dp`` in that spots are found in three blocks at the start of the data set, indexing performed in P1 on this block, then the matrix used to integrate the entire data set in small blocks. This was found during the development of ``xia2`` to give an accurate matrix with minimal computational expense. ``fp3`` however works differently to ``fast_dp`` in the detail as the DIALS programs have a different set of assumptions. The workflow is therefore:

- import the data
- find spots on 0-5°, 45-50° and 90-95° or as close as can be reached
- use this for indexing, determine and refine a UB matrix
- for each 5° block of data:

  * find spots
  * re-index with the known UB matrix
  * refine
  * integrate

- combine all integrate data
- determine symmetry
- scale
- decide resolution limit
- (as necessary) re-scale to determined limit

This is logically *inefficient* since the spot finding on some of the data is performed twice, however it is effective as every 5° block of data may be treated independently however is guaranteed to come out with consistently indexed data (provided that the sample has not slipped excessively). As such, these steps may be performed in parallel across multiple nodes on a cluster, or across multiple cores in a many-core computer. Within each step multiple cores may also be used (e.g. spot finding and integration) allowing for a substantial amount of compute to be applied to the problem. For finely sliced data the spot finding and integration are the most computationally expensive steps, and here these are mitigated in the first instance by limiting the range of data used for determination of the orientation matrix, and in the second by application of substantial parallelism.

Usage
=====

The basic program usage is::

  dials.fp3 \
    max_workers=288 \
    parallelism=drmaa \
    worker_nproc=4 \
    ../Insulin15/Insulin15_3_master.h5

Here we have allowed up to 288 5° blocks to be processed in parallel, each using up to 4 cores (so a total of  up to 1152 cores) using the ``drmaa`` mechanism for cluster job submission. Alternatively within a single machine the following may be used::

  dials.fp3 max_workers=8 worker_nproc=2 ../data_master.h5

which will use up to 16 cores for the processing. For every stage in the processing (spot finding, indexing, refinement, integration &c.) the usual PHIL processing parameters may be passed in as every PHIL scope in DIALS is included in the ``fp3`` scope - for example assigning the unit cell and symmetry in indexing is as easy as setting e.g.::

  known_symmetry.unit_cell=57,57,150,90,90,90 \
  known_symmetry.space_group=P41212

for example, for thaumatin. Clearly including the full PHIL scope for every step provides a *lot* of options - some of which could of course be mutually inconsistent, however in the most general case no options should be necessary.

In the event that e.g. you want to rerun the scaling step with different parameters all that is necessary is to remove the ``scale`` directory and rerun - the results recorded to this point are re-loaded so the processing will pick up at the point of scaling.

Acknowledgements
================

Clearly ``fp3`` would be impossible without the DIALS toolkit. The work has been supported by Diamond Light Source as part of ongoing operations.

License
=======

Copyright (c) 2012-2020 Diamond Light Source.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the
  distribution.
- Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Authors
=======

- Graeme Winter
- Irakli Sikharulidze
