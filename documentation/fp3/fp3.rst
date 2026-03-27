FP3: Faster Parallel Processing Pipeline
========================================

Introduction
------------

``dials.fp3`` is a very simple processing pipeline following the behaviour of ``fast_dp`` - except here using ``dials`` tools for the processing. The workflow is to:

1. Import data
2. Find spots on three 5° wedges of data starting at 0°, 45°, 90° or
   as close as can be reached to there
3. Index these spots, without making any effort to determine the
   symmetry
4. Split the data into 5° blocks, and for each block:
   i. slice down the input experiment to just the block, copy the crystal model from step 3. into the imported experiment
   ii. find spots on the block
   iii. index the data against the known UB matrix
   iv. refine (scan static + scan varying)
   v. integrate
5. Combine the data from all integrated runs
6. Determine the symmetry
7. Scale (no resolution limit)
8. Determine high resolution limit
9. Scale output of step 7 with resolution from step 8

   
