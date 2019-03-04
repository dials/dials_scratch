#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Import a MTZ file as an ExperimentList and reflection_table as per DIALS.

For information on the MTZ format, see
http://www.ccp4.ac.uk/html/mtzformat.html

The DIALS workflow exploits two key object types — 
``dials_array_family_flex_ext.reflection_table`` to represent the reflection
data, and ``dxtbx_model_ext.ExperimentList``, inherited from the CCTBX
Project's dxtbx (diffraction image toolbox), to represent the experiment model.
The former is defined in
``.../dials/array_family/boost_python/flex_reflection_table.cc`` and the latter
in ``.../../cctbx_project/dxtbx/model/experiment_list.h``.  This module allows
one to populate a reflection_table object and an ExperimentList object with
data from a MTZ file, whether merged or unmerged.
"""

from __future__ import absolute_import, division, print_function
