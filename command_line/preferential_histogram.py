from __future__ import division, print_function

def multihist(unmerged_mtz):
    from iotbx.reflection_file_reader import any_reflection_file

    reader = any_reflection_file(unmerged_mtz)
    assert reader.file_type() == 'ccp4_mtz'
    mtz_object = reader.file_content()
    arrays = reader.as_miller_arrays(merge_equivalents=False)

    for ma in arrays:
      if ma.info().labels == ['I', 'SIGI']:
        intensities = ma
      elif ma.info().labels == ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']:
        intensities = ma

    indices = mtz_object.extract_original_index_miller_indices()
    intensities = intensities.customized_copy(
      indices=indices, info=intensities.info())

    intensities = intensities.resolution_filter(d_min=2.0)

    merging = intensities.merge_equivalents()
    multiplicities = merging.redundancies().complete_array(new_data_value=0)
    mult_acentric = multiplicities.select_acentric().data()
    mult_centric = multiplicities.select_centric().data()

if __name__ == '__main__':
    import sys
    multihist(sys.argv[1])
