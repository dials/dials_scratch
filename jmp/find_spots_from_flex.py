from __future__ import print_function

from dials.array_family import flex

def find_spots(image, mask, min_spot_size=1, max_spot_size=1000):
  from dials.algorithms.spot_finding.threshold import XDSThresholdStrategy
  from dials.model.data import PixelList
  from dials.model.data import PixelListLabeller
  threshold_image = XDSThresholdStrategy()

  threshold_mask = threshold_image(image, mask)
  plist = PixelList(0, image, threshold_mask)

  pixel_labeller = PixelListLabeller()

  pixel_labeller.add(plist)

  creator = flex.PixelListShoeboxCreator(
      pixel_labeller,
      0,                   # panel
      0,                   # zrange
      True,                # twod
      min_spot_size,  # min_pixels
      max_spot_size,  # max_pixels
      False)
  shoeboxes = creator.result()

  centroid = shoeboxes.centroid_valid()
  intensity = shoeboxes.summed_intensity()
  observed = flex.observation(shoeboxes.panels(), centroid, intensity)

  return flex.reflection_table(observed, shoeboxes)



if __name__ == '__main__':
  import sys
  from dxtbx.datablock import DataBlockFactory
  db = DataBlockFactory.from_json_file(sys.argv[1])
  iset = db[0].extract_imagesets()[0]
  image = iset.get_raw_data(0)[0]
  mask = iset.get_mask(0)[0]

  reflections = find_spots(image, mask)

  print(len(reflections))
