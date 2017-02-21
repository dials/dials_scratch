from dials_scratch.lbl_feb_2017.apple import Apple
import sys

apple = Apple(sys.argv[1], sys.argv[2])
spots = apple.find_spots()
indexed = apple.index(spots)
print 'Spots found / indexed:', spots.size(), indexed.size()
apple.refine(do_print=True)
hklout = sys.argv[3]
distance_map = apple.render_distance()

mask = apple.get_signal_mask()
background = apple.make_background()
spot = apple.get_background_subtracted_spots()

apple.plot_log_map(spot, 'spot.png')
apple.plot_map(distance_map, 'distance.png')
apple.plot_log_map(background, 'background.png')
apple.plot_map(mask, 'mask.png')

reflections = apple.integrate()
reflections.as_pickle(hklout)
print 'Spots integrated: ', reflections.size()
