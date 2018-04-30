from dials.array_family import flex
import cPickle as pickle

def fix_xy(reflections_in, reflections_out):
    reflections = pickle.load(open(reflections_in, 'r'))

    # validate that input is from 60 panel detector, assumed to be
    # 5 x 12 configuration

    assert flex.max(reflections['panel']) == 59

    delta_x = 487 + 7
    delta_y = 195 + 17

    panel_x = reflections['panel'].as_int() % 5
    panel_y = reflections['panel'].as_int() / 5

    x_offset = delta_x * panel_x
    y_offset = delta_y * panel_y

    # apply fixes
    x, y, z = reflections['xyzobs.px.value'].parts()
    x += x_offset.as_double()
    y += y_offset.as_double()
    reflections['xyzobs.px.value'] = flex.vec3_double(x, y, z)

    x, y, z = reflections['xyzcal.px'].parts()
    x += x_offset.as_double()
    y += y_offset.as_double()
    reflections['xyzcal.px'] = flex.vec3_double(x, y, z)

    # save - should probably do this "properly"
    pickle.dump(reflections, open(reflections_out, 'wb'))

    return

if __name__ == '__main__':
    import sys
    fix_xy(sys.argv[1], sys.argv[2])
