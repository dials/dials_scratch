
def show_image_stack(data, axis=0, interpolation='none', **kwargs):
  '''
  Display a stack of images

  '''
  import matplotlib.pyplot as plt
  from matplotlib.widgets import Slider, Button, RadioButtons

  # check dim
  assert(data.ndim == 3)

  # generate figure
  fig = plt.figure()
  ax = plt.subplot(111)
  fig.subplots_adjust(left=0.25, bottom=0.25)

  # select first image
  s = [slice(0, 1) if i == axis else slice(None) for i in xrange(3)]
  im = data[s].squeeze()

  # display image
  l = ax.imshow(im, interpolation=interpolation, **kwargs)

  # define slider
  ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg='lightgoldenrodyellow')

  slider = Slider(ax, 'Axis %i index' % axis, 0, data.shape[axis] - 1,
                  valinit=0, valfmt='%i')

  def update(val):
    ind = int(slider.val)
    s = [slice(ind, ind + 1) if i == axis else slice(None)
             for i in xrange(3)]
    im = data[s].squeeze()
    l.set_data(im)
    fig.canvas.draw()

  slider.on_changed(update)

  plt.show()


def show_image_stack_multi_view(data, interpolation='none',
                                axis_names=["Axis 0", "Axis 1", "Axis 2"], **kwargs):
  '''
  Display a stack of images

  '''
  import matplotlib.pyplot as plt
  from matplotlib.widgets import Slider, Button, RadioButtons

  # check dim
  assert(data.ndim == 3)

  # generate figure
  fig = plt.figure()
  ax1 = plt.subplot(131)
  ax2 = plt.subplot(132)
  ax3 = plt.subplot(133)
  fig.subplots_adjust(left=0.25, bottom=0.25)

  # select first image
  s1 = [slice(0, 1), slice(None), slice(None)]
  s2 = [slice(None), slice(0, 1), slice(None)]
  s3 = [slice(None), slice(None), slice(0, 1)]
  im1 = data[s1].squeeze()
  im2 = data[s2].squeeze()
  im3 = data[s3].squeeze().transpose()

  ax1.set_ylabel(axis_names[1])
  ax1.set_xlabel(axis_names[2])

  ax2.set_ylabel(axis_names[0])
  ax2.set_xlabel(axis_names[2])

  ax3.set_ylabel(axis_names[1])
  ax3.set_xlabel(axis_names[0])

  # display image
  l = [None] * 3
  l[0] = ax1.imshow(im1, interpolation=interpolation, **kwargs)
  l[1] = ax2.imshow(im2, interpolation=interpolation, **kwargs)
  l[2] = ax3.imshow(im3, interpolation=interpolation, **kwargs)

  # define slider
  ax4 = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg='lightgoldenrodyellow')
  ax5 = fig.add_axes([0.25, 0.06, 0.65, 0.03], axisbg='lightgoldenrodyellow')
  ax6 = fig.add_axes([0.25, 0.02, 0.65, 0.03], axisbg='lightgoldenrodyellow')

  slider1 = Slider(ax4, '%s index' % axis_names[0], 0, data.shape[0] - 1, valinit=0, valfmt='%i')
  slider2 = Slider(ax5, '%s index' % axis_names[1], 0, data.shape[1] - 1, valinit=0, valfmt='%i')
  slider3 = Slider(ax6, '%s index' % axis_names[2], 0, data.shape[2] - 1, valinit=0, valfmt='%i')

  class Update:
    def __init__(self, axis):
      self.axis = axis
    def __call__(self, val):
      ind = int(val)
      s = [slice(ind, ind + 1) if i == self.axis else slice(None)
               for i in xrange(3)]
      im = data[s].squeeze()
      if self.axis == 2:
        im = im.transpose()
      l[self.axis].set_data(im)
      fig.canvas.draw()

  slider1.on_changed(Update(0))
  slider2.on_changed(Update(1))
  slider3.on_changed(Update(2))

  plt.show()

if __name__ == '__main__':

  from dials.array_family import flex
  from math import exp

  data = flex.double(flex.grid(15,15,15))
  for k in range(15):
    for j in range(15):
      for i in range(15):
        A = exp(-(k-7)**2 / (2*1**2))
        B = exp(-(j-7)**2 / (2*2**2))
        C = exp(-(i-7)**2 / (2*3**2))
        data[k,j,i] = A*B*C


  show_image_stack(data.as_numpy_array(), vmax=max(data))
  show_image_stack_multi_view(data.as_numpy_array(), vmax=max(data))
