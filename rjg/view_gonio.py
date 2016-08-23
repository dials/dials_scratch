# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import division
from gltbx import wx_viewer
import copy
import wx
import wxtbx.utils
from gltbx.gl import *
import gltbx
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
import libtbx.phil

help_message = '''

'''

phil_scope= libtbx.phil.parse("""
  angle = None
    .type = float
    .multiple = True
  z_min = None
    .type = float
  z_max = None
    .type = float
  marker_size = 100
    .type = int(value_min=1)
  autospin = False
    .type = bool
""")

def settings():
  return phil_scope.fetch().extract()

class render_3d(object):

  def __init__(self):
    self.reflections = None
    self.goniometer_orig = None

  def load_imageset(self, imageset):
    self.imageset = imageset
    if self.imageset.get_goniometer() is not None:
      from scitbx import matrix
      gonio = self.imageset.get_goniometer()
      setting_rotation = matrix.sqr(gonio.get_setting_rotation())
      axis = setting_rotation * matrix.col(gonio.get_rotation_axis()).elems
      self.viewer.set_rotation_axis(axis)
    self.viewer.set_beam_vector(self.imageset.get_beam().get_s0())

    detector = self.imageset.get_detector()
    beam = self.imageset.get_beam()
    self.set_points()

  def set_goniometer_points(self):
    gonio = self.imageset.get_goniometer()
    scan = self.imageset.get_scan()
    detector = self.imageset.get_detector()

    angle = self.settings.angle
    if angle:
      assert len(angle) == len(gonio.get_angles())
      gonio.set_angles(angle)

    gonio_masker = self.imageset.reader().get_format().get_goniometer_shadow_masker(gonio)
    points = gonio_masker.extrema_at_scan_angle(
      gonio.get_angles()[gonio.get_scan_axis()])

    line_i_seqs = flex.vec2_double(((0,i) for i in range(1, points.size())))
    line_i_seqs += (self.viewer.points.size(), self.viewer.points.size())
    for i_seqs in line_i_seqs:
      self.viewer.line_i_seqs.append([int(i_seq) for i_seq in i_seqs])

    self.viewer.points.extend(points)

    shadow = gonio_masker.project_extrema(
      detector, gonio.get_angles()[gonio.get_scan_axis()])

    for shadow_points, p in zip(shadow, detector):
      n = self.viewer.points.size()
      line_i_seqs = []
      line_colors = {}
      for i in range(shadow_points.size()):
        if i < shadow_points.size() - 1:
          line_i_seqs.append((n+i, n+i+1))
        else:
          line_i_seqs.append((n+i, n))
        line_colors[line_i_seqs[-1]] = (1,1,1)
      self.viewer.line_colors.update(line_colors)
      self.viewer.line_i_seqs.extend(line_i_seqs)
      self.viewer.points.extend(p.get_lab_coord(shadow_points *p.get_pixel_size()[0]))

    #self.viewer.points.extend(projected_points)

  def set_detector_points(self):
    detector = self.imageset.get_detector()
    points = flex.vec3_double()
    line_i_seqs = flex.vec2_double()
    i = 0
    for p in detector:
      image_size = p.get_image_size_mm()
      points.append(p.get_lab_coord((0,0)))
      points.append(p.get_lab_coord((0,image_size[1])))
      points.append(p.get_lab_coord(image_size))
      points.append(p.get_lab_coord((image_size[0],0)))
      line_i_seqs.append((i, i+1))
      line_i_seqs.append((i+1, i+2))
      line_i_seqs.append((i+2, i+3))
      line_i_seqs.append((i+3, i))
      i += 4
    line_i_seqs += (self.viewer.points.size(), self.viewer.points.size())
    self.viewer.points.extend(points)
    for i_seqs in line_i_seqs:
      self.viewer.line_i_seqs.append([int(i_seq) for i_seq in i_seqs])

  def set_points(self):
    # reset point/line lists
    self.viewer.reset()

    try:
      self.set_goniometer_points()
    except Exception:
      pass
    self.set_detector_points()
    self.viewer.update_minimum_covering_sphere()


class ExperimentViewer(wx.Frame, render_3d):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    render_3d.__init__(self)
    self.parent = self.GetParent()
    self.statusbar = self.CreateStatusBar()
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)

    app = wx.GetApp()
    if getattr(app, "settings", None) is not None:
      # XXX copying the initial settings avoids awkward interactions when
      # multiple viewer windows are opened
      self.settings = copy.deepcopy(app.settings)
    else :
      self.settings = settings()

    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND|wx.ALL)
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.Bind(wx.EVT_ACTIVATE, self.OnActive)
    self.viewer.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
    self.viewer.SetFocus()

  def load_imageset(self, imageset):
    render_3d.load_imageset(self, imageset)
    self.settings_panel.add_goniometer_controls(imageset.get_goniometer())

  def OnActive (self, event) :
    if self.IsShown() and type(self.viewer).__name__ != "_wxPyDeadObject":
      self.viewer.Refresh()

  def OnClose (self, event) :
    self.Unbind(wx.EVT_ACTIVATE)
    self.Destroy()
    event.Skip()

  def OnDestroy (self, event) :
    if self.parent is not None:
      self.parent.viewer = None
    event.Skip()

  def OnKeyDown(self, event):
    key = event.GetUnicodeKey()
    if key == wx.WXK_NONE:
      key = event.GetKeyCode()
    dxs = {wx.WXK_LEFT:-1,
           wx.WXK_RIGHT:+1,
           wx.WXK_UP:0,
           wx.WXK_DOWN:0}
    dys = {wx.WXK_LEFT:0,
           wx.WXK_RIGHT:0,
           wx.WXK_UP:+1,
           wx.WXK_DOWN:-1}

    if key in dxs:
      dx = dxs[key]
      dy = dys[key]
      if event.ShiftDown():
        scale = 0.1
      else:
        scale = 1.0
      self.do_Step(dx, dy, scale)

  def do_Step(self, dx, dy, scale):
    v = self.viewer
    rc = v.rotation_center
    glMatrixMode(GL_MODELVIEW)
    gltbx.util.rotate_object_about_eye_x_and_y(
      scale, rc[0], rc[1], rc[2],
      dx, dy, 0, 0)
    v.OnRedraw()

  def create_viewer_panel (self) :
    self.viewer = RLVWindow(settings=self.settings, parent=self, size=(800,600),
      #orthographic=True
      )

  def create_settings_panel (self) :
    self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

  def set_points(self):
    render_3d.set_points(self)

  def update_settings(self, *args, **kwds):
    self.set_points()
    self.viewer.update_settings(*args, **kwds)


class settings_window (wxtbx.utils.SettingsPanel) :
  def __init__ (self, *args, **kwds) :
    wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_CHAR, self.OnChar)

  def OnChar (self, event) :
    self.GetParent().viewer.OnChar(event)

  def add_controls (self) :
    pass
    # d_min control
    #from wx.lib.agw import floatspin

  def add_goniometer_controls(self, goniometer):
    from wx.lib.agw import floatspin

    self.rotation_angle_ctrls = []
    names = goniometer.get_names()
    axes = goniometer.get_axes()
    angles = goniometer.get_angles()
    for name, axis, angle in zip(names, axes, angles):
      ctrl = floatspin.FloatSpin(parent=self, increment=1, digits=1)
      ctrl.SetValue(angle)
      ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
      if wx.VERSION >= (2,9): # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
        ctrl.SetBackgroundColour(self.GetBackgroundColour())
      box = wx.BoxSizer(wx.HORIZONTAL)
      self.panel_sizer.Add(box)
      label = wx.StaticText(self,-1,"%s angle" %name)
      box.Add(ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeSettings, ctrl)
      self.rotation_angle_ctrls.append(ctrl)

  def OnChangeSettings(self, event):
    for i, ctrl in enumerate(self.rotation_angle_ctrls):
      self.settings.angle[i] = ctrl.GetValue()
    self.parent.update_settings()


class RLVWindow(wx_viewer.show_points_and_lines_mixin):

  def __init__(self, settings, *args, **kwds):
    super(RLVWindow, self).__init__(*args, **kwds)
    self.settings = settings
    self.points = flex.vec3_double()
    self.colors = None
    self.rotation_axis = None
    self.beam_vector = None
    self.flag_show_minimum_covering_sphere = False
    self.minimum_covering_sphere = None
    self.field_of_view_y = 20
    if self.settings.autospin:
      self.autospin_allowed = True
      self.yspin = 1
      self.xspin = 1
      self.autospin = True

  def reset(self):
    self.labels = []
    self.points = flex.vec3_double()
    self.line_i_seqs = []
    self.line_colors = {}
    self.spheres = []
    self.labels_display_list = None
    self.points_display_list = None
    self.lines_display_list = None

  def set_points(self, points):
    self.points = points
    self.points_display_list = None
    if self.minimum_covering_sphere is None:
      self.update_minimum_covering_sphere()

  def set_colors(self, colors):
    assert len(colors) == len(self.points)
    self.colors = colors

  def draw_points(self):
    if self.points_display_list is None:
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      glLineWidth(1)
      if self.colors is None:
        self.colors = flex.vec3_double(len(self.points), (1,1,1))
      for point, color in zip(self.points, self.colors):
        self.draw_cross_at(point, color=color)
      self.points_display_list.end()
    self.points_display_list.call()

  def set_rotation_axis(self, axis):
    self.rotation_axis = axis

  def set_beam_vector(self, beam):
    self.beam_vector = beam

  #--- user input and settings
  def update_settings (self) :
    self.points_display_list = None
    self.Refresh()

  def update_minimum_covering_sphere(self):
    self.minimum_covering_sphere = minimum_covering_sphere(
      self.points, epsilon=1e-3)

  def draw_cross_at(self, (x,y,z), color=(1,1,1), f=None):
    if f is None:
      f = 0.01 * self.settings.marker_size
    wx_viewer.show_points_and_lines_mixin.draw_cross_at(
      self, (x,y,z), color=color, f=f)

  def DrawGL(self):
    wx_viewer.show_points_and_lines_mixin.DrawGL(self)
    gonio = self.parent.imageset.get_goniometer()
    beam = self.parent.imageset.get_beam()
    from scitbx import matrix
    R = matrix.identity(3)
    names = reversed(gonio.get_names())
    axes = reversed(gonio.get_axes())
    angles = reversed(gonio.get_angles())
    for name, axis, angle in zip(names, axes, angles):
      axis = R * matrix.col(axis)
      self.draw_axis(axis.elems, name)
      R = axis.axis_and_angle_as_r3_rotation_matrix(angle, deg=True) * R
    self.draw_axis(beam.get_s0(), "beam")

  def draw_axis(self, axis, label):
    if self.minimum_covering_sphere is None:
      self.update_minimum_covering_sphere()
    s = self.minimum_covering_sphere
    scale = max(max(s.box_max()), abs(min(s.box_min())))
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    glDisable(GL_LIGHTING)
    glColor3f(1.0, 1.0, 1.0)
    glLineWidth(1.0)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(axis[0]*scale, axis[1]*scale, axis[2]*scale)
    glEnd()
    glRasterPos3f(0.5+axis[0]*scale, 0.2+axis[1]*scale, 0.2+axis[2]*scale)
    gltbx.fonts.ucs_bitmap_8x13.render_string(label)
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-axis[0]*scale, -axis[1]*scale, -axis[2]*scale)
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  def rotate_view(self, x1, y1, x2, y2, shift_down=False, scale=0.1):
    super(RLVWindow, self).rotate_view(
      x1, y1, x2, y2, shift_down=shift_down, scale=scale)

  def OnLeftUp(self,event):
    self.was_dragged = True
    super(RLVWindow, self).OnLeftUp(event)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_experiments=True,
    #read_reflections=True,
    check_format=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)

  if (len(datablocks) == 0 and len(experiments) == 0):
    parser.print_help()
    exit(0)

  if len(datablocks) == 0 and len(experiments) > 0:
    imagesets = experiments.imagesets()
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  assert len(imagesets) == 1
  imageset = imagesets[0]
  gonio = imageset.get_goniometer()
  if params.angle:
    assert len(params.angle) == len(gonio.get_angles())
  else:
    for angle in gonio.get_angles():
      params.angle.append(angle)

  import wxtbx.app
  a = wxtbx.app.CCTBXApp(0)
  a.settings = params
  f = ExperimentViewer(
    None, -1, "Experiment viewer", size=(1024,768))
  f.load_imageset(imageset)
  f.Show()
  a.SetTopWindow(f)
  #a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
