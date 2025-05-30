from __future__ import division
from __future__ import print_function
import wx
import numpy as np

# set backend before importing pyplot
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def GetBitmap_from_np_array(data2d):
    lc_fig = plt.figure(frameon=False)

    xmax = len(data2d[:, 1])
    ymax = len(data2d[1, :])
    vl_max = np.amax(data2d)
    vl_min = np.amin(data2d)
    d = vl_max - vl_min
    vl_mid_low = vl_min + d / 3.0
    vl_mid_hig = vl_max - d / 3.0

    lc_fig.set_size_inches(xmax * 0.6, ymax * 0.6)
    ax = plt.Axes(lc_fig, [0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    lc_fig.add_axes(ax)
    plt.imshow(np.transpose(data2d), interpolation="nearest", cmap="hot")

    print("xmax =", xmax)
    print("ymax =", ymax)
    for xpos in range(xmax):
        for ypos in range(ymax):
            f_num = data2d[xpos, ypos]
            g = float("{0:.2f}".format(float(f_num)))

            txt_dat = str(g)
            if g < vl_mid_low:
                clr_chr = "yellow"
            elif g > vl_mid_hig:
                clr_chr = "black"
            else:
                clr_chr = "blue"
            plt.annotate(
                txt_dat,
                xy=(xpos - 0.3, ypos + 0.3),
                xycoords="data",
                color=clr_chr,
                size=12.0,
            )

    lc_fig.canvas.draw()
    width, height = lc_fig.canvas.get_width_height()
    np_buf = np.fromstring(lc_fig.canvas.tostring_rgb(), dtype=np.uint8)
    np_buf.shape = (width, height, 3)
    np_buf = np.roll(np_buf, 3, axis=2)
    wx_image = wx.EmptyImage(width, height)
    wx_image.SetData(np_buf)
    wxBitmap = wx_image.ConvertToBitmap()

    plt.close(lc_fig)

    return wxBitmap


def build_np_img(width=64, height=64):
    data2d = np.zeros((width, height), "float")
    print("width, height =", width, height)
    tot_lng = np.sqrt(width**2.0 + height**2.0)
    for x in range(0, width):
        for y in range(0, height):
            dx = float(width - x)
            dy = float(height - y)
            data2d[x, y] = tot_lng - np.sqrt(dx**2.0 + dy**2.0)

    print("data2d.max =", data2d.max())
    return data2d


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, title="Bitmaps")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True


class MyFrame(wx.Frame):
    def __init__(
        self,
        parent,
        id=wx.ID_ANY,
        title="",
        pos=wx.DefaultPosition,
        size=wx.DefaultSize,
        style=wx.DEFAULT_FRAME_STYLE,
        name="MyFrame",
    ):
        super(MyFrame, self).__init__(parent, id, title, pos, size, style, name)
        # Attributes
        self.panel = wx.Panel(self)

        data2d = build_np_img(width=8, height=4)
        print(data2d)
        bitmap = GetBitmap_from_np_array(data2d)

        self.bitmap = wx.StaticBitmap(self.panel, bitmap=bitmap)


if __name__ == "__main__":
    app = MyApp(redirect=False)
    app.MainLoop()
