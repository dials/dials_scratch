from __future__ import annotations

from math import pi
import numpy as np

import libtbx.phil

from libtbx import Auto
from scitbx.array_family import flex
from scitbx.math import minimum_covering_sphere

from dials_scratch.dgw.rlv import Render3d
from dxtbx import flumpy
from collections import namedtuple
from napari.experimental import link_layers

from magicgui import magicgui
import napari

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.reciprocal_lattice.phil_scope

show_rotation_axis = False
  .type = bool
show_beam_vector = False
  .type = bool
show_reciprocal_cell = True
  .type = bool
label_nearest_point = False
  .type = bool
marker_size = Auto
  .type = int(value_min=1)
autospin = False
  .type = bool
model_view_matrix = None
  .type = floats(size=16)
""",
    process_includes=True,
)


@magicgui(auto_call=True, d_min={"label": "high resolution", "step": 0.05})
def rlv_display(
    viewer: napari.Viewer, marker_size: int, d_min: float, z_min: int, z_max: int
):
    for layer in viewer.layers:
        if layer.name.startswith("relps"):
            layer.size[:] = marker_size
            shown = layer.properties["res"] >= d_min
            shown = shown & (layer.properties["z"] >= z_min)
            shown = shown & (layer.properties["z"] <= z_max)
            layer.shown = shown
            layer.refresh()


class ReciprocalLatticeViewer(Render3d):
    def __init__(self, parent, id, title, size, settings=None, *args, **kwds):
        Render3d.__init__(self, settings=settings)

        self.viewer = RLVWindow(
            settings=self.settings,
        )

    def add_to_napari(self, napari_viewer):
        """Add the layers data to a napari viewer"""

        # Get the reciprocal cells for drawing
        cells = self.viewer.draw_cells()

        # Add relps and cell (if present) as points and shapes layers for each id
        for exp_id in sorted(list(set(self.viewer.points_data["id"])), reverse=True):
            sel = self.viewer.points_data["id"] == exp_id

            # Convert points and colors to numpy arrays
            points = flumpy.to_numpy(self.viewer.points.select(sel))
            colors = flumpy.to_numpy(self.viewer.colors.select(sel))

            id_data = self.viewer.points_data["id"].select(sel)
            x, y, z = self.viewer.points_data["xyz"].select(sel).parts()
            panel_data = self.viewer.points_data["panel"].select(sel)
            d_spacing_data = self.viewer.points_data["d_spacing"].select(sel)
            point_properties = {
                "id": flumpy.to_numpy(id_data),
                "x": flumpy.to_numpy(x).round(1),
                "y": flumpy.to_numpy(y).round(1),
                "z": flumpy.to_numpy(z).round(1),
                "panel": flumpy.to_numpy(panel_data),
                "res": flumpy.to_numpy(d_spacing_data).round(3),
            }
            # text = "id:{} panel:{panel} xyz:{x}{y}{z} res:{res}Å"
            if "miller_index" in self.viewer.points_data and exp_id != -1:
                h, k, l = (
                    self.viewer.points_data["miller_index"]
                    .select(sel)
                    .as_vec3_double()
                    .parts()
                )
                h = h.iround()
                k = k.iround()
                l = l.iround()
                point_properties["h"] = flumpy.to_numpy(h)
                point_properties["k"] = flumpy.to_numpy(k)
                point_properties["l"] = flumpy.to_numpy(l)
            #    text += " hkl:{hkl}"
            # Currently not adding the text= to the points as this displays for
            # *every* point. Need a mouseover or tooltip instead. However, the
            # property values are displayed in the status bar, when the relevant
            # layer is selected

            relps_layer = napari_viewer.add_points(
                points,
                properties=point_properties,
                face_color=colors,
                size=self.settings.marker_size,
                name=f"relps id: {exp_id}",
            )
            relps_layer.blending = "translucent_no_depth"

            # Add the cell as a shapes layer, if it exists
            cell = cells.get(exp_id)
            if cell:
                cell_layer = napari_viewer.add_shapes(
                    cell.lines,
                    shape_type="line",
                    edge_width=0.5,
                    edge_color=np.array(cell.colors),
                    name=f"cell id: {exp_id}",
                )

                # Link the visibility of the relps and cell layer
                link_layers([relps_layer, cell_layer], ("visible",))

        # Now add rotation axis. Code extracted from draw_axis
        if self.viewer.minimum_covering_sphere is None:
            self.viewer.update_minimum_covering_sphere()
        s = self.viewer.minimum_covering_sphere
        scale = max(max(s.box_max()), abs(min(s.box_min())))
        axis = self.viewer.rotation_axis
        axis_line = np.array(
            [[0, 0, 0], [axis[0] * scale, axis[1] * scale, axis[2] * scale]]
        )

        axis_layer = napari_viewer.add_shapes(
            [
                axis_line,
            ],
            shape_type="line",
            edge_width=0.5,
            edge_color="white",
            name="axis",
        )

        # Set rotation around the origin
        napari_viewer.camera.center = (0, 0, 0)

        # Add the rlv_display widget and set values and limits
        napari_viewer.window.add_dock_widget(rlv_display, name="rlv display")
        rlv_display.marker_size.value = self.settings.marker_size
        rlv_display.d_min.value = self.settings.d_min
        rlv_display.d_min.min = self.settings.d_min
        rlv_display.z_min.value = self.settings.z_min
        rlv_display.z_min.min = self.settings.z_min
        rlv_display.z_min.max = self.settings.z_max
        rlv_display.z_max.value = self.settings.z_max
        rlv_display.z_max.min = self.settings.z_min
        rlv_display.z_max.max = self.settings.z_max

        return

    def load_models(self, experiments, reflections):
        Render3d.load_models(self, experiments, reflections)
        if self.settings.beam_centre is not None:

            pass
        if self.settings.marker_size is Auto:
            max_radius = max(self.reflections["rlp"].norms())
            volume = 4 / 3 * pi * max_radius ** 3
            density = len(self.reflections) / volume
            # Set marker size depending on relp density, where
            # 1000 < density < 20000 ==> max_size < marker_size < min_size
            # XXX this does not take into account narrow wedges!
            min_size, max_size = 1, 5
            grad = (max_size - min_size) / (20000 - 1000)
            intercept = max_size - 1000 * grad
            marker_size = grad * density + intercept
            marker_size = max(marker_size, min_size)
            marker_size = min(marker_size, max_size)
            self.settings.marker_size = marker_size

    def set_points(self):
        Render3d.set_points(self)

    def update_settings(self, *args, **kwds):
        self.set_beam_centre(self.settings.beam_centre_panel, self.settings.beam_centre)
        self.map_points_to_reciprocal_space()
        self.set_points()
        self.viewer.update_settings(*args, **kwds)


class RLVWindow:
    def __init__(self, settings, *args, **kwds):
        self.settings = settings
        self.points = flex.vec3_double()
        self.colors = None
        self.palette = None
        self.rotation_axis = None
        self.beam_vector = None
        self.recip_latt_vectors = None
        self.recip_crystal_vectors = None
        self.flag_show_minimum_covering_sphere = False
        self.minimum_covering_sphere = None
        self.field_of_view_y = 0.001

    def set_points(self, points):
        self.points = points
        self.points_display_list = None
        if self.minimum_covering_sphere is None:
            self.update_minimum_covering_sphere()

    def set_points_data(self, reflections):
        dstar = reflections["rlp"].norms()
        dstar.set_selected(dstar == 0, 1e-8)
        self.points_data = {
            "panel": reflections["panel"],
            "id": reflections["id"],
            "xyz": reflections["xyzobs.px.value"],
            "d_spacing": 1 / dstar,
        }
        if "miller_index" in reflections:
            self.points_data["miller_index"] = reflections["miller_index"]

    def set_colors(self, colors):
        assert len(colors) == len(self.points)
        self.colors = colors

    def set_palette(self, palette):
        self.palette = palette

    def set_rotation_axis(self, axis):
        self.rotation_axis = axis

    def set_beam_vector(self, beam):
        self.beam_vector = beam

    def set_reciprocal_lattice_vectors(self, vectors_per_crystal):
        self.recip_latt_vectors = vectors_per_crystal

    def set_reciprocal_crystal_vectors(self, vectors_per_crystal):
        self.recip_crystal_vectors = vectors_per_crystal

    def update_minimum_covering_sphere(self):
        n_points = min(1000, self.points.size())
        isel = flex.random_permutation(self.points.size())[:n_points]
        self.minimum_covering_sphere = minimum_covering_sphere(self.points.select(isel))

    def draw_cells(self):
        """Create a dictionary mapping experiment id to lists of lines and
        colours for reciprocal unit cells"""

        result = {}
        CellDrawing = namedtuple("CellDrawing", ["lines", "colors"])
        if self.settings.show_reciprocal_cell:
            # if we don't have one sort of vector we don't have the other either
            vectors = self.recip_latt_vectors
            if self.settings.crystal_frame:
                vectors = self.recip_crystal_vectors

            if vectors:
                for i, axes in enumerate(vectors):
                    if self.settings.experiment_ids:
                        if i not in self.settings.experiment_ids:
                            continue
                    j = (i + 1) % self.palette.size()
                    result[i] = CellDrawing(
                        lines=self.cell_edges(axes), colors=[self.palette[j]] * 12
                    )

        return result

    def cell_edges(self, axes):
        astar, bstar, cstar = axes[0], axes[1], axes[2]
        farpoint = astar + bstar + cstar

        lines = [
            np.array([[0, 0, 0], [*astar.elems]]),
            np.array([[0, 0, 0], [*bstar.elems]]),
            np.array([[0, 0, 0], [*cstar.elems]]),
            np.array([[*astar.elems], [*(astar + bstar).elems]]),
            np.array([[*astar.elems], [*(astar + cstar).elems]]),
            np.array([[*bstar.elems], [*(bstar + astar).elems]]),
            np.array([[*bstar.elems], [*(bstar + cstar).elems]]),
            np.array([[*cstar.elems], [*(cstar + astar).elems]]),
            np.array([[*cstar.elems], [*(cstar + bstar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - astar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - bstar).elems]]),
            np.array([[*farpoint.elems], [*(farpoint - cstar).elems]]),
        ]
        return lines
