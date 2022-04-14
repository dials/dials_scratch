from __future__ import annotations

from math import pi
import numpy as np

import libtbx.phil

from libtbx import Auto
from scitbx.array_family import flex
from scitbx.math import minimum_covering_sphere

from dials_scratch.dgw.rlv import Render3d
from dxtbx import flumpy

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


class ReciprocalLatticeViewer(Render3d):
    def __init__(self, parent, id, title, size, settings=None, *args, **kwds):
        Render3d.__init__(self, settings=settings)

        self.viewer = RLVWindow(
            settings=self.settings,
        )

    def add_to_napari(self, napari_viewer):
        """Add the layers data to a napari viewer"""

        # NB In dials.reciprocal_lattice_viewer, the drawing is handled by an
        # RLVWindow instance, which here is self.viewer.

        # Convert points and colors to numpy arrays
        points = flumpy.to_numpy(self.viewer.points)
        colors = flumpy.to_numpy(self.viewer.colors)

        # Set point labels (could be slow!)
        labels = []
        xyz_data = self.viewer.points_data["xyz"]
        id_data = self.viewer.points_data["id"]
        panel_data = self.viewer.points_data["panel"]
        d_spacing_data = self.viewer.points_data["d_spacing"]
        for xyz, exp_id, panel, d_spacing in zip(
            xyz_data, id_data, panel_data, d_spacing_data
        ):

            label = (
                f"id: {exp_id}; panel: {panel}\n"
                f"xyz: {xyz[0]:.1f} {xyz[1]:.1f} {xyz[2]:.1f}\n"
                f"res: {d_spacing:.2f} Angstrom"
            )
            labels.append(label)
        if "miller_index" in self.viewer.points_data:
            for i, (exp_id, hkl) in enumerate(
                zip(id_data, self.viewer.points_data["miller_index"])
            ):
                if exp_id != -1:
                    labels[i] += f"\nhkl: {hkl}"

        # These labels are displayed on mouseover of the points
        point_properties = {
            "label": labels,
        }

        # Add the points layer. Only a single layer for now, with all points
        points_layer = napari_viewer.add_points(
            points,
            properties=point_properties,
            face_color=colors,
            size=self.settings.marker_size,
        )

        # Now add rotation axis. Code extracted from draw_axis
        if self.viewer.minimum_covering_sphere is None:
            self.viewer.update_minimum_covering_sphere()
        s = self.viewer.minimum_covering_sphere
        scale = max(max(s.box_max()), abs(min(s.box_min())))
        axis = self.viewer.rotation_axis
        line = np.array(
            [[0, 0, 0], [axis[0] * scale, axis[1] * scale, axis[2] * scale]]
        )

        # Struggling to add these shapes individually to separate layers. Let's
        # just accumulate all shapes and then add them all at once
        shapes = [line]
        edge_colors = [(1, 1, 1)]

        axis_layer = napari_viewer.add_shapes(
            shapes,
            shape_type="line",
            edge_width=0.1,
            edge_color=edge_colors,
            name="axis",
        )

        # Add reciprocal cells
        cells, cell_colors = self.viewer.draw_cells()
        # shapes.extend(cells)
        # edge_colors.extend(cell_colors)

        edge_colors = np.array(edge_colors)
        cells_layer = napari_viewer.add_shapes(
            cells,
            shape_type="line",
            edge_width=0.1,
            edge_color=cell_colors,
            name="cells",
        )

        return

    def load_models(self, experiments, reflections):
        Render3d.load_models(self, experiments, reflections)
        if self.settings.beam_centre is not None:

            pass
        if self.settings.marker_size is Auto:
            max_radius = max(self.reflections["rlp"].norms())
            volume = 4 / 3 * pi * max_radius ** 3
            density = len(self.reflections) / volume
            # Set marker size to between 0.05 and 0.5 depending on density, where
            # 1000 < density < 20000 ==> 5 < marker_size < 0.5
            marker_size = (-0.45 / 19000) * density + (0.05 + 9 / 19)
            marker_size = max(marker_size, 0.5)
            marker_size = min(marker_size, 5)
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
        """Create a list of lines corresponding to reciprocal unit cells"""

        lines = []
        colors = []
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
                    colors.extend([self.palette[j]] * 12)
                    lines.extend(self.cell_edges(axes))
        return lines, colors

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
