from __future__ import division
from __future__ import print_function
import sys, os
from libtbx import easy_mp, easy_run
from libtbx.phil import parse
from libtbx.utils import Sorry

"""
Jiffy scipt for doing dials.import and dials.index on a pile of images using easy_mp. Each image should be paired with a strong.pickle like this:
filename.cbf
filename_strong.pickle

Example usage:
libtbx.python mp_index.py path_to_images index.phil
"""

phil_str = """
  output_dir = .
    .type = str
  mp {
    nproc = 1
      .type = int
    method = *multiprocessing
      .type = choice
  }
  image_extension = .cbf
    .type = str
  reference_geometry = None
    .type = str
"""

phil_scope = parse(phil_str)

user_phil = []
root_dirs = []
indexing_phil = None
for arg in sys.argv[1:]:
    if os.path.isdir(arg):
        root_dirs.append(arg)
    elif os.path.isfile(arg):
        assert indexing_phil is None
        indexing_phil = arg
    else:
        try:
            user_phil.append(parse(arg))
        except Exception:
            raise Sorry("Couldn't parse argument %s" % arg)

params = phil_scope.fetch(sources=user_phil).extract()

print("Finding files")

images = []
strongs = []

for root in root_dirs:
    for filename in os.listdir(root):
        if os.path.splitext(filename)[1] != params.image_extension:
            continue
        filepath = os.path.join(root, filename)
        strong_filepath = os.path.join(
            root, os.path.splitext(filename)[0] + "_strong.pickle"
        )
        if not os.path.exists(strong_filepath):
            raise Sorry("Couldn't find spotfinding results for image %s" % filepath)
        images.append(filepath)
        strongs.append(strong_filepath)

print("Found %d images to index" % len(images))


def index(item):
    image, strong = item
    base = os.path.splitext(os.path.basename(image))[0]
    datablock = os.path.join(params.output_dir, base + "_datablock.json")
    command = "dials.import %s output.datablock=%s" % (image, datablock)
    if params.reference_geometry is not None:
        command += " reference_geometry=%s" % params.reference_geometry
    easy_run.fully_buffered(command).raise_if_errors().show_stdout()

    command = "dials.index %s %s output.experiments=%s output.reflections=%s" % (
        datablock,
        strong,
        os.path.join(params.output_dir, base + "_experiments.json"),
        os.path.join(params.output_dir, base + "_indexed.pickle"),
    )
    if indexing_phil is not None:
        command += " %s" % indexing_phil

    easy_run.fully_buffered(command).show_stdout()


easy_mp.parallel_map(
    func=index,
    iterable=zip(images, strongs),
    processes=params.mp.nproc,
    method=params.mp.method,
    preserve_order=False,
    preserve_exception_message=True,
)

print("All done")
