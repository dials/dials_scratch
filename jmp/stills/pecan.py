from __future__ import division
from __future__ import print_function


class PixelLabeller(object):
    """
    Class to do pixel labelling

    """

    def __init__(self, experiments):
        """
        Do the labelling

        """
        from dials.algorithms.spot_prediction import PixelToMillerIndex
        from collections import defaultdict
        from math import floor, sqrt
        from dials.array_family import flex

        # Get the experiment
        experiment = experiments[0]

        # Get the image size
        xsize, ysize = experiment.detector[0].get_image_size()

        # A class to map pixels to miller indices
        transform = PixelToMillerIndex(
            experiment.beam, experiment.detector, experiment.crystal
        )

        # For each pixel, assign to a miller index and also compute the distance
        reflections = defaultdict(list)
        for j in range(ysize):
            for i in range(xsize):
                h = transform.h(0, i, j)
                h0 = tuple(map(lambda x: int(floor(x + 0.5)), h))
                d = sqrt(sum(map(lambda x, y: (x - y) ** 2, h, h0)))
                reflections[h0].append((i, j, d))

        # Initialise arrays
        self._indices = flex.miller_index()
        self._distance = flex.double(flex.grid(ysize, xsize))
        self._label = flex.int(flex.grid(ysize, xsize))
        self._pixels = defaultdict(list)
        for index, (h, pixels) in enumerate(reflections.iteritems()):
            self._indices.append(h)
            for i, j, d in pixels:
                self._distance[j, i] = d
                self._label[j, i] = index
                self._pixels[h].append((i, j))

    def label(self):
        """
        Return an array of labels

        """
        return self._label

    def distance(self):
        """
        Return the distance at each pixel to the miller index

        """
        return self._distance

    def indices(self):
        """
        Return the list of miller indices

        """
        return self._indices

    def pixels(self):
        """
        Return the list of pixels at each reflection

        """
        return self._pixels


class ImageModeller(object):
    """
    Class to do image modelling

    """

    def __init__(self, experiments, mosaicity):
        """
        Do the labelling

        """
        from dials.algorithms.spot_prediction import PixelToMillerIndex
        from collections import defaultdict
        from math import floor, sqrt, pi, exp
        from dials.array_family import flex
        from dials_scratch.jmp.stills import Simulator
        from scitbx import matrix

        # Get the experiment
        experiment = experiments[0]

        # # Get the image size
        # xsize, ysize = experiment.detector[0].get_image_size()

        # # A class to map pixels to miller indices
        # transform = PixelToMillerIndex(
        #   experiment.beam,
        #   experiment.detector,
        #   experiment.crystal)

        data = experiment.imageset.get_raw_data(0)[0]
        mask = experiment.imageset.get_mask(0)[0]

        # For each pixel, assign to a miller index and also compute the distance
        reflections = defaultdict(list)
        simulation = Simulator(
            experiment.beam,
            experiment.detector,
            experiment.crystal,
            mosaicity,
            5,
            data,
            mask,
        )

        reflections = simulation.reflections()
        I = reflections["intensity.sum.value"]
        V = reflections["intensity.sum.variance"]
        S = reflections["intensity.scale"]
        selection = (V > 0) & (S > 0)
        reflections = reflections.select(selection)
        I = reflections["intensity.sum.value"]
        V = reflections["intensity.sum.variance"]
        S = reflections["intensity.scale"]
        selection = I / flex.sqrt(V) > 1
        reflections = reflections.select(selection)

        reflections.as_pickle("integrated.pickle")

        self._distance = simulation.model()
        self._mask = simulation.mask()
        self._labels = simulation.labels()
        self._indices = simulation.indices()
        print(len(self._indices))

        from matplotlib import pylab

        pylab.imshow(self._distance.as_numpy_array(), interpolation="none")
        pylab.show()

        from matplotlib import pylab

        pylab.imshow(self._mask.as_numpy_array(), interpolation="none")
        pylab.show()

        from matplotlib import pylab

        pylab.imshow(self._labels.as_numpy_array(), interpolation="none")
        pylab.show()
        exit(0)
        # for j in range(ysize):
        #   print j
        #   for i in range(xsize):
        #     h = transform.h(0, i, j)
        #     h0 = tuple(map(lambda x: int(floor(x+0.5)), h))
        #     #d = sqrt(sum(map(lambda x,y: (x-y)**2, h, h0)))

        #     A = matrix.col(h)
        #     B = matrix.col(transform.h(0, i+1, j))
        #     C = matrix.col(transform.h(0, i, j+1))
        #     D = matrix.col(transform.h(0, i+1, j+1))
        #     area1 = 0.5 * (B-A).cross(C-A).length()
        #     area2 = 0.5 * (B-D).cross(C-D).length()
        #     area = area1 + area2

        #     sum_f = 0
        #     for jj in range(5):
        #       for ii in range(5):
        #         h = transform.h(0, i+ii/5.0, j+jj/5.0)
        #         d = sqrt(sum(map(lambda x,y: (x-y)**2, h, h0)))
        #         sum_f += (1.0/(sqrt(2*pi)*mosaicity))**3 * exp(-0.5*(d/mosaicity)**2)

        #     I = area * sum_f / (5.0*5.0)
        #     reflections[h0].append((i,j,I))

        # Initialise arrays
        # self._indices = flex.miller_index()
        # self._distance = flex.double(flex.grid(ysize, xsize))
        # self._label = flex.int(flex.grid(ysize, xsize))
        # self._pixels = defaultdict(list)
        # for index, (h, pixels) in enumerate(reflections.iteritems()):
        #   self._indices.append(h)
        #   for i, j, d in pixels:
        #     self._distance[j,i] = d
        #     self._label[j,i] = index
        #     self._pixels[h].append((i,j))

    def label(self):
        """
        Return an array of labels

        """
        return self._label

    def distance(self):
        """
        Return the distance at each pixel to the miller index

        """
        return self._distance

    def indices(self):
        """
        Return the list of miller indices

        """
        return self._indices

    def pixels(self):
        """
        Return the list of pixels at each reflection

        """
        return self._pixels


class ProfileRefiner(object):
    def __init__(self, experiments, bandwidth=0.035, mosaicity=0.01):

        self.experiments = experiments
        self.bandwidth = bandwidth
        self.mosaicity = mosaicity

    def refine(self):

        self.predict()

    def predict(self):

        modeller = ImageModeller(self.experiments, self.mosaicity)

        from dials.array_family import flex

        data = labeller.distance()

        from matplotlib import pylab

        pylab.imshow(data.as_numpy_array(), interpolation="none")
        pylab.show()


if __name__ == "__main__":

    from dxtbx.model.experiment_list import ExperimentListFactory
    import sys

    experiments_filename = sys.argv[1]

    # Read the experiments
    experiments = ExperimentListFactory.from_json_file(experiments_filename)

    # Setup the profile refiner
    refiner = ProfileRefiner(experiments, bandwidth=0.035, mosaicity=0.1)

    # Do the refinement
    refiner.refine()
