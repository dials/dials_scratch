from __future__ import division
from __future__ import print_function


class History(object):
    """
  Store the history of parameter values and likelihood

  """

    def __init__(self, parameter_names):
        """
    Initialise with parameter names

    """
        self.names = parameter_names
        self.values = []
        self.logL = []

    def append(self, parameters, log_likelihood):
        """
    Append a new entry

    """
        assert len(parameters) == len(self.names)
        self.values.append(parameters)
        self.logL.append(log_likelihood)

    def __len__(self):
        """
    Get the number of items

    """
        return len(self.values)

    def __getitem__(self, index):
        """
    Return an item

    """
        return (self.values, self.logL)

    def __iter__(self):
        """
    Iterate through items

    """
        for i in range(len(self)):
            yield self[i]

    def plot(self):
        """
    Plot the history

    """
        from matplotlib import pylab

        fig, (ax1, ax2) = pylab.subplots(2, 1, sharex=True)
        x = list(range(len(self.values)))
        parameters = zip(*self.values)
        for i, p in enumerate(parameters):
            ax1.plot(x, p, label="%s" % self.names[i])
        ax1.legend()
        ax2.plot(x, self.logL)
        ax2.set_xlabel("Iteration")
        ax2.set_ylabel("log(L)")
        ax1.set_ylabel("Parameter value")
        pylab.show()


class ProfileModeller(object):
    """
  A class to do profile modelling

  """

    def __init__(
        self,
        experiments,
        reflections,
        macro_cycles=[10, 100],
        num_integral=5,
        use_mosaic_block_angular_spread=False,
        use_wavelength_spread=False,
    ):
        """
    Init

    """
        from dials_scratch.jmp.profile_modelling import MLTarget3D

        # Save some stuff
        self.experiments = experiments
        self.reflections = reflections
        self.macro_cycles = macro_cycles
        self.num_integral = num_integral
        self.use_mosaic_block_angular_spread = use_mosaic_block_angular_spread
        self.use_wavelength_spread = use_wavelength_spread
        self.parameters = None
        self.simplex = None

        # Get the parameter names
        names = MLTarget3D.parameter_names(
            use_mosaic_block_angular_spread=use_mosaic_block_angular_spread,
            use_wavelength_spread=use_wavelength_spread,
        )

        # Create the history
        self.history = History(list(names))

        # Do the estimation
        for n in self.macro_cycles:
            self.estimate(n)

    def estimate(self, num_reflections):
        """
    Estimate the model parameters

    """
        from scitbx import simplex
        from copy import deepcopy
        from dials.array_family import flex
        from random import uniform

        # Select the reflections
        reflections = self._select_reflections(self.reflections, num_reflections)

        # The number of parameters
        num_parameters = len(self.history.names)

        # Setup the starting simplex
        if self.simplex is None:
            self.simplex = []
            for i in range(num_parameters + 1):
                self.simplex.append(
                    flex.log(
                        flex.double(
                            [uniform(0.0001, 0.01) for j in range(num_parameters)]
                        )
                    )
                )

        class Evaluator(object):
            """
      Evaluator to simplex

      """

            def __init__(
                self,
                history,
                experiment,
                reflections,
                num_integral,
                use_mosaic_block_angular_spread,
                use_wavelength_spread,
            ):
                """
        Initialise

        """
                from dials_scratch.jmp.profile_modelling import MLTarget3D

                self.func = MLTarget3D(
                    experiment,
                    reflections,
                    num_integral=num_integral,
                    use_mosaic_block_angular_spread=use_mosaic_block_angular_spread,
                    use_wavelength_spread=use_wavelength_spread,
                )
                self.count = 1
                self.history = history
                self.logL = None

            def target(self, log_parameters):
                """
        Compute the negative log likelihood

        """
                from dials.array_family import flex

                parameters = flex.exp(log_parameters)

                self.logL = self.func.log_likelihood(parameters)
                logL = flex.sum(self.logL)

                self.count += 1

                print(self.count, list(parameters), logL)

                self.history.append(parameters, logL)

                # Return negative log likelihood
                return -logL

        # Setup the simplex optimizer
        optimizer = simplex.simplex_opt(
            num_parameters,
            matrix=self.simplex,
            evaluator=Evaluator(
                self.history,
                self.experiments[0],
                reflections,
                num_integral=self.num_integral,
                use_mosaic_block_angular_spread=self.use_mosaic_block_angular_spread,
                use_wavelength_spread=self.use_wavelength_spread,
            ),
            tolerance=1e-7,
        )

        # Get the solution
        self.parameters = flex.exp(optimizer.get_solution())

        # Get the final simplex
        self.simplex = optimizer.matrix

        # Save the likelihood for each reflection
        self.log_likelihood = optimizer.evaluator.logL

    def _select_reflections(self, reflections, num):
        """
    Select reflections to use

    """

        def select_used_in_refinement(reflections):
            """
      Select reflections to use

      """
            selection = reflections.get_flags(reflections.flags.used_in_refinement)
            print(
                "Selecting %d/%d strong reflections"
                % (len(reflections), selection.count(True))
            )
            print("")
            return reflections.select(selection)

        def sort_by_intensity(reflections):
            """
      Sort the reflections by intensity

      """
            reflections.sort("intensity.sum.value", reverse=True)
            return reflections

        def select_subset(reflections, num):
            """
      Select most intense reflections

      """
            print("Selecting %d/%d reflections" % (num, len(reflections)))
            return reflections[0:num]

        # Select the strong reflections
        reflections = select_used_in_refinement(reflections)
        reflections = sort_by_intensity(reflections)

        # Select a subset of strong reflections
        return select_subset(reflections, num)

    def display(self, num):
        """
    Display some shoeboxes

    """
        from dials_scratch.jmp.viewer import show_image_stack_multi_view
        from random import sample
        from dials.array_family import flex

        def simulate(experiments, reflections, parameters):
            from dials_scratch.jmp.profile_modelling import MLTarget3D

            func = MLTarget3D(experiments[0], reflections)
            return [func.simulate(i, parameters) for i in range(len(reflections))]

        # Sample from reflections
        reflections = self.reflections.select(
            flex.size_t(sample(range(len(self.reflections)), num))
        )

        # Simulate the reflection profiles from the current model
        simulated = simulate(self.experiments, reflections, self.parameters)

        # Display stuff
        for model, data_sbox in zip(simulated, reflections["shoebox"]):
            data = data_sbox.data
            show_image_stack_multi_view(model.as_numpy_array(), vmax=flex.max(model))
            show_image_stack_multi_view(data.as_numpy_array(), vmax=flex.max(data))


if __name__ == "__main__":

    def read_experiments(filename):
        """
    Read the experiments file

    """
        from dxtbx.model.experiment_list import ExperimentListFactory

        print("Reading experiments from %s" % filename)
        print("")
        return ExperimentListFactory.from_json_file(filename)

    def read_reflections(filename):
        """
    Read the reflections file

    """
        from dials.array_family import flex

        print("Reading reflections from %s" % filename)
        print("")
        return flex.reflection_table.from_pickle(filename)

    def model_profiles(experiments, reflections):
        modeller = ProfileModeller(experiments, reflections)
        return modeller.model

    # Hard code the filenames
    experiments_filename = "/home/upc86896/Data/bag_training/processed_profile/profile_model/experiments.json"
    reflections_filename = "/home/upc86896/Data/bag_training/processed_profile/profile_model/reflections.pickle"

    # Read the experiments and reflections
    experiments = read_experiments(experiments_filename)
    reflections = read_reflections(reflections_filename)

    # Generate the profile model
    modeller = ProfileModeller(
        experiments,
        reflections,
        macro_cycles=[10],
        num_integral=10,
        use_mosaic_block_angular_spread=False,
        use_wavelength_spread=False,
    )
    modeller.history.plot()
    if False:
        modeller.display(5)
    parameters = modeller.parameters

    print("Generated final model:")
    print("")
    print(list(parameters))
    print("")
