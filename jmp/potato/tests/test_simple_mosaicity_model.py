from __future__ import division
from __future__ import print_function
from scitbx import matrix


def tst_simple_mosaicity_model():

    from dials_scratch.jmp.potato.model import SimpleMosaicityModel
    from dials_scratch.jmp.potato.model import SimpleMosaicityParameterisation

    sigma = matrix.sqr((1, 0, 0, 0, 2, 0, 0, 0, 3))

    model = SimpleMosaicityModel(sigma)

    parameterisation = model.parameterisation()

    params = parameterisation.parameters()

    expected = (1.0, 0.0, 1.4142135623730951, 0.0, 0.0, 1.7320508075688772)

    assert all(abs(a - b) < 1e-7 for a, b in zip(params, expected))

    model.compose(parameterisation)

    sigma2 = model.sigma()

    assert all(abs(a - b) < 1e-7 for a, b in zip(sigma, sigma2))

    print("OK")


if __name__ == "__main__":
    tst_simple_mosaicity_model()
