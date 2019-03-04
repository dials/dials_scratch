from __future__ import print_function


if __name__ == "__main__":
    #    from scitbx import matrix
    #    from scitbx.array_family import flex
    #    from iotbx.xds import correction
    #    from cctbx import factor_ev_angstrom, factor_kev_angstrom

    factor_mev_angstrom = (6.6260755 * 2.99792458 / 1.60217733) / 1000.0

    width, height = (2463, 2527)
    wavelength = 0.979500

    #    print "Reading X correction"
    #    handle = correction.reader()
    #    handle.read_file('/home/upc86896/Data/X4_lots_M1S4_1_/GX-CORRECTIONS.cbf')
    #    xcorr = handle.get_correction((height, width))
    #
    #    print "Reading Y correction"
    #    handle = correction.reader()
    #    handle.read_file('/home/upc86896/Data/X4_lots_M1S4_1_/GY-CORRECTIONS.cbf')
    #    ycorr = handle.get_correction((height, width))

    table = [
        [1.00000e-03, 1.570e03, 1.567e03],
        [1.50000e-03, 5.355e02, 5.331e02],
        [1.83890e-03, 3.092e02, 3.070e02],
        [1.83890e-03, 3.192e03, 3.059e03],
        [2.00000e-03, 2.777e03, 2.669e03],
        [3.00000e-03, 9.784e02, 9.516e02],
        [4.00000e-03, 4.529e02, 4.427e02],
        [5.00000e-03, 2.450e02, 2.400e02],
        [6.00000e-03, 1.470e02, 1.439e02],
        [8.00000e-03, 6.468e01, 6.313e01],
        [1.00000e-02, 3.389e01, 3.289e01],
        [1.50000e-02, 1.034e01, 9.794e00],
        [2.00000e-02, 4.464e00, 4.076e00],
        [3.00000e-02, 1.436e00, 1.164e00],
        [4.00000e-02, 7.012e-01, 4.782e-01],
        [5.00000e-02, 4.385e-01, 2.430e-01],
        [6.00000e-02, 3.207e-01, 1.434e-01],
        [8.00000e-02, 2.228e-01, 6.896e-02],
        [1.00000e-01, 1.835e-01, 4.513e-02],
        [1.50000e-01, 1.448e-01, 3.086e-02],
        [2.00000e-01, 1.275e-01, 2.905e-02],
        [3.00000e-01, 1.082e-01, 2.932e-02],
        [4.00000e-01, 9.614e-02, 2.968e-02],
        [5.00000e-01, 8.748e-02, 2.971e-02],
        [6.00000e-01, 8.077e-02, 2.951e-02],
        [8.00000e-01, 7.082e-02, 2.875e-02],
        [1.00000e00, 6.361e-02, 2.778e-02],
        [1.25000e00, 5.688e-02, 2.652e-02],
        [1.50000e00, 5.183e-02, 2.535e-02],
        [2.00000e00, 4.480e-02, 2.345e-02],
        [3.00000e00, 3.678e-02, 2.101e-02],
        [4.00000e00, 3.240e-02, 1.963e-02],
        [5.00000e00, 2.967e-02, 1.878e-02],
        [6.00000e00, 2.788e-02, 1.827e-02],
        [8.00000e00, 2.574e-02, 1.773e-02],
        [1.00000e01, 2.462e-02, 1.753e-02],
        [1.50000e01, 2.352e-02, 1.746e-02],
        [2.00000e01, 2.338e-02, 1.757e-02],
    ]

    from math import log, exp
    import numpy
    from scipy.interpolate import interp1d

    energy, mac, meac = zip(*table)
    energy = numpy.array(energy, dtype=numpy.float32)
    mac = numpy.array(mac, dtype=numpy.float32)

    beam_energy = factor_mev_angstrom / wavelength

    #    from matplotlib import pylab
    #    pylab.plot(energy, mac)
    #    pylab.show()

    print(energy)
    print(mac)

    #    x0 = 0
    #    x1 = 0
    #    x2 = 0
    #
    #    y0 = 0
    #    y1 = 0
    #    y2 = 0

    #    a11 = 2.0 / (x1 - x0)
    #    a12 = 1.0 / (x1 - x0)
    #    a21 = 1.0 / (x1 - x0)
    #    a22 = 2.0 * (1.0 / (x1 - x0) + 1.0 / (x2 - x1))
    #    a23 = 1.0 / (x2 - x1)
    #    a32 = 1.0 / (x2 - x1)
    #    a33 = 2.0 / (x2 - x1)
    #    b1 = 3.0 * (y1 - y0) / (x1 - x0)**2
    #    b2 = 3.0 * ((y1 - y0) / (x1 - x0)**2 + (y2 - y1) / (x2 - x1)**2)
    #    b3 = 3.0 * (y2 - y1) / (x2 - x1)**2

    log_energy = numpy.log(energy)
    log_mac = numpy.log(mac)

    print(log(beam_energy))

    #    from matplotlib import pylab
    #    pylab.plot(energy, log_mac)
    #    pylab.show()

    #    print 1/0

    #    from scipy.interpolate import splrep, splev
    #    tck = splrep(log_energy, log_mac)
    #    print tck
    #    print 1/0
    #    coeff = splev(log(beam_energy), tck, der=0)
    #    f2 = interp1d(numpy.log(energy), numpy.log(mac), kind='linear')
    #    coeff = exp(f2(log(beam_energy)))
    coeff = numpy.exp(numpy.interp(log(beam_energy), numpy.log(energy), numpy.log(mac)))
    #    coeff = 16.9612
    rho = 2.330e00

    print("Energy: {0} Mev".format(beam_energy))
    print("Energy: {0} ev".format(beam_energy * 1000000.0))
    print("Mu/rho: {0} cm^2/g".format(coeff))
    print("Mu:     {0} cm^-1".format(coeff * rho))
    print("rho:    {0} g / cm^3".format(rho))
    from math import exp

    x_arr = []
    p_arr = []

    print("Att Len: {0} microm".format((1.0 / (coeff * rho)) * 10 * 1000))

    for xx in range(0, 100):
        xcm = xx / 1000.0
        xmm = xcm * 10
        p = exp(-coeff * rho * xcm)
        x_arr.append(xmm)
        p_arr.append(p)

#    from matplotlib import pylab
#    pylab.plot(x_arr, p_arr)
#    pylab.axhline(1.0/exp(1.0))
#    pylab.show()
