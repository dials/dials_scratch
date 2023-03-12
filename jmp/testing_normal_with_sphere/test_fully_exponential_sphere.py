from __future__ import print_function

from scitbx import matrix
from math import pi, cos, sin, exp, sqrt, atan2, log, tan, acos


def normal_3d(x, mu, sigma):
    A = 1.0 / sqrt((2 * pi) ** 3 * sigma.determinant())
    B = ((x - mu).transpose() * sigma.inverse() * (x - mu))[0]
    return A * exp(-0.5 * B)


def func_f(theta, phi, mu, s_sq):
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    A = (sin_theta * cos_phi - mu[0]) ** 2 / s_sq[0]
    B = (sin_theta * sin_phi - mu[1]) ** 2 / s_sq[1]
    C = (cos_theta - mu[2]) ** 2 / s_sq[2]
    return 0.5 * (A + B + C)


def func_f_gradient(theta, phi, mu, s_sq):
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    A = (sin_theta * cos_theta * cos_phi**2 - mu[0] * cos_theta * cos_phi) / s_sq[0]
    B = (sin_theta * cos_theta * sin_phi**2 - mu[1] * cos_theta * sin_phi) / s_sq[1]
    C = (cos_theta * sin_theta - mu[2] * sin_theta) / s_sq[2]
    D = (sin_theta**2 * sin_phi * cos_phi - mu[0] * sin_theta * sin_phi) / s_sq[0]
    E = (sin_theta**2 * sin_phi * cos_phi - mu[1] * sin_theta * cos_phi) / s_sq[1]
    return matrix.col((A + B - C, -D + E))


def func_f_hessian(theta, phi, mu, s_sq):
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    cos_2theta = cos(2 * theta)
    sin_2theta = sin(2 * theta)
    sin_2phi = sin(2 * phi)
    cos_2phi = cos(2 * phi)
    A = (cos_2theta * cos_theta**2 + mu[0] * sin_theta * cos_phi) / s_sq[0]
    B = (cos_2theta * sin_theta**2 + mu[1] * sin_theta * sin_phi) / s_sq[1]
    C = (cos_2theta - mu[2] * cos_theta) / s_sq[2]
    D = (0.5 * sin_2theta * sin_2phi - mu[0] * cos_theta * sin_phi) / s_sq[0]
    E = (0.5 * sin_2theta * sin_2phi - mu[1] * cos_theta * cos_phi) / s_sq[1]
    F = (sin_theta**2 * cos_2phi - mu[0] * sin_theta * cos_phi) / s_sq[0]
    G = (sin_theta**2 * cos_2phi + mu[1] * sin_theta * sin_phi) / s_sq[1]
    H1 = A + B - C
    H2 = -D + E
    H3 = -F + G
    return matrix.sqr((H1, H2, H2, H3))


def compute_peak_f(mu, sigma):
    def gradient(theta, phi, mu, s_sq):
        return func_f_gradient(theta, phi, mu, s_sq)

    def hessian(theta, phi, mu, s_sq):
        return func_f_hessian(theta, phi, mu, s_sq)

    phi = atan2(mu[1], mu[0])
    theta = acos(mu[2] / mu.length())

    x0 = matrix.col((theta, phi))

    s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

    while True:

        G = gradient(x0[0], x0[1], mu, s_sq)
        H = hessian(x0[0], x0[1], mu, s_sq)
        x = x0 - H.inverse() * G
        if (x - x0).length() < 1e-7:
            break

        x0 = x

    return x


def compute_peak_fx(mu, sigma):
    def gradient(theta, phi, mu, s_sq):
        G = func_f_gradient(theta, phi, mu, s_sq)
        return G + matrix.col((-1 / tan(theta), tan(phi)))

    def hessian(theta, phi, mu, s_sq):
        H = func_f_hessian(theta, phi, mu, s_sq)
        return H + matrix.sqr((1 / sin(theta) ** 2, 0, 0, 1 / cos(phi) ** 2))

    phi = atan2(mu[1], mu[0])
    theta = acos(mu[2] / mu.length())

    x0 = matrix.col((theta, phi))

    s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

    while True:

        G = gradient(x0[0], x0[1], mu, s_sq)
        H = hessian(x0[0], x0[1], mu, s_sq)
        x = x0 - H.inverse() * G
        if (x - x0).length() < 1e-7:
            break

        x0 = x

    return x


def compute_peak_fy(mu, sigma):
    def gradient(theta, phi, mu, s_sq):
        G = func_f_gradient(theta, phi, mu, s_sq)
        return G + matrix.col((-1 / tan(theta), -1 / tan(phi)))

    def hessian(theta, phi, mu, s_sq):
        H = func_f_hessian(theta, phi, mu, s_sq)
        return H + matrix.sqr((1 / sin(theta) ** 2, 0, 0, 1 / sin(phi) ** 2))

    phi = atan2(mu[1], mu[0])
    theta = acos(mu[2] / mu.length())

    x0 = matrix.col((theta, phi))

    s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

    while True:

        G = gradient(x0[0], x0[1], mu, s_sq)
        H = hessian(x0[0], x0[1], mu, s_sq)
        x = x0 - H.inverse() * G
        if (x - x0).length() < 1e-7:
            break

        x0 = x

    return x


def compute_peak_fz(mu, sigma):
    def gradient(theta, phi, mu, s_sq):
        G = func_f_gradient(theta, phi, mu, s_sq)
        return G + matrix.col((tan(theta), 0))

    def hessian(theta, phi, mu, s_sq):
        H = func_f_hessian(theta, phi, mu, s_sq)
        return H + matrix.sqr((1 / cos(theta) ** 2, 0, 0, 0))

    phi = atan2(mu[1], mu[0])
    theta = acos(mu[2] / mu.length())

    x0 = matrix.col((theta, phi))

    s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))

    while True:

        G = gradient(x0[0], x0[1], mu, s_sq)
        H = hessian(x0[0], x0[1], mu, s_sq)
        x = x0 - H.inverse() * G
        if (x - x0).length() < 1e-7:
            break

        x0 = x

    return x


def compute_mean_estimate(peak_f, peak_fx, peak_fy, peak_fz, mu, sigma):

    sigma_inv = sigma.inverse()

    s1, s2, s3 = sigma[0], sigma[4], sigma[8]

    H_f = matrix.sqr((-1 / s1, 0, 0, 0, -1 / s2, 0, 0, 0, -1 / s3))

    H_fx = H_f + matrix.sqr((-1 / peak_fx[0] ** 2, 0, 0, 0, 0, 0, 0, 0, 0))

    H_fy = H_f + matrix.sqr((0, 0, 0, 0, -1 / peak_fy[1] ** 2, 0, 0, 0, 0))

    H_fz = H_f + matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, -1 / peak_fz[2] ** 2))

    E_f = exp((-0.5 * (peak_f - mu).transpose() * sigma_inv * (peak_f - mu))[0])
    E_fx = peak_fx[0] * exp(
        (-0.5 * (peak_fx - mu).transpose() * sigma_inv * (peak_fx - mu))[0]
    )
    E_fy = peak_fy[1] * exp(
        (-0.5 * (peak_fy - mu).transpose() * sigma_inv * (peak_fy - mu))[0]
    )
    E_fz = peak_fz[2] * exp(
        (-0.5 * (peak_fz - mu).transpose() * sigma_inv * (peak_fz - mu))[0]
    )

    A = (1 / sqrt((-H_f).determinant())) * E_f
    B = (1 / sqrt((-H_fx).determinant())) * E_fz
    C = (1 / sqrt((-H_fy).determinant())) * E_fy
    D = (1 / sqrt((-H_fz).determinant())) * E_fz

    xc = B / A
    yc = C / A
    zc = D / A

    v = matrix.col((xc, yc, zc)).normalize()

    return v


def compute_laplace(mu, sigma):

    theta_f, phi_f = compute_peak_f(mu, sigma)
    theta_fx, phi_fx = compute_peak_fx(mu, sigma)
    theta_fy, phi_fy = compute_peak_fy(mu, sigma)
    theta_fz, phi_fz = compute_peak_fz(mu, sigma)

    # print theta_f, phi_f
    # print theta_fx, phi_fx
    # print theta_fy, phi_fy
    # print theta_fz, phi_fz

    peak_f = matrix.col(
        (sin(theta_f) * cos(phi_f), sin(theta_f) * sin(phi_f), cos(theta_f))
    )

    peak_fx = matrix.col(
        (sin(theta_fx) * cos(phi_fx), sin(theta_fx) * sin(phi_fx), cos(theta_fx))
    )

    peak_fy = matrix.col(
        (sin(theta_fy) * cos(phi_fy), sin(theta_fy) * sin(phi_fy), cos(theta_fy))
    )

    peak_fz = matrix.col(
        (sin(theta_fz) * cos(phi_fz), sin(theta_fz) * sin(phi_fz), cos(theta_fz))
    )

    x = compute_mean_estimate(peak_f, peak_fx, peak_fy, peak_fz, mu, sigma)

    theta = acos(x[2])
    phi = atan2(x[1], x[0])
    print(theta, phi)

    return x
    # mu_phi = atan2(mu[1], mu[0])
    # mu_theta = acos(mu[2] / mu.length())
    # s_sq = matrix.col((sigma[0], sigma[4], sigma[8]))
    # print mu_theta, mu_phi, func_f(mu_theta, mu_phi, mu, s_sq)
    # print theta, phi, func_f(theta, phi, mu, s_sq)


def compute_integrate(mu, sigma):

    I1 = matrix.col((0, 0, 0))
    I2 = 0
    N = 1000
    M = 1000
    T0 = 0
    T1 = pi / 2
    P0 = 0
    P1 = pi / 2
    for j in range(N):
        for i in range(M):
            theta = T0 + j * (T1 - T0) / N
            phi = P0 + i * (P1 - P0) / N

            sin_theta = sin(theta)
            cos_theta = cos(theta)
            sin_phi = sin(phi)
            cos_phi = cos(phi)

            x = matrix.col((sin_theta * cos_phi, sin_theta * sin_phi, cos_theta))

            f = normal_3d(x, mu, sigma)

            I1 += x * f * sin_theta
            I2 += f * sin_theta

    return (I1 / I2).normalize()


if __name__ == "__main__":

    mu = matrix.col((1, 1, 1)).normalize() * 0.95
    sigma = matrix.sqr((0.01, 0, 0, 0, 0.02, 0, 0, 0, 0.03))

    scal = compute_laplace(mu, sigma)

    print(scal.length())

    sest = compute_integrate(mu, sigma)

    theta = acos(sest[2])
    phi = atan2(sest[1], sest[0])
    print(theta, phi)

    print(scal.angle(sest))
