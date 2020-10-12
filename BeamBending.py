# program works for beams with i number of rectangular sections

import numpy as np


def circle_moment_inertia(r):
    Ixx = 1/4 * np.pi * r**4
    Iyy = Ixx
    Ixy = 0
    return Ixx, Iyy, Ixy


def hollow_circle_moment_inertia(r1,r2):
    Ixx = 1/4 * np.pi * (r2 ** 4 - r1 ** 4)
    Iyy = Ixx
    Ixy = 0
    return Ixx, Iyy, Ixy


def rectangle_moment_inertia(b, h):
    Ixx = 1 / 12 * b * h ** 3
    Iyy = 1 / 12 * h * b ** 3
    Ixy = 0
    return Ixx, Iyy, Ixy


def parallel_axis_theorem(a, d):
    """a = area, d = centroid distance from axis"""
    moment = a * d ** 2
    return moment


def get_centroid(A, X, Y):
    """A, X, and Y are lists of areas, and centroid coordinates for i parts"""
    A_sum, AX_sum, AY_sum = 0, 0, 0
    for i in range(0, len(A)):
        A_sum = A_sum + A[i]
        AX_sum = AX_sum + X[i] * A[i]
        AY_sum = AY_sum + Y[i] * A[i]
    X_bar = AX_sum / A_sum
    Y_bar = AY_sum / A_sum
    centroid = [X_bar, Y_bar]
    return centroid


def get_possible_max_bending_coordinates(centroid, B, H, X, Y):
    options = np.zeros((len(B)*4, 2))
    for i in range(0, len(B)):
        # corner 1
        options[i*4][0] = X[i] + B[i]/2 - centroid[0]
        options[i*4][1] = Y[i] + H[i]/2 - centroid[1]

        # corner 2
        options[i*4+1][0] = X[i] + B[i]/2 - centroid[0]
        options[i*4+1][1] = Y[i] - H[i]/2 - centroid[1]

        # corner 3
        options[i*4+2][0] = X[i] - B[i]/2 - centroid[0]
        options[i*4+2][1] = Y[i] + H[i]/2 - centroid[1]

        # corner 4
        options[i*4+3][0] = X[i] - B[i]/2 - centroid[0]
        options[i*4+3][1] = Y[i] - H[i]/2 - centroid[1]
    return options


def get_inertial_moments(centroid, B, H, X, Y):
    Ixx, Iyy, Ixy = 0, 0, 0
    for i in range(0, len(B)):
        X_bar = X[i] - centroid[0]
        Y_bar = Y[i] - centroid[1]
        local_Ixx, local_Iyy, local_Ixy = rectangle_moment_inertia(B[i], H[i])
        Ixx = Ixx + local_Ixx + parallel_axis_theorem(B[i] * H[i], Y_bar)
        Iyy = Iyy + local_Iyy + parallel_axis_theorem(B[i] * H[i], X_bar)
        Ixy = Ixy + B[i] * H[i] * X_bar * Y_bar
    return Ixx, Iyy, Ixy


def get_bending_moments(M, angle):
    angle = np.deg2rad(angle)
    Mx = M * np.cos(angle)
    My = -M * np.sin(angle)
    return Mx, My


def get_sigma_z(Mx, My, Ixx, Iyy, Ixy):
    x_coef = (My * Ixx - Mx * Ixy) / (Ixx * Iyy - Ixy ** 2)
    y_coef = (Mx * Iyy - My * Ixy) / (Ixx * Iyy - Ixy ** 2)
    sigma_z = [x_coef, y_coef]
    return sigma_z


def get_sigma_z_max(sigma_z, coordinates):
    sigma_z_max, x_max, y_max = 0, 0, 0
    for i in range(0, int(np.size(coordinates)/2)):
        new = abs(sigma_z[0] * coordinates[i, 0] + sigma_z[1] * coordinates[i, 1])
        if abs(new) > sigma_z_max:
            sigma_z_max = new
            x_max = coordinates[i, 0]
            y_max = coordinates[i, 1]
    return sigma_z_max, x_max, y_max


if __name__ == '__main__':

    # unsymmetric section
    M, angle = 70.7 * 10 ** 3, -45
    A = [16, 16]
    X = [1, 4]
    Y = [6, 1]
    B = [2, 8]
    H = [8, 2]
    centroid = get_centroid(A, X, Y)
    options = get_possible_max_bending_coordinates(centroid, B, H, X, Y)
    Ixx, Iyy, Ixy = get_inertial_moments(centroid, B, H, X, Y)
    Mx, My = get_bending_moments(M, angle)
    sigma_z = get_sigma_z(Mx, My, Ixx, Iyy, Ixy)
    sigma_z_max, x_max, y_max = get_sigma_z_max(sigma_z, options)
    print(f"Centroid: {centroid[0], centroid[1]}")
    print(f"Ixx = {round(Ixx, 2)}, Iyy = {round(Iyy,2)}, Ixy = {round(Ixy,2)}")
    print(f"Mx = {round(Mx)}, My = {round(My)}")
    print(f"sigma_z = {round(sigma_z[0], 2)}x + {round(sigma_z[1], 2)}y")
    print(f"sigma_z_max = {sigma_z_max}")
    print(f"Location of max bending: {x_max,y_max}")