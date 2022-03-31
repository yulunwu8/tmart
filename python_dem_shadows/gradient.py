import numpy as np


def gradient(grid, length_x, length_y=None):
    """
    Calculate the numerical gradient of a matrix in X, Y and Z directions.

    :param grid: Matrix
    :param length_x: Length between two columns
    :param length_y: Length between two rows
    :return:
    """
    if length_y is None:
        length_y = length_x

    assert len(grid.shape) == 2, "Grid should be a matrix."

    grad = np.empty((*grid.shape, 3))
    grad[:] = np.nan
    grad[:-1, :-1, 0] = 0.5 * length_y * (
        grid[:-1, :-1] - grid[:-1, 1:] + grid[1:, :-1] - grid[1:, 1:]
    )
    grad[:-1, :-1, 1] = 0.5 * length_x * (
        grid[:-1, :-1] + grid[:-1, 1:] - grid[1:, :-1] - grid[1:, 1:]
    )
    grad[:-1, :-1, 2] = length_x * length_y

    # Copy last row and column
    grad[-1, :, :] = grad[-2, :, :]
    grad[:, -1, :] = grad[:, -2, :]

    area = np.sqrt(
        grad[:, :, 0] ** 2 +
        grad[:, :, 1] ** 2 +
        grad[:, :, 2] ** 2
    )
    for i in range(3):
        grad[:, :, i] /= area
    return grad


def check_gradient(grad):
    assert len(grad.shape) == 3 and grad.shape[2] == 3, \
        "Gradient should be a tensor with 3 layers."


def aspect(grad, degrees=False):
    """
    Calculate the elevation aspect angle given the gradient.

    Aspect is the direction a slope is facing to.

    :param grad: Tensor representing the X,Y,Z gradient
    :param degrees: Output in degrees or radians
    :return: Matrix with aspect per grid cell.
    """
    check_gradient(grad)

    y_grad = grad[:, :, 1]
    x_grad = grad[:, :, 0]
    asp = np.arctan2(y_grad, x_grad) + (np.pi / 2)
    asp[asp < 0] += (2 * np.pi)

    if degrees:
        asp = np.rad2deg(asp)
    return asp


def slope(grad, degrees=False):
    """
    Calculate the slope inclination angle given the gradient.
    :param grad: Tensor representing the X,Y,Z gradient
    :param degrees:
    :return:
    """
    check_gradient(grad)

    sl = np.arccos(grad[:, :, 2])
    if degrees:
        sl = np.rad2deg(sl)
    return sl


def normal_vector(slope_deg, aspect_deg):
    """
    Calculate the unit vector normal to the surface defined by slope and aspect.
    :param slope_deg: slope inclination in degrees
    :param aspect_deg: slope aspect in degrees
    :return: 3-dim unit normal vector
    """
    slope_rad = np.deg2rad(slope_deg)
    aspect_rad = np.deg2rad(aspect_deg)

    nvx = np.sin(aspect_rad) * np.sin(slope_rad)
    nvy = -np.cos(aspect_rad) * np.sin(slope_rad)
    nvz = np.cos(slope_rad)
    return np.array([nvx, nvy, nvz])


def hill_shade(grad, sun_vector):
    """
    Compute the intensity of illumination on a surface given the sun position.
    :param grad:
    :param sun_vector:
    :return:
    """
    check_gradient(grad)

    hsh = (
        grad[:, :, 0] * sun_vector[0] +
        grad[:, :, 1] * sun_vector[1] +
        grad[:, :, 2] * sun_vector[2]
    )
    # Remove negative incidence angles - indicators for self-shading.
    hsh = (hsh + abs(hsh)) / 2.

    return hsh
