# Using SPR distribution data, integrate numerically over surface of detector
# Note: does not account for horizontal dist. of beam


import json
import numpy as np


def solid_angle(theta):
    d_theta = np.pi / 500
    d_phi = np.pi / 30 / 500
    return np.abs(np.sin(theta + d_theta / 2) - np.sin(theta - d_theta / 2)) * d_phi


def is_point_on_square(max_dist, theta_square, theta_point): # don't worry about phi it should always be fine
  return np.abs(np.tan(theta_square - theta_point)) <= max_dist


if __name__ == '__main__':

    # print(solid_angle(np.pi / 4))

    filepath = input('Enter data file path: ')

    with open(filepath, 'r') as f:
        data = json.loads(f.read())

    theta_square = float(input('Enter theta of center of detector: ')) * np.pi / 180
    detector_center_phi = 0
    max_dist = float(input('Enter ratio of detector width and dist from grating to detector: '))

    phis, thetas = np.meshgrid(np.linspace(-np.pi / 60, np.pi / 60, 501), np.linspace(0, np.pi, 501))

    result = 0

    for i in range(501):
        for j in range(501):
            theta_point = i * np.pi / 500

            if is_point_on_square(max_dist, theta_square, theta_point):
              result += data[i][j] * solid_angle(theta_point)

    print(f'Expected photons detected per second per electron: {result}')
