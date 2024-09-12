import numpy as np

def generate_fibonacci_hemisphere(num_points):
    """Generates points on a hemisphere using the Fibonacci lattice method.

    Args:
        num_points (int): The number of points to generate on the hemisphere.

    Returns:
        numpy.ndarray: An array of shape (num_points, 3) containing the (x, y, z) coordinates of the points on the hemisphere.
    """

    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians

    for i in range(num_points):
        y = 1 - (i / float(num_points - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)             # radius at y

        theta = phi * i                         # Golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        # Since it's a hemisphere, ensure y >= 0
        if y < 0:
            y = -y
            x = -x
            z = -z

        points.append((x, y, z))  # Append as tuple

    return points

