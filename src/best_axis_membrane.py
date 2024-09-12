import numpy as np

# def max_sub_array_sum(arr):
#     """Find the maximum sum of a contiguous subarray.

#     Args:
#         arr: List or numpy array of numbers.

#     Returns:
#         tuple: (start_index, end_index, max_sum) where max_sum is the maximum sum
#                of the subarray found.
#     """
#     n = len(arr)
#     max_sum = float('-inf')
#     start_index = end_index = 0
#     temp_start = 0

#     # Iterate over all possible starting points
#     for i in range(n):
#         current_sum = 0
#         # Iterate over all possible ending points
#         for j in range(i, n):
#             current_sum += arr[j]
#             # Update the maximum sum if the current sum is greater
#             if current_sum > max_sum:
#                 max_sum = current_sum
#                 start_index = i
#                 end_index = j

#     return (start_index, end_index + 1, max_sum)

def max_sub_array_sum(arr):
    """Find the maximum sum of a contiguous subarray.

    Args:
        arr: List or numpy array of numbers.

    Returns:
        tuple: (start_index, end_index, max_sum) where max_sum is the maximum sum
               of the subarray found. Returns (None, None, None) for empty arrays.
    """
    n = len(arr)
    if n == 0:
        return (None, None, None)

    max_sum = float('-inf')
    start_index = end_index = 0

    # Iterate over all possible starting points
    for i in range(n):
        current_sum = 0
        # Iterate over all possible ending points
        for j in range(i, n):
            current_sum += arr[j]
            # Update the maximum sum if the current sum is greater
            if current_sum > max_sum:
                max_sum = current_sum
                start_index = i
                end_index = j

    return (start_index, end_index + 1, max_sum)


def get_best_results(processed_lines):
    """Find the best result among the processed lines.

    Args:
        processed_lines (list): A list of dictionaries where each dictionary represents a processed line
                                and contains hydrophobicity and other related data.

    Returns:
        list: A list containing the best line data, and the start and end indices of the best subarray,
              along with the sum of hydrophobicity and metadata (nb_steps, shortest_distance).
    """
    lines = []
    best_slices = None
    nb_steps = None
    shortest_distance = None

    # Flatten the processed_lines into a single list
    for index in processed_lines:
        # Ensure index is iterable, handle any edge cases
        if isinstance(index, dict):
            lines.append(index)
        else:
            print(f"Warning: Skipping invalid entry: {index}")

    if not lines:
        raise ValueError("No valid lines found in processed_lines.")

    # Find the best line based on the second value in the "axis_average_hydrophobicity" tuple
    best_line = max(lines, key=lambda line: line["axis_average_hydrophobicity"][1])

    # Find the best slices, nb_steps, and shortest_distance based on the best line
    for line in lines:
        if line["axis_average_hydrophobicity"] == best_line["axis_average_hydrophobicity"]:
            best_slices = line["slice_hydrophobicity"]
            nb_steps = line["total_slices"]
            shortest_distance = line["min_distance"]

    # Find the best subarray (slice) of hydrophobicity values
    start_index, best_index, best = max_sub_array_sum(best_slices)

    # Return the best line information along with subarray data and metadata
    return [best_line, start_index, best_index, best, nb_steps, shortest_distance]


# def generate_membranes(processed_lines, best_results, resolution):
#     """Generate points in the space to simulate the membranes.

#     Args:
#         processed_lines: List containing lines with slice infos and hydrophobicity values.
#         best_results: List with results including plane normal, distances, etc.
#         resolution: Integer setting the step of sliding.

#     Returns:
#         tuple: (points_membrane_1, points_membrane_2) as numpy arrays.
#     """
#     # Unpack results
#     plane_normal = best_results[0]["axis_average_hydrophobicity"][0]
#     shortest_distance = best_results[5]
#     start_index = best_results[1]
#     best_index = best_results[2]

#     # Calculate distances for the two membranes
#     dist_m1 = resolution * (start_index + 1)
#     dist_m2 = resolution * (best_index + 1)

#     # Create points for each membrane
#     def create_plane_points(distance):
#         # Define the points on the plane
#         point_on_plane = plane_normal * (shortest_distance + distance)
#         d = - (point_on_plane.x * plane_normal.x + point_on_plane.y * plane_normal.y + point_on_plane.z * plane_normal.z)

#         # Create a grid of points
#         grid_size = 20
#         X, Y = np.meshgrid(np.linspace(-10, 10, grid_size), np.linspace(-10, 10, grid_size))
#         X_flat, Y_flat = X.ravel(), Y.ravel()

#         # Calculate Z coordinates for the grid points
#         Z_flat = (-plane_normal.x * X_flat - plane_normal.y * Y_flat - d) / plane_normal.z

#         # Combine X, Y, Z into an array of points
#         return np.vstack((X_flat, Y_flat, Z_flat)).T

#     # Generate points for both membranes
#     points_membrane_1 = create_plane_points(dist_m1)
#     points_membrane_2 = create_plane_points(dist_m2)

#     print(f"Distance to membrane 1: {dist_m1}")
#     print(f"Distance to membrane 2: {dist_m2}")

#     return points_membrane_1, points_membrane_2

def generate_membranes(processed_lines, best_results, resolution):
    """Generate points in the space to simulate the membranes.

    Args:
        processed_lines: List containing lines with slice infos and hydrophobicity values.
        best_results: List with results including plane normal, distances, etc.
        resolution: Integer setting the step of sliding.

    Returns:
        tuple: (points_membrane_1, points_membrane_2) as numpy arrays.
    """
    # Unpack results
    best_line = best_results[0]
    start_index = best_results[1]
    best_index = best_results[2]
    shortest_distance = best_results[5]
    
    # Extract plane_normal from best_line
    plane_normal = best_line["axis_average_hydrophobicity"][0]
    
    # Calculate distances for the two membranes
    dist_m1 = resolution * start_index  # Adjusted to avoid double +1
    dist_m2 = resolution * best_index    # Adjusted to avoid double +1

    # Create points for each membrane
    def create_plane_points(distance):
        # Define the points on the plane
        point_on_plane = plane_normal * (shortest_distance + distance)
        # Correct plane equation: d = - (n . p)
        d = - (plane_normal.x * point_on_plane.x + plane_normal.y * point_on_plane.y + plane_normal.z * point_on_plane.z)

        # Handle division by zero for Z coordinate
        if plane_normal.z == 0:
            raise ValueError("Plane normal's z-component is zero, leading to division by zero in Z coordinate calculation.")

        # Create a grid of points
        grid_size = 20
        X, Y = np.meshgrid(np.linspace(-10, 10, grid_size), np.linspace(-10, 10, grid_size))
        X_flat, Y_flat = X.ravel(), Y.ravel()

        # Calculate Z coordinates for the grid points
        Z_flat = (-plane_normal.x * X_flat - plane_normal.y * Y_flat - d) / plane_normal.z

        # Combine X, Y, Z into an array of points
        return np.vstack((X_flat, Y_flat, Z_flat)).T

    try:
        # Generate points for both membranes
        points_membrane_1 = create_plane_points(dist_m1)
        points_membrane_2 = create_plane_points(dist_m2)
    except ValueError as e:
        print(f"Error generating membrane points: {e}")
        return None, None

    print(f"Distance to membrane 1: {dist_m1}")
    print(f"Distance to membrane 2: {dist_m2}")

    return points_membrane_1, points_membrane_2
