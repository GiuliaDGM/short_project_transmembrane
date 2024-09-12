import numpy as np
from numbers import Number
import math


class Vector3D:
    ##float^3->Vector3D
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    
    ##Vector3D->str
    def __str__(self):
        return "[{:.2f}, {:.2f}, {:.2f}]".format(self.x, self.y, self.z)
    
    ##Vector3D->str
    def __repr__(self):
        return self.__str__()
    
    ##Vector3D->Vector3D
    def __add__(self, other):
        if isinstance(other, Vector3D):
            return Vector3D(self.x + other.x, self.y + other.y, self.z + other.z)
        elif isinstance(other, Number):
            return Vector3D(self.x + other, self.y + other, self.z + other)
        else:
            raise TypeError("Unsupported type for addition")
    
    ##Vector3D->Vector3D
    def __sub__(self, other):
        if isinstance(other, Vector3D):
            return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)
        elif isinstance(other, Number):
            return Vector3D(self.x - other, self.y - other, self.z - other)
        else:
            raise TypeError("Unsupported type for subtraction")
    
    ##Vector3D->Vector3D
    def __mul__(self, value):
        if isinstance(value, Vector3D):
            return Vector3D(self.x * value.x, self.y * value.y, self.z * value.z)
        elif isinstance(value, Number):
            return Vector3D(self.x * value, self.y * value, self.z * value)
        else:
            raise TypeError("Unsupported type for multiplication")
    
    ##Vector3D->Vector3D
    def __rmul__(self, value):
        return self.__mul__(value)
    
    ##Vector3D->Vector3D
    def __truediv__(self, value):
        if isinstance(value, Vector3D):
            return Vector3D(self.x / value.x, self.y / value.y, self.z / value.z)
        elif isinstance(value, Number):
            return Vector3D(self.x / value, self.y / value, self.z / value)
        else:
            raise TypeError("Unsupported type for division")
    
    ##Vector3D->Vector3D
    def __neg__(self):
        return Vector3D(-self.x, -self.y, -self.z)
    
    ##Vector3D->float
    def norm(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def normalized(self):
        """
        Return a normalized (unit length) version of the vector.
        
        Returns:
            Vector3D: A new Vector3D instance with the same direction but unit length.
        
        Raises:
            ValueError: If the vector is zero-length.
        """
        norm = self.norm()
        if norm == 0:
            raise ValueError("Cannot normalize a zero vector.")
        return Vector3D(self.x / norm, self.y / norm, self.z / norm)

    
    def dist_to_plane(self, vector_plane_normal, d=0):
        """
        Calculate the distance from this vector (point) to the plane defined by its normal vector.
        Helps in determining which residues are included in each slice.
        (Positive or negative distance from membrane plane = Inside or Outside)
        
        Args:
            vector_plane_normal (Vector3D): The normal vector (A, B, C) of the plane.
            d (float): The D value in the plane equation Ax + By + Cz + D = 0.
                    Default is 0 if the plane passes through the origin.
        
        Returns:
            float: The distance from the point to the plane.
        """
        if not isinstance(vector_plane_normal, Vector3D):
            raise TypeError("plane_normal must be a Vector3D instance")
        
        # Calculate the dot product of the point (self) with the plane normal vector (A, B, C)
        numerator = np.abs(vector_plane_normal.x * self.x + vector_plane_normal.y * self.y + vector_plane_normal.z * self.z + d)
        
        # Calculate the magnitude of the plane normal vector (sqrt(A^2 + B^2 + C^2))
        denominator = np.sqrt(vector_plane_normal.x ** 2 + vector_plane_normal.y ** 2 + vector_plane_normal.z ** 2)
        
        # Return the distance
        return numerator / denominator
