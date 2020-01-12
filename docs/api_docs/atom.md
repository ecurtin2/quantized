# quantized.atom Module

This module contains the atom class



## Classes

---

###Atom

```
Atom(element, x, y, z)
```
Atom class containing coordinates, basis and mass.
The atom will generally not be used in isolation, but will
likely be part of a molecule. This is expected to be used
as a structured container for atomic information, but
it does contain logic to transform atoms in space.


#### Fields

 **element (Element)** The element

 **x (float)** The x coordinate of the atom

 **y (float)** The y coordinate of the atom

 **z (float)** The z coordinate of the atom

####Methods

#####angle\_to\_xy\_plane
```
Atom.angle_to_xy_plane(self) -> float
```


Angle, in radians, between this atom's coordinate vector and the xy plane
#####angle\_to\_xz\_plane
```
Atom.angle_to_xz_plane(self) -> float
```


Angle, in radians, between this atom's coordinate vector and the xz plane
#####angle\_to\_yz\_plane
```
Atom.angle_to_yz_plane(self) -> float
```


Angle, in radians, between this atom's coordinate vector and the yz plane
#####distance
```
Atom.distance(self, other: 'Atom') -> float
```


Determine the distance between this atom and another atom.
#####flipped\_x
```
Atom.flipped_x(self) -> Atom
```


Return an equivalent atom, but the x coordinate is the opposite
#####rotated
```
Atom.rotated(self, r: 'np.ndarray') -> Atom
```


Return an equivalent atom rotated by the given rotation matrix

The matrix must have shape (3, 3)
#####rotated\_about\_x
```
Atom.rotated_about_x(self, angle: 'float') -> 'Atom'
```


Return an equivalent atom, rotated by `angle` radians about the x axis
#####rotated\_about\_y
```
Atom.rotated_about_y(self, angle: 'float') -> 'Atom'
```


Return an equivalent atom, rotated by `angle` radians about the y axis
#####rotated\_about\_z
```
Atom.rotated_about_z(self, angle: 'float') -> 'Atom'
```


Return an equivalent atom, rotated by `angle` radians about the z axis
#####rotation\_matrix\_to
```
Atom.rotation_matrix_to(self, x: 'float', y: 'float', z: 'float') -> np.ndarray
```


Get a matrix that would rotate this atom to align with the given coordinates

https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
#####scaled
```
Atom.scaled(self, factor: 'float') -> Atom
```


Return an equivalent atom with all coordinates scaled by some factor
#####translated
```
Atom.translated(self, x=0.0, y=0.0, z=0.0) -> Atom
```


Return an equivalent atom translated in the direction given
#####with\_coords
```
Atom.with_coords(self, x: 'float', y: 'float', z: 'float') -> 'Atom'
```


Return an equivalent atom at these coordinates



####Properties

#####coords


Three dimensional array of coordinates, [x, y, z]
#####mass


The mass of the atom
#####normalized\_coords


Return a unit vector pointing towards the atom







####Dunder Methods

#####\_\_eq\_\_
```
Atom.__eq__(self, other: "'Atom'") -> bool
```


Return True if the other atom is the same element, and very close

