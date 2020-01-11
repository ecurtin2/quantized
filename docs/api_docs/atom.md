# quantized.atom Module

This module contains the atom class



## Classes

---

###Atom

```
Atom(element, x, y, z)
```
Atom class containing coordinates, basis and mass.

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


None
#####angle\_to\_xz\_plane
```
Atom.angle_to_xz_plane(self) -> float
```


None
#####angle\_to\_yz\_plane
```
Atom.angle_to_yz_plane(self) -> float
```


None
#####distance
```
Atom.distance(self, other: 'Atom') -> float
```


None
#####flipped\_x
```
Atom.flipped_x(self) -> Atom
```


None
#####rotated
```
Atom.rotated(self, r: 'np.ndarray') -> Atom
```


None
#####rotated\_about\_x
```
Atom.rotated_about_x(self, angle: 'float') -> 'Atom'
```


None
#####rotated\_about\_y
```
Atom.rotated_about_y(self, angle: 'float') -> 'Atom'
```


None
#####rotated\_about\_z
```
Atom.rotated_about_z(self, angle: 'float') -> 'Atom'
```


None
#####rotation\_matrix\_to
```
Atom.rotation_matrix_to(
    self,
    x: 'float',
    y: 'float',
    z: 'float'
) -> <class 'inspect._empty'>
```


https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
#####scaled
```
Atom.scaled(self, factor: 'float') -> Atom
```


None
#####translated
```
Atom.translated(self, x=0.0, y=0.0, z=0.0) -> Atom
```


None
#####with\_coords
```
Atom.with_coords(self, x: 'float', y: 'float', z: 'float') -> 'Atom'
```


None



####Properties

#####coords


Three dimensional array of coordinates, [x, y, z]
#####mass


The mass of the atom
#####normalized\_coords


None







## Functions

----

