# quantized.molecule Module



## Classes

---

###Molecule

```
Molecule(atoms: List[quantized.atom.Atom])
```
#### Fields

 **atoms (Atom)** 

####Methods

#####com\_as\_origin
```
Molecule.com_as_origin(self) -> Molecule
```


None
#####flipped\_x
```
Molecule.flipped_x(self) -> Molecule
```


None
#####map
```
Molecule.map(self, f: Callable) -> Molecule
```


None
#####rotated
```
Molecule.rotated(self, r: <built-in function array>) -> Molecule
```


None
#####rotated\_about\_x
```
Molecule.rotated_about_x(self, angle: float) -> Molecule
```


None
#####rotated\_about\_y
```
Molecule.rotated_about_y(self, angle: float) -> Molecule
```


None
#####rotated\_about\_z
```
Molecule.rotated_about_z(self, angle: float) -> Molecule
```


None
#####scaled
```
Molecule.scaled(self, factor: float) -> Molecule
```


None
#####sorted
```
Molecule.sorted(self, atomic_key: Callable) -> Molecule
```


None
#####translated
```
Molecule.translated(
    self,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0
) -> Molecule
```


None
#####with\_atom\_aligned\_to
```
Molecule.with_atom_aligned_to(
    self,
    atom: quantized.atom.Atom,
    x: float,
    y: float,
    z: float
) -> Molecule
```


None



####Properties

#####R


Calculate the distances for each atom-atom pair.
#####center\_of\_mass


Determine the center of mass of the molecule.
#####coords


None
#####mass


None



####Static Methods

#####from\_xyz
```
Molecule.from_xyz(xyz: str) -> Molecule
```


Create a molecule from an xyz-file formatted string





## Functions

----

