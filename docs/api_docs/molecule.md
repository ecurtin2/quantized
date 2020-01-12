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



#####flipped\_x
```
Molecule.flipped_x(self) -> Molecule
```



#####map
```
Molecule.map(self, f: Callable) -> Molecule
```



#####rotated
```
Molecule.rotated(self, r: <built-in function array>) -> Molecule
```



#####rotated\_about\_x
```
Molecule.rotated_about_x(self, angle: float) -> Molecule
```



#####rotated\_about\_y
```
Molecule.rotated_about_y(self, angle: float) -> Molecule
```



#####rotated\_about\_z
```
Molecule.rotated_about_z(self, angle: float) -> Molecule
```



#####scaled
```
Molecule.scaled(self, factor: float) -> Molecule
```



#####sorted
```
Molecule.sorted(self, atomic_key: Callable) -> Molecule
```



#####translated
```
Molecule.translated(
    self,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0
) -> Molecule
```



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






####Properties

#####R


Calculate the distances for each atom-atom pair.
#####center\_of\_mass


Determine the center of mass of the molecule.
#####coords



#####mass








####Static Methods

#####from\_xyz
```
Molecule.from_xyz(xyz: str) -> Molecule
```


Create a molecule from an xyz-file formatted string



