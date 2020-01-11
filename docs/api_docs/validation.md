# quantized.validation Module



## Classes

---

###Range

```
Range(min, max)
```
Attrs Validator to ensure values between a min and a max value.


#### Fields

 **min (None)** The minimum allowed value

> Constraints:  <function Range.check at 0x7fb2288057b8>

 **max (None)** The maximum allowed value

####Methods

#####check
```
Range.check(self, attribute, value) -> <class 'inspect._empty'>
```


None









## Functions

----

###not_nan

```
not_nan(
    instance: Any,
    attribute: attr._make.Attribute,
    value: float
) -> <class 'inspect._empty'>
```


None

###not_inf

```
not_inf(
    instance: Any,
    attribute: attr._make.Attribute,
    value: float
) -> <class 'inspect._empty'>
```


None