# quantized.config Module



## Classes

---

###Config

```
Config(
    harmonic_oscillator_max_n=50,
    small_number=1e-08,
    large_number=1000.0,
    float_tol=1e-06,
    enable_progressbar=False,
    cache_dir='~/.quantized/cache',
    joblib_verbosity=0
)
```
#### Fields

 **harmonic_oscillator_max_n (int)** 

> Default: 50

 **small_number (float)** 

> Default: 1e-08

 **large_number (float)** 

> Default: 1000.0

 **float_tol (float)** 

> Default: 1e-06

 **enable_progressbar (bool)** 

> Default: False

 **cache_dir (Path)** 

> Default: ~/.quantized/cache

 **joblib_verbosity (int)** 

> Default: 0









## Functions

----

###to_bool

```
to_bool(x: Union[str, bool]) -> <class 'bool'>
```


None