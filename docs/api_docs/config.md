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
    cache_dir='.quantized/cache',
    joblib_verbosity=0
)
```
Class to hold configuration for quantized

#### Fields

 **harmonic_oscillator_max_n (int)** Limit of n for the harmonic oscillator

> Default: 50

 **small_number (float)** A value to be used to ensure quantities are nonzero

> Default: 1e-08

 **large_number (float)** A cutoff value that's used to guard against errors, where a number is not normally expected to be large.

> Default: 1000.0

 **float_tol (float)** Two floats are considered equal if they are within this value.

> Default: 1e-06

 **enable_progressbar (bool)** If True, certain methods will print a progress bar to the screen

> Default: False

 **cache_dir (Path)** The directory where cached objects are stored

> Default: .quantized/cache

 **joblib_verbosity (int)** Verbosity level for joblib's cache

> Default: 0











## Functions

----

###to_bool

```
to_bool(x: Union[str, bool]) -> <class 'bool'>
```


Convert x to a bool

For string arguments, "true" and "yes" are considered
Truthy, and are compared in a case-insensitive manner.