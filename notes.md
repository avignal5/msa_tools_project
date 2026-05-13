# Notes on the msa_tools module
## \_\_init\_\_.py
* This file turns the directory msa_tools into a Python package. 
```python
from .alignment import MultipleSequenceAlignment
```
* the . means relative import: import from the same package directory.
* \_\_all\_\_ : list of what will be publicly available when:
```python
from mypackage import *
```
## io.py
* for reading a multiple alignment
* keeping parsing separate from the class will allow to add other formats : read_clustal(), etc.
## alignment.py
### def \_\_init\_\_
* Defines a class MultipleSequenceAlignment
```python
def __init__(self, sequences: Dict[str, str], reference: str = None):
```
* Automatically called when creating a new object.
  * self : Refers to the instance being created and Allows to store data on the object
  * sequences: Dict[str, str] : sequences should be a dictionary, where keys and values are strings
  * reference: str = None : optional, default = None, string if provided
### def \_\repr\_\_(self)
```python
def __repr__(self):
    return f"<MSA: {len(self.sequences)} sequences, length={self.length}>"
```
* \_\_repr\_\_ controls what is shown when print in the interpreter or inspect in debugger or type the object name in notebook
### def _build_coordinate_map(self):
*ls -l


 The leading underscore signals that the function (or variable) is intended to be "private" : not accessible directly.


