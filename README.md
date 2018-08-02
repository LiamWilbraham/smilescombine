## `smilesenumerator` :computer:
Python library for combining smiles 'skeletons' with functional group substituents in a combinatorial fashion.

`smilesenumerator` provides molecular building blocks (in the form of SMILES strings) that can be used in combination 
supramolecular structure-building software such as `stk`. Alternatively, the bulding blocks can be used to generate molecular
libraries in their own right.

## Functionality
Using `polyhts` begins by initialising a class `Skeleton`, which accepts a SMILES string of the molecular skeleton we will
combine with functional groups and a list of those functional groups (`substituents`, represented by SMILES).

```python
substituents = ['(N(C)C)', '(N)', '(OC)', '(O)', '(S)', '(C)', '(F)', '(Cl)', '(CC)', '(C=O)', '(C(=O)OC)']

skeleton = Skeleton('c1(Br)c(Br)c(Br)c(Br)c(Br)c1(Br)', substituents)
```

We can then use the `combine_substituents()` method of `Skeleton` to generate our substituted molecules. 
```python
combinations = skeleton.combine_substituents()
```
Now, the object `combinations` is a list of SMILES strings containing all possible positional and compositonal
combinations of our molecular skeleton and the substituent SMILES supplied in `substituents`.

#### Defining a molecular skeleton SMILES
As shown above, the first argument to the class `Skeleton` is our molecular skeleton of interest. Acceptable
skeleton SMILES contain dummy atoms (in this case `(Br)`) in locations where functional groups may be substituted.



