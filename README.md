  ## `smilescombine` :computer:
Python library for combining smiles 'skeletons' with functional groups &
substituents in a combinatorial fashion.

`smilescombine` provides molecular building blocks (in the form of SMILES strings)
that can be used in combination with Supramolecular structure-building software such
as `stk`. Alternatively, the bulding blocks can be used to generate
structured molecular libraries in their own right.

## Functionality
Using `smilescombine` begins by initialising a class `Combiner`, which accepts a
SMILES string of the molecular skeleton we will combine with functional groups
and a list of those functional groups (`substituents`, represented by SMILES).

```python
substituents = ['(N(C)C)', '(N)', '(OC)', '(O)', '(CC)', '(C=O)', '(C(=O)OC)']

skeleton = Combiner('c1ccccc1', substituents, nmax=4, nconnect=0, auto_placement=True)
skeleton.combine_substituents()
```

We can then use the `combine_substituents()` method of `Combiner` to generate our
molecular library for a this skeleton. With the above arguments, `combine_substituents()`
will allow a maximum of 4 substitutions and place them automatically on all
accessible aromatic carbon atoms within the skeleton.

Now, the attribute `combinations` of the instance of the Combiner class `skeleton`
is a list of SMILES strings containing all possible positional and compositonal
combinations of our molecular skeleton and the substituent SMILES supplied by
`substituents`.

#### Defining a molecular skeleton SMILES
If automatic placement of substituents is accectable, we can define an ordinary
SMILES string. If we wish to specify where substitutions may take place, we can
specify `auto_placement=False` within `Combiner` and supply a custom skeleton
SMILES string:

```python
skeleton = Combiner('c1c(Br)cc(Br)cc1', substituents, nmax=4, nconnect=0, auto_placement=True)
```
Where allowed substitution sites are indicated by '(Br)' within the string.

## Installation & Requirements

Simply clone the repository

#### rdkit
smilescombine relies on rdkit, which can be installed via conda (recommended)
```
conda install -c rdkit rdkit
```
