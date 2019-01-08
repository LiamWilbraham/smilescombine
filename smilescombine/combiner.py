# coding: utf-8

import stk
import rdkit, rdkit.Chem as rdkit
from rdkit.Chem import rdMolDescriptors
import itertools

class SpecificationError(Exception):
    def __init__(self, message):
        self.message = message


class Combiner:
    """
    Base Combiner class. Used to combine an aromatic molecular skeleton
    combinatorially with a list of substituents. Substituents are placed
    on arbitrary molecular skeletons. Substitution sites may either be
    specified randomly or identified automatically (i.e. all available aromatic
    carbon atoms).

    Attributes
    ----------

    skeleton : :class:`str`
        SMILES string representing molecular skeleton onto which substituent
        groups will be placed. Two methods can be used to supply skeleton
        SMILES:
            1. An ordinary SMILES string (e.g. c1ccccc1) if smilescombine is
            to automatically place substituents on available aromatic carbons.

            2. A modified SMILES string where substitution positions can be
            specified manually, for example, indicated by Bromine atoms
            (e.g. c1c(Br)cc(Br)cc1). The atom type indicating allowed
            substitution positions should be the same as connect_atom (below).

    substituents : :class:`list`
        A list of allowed substituents. Each substituent should be enclosed
        within parentheses (e.g. ['(N(C)C)', '(N)', '(OC)', '(O)' '(S)']).

    nmax : :class:`int` (default = ``2``)
        Maxumum number of substitutions allowed for a given skeleton.

    nconnect : :class:`int` (default = ``2``)
        Number of 'connecting' atoms to be left on each substituted skeleton.
        This is indended to be used in conjunction with STK to produce libraries
        of supramolecules (https://github.com/lukasturcani/stk).

    connect_atom : :class:`str` (default = ``Br``)
        Atom type to be used both to indicate user-defined substitution positions
        and for use as connection points for the definition of STK StructUnit
        objects.

    auto_placement : :class:`bool` (default = ``True``)
        Specified whether substitution positions are decided automatically or
        by the user.

    Combiner can be used to produce standalone substituted molecules, but is
    also intended for use in conjunction with STK (https://github.com/lukasturcani/stk)
    to construct structured libraries of supramolecules. Therefore, the number of
    required connection atoms and their labels (required by STK) may also be
    specified.
    """

    def __init__(self, skeleton, substituents, nmax=2, nconnect=0,
                 connect_atom='Br', auto_placement=True):

        self.skeleton_smiles = skeleton
        self.substituents = self.assign_ring_order(skeleton, substituents)
        self.nmax = nmax
        self.nconnect = nconnect
        self.connect_atom = connect_atom
        self.auto_placement = auto_placement
        self.combinations = []


    def combine_substituents(self, filename):

        """
        Generates all possible unique structures formed from the skeleton &
        each of the substituents.

        Arguments
        ---------

        filename : :class:`str` Path to and name of file into which
        SMILES combinations will be written.
        """

        template = self.get_skeleton_template()

        for smiles in self.get_substituent_permutations(template):
            if smiles not in self.combinations:
                self.combinations.append(smiles)

        self.combinations = sorted(self.combinations, reverse=True)
        self.n_combinations = len(self.combinations)

        print('Skeleton SMILES:', self.skeleton_smiles)
        print('Number of vacant sites:', self.vacant_sites)
        print('Numer of unique substituent permutations:', self.n_combinations, '\n')

        self.write_smiles(filename, self.combinations)


    def get_substituent_permutations(self, template):

        """
        Generator that yields all combinations of user-specified substituents.

        Arguments
        ---------

        template : :class:`str` SMILES string representing aromatic skeleton and
            available substitution sites.

        Yields
        -------

        smiles : :class:`str` permutation of a given subset of substituents.

        """

        if self.nmax is not None:
            if self.nmax >=self.vacant_sites:
                vacancies = self.vacant_sites - self.nconnect
            else:
                vacancies = self.nmax - self.nconnect
        else:
            vacancies = self.vacant_sites - self.nconnect

        for i in range(vacancies+1):
            combinations = itertools.product(self.substituents, repeat=i)

            combinations = [(list(i) + ['('+self.connect_atom+')']*self.nconnect) for i in combinations]
            combinations = [i+['']*(self.vacant_sites - len(i)) for i in combinations]
            #combinations = self.assign_ring_order(combinations)

            for combination in combinations:
                for permutation in set(itertools.permutations(combination)):
                    smiles = rdkit.MolToSmiles(rdkit.MolFromSmiles(template.format(*permutation)), canonical=True)
                    yield smiles


    def get_skeleton_template(self):

        """
        Converts skeleton SMILES string into template where possible
        substitution sites are identified.

        Returns
        -------

        template : :class:`str`
            Pseudo-SMILES string with possible substitution sites indicated
            by parentheses '{}'.
        """

        if self.auto_placement:
            mol_h = rdkit.MolFromSmiles(self.skeleton_smiles)
            rdkit.AddHs(mol_h)
            template = rdkit.MolToSmiles(mol_h, allHsExplicit=True)
            template = template.replace('[cH]', 'c{}').replace('[c]', 'c')
        else:
            template = self.skeleton_smiles.replace('(Br)', '{}')

        self.vacant_sites = template.count('{}')

        if self.nconnect > self.vacant_sites:
            raise SpecificationError(
                "Number of connections cannot be greater than the number of possible substitution sites.")
        if self.nmax is not None:
            if self.nconnect > self.nmax:
                raise SpecificationError(
                    "Number of connections cannot be greater than the maximum number of allowed substitutions.")

        return template


    def assign_ring_order(self, skeleton, substituents):

        """
        Assures that ring numbering in substituents is compatible with the
        number of rings present in the skeleton.

        Arguments
        ---------

        skeleton : :class:`str`
            SMILES string representing molecular skeleton onto which substituent
            groups will be placed.

        substituents : :class:`list`
            A list of allowed substituents, represented by SMILES strings.

        Returns
        -------

        substituents : :class:`list`
            The list of allowed substituents, still represented by SMILES
            strings, with their ring open/close numbering adjusted to be
            compatible with the number of rings present in the skeleton.

        """

        n = rdMolDescriptors.CalcNumAromaticRings(rdkit.MolFromSmiles(skeleton))

        for i, item in enumerate(substituents):
            rings = rdMolDescriptors.CalcNumAromaticRings(rdkit.MolFromSmiles(item[1:-1]))
            if rings > 0:
                for j in reversed(range(rings+1)):
                    item = item.replace(str(j), str(j+n))
                substituents[i] = item

        return substituents


    def write_smiles(self, filename, combinations):

        """
        Writes SMILES to *.csv file
        """

        with open(filename, 'w') as f:
            for smi in combinations:
                f.write(smi + '\n')


    def __str__(self):
        string = 'Skeleton SMILES: ' + self.skeleton_smiles + '\n'
        string += 'Substituents: ' + str(self.substituents) + '\n'
        string += 'Max number of substitutions: ' + str(self.nmax) + '\n'
        string += 'Possible substitution sites: ' + str(self.vacant_sites) + '\n'
        string += 'Number of unique combinations: ' + str(len(self.combinations)) + '\n'
        if self.nconnect > 0:
            string += 'Connection points: ' + str(self.nconnect) + '\n'
        return string


    def __repr__(self):
        return str(self)
