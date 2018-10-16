# coding: utf-8

import stk
import rdkit, rdkit.Chem as rdkit
import itertools

class Combiner:

    """
    Base Skeleton class. Once initialized, can be used to combine a skeleton
    combinatorially with substituents and write the resulting substituted
    molecules in *.mol file format.
    """

    def __init__(self, skeleton, substituents, nmax=None, nconnect=1,
                 connect_atom='Br', auto_placement=False):

        self.skeleton_smiles = skeleton
        self.substituents = substituents
        self.nmax = nmax
        self.nconnect = nconnect
        self.connect_atom = connect_atom
        self.auto_placement = auto_placement


    def combine_substituents(self):

        """
        Get all possible unique structures formed from the skeleton &
        each of the substituents, retaining two Br functional groups
        to be used later when forming supramolecular structure.

        Arguments
        ---------

        skeleton : `str` SMILES string describing skeleton onto which
            substituents will be places. See README.md for how
            skeleton SMILES should be specified.

        Returns
        -------

        canonical_smiles : `list` All unique SMILES obtained
             by placing maximum of two substituents on a molecular
             skeleton.

        """

        smiles = []
        if self.auto_placement:
            template = self.skeleton_smiles.replace('c', 'c{}')
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

        unique_smiles = []
        for smiles in self.get_substituent_permutations(template):
            if smiles not in unique_smiles:
                unique_smiles.append(smiles)
        unique_smiles = sorted(unique_smiles)

        print('Skeleton:', self.skeleton_smiles)
        print('Number of vacant sites:', self.vacant_sites)
        print('Numer of unique substituent permutations:', len(unique_smiles), '\n')

        self.unique_combinations = len(unique_smiles)

        return unique_smiles


    def get_substituent_permutations(self, template):

        """
        Finds all combinations of user-specified substituents. A maximum
        of two substituents may be selected at once.

        Arguments
        ---------

        substituents : `list` user-specified SMILES strings representing substituents
            to be combined with molecular skeleton.

        vacant_sites : `int` numer of sites that can be substituted onto for a given
            molecular skeleton.

        Returns
        -------

        permutations : `list` All permutations of two substituents plus
            two bromine functional groups to be used to build supramolecules

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
                for permutation in list(itertools.permutations(combination)):
                    smiles = rdkit.MolToSmiles(rdkit.MolFromSmiles(template.format(*permutation)), canonical=True)
                    yield smiles


    def assign_ring_order(self, sub_combinations):

        """
        Numerically labels 'opening' and 'closing' of aromatic rings to facilitate
        canonicalisation of smiles strings. Since any rings in the skeleton will
        already be numerically labelled, the number of rings in the skeleton is
        counted and substituent rings are labelled subsequently.
        """

        for combination in sub_combinations:
            m = rdkit.MolFromSmiles(self.skeleton_smiles)
            ring_num = 1+ m.GetRingInfo().NumRings()
            for index, smiles in enumerate(combination):
                combination[index] = smiles.replace('x', str(ring_num))
                ring_num += 1

        return sub_combinations


    def get_embedded_structures(permutations):

        """
        Embeds structures and writes them in *.mol format
        """

        pass



class SpecificationError(Exception):
    def __init__(self, message):
        self.message = message
