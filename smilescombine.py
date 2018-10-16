# coding: utf-8

import rdkit, rdkit.Chem as rdkit
import itertools

class Combiner:

    """
    Base Skeleton class. Once initialized, can be used to combine a skeleton
    combinatorially with substituents and write the resulting substituted
    molecules in *.mol file format.
    """

    def __init__(self, skeleton, substituents, nmax=None, nconnect=2, autoplacement=False, connect_atom='Br'):

        self.skeleton_smiles = skeleton
        self.substituents = [''] + substituents
        self.nmax = nmax
        self.nconnect = nconnect
        self.autoplacement = autoplacement
        self.connect_atom = connect_type


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
        template = self.skeleton_smiles.replace(self.connect_atom, '{}')
        self.vacant_sites = template.count('{}')

        if self.nconnect > self.vacant_sites:
            raise(SpecificationError
                "Number of connections cannot be greater than the number of vacant sites.")
        elif self.nconnect > self.nmax:
            raise(SpecificationError
                "Number of connections cannot be greater than maximum number of allowed substitutions.")

        perms = self.get_substituent_permutations()
        smiles = [template.format(*perm) for perm in perms]
        canonical_smiles = [rdkit.MolToSmiles(rdkit.MolFromSmiles(smi), canonical=True) for smi in smiles]
        canonical_smiles = remove_duplicates(canonical_smiles)

        print('Skeleton:', self.skeleton_smiles)
        print('Number of vacant sites:', self.vacant_sites)
        print('Numer of unique substituent permutations:', len(canonical_smiles), '\n')

        self.unique_combinations = len(canonical_smiles)

        return canonical_smiles


    def get_substituent_permutations(self):

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
            if self.nmax >= self.vacant_sites:
                vacancies = self.vacant_sites - self.nconnect
            else:
                vacancies = self.nmax - self.nconnect
        else:
            vacancies = self.vacant_sites - self.nconnect

        permutations = []
        for i in range(vacancies+1):
            combinations = itertools.product(substituents, repeat=i)

            combinations = [(list(i) + ['Br']*nconnect) for i in combinations]
            #combinations = self.assign_ring_order(combinations)

            for combination in combinations:
                for permuation in list(itertools.permutations(
                    combination+['']*(self.vacant_sites - len(combination)), self.vacant_sites)):
                    permutations.append(tuple(permuation))

        permutations = sorted(set(permutations))
        return permutations


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
