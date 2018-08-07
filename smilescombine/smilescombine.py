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

    def __init__(self, skeleton, substituents):

        self.skeleton_smiles = skeleton
        self.substituents = [''] + substituents


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
        template = self.skeleton_smiles.replace('(Br)', '{}')
        self.vacant_sites = template.count('{}')
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

        if self.vacant_sites >= 4:
            sub_combinations = list(itertools.combinations(self.substituents, 2))
            for i in self.substituents:
                sub_combinations.append([i, i])

        elif self.vacant_sites == 3:
            sub_combinations = list(itertools.combinations(self.substituents, 1))

        else:
            sub_combinations = list(itertools.combinations(self.substituents, 0))

        sub_combinations = [(list(i) + ['(Br)', '(Br)']) for i in sub_combinations]
        sub_combinations = self.assign_ring_order(sub_combinations)

        permutations = []
        for combination in sub_combinations:
            for permuation in list(itertools.permutations(
                combination+['']*(self.vacant_sites - len(combination)), self.vacant_sites)):
                permutations.append(list(permuation))

        permutations = remove_duplicates(permutations)

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


def remove_duplicates(x):
    # simple tool for removing duplicates in lists
    x_unique = []
    for item in x:
        if item not in x_unique:
            x_unique.append(item)
    return x_unique


# #### Notes
# * problems with Si and Carboxylic Acid substituents : cannot canonicalize
