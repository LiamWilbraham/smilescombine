# coding: utf-8

import stk
import rdkit, rdkit.Chem as rdkit
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

    Combiner can be used to produce standalone substituted molecules, but is
    also intended for use in conjunction with STK (https://github.com/lukasturcani/stk)
    to construct structured libraries of supramolecules. Therefore, the number of
    required connection atoms and their labels (required by STK) may also be specified.
    """

    def __init__(self, skeleton, substituents, nmax=None, nconnect=0,
                 connect_atom='Br', auto_placement=True):

        self.skeleton_smiles = skeleton
        self.substituents = substituents
        self.nmax = nmax
        self.nconnect = nconnect
        self.connect_atom = connect_atom
        self.auto_placement = auto_placement
        self.combinations = []


    def combine_substituents(self):

        """
        Generates all possible unique structures formed from the skeleton &
        each of the substituents.
        """

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

        for smiles in self.get_substituent_permutations(template):
            if smiles not in self.combinations:
                self.combinations.append(smiles)
        self.combinations = sorted(self.combinations, reverse=True)
        self.n_combinations = len(self.combinations)

        print('Skeleton SMILES:', self.skeleton_smiles)
        print('Number of vacant sites:', self.vacant_sites)
        print('Numer of unique substituent permutations:', self.n_combinations, '\n')


    def get_substituent_permutations(self, template):

        """
        Generator that yields all combinations of user-specified substituents.

        Arguments
        ---------

        template : `str` SMILES string representing aromatic skeleton and
            available substitution sites.

        Yields
        -------

        smiles : `str` permutation of a given subset of substituents.

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
