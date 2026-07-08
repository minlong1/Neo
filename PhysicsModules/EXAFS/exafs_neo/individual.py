from larch.xafs import feffdat

import numpy as np

from Solvers.core import GeneRange, ParameterSpace
from Solvers.core import Individual as SolverIndividual

from typing import Dict, List, AnyStr

"""
Construct individuals for the GA.

The EXAFS individual is a Solvers individual over a flat genome:

    genes = [e0, s02_0, sigma2_0, deltaR_0, s02_1, sigma2_1, deltaR_1, ...]

with one shared e0 gene followed by three genes per FEFF path. The
path-oriented accessors (get_path/set_path/get_e0/...) preserve the
historical EXAFS API on top of that genome.
"""


def build_parameter_space(pathrange_pars: List) -> ParameterSpace:
    """Build the flat-genome ParameterSpace from per-path Pathrange_limits."""
    gene_ranges = [GeneRange(pathrange_pars[0].getrange_E0(), name="e0")]
    for i, pathrange in enumerate(pathrange_pars):
        gene_ranges.append(GeneRange(pathrange.getrange_S02(), name=f"path{i}_s02"))
        gene_ranges.append(GeneRange(pathrange.getrange_Sigma2(), name=f"path{i}_sigma2"))
        gene_ranges.append(GeneRange(pathrange.getrange_DeltaR(), name=f"path{i}_deltaR"))
    return ParameterSpace(gene_ranges)


class Individual(SolverIndividual):
    """A Solvers Individual over the shared-e0 genome (see module
    docstring); adds get_path/set_path/get_e0/... accessors so EXAFS
    fitness/path-optimization code can keep addressing genes by path
    instead of by raw index."""

    def __init__(
        self,
        npaths: int,
        pathDictionary: Dict,
        pathrange_Dict: Dict,
        pathlists: List,
        e0: float,
        pathName: AnyStr,
    ):
        """

        :param npaths:
        :param pathDictionary:
        :param pathrange_Dict:
        :param pathlists:
        :param e0:
        :param pathName:
        """

        self.npaths = npaths
        self.path_lists = pathlists
        self.pathname = pathName
        self.pathrange_Dict = pathrange_Dict
        self.pathDictionary = pathDictionary

        space = build_parameter_space(list(pathrange_Dict)[:npaths])
        genes = space.sample()
        genes[0] = e0
        super().__init__(space, genes)

    def _path_slice(self, i: int) -> slice:
        return slice(1 + 3 * i, 4 + 3 * i)

    def get(self) -> List:
        """Get all vars

        Returns:
            list of (npaths,4) 2D list
        """
        return [self.get_path(i) for i in range(self.npaths)]

    def get_var(self) -> List:
        """Get all parameters except e0

        Returns:
            list of (npaths,3) 2D list
        """
        return [list(self.genes[self._path_slice(i)]) for i in range(self.npaths)]

    def get_e0(self) -> float:
        """Get e0 of the individual
        Returns:
            int: e0 value
        """
        return float(self.genes[0])

    def get_path(self, i):
        s02, sigma2, deltaR = self.genes[self._path_slice(i)]
        return [s02, self.get_e0(), sigma2, deltaR]

    def verbose(self) -> None:
        """
        Print out the populations
        """
        for i in range(self.npaths):
            s02, e0, sigma2, deltaR = self.get_path(i)
            print(s02, e0, sigma2, deltaR)

    def set_path(self, i: int, s02: float, sigma2: float, deltaR: float) -> None:
        self.genes[self._path_slice(i)] = [s02, sigma2, deltaR]

    def set_e0(self, e0: float) -> None:
        self.genes[0] = e0

    def mutate_paths(self, chance: float) -> None:
        """Resample path genes; `chance` keeps the historical percent scale
        (a gene mutates when random()*100 < chance)."""
        for gene_idx in range(1, len(self.genes)):
            if np.random.random() * 100 < chance:
                self.genes[gene_idx] = self.space.sample_gene(gene_idx)

    def mutate(self, chance: float) -> int:
        """Solvers-facing mutation hook; delegates to the EXAFS path mutation
        so generic mutators preserve historical EXAFS behavior (shared e0
        gene untouched, percent-scaled chance)."""
        self.mutate_paths(chance)
        return 0

    def verbose_yTotal(self, interval_k: List) -> List:
        yTotal = [0] * 401
        for i in range(self.npaths):
            path = self.pathDictionary.get(self.pathname[i])
            path_pars = self.get_path(i)
            path.s02 = path_pars[0]
            path.e0 = path_pars[1]
            path.sigma2 = path_pars[2]
            path.deltar = path_pars[3]
            feffdat.path2chi(path)
            y = path.chi
            for k in interval_k:
                yTotal[int(k)] += y[int(k)]

        return yTotal

    def __len__(self) -> int:
        """Return the number of independent parameters in the individual."""

        return int((3 * self.npaths) + 1)

    def __str__(self) -> AnyStr:
        return f"Individual with values of {self.get()}"
