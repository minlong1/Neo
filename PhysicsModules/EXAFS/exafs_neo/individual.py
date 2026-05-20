from larch.xafs import feffdat

from PhysicsModules.EXAFS.exafs_neo.pathObj import PathObject
from typing import Dict, List, AnyStr

"""
Construct individuals for the GA
"""


class Individual:
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
        self.Population = []
        self.pathrange_Dict = pathrange_Dict
        self.pathDictionary = pathDictionary

        for this_path in self.pathrange_Dict:
            self.Population.append(PathObject(this_path, e0))

    def get(self) -> List:
        """Get all vars

        Returns:
            list of (npaths,4) 2D list
        """
        population = []
        for i in range(self.npaths):
            population.append(self.Population[i].get())
        return population

    def get_var(self) -> List:
        """Get all parameters except e0

        Returns:
            list of (npaths,3) 2D list
        """
        population = []
        for i in range(self.npaths):
            population.append(self.Population[i].get_var())
        return population

    def get_e0(self) -> float:
        """Get e0 of the individual
        Returns:
            int: e0 value
        """
        return self.Population[0].get_e0()

    def get_path(self, i):
        return self.Population[i].get()

    def verbose(self) -> None:
        """
        Print out the populations
        """
        for i in range(self.npaths):
            self.Population[i].verbose()

    def set_path(self, i: int, s02: float, sigma2: float, deltaR: float) -> None:
        self.Population[i].set(s02, sigma2, deltaR)

    def set_e0(self, e0: float) -> None:
        for i in range(self.npaths):
            self.Population[i].set_e0(e0)

    def mutate_paths(self, chance: float) -> None:
        for path in self.Population:
            path.mutate(chance)

    def verbose_yTotal(self, interval_k: List) -> List:
        yTotal = [0] * 401
        for i in range(self.npaths):
            path = self.pathDictionary.get(self.pathname[i])
            Individual = self.get()
            path.e0 = Individual[i][1]
            path.s02 = Individual[i][0]
            path.sigma2 = Individual[i][2]
            path.deltar = Individual[i][3]
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
