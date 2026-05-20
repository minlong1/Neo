import random

import attr
import numpy as np
from PhysicsModules.EXAFS.exafs_neo.pathrange import Pathrange_limits

"""
Author: Andy Lau

"""


@attr.s(slots=True)
class PathObject:
    """
    Paths Objects of a specific paths
    """

    #
    pathrange_obj: Pathrange_limits = attr.ib()
    e0: float = attr.ib()
    S02_range = attr.ib(default=list)
    E0_range = attr.ib(default=list)
    Sigma2_range = attr.ib(default=list)
    DeltaR_range = attr.ib(default=list)

    s02 = attr.ib(default=0)
    sigma2 = attr.ib(default=0)
    deltaR = attr.ib(default=0)

    def __attrs_post_init__(self):
        self.S02_range = self.pathrange_obj.getrange_S02()
        self.E0_range = self.pathrange_obj.getrange_E0()
        self.Sigma2_range = self.pathrange_obj.getrange_Sigma2()
        self.DeltaR_range = self.pathrange_obj.getrange_DeltaR()

        self.s02 = np.random.choice(self.S02_range)
        self.sigma2 = np.random.choice(self.Sigma2_range)
        self.deltaR = np.random.choice(self.DeltaR_range)

    def get(self):
        return [self.s02, self.e0, self.sigma2, self.deltaR]

    def get_var(self):
        return [self.s02, self.sigma2, self.deltaR]

    def verbose(self):
        print(self.s02, self.e0, self.sigma2, self.deltaR)

    def get_e0(self):
        return self.e0

    # -----------------
    def set_s02(self, s02):
        self.s02 = s02

    def set_e0(self, e0):
        """_summary_

        TODO: implement a checker for this checker for this
        Args:
            e0 (_type_): _description_
        """
        np.clip(e0, self.E0_range[0], self.E0_range[1])
        self.e0 = e0

    def set_sigma2(self, sigma2):
        self.sigma2 = sigma2

    def set_deltaR(self, deltaR):
        self.deltaR = deltaR

    def set(self, s02, sigma2, deltaR):
        self.set_s02(s02)
        self.set_sigma2(sigma2)
        self.set_deltaR(deltaR)

    def __mutate_s02(self, chance):
        if random.random() * 100 < chance:
            self.s02 = np.random.choice(self.pathrange_obj.getrange_S02())

    def __mutate_sigma2(self, chance):
        if random.random() * 100 < chance:
            self.sigma2 = np.random.choice(self.pathrange_obj.getrange_Sigma2())

    def __mutate_deltaR(self, chance):
        if random.random() * 100 < chance:
            self.deltaR = np.random.choice(self.pathrange_obj.getrange_DeltaR())

    def mutate(self, chance):
        """
        Mutated each of the parameters

        :param float chance: mutation chance
        :return:
        """
        self.__mutate_s02(chance)
        self.__mutate_sigma2(chance)
        self.__mutate_deltaR(chance)

    def __str__(self):
        return f"PathObject: s02: {np.round(self.s02, 2)}, e0: {np.round(self.e0, 2)}, sigma2: {np.round(self.sigma2, 3)}, deltaR: {np.round(self.deltaR, 4)}"


if __name__ == "__main__":
    e0 = 0.37
    pathrange_Obj = Pathrange_limits(0)
    pathObj = PathObject(pathrange_Obj, 0.37)
    print(pathObj)
