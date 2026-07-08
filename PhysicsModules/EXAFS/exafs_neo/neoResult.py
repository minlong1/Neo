import numpy as np
from attrs import define, field
import pickle

from matplotlib import pyplot as plt

from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger
from PhysicsModules.EXAFS.exafs_neo.individual import Individual
from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations


@define
class NeoResult:
    """Per-generation history (best fitness, reduced chi-square) plus the
    overall best individual; `collect()` is called once per generation from
    `ExafsNeo.run()`."""

    exafs_pars: NeoPars = None
    best_individual: Individual = None
    historyBest: list[float] = field(factory=list)
    historyBestChiR: list[float] = field(factory=list)
    genBest: list[float] = field(factory=list)
    crossover_scorelist: list[float] = field(factory=list)
    mutation_scorelist: list[float] = field(factory=list)
    logger: NeoLogger = None

    def __str__(self) -> str:
        if self.best_individual is None:
            return "Best Individual: None"
        else:
            return f"Best Individual: {self.best_individual}"

    def save(self, filename: str) -> None:
        """
        Save the Neo Result file
        """

        try:
            with open(filename, "wb") as file:
                pickle.dump(self, file)
        except Exception as e:
            self.logger.print(f"Error saving file: {e}")

    @classmethod
    def load(cls, filename: str) -> None:
        """
        Load the result file back into it
        """
        try:
            with open(filename, "rb") as file:
                obj = pickle.load(file)
            print(f"Object loaded from {filename}")
            return obj
        except Exception as e:
            print(f"Error loading object: {e}")
            return None

    def initialize(self, exafs_pars) -> None:
        self.exafs_pars = exafs_pars

    def collect(self, neo_population: NeoPopulations, exafs_pars: NeoPars) -> None:
        """
        Collecting all result into a single area for post-processing.
        """
        self.exafs_pars = exafs_pars
        self.best_individual = neo_population.population_sorted[0][0]
        self.historyBest.append(exafs_pars.bestFitPars.globBestVal)
        self.genBest.append(exafs_pars.bestFitPars.currBestVal)
        global_r = exafs_pars.bestFitPars.globBestVal / (
            len(exafs_pars.exafsPars.intervalK) - 3 * exafs_pars.exafsPars.npath + 1
        )
        self.historyBestChiR.append(global_r)

    def plot_fitness(self) -> None:
        """
        Plot the fitness function
        :return:
        """
        x = np.arange(0, len(self.historyBest), 1)

        fig, ax1 = plt.subplots()
        ax1.plot(x, self.historyBest, "Fitness")
        ax1.set_xlabel("Generation")
        ax1.set_ylabel("Fitness", color="b")
        # ax1.tick_params('y', colors='b')
        ax1.set_xticks(x)

        ax2 = ax1.twinx()
        ax2.plot(x, self.historyBestChiR, "Fitness Reduced Chi2")
        ax2.set_ylabel("Fitness Reduced Chi2", color="r")

        # locs, labels = plt.xticks()
        # plt.plot(x, self.historyBest, label='Fitness ')
        # # plt.plot(x,self.historyBestChiR,label='Fitness Reduced Chi2')
        # plt.xlabel('Generation')
        # plt.ylabel('Fitness')
        # plt.xticks(x)


if __name__ == "__main__":
    x = "test"
    neo_pars = NeoPars()
    neo_result = NeoResult()

    neo_result.save("Test.pkl")

    neo_result = NeoResult.load("Test.pkl")
    print(neo_result)
