"""
Generic population container, ported from the EXAFS NeoPopulations but
scoring through an OptimizationProblem instead of a hard-wired fitness.
"""

import operator
from typing import List

from Solvers.core.individual import Individual
from Solvers.core.problem import OptimizationProblem
from Solvers.core.run_state import RunState


class Population:
    def __init__(self, problem: OptimizationProblem, state: RunState, n_pops: int = 100):
        self.problem = problem
        self.state = state
        self.n_pops = n_pops

        self.population: List[Individual] = []
        self.population_sorted: List = []  # list of (individual, score), best first
        self.next_population: List[Individual] = []

    def generate_individual(self) -> Individual:
        return Individual(self.problem.space, self.problem.sample_genes())

    def initialize_populations(self) -> None:
        self.population = [self.generate_individual() for _ in range(self.n_pops)]
        self.eval_population()

    def eval_population(self, replace: bool = True, sorting: bool = True) -> List[float]:
        scores = []
        performance = {}
        for individual in self.population:
            score = self.problem.fitness(individual.genes)
            scores.append(score)
            performance[individual] = score

        if sorting:
            self.population_sorted = sorted(
                performance.items(), key=operator.itemgetter(1)
            )
        if replace and self.population_sorted:
            best_individual, best_value = self.population_sorted[0]
            self.state.update_best(best_individual, best_value)
        return scores

    def best(self) -> Individual:
        if not self.population_sorted:
            raise ValueError("Population has not been evaluated yet")
        return self.population_sorted[0][0]

    def __len__(self) -> int:
        return len(self.population)

    def __getitem__(self, item):
        return self.population_sorted[item]
