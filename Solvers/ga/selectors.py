"""
GA selection operators, ported from EXAFS neoSelector.py.

Numeric option IDs are preserved so physics-module input files keep their
meaning: 0 = roulette wheel, 1 = tournament.
"""

from Solvers.core.population import Population


class SelectorBase:
    selector_type = None
    selector_operator = "Base"

    def __init__(self, n_pops: int, n_best_sample: float = 0.3, n_lucky_sample: float = 0.2, logger=None):
        self.logger = logger
        self.n_pops = n_pops
        self.nBest_Percent = n_best_sample
        self.nLucky_Percent = n_lucky_sample

        self.nBest = int(self.nBest_Percent * self.n_pops)
        self.nLucky = int(self.nLucky_Percent * self.n_pops)
        self.nCross = int(self.n_pops - self.nBest - self.nLucky)

    def select(self, pops: Population):
        pass

    def __str__(self):
        return f"Top Percentage: {100 * self.nBest_Percent}%, Lucky: {100 * self.nLucky_Percent}%"


class RouletteWheelSelector(SelectorBase):
    selector_type = 0
    selector_operator = "Roulette Wheel"

    def select(self, pops: Population):
        next_population = []
        for i in range(self.nBest):
            next_population.append(pops.population_sorted[i][0])

        for _ in range(self.nLucky):
            next_population.append(pops.generate_individual())

        pops.next_population = next_population


class TournamentSelector(SelectorBase):
    # TODO: implement tournament selection (stub carried over from neoSelector.py)
    selector_type = 1
    selector_operator = "Tournament"

    def select(self, pops: Population):
        pass


_SELECTORS = {
    0: RouletteWheelSelector,
    1: TournamentSelector,
}


def create_selector(sel_opt: int, n_pops: int, n_best_sample: float = 0.3,
                    n_lucky_sample: float = 0.2, logger=None) -> SelectorBase:
    if sel_opt not in _SELECTORS:
        raise ValueError(f"Invalid selector type: {sel_opt}")
    return _SELECTORS[sel_opt](n_pops, n_best_sample, n_lucky_sample, logger=logger)
