def neocrossover_int2str(crossover_type):
    """
    Convert integer to string
    :param int crossover_type:
    :return str:
    """
    if crossover_type == 0:
        return "Uniform Crossover"
    elif crossover_type == 1:
        return "Single Point Crossover"
    elif crossover_type == 2:
        return "Dual Point Crossover"
    elif crossover_type == 3:
        return "Arithmetic Crossover"
    elif crossover_type == 4:
        return "Or Crossover"
    elif crossover_type == 5:
        return "Average Crossover"
    else:
        return "Unknown Crossover"


def neocrossover_str2int(crossover_type):
    """
    Convert string to integer
    :param str crossover_type:
    :return int:
    """
    if crossover_type == "Uniform Crossover":
        return 0
    elif crossover_type == "Single Point Crossover":
        return 1
    elif crossover_type == "Dual Point Crossover":
        return 2
    elif crossover_type == "Arithmetic Crossover":
        return 3
    elif crossover_type == "Or Crossover":
        return 4
    elif crossover_type == "Average Crossover":
        return 5
    else:
        return -1


def neomutator_int2str(mutator_type):
    """
    Convert integer to string
    :param int mutator_type:
    :return str:
    """
    if mutator_type == 0:
        return "Mutate Per Individual"
    elif mutator_type == 1:
        return "Mutate Per Path"
    elif mutator_type == 2:
        return "Mutate Per Trait"
    elif mutator_type == 3:
        return "Mutate Metropolis"
    elif mutator_type == 4:
        return "Mutate Bounded Per Range"
    elif mutator_type == 5:
        return "Mutate Differential Evolution"
    else:
        return "Unknown Mutator"


def neomutator_str2int(mutator_type):
    """
    Convert string to integer
    :param str mutator_type:
    :return int:
    """
    if mutator_type == "Mutate Per Individual":
        return 0
    elif mutator_type == "Mutate Per Path":
        return 1
    elif mutator_type == "Mutate Per Trait":
        return 2
    elif mutator_type == "Mutate Metropolis":
        return 3
    elif mutator_type == "Mutate Bounded Per Range":
        return 4
    elif mutator_type == "Mutate Differential Evolution":
        return 5
    else:
        return -1


def neoselector_int2str(selector_type):
    """
    Convert integer to string
    :param int selector_type:
    :return str:
    """
    if selector_type == 0:
        return "Roulette Wheel"
    elif selector_type == 1:
        return "Tournament"
    else:
        return "Unknown Selector"


def neoselector_str2int(selector_type):
    """
    Convert string to integer
    :param str selector_type:
    :return int:
    """
    if selector_type == "Roulette Wheel":
        return 0
    elif selector_type == "Tournament":
        return 1
    else:
        return -1


def neosolver_int2str(solver_type):
    """
    Convert integer to string
    :param int solver_type:
    :return str:
    """
    if solver_type == 0:
        return "Genetic Algorithm"
    elif solver_type == 1:
        return "Genetic Algorithm with Rechenberg"
    elif solver_type == 2:
        return "Differential Evolution"
    else:
        return "Unknown Solver"


def neosolver_str2int(solver_type):
    """
    Convert string to integer
    :param str solver_type:
    :return int:
    """
    if solver_type == "Genetic Algorithm":
        return 0
    elif solver_type == "Genetic Algorithm with Rechenberg":
        return 1
    elif solver_type == "Differential Evolution":
        return 2
    else:
        return -1
