"""This module provies exact stochastic simulation via the Gillespie
algorithm.

Data structures:

A simulation is a list of reactions + an initial state

A reaction is a list of products, a list of reactants and a reaction
rate constant.

Output:

A list of lists of the form [[time, species1, species2,...,speciesN]]

Must optionally incorporate forcing function, i.e. species
concentration is known at all times

Usage:


"""
