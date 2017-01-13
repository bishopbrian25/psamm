# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Implementation of Regulatory on/off minimization."""

import logging
from itertools import product

from six import iteritems

from psamm.lpsolver import lp
from psamm.lpsolver.lp import VariableType

logger = logging.getLogger(__name__)


class ROOMError(Exception):
    """Error indicating an error solving ROOM."""


class ConstraintGroup(object):
    def __init__(self, room, *args):
        self._room = room
        self._constrs = []
        if len(args) > 0:
            self.add(*args)

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        self.delete()

    def add(self, *args):
        self._constrs.extend(self._room._prob.add_linear_constraints(*args))

    def delete(self):
        self._room._remove_constr.extend(self._constrs)


class ROOMProblem(object):
    def __init__(self, model, solver):
        self._prob = solver.create_problem()
        self._model = model

        self._remove_constr = []

        self._v_wt = v_wt = self._prob.namespace()
        self._v = v = self._prob.namespace()
        self._y = y = self._prob.namespace()
        self._z = None

        # Define flux variables
        for reaction_id in self._model.reactions:
            lower, upper = self._model.limits[reaction_id]
            v_wt.define([reaction_id], lower=lower, upper=upper)
            v.define([reaction_id], lower=lower, upper=upper)
            y.define([reaction_id], lower=0, upper=1)

        # Define constraints
        mass_balance = self.constraints()
        massbalance_lhs = {
            spec: 0 for spec in product(model.compounds, ('wt', 'mod'))}
        for (compound, reaction_id), value in iteritems(self._model.matrix):
            massbalance_lhs[compound, 'wt'] += v_wt(reaction_id) * value
            massbalance_lhs[compound, 'mod'] += v(reaction_id) * value
        for compound, lhs in iteritems(massbalance_lhs):
            mass_balance.add(lhs == 0)



    @property
    # Returns the problem that we are working on
    def prob(self):
        return self._prob

    # Creates a constraints object
    def constraints(self, *args):
        return ConstraintGroup(self, *args)

    # Returns a generator of all the non exchange reactions in the model
    def _adjustment_reactions(self):
        for reaction_id in self._model.reactions:
            if not self._model.is_exchange(reaction_id):
                yield reaction_id
    def _reactions(self):
        for reaction_id in self._model.reactions:
            yield reaction_id

    # Solve the linear programming problem
    # Sense tells the solver to minimize or maximize the result
    def _solve(self, sense=None):
        # Remove temporary constraints
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()

        # Try to solve the problem
        try:
            return self._prob.solve(sense=sense)
        # Reset the remove constraints list before the next solve
        finally:
            self._remove_constr = []

    # Solves for the standard FBA of the problem
    # Objective is the biomass
    def _solve_fba(self, objective):
        self._prob.set_objective(self._v_wt(objective))

        # Solve and store the result
        result = self._solve(lp.ObjectiveSense.Maximize)

        # If no solution was found an error is raised and told to the user
        if not result:
            raise ROOMError('Unable to solve initial FBA: {}'.format(
                result.status))

        # Return the result object
        return result

    # Returns all the flux values for each of the reactions. This is used in
    # ROOM LP2 and QLP2 to minimize the change in the flux flow following a
    # gene deletion.
    def get_fba_flux(self, objective):
        # Get the result of the wild type FBA problem
        flux_result = self._solve_fba(objective)

        # Creat a dictionary of all the reaction ids and their fluxes.
        # {'reaction1': 1000, 'reaction2': 2000, ...}
        fba_fluxes = {}

        # Place all the flux values in a dictionary
        for key in self._model.reactions:
            fba_fluxes[key] = flux_result.get_value(self._v_wt(key))

        # Return the dictionary with all of the fluxes associated with the
        # reaction id
        # {'reaction1': 1000, 'reaction2': 2000, ...}
        return fba_fluxes

    # Solves for the FBA biomass. Uses the result returned by _solve_fba and
    # finds the value of the objective (biomass) function.
    def get_fba_biomass(self, objective):
        # Run FBA and store the result
        flux_result = self._solve_fba(objective)
        # Pull out the biomass value
        return flux_result.get_value(self._v_wt(objective))

    """
    biomass_constraint:
        True -> Add a constraint for the biomass
        False -> Don't add any constraints for the biomass
    """
    def minimize_room(self, wt_fluxes, biomass, biomass_constraint):
        # These are all of our non eachange reactions
        reactions = set(self._adjustment_reactions())

        # Create a constraints object
        constr = self.constraints()

        # Allows for some small flux changes without penalty
        d = 0.03
        e = 0.001

        v = self._v # Knockout variables
        y = self._y # Significance variables

        obj_exp = 0

        for reaction_id in self._adjustment_reactions():
            if reaction_id == biomass and biomass_constraint == True:
                continue
            # Upper and lower threshold for significance
            wt_u = wt_fluxes[reaction_id] + (d * abs(wt_fluxes[reaction_id])) + e
            wt_l = wt_fluxes[reaction_id] - (d * abs(wt_fluxes[reaction_id])) - e

            # Establish the upper and lower kinetic bounds to the reactions
            lower, upper = self._model.limits[reaction_id]
            # print(float(upper) - wt_u)
            # Add contraint that establishes the value of y
            constr.add(
                v(reaction_id) + (-1 * y(reaction_id) * (float(upper) - wt_u)) <= wt_u,
                v(reaction_id) + (-1 * y(reaction_id) * (float(lower) - wt_l)) >= wt_l)

            # Objective function for the sum of all signigicant flux changes
            # obj_exp += y(reaction_id)

        # self._prob.set_objective(obj_exp)
        self._prob.set_objective(y.sum(self._adjustment_reactions()))
        # Minimize the number of big flux changes. Biology favors a bunch
        # of small changes to flux distribution.
        result = self._solve(lp.ObjectiveSense.Minimize)

        # Need to make sure we catch any problems that occur in the solver
        if not result:
            raise ROOMError('Unable to solve ROOM: {}'.format(
                result.status))

        # Go back and manually delete all the constaints
        constr.delete()

    # Get the flux of the knockout model
    def get_flux(self, reaction):
        return self._prob.result.get_value(self._v(reaction))

    # Get the variable object for a specific reaction
    def get_flux_var(self, reaction):
        return self._v(reaction)
