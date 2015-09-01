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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

from ..command import Command, SolverCommandMixin, CommandError
from .. import fluxanalysis

from six.moves import range


class RobustnessCommand(SolverCommandMixin, Command):
    """Run robustness analysis on the model.

    Given a reaction to maximize and a reaction to vary,
    the robustness analysis will run FBA while fixing the
    reaction to vary at each iteration. The reaction will
    be fixed at the specified number of steps between the
    minimum and maximum flux value specified in the model.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--steps', metavar='N', type=int, default=10,
            help='Number of flux value steps for varying reaction')
        parser.add_argument(
            '--minimum', metavar='V', type=float,
            help='Minumum flux value of varying reacton')
        parser.add_argument(
            '--maximum', metavar='V', type=float,
            help='Maximum flux value of varying reacton')
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FBA',
            action='store_true')
        parser.add_argument(
            '--reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument('varying', help='Reaction to vary')
        super(RobustnessCommand, cls).init_parser(parser)

    def run(self):
        """Run flux analysis command"""

        def run_fba_fmin(model, reaction):
            fba = fluxanalysis.FluxBalanceProblem(model, solver=solver)
            fba.maximize(reaction)
            optimum = fba.get_flux(reaction)
            return fluxanalysis.flux_minimization(
                model, {reaction: optimum}, solver=solver)

        def run_tfba(model, reaction):
            return fluxanalysis.flux_balance(
                model, reaction, tfba=True, solver=solver)

        if self._args.no_tfba:
            run_fba = run_fba_fmin
            solver = self._get_solver()
        else:
            run_fba = run_tfba
            solver = self._get_solver(integer=True)

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        varying_reaction = self._args.varying
        if not self._mm.has_reaction(varying_reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                varying_reaction))

        steps = self._args.steps
        if steps <= 0:
            raise CommandError('Invalid number of steps: {}\n'.format(steps))

        # Run FBA on model at different fixed flux values
        flux_min = self._mm.limits[varying_reaction].lower
        flux_max = self._mm.limits[varying_reaction].upper
        if self._args.minimum is not None:
            flux_min = self._args.minimum
        if self._args.maximum is not None:
            flux_max = self._args.maximum

        if flux_min > flux_max:
            raise CommandError('Invalid flux range: {}, {}\n'.format(
                flux_min, flux_max))

        for i in range(steps):
            fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
            test_model = self._mm.copy()
            test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

            try:
                for other_reaction, flux in run_fba(test_model, reaction):
                    print('{}\t{}\t{}'.format(
                        other_reaction, fixed_flux, flux))
            except fluxanalysis.FluxBalanceError:
                pass
