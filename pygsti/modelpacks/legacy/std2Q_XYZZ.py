#***************************************************************************************************
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************
"""
Variables for working with the 2-qubit model containing the gates
I*X(pi/2), I*Y(pi/2), X(pi/2)*I, Y(pi/2)*I, and ZZ.
"""

import sys as _sys
import numpy as np

from ...circuits import circuitconstruction as _strc
from ...models import modelconstruction as _setc
from .. import stdtarget as _stdtarget
from ... import change_basis, unitary_to_process_mx

description = "I*X(pi/2), I*Y(pi/2), X(pi/2)*I, Y(pi/2)*I, and zz gates"

gates = ['Gix', 'Giy', 'Gxi', 'Gyi', 'Gzz']

fiducials16 = _strc.to_circuits(
    [(), ('Gix',), ('Giy',), ('Gix', 'Gix'),
     ('Gxi',), ('Gxi', 'Gix'), ('Gxi', 'Giy'), ('Gxi', 'Gix', 'Gix'),
     ('Gyi',), ('Gyi', 'Gix'), ('Gyi', 'Giy'), ('Gyi', 'Gix', 'Gix'),
     ('Gxi', 'Gxi'), ('Gxi', 'Gxi', 'Gix'), ('Gxi', 'Gxi', 'Giy'), ('Gxi', 'Gxi', 'Gix', 'Gix')], line_labels=('*',))

fiducials36 = _strc.to_circuits(
    [(),
     ('Gix',),
     ('Giy',),
     ('Gix', 'Gix'),
     ('Gix', 'Gix', 'Gix'),
     ('Giy', 'Giy', 'Giy'),
     ('Gxi',),
     ('Gxi', 'Gix'),
     ('Gxi', 'Giy'),
     ('Gxi', 'Gix', 'Gix'),
     ('Gxi', 'Gix', 'Gix', 'Gix'),
     ('Gxi', 'Giy', 'Giy', 'Giy'),
     ('Gyi',),
     ('Gyi', 'Gix'),
     ('Gyi', 'Giy'),
     ('Gyi', 'Gix', 'Gix'),
     ('Gyi', 'Gix', 'Gix', 'Gix'),
     ('Gyi', 'Giy', 'Giy', 'Giy'),
     ('Gxi', 'Gxi'),
     ('Gxi', 'Gxi', 'Gix'),
     ('Gxi', 'Gxi', 'Giy'),
     ('Gxi', 'Gxi', 'Gix', 'Gix'),
     ('Gxi', 'Gxi', 'Gix', 'Gix', 'Gix'),
     ('Gxi', 'Gxi', 'Giy', 'Giy', 'Giy'),
     ('Gxi', 'Gxi', 'Gxi'),
     ('Gxi', 'Gxi', 'Gxi', 'Gix'),
     ('Gxi', 'Gxi', 'Gxi', 'Giy'),
     ('Gxi', 'Gxi', 'Gxi', 'Gix', 'Gix'),
     ('Gxi', 'Gxi', 'Gxi', 'Gix', 'Gix', 'Gix'),
     ('Gxi', 'Gxi', 'Gxi', 'Giy', 'Giy', 'Giy'),
     ('Gyi', 'Gyi', 'Gyi'),
     ('Gyi', 'Gyi', 'Gyi', 'Gix'),
     ('Gyi', 'Gyi', 'Gyi', 'Giy'),
     ('Gyi', 'Gyi', 'Gyi', 'Gix', 'Gix'),
     ('Gyi', 'Gyi', 'Gyi', 'Gix', 'Gix', 'Gix'),
     ('Gyi', 'Gyi', 'Gyi', 'Giy', 'Giy', 'Giy')], line_labels=('*',))

fiducials = fiducials16
prepStrs = fiducials16

effectStrs = _strc.to_circuits(
    [(), ('Gix',), ('Giy',),
     ('Gix', 'Gix'), ('Gxi',),
     ('Gyi',), ('Gxi', 'Gxi'),
     ('Gxi', 'Gix'), ('Gxi', 'Giy'),
     ('Gyi', 'Gix'), ('Gyi', 'Giy')], line_labels=('*',))

legacy_effectStrs = _strc.to_circuits(
    [(), ('Gix',), ('Giy',), ('Gxi',), ('Gyi',),
     ('Gix', 'Gxi'), ('Gxi', 'Giy'), ('Gyi', 'Gix'),
     ('Gyi', 'Giy'), ('Gxi', 'Gxi')], line_labels=('*',))

germs = _strc.to_circuits(
    [('Gxi',),
     ('Gyi',),
     ('Gix',),
     ('Giy',),
     ('Gzz',),
     ('Gxi', 'Gyi'),
     ('Gix', 'Giy'),
     ('Giy', 'Gyi'),
     ('Gix', 'Gxi'),
     ('Gix', 'Gyi'),
     ('Giy', 'Gxi'),
     ('Giy', 'Gzz'),
     ('Gyi', 'Gzz'),
     ('Gxi', 'Gzz'),
     ('Gix', 'Gzz'),
     ('Gxi', 'Gxi', 'Gyi'),
     ('Gix', 'Gix', 'Giy'),
     ('Gxi', 'Gyi', 'Gyi'),
     ('Gix', 'Giy', 'Giy'),
     ('Giy', 'Gxi', 'Gxi'),
     ('Giy', 'Gxi', 'Gyi'),
     ('Gix', 'Gxi', 'Giy'),
     ('Gix', 'Gyi', 'Gxi'),
     ('Gix', 'Gyi', 'Giy'),
     ('Gix', 'Giy', 'Gyi'),
     ('Gix', 'Giy', 'Gxi'),
     ('Giy', 'Gyi', 'Gxi'),
     ('Giy', 'Gzz', 'Gxi'),
     ('Gix', 'Gix', 'Gzz'),
     ('Gxi', 'Gzz', 'Gzz'),
     ('Giy', 'Gyi', 'Gzz'),
     ('Gyi', 'Gzz', 'Gzz'),
     ('Gix', 'Giy', 'Gzz'),
     ('Giy', 'Gzz', 'Gzz'),
     ('Gzz', 'Gix', 'Gxi', 'Gxi'),
     ('Gyi', 'Gix', 'Gxi', 'Giy'),
     ('Gix', 'Giy', 'Gxi', 'Gyi'),
     ('Gix', 'Gix', 'Gix', 'Giy'),
     ('Gxi', 'Gyi', 'Gyi', 'Gyi'),
     ('Gyi', 'Gyi', 'Giy', 'Gyi'),
     ('Gyi', 'Gix', 'Gix', 'Gix'),
     ('Gxi', 'Gyi', 'Gix', 'Gix'),
     ('Gzz', 'Gyi', 'Gzz', 'Gxi'),
     ('Gzz', 'Gix', 'Gzz', 'Giy'),
     ('Gix', 'Gxi', 'Gzz', 'Giy'),
     ('Gix', 'Giy', 'Gxi', 'Gzz'),
     ('Giy', 'Gyi', 'Gxi', 'Gxi', 'Giy'),
     ('Gxi', 'Gxi', 'Giy', 'Gyi', 'Giy'),
     ('Giy', 'Gix', 'Gxi', 'Gix', 'Gxi'),
     ('Gyi', 'Giy', 'Gyi', 'Gix', 'Gix'),
     ('Giy', 'Gxi', 'Gix', 'Giy', 'Gyi'),
     ('Giy', 'Giy', 'Gxi', 'Gyi', 'Gxi'),
     ('Gyi', 'Gzz', 'Giy', 'Giy', 'Gxi'),
     ('Gxi', 'Gix', 'Giy', 'Gxi', 'Giy', 'Gyi'),
     ('Gxi', 'Giy', 'Gix', 'Gyi', 'Gix', 'Gix'),
     ('Gyi', 'Giy', 'Gxi', 'Gyi', 'Gxi', 'Gzz'),
     ('Gxi', 'Gxi', 'Gyi', 'Gxi', 'Gyi', 'Gyi'),
     ('Gix', 'Gix', 'Giy', 'Gix', 'Giy', 'Giy'),
     ('Gyi', 'Gxi', 'Gix', 'Giy', 'Gxi', 'Gix'),
     ('Gyi', 'Gxi', 'Gix', 'Gxi', 'Gix', 'Giy'),
     ('Gxi', 'Gix', 'Giy', 'Giy', 'Gxi', 'Gyi'),
     ('Gix', 'Giy', 'Giy', 'Gix', 'Gxi', 'Gxi'),
     ('Gyi', 'Giy', 'Gxi', 'Giy', 'Giy', 'Giy'),
     ('Gyi', 'Gyi', 'Gyi', 'Giy', 'Gyi', 'Gix'),
     ('Giy', 'Giy', 'Gxi', 'Giy', 'Gix', 'Giy'),
     ('Giy', 'Gix', 'Gyi', 'Gyi', 'Gix', 'Gxi', 'Giy'),
     ('Gyi', 'Gxi', 'Giy', 'Gxi', 'Gix', 'Gxi', 'Gyi', 'Giy'),
     ('Gix', 'Gix', 'Gyi', 'Gxi', 'Giy', 'Gxi', 'Giy', 'Gyi')
     ], line_labels=('*',))

germs_lite = _strc.to_circuits(
    [('Gxi',),
     ('Gyi',),
     ('Gix',),
     ('Giy',),
     ('Gzz',),
     ('Gxi', 'Gyi'),
     ('Gix', 'Giy'),
     ('Gxi', 'Gxi', 'Gyi'),
     ('Gix', 'Gix', 'Giy'),
     ('Gzz', 'Gix', 'Gxi', 'Gxi'),
     ('Gxi', 'Gix', 'Giy', 'Gxi', 'Giy', 'Gyi'),
     ('Gxi', 'Giy', 'Gix', 'Gyi', 'Gix', 'Gix'),
     ('Gyi', 'Giy', 'Gxi', 'Gyi', 'Gxi', 'Gzz'),
     ('Gyi', 'Gxi', 'Giy', 'Gxi', 'Gix', 'Gxi', 'Gyi', 'Giy')
     ], line_labels=('*',))


legacy_germs = _strc.to_circuits(
    [('Gxi',), ('Gyi',), ('Gix',), ('Giy',), ('Gzz',),
     ('Gxi', 'Gyi'), ('Gix', 'Giy'), ('Giy', 'Gyi'), ('Gix', 'Gyi'),
     ('Gyi', 'Gzz'), ('Giy', 'Gzz'),
     ('Gxi', 'Gzz', 'Gzz'),
     ('Giy', 'Gxi', 'Gzz'),
     ('Giy', 'Gzz', 'Gyi'),
     ('Giy', 'Gyi', 'Gzz'),
     ('Gix', 'Gxi', 'Gzz'),
     ('Giy', 'Giy', 'Gzz'),
     ('Giy', 'Gzz', 'Gxi'),
     ('Gix', 'Giy', 'Gzz'),
     ('Giy', 'Gxi', 'Gyi'),
     ('Gix', 'Giy', 'Gyi'),
     ('Gyi', 'Gyi', 'Gyi', 'Gxi'),
     ('Giy', 'Giy', 'Giy', 'Gix'),
     ('Gxi', 'Gyi', 'Gix', 'Giy'),
     ('Gzz', 'Gix', 'Gyi', 'Gyi'),
     ('Gzz', 'Gix', 'Gix', 'Gzz'),
     ('Gxi', 'Gzz', 'Gyi', 'Gyi'),
     ('Gyi', 'Gyi', 'Gyi', 'Gix'),
     ('Gix', 'Gix', 'Giy', 'Gzz', 'Gzz'),
     ('Gzz', 'Giy', 'Giy', 'Gix', 'Giy'),
     ('Gyi', 'Gzz', 'Gix', 'Giy', 'Gyi'),
     ('Giy', 'Gxi', 'Gzz', 'Gxi', 'Gzz'),
     ('Gyi', 'Gzz', 'Gxi', 'Gzz', 'Gxi'),
     ('Gyi', 'Gxi', 'Gyi', 'Gxi', 'Gxi', 'Gxi'),
     ('Gyi', 'Gxi', 'Gyi', 'Gyi', 'Gxi', 'Gxi'),
     ('Gyi', 'Gyi', 'Gyi', 'Gxi', 'Gyi', 'Gxi'),
     ('Gxi', 'Gxi', 'Gyi', 'Gxi', 'Gyi', 'Gyi'),
     ('Giy', 'Gix', 'Giy', 'Gix', 'Gix', 'Gix'),
     ('Giy', 'Gix', 'Giy', 'Giy', 'Gix', 'Gix'),
     ('Giy', 'Giy', 'Giy', 'Gix', 'Giy', 'Gix'),
     ('Gix', 'Gix', 'Giy', 'Gix', 'Giy', 'Giy'),
     ('Gzz', 'Gyi', 'Giy', 'Gxi', 'Gix', 'Gzz'),
     ('Gxi', 'Giy', 'Gxi', 'Gzz', 'Gyi', 'Gix'),
     ('Gxi', 'Giy', 'Giy', 'Giy', 'Gzz', 'Gxi'),
     ('Gzz', 'Gxi', 'Gzz', 'Gxi', 'Giy', 'Gix'),
     ('Gyi', 'Gix', 'Gyi', 'Gix', 'Gxi', 'Gxi'),
     ('Gix', 'Gzz', 'Gxi', 'Gix', 'Gxi', 'Gzz'),
     ('Gxi', 'Giy', 'Gyi', 'Gxi', 'Gzz', 'Gzz'),
     ('Gix', 'Gix', 'Giy', 'Gzz', 'Giy', 'Gzz', 'Gxi'),
     ('Giy', 'Gxi', 'Gzz', 'Gix', 'Gix', 'Giy', 'Giy'),
     ('Gxi', 'Gzz', 'Giy', 'Gyi', 'Gxi', 'Gix', 'Giy'),
     ('Gzz', 'Gzz', 'Gix', 'Gxi', 'Giy', 'Gxi', 'Gxi'),
     ('Gxi', 'Gix', 'Giy', 'Gyi', 'Gix', 'Gix', 'Gix'),
     ('Gxi', 'Gix', 'Gyi', 'Gix', 'Gyi', 'Giy', 'Gyi'),
     ('Gix', 'Gix', 'Gix', 'Gix', 'Gxi', 'Gxi', 'Gyi'),
     ('Giy', 'Gzz', 'Gxi', 'Gyi', 'Gyi', 'Gzz', 'Gix', 'Gzz'),
     ('Gxi', 'Gyi', 'Gxi', 'Giy', 'Gxi', 'Giy', 'Gix', 'Giy'),
     ('Giy', 'Giy', 'Gyi', 'Gix', 'Gzz', 'Gxi', 'Gyi', 'Gyi'),
     ('Gxi', 'Gix', 'Gzz', 'Gyi', 'Gix', 'Gzz', 'Gix', 'Giy'),
     ('Gix', 'Gxi', 'Gxi', 'Giy', 'Gxi', 'Gyi', 'Gix', 'Gzz'),
     ('Gix', 'Gix', 'Gyi', 'Gxi', 'Giy', 'Gix', 'Gzz', 'Gyi'),
     ('Gix', 'Giy', 'Gix', 'Gxi', 'Gix', 'Giy', 'Gxi', 'Gxi'),
     ('Giy', 'Gix', 'Gzz', 'Gxi', 'Gzz', 'Gxi', 'Gzz', 'Gyi'),
     ('Gxi', 'Giy', 'Gix', 'Gix', 'Gxi', 'Giy', 'Gxi', 'Gzz'),
     ('Gyi', 'Gyi', 'Gyi', 'Gyi', 'Gix', 'Giy', 'Gix', 'Gyi')
     ])


#Construct the target model
_target_model = _setc.create_explicit_model_from_expressions( 
            [('Q0','Q1')],['Gix','Giy','Gxi','Gyi'], 
            [ "X(pi/2,Q1)", "Y(pi/2,Q1)", "X(pi/2,Q0)", "Y(pi/2,Q0)" ],
            effect_labels=['00','01','10','11'], effect_expressions=["0","1","2","3"])

#Unitary in acting on the state-space { |A>, |B>, |C>, |D> } == { |00>, |01>, |10>, |11> }.
# This unitary rotates the second qubit by pi/2 in either the (+) or (-) direction based on 
# the state of the first qubit.
_myUnitary = 1./np.sqrt(2) * np.array([[1-1j,0,0,0],
                                      [0,1+1j,0,0],
                                      [0,0,1+1j,0],
                                      [0,0,0,1-1j]])

#Convert this unitary into a "superoperator", which acts on the 
# space of vectorized density matrices instead of just the state space.
# These superoperators are what GST calls "gates".
_mySuperOp_stdbasis = unitary_to_process_mx(_myUnitary)

#After the call to unitary_to_process_mx, the superoperator is a complex matrix
# in the "standard" or "matrix unit" basis given by { |A><A|, |A><B|, etc }.
# For use in GST, we want to work with a *real* matrix in either the 
# Gell-Mann or Pauli-product basis. Here we choose the Pauli-product basis,
# which is typically more intuitive when working with 2 qubits.
_mySuperOp_ppbasis = change_basis(_mySuperOp_stdbasis, "std", "pp")

#The resulting superoperator in the Pauli-product basis is exactly
# what goes into the Model object, which can be set using 
# dictionary syntax.  The line below names our two-qubit gate 'Gtq'
_target_model['Gzz'] = _mySuperOp_ppbasis

_gscache = {("full", "auto"): _target_model}


def processor_spec():
    return target_model('static').create_processor_spec(None)


def target_model(parameterization_type="full", sim_type="auto"):
    """
    Returns a copy of the target model in the given parameterization.

    Parameters
    ----------
    parameterization_type : {"TP", "CPTP", "H+S", "S", ... }
        The gate and SPAM vector parameterization type. See
        :function:`Model.set_all_parameterizations` for all allowed values.

    sim_type : {"auto", "matrix", "map", "termorder:X" }
        The simulator type to be used for model calculations (leave as
        "auto" if you're not sure what this is).

    Returns
    -------
    Model
    """
    return _stdtarget._copy_target(_sys.modules[__name__], parameterization_type,
                                   sim_type, _gscache)


#Wrong zz (bad 1Q phase factor)
legacy_gs_target = _target_model


global_fidPairs = [
    (0, 1), (0, 5), (1, 3), (2, 4), (3, 8), (5, 5), (7, 0), (9, 3),
    (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6), (15, 0),
    (15, 5)]

pergerm_fidPairsDict = {
    ('Gyi',): [
        (3, 1), (4, 1), (4, 2), (5, 0), (5, 1), (5, 7), (6, 0),
        (6, 8), (7, 2), (7, 4), (7, 9), (8, 0), (8, 7), (9, 2),
        (9, 3), (10, 9), (10, 10), (14, 7), (14, 9), (15, 10)],
    ('Gix',): [
        (0, 5), (1, 0), (1, 1), (2, 2), (2, 5), (2, 9), (3, 3),
        (3, 4), (3, 8), (4, 0), (4, 2), (4, 7), (4, 8), (4, 10),
        (5, 0), (5, 1), (5, 2), (5, 6), (5, 8), (6, 7), (6, 8),
        (6, 9), (7, 0), (7, 4), (8, 5), (8, 9), (9, 5), (10, 8),
        (10, 10), (12, 2), (12, 4), (12, 7), (13, 2), (13, 3),
        (13, 9), (14, 0), (14, 5), (14, 6), (15, 5), (15, 8),
        (15, 9)],
    ('Giy',): [
        (0, 0), (0, 7), (1, 1), (3, 5), (3, 6), (4, 2), (4, 4),
        (4, 5), (5, 3), (5, 7), (7, 1), (7, 8), (8, 5), (9, 4),
        (9, 5), (9, 9), (10, 5), (11, 5), (11, 6), (11, 8), (11, 10),
        (12, 0), (12, 3), (13, 10), (14, 0), (14, 5), (14, 6),
        (14, 7), (15, 0), (15, 6), (15, 9)],
    ('Gzz',): [
        (0, 7), (1, 0), (2, 4), (4, 3), (5, 5), (5, 10), (6, 1),
        (6, 10), (7, 8), (8, 2), (8, 5), (8, 6), (9, 0), (9, 4),
        (9, 6), (9, 9), (10, 2), (10, 4), (10, 5), (11, 5), (12, 2),
        (13, 7), (14, 1), (14, 2), (14, 3), (14, 6), (15, 8),
        (15, 10)],
    ('Gxi',): [
        (0, 7), (1, 1), (1, 7), (2, 7), (3, 3), (4, 9), (5, 4),
        (7, 2), (7, 10), (8, 2), (9, 2), (9, 8), (9, 9), (10, 1),
        (10, 10), (11, 2), (11, 5), (11, 6), (13, 2), (14, 7),
        (15, 2), (15, 3)],
    ('Giy', 'Gyi'): [
        (0, 6), (0, 8), (0, 10), (1, 0), (1, 1), (1, 3), (2, 9),
        (3, 8), (4, 4), (4, 7), (5, 7), (6, 1), (7, 0), (7, 8),
        (9, 10), (10, 5), (11, 5), (12, 5), (12, 6), (14, 0),
        (15, 0), (15, 6), (15, 8)],
    ('Gix', 'Gxi'): [
        (0, 0), (1, 5), (2, 4), (3, 3), (3, 5), (5, 2), (6, 1),
        (6, 8), (6, 10), (8, 6), (10, 2), (10, 8), (10, 10),
        (11, 8), (12, 1), (13, 1), (13, 4), (13, 6), (13, 10),
        (14, 8), (15, 3)],
    ('Giy', 'Gzz'): [
        (0, 2), (1, 0), (1, 4), (1, 9), (3, 10), (4, 3), (5, 7),
        (7, 4), (7, 7), (7, 8), (8, 7), (8, 9), (9, 2), (9, 6),
        (10, 3), (14, 10), (15, 4)],
    ('Gxi', 'Gyi'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Giy', 'Gxi'): [
        (1, 1), (2, 8), (3, 0), (3, 2), (3, 6), (4, 7), (7, 2),
        (8, 6), (9, 1), (9, 7), (9, 9), (10, 2), (10, 10), (11, 8),
        (12, 6), (13, 2), (13, 7), (14, 2), (15, 5)],
    ('Gix', 'Gzz'): [
        (2, 2), (3, 1), (4, 5), (4, 7), (4, 8), (5, 3), (5, 4),
        (8, 5), (9, 4), (9, 6), (10, 10), (12, 0), (13, 0), (14, 5),
        (15, 1), (15, 7), (15, 10)],
    ('Gyi', 'Gzz'): [
        (0, 2), (1, 0), (1, 4), (1, 9), (3, 10), (4, 3), (5, 7),
        (7, 4), (7, 7), (7, 8), (8, 7), (8, 9), (9, 2), (9, 6),
        (10, 3), (14, 10), (15, 4)],
    ('Gxi', 'Gzz'): [
        (0, 1), (0, 5), (1, 3), (2, 4), (2, 10), (3, 8), (5, 5),
        (7, 0), (9, 3), (9, 9), (9, 10), (10, 8), (12, 2), (12, 6),
        (14, 6), (15, 0), (15, 5)],
    ('Gix', 'Giy'): [
        (1, 0), (1, 10), (4, 0), (4, 4), (4, 7), (4, 8), (5, 5),
        (7, 6), (8, 9), (9, 9), (10, 2), (10, 8), (11, 10), (12, 6),
        (12, 9), (13, 9), (15, 1)],
    ('Gix', 'Gyi'): [
        (0, 5), (0, 9), (1, 6), (3, 1), (3, 2), (5, 0), (5, 4),
        (6, 0), (6, 8), (9, 7), (10, 9), (11, 1), (11, 4), (14, 4),
        (14, 9), (15, 5), (15, 7)],
    ('Gix', 'Gxi', 'Giy'): [
        (0, 6), (3, 0), (5, 0), (6, 7), (7, 1), (8, 3), (9, 9),
        (10, 4), (10, 9), (12, 9), (13, 2), (14, 5), (14, 8),
        (14, 10), (15, 6)],
    ('Giy', 'Gxi', 'Gyi'): [
        (0, 9), (1, 1), (1, 9), (2, 7), (3, 4), (4, 4), (4, 10),
        (6, 0), (6, 3), (7, 0), (9, 4), (11, 5), (12, 4), (13, 7),
        (14, 0)],
    ('Giy', 'Gyi', 'Gxi'): [
        (0, 9), (1, 1), (1, 9), (2, 7), (3, 4), (4, 4), (4, 10),
        (6, 0), (6, 3), (7, 0), (9, 4), (11, 5), (12, 4), (13, 7),
        (14, 0)],
    ('Gix', 'Gix', 'Gzz'): [
        (0, 6), (0, 7), (2, 2), (2, 3), (2, 5), (2, 9), (4, 4),
        (4, 9), (6, 1), (8, 10), (12, 3), (12, 9), (13, 3), (14, 10),
        (15, 3)],
    ('Gix', 'Giy', 'Gxi'): [
        (0, 6), (3, 0), (5, 0), (6, 7), (7, 1), (8, 3), (9, 9),
        (10, 4), (10, 9), (12, 9), (13, 2), (14, 5), (14, 8),
        (14, 10), (15, 6)],
    ('Gxi', 'Gyi', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gix', 'Giy', 'Gzz'): [
        (0, 3), (1, 0), (1, 4), (3, 10), (4, 3), (5, 7), (7, 2),
        (7, 4), (7, 7), (7, 8), (8, 1), (8, 5), (8, 7), (8, 9),
        (9, 2), (9, 6), (10, 3), (14, 10), (15, 4)],
    ('Gix', 'Giy', 'Gyi'): [
        (0, 1), (4, 2), (4, 7), (6, 7), (8, 3), (9, 5), (9, 7),
        (10, 0), (10, 4), (10, 5), (11, 2), (11, 9), (14, 6),
        (14, 8), (15, 3)],
    ('Giy', 'Gyi', 'Gzz'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (7, 8), (9, 3), (9, 9), (9, 10),
        (10, 8), (12, 2), (12, 6), (14, 6), (15, 0), (15, 4),
        (15, 5)],
    ('Giy', 'Gzz', 'Gxi'): [
        (0, 1), (0, 5), (1, 3), (2, 4), (2, 10), (3, 8), (5, 5),
        (7, 0), (9, 3), (9, 9), (9, 10), (10, 8), (12, 2), (12, 6),
        (14, 6), (15, 0), (15, 5)],
    ('Gyi', 'Gzz', 'Gzz'): [
        (3, 1), (4, 1), (4, 2), (5, 0), (5, 1), (5, 7), (6, 0),
        (6, 8), (7, 2), (7, 4), (7, 9), (8, 0), (8, 7), (9, 2),
        (9, 3), (10, 9), (10, 10), (14, 7), (14, 9), (15, 10)],
    ('Gix', 'Gix', 'Giy'): [
        (0, 0), (0, 6), (1, 0), (1, 10), (4, 0), (4, 4), (4, 7),
        (4, 8), (5, 5), (6, 7), (7, 6), (8, 9), (9, 9), (10, 2),
        (10, 8), (11, 10), (12, 6), (12, 9), (13, 1), (13, 9),
        (15, 1)],
    ('Giy', 'Gxi', 'Gxi'): [
        (1, 7), (2, 2), (4, 8), (7, 2), (7, 10), (8, 6), (9, 8),
        (9, 9), (10, 1), (11, 4), (11, 9), (12, 8), (12, 9),
        (13, 0), (13, 1), (13, 9)],
    ('Gix', 'Giy', 'Giy'): [
        (0, 4), (0, 5), (0, 7), (1, 1), (1, 6), (2, 3), (4, 10),
        (5, 4), (6, 8), (7, 4), (7, 10), (8, 8), (8, 9), (10, 5),
        (11, 5), (11, 6), (11, 9), (13, 10), (14, 1), (14, 9)],
    ('Giy', 'Gzz', 'Gzz'): [
        (0, 0), (0, 7), (1, 1), (3, 5), (3, 6), (4, 2), (4, 4),
        (4, 5), (5, 3), (5, 7), (7, 1), (7, 8), (8, 5), (9, 4),
        (9, 5), (9, 9), (10, 5), (11, 5), (11, 6), (11, 8), (11, 10),
        (12, 0), (12, 3), (13, 10), (14, 0), (14, 5), (14, 6),
        (14, 7), (15, 0), (15, 6), (15, 9)],
    ('Gix', 'Gyi', 'Gxi'): [
        (1, 10), (2, 10), (4, 8), (5, 5), (5, 6), (6, 10), (7, 0),
        (7, 5), (7, 6), (7, 8), (8, 5), (12, 5), (13, 0), (13, 2),
        (14, 1)],
    ('Gxi', 'Gxi', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gxi', 'Gzz', 'Gzz'): [
        (0, 7), (1, 1), (1, 7), (2, 7), (3, 3), (4, 9), (5, 4),
        (7, 2), (7, 10), (8, 2), (9, 2), (9, 8), (9, 9), (10, 1),
        (10, 10), (11, 2), (11, 5), (11, 6), (13, 2), (14, 7),
        (15, 2), (15, 3)],
    ('Gix', 'Gyi', 'Giy'): [
        (3, 0), (4, 4), (5, 1), (5, 8), (6, 5), (7, 3), (8, 6),
        (8, 7), (9, 5), (10, 3), (11, 4), (14, 0), (14, 6), (14, 9),
        (15, 5)],
    ('Gxi', 'Gyi', 'Gix', 'Gix'): [
        (1, 10), (2, 10), (4, 8), (5, 5), (5, 6), (6, 10), (7, 0),
        (7, 5), (7, 6), (7, 8), (8, 5), (12, 5), (13, 0), (13, 2),
        (14, 1)],
    ('Gix', 'Giy', 'Gxi', 'Gzz'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gix', 'Giy', 'Gxi', 'Gyi'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gyi', 'Gix', 'Gix', 'Gix'): [
        (0, 5), (0, 9), (1, 6), (3, 1), (3, 2), (5, 0), (5, 4),
        (6, 0), (6, 8), (9, 7), (10, 9), (11, 1), (11, 4), (14, 4),
        (14, 9), (15, 5), (15, 7)],
    ('Gxi', 'Gyi', 'Gyi', 'Gyi'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gyi', 'Gyi', 'Giy', 'Gyi'): [
        (0, 2), (1, 1), (1, 4), (2, 1), (2, 10), (3, 10), (4, 0),
        (5, 3), (5, 7), (6, 4), (6, 10), (8, 2), (8, 3), (9, 0),
        (10, 8), (11, 1), (11, 7), (13, 1), (13, 8)],
    ('Gzz', 'Gyi', 'Gzz', 'Gxi'): [
        (0, 0), (1, 0), (1, 10), (4, 0), (4, 4), (5, 10), (7, 9),
        (8, 9), (9, 0), (9, 9), (10, 3), (10, 4), (10, 8), (11, 5),
        (12, 5), (12, 9), (13, 1), (13, 9), (15, 1)],
    ('Gix', 'Gix', 'Gix', 'Giy'): [
        (1, 0), (1, 10), (4, 0), (4, 4), (4, 7), (4, 8), (5, 5),
        (7, 6), (8, 9), (9, 9), (10, 2), (10, 8), (11, 10), (12, 6),
        (12, 9), (13, 9), (15, 1)],
    ('Gyi', 'Gix', 'Gxi', 'Giy'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gzz', 'Gix', 'Gzz', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gzz', 'Gix', 'Gxi', 'Gxi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gix', 'Gxi', 'Gzz', 'Giy'): [
        (0, 1), (0, 2), (0, 5), (0, 10), (1, 3), (1, 7), (2, 1),
        (2, 7), (3, 6), (5, 1), (6, 10), (7, 0), (7, 6), (9, 0),
        (9, 2), (10, 4), (11, 7), (11, 9), (13, 7), (14, 3),
        (14, 7), (15, 5), (15, 10)],
    ('Giy', 'Giy', 'Gxi', 'Gyi', 'Gxi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Giy', 'Gyi', 'Gix', 'Gix'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Gzz', 'Giy', 'Giy', 'Gxi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gxi', 'Gxi', 'Giy', 'Gyi', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Giy', 'Gyi', 'Gxi', 'Gxi', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Giy', 'Gix', 'Gxi', 'Gix', 'Gxi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Giy', 'Gxi', 'Gix', 'Giy', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gix', 'Giy', 'Giy', 'Gix', 'Gxi', 'Gxi'): [
        (1, 1), (2, 5), (4, 3), (5, 5), (6, 3), (7, 1), (10, 2),
        (10, 5), (11, 2), (11, 5), (12, 7), (12, 10), (13, 0),
        (13, 4), (14, 5)],
    ('Gxi', 'Gix', 'Giy', 'Gxi', 'Giy', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gxi', 'Gxi', 'Gyi', 'Gxi', 'Gyi', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Giy', 'Gxi', 'Giy', 'Giy', 'Giy'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gix', 'Gix', 'Giy', 'Gix', 'Giy', 'Giy'): [
        (1, 0), (1, 10), (4, 0), (4, 4), (4, 7), (4, 8), (5, 5),
        (7, 6), (8, 9), (9, 9), (10, 2), (10, 8), (11, 10), (12, 6),
        (12, 9), (13, 9), (15, 1)],
    ('Giy', 'Giy', 'Gxi', 'Giy', 'Gix', 'Giy'): [
        (0, 4), (0, 6), (1, 1), (2, 2), (4, 1), (4, 3), (5, 1),
        (5, 3), (6, 10), (8, 2), (8, 8), (9, 4), (10, 7), (12, 1),
        (13, 2), (15, 6), (15, 9)],
    ('Gxi', 'Gix', 'Giy', 'Giy', 'Gxi', 'Gyi'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Gyi', 'Gyi', 'Giy', 'Gyi', 'Gix'): [
        (0, 3), (1, 0), (1, 4), (3, 10), (4, 3), (5, 7), (7, 2),
        (7, 4), (7, 7), (7, 8), (8, 1), (8, 5), (8, 7), (8, 9),
        (9, 2), (9, 6), (10, 3), (14, 10), (15, 4)],
    ('Gxi', 'Giy', 'Gix', 'Gyi', 'Gix', 'Gix'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gyi', 'Giy', 'Gxi', 'Gyi', 'Gxi', 'Gzz'): [
        (0, 1), (0, 2), (0, 5), (1, 3), (1, 9), (2, 4), (2, 10),
        (3, 8), (5, 5), (7, 0), (9, 3), (9, 9), (9, 10), (10, 8),
        (12, 2), (12, 6), (14, 6), (15, 0), (15, 5)],
    ('Gyi', 'Gxi', 'Gix', 'Gxi', 'Gix', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Gxi', 'Gix', 'Giy', 'Gxi', 'Gix'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Giy', 'Gix', 'Gyi', 'Gyi', 'Gix', 'Gxi', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gyi', 'Gxi', 'Giy', 'Gxi', 'Gix', 'Gxi', 'Gyi', 'Giy'): [
        (0, 1), (0, 5), (1, 3), (3, 8), (5, 5), (7, 0), (9, 3),
        (9, 9), (9, 10), (10, 8), (12, 2), (12, 6), (14, 6),
        (15, 0), (15, 5)],
    ('Gix', 'Gix', 'Gyi', 'Gxi', 'Giy', 'Gxi', 'Giy', 'Gyi'): [
        (1, 1), (2, 5), (4, 3), (5, 5), (6, 3), (7, 1), (10, 2),
        (10, 5), (11, 2), (11, 5), (12, 7), (12, 10), (13, 0),
        (13, 4), (14, 5)],
}
