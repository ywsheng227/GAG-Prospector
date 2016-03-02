"""
The MIT License (MIT)

Copyright (c) 2015 Yuewei Sheng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
commercial use

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
__author__ = 'Yuewei Sheng, PhD'

from __future__ import print_function
import math
import logging
from collections import namedtuple

ATOM_MASS = {'C':12, 'H':1.007825, 'N':14.003074, 'O':15.9949146,
             'S':31.972071, 'Na+':22.989222, 'K+':38.963158, 'Ca2+':39.961494,
             'H+':1.007276, 'Li+':7.015455} 
                                                                     
GAG = namedtuple('GAG', 'dHexA HexA HexN HexNAc Mann Ac SO3 NH4 HOHloss Na K Ca Li')
Formula = namedtuple('Formula', 'C H N O S Na K Ca Li proton')

def gag_to_formula(gag):
    # calculate the num of atoms from GAG formula
    num_C = 6 * (gag.dHexA + gag.HexA + gag.HexN) + 8 * gag.HexNAc \
            + 2 * gag.Ac + 6 * gag.Mann
    if gag.Mann != 0:
        num_H = 6 * gag.dHexA + 8 * gag.HexA + 11 *                                                                                                                                                                                                                                                                        gag.HexN + 13 * gag.HexNAc \
                + 2 * gag.Ac + 3 * gag.NH4 + 11 * gag.Mann + 1 \
                - 2 * gag.HOHloss 
    else:
        num_H = 6 * gag.dHexA + 8 * gag.HexA + 11 * gag.HexN + 13 * gag.HexNAc \
                + 2 * gag.Ac + 3 * gag.NH4 + 2 \
                - 2 * gag.HOHloss 
    num_proton = - gag.Na - gag.K - gag.Li - gag.Ca * 2
    num_N = gag.HexN + gag.HexNAc + gag.NH4 
    num_O = 5 * gag.dHexA + 6 * gag.HexA + 4 * gag.HexN + 5 * gag.HexNAc \
            + 3 * gag.SO3 + gag.Ac + 4 * gag.Mann + 1 - gag.HOHloss 
    num_S = gag.SO3
    num_Na = gag.Na
    num_K = gag.K
    num_Ca = gag.Ca
    num_Li = gag.Li
    return Formula(num_C, num_H, num_N, num_O, num_S, num_Na, num_K, num_Ca, num_Li, num_proton)

def calculate_MW(f):
    # from a chemical formula, calculate MW
    return f.C * ATOM_MASS['C'] + f.H * ATOM_MASS['H'] + f.N * ATOM_MASS['N'] + \
            f.O * ATOM_MASS['O'] + f.S * ATOM_MASS['S'] + f.Na * ATOM_MASS['Na+'] + \
            f.K * ATOM_MASS['K+'] + f.Ca * ATOM_MASS['Ca2+'] + f.Li * ATOM_MASS['Li+'] + \
            f.proton * ATOM_MASS['H+']

class GagProspector:
    def __init__(self, min_charge, max_charge, min_mz, max_mz, dHexA, HexA, HexN,
                 HexNAc, Mann, Ac=0, SO3=0, NH4=0, HOHloss=0, 
                 Na=0, K=0, Ca=0, Li=0, decon_type='normal'):
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.dHexA = dHexA
        self.HexA = HexA
        self.HexN = HexN
        self.HexNAc = HexNAc
        self.Mann = Mann
        self.Ac = Ac
        self.SO3 = SO3
        self.NH4 = NH4
        self.HOHloss = HOHloss
        self.Na = Na
        self.K = K
        self.Ca = Ca
        self.Li = Li
        self.decon_type = decon_type
                 
    def calculate_mz(self, gag):
        # calculate m/z list for a single gag
        f = gag_to_formula(gag)
        mw = calculate_MW(f)
        charge, ans = self.min_charge, []
        while charge <= self.max_charge:
            mz = mw / charge - ATOM_MASS['H+']
            if mz >= self.min_mz and mz <= self.max_mz:
                ans.append([mz, charge, gag])
            if mz < self.min_mz: break
            charge += 1
        return ans

    def calculate_mass(self, gag):
        #calculate mass for a single gag
        f = gag_to_formula(gag)
        return calculate_MW(f)

    def build_database(self):
        db, num_gag = [], 0
        for num_Ac in range(min(self.HexN, self.Ac) + 1):
            for num_SO3 in range(1, self.SO3 + 1):
                    for num_NH4 in range(self.NH4 + 1):
                        for num_Na in range(self.Na + 1):
                            for num_K in range(self.K + 1):
                                for num_Ca in range(self.Ca + 1):
                                    for num_Li in range(self.Li + 1):
                                        for num_HOHloss in range(self.HOHloss + 1): 
                                            gag = GAG(self.dHexA, self.HexA, self.HexN, self.HexNAc, 
                                                    self.Mann, num_Ac, num_SO3, num_NH4, num_HOHloss, 
                                                    num_Na, num_K, num_Ca, num_Li)
                                            num_gag += 1
                                            if self.decon_type == 'normal':
                                                db += [(self.calculate_mass(gag), gag)]
                                            else:
                                                db += self.calculate_mz(gag)
        if self.decon_type == 'normal':
            db.sort(key = lambda x: x[0])
        return db, num_gag


