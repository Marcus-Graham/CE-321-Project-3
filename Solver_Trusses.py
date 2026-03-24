#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads
def ComputeMemberForces(bars):
    # COMPLETE THIS FUNCTION
    # Compute member forces for all bars using equation 14-23
    for bar in bars:
        A=bar.A
        E=bar.E
        L=bar.Length()
        lambdax, lambday=bar.LambdaTerms()
        u1x=bar.init_node.xdisp
        u1y=bar.init_node.ydisp
        u2x=bar.end_node.xdisp
        u2y=bar.end_node.ydisp
        u=np.array([u1x,u1y,u2x,u2y])
        T=np.array([-lambdax,-lambday,lambdax,lambday])
        axial_force=(A*E/L)*np.dot(T,u)
        bar.axial_load=axial_force
    
# compute the normal stresses
def ComputeNormalStresses(bars):
    # COMPLETE THIS FUNCTION
    # Compute normal stress for all bars
    for bar in bars:
        if np.isnan(bar.axial_load):
            raise Exception("Axial load needs to be computed before stress")
        bar.normal_stress=bar.axial_load/bar.A

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    # COMPLETE THIS FUNCTION
    # Compute critical buckling load for all bars
    for bar in bars:
        E=bar.E
        L=bar.Length()
        I=bar.Iu if hasattr(bar, 'Iu') and bar.Iu != 0 else bar.It
        K=1.0
        if I==0:
            bar.buckling_load=np.nan
            continue
        Pcr=(np.pi**2*E*I)/((K*L)**2)
        bar.buckling_load=Pcr
