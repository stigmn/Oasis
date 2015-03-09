__author__ = 'Joakim Boe <joakim.bo@mn.uio.no>'
__date__ = '2015-02-28'
__copyright__ = 'Copyright (C) 2015 ' + __author__
__license__  = 'GNU Lesser GPL version 3 or any later version'

from dolfin import TrialFunction, TestFunction, dx 
from DynamicModules import tophatfilter, lagrange_average, compute_Lij,\
        compute_Mij, compute_Leonard, update_mixedLESSource, dyn_u_ops
import DynamicLagrangian
import numpy as np

__all__ = ['les_setup', 'les_update']

def les_setup(u_, mesh, dt, krylov_solvers, V, assemble_matrix, CG1Function, nut_krylov_solver, 
        bcs, u_components, MixedDynamicLagrangian, DynamicSmagorinsky, **NS_namespace):
    """
    Set up for solving the mixed scale similar Germano Dynamic LES model applying
    Lagrangian Averaging. The implementation is based on the work of
    Vreman et.al. 1994, "On the Formulation of the Dynamic Mixed Subgrid-scale
    Model" and their so called DMM2 model. 
    
    Results showed that both the DMM1 and DMM2 models performed better than 
    the Germano Dynamic Model. However an inconsistency present in DMM1 is 
    fixed in DMM2, resulting in even better results.

    Cs**2 = avg((Lij-Hij)Mij)/avg(MijMij)
    
    where

    Lij = F(uiuj) - F(ui)F(uj)
    Mij = 2*delta**2*(F(|S|Sij)-alpha**2*F(|S|)F(Sij))
    Hij = F(G(F(ui)F(uj)))-F(G(F(ui)))F(G(F(uj))) - F(G(uiuj)) + F(G(ui)G(uj))
    
    and the Leonard tensor is

    Lij_L = dev(G(uiuj)-G(ui)G(uj))
    
    SGS stress modeled as

    tau_ij = Lij_L - 2*Cs**2*delta**2*|S|*Sij

    - Test filter F = 1 iteration on filter
    - Grid filter G = 0.75 iteration on filter
    """
    
    # The setup is 99% equal to DynamicLagrangian, hence use its les_setup
    dyn_dict = DynamicLagrangian.les_setup(**vars())
    
    # Set up functions for scale similarity tensor Hij
    Hij = [dyn_dict["dummy"].copy() for i in range(dyn_dict["tensdim"])]
    mixedmats = [assemble_matrix(TrialFunction(dyn_dict["CG1"]).dx(i)*TestFunction(V)*dx)
        for i in range(dyn_dict["dim"])]
    
    dummy2 = dyn_dict["dummy"].copy()
    dummy2.zero()
    
    if MixedDynamicLagrangian["model"] == "DMM2":
        from DynamicModules import compute_Hij
    elif MixedDynamicLagrangian["model"] == "DMM1":
        from DynamicModules import compute_Hij_DMM1 as compute_Hij
    
    # For check if Cs has been computed once
    Cs_bool = [False]

    dyn_dict.update(Hij=Hij, mixedmats=mixedmats, compute_Hij=compute_Hij,
            dummy2=dummy2, Cs_bool=Cs_bool)

    return dyn_dict

def les_update(u_ab, nut_, nut_form, dt, CG1, delta, tstep, u_components, V,
            DynamicSmagorinsky, Cs, u_CG1, u_filtered, Lij, Mij, Hij, lag_dt,
            JLM, JMM, dim, tensdim, G_matr, G_under, ll, mixedLESSource,
            dummy, uiuj_pairs, Sijmats, Sijcomps, Sijfcomps, delta_CG1_sq, 
            mixedmats, Sij_sol, compute_Hij, bcs_u_CG1, dummy2, vdegree,
            Cs_bool, **NS_namespace):

    # Check if Cs is to be computed, if not update nut_, mixedLESSource and break
    if tstep%DynamicSmagorinsky["Cs_comp_step"] != 0 and Cs:
        # Update nut_
        nut_()
        
        # Only update mixed source if Cs has been computed once; if not
        # the Leonard term may cause instabilities
        if Cs_bool[0] == True:
            # Update CG1-velocity
            if vdegree == 1:
                for i in range(dim):
                    u_CG1[i].vector().zero()
                    # Assign vector
                    u_CG1[i].vector().axpy(1.0, u_ab[i].vector())
            else:
                for i in range(dim):
                    ll.interpolate(u_CG1[i], u_ab[i])
            # Update mixedLESSource
            update_mixedLESSource(**vars())

        # Break les_update function
        return

    # All velocity components must be interpolated to CG1 then filtered
    dyn_u_ops(**vars())

    # Compute Lij applying dynamic modules function
    compute_Lij(u=u_CG1, uf=u_filtered, **vars())

    # Compute Mij applying dynamic modules function
    alpha = 2.0
    magS = compute_Mij(alphaval=alpha, u_nf=u_CG1, u_f=u_filtered, **vars())

    # Compute Hij
    compute_Hij(u=u_CG1, uf=u_filtered, **vars())

    # Compute Aij = Lij-Hij and add to Lij
    [Lij[i].axpy(-1.0, Hij[i]) for i in xrange(tensdim)]

    # Lagrange average (Lij-Hij) and Mij
    lagrange_average(J1=JLM, J2=JMM, Aij=Lij, Bij=Mij, **vars())

    # Update Cs = JLM/JMM and filter/smooth
    """
    Important that the term in nut_form is Cs and not Cs**2
    since Cs here is stored as JLM/JMM.
    """
    Cs.vector().set_local((JLM.vector().array()/JMM.vector().array()).clip(max=0.09))
    Cs.vector().apply("insert")
    # Filter Cs twice
    [tophatfilter(unfiltered=Cs.vector(), filtered=Cs.vector(), **vars()) for i in xrange(2)]

    # Update nut_
    nut_.vector().zero()
    nut_.vector().axpy(1.0, Cs.vector() * delta_CG1_sq * magS)
    [bc.apply(nut_.vector()) for bc in nut_.bcs]

    # Update MixedSources for rhs NS
    update_mixedLESSource(**vars())

    Cs_bool[0] = True
