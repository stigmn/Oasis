from ..NSfracStep import *


NS_parameters.update(
    nu = 0.0035,
    T  = 1.0,
    dt = 0.01,
    plot_interval = 20,
    print_intermediate_info = 10,
    mesh_path = "/home/stigmn/first-test/carotid2/carotid_mesh.xml",
    use_krylov_solvers = True)

set_log_active(False)

def mesh(mesh_path, **NS_namespace):
    return Mesh(mesh_path)

inlet = Expression("0.1")
inlety = Constant(0.0)
inletz = Constant(0.0)
def create_bcs(V, **NS_namespace):
    bc01=DirichletBC(V,inletz,1)
    bc02=DirichletBC(V,inlety,1)
    bc03=DirichletBC(V,inlet,1)
    bc04=DirichletBC(V,inletz,0)
    return dict(u0 = [bc01,bc04],
                u1 = [ bc02,bc04],
                u2 = [bc03, bc04],
                p  = [])

def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:    
        [bc.apply(x_2[ui]) for bc in bcs[ui]]

def pre_solve_hook(mesh, velocity_degree, **NS_namespace):
    Vv = VectorFunctionSpace(mesh, 'CG', velocity_degree)
    return dict(uv=Function(Vv))

def temporal_hook(q_, tstep, u_, uv, p_, plot_interval, **NS_namespace):
    if tstep % plot_interval == 0:
        assign(uv.sub(0), u_[0])
        assign(uv.sub(1), u_[1])
        assign(uv.sub(2), u_[2])
        file = File("./results/Bifurcation/u%d.pvd"%(tstep) )
        file << uv
