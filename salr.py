import hoomd
from hoomd import *
from hoomd import md
from os import *
import random
import numpy as np
import math as ma

n1 = 25; n2 = 25; n3 = 25; N = n1*n2*n3

dt        = 0.005
T_start   = 0.6
Temp      = 0.001
Press     = 0.5
npt_steps = 1e8

def thermalize(system, kT):
    npart = len(system.particles)
    px = 0
    py = 0
    pz = 0
    for p in system.particles:
        # thermalize velocity distribution
        mass = p.mass
        vx = random.gauss(0,np.sqrt(1.0*kT/mass))
        vy = random.gauss(0,np.sqrt(1.0*kT/mass))
        vz = random.gauss(0,np.sqrt(1.0*kT/mass))
        p.velocity = (vx,vy,vz)
        # sum the total system momentum
        px += mass*vx
        py += mass*vy
        pz += mass*vz

    # compute average momentum
    px /= npart
    py /= npart
    pz /= npart

    # subtract that average momentum from each particle
    for p in system.particles:
        mass = p.mass
        v = p.velocity
        p.velocity = (v[0] - px/mass, v[1] - py/mass, v[2] - pz/mass)

# User-defined potential
def LJ12_6shift(r, rmin, rmax, epsilon, sigma):
    """12-6 shifted Lennard-Jones potential"""
    if (rmin < r < rmax):
        V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6 - (sigma / rmax)**12 + (sigma / rmax)**6);
        F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
    else:
        V = 0.0
        F = 0.0
    return (V, F)

def Yukawashift(r, rmin, rmax, epsilon, kappa):
    """Shifted Yukawa potential"""
    if (rmin < r < rmax):
        V = epsilon * ( ma.exp( -kappa * r) / r - ma.exp( -kappa * rmax) / rmax)
        F = epsilon * ma.exp( -kappa * r) * (kappa  + 1 / r) / r
    else:
        V = 0.0
        F = 0.0
    return (V, F)

def SALR(r, rmin, rmax, epsLJ, sigma, epsY, kappa):
    """Pekalski & Santos short-ranged attraction long-ranged repulsion potential"""
    if (rmin < r < rmax):
        LJ = LJ12_6shift(r, rmin, rmax, epsLJ, sigma)
        Yukawa = Yukawashift(r, rmin, rmax, epsY, kappa)
        V = LJ[0] + Yukawa[0]
        F = LJ[1] + Yukawa[1]
    else:
        V = 0.0
        F = 0.0
    return (V, F)

# pull from a gaussian distribution of velocities

# Random initial configuration
s1 = hoomd.context.initialize("--notice-level=1 --mode=gpu");
system = hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=1.5), n=[n1,n2,n3]);

# Interaction potential
nl = hoomd.md.nlist.cell();
all = group.all();
tab_potential = md.pair.table(width=10000, nlist=nl)
tab_potential.pair_coeff.set('A', 'A', func=SALR, rmin=0.5, rmax=7.0, coeff=dict(epsLJ=1.0, sigma=1.0, epsY=0.50, kappa=0.5))
thermalize(system, T_start)

hoomd.md.integrate.mode_standard(dt=dt);
hoomd.md.integrate.npt(group=all, kT=variant.linear_interp(points = [(0,T_start),(npt_steps,Temp)]), tau=0.5, tauP=0.5, P=Press, couple='xyz');
dump.gsd(filename=f'dump_npt_p{Press:.2f}_T{Temp:.2f}.gsd',period=1e4,truncate=False,group=group.all(),phase=0, dynamic=['momentum'])

# Collect output
analyze.log(filename=f'analiza_npt_p{Press:.2f}_T{Temp:.2f}.log',quantities=['potential_energy','kinetic_energy', 'temperature', 'pressure','volume','lx','ly','lz'],period=1e4,overwrite=True)

# Run
run(npt_steps)




