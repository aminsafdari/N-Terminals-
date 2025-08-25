import sys
from hoomd import md
import pandas as pd
import gsd.hoomd
import hoomd
import itertools
import numpy
import matplotlib.pyplot as plt
import math
import matplotlib.cm
from matplotlib.colors import Normalize
import numpy as np
from random import seed
from random import random
from hoomd import _hoomd
from hoomd.md import _md
from hoomd.md import force
import time
# Start timing the script
t0 = time.time()

# Set simulation parameters
param = {
    'r_ver': 0.1,       # Vertex radius
    'r_wall': 0.2,      # Wall radius
    'r_terminal': 0.05, # Terminal bead radius
    'r_gen': 0.05       # Genome bead radius
}



# ------------------------
# INPUT PARAMETERS
# ------------------------

# Read input file name and scale factor from command line
input_file_genome = 'length' + str(sys.argv[1]) + '_0_05.gsd'
vertex_file_name = '3_vertexs_1_4.csv'
scale_2 = float(sys.argv[2])
end_terminal_length = 5  # Length of the end terminal beads

print('input_file_genome',input_file_genome)
v = 1
total_time = 50000000
# total_time = 1
# Construct output file name using sanitized parts of the input

e1 = input_file_genome.replace('.','')
e2 = vertex_file_namge.replace('.','')
print('scale_2',scale_2)
output_file_genome =  f"genome_{e1[0:len(e1)-3]}_inside_r_capsid_{e2[0:len(e2)-3]}_NT_Length_{end_terminal_length}_epsilon_is_{scale_2*0.1:.5f}"
output_file_genome = output_file_genome.replace('.','_')
output_file_genome += '.gsd'


# ------------------------
# FUNCTION DEFINITIONS
# ------------------------

# Function to count how many polymer beads are inside the capsid
def count_inside(snap):
    size_p = len(index_p)
    size_v = len(index_c)
    max_r = 0
    dis_1 = 0
    centeral_mass = [0,0,0]
    count_p = 0
    
    for i in range(size_v):
        centeral_mass[0] += snap.particles.position[index_c[i]][0]
        centeral_mass[1] += snap.particles.position[index_c[i]][1]
        centeral_mass[2] += snap.particles.position[index_c[i]][2]
        
    centeral_mass[0] = centeral_mass[0] / size_v
    centeral_mass[1] = centeral_mass[1] / size_v
    centeral_mass[2] = centeral_mass[2] / size_v

    for i in range(size_v):
        dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_c[i]][0]
                 ,snap.particles.position[index_c[i]][1],snap.particles.position[index_c[i]][2])
        if (dis_1 > max_r):
            max_r = dis_1
            
    for j in range(size_p):
        dis_1 = distance(centeral_mass[0],centeral_mass[1],centeral_mass[2],snap.particles.position[index_p[j]][0]
             ,snap.particles.position[index_p[j]][1],snap.particles.position[index_p[j]][2])
        if (dis_1 < max_r ):
            count_p +=1 
    CGer = '\033[32m'
    CEND = '\033[0m'
    print(CGer+'number of polymer inside the capsid',count_p,CEND)

def density_fun(snap):
    size_p = len(index_p)
    size_v = len(index_c)
    delta_r = 0.04  # Bin width for radial histogram
    max_r = 0
    dis_1 = 0
    centeral_mass = [0, 0, 0]

    # Calculate center of mass of capsid
    for i in range(size_v):
        centeral_mass[0] += snap.particles.position[index_c[i]][0]
        centeral_mass[1] += snap.particles.position[index_c[i]][1]
        centeral_mass[2] += snap.particles.position[index_c[i]][2]
    centeral_mass = [x / size_v for x in centeral_mass]

    # Determine maximum radius of capsid from center
    for i in range(size_v):
        dis_1 = distance(centeral_mass[0], centeral_mass[1], centeral_mass[2],
                        snap.particles.position[index_c[i]][0],
                        snap.particles.position[index_c[i]][1],
                        snap.particles.position[index_c[i]][2])
        if dis_1 > max_r:
            max_r = dis_1

    step = math.floor(max_r / delta_r)  # Number of bins
    density = np.zeros(step)
    rad = np.zeros(step)

    # Loop over bins and count polymer beads in each radial shell
    for i in range(step):
        for j in range(size_p):
            dis_1 = distance(centeral_mass[0], centeral_mass[1], centeral_mass[2],
                            snap.particles.position[index_p[j]][0],
                            snap.particles.position[index_p[j]][1],
                            snap.particles.position[index_p[j]][2])
            if (i * delta_r <= dis_1 < (i + 1) * delta_r):
                density[i] += 1
        rad[i] = (i + 0.5) * delta_r  # Midpoint of each bin

    return rad, density

# Utility function to calculate Euclidean distance between two 3D points
def distance(x1,y1,z1,x2,y2,z2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

def move_to_origin(snap):
    """
    Translate all particles in the snapshot so that the system's center of mass (COM)
    is moved to the origin (0, 0, 0).
    
    Parameters:
        snap: A HOOMD snapshot object with particle positions.
    """
    # Get the number of particles
    size = len(snap.particles.position)
    
    # Compute the center of mass
    center_of_mass = [0.0, 0.0, 0.0]
    for i in range(size):
        for j in range(3):  # x, y, z components
            center_of_mass[j] += snap.particles.position[i][j]
    
    center_of_mass = [coord / size for coord in center_of_mass]

    # Shift all particles so the center of mass is at the origin
    for i in range(size):
        for j in range(3):  # x, y, z components
            snap.particles.position[i][j] -= center_of_mass[j]

def create_shell():
    """
    Constructs a shell structure based on triangular mesh data.
    Returns vertex, line, triangle data, generated wall vertices,
    and the center of mass of the capsid shell.
    """

    # Load geometry data
    vertex = pd.read_csv(vertex_file_namge)
    line = pd.read_csv('3_lines.csv')
    triangle = pd.read_csv('3_triangle.csv')

    # Initialize DataFrame to store wall vertices
    wall_vertex = pd.DataFrame(columns=['x', 'y', 'z', 't'])

    # Compute center of mass of the vertices
    centeral_mass = [
        vertex['x'].mean(),
        vertex['y'].mean(),
        vertex['z'].mean()
    ]

    # Diameter used for spacing wall mesh points
    diameter_of_m = 1 * param['r_wall']

    # Loop over each triangle to interpolate points
    for i in range(len(triangle) - 1):
        # Estimate number of divisions between two triangle vertices
        v1 = triangle.tv1[i] - 1
        v2 = triangle.tv2[i] - 1
        v3 = triangle.tv3[i] - 1

        n_div_1 = int(distance(vertex.x[v1], vertex.y[v1], vertex.z[v1],
                               vertex.x[v3], vertex.y[v3], vertex.z[v3]) / diameter_of_m)

        for j in range(n_div_1):
            # Interpolate along v1 to v3
            alpha = j / max(n_div_1, 1)
            p1x = vertex.x[v1] + alpha * (vertex.x[v3] - vertex.x[v1])
            p1y = vertex.y[v1] + alpha * (vertex.y[v3] - vertex.y[v1])
            p1z = vertex.z[v1] + alpha * (vertex.z[v3] - vertex.z[v1])

            # Interpolate along v2 to v3
            p2x = vertex.x[v2] + alpha * (vertex.x[v3] - vertex.x[v2])
            p2y = vertex.y[v2] + alpha * (vertex.y[v3] - vertex.y[v2])
            p2z = vertex.z[v2] + alpha * (vertex.z[v3] - vertex.z[v2])

            # Number of divisions between p1 and p2
            n_div_2 = int(distance(p1x, p1y, p1z, p2x, p2y, p2z) / diameter_of_m)

            for k in range(n_div_2):
                beta = k / max(n_div_2, 1)
                wall_x = p1x + beta * (p2x - p1x)
                wall_y = p1y + beta * (p2y - p1y)
                wall_z = p1z + beta * (p2z - p1z)

                # Append the interpolated wall vertex
                wall_vertex.loc[len(wall_vertex)] = [wall_x, wall_y, wall_z, i]

    # Return full shell data
    return shell(vertex, line, triangle, wall_vertex, centeral_mass)

class Shell:
    """
    A class to represent the capsid shell structure, including geometry,
    particle types, bonds, and additional metadata used in coarse-grained modeling.
    """

    def __init__(self, vertex, line, triangle, wall_vertex, centeral_mass):
        """
        Initializes a Shell object with geometric and structural data.

        Parameters:
        - vertex: DataFrame of vertex positions (x, y, z)
        - line: DataFrame of bonds/edges (connectivity)
        - triangle: DataFrame of triangle connectivity
        - wall_vertex: DataFrame of wall vertices interpolated across the shell
        - centeral_mass: 3-element list of the shell's center of mass
        """

        # Shell geometry and physical parameters (some hard-coded)
        self.rp = param['r_ver']              # Radius of each vertex bead
        self.R = 1.1                          # (not from param) Total shell radius
        self.l0 = 1.0                         # Equilibrium bond length
        self.theta = 0.644                    # Preferred bond angle
        self.theta_upper = 40.0               # Upper angular limit (degrees)
        self.theta_lower = 15.0               # Lower angular limit (degrees)
        self.gamma = 120                      # Dihedral angle preference
        self.sigma = 0.0001                   # Interaction strength or LJ sigma
        self.seed = [1]                       # RNG seed list
        self.typeid = [1]                     # Type ID for shell particles
        self.bodyid = [-1, -2, -3]            # Body ID for rigid group handling
        self.bondid = [1, 2]                  # Bond types used
        self.dihedralid = [0]                 # Dihedral type IDs

        # Assign geometric and topological structures
        self.vertex = vertex
        self.line = line
        self.triangle = triangle
        self.wall_vertex = wall_vertex
        self.centeral_mass = centeral_mass

        # Shell topology and particle tracking
        self.n = len(self.vertex)                             # Number of main shell vertices
        self.n_wall = len(self.wall_vertex)                   # Number of wall vertices
        self.particles = self.vertex[['x', 'y', 'z']].values  # Nx3 array of vertex positions
        self.pids = [1] * self.n                              # Particle type IDs for shell
        self.particles_wall = self.wall_vertex[['x', 'y', 'z']].values
        self.pids_wall = [2] * self.n_wall                    # Type ID 2 for wall particles
        self.pids_spik = [3] * (self.n * end_terminal_length) # N-terminals particle IDs
        self.bonds = np.int_(self.line[['lv1', 'lv2']])       # Bond list for shell edges
        self.bonds_spik = self.n * end_terminal_length        # Number of N-terminals bonds
        self.dihedrals = []                                   # Empty dihedral list (to be filled later)

    def update_shell_info(self):
        """
        Updates shell properties such as number of particles, bond list,
        and calls methods to update PIDs and dihedrals.
        """
        self.n = len(self.vertex)
        self.update_pids()
        self.particles = np.array(self.vertex[['x', 'y', 'z']])
        self.bonds = np.int_(self.line[['v0', 'v1']])
        self.update_dihedrals()

def initial2():
    """
    Initializes a HOOMD snapshot from a specific frame (frame 249) of a GSD file.
    The snapshot includes particle and bond information, with added metadata for types.
    
    Returns:
        snap (hoomd.Snapshot): Configured snapshot ready for further simulation setup.
    """

    # Open the GSD trajectory file in read-only binary mode
    f = gsd.hoomd.open(name=input_file_genome, mode='rb')

    # Create a blank HOOMD snapshot object
    snap = hoomd.Snapshot()

    # Set number of particles from frame 249 (last frame)
    snap.particles.N = f[249].particles.N

    # Copy particle information (position, typeid, body, diameter) from frame 249
    for i in range(snap.particles.N):
        snap.particles.position[i] = f[249].particles.position[i]
        snap.particles.typeid[i] = f[249].particles.typeid[i]
        snap.particles.body[i] = f[249].particles.body[i]
        snap.particles.diameter[i] = f[249].particles.diameter[i]

    # Set number of bonds and initialize bond types
    snap.bonds.N = f[249].bonds.N
    bond_type = []

    # Copy bond group and typeid from frame 249, assign temporary name "G" for all bond types
    for i in range(snap.bonds.N):
        bond_type.append("G")
        snap.bonds.group[i] = f[249].bonds.group[i]
        snap.bonds.typeid[i] = f[249].bonds.typeid[i]

    # Set the bond type names
    snap.bonds.types = bond_type

    # Move system so center of mass is at the origin
    move_to_origin(snap)

    # Define custom dihedral, particle, bond, and angle types
    snap.dihedrals.types = ['capsomer']
    snap.particles.types = ['P', 'X', 'W', 'S']
    snap.bonds.types = ['G', 'shell', 'spik']
    snap.angles.types = ['ds']

    return snap

 def assign_snap_particles(snap, obj, cover):
    """
    Assigns particle and bond data from a given object into a HOOMD snapshot.
    """
    
    # Determine starting particle index
    if cover:
        # Find first index where the typeid matches any in obj.typeid
        for i, pid in enumerate(snap.particles.typeid):
            if pid in obj.typeid:
                pstart = i
                break
    else:
        # If not covering existing particles, append to end
        pstart = snap.particles.N

    # Compute end index of main + terminal particles
    pend = pstart + obj.n + obj.n * end_terminal_length

    # Update the total number of particles to include main, terminal, and wall particles
    snap.particles.N = pend + obj.n_wall

    # Assign main particles
    for i, (cor, pid) in enumerate(zip(obj.particles, obj.pids), start=pstart):
        snap.particles.position[i] = cor
        snap.particles.typeid[i] = pid
        snap.particles.body[i] = obj.bodyid[0]
        snap.particles.diameter[i] = 2 * obj.rp

    # Initialize bond starting index and update total number of bonds
    bstart_bond = snap.bonds.N
    snap.bonds.N = bstart_bond + end_terminal_length * obj.n

    count = 0  # Track number of spike regions

    # Assign spike/terminal particles and bonds
    for i, (cor, pid) in enumerate(zip(obj.particles, obj.pids), start=pstart + obj.n):
        # Assign diameters for spike particles
        for j in range(end_terminal_length):
            snap.particles.diameter[pstart + obj.n + end_terminal_length * count + j] = 2 * param['r_terminal']

        # Compute direction vector from center of mass to current particle
        dx = -cor[0] + obj.centeral_mass[0]
        dy = -cor[1] + obj.centeral_mass[1]
        dz = -cor[2] + obj.centeral_mass[2]
        norm = math.sqrt(dx**2 + dy**2 + dz**2)

        # Assign spike particle positions along vector from core to center of mass
        for j in range(end_terminal_length):
            scale = 2 * param['r_terminal'] * j / norm + (param['r_ver'] + param['r_terminal']) / norm
            snap.particles.position[pstart + obj.n + end_terminal_length * count + j] = [
                cor[0] + dx * scale,
                cor[1] + dy * scale,
                cor[2] + dz * scale
            ]

        # Assign spike particle types and body ID
        for j in range(end_terminal_length):
            idx = pstart + obj.n + end_terminal_length * count + j
            snap.particles.typeid[idx] = obj.pids_spik[end_terminal_length * (i - pstart - obj.n) + j]
            snap.particles.body[idx] = obj.bodyid[2]

        # Assign bonds between core and spike particles and among spike particles
        if count == 0:
            snap.bonds.group[bstart_bond + end_terminal_length * count] = [i - obj.n, i + end_terminal_length * count]
        else:
            snap.bonds.group[bstart_bond + end_terminal_length * count] = [i - obj.n, i + (end_terminal_length - 1) * count]
        snap.bonds.typeid[bstart_bond + end_terminal_length * count] = obj.bondid[1]

        for j in range(end_terminal_length - 1):
            idx = bstart_bond + end_terminal_length * count + j + 1
            snap.bonds.group[idx] = [
                i + (end_terminal_length - 1) * count + j,
                i + (end_terminal_length - 1) * count + j + 1
            ]
            snap.bonds.typeid[idx] = obj.bondid[1]

        count += 1

    # Assign wall particles
    for i, (cor_wall, pid_wall) in enumerate(zip(obj.particles_wall, obj.pids_wall), start=pstart + (end_terminal_length + 1) * obj.n):
        snap.particles.diameter[i] = 0.4 * param['r_wall']
        snap.particles.position[i] = cor_wall
        snap.particles.typeid[i] = pid_wall
        snap.particles.body[i] = obj.bodyid[2]

    return pstart, pend
      
def assign_snap_bonds(snap, obj, pstart, cover):
    """
    Assigns bond information from an object to a HOOMD snapshot.
    """
    if cover:
        # Find the first bond index matching the desired typeid
        bstart = int(np.nonzero(snap.bonds.typeid == obj.bondid)[0][0])
    else:
        # Start appending at the end
        bstart = snap.bonds.N

    # Resize bond array to hold new bonds
    snap.bonds.N = bstart + len(obj.bonds)

    # Assign each bond
    for i, bond in enumerate(obj.bonds, start=bstart):
        snap.bonds.group[i] = pstart - 1 + bond  # Adjust bond indices relative to snapshot
        snap.bonds.typeid[i] = obj.bondid[0]     # Assign bond type


def assign_snap_dihedrals(snap, obj, pstart, cover):
    """
    Assigns dihedral information from an object to a HOOMD snapshot.
    """
    if cover:
        dstart = int(np.nonzero(snap.dihedrals.typeid == obj.dihedralid)[0][0])
    else:
        dstart = snap.dihedrals.N

    # Resize dihedral array to hold new dihedrals
    snap.dihedrals.N = dstart + len(obj.dihedrals)

    # Assign each dihedral
    for i, dihedral in enumerate(obj.dihedrals, start=dstart):
        snap.dihedrals.group[i] = pstart - 1 + dihedral  # Adjust dihedral indices
        snap.dihedrals.typeid[i] = obj.dihedralid


def assign_obj_to_snap(snap, obj, cover=True, particle=True, bond=True, dihedral=False):
    """
    Assigns an object's particle, bond, and dihedral data to a snapshot.
    """
    if particle:
        pstart, pend = assign_snap_particles(snap, obj, cover)
    if bond:
        assign_snap_bonds(snap, obj, pstart, cover)
    if dihedral:
        assign_snap_dihedrals(snap, obj, pstart, cover)
    return pstart, pend


def assign_shell(snap):
    """
    Creates and assigns a shell structure to the snapshot.
    """
    shell = create_shell()
    shstart, shend = assign_obj_to_snap(snap, shell, cover=False, dihedral=True)
    return shell


def dihedral_harmonic(theta, kappa, theta0):
    """
    Computes harmonic dihedral potential energy and force.
    """
    V = 0.5 * kappa * (1 - np.cos(theta - theta0))
    F = -0.5 * kappa * np.sin(theta - theta0)
    return V, F

def main():
    import hoomd
    import hoomd.md
    import numpy as np
    import time

    global param, output_file_genome, total_time
    t0 = time.time()

    number_of_bid = 0

    # Initialize GPU device and simulation
    hoomd.device.GPU()
    sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=1)

    # Initialize snapshot from GSD
    snap_1 = initial2()
    number_of_bid = snap_1.particles.N

    # Assign shell particles and topology
    object_of_capsid = assign_shell(snap_1)

    # Set simulation box
    BOX = hoomd.Box(Lx=100, Ly=100, Lz=100)
    snap_1.configuration.box = BOX

    print('object_of_capsid.n_wall:', object_of_capsid.n_wall)
    print("Total number of bonds:", snap_1.bonds.N)

    # Create simulation state from snapshot
    sim.create_state_from_snapshot(snap_1)
    new_snap = sim.state.get_snapshot()

    # Categorize particles by typeid
    index_p, index_c, index_w, index_s = [], [], [], []
    for i, pid in enumerate(new_snap.particles.typeid):
        if pid == 0: index_p.append(i)
        elif pid == 1: index_c.append(i)
        elif pid == 2: index_w.append(i)
        elif pid == 3: index_s.append(i)

    # Define harmonic bond potentials
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['shell'] = dict(k=100 * np.sqrt(10), r0=2.1)
    harmonic.params['G'] = dict(k=10000, r0=2 * param['r_gen'])
    harmonic.params['spik'] = dict(k=100.0, r0=2 * param['r_ver'])

    # Define Lennard-Jones (LJ) interactions
    scale = 1

    nl = hoomd.md.nlist.Cell()
    lj = hoomd.md.pair.LJ(nlist=nl)

    # Define interacting LJ pairs
    lj.params[('P', 'P')] = {'sigma': 2 * param['r_gen'] * 1 ** (-1 / 6), 'epsilon': 0.01}
    lj.r_cut[('P', 'P')] = 2 * param['r_gen']

    lj.params[('P', 'W')] = {'sigma': ((param['r_gen'] + param['r_wall']) / 1) * 2 ** (-1 / 6), 'epsilon': scale * 0.1}
    lj.r_cut[('P', 'W')] = param['r_gen'] + param['r_wall']

    lj.params[('W', 'S')] = {'sigma': (param['r_wall'] + param['r_terminal']) * 2 ** (-1 / 6), 'epsilon': scale * 0.1}
    lj.r_cut[('W', 'S')] = param['r_wall'] + param['r_terminal']

    lj.params[('X', 'S')] = {'sigma': (param['r_ver'] + param['r_terminal']) * 2 ** (-1 / 6), 'epsilon': scale * 0.1}
    lj.r_cut[('X', 'S')] = param['r_ver'] + param['r_terminal']

    lj.params[('P', 'X')] = {'sigma': (param['r_gen'] + param['r_ver']) * 2 ** (-1 / 6), 'epsilon': scale * 0.1}
    lj.r_cut[('P', 'X')] = param['r_gen'] + param['r_ver']

    lj.params[('P', 'S')] = {'sigma': (param['r_gen'] + param['r_terminal']) * 2 ** (-1 / 6), 'epsilon': scale_2 * 0.1}
    lj.r_cut[('P', 'S')] = 5

    # Suppress LJ interactions between some types
    zero_pairs = [('S', 'S'), ('X', 'W'), ('W', 'W'), ('X', 'X')]
    for pair in zero_pairs:
        lj.params[pair] = {'sigma': 0, 'epsilon': 0}
        lj.r_cut[pair] = 0.01

    # Filter groups for potential methods (optional, here for clarity)
    group_all = hoomd.filter.All()
    group_p = hoomd.filter.Type('P')

    # Evaluate particle density (custom function assumed)
    new_snap = sim.state.get_snapshot()
    density_fun(new_snap)

    # Set up GSD writer
    gsd_writer = hoomd.write.GSD(
        trigger=hoomd.trigger.Periodic(50000),
        filename=output_file_genome,
        filter=hoomd.filter.All(),
        mode='wb'
    )
    sim.operations.writers.append(gsd_writer)

    # Apply Langevin dynamics to particles of type 'P'
    langevin = hoomd.md.methods.Langevin(filter=group_p, kT=0.1, alpha=0.2)

    # Create integrator
    integrator = hoomd.md.Integrator(dt=0.001, forces=[harmonic, lj], methods=[langevin])
    sim.operations.integrator = integrator

    # Run the simulation
    sim.run(total_time)

    # Final snapshot analysis
    new_snap = sim.state.get_snapshot()
    count_inside(new_snap)

    # Timing
    t1 = time.time() - t0
    print("Time elapsed:", t1)
    print("Timesteps per second:", sim.tps)


# To execute:
if __name__ == "__main__":
    main()
