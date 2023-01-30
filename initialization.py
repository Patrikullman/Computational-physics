""" Methods for generating new random starting candidates. """
import numpy as np
from numpy import random
from ase import Atoms, Atom


def get_random_direction():
    """ Returns a random vector with norm 1. """
    phi, theta = random.rand(2) * np.pi
    phi *= 2
    direction = np.array([np.sin(theta) * np.cos(phi),
                          np.sin(theta) * np.sin(phi), np.cos(theta)])
    return direction


def min_distance_fulfilled(
        atom,
        atomic_numbers,
        positions,
        bond_lengths,
        min_distance_ratio=1):
    """ Checks if the distance between atom and atoms described by
        atomic_numbers and positions is large enough in regards to the provided
        bond lengths multiplied by min_distance_ratio.

        Parameters:

        atom : Atom object containing a position and atom number.
        atomic_numbers : The atomic numbers of the atoms to check distance to.
        positions : The positions of the atoms to check distance to.
        bond_lengths : A dictionary containing the bond lengths
        min_distance_ratio : ratio of the bond length allowed between atoms.
    """
    for number, position in zip(atomic_numbers, positions):
        distance = np.linalg.norm(position - atom.position)
        if distance < min_distance_ratio * bond_lengths[(atom.number, number)]:
            return False

    return True


class StartGenerator(object):
    """ Class used to generate random starting candidates.
        The candidates are generated by iteratively trying to add
        a neighbour to a random atom in the cluster. The bond length
        of the atom and the neighbour is fixed. Minimum distance to
        other atoms can be adjusted when calling get_new_candidate()
        (default is bond_length).

        Parameters:

        slab : The atoms object describing the super cell to
        optimize within.
        atom_numbers : A list of the atomic numbers that needs
        to be optimized.
        bond_lengths : A dictionary containing the bond lengths
        between all different atom types used in the cluster.
    """

    def __init__(self, slab, atom_numbers,
                 bond_lengths):
        # the generator does not work if the slab contains atoms right now.
        # The slab is only used to get cell and boundary conditions.
        assert(len(slab) == 0)
        self.slab = slab
        self.atom_numbers = atom_numbers
        self.blmin = bond_lengths

    def get_new_candidate(
            self,
            max_tries_to_place_atom=10,
            min_distance_ratio=1.0):
        """ Returns a new candidate (None if unsuccessful).

            Parameters:

            max_tries_to_place_atom : how many random angles to try before
            switching to a different neighbour.
            min_distance_ratio : describes the minimum distance to non-neighbours expressed
            in ratio of the bond length

            Returns:

            Atoms object if successful, else None.
        """
        atom_numbers = self.atom_numbers

        # The ordering is shuffled so different atom
        # types are added in different order for different
        # candidates.
        random.shuffle(atom_numbers)

        # Create a cluster containing the first atom (positioned at 0,0,0)
        # use same properties as for the slab. As the cluster will be
        # centered in the cell at the end the position of the first atom
        # is not important.
        cell = self.slab.get_cell()
        pbc = self.slab.get_pbc()
        cluster = Atoms(numbers=[atom_numbers[0]],
                        pbc=pbc, cell=cell)

        # Add remaining atoms
        for atomic_number in atom_numbers[1:]:
            success = self.add_atom(atomic_number, cluster,
                                    max_tries=max_tries_to_place_atom,
                                    min_distance_ratio=min_distance_ratio)
            # If a new atom could not be placed anywhere, return None.
            # Construction of the candidate failed.
            if not success:
                return None

        # Sort to get atom numbers in correct order, positions are reordered as
        # appropriate
        orderd_numbers, orderd_positions = zip(*sorted(
            zip(cluster.numbers, cluster.positions),
            key=lambda x: x[0]))
        cluster.set_atomic_numbers(orderd_numbers)
        cluster.set_positions(orderd_positions)
        cluster.center()

        return self.slab + cluster

    def add_atom(
            self,
            atomic_number,
            cluster,
            max_tries=10,
            min_distance_ratio=1.0):
        """ Adds an atom to the cluster by iterating through all cluster atoms
            in random order, trying to place it at bond_length's distance from
            this atom and at min_distance_ratio x bond_length distance from all
            other atoms using the place_neighbour method. If it succeeds the
            atom is added to the cluster and True is returned, else it returns
            False.

            Parameters:

            atomic_number : {int, str} Atomic number of the atom to add
            cluster : Atoms object to add an atom to.
            max_tries : {positive int} The number of random directions to try
            and add the new atom at every neighbour before switching to the next one.
            min_distance_ratio : {positive float} Ratio of bond_length allowed
            as minimum distance to other atoms (not the one it is placed by).
            Only use values less than 2.

            Returns:

            True if successful (an atom has been added), else false.
        """
        # Try to place the new atom at every atom present in cluster (random order), return if
        # successful.
        for neighbour_index in random.permutation(
                cluster.get_number_of_atoms()):
            new_atom = self.place_neighbour(
                atomic_number,
                cluster,
                neighbour_index,
                max_tries=max_tries,
                min_distance_ratio=min_distance_ratio)
            if new_atom is not None:
                cluster += new_atom
                return True

        # Did not manage to place the atom.
        # Nothing is added to the cluster.
        return False

    def place_neighbour(
            self,
            atomic_number,
            cluster,
            neighbour_index,
            max_tries=10,
            min_distance_ratio=1.0):
        """ Tries to add the new atom to the cluster by placing it at
            bond_length's distance from the atom at neighbour_index and at
            min_distance_ratio x bond_length distance from all other atoms. If
            it succeeds a Atom object is return, else None.

            Parameters:

            atomic_number : {int, str} Atomic number of the new atom
            cluster : Atoms object to add an atom to.
            max_tries : {positive int} The number of random directions to try
            and add the new atom at every neighbour before switching to the next one.
            min_distance_ratio : {positive float} Ratio of bond_length allowed
            as minimum distance to other atoms (not the one it is placed by).
            Only use value less than 2.

            Returns:

            Atom object if successful, else None.
        """
        new_atom = Atom(symbol=atomic_number)
        counter = 0
        found_position = False
        while not found_position and counter < max_tries:
            counter += 1
            # Start at neighbour position
            new_atom.position = cluster.positions[neighbour_index].copy()
            # Add random vector with norm of the provided bond_length
            bond_length = self.blmin[(
                atomic_number, cluster.numbers[neighbour_index])]
            new_atom.position += get_random_direction() * bond_length

            # Check if other distances are large enough
            atoms_to_check = np.array([True] * cluster.get_number_of_atoms())
            # neighbour distance may be smaller than other distances
            atoms_to_check[neighbour_index] = False
            found_position = min_distance_fulfilled(
                new_atom,
                cluster.numbers[atoms_to_check],
                cluster.positions[atoms_to_check],
                self.blmin,
                min_distance_ratio)

        if counter == max_tries:
            # Could not place atom
            return None

        return new_atom


if __name__ == '__main__':
    from ase.ga.utilities import closest_distances_generator
    from ase.ga.data import PrepareDB
    from ase.calculators.lj import LennardJones
    from ase.optimize import BFGS
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('natoms', type=int,
                        help='Number of atoms in the cluster')
    defval = 20
    parser.add_argument('-p', '--population', type=int, default=defval,
                        help='Size of starting population (default: {})'.format(defval))
    parser.add_argument('-v', '--visualize', action='store_true',
                        help='Visualize the population with the ase gui')
    args = parser.parse_args()

    db_file = 'gadb.db'
    calc = LennardJones(epsilon=0.46252798, sigma=2.57488462)

    # Make random Na clusters.
    # All atoms are at least one covalent radius apart.
    atom_numbers = [11] * args.natoms
    cd = closest_distances_generator(atom_numbers=atom_numbers,
                                     ratio_of_covalent_radii=0.7)
    sg = StartGenerator(slab=Atoms(),  # pbc from this object are used
                        atom_numbers=atom_numbers,
                        bond_lengths=cd)
    sg.slab.set_cell([16, 16, 16])

    # generate the starting population
    population_size = args.population
    starting_population = [sg.get_new_candidate() for i in range(population_size)]

    for cluster in starting_population:
        # Relax the cluster with a homemade LJ-potential
        cluster.set_calculator(calc)
        dyn = BFGS(cluster, trajectory=None, logfile='qn.log')
        dyn.run(fmax=0.05, steps=1000)

    if args.visualize:
        from ase.visualize import view
        view(starting_population)

    # create the database to store information in
    d = PrepareDB(db_file_name=db_file,
                  simulation_cell=sg.slab,
                  stoichiometry=atom_numbers)

    for a in starting_population:
        d.add_unrelaxed_candidate(a)

