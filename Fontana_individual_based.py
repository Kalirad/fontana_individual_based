
import numpy as np
import numpy.random as rnd
import pickle
import copy


class Genotype(object):
    """
    A sequence.

    **Attributes:**

    * **n_alleles** (*int* between 2 and 10) -- number of alleles per locus (size of the alphabet)
    * **n_loci** -- number of loci, given by ``len(seq)``
    * **seq** (*str* of integers in ``range(n_alleles)``) -- sequence of length **n_loci**
    """

    def __init__(self, seq, n_alleles):
        """
        Initialize a Genotype object with an input **seq** and **n_alleles**.

        :param seq: sequence of length **n_loci**
        :type seq: str
        :param n_alleles: number of alleles per locus
        :type n_alleles: int
        :return: Genotype object

        >>> x = Genotype("00101", 3)
        >>> print x.seq, x.n_loci, x.n_alleles, x.int_repr
        00101 5 3 10
        """
        assert n_alleles > 1 and n_alleles < 11
        assert type(seq) == str
        for allele in seq:
            assert int(allele) < n_alleles
        self.seq = seq
        self.n_alleles = n_alleles

    @property
    def n_loci(self):
        return len(self.seq)


    def mutate(self, u):
        """
        Muate a sequence with the probability of u.
        """
        assert 0 <= u <= 1
        if rnd.rand() <= u:
            mut = list(self.seq)
            index = rnd.randint(0, self.n_loci)
            while True:
                mutation = rnd.randint(0, self.n_alleles)
                # Ensures that mutation is always different from the original sequence
                if mutation != mut[index]:
                    break
            mut[index] = str(mutation)
            self.seq = "".join(mut)
            return self
    
    @staticmethod
    def cross_over(seq1, seq2, n_alleles, n_loci):
        """
        Peform cross over between two sequences.
        Recombination happens at a random position, randint drawn from n_loci.
        One of the two recombinant sequences is being chosen by chance (50:50).
        """
        pos = rnd.randint(1, n_loci)
        rec1 = seq1[0:pos] + seq2[pos:]
        rec2 = seq2[0:pos] + seq1[pos:]
        if rnd.randint(0,2) == 0:
            return Genotype(rec1, n_alleles)
        else:
            return Genotype(rec2, n_alleles)
        
    
    def generate_offspring(self, other, u, r):
        """
        Generate offspring using recombination and then mutation.
        """
        assert type(self) == Genotype
        assert type(other) == Genotype
        assert self.n_loci == other.n_loci
        assert self.n_alleles == other.n_alleles
        assert r >= 0
        assert r <= 0.5
        
        offspring = copy.deepcopy(self)
        if rnd.rand() <= r:
            rec = Genotype.cross_over(self.seq, other.seq, self.n_alleles, self.n_loci)
            offspring.seq = rec.seq
        else:
            if rnd.randint(0,2) == 0:
                offspring.seq = other.seq
        offspring.mutate(u)
        return offspring



class Population(object):
    
    def __init__(self, N, r, u, starting_seq, network_list):
        """
        Initiate the population object as a list of Genotype objects with seq "000" and 4 loci.
        """
        assert type(network_list) == list
        assert type(starting_seq) == str
        self.population = []
        for i in range(N):
            self.population.append(Genotype(starting_seq, 4))
        self.N = N
        self.r = r
        self.u = u
        self.network_list = network_list
        
    def is_in_network(self, offspring):
        """
        Check to see whether the genotype is in the network
        """
        assert type(offspring) == Genotype
        for i in range(len(self.network_list)):
            if self.network_list[i] == offspring.seq:
                return True
            
            
    def get_next_generation(self):
        """
        Create the next generation. The events occur in the following order: 
        (i) Picking two individulas from the population at random
        (ii) recombining and mutate the individulas proportional to their r and u
        (iii) See whether the offpring is a part of the network and if so add it to the next generation 
        """
        next_generation = copy.deepcopy(self)
        next_generation.population = []
        while len(next_generation.population) < self.N:
            #(i)Picking
            ind_1 = self.population[randint(0, self.N)]
            ind_2 = self.population[randint(0, self.N)]
            #(ii) Generating the offspring
            offspring = ind_1.generate_offspring(ind_2, self.u, self.r)
            for i in range(len(self.network_list)):
                #(iii) Check whether the offspring is in the network
                if self.is_in_network(offspring) == True:
                    next_generation.population.append(offspring)
        return next_generation 
                

    def get_stats(self):
        """
        Count the number of individulas at each node in the network.
        It creats a list of sequences (step(i)). Then counts the number of each sequence in the list and adds the sequence and its count
        to the self.stats dictionary (step(ii)).
        """
        self.stats = {}
        #(i)
        sequences = []
        for i in range(len(self.population)):
            sequences.append(self.population[i].seq)
        #(ii)
        for i in self.network_list:
            self.stats[i] = sequences.count(i)
        genotype_count = (self.stats['001'] + self.stats['101'] + self.stats['201']) - (self.stats['010'] + self.stats['011'] + self.stats['012'])
        self.stats['genotype_count'] = genotype_count
        self.stats['<m>'] = np.divide(float(genotype_count),float(self.N))
        return self.stats
            


class Evolution(object):
    """
    This object manipulates the object population to simulate over generations.
    """
    def __init__(self, population, n_generations, period = 1, verbose = False, name = "simulation"):
        gen = 0
        population.get_stats()
        self.curr_population = population
        self.evolution = {0: population.stats}
        for i in range(n_generations):
            next_generation = self.curr_population.get_next_generation()
            self.curr_population = next_generation
            gen += 1
            if gen % period == 0:
                next_generation.get_stats()
                self.evolution.update({gen: next_generation.stats})
                file = open(name + "_" + "r_" + str(next_generation.r) + "_u_" + str(next_generation.u) + "_" + str(gen) + ".pickle", "w")
                pickle.dump(next_generation, file)
                file.close()
                if verbose:
                    print gen


