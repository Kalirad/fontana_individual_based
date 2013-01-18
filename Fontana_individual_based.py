
import numpy as np
import numpy.random as rnd
import pickle
import copy


class Population(object):

    def __init__(self, N, r, u, n_alleles, starting_seq):
        """
        Initiate the population object as a list of Genotype objects with seq "000" and 4 loci.

        :param N: population size
        :param r: recombination rate
        :param u: mutation rate
        :param n_alleles: number of alleles
        :param starting_seq: the sequence of population at the start of simulation
        :type starting_seq: string
        """
        assert type(starting_seq) == str
        assert 0 <= u <= 1
        assert r >= 0
        assert r <= 0.5
        self.population = []
        self.N = N
        self.r = r
        self.u = u
        self.n_alleles = n_alleles
        self.n_loci = len(starting_seq)
        for i in range(N):
            self.population.append(starting_seq)

    def add_network_to_population(self, lst):
        """
        Add a list of sequences to the network.

        :param lst: list of sequences in the network
        :type lst: list
        """
        assert type(lst) == list
        self.network_list = lst


    def contain_sequence(self, offspring):
        """
        Check to see whether the genotype is in the network
        """
        assert type(offspring) == str
        if offspring in self.network_list:
            return True


    def define_modules(self, module_1, module_2):
        """
        Define the modules in the network and make sure that they belong to the network.

        :param module_1: a list of sequences forming the first module in the network.
        :param module_2: a list of sequences forming the second module in the network.
        """
        assert type(module_1) == list
        assert type(module_2) == list
        assert len(module_1 + module_2) == len(self.network_list)

        # Make sure that modules belong to the network
        for i in range(len(module_1)):
            assert module_1[i] in self.network_list
        for i in range(len(module_2)):
            assert module_2[i] in self.network_list
        self.module_1 = module_1
        self.module_2 = module_2


    '''
    def define_module(self, *args):
        
        self.modules = {}
        i = 1
        for arg in args:
            assert type(arg) == list
            self.module[i] = arg
            i += 1
    '''

    @staticmethod
    def cross_over(seq1, seq2, n_alleles):
        """
        Peform cross over between two sequences.
        Recombination happens at a random position, randint drawn from n_loci.
        One of the two recombinant sequences is being chosen by chance (50:50).

        :param seq1:
        :param seq2:
        :param n_alleles:
        :return rec: the product of recombination
        """
        n_loci = len(seq1)
        pos = rnd.randint(1, n_loci)
        if rnd.randint(0,2) == 0:
            rec = seq1[0:pos] + seq2[pos:]
            return rec
        else:
            rec = seq2[0:pos] + seq1[pos:]
            return rec


    @staticmethod
    def mutate(seq, n_alleles):
        """
        Mutate "seq" by the probability of "u".

        :param u: mutation rate
        :param n_alleles: number of alleles
        :param seq:
        :type seq: string
        :return mutant: the product of mutation
        """
        assert type(seq) == str
        n_loci = len(seq1)    
        #Determihe soon-too-be-mutated locus at random
        pos = rnd.randint(0, n_loci)
        #Draw the mutation from available alleles
        mut = rnd.randint(0, n_alleles)
        #Ensure that the mutation is different from the current allele at locus of interest
        if mut == int(seq[pos]):
            #If the mutant allele is the same as current allele, draw a new mutation at random from availible alleles untill the condition is met.
            while mut == int(seq[pos]):
                mut = rnd.randint(0, n_alleles)
            else:
                mutant = seq[0:pos] + str(mut) + seq[pos+1:]
                return mutant
        else:
            mutant = seq[0:pos] + str(mut) + seq[pos+1:]
            return mutant
     

    @staticmethod
    def generate_offspring(seq1, seq2, n_alleles, u, r):
        """
        Generate offspring using recombination and then mutation.

        :return offspring: sequences which has undergone recombination and mutation.
        """
        assert len(seq1) == len(seq2)
        if rnd.rand() <= r:
            rec = Population.cross_over(seq1, seq2, n_alleles)
            if rnd.rand() <= u:
                offspring = Population.mutate(n_alleles, rec)
                return offspring
            else:
                offspring = rec
                return offspring
        else:
            if rnd.randint(0,2) == 0:
                rec = seq1
            else:
                rec = seq2
        if rnd.rand() <= u:
            offspring = Population.mutate(n_alleles, rec)
            return offspring
        else:
            offspring = rec
            return offspring

    def get_next_generation(self):
        """
        Create the next generation. The events occur in the following order:
        
        (i) Picking two individulas from the population at random
        (ii) recombining and mutate the individulas proportional to their r and u
        (iii) See whether the offpring is a part of the network and if so add it to the next generation

        :return next_generation: a list of sequences comprising the next generation
        """
        next_generation = copy.deepcopy(self)
        next_generation.population = []
        while len(next_generation.population) < self.N:
            #(i)Picking
            ind_1 = self.population[rnd.randint(0, self.N)]
            ind_2 = self.population[rnd.randint(0, self.N)]
            #(ii) Generating the offspring
            offspring = Population.generate_offspring(ind_1, ind_2, self.n_alleles, self.u, self.r)
            for i in range(len(self.network_list)):
                #(iii) Check whether the offspring is in the network
                if self.contain_sequence(offspring):
                    next_generation.population.append(offspring)
        return next_generation


    def get_stats(self):
        """
        Count the number of individulas at each node in the network.
        """
        self.stats = {}

        for i in self.network_list:
            self.stats[i] = self.population.count(i)

        # Count the number of individuals in the first module
        module_1_genotype_count = 0
        for i in range(len(self.module_1)):
            module_1_genotype_count += self.stats[self.module_1[i]]

        # Count the number of individuals in the second module
        module_2_genotype_count = 0
        for i in range(len(self.module_2)):
            module_2_genotype_count += self.stats[self.module_2[i]]

        # colculate the genotype count by subtractiing the #individuals in modules form each other
        genotype_count = module_1_genotype_count - module_2_genotype_count
        #self.stats['genotype_count'] = genotype_count
        #self.stats['<m>'] = np.divide(float(genotype_count),float(self.N))
        return self.stats



class Evolution(object):
    """
    This object manipulates the object population to simulate over generations.
    """
    def __init__(self, population, n_generations, period = 1, verbose = False, name = "simulation"):
        """
        Save the population objects in the course of simulation.

        :param population: starting population
        :param n_generations: the number of generations
        :param period: period of pickling
        """
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

