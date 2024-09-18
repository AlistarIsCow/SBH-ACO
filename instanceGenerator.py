import math
import random

from graph import Graph


class Instance:
    def __init__(self, data):
        random.seed()

        if len(data.keys()) == 4:
            self.dna = data["dna"]
            self.spectrum_with_errors = data["spectrum"]
            self.oligomer_length = len(data["spectrum"][0])
            self.dna_length = len(self.dna)
        else:
            self.dna_length = data["dna_length"]
            self.oligomer_length = data["oligomer_length"]
            self.nucleobases = data["nucleobases"]
            self.negatives = data["negatives"]
            self.positives = data["positives"]

            self.dna = ""
            self.spectrum = []
            self.spectrum_with_errors = []

            self.__generate_dna()
            self.__generate_spectrum()
            self.__add_errors_to_spectrum()

        self.save_graph = data["save_graph"]
        self.console = data["console"]
        self.graph = None
        self.__generate_graph()
        if self.console:
            print(self)

    def __str__(self):
        return "\nInstance:\n" + \
               "\tDNA: {}\n".format(self.dna) + \
               "\tSpectrum: {}\n".format(self.spectrum) + \
               "\tSpectrumWithErrors: {}\n".format(self.spectrum_with_errors)

    def __generate_dna(self):
        for i in range(self.dna_length):
            self.dna += random.choice(self.nucleobases)

    def __generate_spectrum(self):
        for i in range(0, len(self.dna) - self.oligomer_length + 1):
            oligomer = self.dna[i:i + self.oligomer_length]
            if oligomer not in self.spectrum:
                self.spectrum.append(oligomer)
        random.shuffle(self.spectrum)

    def __add_errors_to_spectrum(self):
        self.spectrum_with_errors = self.spectrum.copy()
        length = len(self.spectrum)
        nr_neg = math.ceil(length * self.negatives)
        nr_pos = math.ceil(length * self.positives)
        to_pop = []
        for i in range(nr_neg):
            iteration = 0
            while iteration < 5:
                iteration += 1
                pop = random.randint(0, len(self.spectrum_with_errors) - 1)
                if pop not in to_pop:
                    to_pop.append(pop)
                    break
        to_pop.sort(reverse=True)
        to_add = []
        for i in range(nr_pos):
            iteration = 0
            while iteration < 5:
                iteration += 1
                oligomer = "".join([random.choice(self.nucleobases) for x in range(self.oligomer_length)])
                if oligomer not in self.spectrum_with_errors and oligomer not in to_add:
                    to_add.append(oligomer)
                    break

        for pop in to_pop:
            self.spectrum_with_errors.pop(pop)
        for add in to_add:
            self.spectrum_with_errors.insert(random.randint(0, len(self.spectrum_with_errors) - 1), add)

        if self.dna[:self.oligomer_length] not in self.spectrum_with_errors:
            self.spectrum_with_errors.insert(random.randint(0, len(self.spectrum_with_errors) - 1), self.dna[:self.oligomer_length])

    def __generate_graph(self):
        self.graph = Graph(self.spectrum_with_errors, self.oligomer_length, self.save_graph)
