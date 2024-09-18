from numpy.random import choice
from operator import itemgetter
from itertools import zip_longest
from pathlib import Path
from multiprocessing import Process, Queue, cpu_count


# Concatenate given strings by their overlap weight
def concat(a, b, weight):
    return a + b[weight:]


# Calculate difference between given strings
def hamming_with_padding(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip_longest(str1, str2))


class ACOSearch:
    def __init__(self,
                 instance,
                 repeats=5,
                 max_cycles=100,
                 ants=30,
                 chosen_ants=10,
                 alpha=1,
                 beta=1,
                 repetition_factor=0.2,
                 evaporation_rate=0.5,
                 pheromone_multiplication_value=2,
                 console=False,
                 logs=False,
                 save=False):
        self.instance = instance
        self.dna = instance.dna
        self.dna_length = instance.dna_length
        self.graph = instance.graph
        self.starter = instance.dna[:instance.oligomer_length]
        self.oligomer_length = instance.oligomer_length

        self.max_repeats = repeats
        self.max_cycles = max_cycles
        self.ants = ants
        self.chosen_ants = chosen_ants
        self.alpha = alpha
        self.beta = beta
        self.repetition_factor = repetition_factor
        self.evaporation_rate = evaporation_rate
        self.pheromone_multiplication_value = pheromone_multiplication_value
        self.result = None
        self.to_save = {}

        self.console = console  # print result
        self.logs = console and logs  # print iteration status
        self.save = save  # save to given location
        self.processes = cpu_count()

        self.keys = ["id", "oligomer", "overlap", "full_sequence", "used_oligomers", "arcs", "probability"]
        self.__aco()
        if self.save is not False:
            self.__save()

    def __str__(self):
        return "Sequence: {}\nQuality: {}\nFitness: {}%".format(
            self.result["sequence"],
            self.result["quality"],
            self.result["fitness"])

    def __save(self):
        Path("./data/output/" + '/'.join(self.save.split("/")[2:-1])).mkdir(parents=True, exist_ok=True)
        f = open("./data/output/" + '/'.join(self.save.split("/")[2:]) + ".txt", "w")
        for key, part in self.to_save.items():
            f.write("{} {} {}\n".format(key + 1, part["quality"], part["fitness"]))
        f.write("{}\n".format(self.result["sequence"]))
        f.write("{}\n".format(self.result["quality"]))
        f.write("{}\n".format(self.result["fitness"]))
        f.close()

    def __aco(self):
        repeats = 0

        for i in range(self.max_cycles):
            if self.console:
                print("Iteration: " + str(i + 1) + ", repeats: " + str(repeats))

            # Calculate number of tasks per process
            processes = self.processes
            tasks = self.ants
            k = 0
            while tasks > processes:
                k += 1
                tasks -= processes
            tasks_per_process = [k] * processes
            if tasks != 0:
                for j in range(tasks):
                    tasks_per_process[j] += 1

            qq = Queue()
            ant_processes = []
            for j in range(processes):
                ant_processes.append(Process(target=self.ant, args=(tasks_per_process[j], qq,)))
                ant_processes[j].daemon = True
                ant_processes[j].start()

            iteration_results = []
            finished_tasks = 0
            while True:
                iteration_results.append(qq.get())
                finished_tasks += 1
                if finished_tasks == self.ants:
                    break

            [x.join() for x in ant_processes]
            sorted_results = sorted(iteration_results, key=itemgetter("quality"), reverse=True)
            
            if self.logs:
                for nr, ant in enumerate(sorted_results):
                    print("{}-ant quality: {}".format(nr + 1, ant["quality"]))
                print("")

            if self.result is None or sorted_results[0]["quality"] > self.result["quality"]:
                self.result = sorted_results[0]
                repeats = 0
                if self.logs:
                    print("NEW BEST RESULT:\n{}\n".format(self))
            else:
                repeats += 1
                if repeats > self.max_repeats:
                    break

            # Pheromone evaporation
            for j in range(len(self.graph.weights)):
                self.graph.pheromones[j] *= 1 - self.evaporation_rate

            # New pheromones
            for result in sorted_results[:self.chosen_ants]:
                for arc in result["arcs"]:
                    self.graph.pheromones[arc] += result["quality"] * self.pheromone_multiplication_value

            if self.save:
                self.to_save.update({i: self.result.copy()})

        if self.console:
            print("\nRESULT:\n{}".format(self))

    def ant(self, repetitions, qq):
        for i in range(repetitions):
            begin = self.graph.get_vertex_of_label(self.starter, "left")
            initial_candidate = dict(zip(self.keys[:-1], [
                begin["id"],
                begin["label"],
                0,
                begin["label"],
                0,
                [],
                0
            ]))

            # INITIATE CURRENT CANDIDATE
            current_candidate = initial_candidate.copy()

            used_oligomers = []
            while not self.__searching_stopping_condition(current_candidate["full_sequence"]):
                # NEIGHTBORHOOD
                neighborhood = self.__get_neighborhood(current_candidate, used_oligomers)
                current_candidate = choice(neighborhood, 1, p=[neighbor["probability"] for neighbor in neighborhood])[0]
                used_oligomers.append(current_candidate["id"])

            qq.put({
                "sequence": current_candidate["full_sequence"],
                "quality": current_candidate["used_oligomers"] / (
                        self.instance.dna_length - self.instance.oligomer_length + 1),
                "fitness": (self.dna_length - hamming_with_padding(self.dna, current_candidate["full_sequence"])) / self.dna_length * 100,
                "arcs": current_candidate["arcs"]
                })

    def __searching_stopping_condition(self, sequence):
        return len(sequence) >= self.dna_length

    # Get neighborhood of the current candidate
    def __get_neighborhood(self, current_candidate, used_oligomers):
        neighbors = self.graph.get_neighbors_by_arcs_out_of_idx(current_candidate["id"])
        if len(neighbors) == 0:
            for idx, label in enumerate(self.instance.spectrum):
                neighbors.append({"id": idx, "label": label})

        neighborhood = []
        weights_sum = 0
        for neighbor in neighbors:
            neighbor_possibilities = self.__get_next_neighbor_info(current_candidate, neighbor, used_oligomers)
            for candidate in neighbor_possibilities:
                weights_sum += candidate["probability"]
                neighborhood.append(candidate)

        for neighbor in neighborhood:
            neighbor["probability"] = float(neighbor["probability"]) / weights_sum

        return neighborhood

    # Get next neighbor info of the current candidate
    def __get_next_neighbor_info(self, current_candidate, neighbor, used_oligomers):
        neighbor_possibilities = []
        preffered_arcs = self.graph.find_weights(current_candidate["id"], neighbor["id"])

        arcs = {
            "weights": preffered_arcs["weights"],
            "pheromones": preffered_arcs["pheromones"],
            "arc_ids": preffered_arcs["arc_ids"]
        }

        if len(arcs) == 0:
            arcs = {"weights": [], "pheromones": [], "arc_ids": []}
        for overlap, pheromone, arc in zip(arcs["weights"], arcs["pheromones"], arcs["arc_ids"]):
            oligomer = neighbor["label"]
            full_sequence = concat(current_candidate["full_sequence"], oligomer, overlap)
            overlap_value = overlap ** self.beta
            new_arcs = list(current_candidate["arcs"])
            new_arcs.append(arc)
            candidate = dict(zip(self.keys, [
                neighbor["id"],
                oligomer,
                overlap,
                full_sequence,
                current_candidate["used_oligomers"] + 1,
                new_arcs,
                (pheromone ** self.alpha * overlap_value) if pheromone != 0 else (0.5 ** self.alpha * overlap_value)
            ]))

            if candidate["id"] in used_oligomers:
                candidate["probability"] = candidate["probability"] ** self.repetition_factor

            neighbor_possibilities.append(candidate)

        return neighbor_possibilities
