import copy
from instanceGenerator import Instance
from ACOSearch import ACOSearch
from os import listdir
from os.path import isfile, join
import time


def process(instance, sequence_file, ants, evaporation_rate, alpha, beta):
    print("Processing: {} ants-{} evaporation_rate-{} alpha-{} beta-{} ".format(
        sequence_file.split("/")[-1].split('.')[0], ants, evaporation_rate, alpha, beta), end="")
    ACOSearch(
        instance=instance,
        repeats=5,
        max_cycles=50,
        ants=ants,
        chosen_ants=10,
        alpha=alpha,
        beta=beta,
        repetition_factor=0.3,
        evaporation_rate=evaporation_rate,
        pheromone_multiplication_value=2,
        console=False,
        logs=False,
        save=sequence_file.split('.')[1] + "-ants_{}-evrate_{:.0f}-alpha_{}-beta_{}".format(ants,
                                                                                            evaporation_rate * 100,
                                                                                            alpha, beta)
    )
    print("DONE")


def read_data(sequence_path, spectrum_path):
    sequence_files = sequence_path
    spectrum_files = spectrum_path
    sq = open(sequence_files, "r")
    sp = open(spectrum_files, "r")
    dna = sq.readline()
    spectrum = []
    for i, oligomer in enumerate(sp):
        if i > 1:
            spectrum.append(oligomer[:-1])
    spectrum.insert(0, dna[:len(spectrum[0])])
    instance = Instance({
        "dna": dna,
        "spectrum": spectrum,
        "save_graph": False,
        "console": False
    })

    process(copy.deepcopy(instance), sequence_files, 40, 0.3, 1, 7)


def generate_data():
    instance = Instance({
        "dna_length": 500,
        "oligomer_length": 7,
        "nucleobases": ("A", "C", "G", "T"),
        "negatives": 0.1,
        "positives": 0.1,
        "save_graph": True,
        "console": True
    })
    ACOSearch(
        instance=instance,
        repeats=5,
        max_cycles=100,
        ants=40,
        chosen_ants=10,
        alpha=1,
        beta=7,
        repetition_factor=0.3,
        evaporation_rate=0.3,
        pheromone_multiplication_value=2,
        console=True,
        logs=True,
        save=False
    )

def main():
    start = time.time()
    data = [["./data/500_5.seq", "./data/500_5"],
            ["./data/500_10.seq", "./data/500_10"],
            ["./data/500_20.seq", "./data/500_20"]]
    for files in data:
        read_data(files[0], files[1])
    #generate_data()
    print("Total time: {}".format(time.time() - start))


if __name__ == "__main__":
    main()
