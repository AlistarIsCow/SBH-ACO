import igraph as ig


class Graph:
    def __init__(self, spectrum, initial_oligomer_length, save=False):
        self.spectrum = spectrum
        self.initial_oligomer_length = initial_oligomer_length
        self.vertices_count = len(spectrum)
        self.save = save

        self.arcs = {}
        self.weights = []
        self.pheromones = []
        self.graph = None

        self.__generate_graph()
        self.__build_graph()
        self.__save_graph()

    def __generate_graph(self):
        k = 0
        for i in range(self.vertices_count):
            for j in range(self.vertices_count):
                weights = self.__calc_weights(self.spectrum[i], self.spectrum[j])
                for weight in weights:
                    if (i, j) in self.arcs:
                        self.arcs[(i, j)].append(k)
                    else:
                        self.arcs[(i, j)] = [k]
                    k += 1
                    self.weights.append(weight)
                    self.pheromones.append(0.0)

    def __build_graph(self):
        self.graph = ig.Graph(n=self.vertices_count, edges=self.arcs, edge_attrs={'label': self.weights,
            'pheromone': self.pheromones}, vertex_attrs={'label': self.spectrum}, directed=True)

    def __calc_weights(self, oligo1, oligo2):
        weights = []
        length1 = len(oligo1)
        length2 = len(oligo2)
        diff = length1 - length2
        begin = length1 - self.initial_oligomer_length + 1
        for i in range(begin, length1):
            for j in range(len(oligo1[i:])):
                if oligo1[i:][j] != oligo2[j]:
                    break
                if j == length1 - i - 1:
                    # Number of overlaps
                    weights.append(min(length1, length2) - i + (diff if diff > 0 else 0))
        return weights

    def __save_graph(self):
        if not self.save:
            return
        layout = self.graph.layout("kk")
        out = ig.plot(self.graph, layout=layout, margin=100)
        out.save("graph.png")

    # Get vertex by pattern and edge of label
    def get_vertex_of_label(self, pattern, side):
        for i, label in enumerate(self.spectrum):
            if side == "left" and pattern == label[:len(pattern)]:
                return {"id": i, "label": label}
            elif side == "right" and pattern == label[-len(pattern):]:
                return {"id": i, "label": label}
        return None

    # Get label by vertex index
    def get_vertex_of_idx(self, idx):
        return {"id": idx, "label": self.spectrum[idx]}

    # Get neighborhood of vertex
    def get_neighbors_of_idx(self, idx):
        result = []
        for neighbor in self.graph.neighbors(idx):
            result.append(self.get_vertex_of_idx(neighbor))
        return result

    # Get neighborhood of vertex by IN arcs
    def get_neighbors_by_arcs_in_of_idx(self, idx):
        result = []
        for neighbor in self.graph.neighbors(idx, mode="in"):
            result.append(self.get_vertex_of_idx(neighbor))
        return result

    # Get neighborhood of vertex by OUT arcs
    def get_neighbors_by_arcs_out_of_idx(self, idx):
        try:
            result = []
            for neighbor in self.graph.neighbors(idx, mode="out"):
                result.append(self.get_vertex_of_idx(neighbor))
            return result
        except:
            return []

    # Get all arcs of vertex
    def get_incident_arcs_of_idx(self, idx):
        result = []
        for neighbor in self.graph.vs[idx].incident():
            result.append(self.get_vertex_of_idx(neighbor))
        return result

    # Get all arcs of vertex by IN ars
    def get_incident_arcs_by_arcs_in(self, idx):
        result = []
        for neighbor in self.graph.vs[idx].in_arcs():
            result.append(self.get_vertex_of_idx(neighbor))
        return result

    # Get all arcs of vertex by OUT arcs
    def get_incident_arcs_by_arcs_out(self, idx):
        result = []
        for neighbor in self.graph.vs[idx].out_arcs():
            result.append(self.get_vertex_of_idx(neighbor))
        return result

    # Find all weights from vertex A to vertex B
    def find_weights(self, idx_A, idx_B):
        weights = []
        pheromones = []
        ids = []
        if (idx_A, idx_B) in self.arcs:
            for idx in self.arcs[(idx_A, idx_B)]:
                ids.append(idx)
                weights.append(self.weights[idx])
                pheromones.append(self.pheromones[idx])

        return {"weights": weights, "pheromones": pheromones, "arc_ids": ids}