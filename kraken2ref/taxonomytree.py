import logging
from kraken2ref.taxonlevel import TaxonLevel, find_parent

class TaxonomyTree:
    """Class to encapsulate taxonomy tree and handle related exceptions/standard behaviours
    """
    def __init__(self, nodes: list = None, tree = None):
        """Initialiser. Essentially populates data into instance using methods.

        Args:
            nodes (list, optional): Node list to build tree from. Defaults to None.
            tree (dict, optional): Previous TaxonomyTree.graph object to initialise from. Defaults to None.

        """
        if nodes == None and tree == None:
            raise AttributeError("No data found! Please provide either a list of indexed nodes, or a prebuilt tree.")

        ## ignore tree if both nodes and tree given
        if nodes and tree:
            tree = None

        if nodes:
            self.nodes = sorted(nodes)
            self.graph = self.build_from(self.nodes)

        if tree:
            self.nodes = sorted(list(tree.keys()))
            self.graph = tree

        self.root = self.nodes[0]
        self.root_idx = self.root[0]

        self.max_lvl = max([TaxonLevel(t) for (i, t) in self.nodes]).lvl

        self.leaf_lvls = sorted(set([t for (i, t) in self.graph.keys() if len(self.graph[(i, t)]) == 0]))
        self.subterminal_lvls = [(TaxonLevel(t) - 1).lvl for t in self.leaf_lvls]
        self.complexity = self.get_tree_complexity(self.graph)

        self.subterminal_nodes = sorted([node for node in self.nodes if node[1] in self.subterminal_lvls and len(self.graph[node]) > 0])
        self.leaf_nodes = sorted([node for node in self.nodes if node[1] in self.leaf_lvls and len(self.graph[node]) == 0])
        if self.complexity > 0:
            self.subgraphs = self.decompose_tree(self.graph, [])
        else:
            self.subgraphs = None

    def build_from(self, node_list: list):
        """Function to build a graph representation from a list of indexed nodes.

        Args:
            node_list (list(tuples)): List of nodes; each node is represented as a tuple,
                    where:
                        node[0] = index of that node in the original data
                        node[1] = taxon level of that node ("S"/"S1"/etc)

        Returns:
            graph (dict): A dictionary representation of the graph, where:
                        key = Source node
                        value = Target node
                        (Nodes are represented as described above)
        """

        ## initialise required data structures for graph building
        graph = {}

        plain_levels = []  ## e.g. ["S", "S1", "S2", "S3", "S3"]

        ## levels_as_taxa is used to backtrack
        ## good to have it computed once as will probably be needed often
        levels_as_taxa = []  ## same as above except as list of TaxonLevel objects

        ## iterate over node_list list
        ## (looks like: [(11,"S"), (12,"S1)", (13,"S2"), (14,"S3"), (15,"S3")] )
        ## to populate data structures
        for indexed_node in node_list:
            graph[indexed_node] = []
            plain_levels.append(indexed_node[1])
            levels_as_taxa.append(TaxonLevel(indexed_node[1]))

        ## set up counters for list trawling
        start = 0

        ## trawl over list to find edges in graph
        while start < len(plain_levels)-1:
            end = start + 1
            left_node = TaxonLevel(plain_levels[start])
            right_node = TaxonLevel(plain_levels[end])

            ## if left node less than right node: that's an edge
            if left_node < right_node:
                ## generalise this, do not depend on offset
                graph[node_list[start]].append(node_list[end])

            ## if left node equals right node: oops! find the closest parent
            if left_node >= right_node:
                sublist = levels_as_taxa[:end]
                parent = right_node - 1
                parent_idx = find_parent(sublist, parent)
                ## generalise this, do not depend on offset
                graph[node_list[parent_idx]].append(node_list[end])
            start += 1

        return graph

    def get_tree_complexity(self, graph: dict):
        """Get a measure of graph "complexity" based on multiplicity of nodes at each level.
            Complexity is:
                0 IF graph has multiplicity at only the leaf level
                1 IF graph has multiplicity at multiple levels AND ALL PATHS END AT TAXON LEVEL
                2 IF graph has multiplicity at multiple levels AND THERE ARE MULTIPLE TERMINAL TAXON LEVELS


        Args:
            graph (dict): A dictionary representation of the graph, where:
                        key = Source node
                        value = Target node
                        (Nodes are represented as described elsewhere)

        Returns:
            complexity (int): Integer description of graph complexity.
        """
        terminal_levels = set()
        terminal_nodes = []
        for k, v in graph.items():
            if len(v) == 0:
                terminal_nodes.append(k)
                terminal_levels.add(k[1])
        if len(terminal_levels) > 1:
            complexity = 2
        else:
            breadth_indicator_level = (TaxonLevel(list(terminal_levels)[0]) - 1).lvl
            indicator_nodes = [node for node in graph.keys() if node[1] == breadth_indicator_level]
            if len(indicator_nodes) > 1:
                complexity = 1
            else:
                complexity = 0

        return complexity

    def find_all_paths(self, graph: dict, source: tuple, target: tuple, path: list =[]):
        """Generic function to find all paths from given source node to given target node.

        Args:
            graph (dict): A dictionary representation of the graph, where:
                            key = Source node
                            value = Target node
            source (hashable object): Source node from which to start
            target (hashable object): Target node to which paths should lead
            path (list, optional): A known path from given source to given target [Default = []]

        Returns:
            list(lists): A list of paths (where each path is a list)
        """
        ## initialise path by adding source
        path = path + [source]

        ## handle edge cases
        if source == target:
            return [path]
        if source not in graph:
            return []

        ## paths is a list in case of multiple paths from source to target
        paths = []

        ## iterate over graph dict (i.e. edges)
        for node in graph[source]:
            if node not in path:
                ## recursively find all "sub-paths" and update paths list
                newpaths = self.find_all_paths(graph, node, target, path)
                for newpath in newpaths:
                    paths.append(newpath)

        return paths

    def split_tree(self, graph: dict, split_at: str):
        """Split a given graph at specified level. For example, splitting the following graph at "S1":
                                    (0, S)
                                    |
                                    |
                            (1, S1)--------(2, S1)   <----- Splitting at this level
                                |             |
                                |             |
                        (3, S2)---(4, S2)  (5, S2)

            Gives:
                    (0, S)                  (0, S)
                    |                       |
                    |                       |
                    (1, S1)         AND     (2, S1)
                    |                       |
                    |                       |
            (3, S2)---(4, S2)            (5, S2)


        Args:
            graph (dict): A dictionary representation of the graph, where:
                        key = Source node
                        value = Target node
                        (Nodes are represented as described elsewhere)
            split_at (str): Taxon level to split graph at

        Returns:
            subgraphs (list): List of subgraphs found after splitting intput graph at given level
        """
        subgraphs = []
        root = sorted(list(graph.keys()))[0]

        split_points = [k for k in graph.keys() if k[1] == split_at]
        if len(split_points) <= 1:
            logging.debug(msg=f"Can't split at {split_at}")
            return graph

        for split_pt in split_points:
            path_to_split_pt = self.find_all_paths(graph, root, split_pt)[0]
            if len(graph[split_pt]) != 0:
                post_split_leaves = graph[split_pt]
                path_to_split_pt.extend(post_split_leaves)
                for node in post_split_leaves:
                    if len(graph[node]) != 0:
                        path_to_split_pt.extend(graph[node])
            final_path = sorted(path_to_split_pt)

            subgraphs.append(self.build_from(final_path))

        return subgraphs

    def decompose_tree(self, graph: dict, known_trees: list = []):
        """Split a complex graph into simple graphs.
            Here, a simple graph is defined as one where there is only one node at all levels except the lowest one.

        Args:
            graph (dict): A dictionary representation of the graph, where:
                        key = Source node
                        value = Target node
                        (Nodes are represented as described elsewhere)
            known_subgraphs (list, optional): List of known subgraphs to output. Defaults to [].

        Returns:
            subgraphs (list(dicts)): List of simple subgraphs decomposed from input graph
        """
        ## collect known subgraphs
        subgraphs = known_trees

        ## get a list of taxon levels and find the highest level with count > 1
        levels_list = [k[1] for k in graph.keys()]
        for i in sorted(levels_list):
            if levels_list.count(i) > 1:
                multiplicity_level = i
                break

        ## dump graphs split at level identified above into temporary list
        ## for each subgraph check if it is simple, if not, call this function on it again
        tmp_subgraphs = self.split_tree(graph, multiplicity_level)

        for g in tmp_subgraphs:
            if self.get_tree_complexity(g) == 0:
                subgraphs.append(g)
            else:
                logging.debug(f"Decomposing {g.keys()} further")
                ## recursion is cool!
                ## (but only when it works)
                cmp_subgraphs = self.decompose_tree(g, subgraphs)

        return subgraphs

