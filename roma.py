import networkx as nx
import matplotlib.pyplot as plt

class RomaPuzzle:
    def __init__(self, n, roma_cell, arrows, boxes):
        """
        Initialize a Roma Puzzle instance.
        
        Parameters:
        - n: Size of the grid (n x n)
        - roma_cell: Tuple (i, j) for the position of the Roma-cell ◦
        - arrows: Dictionary with keys as (i, j) tuples and values as directions ('↑', '→', '↓', '←')
        - boxes: Dictionary mapping (i, j) tuples to box identifiers
        """
        self.n = n  # Grid size n x n
        self.grid = [[' ' for _ in range(n)] for _ in range(n)]  # Empty grid
        self.boxes = boxes  # Box assignments
        self.roma_cell = roma_cell # Roma cell
        self.grid[roma_cell[0]][roma_cell[1]] = '◦' # Place the Roma-cell in the grid
        for (i, j), arrow in arrows.items(): # Place arrows in the grid
            self.grid[i][j] = arrow

        # Board condition
        # 1st -> Board size (n * n) = number of cells
        size = len([cell for cell in boxes.keys()])
        if n * n != size:
            raise ValueError("Board shape n * n != number of cells")

        # Boxes condition
        # Initialize condition
        # 1st -> Check boxes condition boxes are made up of adjacent cells
        boxes_inversed = {}
        for i in range(self.n):
            for j in range(self.n):
                if (i, j) in boxes:
                    box_id = boxes[(i, j)]
                    if box_id not in boxes_inversed:
                        boxes_inversed[box_id] = [(i, j)]
                    else:
                        boxes_inversed[box_id].append((i, j))
        for box_id, cells in boxes_inversed.items():
            if not self._are_cells_adjacent(cells):
                raise ValueError(f"Cells in box {box_id} are not adjacent.")
        self.boxes_inversed = boxes_inversed

        # 2nd -> Check boxes condition maximum size 4
        for box_id, box in boxes_inversed.items():
            if len(box) > 4:
                raise ValueError(f"Number of cells in box {box_id} exceeds 4.")

        # 3rd -> Check that roma cell is in her own box
        for box_id, box in boxes_inversed.items():
            if roma_cell in box and len(box) > 1:
                raise ValueError(f"Roma cell must be in her own box")

        # 4th -> Check that a cell appears in one box
        cells = []
        for box in boxes_inversed.values():
            cells += box
        for cell in cells:
            cpt = 0
            for box in boxes_inversed.values():
                if cell in box:
                    cpt += 1
            if cpt != 1:
                raise ValueError(f"{cell} does not appear in one box or appear in more than 2 boxes.")

        # Condition to be checked
        # 5th -> Check that an arrow doesn't appear more than one time in boxes
        if not self.box_condition():
            raise ValueError("The set of default arrows is invalid")
        
        # Graph condition
        self.adj_matrix = self.adjacency_matrix()
        self.graph = self.generate_graph()

        # 1st -> Check that graph is acyclic
        if not nx.is_directed_acyclic_graph(self.graph):
            raise ValueError("The directed graph representing Roma puzzle contains a cycle.")
        
        # 2nd -> Roma cell has an out-degree of output of 0
        if self.graph.out_degree(self.roma_cell) != 0:
            raise ValueError("Roma cell must have an out-degree of 0.")

    def display_board(self):
        """Display the grid using matplotlib with arrows and boxes highlighted."""
        _, ax = plt.subplots()

        # Set grid limits
        ax.set_xlim(0, self.n)
        ax.set_ylim(0, self.n)

        # Draw grid lines
        for i in range(self.n + 1):
            ax.plot([i, i], [0, self.n], color='black', lw=1)
            ax.plot([0, self.n], [i, i], color='black', lw=1)

        # Highlight box borders with thicker lines
        for i in range(self.n):
            for j in range(self.n):
                current_box = self.boxes[(i, j)]

                # Check the right border (between boxes)
                if j < self.n - 1 and self.boxes[(i, j + 1)] != current_box:
                    ax.plot([j + 1, j + 1], [self.n - i - 1, self.n - i], color='black', lw=3)

                # Check the bottom border (between boxes)
                if i < self.n - 1 and self.boxes[(i + 1, j)] != current_box:
                    ax.plot([j, j + 1], [self.n - i - 1, self.n - i - 1], color='black', lw=3)

        # Add arrows and Roma-cell to the grid
        for i in range(self.n):
            for j in range(self.n):
                symbol = self.grid[i][j]
                if symbol != ' ':
                    ax.text(j + 0.5, self.n - i - 0.5, symbol, fontsize=20, ha='center', va='center')

        # Remove ticks and labels
        ax.set_xticks([])
        ax.set_yticks([])

        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def box_condition(self):
        """Check the box condition: No box contains more than one arrow of the same direction."""
        box_arrows = {}
        for i in range(self.n):
            for j in range(self.n):
                box_id = self.boxes[(i, j)]
                if box_id not in box_arrows:
                    box_arrows[box_id] = set()
                if self.grid[i][j] in {'↑', '→', '↓', '←'}:
                    if self.grid[i][j] in box_arrows[box_id]:
                        return False
                    box_arrows[box_id].add(self.grid[i][j])
        return True
    
    def _are_cells_adjacent(self, cells):
        """
        Check if a set of cells are adjacent. Cells are adjacent if they are next to each other
        horizontally or vertically (not diagonally).
        """
        if len(cells) == 1:
            return True
        else:
            for (i, j) in cells:
                has_neighbor = False
                for (x, y) in cells:
                    if (i == x and abs(j - y) == 1) or (j == y and abs(i - x) == 1):
                        has_neighbor = True
                        break
                if not has_neighbor:
                    return False
            return True
        
    def adjacency_matrix(self):
        """Generate an adjacency matrix from the Roma board."""
        # Initialisation of the adjacency matrix (n*n x n*n) (directed graph)
        matrix_size = self.n * self.n
        adjacency_matrix = [[0 for _ in range(matrix_size)] for _ in range(matrix_size)]
        
        # Helper to convert position (i, j) into an index in the adjacency matrix
        def position_to_index(i, j):
            return i * self.n + j

        # Go through each cell on the board
        for i in range(self.n):
            for j in range(self.n):
                current_cell = self.grid[i][j]
                current_index = position_to_index(i, j)

                # Check directions and add edges accordingly
                if current_cell == '→' and j + 1 < self.n:
                    right_index = position_to_index(i, j + 1)
                    adjacency_matrix[current_index][right_index] = 1
                elif current_cell == '←' and j - 1 >= 0:
                    left_index = position_to_index(i, j - 1)
                    adjacency_matrix[current_index][left_index] = 1
                elif current_cell == '↓' and i + 1 < self.n:
                    down_index = position_to_index(i + 1, j)
                    adjacency_matrix[current_index][down_index] = 1
                elif current_cell == '↑' and i - 1 >= 0:
                    up_index = position_to_index(i - 1, j)
                    adjacency_matrix[current_index][up_index] = 1

        return adjacency_matrix
    
    def generate_graph(self):
        """Generates the directed graph representing the incidence matrix """
        # Construct the graph from the adjacency matrix
        adj_matrix = self.adjacency_matrix()
        G = nx.DiGraph()  # Oriented graph

        # Ajouter les nœuds
        for i in range(self.n):
            for j in range(self.n):
                G.add_node((i, j))  # Add a node for each cell

        # Add directed edges according to the adjacency matrix
        for i in range(self.n * self.n):
            for j in range(self.n * self.n):
                if adj_matrix[i][j] == 1:
                    # Convert matrix indices into (i, j) coordinates
                    node_i = (i // self.n, i % self.n)
                    node_j = (j // self.n, j % self.n)
                    G.add_edge(node_i, node_j)
        
        return G

    def display_graph(self):
        """Display the graph constructed from the adjacency matrix."""
        # Positioning nodes on a grid
        pos = {(i, j): (j, -i) for i in range(self.n) for j in range(self.n)}

        # Drawing the graph with NetworkX
        plt.figure(figsize=(8, 8))
        nx.draw(self.graph, pos, with_labels=True, node_size=500, node_color="skyblue", font_size=10, font_weight="bold", arrows=True, arrowstyle='->', arrowsize=20)
        plt.show()
    
    def update_arrows(self, new_arrows):
        """
        Update the arrows on the board with new directions and validate the board after modification.

        Parameters:
        - new_arrows: Dictionary with keys as (i, j) tuples and values as directions ('↑', '→', '↓', '←')
        """
        # Step 1: Update the arrows on the grid
        for (i, j), new_arrow in new_arrows.items():
            self.grid[i][j] = new_arrow

        # Step 2: Check that arrows respect the box condition
        if not self.box_condition():
            raise ValueError("The added arrows do not fit the box condition (duplicate arrows in a box).")
        
        # Step 3: Recalculate the adjacency matrix and graph
        self.adj_matrix = self.adjacency_matrix()  # Update the adjacency matrix
        self.graph = self.generate_graph()  # Generate the new graph
        
        # Step 4: Check if the graph is acyclic
        if not nx.is_directed_acyclic_graph(self.graph):
            raise ValueError("The directed graph representing Roma puzzle contains a cycle.")
        
        # Step 5: Check if the Roma cell has an out-degree of 0
        if self.graph.out_degree(self.roma_cell) != 0:
            raise ValueError("Roma cell must have an out-degree of 0.")