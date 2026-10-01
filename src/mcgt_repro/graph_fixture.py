"""Twelve-vertex grid, its oriented boundary, and the four declared cycles."""

from fractions import Fraction

from mcgt_repro.exact_rank import nullspace, rank_rows


EDGES = [
    (0, 1), (1, 2), (2, 3),
    (4, 5), (5, 6), (6, 7),
    (8, 9), (9, 10), (10, 11),
    (0, 4), (1, 5), (2, 6), (3, 7),
    (4, 8), (5, 9), (6, 10), (7, 11),
]

TREE_EDGES = [
    (0, 4), (1, 5), (2, 6), (3, 7),
    (4, 8), (5, 9), (6, 10), (7, 11),
    (0, 1), (1, 2), (2, 3),
]

DECLARED_WALKS = [
    [(0, 1), (1, 5), (5, 4), (4, 0)],
    [(1, 2), (2, 6), (6, 5), (5, 1)],
    [(4, 5), (5, 9), (9, 8), (8, 4)],
    [(5, 6), (6, 10), (10, 9), (9, 5)],
]

N_VERTICES = 12


def edge_index():
    return {pair: index for index, pair in enumerate(EDGES)}


def tree_indices():
    lookup = edge_index()
    return [lookup[pair] for pair in TREE_EDGES]


def _signed_row(walk):
    lookup = edge_index()
    row = [Fraction(0)] * len(EDGES)
    for tail, head in walk:
        if (tail, head) in lookup:
            row[lookup[(tail, head)]] += Fraction(1)
        elif (head, tail) in lookup:
            row[lookup[(head, tail)]] += Fraction(-1)
        else:
            raise ValueError("walk step is not a free edge")
    return row


def declared_cycle_rows():
    return [_signed_row(walk) for walk in DECLARED_WALKS]


def boundary_matrix():
    """Rows of the oriented boundary D: edge coordinates to vertex coordinates."""
    rows = []
    for vertex in range(N_VERTICES):
        row = [Fraction(0)] * len(EDGES)
        for index, (tail, head) in enumerate(EDGES):
            if head == vertex:
                row[index] += Fraction(1)
            if tail == vertex:
                row[index] -= Fraction(1)
        rows.append(row)
    return rows


def component_count():
    parent = list(range(N_VERTICES))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for tail, head in EDGES:
        left, right = find(tail), find(head)
        if left != right:
            parent[left] = right
    return len({find(vertex) for vertex in range(N_VERTICES)})


def tree_is_spanning_forest():
    parent = list(range(N_VERTICES))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for tail, head in TREE_EDGES:
        left, right = find(tail), find(head)
        if left == right:
            return False
        parent[left] = right
    components = len({find(vertex) for vertex in range(N_VERTICES)})
    return components == 1 and len(TREE_EDGES) == N_VERTICES - 1


def cycle_space_basis():
    """A basis of ker D. Each vector has one coordinate per free edge."""
    return nullspace(boundary_matrix(), len(EDGES))


def rows_in_kernel(rows, operator_rows):
    images_are_zero = []
    for row in rows:
        image = [
            sum(left * right for left, right in zip(operator_row, row))
            for operator_row in operator_rows
        ]
        images_are_zero.append(all(entry == 0 for entry in image))
    return all(images_are_zero)


def cycle_counts():
    boundary = boundary_matrix()
    basis = cycle_space_basis()
    declared = declared_cycle_rows()
    return {
        "component_count": component_count(),
        "tree_is_spanning_forest": tree_is_spanning_forest(),
        "boundary_rank": rank_rows(boundary),
        "q_max": rank_rows(basis),
        "q": rank_rows(declared),
        "declared_rows_lie_in_kernel": rows_in_kernel(declared, boundary),
        "cycle_basis_lies_in_kernel": rows_in_kernel(basis, boundary),
    }
