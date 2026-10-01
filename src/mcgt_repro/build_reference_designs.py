"""Coordinate rows of Designs A and B at the normalized flat reference.

Design A reads every vertex Hermitian block and the skew-Hermitian block of
each spanning-tree transport. Design B reads the Hermitian block at one root
and the full transport block of each spanning-tree edge.
"""

from fractions import Fraction


def _edge_offset(n_vertices, sector, edge_index):
    return n_vertices * (sector * sector) + edge_index * (2 * sector * sector)


def design_a_rows(n_vertices, n_edges, tree_indices, sector):
    metric_width = sector * sector
    domain = n_vertices * metric_width + n_edges * (2 * metric_width)
    rows = []
    for vertex in range(n_vertices):
        for coordinate in range(metric_width):
            row = [Fraction(0)] * domain
            row[vertex * metric_width + coordinate] = Fraction(1)
            rows.append(row)
    for edge in tree_indices:
        base = _edge_offset(n_vertices, sector, edge) + metric_width
        for coordinate in range(metric_width):
            row = [Fraction(0)] * domain
            row[base + coordinate] = Fraction(1)
            rows.append(row)
    return rows


def design_b_rows(n_vertices, n_edges, tree_indices, sector, root=0):
    metric_width = sector * sector
    transport_width = 2 * metric_width
    domain = n_vertices * metric_width + n_edges * transport_width
    rows = []
    for coordinate in range(metric_width):
        row = [Fraction(0)] * domain
        row[root * metric_width + coordinate] = Fraction(1)
        rows.append(row)
    for edge in tree_indices:
        base = _edge_offset(n_vertices, sector, edge)
        for coordinate in range(transport_width):
            row = [Fraction(0)] * domain
            row[base + coordinate] = Fraction(1)
            rows.append(row)
    return rows
