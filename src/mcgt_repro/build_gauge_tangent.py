"""Infinitesimal gauge action at the normalized flat reference.

A generator A_v acts by h_v = -(A_v* + A_v) and x_e = A_head - A_tail.
"""

from mcgt_repro.real_matrix_basis import (
    add,
    adjoint,
    gl_basis,
    gl_coords,
    herm_coords,
    scale,
    sub,
    zeros,
)


def gauge_tangent_columns(edges, sector, n_vertices):
    n_edges = len(edges)
    metric_width = sector * sector
    transport_width = 2 * metric_width
    expected = n_vertices * metric_width + n_edges * transport_width
    columns = []
    for vertex in range(n_vertices):
        for generator in gl_basis(sector):
            metrics = [zeros(sector) for _ in range(n_vertices)]
            metrics[vertex] = scale(add(adjoint(generator), generator), -1)
            transports = [zeros(sector) for _ in range(n_edges)]
            for index, (tail, head) in enumerate(edges):
                increment = zeros(sector)
                if head == vertex:
                    increment = add(increment, generator)
                if tail == vertex:
                    increment = sub(increment, generator)
                transports[index] = increment
            coordinates = []
            for metric in metrics:
                coordinates.extend(herm_coords(metric))
            for transport in transports:
                coordinates.extend(gl_coords(transport))
            if len(coordinates) != expected:
                raise RuntimeError("gauge coordinate length mismatch")
            columns.append(coordinates)
    return columns
