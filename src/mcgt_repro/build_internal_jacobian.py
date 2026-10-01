"""Internal Jacobian at the normalized flat reference g_v = I, L_e = I."""

from mcgt_repro.real_matrix_basis import (
    add,
    adjoint,
    gl_basis,
    gl_coords,
    herm_basis,
    herm_coords,
    scale,
    sub,
    zeros,
)


def domain_dimension(n_vertices, n_edges, sector):
    return n_vertices * (sector * sector) + n_edges * (2 * sector * sector)


def _metric_block(h_source, h_target, x_edge):
    """Hermitian coordinates of d(L* g_t L - g_s) at the identity."""
    residual = add(adjoint(x_edge), h_target)
    residual = add(residual, x_edge)
    residual = sub(residual, h_source)
    return herm_coords(residual)


def internal_jacobian_rows(edges, declared_rows, sector, n_vertices):
    """Rows of J_int. Columns run over the real basis of the state space."""
    n_edges = len(edges)
    metric_basis = herm_basis(sector)
    transport_basis = gl_basis(sector)
    columns = []

    def pack(metrics, transports):
        coordinates = []
        for index, (source, target) in enumerate(edges):
            coordinates.extend(
                _metric_block(metrics[source], metrics[target], transports[index])
            )
        for signed in declared_rows:
            total = zeros(sector)
            for index, coefficient in enumerate(signed):
                if coefficient != 0:
                    total = add(total, scale(transports[index], coefficient))
            coordinates.extend(gl_coords(total))
        return coordinates

    for vertex in range(n_vertices):
        for basis_matrix in metric_basis:
            metrics = [zeros(sector) for _ in range(n_vertices)]
            metrics[vertex] = basis_matrix
            transports = [zeros(sector) for _ in range(n_edges)]
            columns.append(pack(metrics, transports))
    for edge in range(n_edges):
        for basis_matrix in transport_basis:
            metrics = [zeros(sector) for _ in range(n_vertices)]
            transports = [zeros(sector) for _ in range(n_edges)]
            transports[edge] = basis_matrix
            columns.append(pack(metrics, transports))

    rows = [list(row) for row in zip(*columns)]
    return rows, domain_dimension(n_vertices, n_edges, sector)
