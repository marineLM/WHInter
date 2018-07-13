import numpy as np
import scipy.sparse as sp
import csv


def SPPformat_to_csc(file):
    """
    Read data in SPP format,
    and return a csc (compressed sparse column) matrix.

    Parameters
    ----------
    file: str
        Path to the SPP file.

    Returns
    -------
    A dictionary with entries:
        X: csc matrix
        y: response vector
    """
    f = open(file)
    indices = np.array([])
    indptr = np.array([0])
    y = np.array([])
    tot = 0

    for i, row in enumerate(csv.reader(f, delimiter="\t")):
        y = np.append(y, float(row[0]))
        row = [int(s.split(":")[0]) - 1 for s in row[1:]]
        tot += len(row)
        indices = np.append(indices, row)
        indptr = np.append(indptr, tot)

    X_csc = sp.csr_matrix((np.ones(tot), indices, indptr), dtype=int)
    X_csc = sp.csc_matrix(X_csc)

    return {'X': X_csc, 'y': y}
