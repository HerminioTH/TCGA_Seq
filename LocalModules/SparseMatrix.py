import numpy as np
from scipy.sparse import coo_matrix

class SparseMatrixCOO(object):
    def __init__(self, nVertices, nDOF):
        self.n = nVertices
        self.shape = (self.n*nDOF, self.n*nDOF)
        self.rows = []
        self.cols = []
        self.data = []

    def addValueToMatrix(self, row, col, value):
        self.rows.append(row)
        self.cols.append(col)
        self.data.append(value)

    def initialize(self):
        pass

    def build(self):
        self.matrix = coo_matrix((self.data, (self.rows, self.cols)), shape=self.shape)
