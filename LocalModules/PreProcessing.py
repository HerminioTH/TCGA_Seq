import PyEFVLib
import sys
sys.path.append("../LocalModules")
import pickle
import json
import numpy as np
from time import perf_counter
import matplotlib.pylab as plt
from scipy.sparse import bmat, coo_matrix
from BiotModelHeterogeneous import *
from ControllersHandler import TimeHandler
from Utils import *
from itertools import product

#	| K      		 L        ||Uᵏ⁺¹|   |0   0||U⁰|   |0      0   ||Uᵏ|   | bu  |
#	|        		          ||    | = |     ||  | + |           ||  | + |     |
#	| 0  (A + (1/δ)*B + Δt*H) ||Pᵏ⁺¹|   |Q   A||P⁰|   |Q   (1/δ)*B||Pᵏ|   |Δt*bp|

class PreProcessing(object):
	def __init__(self, settings_path, show_time=False):
		self.show_time = show_time
		self.settings = readJsonData(settings_path)
		self.buildProblemData()
		self.writeProperties()
		self.buildMatrices()
		self.buildVector()
		self.saveLinearSystem()

	def buildProblemData(self):
		tic = perf_counter()
		self.problemData = PyEFVLib.ProblemData(
			meshFilePath = self.settings["Grid"],
			outputFilePath = self.settings["Output"],
			propertyData = PyEFVLib.PropertyData(self.settings["Properties"]),
			boundaryConditions = PyEFVLib.BoundaryConditions(self.settings["BoundaryConditions"]),
		)
		self.n = self.problemData.grid.numberOfVertices
		self.d = self.problemData.grid.dimension
		toc = perf_counter()
		if self.show_time: print("ProblemData: %.5f"%(toc-tic))

	def writeProperties(self):
		tic = perf_counter()
		for region, prop_name in product(self.problemData.grid.regions, ["nu", "G", "cs", "phi", "kx", "ky", "kz", "rhos"]):
			np.random.seed(self.settings["Properties"][region.name][prop_name]["seed"])
			delta 		= self.settings["Properties"][region.name][prop_name]["delta"]
			prop_value 	= self.settings["Properties"][region.name][prop_name]["value"]
			self.settings["Properties"][region.name][prop_name]["value"] = np.random.uniform(prop_value - delta, prop_value + delta, size=self.problemData.grid.elements.size)
		toc = perf_counter()
		if self.show_time: print("writeProperties: %.5f"%(toc-tic))

	def buildMatrices(self):
		tic = perf_counter()

		ls = assemble_K_L_Q_H_heterogeneous(self.problemData, self.settings["Properties"])
		apply_dirichlet_displacement_matrix(self.problemData, ls, diagonal_value=1.0)
		apply_dirichlet_pressure_matrix(self.problemData, ls, diagonal_value=0.0, shift=3)
		self.mat_K, self.mat_L, self.mat_Q, self.mat_H = split_full_matrix(ls.matrix, self.n)

		ls_A = assemble_A_heterogeneous(self.problemData, self.settings["Properties"])
		apply_dirichlet_pressure_matrix(self.problemData, ls_A, diagonal_value=1.0)
		self.mat_A = ls_A.matrix

		apply_dirichlet_pressure_matrix(self.problemData, ls_A, diagonal_value=0.0)
		self.mat_A_rhs = ls_A.matrix.copy()

		ls_B = assemble_B_heterogeneous(self.problemData, self.settings["Properties"])
		apply_dirichlet_pressure_matrix(self.problemData, ls_B, diagonal_value=0.0)
		self.mat_B = ls_B.matrix

		toc = perf_counter()
		if self.show_time: print("buildMatrices: %.5f"%(toc-tic))

	def buildVector(self):
		tic = perf_counter()

		self.vec_b = assemble_bu_bp_heterogeneous(self.problemData, self.settings["Properties"])
		apply_neumann_displacements_vector(self.problemData, self.vec_b)
		apply_dirichlet_displacements_vector(self.problemData, self.vec_b)
		apply_neumann_pressure_vector(self.problemData, self.vec_b, shift=3)
		apply_dirichlet_pressure_vector(self.problemData, self.vec_b, shift=3)

		toc = perf_counter()
		if self.show_time: print("buildVector: %.5f"%(toc-tic))

	def saveLinearSystem(self):
		tic = perf_counter()
		savePickleData(self.settings["Matrices"], "/mat_K", 	self.mat_K)
		savePickleData(self.settings["Matrices"], "/mat_L", 	self.mat_L)
		savePickleData(self.settings["Matrices"], "/mat_Q", 	self.mat_Q)
		savePickleData(self.settings["Matrices"], "/mat_H", 	self.mat_H)
		savePickleData(self.settings["Matrices"], "/mat_B", 	self.mat_B)
		savePickleData(self.settings["Matrices"], "/mat_A", 	self.mat_A)
		savePickleData(self.settings["Matrices"], "/mat_A_rhs", self.mat_A_rhs)
		savePickleData(self.settings["Matrices"], "/vec_b", 	self.vec_b)
		toc = perf_counter()
		if self.show_time: print("saveLinearSystem: %.5f"%(toc-tic))


if __name__ == '__main__':
	preproc = PreProcessing("./Settings/Terzaghi/settings.json")
