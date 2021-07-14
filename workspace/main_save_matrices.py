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

def create_element_props(prop):
	global grid
	delta = prop["delta"]
	prop_value = prop["value"]
	np.random.seed( prop["seed"] )
	poissons = np.random.uniform(prop_value - delta, prop_value + delta, size=grid.elements.size)
	prop["value"] = poissons


def main():
	case = "Terzaghi"
	settings = readJsonData("./Settings/" + case + "/settings.json")


	props = settings["Properties"]

	meshFilePaths = [settings["Grid"]["Paths"][i] for i in settings["Grid"]["Indices"]]
	# print(meshFilePaths)

	for meshFilePath in meshFilePaths:
		outputFilePath = "../temp/" + case + meshFilePath[meshFilePath.rfind("/"):-4]
		print("outputFilePath")
		print(outputFilePath)

		tic = perf_counter()
		problemData = PyEFVLib.ProblemData(
			meshFilePath = meshFilePath,
			outputFilePath = outputFilePath,
			propertyData = PyEFVLib.PropertyData(settings["Properties"]),
			boundaryConditions = PyEFVLib.BoundaryConditions(settings["BoundaryConditions"]),
		)
		toc = perf_counter()
		print("ProblemData: %.5f"%(toc-tic))

		print(problemData.libraryPath)

		# ------------------ BUILD GRID -----------------------
		global grid
		grid = problemData.grid
		n 	 = grid.numberOfVertices
		d	 = grid.dimension
		# -----------------------------------------------------

		# ------------------ PROPERTIES -----------------------
		for region, prop in product(grid.regions, ["nu", "G", "cs", "phi", "kx", "ky", "kz", "rhos"]):
			create_element_props(settings["Properties"][region.name][prop])
		print(settings["Properties"])
		# -----------------------------------------------------

		# ------------ ASSEMBLE LINEAR SYSTEM -----------------
		tic = perf_counter()

		ls = assemble_K_L_Q_H_heterogeneous(problemData, props)
		apply_dirichlet_displacement_matrix(problemData, ls, diagonal_value=1.0)
		apply_dirichlet_pressure_matrix(problemData, ls, diagonal_value=0.0, shift=3)
		mat_K, mat_L, mat_Q, mat_H = split_full_matrix(ls.matrix, n)


		ls_A = assemble_A_heterogeneous(problemData, props)
		apply_dirichlet_pressure_matrix(problemData, ls_A, diagonal_value=1.0)
		mat_A = ls_A.matrix

		ls_B = assemble_B_heterogeneous(problemData, props)
		apply_dirichlet_pressure_matrix(problemData, ls_B, diagonal_value=0.0)
		mat_B = ls_B.matrix

		vec_b = assemble_bu_bp_heterogeneous(problemData, props)
		apply_neumann_displacements_vector(problemData, vec_b)
		apply_dirichlet_displacements_vector(problemData, vec_b)
		apply_neumann_pressure_vector(problemData, vec_b, shift=3)
		apply_dirichlet_pressure_vector(problemData, vec_b, shift=3)

		toc = perf_counter()
		print("Assembly time: %.5f"%(toc-tic))
		# -----------------------------------------------------


if __name__ == '__main__':
	main()
