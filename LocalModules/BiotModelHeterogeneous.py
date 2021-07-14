import sys
sys.path.append("/home/herminio/git_remote_repos/PyEFVLib")
import PyEFVLib
import numpy as np
from SparseMatrix import SparseMatrixCOO
from scipy.sparse import csc_matrix, csr_matrix
from BiotModelEquations import *

#	| K      		 L        ||Uᵏ⁺¹|   |0   0||U⁰|   |0      0   ||Uᵏ|   | bu  |
#	|        		          ||    | = |     ||  | + |           ||  | + |     |
#	| 0  (A + (1/δ)*B + Δt*H) ||Pᵏ⁺¹|   |Q   A||P⁰|   |Q   (1/δ)*B||Pᵏ|   |Δt*bp|

def assemble_K_L_Q_H_heterogeneous2(problemData, props):
	propertyData = problemData.propertyData
	grid = problemData.grid
	n 	 = grid.numberOfVertices
	d	 = grid.dimension
	nDOF = 4

	ls = SparseMatrixCOO(n, nDOF)

	for region in grid.regions:
		cf   = props[region.name]["cf"]
		mu   = props[region.name]["mu"]
		rhof = props[region.name]["rhof"]
		for element in region.elements:
			nu   = props[region.name]["nu"]["value"][element.handle]
			G    = props[region.name]["G"]["value"][element.handle]
			cs   = props[region.name]["cs"]["value"][element.handle]
			phi  = props[region.name]["phi"]["value"][element.handle]
			kx   = props[region.name]["kx"]["value"][element.handle]
			ky   = props[region.name]["ky"]["value"][element.handle]
			kz   = props[region.name]["kz"]["value"][element.handle]
			rhos = props[region.name]["rhos"]["value"][element.handle]

			k = np.array([[kx, 0, 0], [0, ky, 0], [0, 0, kz]])

			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([	[2*G+lame,	lame,		lame,		0,	0,	0],
							[lame,		2*G+lame,	lame,	 	0, 	0, 	0],
							[lame, 		lame, 		2*G+lame,	0,	0, 	0],
							[0, 		0, 			0, 			G, 	0, 	0],
							[0, 		0, 			0, 			0, 	G, 	0],
							[0, 		0, 			0, 			0, 	0, 	G]])
			rho = phi*rhof + (1 - phi)*rhos
			K = 2*G*(1 + nu)/(3*(1-2*nu))
			cb = 1/K
			alpha = 1 - cs/cb
			if alpha < phi:
				raise Exception("Inconsistent properties: alpha = %f < phi = %f"%(alpha,phi))
			S = phi*cf + (alpha - phi)*cs

			for face in element.innerFaces:
				m = element.vertices.size
				area = face.area.getCoordinates()[:d]
				transposedVoigtArea = get_transposed_Voigt_area(face)
				shapeFunctions = face.getShapeFunctions()
				backwardHandle, forwardHandle = face.getNeighborVerticesHandles()

				# Matrix H -> -k * mu * grad(p)
				matrixCoefficients = -mu*np.matmul(area.T, np.matmul(k, face.globalDerivatives))
				for local, vertex in enumerate(element.vertices):
					ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(d)*n, 1)
					ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(d)*n, 1)

				# Matrix K -> Ce * grad_s(u)
				voigtGradientOperator = get_Voigt_gradient_operator(face.globalDerivatives)
				matrixCoefficients = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, Ce, voigtGradientOperator)
				for i in range(d):
					for j in range(d):
						for local, vertex in enumerate(element.vertices):
							ls.addValueToMatrix(backwardHandle+(i+0)*n, vertex.handle+(j+0)*n, 1)
							ls.addValueToMatrix(forwardHandle+(i+0)*n, vertex.handle+(j+0)*n, 1)

				identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)])
				matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)
				for coord in range(d):
					for local, vertex in enumerate(element.vertices):
						# Matrix L -> -alpha * grad(p)
						ls.addValueToMatrix(backwardHandle+(coord+0)*n, vertex.handle+(d)*n, 1)
						ls.addValueToMatrix(forwardHandle+(coord+0)*n, vertex.handle+(d)*n, 1)

						# Matrix Q -> alpha * u
						ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(coord+0)*n, 1)
						ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(coord+0)*n, 1)

	# Matrix Q -> alpha * u
	for facet in grid.facets:
		for outerFace in facet.outerFaces:
			area = outerFace.area.getCoordinates()[:d]
			shapeFunctions = outerFace.getShapeFunctions()
			for coord in range(d):
				for local, vertex in enumerate(facet.vertices):
					ls.addValueToMatrix(outerFace.vertex.handle+(d)*n, vertex.handle+(coord+0)*n, 1)
	ls.build()
	return ls


def assemble_K_L_Q_H_heterogeneous(problemData, props):
	propertyData = problemData.propertyData
	grid = problemData.grid
	n 	 = grid.numberOfVertices
	d	 = grid.dimension
	nDOF = 4

	ls = SparseMatrixCOO(n, nDOF)

	for region in grid.regions:
		cf   = props[region.name]["cf"]
		mu   = props[region.name]["mu"]
		rhof = props[region.name]["rhof"]
		for element in region.elements:
			nu   = props[region.name]["nu"]["value"][element.handle]
			G    = props[region.name]["G"]["value"][element.handle]
			cs   = props[region.name]["cs"]["value"][element.handle]
			phi  = props[region.name]["phi"]["value"][element.handle]
			kx   = props[region.name]["kx"]["value"][element.handle]
			ky   = props[region.name]["ky"]["value"][element.handle]
			kz   = props[region.name]["kz"]["value"][element.handle]
			rhos = props[region.name]["rhos"]["value"][element.handle]

			k = np.array([[kx, 0, 0], [0, ky, 0], [0, 0, kz]])

			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([	[2*G+lame,	lame,		lame,		0,	0,	0],
							[lame,		2*G+lame,	lame,	 	0, 	0, 	0],
							[lame, 		lame, 		2*G+lame,	0,	0, 	0],
							[0, 		0, 			0, 			G, 	0, 	0],
							[0, 		0, 			0, 			0, 	G, 	0],
							[0, 		0, 			0, 			0, 	0, 	G]])
			rho = phi*rhof + (1 - phi)*rhos
			K = 2*G*(1 + nu)/(3*(1-2*nu))
			cb = 1/K
			alpha = 1 - cs/cb
			if alpha < phi:
				raise Exception("Inconsistent properties: alpha = %f < phi = %f"%(alpha,phi))
			S = phi*cf + (alpha - phi)*cs

			for face in element.innerFaces:
				m = element.vertices.size
				area = face.area.getCoordinates()[:d]
				transposedVoigtArea = get_transposed_Voigt_area(face)
				shapeFunctions = face.getShapeFunctions()
				backwardHandle, forwardHandle = face.getNeighborVerticesHandles()

				# Matrix H -> -k * mu * grad(p)
				matrixCoefficients = -(1/mu)*np.matmul(area.T, np.matmul(k, face.globalDerivatives))
				for local, vertex in enumerate(element.vertices):
					ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(d)*n, matrixCoefficients[local])
					ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(d)*n, -matrixCoefficients[local])

				# Matrix K -> Ce * grad_s(u)
				voigtGradientOperator = get_Voigt_gradient_operator(face.globalDerivatives)
				matrixCoefficients = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, Ce, voigtGradientOperator)
				for i in range(d):
					for j in range(d):
						for local, vertex in enumerate(element.vertices):
							ls.addValueToMatrix(backwardHandle+(i+0)*n, vertex.handle+(j+0)*n, matrixCoefficients[i][j][local])
							ls.addValueToMatrix(forwardHandle+(i+0)*n, vertex.handle+(j+0)*n, -matrixCoefficients[i][j][local])

				identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)])
				matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)
				for coord in range(d):
					for local, vertex in enumerate(element.vertices):
						# Matrix L -> -alpha * grad(p)
						ls.addValueToMatrix(backwardHandle+(coord+0)*n, vertex.handle+(d)*n, matrixCoefficients[coord][local])
						ls.addValueToMatrix(forwardHandle+(coord+0)*n, vertex.handle+(d)*n, -matrixCoefficients[coord][local])

						# Matrix Q -> alpha * u
						ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(coord+0)*n, alpha*shapeFunctions[local]*area[coord])
						ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(coord+0)*n, -alpha*shapeFunctions[local]*area[coord])

	# # Matrix Q -> alpha * u
	for facet in grid.facets:
		elem = facet.element
		for outerFace in facet.outerFaces:
			area = outerFace.area.getCoordinates()[:d]
			shapeFunctions = outerFace.getShapeFunctions()
			for coord in range(d):
				for local, vertex in enumerate(elem.vertices):
					ls.addValueToMatrix(outerFace.vertex.handle+(d)*n, vertex.handle+(coord+0)*n, alpha * shapeFunctions[local] * area[coord])

	ls.build()
	return ls


def assemble_A_heterogeneous(problemData, props):
	propertyData = problemData.propertyData
	grid = problemData.grid
	ls_A = SparseMatrixCOO(grid.numberOfVertices, 1)
	for region in grid.regions:
		cf   = props[region.name]["cf"]
		rhof = props[region.name]["rhof"]
		for element in region.elements:
			nu   = props[region.name]["nu"]["value"][element.handle]
			G    = props[region.name]["G"]["value"][element.handle]
			cs   = props[region.name]["cs"]["value"][element.handle]
			phi  = props[region.name]["phi"]["value"][element.handle]
			rhos = props[region.name]["rhos"]["value"][element.handle]
			K = 2*G*(1 + nu)/3*(1-2*nu)
			alpha = 1 - cs*K
			S = phi*cf + (alpha - phi)*cs
			for vertex in element.vertices:
				ls_A.addValueToMatrix(vertex.handle, vertex.handle, vertex.getSubElementVolume(element) * S)
	ls_A.build()
	return ls_A

def assemble_B_heterogeneous(problemData, props):
	propertyData = problemData.propertyData
	grid = problemData.grid
	ls_B = SparseMatrixCOO(grid.numberOfVertices, 1)
	for region in grid.regions:
		cf   = props[region.name]["cf"]
		rhof = props[region.name]["rhof"]
		for element in region.elements:
			nu = props[region.name]["nu"]["value"][element.handle]
			G  = props[region.name]["G"]["value"][element.handle]
			cs = props[region.name]["cs"]["value"][element.handle]
			K  = 2*G*(1 + nu)/3*(1-2*nu)
			alpha = 1 - cs*K
			for vertex in element.vertices:
				ls_B.addValueToMatrix(vertex.handle, vertex.handle, vertex.getSubElementVolume(element)*alpha**2/K)
	ls_B.build()
	return ls_B

def assemble_bu_bp_heterogeneous(problemData, props):
	propertyData = problemData.propertyData
	grid = problemData.grid
	n 	 = grid.numberOfVertices
	d	 = grid.dimension
	nDOF = 4
	b = np.zeros((1+d)*n)

	g = np.array([0.0, 0.0, 0.0])[:d]
	for region in grid.regions:
		cf   = props[region.name]["cf"]
		mu   = props[region.name]["mu"]
		rhof = props[region.name]["rhof"]
		for element in region.elements:
			nu   = props[region.name]["nu"]["value"][element.handle]
			G    = props[region.name]["G"]["value"][element.handle]
			cs   = props[region.name]["cs"]["value"][element.handle]
			phi  = props[region.name]["phi"]["value"][element.handle]
			kx   = props[region.name]["kx"]["value"][element.handle]
			ky   = props[region.name]["ky"]["value"][element.handle]
			kz   = props[region.name]["kz"]["value"][element.handle]
			rhos = props[region.name]["rhos"]["value"][element.handle]

			k = np.array([[kx, 0, 0], [0, ky, 0], [0, 0, kz]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			alpha = 1 - cs*K
			S = (phi * cf + (alpha-phi) * cs)
			for vertex in element.vertices:
				# (-1) * rho * g
				for coord in range(d):
					b[vertex.handle+coord*n] += vertex.getSubElementVolume(element) * (-1) * rho * g[coord]

			for innerFace in element.innerFaces:
				area = innerFace.area.getCoordinates()[:d]
				shapeFunctions = innerFace.getShapeFunctions()
				backwardHandle, forwardHandle = innerFace.getNeighborVerticesHandles()

				# (-1) * k * mu * rho * g
				coefficient = -mu*np.matmul(area.T, np.matmul(k,rho*g))
				b[backwardHandle+(d)*n] += coefficient
				b[forwardHandle+(d)*n]  -= coefficient
	return b
