import numpy as np
from SparseMatrix import SparseMatrixCOO
from scipy.sparse import csc_matrix, csr_matrix

#	| K      L 	    ||U|   |0   0||U⁰|   | bu  |
#	|               || | = |     ||  | + |     |
#	| Q  (A + Δt*H) ||P|   |Q   A||P⁰|   |Δt*bp|


def mat_zero_row(ls, row, diagonal_value):
	if not isinstance(ls.matrix, csr_matrix):
		ls.matrix = ls.matrix.tocsr()
	pos1 = ls.matrix.indptr[row]
	pos2 = ls.matrix.indptr[row+1]
	ls.matrix.data[ls.matrix.indptr[row]:ls.matrix.indptr[row+1]] = 0.0
	index = pos1 + np.where(ls.matrix.indices[pos1:pos2] == row)[0][0]
	ls.matrix.data[index] = diagonal_value


def mat_zero_sym(ls, row, diagonal_value):
	def apply(mat, value, matrix_format):
		if matrix_format == "columns":
			if not isinstance(mat, csc_matrix):
				mat = mat.tocsc()
		elif matrix_format == "rows":
			if not isinstance(mat, csr_matrix):
				mat = mat.tocsr()
		else:
			raise NameError("Choose between rows or columns")
		pos1 = mat.indptr[row]
		pos2 = mat.indptr[row+1]
		mat.data[mat.indptr[row]:mat.indptr[row+1]] = 0.0
		index = pos1 + np.where(mat.indices[pos1:pos2] == row)[0][0]
		mat.data[index] = value
		return mat
	ls.matrix = apply(ls.matrix, diagonal_value, "columns")
	ls.matrix = apply(ls.matrix, diagonal_value, "rows")



def get_transposed_Voigt_area(innerFace):
	Sx, Sy, Sz = innerFace.area.getCoordinates()
	return np.array([	[Sx,	0,		0,		Sy,		0,		Sz 	],
						[0,		Sy,		0,		Sx,		Sz,		0	],
						[0,		0,		Sz,		0,		Sy,		Sx 	]])

def get_Voigt_gradient_operator(globalDerivatives):
	Nx,Ny,Nz = globalDerivatives
	zero=np.zeros(Nx.size)
	return np.array([	[Nx,	zero,	zero],
						[zero,	Ny,		zero],
						[zero,	zero,	Nz	],
						[Ny,	Nx,		zero],
						[zero,	Nz,		Ny 	],
						[Nz,	zero,	Nx	]])

def split_full_matrix(mat_FULL, n):
	'''
		mat_FULL - coo_matrix
		n - number of vertices
	'''
	aux_FULL = mat_FULL.tocsr()
	upperBlock = aux_FULL[:3*n,:].tocsc()
	lowerBlock = aux_FULL[3*n:,:].tocsc()

	Block11 = upperBlock[:,:3*n]
	Block12 = upperBlock[:,3*n:]
	Block21 = lowerBlock[:,:3*n]
	Block22 = lowerBlock[:,3*n:]

	return Block11, Block12, Block21, Block22

def assemble_K_L_Q_H(problemData):
	propertyData = problemData.propertyData
	grid = problemData.grid
	n 	 = grid.numberOfVertices
	d	 = grid.dimension
	nDOF = 4

	ls = SparseMatrixCOO(n, nDOF)

	for region in grid.regions:
		nu   = propertyData.get(region.handle, "nu")
		G    = propertyData.get(region.handle, "G")
		cs   = propertyData.get(region.handle, "cs")
		phi  = propertyData.get(region.handle, "phi")
		k    = propertyData.get(region.handle, "kx")
		cf   = propertyData.get(region.handle, "cf")
		mu   = propertyData.get(region.handle, "mu")
		rhos = propertyData.get(region.handle, "rhos")
		rhof = propertyData.get(region.handle, "rhof")

		g = np.array([0.0, 0.0, 0.0])[:d]
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
		S = phi*cf + (alpha - phi)*cs

		for element in region.elements:
			for face in element.innerFaces:
				m = element.vertices.size
				area = face.area.getCoordinates()[:d]
				transposedVoigtArea = get_transposed_Voigt_area(face)
				shapeFunctions = face.getShapeFunctions()
				backwardHandle, forwardHandle = face.getNeighborVerticesHandles()

				# (-1) * k * mu * grad(p)
				matrixCoefficients = -k*mu*np.matmul(area.T, face.globalDerivatives)
				for local, vertex in enumerate(element.vertices):
					ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(d)*n, matrixCoefficients[local])
					ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(d)*n, -matrixCoefficients[local])

				# Ce * grad_s(u)
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
						# (-1) * alpha * grad(p)
						ls.addValueToMatrix(backwardHandle+(coord+0)*n, vertex.handle+(d)*n, matrixCoefficients[coord][local])
						ls.addValueToMatrix(forwardHandle+(coord+0)*n, vertex.handle+(d)*n, -matrixCoefficients[coord][local])

						# alpha * u
						ls.addValueToMatrix(backwardHandle+(d)*n, vertex.handle+(coord+0)*n, alpha*shapeFunctions[local]*area[coord])
						ls.addValueToMatrix(forwardHandle+(d)*n, vertex.handle+(coord+0)*n, -alpha*shapeFunctions[local]*area[coord])

	# alpha * u
	for facet in grid.facets:
		for outerFace in facet.outerFaces:
			area = outerFace.area.getCoordinates()[:d]
			shapeFunctions = outerFace.getShapeFunctions()
			for coord in range(d):
				for local, vertex in enumerate(facet.vertices):
					ls.addValueToMatrix(outerFace.vertex.handle+(d)*n, vertex.handle+(coord+0)*n, alpha * shapeFunctions[local] * area[coord])
	ls.build()
	return ls




def assemble_A(problemData):
	propertyData = problemData.propertyData
	grid = problemData.grid
	ls_A = SparseMatrixCOO(grid.numberOfVertices, 1)
	for region in grid.regions:
		nu  = propertyData.get(region.handle, "nu")
		G   = propertyData.get(region.handle, "G")
		cs  = propertyData.get(region.handle, "cs")
		phi = propertyData.get(region.handle, "phi")
		cf  = propertyData.get(region.handle, "cf")
		K = 2*G*(1 + nu)/3*(1-2*nu)
		alpha = 1 - cs*K
		S = phi*cf + (alpha - phi)*cs
		for element in region.elements:
			for vertex in element.vertices:
				ls_A.addValueToMatrix(vertex.handle, vertex.handle, vertex.getSubElementVolume(element) * S)
	ls_A.build()
	return ls_A

def assemble_B(problemData):
	propertyData = problemData.propertyData
	grid = problemData.grid
	ls_B = SparseMatrixCOO(grid.numberOfVertices, 1)
	for region in grid.regions:
		nu  = propertyData.get(region.handle, "nu")
		G   = propertyData.get(region.handle, "G")
		cs  = propertyData.get(region.handle, "cs")
		K = 2*G*(1 + nu)/3*(1-2*nu)
		alpha = 1 - cs*K
		for element in region.elements:
			for vertex in element.vertices:
				# ls_B.addValueToMatrix(vertex.handle, vertex.handle, 0.0)
				ls_B.addValueToMatrix(vertex.handle, vertex.handle, vertex.getSubElementVolume(element)*alpha**2/K)
	ls_B.build()
	return ls_B


def assemble_bu_bp(problemData):
	propertyData = problemData.propertyData
	grid = problemData.grid
	n 	 = grid.numberOfVertices
	d	 = grid.dimension
	nDOF = 4
	b = np.zeros((1+d)*n)

	for region in grid.regions:
		nu   = propertyData.get(region.handle, "nu")
		G    = propertyData.get(region.handle, "G")
		cs   = propertyData.get(region.handle, "cs")
		phi  = propertyData.get(region.handle, "phi")
		k    = propertyData.get(region.handle, "k")
		cf   = propertyData.get(region.handle, "cf")
		mu   = propertyData.get(region.handle, "mu")
		rhos = propertyData.get(region.handle, "rhos")
		rhof = propertyData.get(region.handle, "rhof")

		g = np.array([0.0, 0.0, 0.0])[:d]
		rho = phi * rhof + (1-phi) * rhos
		K = 2*G*(1 + nu) / 3*(1-2*nu)
		cb = 1 / K
		alpha = 1 - cs / cb
		S = (phi * cf + (alpha-phi) * cs)

		for element in region.elements:
			for vertex in element.vertices:
				# (-1) * rho * g
				for coord in range(d):
					b[vertex.handle+coord*n] += vertex.getSubElementVolume(element) * (-1) * rho * g[coord]

			for innerFace in element.innerFaces:
				area = innerFace.area.getCoordinates()[:d]
				shapeFunctions = innerFace.getShapeFunctions()
				backwardHandle, forwardHandle = innerFace.getNeighborVerticesHandles()

				# (-1) * k * mu * rho * g
				coefficient = np.matmul(area.T, (-1)*k*mu*rho*g)
				b[backwardHandle+(d)*n] += coefficient
				b[forwardHandle+(d)*n]  -= coefficient
	return b





def apply_dirichlet_pressure_matrix(problemData, ls, diagonal_value=1.0, shift=0):
	n = problemData.grid.numberOfVertices
	if not isinstance(ls.matrix, csr_matrix):
		ls.matrix = ls.matrix.tocsr()
	for bCondition in problemData.dirichletBoundaries["p"]:
		for vertex in bCondition.boundary.vertices:
			# mat_zero_row(ls, vertex.handle + shift*n, diagonal_value)
			mat_zero_sym(ls, vertex.handle + shift*n, diagonal_value)

def apply_dirichlet_displacement_matrix(problemData, ls, diagonal_value=1.0):
	n = problemData.grid.numberOfVertices
	if not isinstance(ls.matrix, csr_matrix):
		ls.matrix = ls.matrix.tocsr()

	for bCondition in problemData.dirichletBoundaries["u_x"]:
		for vertex in bCondition.boundary.vertices:
			mat_zero_row(ls, vertex.handle+(0)*n, diagonal_value)
			# mat_zero_sym(ls, vertex.handle+(0)*n, diagonal_value)

	for bCondition in problemData.dirichletBoundaries["u_y"]:
		for vertex in bCondition.boundary.vertices:
			mat_zero_row(ls, vertex.handle+(1)*n, diagonal_value)
			# mat_zero_sym(ls, vertex.handle+(1)*n, diagonal_value)

	for bCondition in problemData.dirichletBoundaries["u_z"]:
		for vertex in bCondition.boundary.vertices:
			mat_zero_row(ls, vertex.handle+(2)*n, diagonal_value)
			# mat_zero_sym(ls, vertex.handle+(2)*n, diagonal_value)

def apply_neumann_displacements_vector(problemData, vector_b):
	n = problemData.grid.numberOfVertices
	for bCondition in problemData.neumannBoundaries["u_x"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				vector_b[outerFace.vertex.handle+(0)*n] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	for bCondition in problemData.neumannBoundaries["u_y"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				vector_b[outerFace.vertex.handle+(1)*n] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	for bCondition in problemData.neumannBoundaries["u_z"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				vector_b[outerFace.vertex.handle+(2)*n] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

def apply_dirichlet_displacements_vector(problemData, vector_b):
	n = problemData.grid.numberOfVertices
	for bCondition in problemData.dirichletBoundaries["u_x"]:
		for vertex in bCondition.boundary.vertices:
			vector_b[vertex.handle+(0)*n] = bCondition.getValue(vertex.handle)

	for bCondition in problemData.dirichletBoundaries["u_y"]:
		for vertex in bCondition.boundary.vertices:
			vector_b[vertex.handle+(1)*n] = bCondition.getValue(vertex.handle)

	for bCondition in problemData.dirichletBoundaries["u_z"]:
		for vertex in bCondition.boundary.vertices:
			vector_b[vertex.handle+(2)*n] = bCondition.getValue(vertex.handle)


def apply_neumann_pressure_vector(problemData, vector_b, shift=0):
	n = problemData.grid.numberOfVertices
	for bCondition in problemData.neumannBoundaries["p"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				vector_b[outerFace.vertex.handle+shift*n] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

def apply_dirichlet_pressure_vector(problemData, vector_b, shift=0):
	n = problemData.grid.numberOfVertices
	for bCondition in problemData.dirichletBoundaries["p"]:
		for vertex in bCondition.boundary.vertices:
			vector_b[vertex.handle+shift*n] = bCondition.getValue(vertex.handle)


def updateIndependent(problemData, vec_b, mat_Q, mat_A, oldUField, oldPField, timeStep):
	n = problemData.grid.numberOfVertices
	d = problemData.grid.dimension
	independent = vec_b.copy()
	independent[3*n:] *= timeStep
	independent[3*n:] += mat_Q*oldUField + mat_A*oldPField
	apply_dirichlet_pressure_vector(problemData, independent, shift=3)
	return independent

def buildIndependent(vec_b, mat_Q, mat_A_rhs, oldUField, oldPField, timeStep):
	'''
		Be careful, only Dirichlet equal to ZERO or Neumann equal do ZERO is considered for
		pressure.
	'''
	n = mat_A_rhs.shape[0]
	independent = vec_b.copy()
	independent[3*n:] = timeStep*independent[3*n:] + mat_Q*oldUField + mat_A_rhs*oldPField
	return independent
