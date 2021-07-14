import json
import numpy as np
import pickle
import os

def computeRate(error):
    # return (np.log10(error[0]) - np.log10(max(1e-200, error[-1])))/len(error)
    return -np.log10(error[-1])/len(error)

def readJsonData(data_file):
    with open(data_file, "r") as jsonFile:
        data = json.load(jsonFile)
    return data

def saveDataToJson(data_file, data):
	with open(data_file, "w") as jsonFile:
		json.dump(data, jsonFile, indent=3)

def loadMatrix(fileName):
	mat = np.loadtxt(fileName, delimiter=" ")
	return coo_matrix((mat[:,2].astype(float), (mat[:,0].astype(int), mat[:,1].astype(int))))

def loadVector(fileName):
	return np.loadtxt(fileName)

def loadPickleData(fileName):
	with open(fileName, "rb") as f:
		data_pickle = pickle.load(f)
	return data_pickle

def savePickleData(directory, fileName, matrix):
	if not os.path.exists(directory):
		os.makedirs(directory)
	with open(directory + fileName, "wb") as f:
		pickle.dump(matrix, f)

def generateRandom(grid, props, material_name, prop_name):
	prop  = props[material_name][prop_name]["value"]
	delta = props[material_name][prop_name]["delta"]
	np.random.seed(props[material_name][prop_name]["seed"])
	props[material_name][prop_name]["value"] = np.random.uniform(prop - delta, prop + delta, size=grid.elements.size)


def computeCouplingStrength(material_props):
	G = material_props["G"]["value"]
	nu = material_props["nu"]["value"]
	cs = material_props["cs"]["value"]
	cf = material_props["cf"]
	phi = material_props["phi"]["value"]
	K = 2*G*(1 + nu)/(3 - 6*nu)
	alpha = 1 - cs*K
	if alpha < phi:
		raise Exception(f"Inconsistent properties: alpha = {alpha} < phi = {phi}")
	S = phi*cf + (alpha - phi)*cs
	return alpha*alpha/S/K

def computeConsolidationCoefficient(material_props):
	mu    = material_props["mu"]
	cf    = material_props["cf"]
	cs    = material_props["cs"]["value"]
	kx    = material_props["kx"]["value"]
	ky    = material_props["ky"]["value"]
	kz    = material_props["kz"]["value"]
	nu    = material_props["nu"]["value"]
	phi   = material_props["phi"]["value"]
	G     = material_props["G"]["value"]

	K     = 2*G*(1 + nu)/(3 - 6*nu)
	alpha = 1 - cs*K
	M 	  = 2*G*(1 - nu)/(1 - 2*nu)
	Q 	  = 1./(cf*phi + cs*(alpha - phi))
	k     = (kx + ky + kz)/3
	return k*Q*M/(mu*(M + alpha*alpha*Q))

def computeTime(material_props, h, t_dimensionless=0.5):
	c = computeConsolidationCoefficient(material_props)
	return h*h*t_dimensionless/c
