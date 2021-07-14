import PyEFVLib
import sys
sys.path.append("../../../LocalModules")
import pickle
import json
import numpy as np
from time import perf_counter
import matplotlib.pylab as plt
from scipy.sparse import bmat, coo_matrix
from BiotModelEquations import *
from ControllersHandler import TimeHandler
from PreProcessing import *


def main():
	preproc = PreProcessing("settings.json")

if __name__ == '__main__':
	main()
