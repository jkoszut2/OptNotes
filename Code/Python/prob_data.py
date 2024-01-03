# Problem Data

import numpy as np


# PROBLEM SETUP
A = np.array([[1, 1, 0]])
b = np.array([[5]])
c = np.array([[0],
			  [0],
			  [1]])
# G = np.array([[ 0, -1,  0],
#               [ 0,  0, -1],
#               [-1,  0,  0],
#               [ 0, -1,  0]])
G = np.array([[ 0, -1,  0],
              [ 0,  0, -1.1],
              [-1.2,  0,  0],
              [ 0, -1.3,  0]])
h = np.array([[-2],
			  [ 0],
			  [ 0],
			  [ 0]])

cone_orders = np.array([1, 3]) # order of each cone
cone_exprs = np.array([1, 1]) # number of expressions in each cone
cone_vars = cone_orders*cone_exprs # number of variables in each cone
# Examples:
# 2*x1 + 2*x2 >= 0 is of order 1 and has 1 expression
# x1>=0, x2>=0 is of order 1 and has 2 expressions
# ||x||_2 >=||[y,z]^T||_2 is of order 3 and has 1 expression

maxit = 20 # maximum number of iterations
