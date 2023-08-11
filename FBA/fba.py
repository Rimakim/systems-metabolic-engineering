import cobra
from cobra.io import read_sbml_model
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.array import constraint_matrices
import numpy as np
import scipy.sparse as sp 
from tqdm import tqdm
from gurobipy import GRB 
from gurobipy import *



class FBA:
    def __init__(self, model, objective):
        
        self.model = model
        self.objective = self.model.reactions.get_by_id(objective)
        self.M = len(self.model.reactions)
        self.N = len(self.model.metabolites)
        self.LP = Model()
        
    def make_valuable(self):
        v = self.LP.addMVar(shape = 2*self.M, name = 'v')
        return v
    
    def make_A_mat(self):
        S = create_stoichiometric_matrix(self.model)       
        S_mat = np.concatenate((S, -S), axis = 1)          # shape = N * 2M
        
        v_eye = np.eye(self.M)                             
        V_eye = np.concatenate((v_eye, -v_eye), axis = 1)  # shape = M * 2M
        
        # shape = (2N + 2M) * 2M
        A_matrix = np.concatenate((S_mat, -S_mat, V_eye, -V_eye), axis = 0)   
        
        return A_matrix
    
    def make_b_mat(self):
        rhs_S = np.array([0 for i in range(self.N)])
        rhs_U = np.array([rxn.upper_bound for rxn in self.model.reactions])
        rhs_L = np.array([rxn.lower_bound for rxn in self.model.reactions])
        
        b_matrix = np.concatenate((rhs_S, -rhs_S, rhs_U, -rhs_L), axis = 0)
        
        return b_matrix
        
    def make_c_t_mat(self):
        ob_idx = 0  
        for i, rxn in enumerate(model.reactions):
            if rxn == self.objective:
                ob_idx = i 
        
        c_t = np.array([])
        for i in range(2*self.M):
            if i == ob_idx:
                c_t = np.append(c_t, 1)

            elif i == self.M + ob_idx:
                c_t = np.append(c_t, -1)

            else:
                c_t = np.append(c_t, 0)
        
        return c_t
    
    def optimize(self):
        x = self.make_valuable()
        A = self.make_A_mat()
        b = self.make_b_mat()
        c_t = self.make_c_t_mat()
        
        self.LP.setObjective(c_t @ x, GRB.MAXIMIZE)
        self.LP.addConstr(A @ x <= b, name = 'primal_constraints')
        
        return self.LP.optimize()
    
    def dual_optimize(self):
        y = self.LP.addMVar(shape = 2*(self.N +self.M), name = 'y')
        A_t = self.make_A_mat().T
        b_t = self.make_b_mat().T
        c = self.make_c_t_mat().T
        
        self.LP.setObjective(b_t @ y, GRB.MINIMIZE)
        self.LP.addConstr(A_t @ y >= c, name = 'dual_constraints')
        
        return self.LP.optimize()
