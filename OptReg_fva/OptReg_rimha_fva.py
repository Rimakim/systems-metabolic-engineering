import os
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.flux_analysis import add_loopless

import numpy as np
import pandas as pd
from gurobipy import *
from tqdm import tqdm
from utils import argument_parser




if __name__=='__main__':
    
    parser = argument_parser()
    options = parser.parse_args()

    model_path = options.model_path
    strain = options.strain
    output_path = options.output_path
    biomass_lb = options.biomass_lb
    target_rxn = options.target_rxn
    max_manipulation = options.max_manipulation
    regulation_strength = options.regulation_strength
    eps = options.eps
    oxygen_condition = options.oxygen_condition
    fraction_of_optimum = options.fraction_of_optimum
    loopless = options.loopless
    
    iter_num = options.iter_num
    cpu_num = options.cpu_num
    default_bound = options.default_bound





    
    class OptReg:
    
        def __init__(self, model_path, strain, biomass_lb, target_rxn, L, C, eps, oxygen_condition, fraction_of_optimum, iter_num):
    
            self.model = read_sbml_model(model_path)
            self.strain = strain
            self.target_rxn = target_rxn
            self.biomass_lb = biomass_lb
            self.L = L
            self.C = C
            self.eps = eps
            self.oxygen_condition = oxygen_condition
            self.f_opt = fraction_of_optimum
            self.iter_num = iter_num
    
            # biomass, target
            self.biomass_rxn = None
    
            for rxn in self.model.reactions:
                if rxn.objective_coefficient == 1:
                    self.biomass_rxn = rxn.id

            # oxygen condition
            if self.oxygen_condition  == 'anaerobic':
                self.model.reactions.EX_o2_e.lower_bound = 0
            else:
                pass

            
            # Possible execretion condition
            '''
            Originated from Kamp et al. Nature Comm. (2017)
            
            <E.coli(iML1515)>
            ethanol, lactate, formate, succinate, hydrogen, methanol
            
            <S.cerevisiae(iMM904)>
            ethanol, glycerol, pyruvate, acetate, succinate
            
            <A.niger>
            gluconate, citrate, oxalate, malate, succinate, erythritol
            
            <C.glutamicum>
            glutamate, succinate, lysine, lactate, acetate, alanine, isoleucine, glycine
            '''

            if self.strain == 'e_coli':
                
                pos_ex_list = [self.target_rxn, 'EX_ac_e', 'EX_co2_e', 'EX_etoh_e', 
                               'EX_for_e', 'EX_h_e', 'EX_lac__L_e', 
                               'EX_succ_e', 'EX_meoh_e', 'EX_h2o_e'
                              ]

                # execretion flux = 0 if metabolite is organic chemocals 
                for rxn in self.model.reactions:
                    
                    if 'EX_' not in rxn.id:
                        continue
                    if rxn.upper_bound <= 0:
                        continue
    
                    met = self.model.metabolites.get_by_id(rxn.id[3:])
                    formula = met.formula
                    
                    if 'C' in formula and 'H' in formula and 'O' in formula: # organic chemicals
                        rxn.upper_bound = 0
    
                
                for rxn_id in pos_ex_list:
                    self.model.reactions.get_by_id(rxn_id).upper_bound = 1000

            else:
                pass

        def essential_reactions(self):

            essential_list = []
            
            with self.model as m:
                for rxn in m.reactions:
                    m.reactions.get_by_id(rxn.id).knock_out()
                    sol = m.slim_optimize()

                    if sol < self.biomass_lb:
                        essential_list.append(rxn.id)

            return essential_list
                    

        def calculate_cellular_state(self, loopless):
    
            self.rxn_ub = {}
            self.rxn_lb = {}
            self.vo_u = {}
            self.vo_l = {}

            self.loopless = loopless
    
            # calculation cellular max, min with FBA
            print('finding cellular flux max, min, steady-state value...')

            
            if self.loopless == 'T':
                add_loopless(self.model)
            else:
                pass
                

            for rxn in tqdm(self.model.reactions):
        
                with self.model as m:
                    m.objective = rxn.id

                    m.objective_direction = 'max'
                    ub = m.slim_optimize()
                    self.rxn_ub[rxn.id] = ub
                    
                    m.objective_direction = 'min'
                    lb = m.slim_optimize()
                    self.rxn_lb[rxn.id] = lb
                    
            
            # calculation steady-state range with FVA
            
            with self.model as m:
                
                m.objective = self.biomass_rxn

                if self.loopless == 'T':
                    fva_sol = fva(m, fraction_of_optimum=self.f_opt, loopless=True)
                else:
                    fva_sol = fva(m, fraction_of_optimum=self.f_opt, loopless=False)

            
            for rxn in self.model.reactions:
                self.vo_u[rxn.id] = fva_sol.loc[rxn.id, 'maximum']
                self.vo_l[rxn.id] = fva_sol.loc[rxn.id, 'minimum']

            
            # reversible/irreversible reactions
            
            self.for_list = []       # only forward (flux > 0)
            self.rev_list = []       # only reverse (flux < 0)
            self.for_rev_list = []   # reversible
    
            for rxn in self.model.reactions:
                if self.rxn_ub[rxn.id] > 0 and self.rxn_lb[rxn.id] >= 0:
                    self.for_list.append(rxn.id)
            
                elif self.rxn_ub[rxn.id] <= 0 and self.rxn_lb[rxn.id] < 0:
                    self.rev_list.append(rxn.id)
            
                else:
                    self.for_rev_list.append(rxn.id) 


            return self.rxn_ub, self.rxn_lb, self.vo_u, self.vo_l
    
    
        def model_setting(self, cpu_num, default_bound):
        
            self.cpu_num = cpu_num
            self.ub_default = default_bound
            self.lb_default = -1*default_bound
    
            self.m = Model()
            self.m.reset()
            self.m.params.Threads = cpu_num
    
            self.m.params.FeasibilityTol = 1e-9
            self.m.params.OptimalityTol = 1e-9
            self.m.params.BarConvTol = 1e-12
            self.m.params.IntFeasTol = 1e-9
            
            self.m.update()
    
    
        def variable_setting(self):
            
            self.v = {}
    
            '''
            <primal variable addition>
            
            only forward: rxn name_for
            only reverse: rxn name_rev
            reversible: both are included
            '''
                
            for rxn in self.model.reactions:
            
                if rxn.id in self.for_list:
                    self.v[f'{rxn.id}_for'] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} for')
            
                elif rxn.id in self.rev_list:
                    self.v[f'{rxn.id}_rev'] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} rev')
            
                else:
                    self.v[f'{rxn.id}_for'] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} for')
                    self.v[f'{rxn.id}_rev'] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} rev')
            
            self.m.update()
    
    
            '''
            <max, min values for primary(flux) variables>
            all flux values >= 0
            '''
            
            self.v_ub = {}  # upper bounds
            self.v_lb = {}  # lower bounds
            self.var_o_u = {} # cellular-state up
            self.var_o_l = {} # cellular-state low
            
            for var in self.v.keys():

                if '_for' in var:
                    self.v_lb[var] = max(0, self.rxn_lb[var[:-4]]) 
                    self.var_o_l[var] = max(0, self.vo_l[var[:-4]])  
                    self.var_o_u[var] = max(0, self.vo_u[var[:-4]])
                    self.v_ub[var] = max(0, self.rxn_ub[var[:-4]])

                else:
                    self.v_lb[var] = abs(min(0, self.rxn_ub[var[:-4]])) 
                    self.var_o_l[var] = abs(min(0, self.vo_u[var[:-4]]))  
                    self.var_o_u[var] = abs(min(0, self.vo_l[var[:-4]]))
                    self.v_ub[var] = abs(min(0, self.rxn_lb[var[:-4]]))     
               
            
            # dual, binary variables
            self.ds = {}
            self.qku = {}
            self.qkl = {}
            self.quu = {}
            self.qul = {}
            self.qdu = {}
            self.qdl = {}
            self.d_bio_atp = {}
            
            self.yk = {}
            self.yu = {}
            self.yd = {}

            self.zku = {}
            self.zkl = {}
            self.zul = {}
            self.zdu = {}

            # for reversible reactions
            self.yuf = {} 
            self.yur = {}
            self.ydf = {}
            self.ydr = {}
            
            for rxn in self.v.keys():

                self.qku[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO ub dual')
                self.qkl[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO lb dual')
                
                self.quu[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} up ub dual')
                self.qul[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} up lb dual')
                    
                self.qdu[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} down ub dual')
                self.qdl[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} down lb dual')

                self.zku[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO ub z')
                self.zkl[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO lb z')
                self.zul[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} up lb z')
                self.zdu[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} down ub z')

                
                if rxn[:-4] in self.for_rev_list:

                    if '_for' in rxn:
                        self.yuf[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')
                        self.ydf[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')
                    else:
                        self.yur[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')
                        self.ydr[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')

                else:
                    self.yu[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')     # on = 0, off = 1
                    self.yd[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')     # on = 0, off = 1

                
                if rxn[:-4] == self.biomass_rxn or rxn[:-4] == 'ATPM':
                    self.d_bio_atp[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} min dual')

            
            # binary variables for KO
            for rxn in self.for_list:
                self.yk[f'{rxn}_for'] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} KO on/off')     # on = 0, off = 1     
            for rxn in self.rev_list:
                self.yk[f'{rxn}_rev'] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} KO on/off')     # on = 0, off = 1
            for rxn in self.for_rev_list:
                self.yk[f'{rxn}'] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} KO on/off')     # on = 0, off = 1         

            
            for met in self.model.metabolites:
                self.ds[met.id] = self.m.addVar(lb = self.lb_default, ub = self.ub_default, name = f'{met.id} stoich')
    
            self.m.update()
    
        
        def make_S_matrix(self):
            # S_matrix (forward, reverse reaction 반영)
            # S*v_forward - S*v_reverse = 0
            S_matrix_v = {}

            
            for var in tqdm(self.v.keys()):
                
                rxn = self.model.reactions.get_by_id(var[:-4])
               
                if '_for' in var:   
                    
                    reactant_list = list(rxn.reactants)
                    product_list = list(rxn.products)
                    
                    reactant_coff_list = list(rxn.get_coefficients(reactant_list))
                    product_coff_list = list(rxn.get_coefficients(product_list))
                    
                    met_list = reactant_list + product_list
                    coff_list = reactant_coff_list + product_coff_list
                    
                    for met, coff in zip(met_list, coff_list):
                        S_matrix_v[(var, met.id)] = coff
                
                else:
                    reactant_list = list(rxn.products)
                    product_list = list(rxn.reactants)
                    
                    reactant_coff_list = [-1*coff for coff in rxn.get_coefficients(reactant_list)]
                    product_coff_list = [-1*coff for coff in rxn.get_coefficients(product_list)]
                                             
                    met_list = reactant_list + product_list
                    coff_list = reactant_coff_list + product_coff_list
                    
                                         
                    for met, coff in zip(met_list, coff_list):
                        S_matrix_v[(var, met.id)] = coff
                        
            
            pair_v, self.coff_value_v = multidict(S_matrix_v)
            self.pair_v = tuplelist(pair_v)
    
        
        def add_constraints(self):
    
            '''
            <primal constraints>
            1. Stoichiometric constraints 
            2. KO constraints
            3. UP_regulation constraints
            4. DOWN_regulation constraints
            5. minimal biomass, ATPM constraints
            '''

    
            
            # 1. Stoichiometric constraints
            
            for met in tqdm(self.model.metabolites):

                self.m.addConstr(quicksum( self.v[val_id]*self.coff_value_v[val_id, met_id] for val_id, met_id in self.pair_v.select('*', met.id) ) == 0,
                           name = f'{met.id} mass balance v')
            
            self.m.update()
    
    
            # KO, UP, DOWN constraints
            for var in self.v.keys():
    
                # 2. KO constraints (reversible reactions are Knockouted together)
                if var[:-4] in self.for_rev_list:
                    self.m.addConstr(self.v[var] >= self.v_lb[var] * self.yk[var[:-4]])
                    self.m.addConstr(self.v[var] <= self.v_ub[var] * self.yk[var[:-4]])
                else:
                    self.m.addConstr(self.v[var] >= self.v_lb[var] * self.yk[var])
                    self.m.addConstr(self.v[var] <= self.v_ub[var] * self.yk[var])

            
            for var in self.v.keys():

                if var[:-4] in self.for_rev_list:
                    if '_for' in var:
                        # up_f
                        self.m.addConstr(self.v[var] >= ( self.var_o_u[var]*(1-self.C) + self.v_ub[var]*self.C ) * ( 1-self.yuf[var] ) + ( self.v_lb[var] * self.yuf[var] ))
                        self.m.addConstr(self.v[var] <= self.v_ub[var])
                    
                        # dw_f
                        self.m.addConstr(self.v[var] >= self.v_lb[var])
                        self.m.addConstr(self.v[var] <= ( self.var_o_l[var]*(1-self.C) + self.v_lb[var]*self.C ) * ( 1-self.ydf[var] ) + ( self.v_ub[var] * self.ydf[var] ))

                        # reverse reaction switch
                        self.m.addConstr(self.v[f'{var[:-4]}_rev'] >= self.v_lb[f'{var[:-4]}_rev'] * self.yuf[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_rev'] <= self.v_ub[f'{var[:-4]}_rev'] * self.yuf[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_rev'] >= self.v_lb[f'{var[:-4]}_rev'] * self.ydf[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_rev'] <= self.v_ub[f'{var[:-4]}_rev'] * self.ydf[var])
                        
                    else:
                        # up_r
                        self.m.addConstr(self.v[var] >= ( self.var_o_u[var]*(1-self.C) + self.v_ub[var]*self.C ) * ( 1-self.yur[var] ) + ( self.v_lb[var] * self.yur[var] ))
                        self.m.addConstr(self.v[var] <= self.v_ub[var])

                        # dw_r
                        self.m.addConstr(self.v[var] >= self.v_lb[var])
                        self.m.addConstr(self.v[var] <= ( self.var_o_l[var]*(1-self.C) + self.v_lb[var]*self.C ) * ( 1-self.ydr[var] ) + ( self.v_ub[var] * self.ydr[var] ))

                        # reverse reaction switch
                        self.m.addConstr(self.v[f'{var[:-4]}_for'] >= self.v_lb[f'{var[:-4]}_for'] * self.yur[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_for'] <= self.v_ub[f'{var[:-4]}_for'] * self.yur[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_for'] >= self.v_lb[f'{var[:-4]}_for'] * self.ydr[var])
                        self.m.addConstr(self.v[f'{var[:-4]}_for'] <= self.v_ub[f'{var[:-4]}_for'] * self.ydr[var])
                

                
                else:
                    # 3. Up regulation constraints
                    self.m.addConstr(self.v[var] >= ( self.var_o_u[var]*(1-self.C) + self.v_ub[var]*self.C ) * ( 1-self.yu[var] ) + ( self.v_lb[var] * self.yu[var] ))
                    self.m.addConstr(self.v[var] <= self.v_ub[var])
                
                    # 4. Down regulation constraints
                    self.m.addConstr(self.v[var] >= self.v_lb[var])
                    self.m.addConstr(self.v[var] <= ( self.var_o_l[var]*(1-self.C) + self.v_lb[var]*self.C ) * ( 1-self.yd[var] ) + ( self.v_ub[var] * self.yd[var] ))
            
            self.m.update()
    
            
            # 5. minimal biomass, atp_maint
            self.m.addConstr(self.v[f'{self.biomass_rxn}_for'] >= self.biomass_lb)
            self.m.addConstr(self.v['ATPM_for'] >= self.v_lb['ATPM_for'])
            
            self.m.update()
    
    
            '''
            <Dual constraints> 
            '''
            for var in tqdm(self.v.keys()):
                
                if var[:-4] == self.biomass_rxn:
                    self.m.addConstr(quicksum(self.coff_value_v[var_id, met_id] * self.ds[met_id] for var_id, met_id in self.pair_v.select(var, '*')) 
                                + self.qku[var] - self.qkl[var] + self.quu[var] - self.qul[var] + self.qdu[var] - self.qdl[var] - self.d_bio_atp[var] >= 1 - self.eps)
                
                elif var[:-4] == 'ATPM':
                    self.m.addConstr(quicksum(self.coff_value_v[val_id, met_id] * self.ds[met_id] for val_id, met_id in self.pair_v.select(var, '*')) 
                                + self.qku[var] - self.qkl[var] + self.quu[var] - self.qul[var] + self.qdu[var] - self.qdl[var] - self.d_bio_atp[var] >= -1*self.eps)
                
                else:
                    self.m.addConstr(quicksum(self.coff_value_v[val_id, met_id] * self.ds[met_id] for val_id, met_id in self.pair_v.select(var, '*')) 
                                + self.qku[var] - self.qkl[var] + self.quu[var] - self.qul[var] + self.qdu[var] - self.qdl[var]  >= -1*self.eps)
                
                
            self.m.update()    
    
    
            '''
            <constraints for z variables for linearization (z = q*y)>
            '''
            for var in self.v.keys():
    
                if var[:-4] in self.for_rev_list:
                    self.m.addConstr(self.zku[var] >= 0)
                    self.m.addConstr(self.zku[var] <= self.ub_default * self.yk[var[:-4]])
                    self.m.addConstr(self.zku[var] >= self.qku[var] - self.ub_default*(1 - self.yk[var[:-4]]))
                    self.m.addConstr(self.zku[var] <= self.qku[var])
                    
                    self.m.addConstr(self.zkl[var] >= 0)
                    self.m.addConstr(self.zkl[var] <= self.ub_default * self.yk[var[:-4]])
                    self.m.addConstr(self.zkl[var] >= self.qkl[var] - self.ub_default*(1 - self.yk[var[:-4]]))
                    self.m.addConstr(self.zkl[var] <= self.qkl[var])

                    if '_for' in var:
                        self.m.addConstr(self.zul[var] >= 0)
                        self.m.addConstr(self.zul[var] <= self.ub_default * self.yuf[var])
                        self.m.addConstr(self.zul[var] >= self.qul[var] - self.ub_default*(1 - self.yuf[var]))
                        self.m.addConstr(self.zul[var] <= self.qul[var])
                        
                        self.m.addConstr(self.zdu[var] >= 0)
                        self.m.addConstr(self.zdu[var] <= self.ub_default * self.ydf[var])
                        self.m.addConstr(self.zdu[var] >= self.qdu[var] - self.ub_default*(1 - self.ydf[var]))
                        self.m.addConstr(self.zdu[var] <= self.qdu[var])

                    else:
                        self.m.addConstr(self.zul[var] >= 0)
                        self.m.addConstr(self.zul[var] <= self.ub_default * self.yur[var])
                        self.m.addConstr(self.zul[var] >= self.qul[var] - self.ub_default*(1 - self.yur[var]))
                        self.m.addConstr(self.zul[var] <= self.qul[var])
                        
                        self.m.addConstr(self.zdu[var] >= 0)
                        self.m.addConstr(self.zdu[var] <= self.ub_default * self.ydr[var])
                        self.m.addConstr(self.zdu[var] >= self.qdu[var] - self.ub_default*(1 - self.ydr[var]))
                        self.m.addConstr(self.zdu[var] <= self.qdu[var])
                        
                else:
                    self.m.addConstr(self.zku[var] >= 0)
                    self.m.addConstr(self.zku[var] <= self.ub_default * self.yk[var])
                    self.m.addConstr(self.zku[var] >= self.qku[var] - self.ub_default*(1 - self.yk[var]))
                    self.m.addConstr(self.zku[var] <= self.qku[var])
                    
                    self.m.addConstr(self.zkl[var] >= 0)
                    self.m.addConstr(self.zkl[var] <= self.ub_default * self.yk[var])
                    self.m.addConstr(self.zkl[var] >= self.qkl[var] - self.ub_default*(1 - self.yk[var]))
                    self.m.addConstr(self.zkl[var] <= self.qkl[var])
                
                    self.m.addConstr(self.zul[var] >= 0)
                    self.m.addConstr(self.zul[var] <= self.ub_default * self.yu[var])
                    self.m.addConstr(self.zul[var] >= self.qul[var] - self.ub_default*(1 - self.yu[var]))
                    self.m.addConstr(self.zul[var] <= self.qul[var])
                    
                    self.m.addConstr(self.zdu[var] >= 0)
                    self.m.addConstr(self.zdu[var] <= self.ub_default * self.yd[var])
                    self.m.addConstr(self.zdu[var] >= self.qdu[var] - self.ub_default*(1 - self.yd[var]))
                    self.m.addConstr(self.zdu[var] <= self.qdu[var])

            self.m.update()

            
            '''
            <Constraints for strong duality>
            primal obj function(max) = dual obj function(min)
            '''
            
            self.m.addConstr(
                self.v[f'{self.biomass_rxn}_for'] - quicksum(self.v[var] for var in self.v.keys())*self.eps
                == -1*( ( self.biomass_lb * self.d_bio_atp[f'{self.biomass_rxn}_for'] ) + ( self.v_lb['ATPM_for'] * self.d_bio_atp['ATPM_for'] ) )
            + quicksum( ( (self.zku[var]*self.v_ub[var]) - (self.zkl[var]*self.v_lb[var]) ) for var in self.v.keys() )
            + quicksum( self.quu[var]*self.v_ub[var] - self.zul[var]*self.v_lb[var] + ( self.var_o_u[var]*(1-self.C) + self.v_ub[var]*self.C )*( self.zul[var] - self.qul[var] ) for var in self.v.keys() ) 
            + quicksum( self.zdu[var]*self.v_ub[var] - self.qdl[var]*self.v_lb[var] + ( self.var_o_l[var]*(1-self.C) + self.v_lb[var]*self.C )*( self.qdu[var] - self.zdu[var] ) for var in self.v.keys() ) 
            )

            
            self.m.update()
    
    
            '''
            <Constraints for binary variables>
            
            1. Each reaction can have only one(or zero) regulation state.
            2. Maximum regulations
            3. In reversible reactions, forward and reverse reactions cannot be up or down-regulated at the same time.
            4. Reactions that are NOT linked with GPR cannot be manipulated. (yk, yu, yd values are always 1)
            5. Do not KO Essential genes: make biomass < min_biomass when Knockouted 
            '''
    
            # 1. Each reaction can have only one(or zero) regulation state.

            for var in self.v.keys():
                if var[:-4] in self.for_rev_list:
                    var = var[:-4]
                    self.m.addConstr((1-self.yk[var]) + (1-self.yuf[f'{var}_for']) + (1-self.ydf[f'{var}_for']) + (1-self.yur[f'{var}_rev']) + (1-self.ydr[f'{var}_rev']) <= 1.1)
                else:
                    self.m.addConstr((1-self.yk[var]) + (1-self.yu[var]) + (1-self.yd[var]) <= 1.1)
                    
            self.m.update()
            
            # 2. Maximum regulations
                
            rev = quicksum((1-self.yk[rxn]) + (1-self.yuf[f'{rxn}_for']) + (1-self.ydf[f'{rxn}_for'])
                           + (1-self.yur[f'{rxn}_rev']) + (1-self.ydr[f'{rxn}_rev']) for rxn in self.for_rev_list)
            
            irrev = quicksum((1-self.yk[var]) + (1-self.yu[var]) + (1-self.yd[var]) for var in self.v.keys() if var[:-4] not in self.for_rev_list)
            
            self.m.addConstr(rev + irrev <= self.L + 0.1)
               
            self.m.update()
            
            # 4. Reactions that are NOT linked with gene_reaction_rule cannot be manipulated. (yk, yu, yd values are always 1) 
    
            for val in self.v.keys():
                gpr = self.model.reactions.get_by_id(val[:-4]).gene_reaction_rule
                
                if gpr == '' or 's' in gpr or 'ex' in val:    # not linked GPR or spontaneous reactions    
            
                    if val[:-4] in self.for_rev_list:
                        val = val[:-4]
                        self.m.addConstr(self.yk[val] >= 0.9)
                        self.m.addConstr(self.yuf[f'{val}_for'] >= 0.9)
                        self.m.addConstr(self.ydf[f'{val}_for'] >= 0.9)
                        self.m.addConstr(self.yur[f'{val}_rev'] >= 0.9)
                        self.m.addConstr(self.ydr[f'{val}_rev'] >= 0.9)
                    else:
                        self.m.addConstr(self.yk[val] >= 0.9)
                        self.m.addConstr(self.yu[val] >= 0.9)
                        self.m.addConstr(self.yd[val] >= 0.9)

            
            # 5. Do not KO Essential genes: make biomass < min_biomass when Knockouted

            essential_list = self.essential_reactions()

            for val in self.v.keys():
                if val[:-4] in essential_list:

                    if val[:-4] in self.for_rev_list:
                        self.m.addConstr(self.yk[val[:-4]] >= 0.9)
                    else:
                        self.m.addConstr(self.yk[val] >= 0.9)
                    
            
            self.m.update()
    
    
            '''
            <Objective function>
            Maximize target chemical
            '''
            if self.target_rxn in self.for_rev_list:
                self.m.setObjective(self.v[f'{self.target_rxn}_for'] - self.v[f'{self.target_rxn}_rev'], GRB.MAXIMIZE)

            else:
                self.m.setObjective(self.v[f'{self.target_rxn}_for'], GRB.MAXIMIZE)

            self.m.update()
    
        def return_results(self):
    
            manipulation_target_iter = []
            optimized_biomass_iter = []
            optimized_target_iter = []
    
            # finding multiple manipulation target sets 
            for i in range(self.iter_num):  
                
                print()
                print('iterations:', i+1)
                
                self.m.optimize()
                opt_status = self.m.status
                
                if opt_status == 2:
    
                    manipulation_target = []
                    manipulation_state = []
                    manipulation = []
    
                    for rxn_id, yk_value in self.yk.items():
                        if yk_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_KO')
                    
                    for rxn_id, yu_value in self.yu.items():
                        if yu_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_UP')
                    
                    for rxn_id, yd_value in self.yd.items():
                        if yd_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_DW')

                    for rxn_id, yuf_value in self.yuf.items():
                        if yuf_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_DW')

                    for rxn_id, ydf_value in self.ydf.items():
                        if ydf_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_DW')

                    for rxn_id, yur_value in self.yur.items():
                        if yur_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_DW')

                    for rxn_id, ydr_value in self.ydr.items():
                        if ydr_value.x < 0.1:
                            manipulation_target.append(rxn_id)
                            manipulation_state.append('_DW')
    
                    
                    for tgt, ste in zip(manipulation_target, manipulation_state):
                        manipulation.append(tgt+ste)
    
                    manipulation_target_iter.append(manipulation)
                    optimized_biomass_iter.append(self.v[f'{self.biomass_rxn}_for'].x)
                    optimized_target_iter.append(self.v[f'{self.target_rxn}_for'].x)
                    
                    # manipulated reactions in previous step
                    for rxn_id in manipulation_target:

                        if '_for' not in rxn_id and '_rev' not in rxn_id:    # KO of reversible reaction 

                            self.m.addConstr(self.yk[rxn_id] >= 0.9)
                            self.m.addConstr(self.yuf[f'{rxn_id}_for'] >= 0.9)
                            self.m.addConstr(self.ydf[f'{rxn_id}_for'] >= 0.9)
                            self.m.addConstr(self.yur[f'{rxn_id}_rev'] >= 0.9)
                            self.m.addConstr(self.ydr[f'{rxn_id}_rev'] >= 0.9)

                        else:                                                                  
                            if rxn_id[:-4] in self.for_rev_list:             # UP, DW of reversible reaction
                                rxn_id = rxn_id[:-4]
                                self.m.addConstr(self.yk[rxn_id] >= 0.9)
                                self.m.addConstr(self.yuf[f'{rxn_id}_for'] >= 0.9)
                                self.m.addConstr(self.ydf[f'{rxn_id}_for'] >= 0.9)
                                self.m.addConstr(self.yur[f'{rxn_id}_rev'] >= 0.9)
                                self.m.addConstr(self.ydr[f'{rxn_id}_rev'] >= 0.9)
                            
                            else:                                            # KO, UP, DW of irreversible reactions
                                self.m.addConstr(self.yk[rxn_id] >= 0.9)
                                self.m.addConstr(self.yu[rxn_id] >= 0.9)
                                self.m.addConstr(self.yd[rxn_id] >= 0.9)

                    self.m.update()
                        
                   
                else:
                    break
                    
            return manipulation_target_iter, optimized_biomass_iter, optimized_target_iter, opt_status


##############################################################################################################################################
    
    '''
    <Run OptReg and Output results>
    '''
    
    optreg = OptReg(model_path = model_path,
                    strain = strain,
                    biomass_lb = biomass_lb,
                    target_rxn = target_rxn,
                    L = max_manipulation,
                    C = regulation_strength,
                    eps = eps, 
                    oxygen_condition = oxygen_condition,
                    fraction_of_optimum = fraction_of_optimum,
                    iter_num = iter_num,
                   )
    
    rxn_ub, rxn_lb, vo_u, vo_l = optreg.calculate_cellular_state(loopless = loopless)
    
    optreg.model_setting(cpu_num = cpu_num, default_bound = default_bound)
    optreg.variable_setting()
    optreg.make_S_matrix()
    optreg.add_constraints()
    
    manipulation_target_iter, optimized_biomass_iter, optimized_target_iter, opt_status = optreg.return_results()
    
    
    df = pd.DataFrame()
    model = read_sbml_model(model_path)
    
    if len(manipulation_target_iter) < iter_num:
            manipulation_target_iter = manipulation_target_iter + (iter_num - len(manipulation_target_iter)) * ['-']
        
    for i in range(iter_num):
        if manipulation_target_iter[i] == '-':
            pass
            
        else:
            biomass_value = optimized_biomass_iter[i]
            objective_value = optimized_target_iter[i]
            manipulation_targets = manipulation_target_iter[i]
            manipulation_target_names = [] 
            manipulation_target_GPRs = []

            for tgt in manipulation_targets:

                if '_for' in tgt:
                    manipulation_target_names.append(model.reactions.get_by_id(tgt[:-7]).name)
                    manipulation_target_GPRs.append(model.reactions.get_by_id(tgt[:-7]).gene_reaction_rule)

                elif '_rev' in tgt:
                    manipulation_target_names.append(model.reactions.get_by_id(tgt[:-7]).name)
                    manipulation_target_GPRs.append(model.reactions.get_by_id(tgt[:-7]).gene_reaction_rule)

                else: 
                    manipulation_target_names.append(model.reactions.get_by_id(tgt[:-3]).name)
                    manipulation_target_GPRs.append(model.reactions.get_by_id(tgt[:-3]).gene_reaction_rule)
                            
            df[f'step{i+1}'] = [biomass_value, 
                                objective_value, 
                                manipulation_targets, 
                                manipulation_target_names, 
                                manipulation_target_GPRs
                               ]



    # Making result csv file
    
    df.index = ['biomass', 'objective', 'manipulation targets', 'manipulation target names', 'manipulation target GPRs']
    df.to_csv(f'./{output_path}.csv', encoding = 'cp949')

    df_cell = pd.DataFrame()

    for rxn in model.reactions:
        df_cell[rxn.id] = [rxn_lb[rxn.id], vo_l[rxn.id], vo_u[rxn.id], rxn_ub[rxn.id]]

    df_cell.index = ['V_min', 'Vo_lb', 'Vo_ub', 'V_max']
    
    df_cell = df_cell.transpose()
    df_cell.to_csv(f'./{output_path}_cellular_state.csv', encoding = 'cp949')

    print()
    print(df)
    print()
    print('status:', opt_status)


