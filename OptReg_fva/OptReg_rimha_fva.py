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
    
        def __init__(self, model_path, biomass_lb, target_rxn, L, C, eps, oxygen_condition, fraction_of_optimum, iter_num):
    
            self.model = read_sbml_model(model_path)
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

            # oxygen_condition
            if self.oxygen_condition  == 'anaerobic':
                self.model.reactions.EX_o2_e.lower_bound = 0
            else:
                pass
    
        
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
                    
                    ub = m.optimize(objective_sense = 'maximize')
                    self.rxn_ub[rxn.id] = ub.objective_value
        
                    lb = m.optimize(objective_sense = 'minimize')
                    self.rxn_lb[rxn.id] = lb.objective_value

            
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
                if self.rxn_ub[rxn.id] >= 0 and self.rxn_lb[rxn.id] >= 0:
                    self.for_list.append(rxn.id)
            
                elif self.rxn_ub[rxn.id] <= 0 and self.rxn_lb[rxn.id] <= 0:
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
            
            self.flux = {}
            self.v = {}
    
            '''
            <primal variable addition>
            
            only forward: rxn name_for
            only reverse: rxn name_rev
            reversible: both are included
            '''
            # true flux value
            for rxn in self.model.reactions:
                self.flux[rxn.id] = self.m.addVar(lb = self.rxn_lb[rxn.id], ub = self.rxn_ub[rxn.id], name = f'{rxn.id} flux')
                
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
            self.var_o_u = {} # steady-state up
            self.var_o_l = {} # steady-state low
            
            for var in self.v.keys():
                
                if var[:-4] in self.for_list:
                    self.v_ub[var] = self.rxn_ub[var[:-4]]
                    self.v_lb[var] = self.rxn_lb[var[:-4]]
                    self.var_o_u[var] = self.vo_u[var[:-4]]
                    self.var_o_l[var] = self.vo_l[var[:-4]]
                
                elif var[:-4] in self.rev_list:
                    self.v_ub[var] = abs(self.rxn_lb[var[:-4]])
                    self.v_lb[var] = abs(self.rxn_ub[var[:-4]])
                    self.var_o_u[var] = abs(self.vo_l[var[:-4]])
                    self.var_o_l[var] = abs(self.vo_u[var[:-4]])

                # for reversible reactions
                else:
                    # forward
                    if '_for' in var:
                        self.v_ub[var] = self.rxn_ub[var[:-4]]
                        self.v_lb[var] = 0
                        
                        if self.vo_u[var[:-4]] >= 0 and self.vo_l[var[:-4]] >= 0:
                            self.var_o_u[var] = self.vo_u[var[:-4]]
                            self.var_o_l[var] = self.vo_l[var[:-4]]

                        elif self.vo_u[var[:-4]] >= 0 and self.vo_l[var[:-4]] < 0:
                            self.var_o_u[var] = self.vo_u[var[:-4]]
                            self.var_o_l[var] = 0 
                            
                        else:
                            self.var_o_u[var] = 0
                            self.var_o_l[var] = 0

                    # reverse
                    else:
                        self.v_ub[var] = abs(self.rxn_lb[var[:-4]])
                        self.v_lb[var] = 0   
                        
                        if self.vo_u[var[:-4]] <= 0 and self.vo_l[var[:-4]] <= 0:
                            self.var_o_u[var] = abs(self.vo_l[var[:-4]])
                            self.var_o_l[var] = abs(self.vo_u[var[:-4]])

                        elif self.vo_u[var[:-4]] > 0 and self.vo_l[var[:-4]] <= 0:
                            self.var_o_u[var] = abs(self.vo_l[var[:-4]]) 
                            self.var_o_l[var] = 0
                            
                        else:
                            self.var_o_u[var] = 0
                            self.var_o_l[var] = 0
    
            
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
            
            for rxn in self.v.keys():
            
                self.qku[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO ub dual')
                self.qkl[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} KO lb dual')
            
                self.quu[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} up ub dual')
                self.qul[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} up lb dual')
                
                self.qdu[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} down ub dual')
                self.qdl[rxn] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn} down lb dual')
            
                self.yu[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')     # on = 0, off = 1
                self.yd[rxn] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn} up on/off')     # on = 0, off = 1
               
                
                self.zku[rxn] = self.m.addVar(lb = -1*self.ub_default, ub = self.ub_default, name = f'{rxn} KO ub z')
                self.zkl[rxn] = self.m.addVar(lb = -1*self.ub_default, ub = self.ub_default, name = f'{rxn} KO lb z')
                self.zul[rxn] = self.m.addVar(lb = -1*self.ub_default, ub = self.ub_default, name = f'{rxn} up lb z')
                self.zdu[rxn] = self.m.addVar(lb = -1*self.ub_default, ub = self.ub_default, name = f'{rxn} down ub z')
            
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
    
            S_matrix = {}
            
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
                        S_matrix[(var, met.id)] = coff
                
                else:
                    reactant_list = list(rxn.products)
                    product_list = list(rxn.reactants)
                    
                    reactant_coff_list = [-1*coff for coff in rxn.get_coefficients(reactant_list)]
                    product_coff_list = [-1*coff for coff in rxn.get_coefficients(product_list)]
                                             
                    met_list = reactant_list + product_list
                    coff_list = reactant_coff_list + product_coff_list
                    
                                         
                    for met, coff in zip(met_list, coff_list):
                        S_matrix[(var, met.id)] = coff
                        
            
            pair, self.coff_value = multidict(S_matrix)
            self.pair = tuplelist(pair)
    
        
        def add_constraints(self):
    
            '''
            <primal constraints>
            0. flux = vf-vr
            1. Stoichiometric constraints 
            2. KO constraints
            3. UP_regulation constraints
            4. DOWN_regulation constraints
            5. minimal biomass, ATPM constraints
            '''

            # 0. flux constraints
            
            for rxn in self.flux.keys():

                if rxn in self.for_list:
                    self.m.addConstr(self.flux[rxn] - self.v[f'{rxn}_for'] == 0)

                elif rxn in self.rev_list:
                    self.m.addConstr(self.flux[rxn] + self.v[f'{rxn}_rev'] == 0)

                else:
                    self.m.addConstr(self.flux[rxn] - self.v[f'{rxn}_for'] + self.v[f'{rxn}_rev'] == 0)
            
            
            # 1. Stoichiometric constraints
            
            for met in tqdm(self.model.metabolites):
                self.m.addConstr(quicksum( self.v[val_id]*self.coff_value[val_id, met_id] for val_id, met_id in self.pair.select('*', met.id) ) == 0,
                           name = f'{met.id} mass balance')
            
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
                    self.m.addConstr(quicksum(self.coff_value[var_id, met_id] * self.ds[met_id] for var_id, met_id in self.pair.select(var, '*')) 
                                + self.qku[var] - self.qkl[var] + self.quu[var] - self.qul[var] + self.qdu[var] - self.qdl[var] - self.d_bio_atp[var] >= 1 - self.eps)
                
                elif var[:-4] == 'ATPM':
                    self.m.addConstr(quicksum(self.coff_value[val_id, met_id] * self.ds[met_id] for val_id, met_id in self.pair.select(var, '*')) 
                                + self.qku[var] - self.qkl[var] + self.quu[var] - self.qul[var] + self.qdu[var] - self.qdl[var] - self.d_bio_atp[var] >= -1*self.eps)
                
                else:
                    self.m.addConstr(quicksum(self.coff_value[val_id, met_id] * self.ds[met_id] for val_id, met_id in self.pair.select(var, '*')) 
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
            primal_obj = self.v[f'{self.biomass_rxn}_for'] - quicksum(self.v[var] for var in self.v.keys())*self.eps
            
            dual_obj = (-1*( ( self.biomass_lb * self.d_bio_atp[f'{self.biomass_rxn}_for'] ) + ( self.v_lb['ATPM_for'] * self.d_bio_atp['ATPM_for'] ) )
                        + quicksum( ( (self.zku[var]*self.v_ub[var]) - (self.zkl[var]*self.v_lb[var]) ) for var in self.v.keys() )
                        + quicksum( self.quu[var]*self.v_ub[var] - self.zul[var]*self.v_lb[var] + ( self.var_o_u[var]*(1-self.C) + self.v_ub[var]*self.C )*( self.zul[var] - self.qul[var] ) for var in self.v.keys() )
                        + quicksum( self.zdu[var]*self.v_ub[var] - self.qdl[var]*self.v_lb[var] + ( self.var_o_l[var]*(1-self.C) + self.v_lb[var]*self.C )*( self.qdu[var] - self.zdu[var] ) for var in self.v.keys() )
                       )
            
            self.m.addConstr(primal_obj == dual_obj)
            
            self.m.update()
    
    
            '''
            <Constraints for binary variables>
            
            1. Each reaction can have only one(or zero) regulation state.
            2. Maximum regulations
            3. In reversible reactions, forward and reverse reactions cannot be up or down-regulated at the same time.
            4. Reactions that are NOT linked with GPR cannot be manipulated. (yk, yu, yd values are always 1) 
            '''
    
            # 1. Each reaction can have only one(or zero) regulation state.
    
            for var in self.v.keys():
                if var[:-4] in self.for_rev_list:
                    self.m.addConstr((1-self.yk[var[:-4]]) + (1-self.yu[var]) + (1-self.yd[var]) <= 1)
                else:
                    self.m.addConstr((1-self.yk[var]) + (1-self.yu[var]) + (1-self.yd[var]) <= 1)
                    
            
            # 2. Maximum regulations
            
            rev = quicksum((1-self.yk[var[:-4]]) + (1-self.yu[var]) + (1-self.yd[var]) for var in self.v.keys() if var[:-4] in self.for_rev_list)
            irrev = quicksum((1-self.yk[var]) + (1-self.yu[var]) + (1-self.yd[var]) for var in self.v.keys() if var[:-4] not in self.for_rev_list)
            
            self.m.addConstr(rev + irrev <= self.L + 0.1)
            
            # 3. In reversible reactions, forward and reverse reactions cannot be up or down-regulated at the same time.
            for rxn in self.for_rev_list:
                
                self.m.addConstr(self.yu[f'{rxn}_for'] + self.yu[f'{rxn}_rev'] + self.yd[f'{rxn}_for'] + self.yd[f'{rxn}_rev'] >= 2.9)
               
                
            # 4. Reactions that are NOT linked with gene_reaction_rule cannot be manipulated. (yk, yu, yd values are always 1) 
    
            for val in self.v.keys():
                gpr = self.model.reactions.get_by_id(val[:-4]).gene_reaction_rule
                
                if gpr == '' or 's' in gpr or 'tex' in val:    # not linked GPR or spontaneous reactions    
                    self.m.addConstr(self.yu[val] == 1)
                    self.m.addConstr(self.yd[val] == 1)
            
                    if val[:-4] in self.for_rev_list:
                        self.m.addConstr(self.yk[val[:-4]] == 1)
                    else:
                        self.m.addConstr(self.yk[val] == 1)

                
                        
            
            self.m.update()
    
    
            '''
            <Objective function>
            Maximize target chemical
            '''
            
            self.m.setObjective(self.flux[self.target_rxn], GRB.MAXIMIZE)

            '''
            <Objective function>
            Maximize target chemical
            '''
            
            #if self.target_rxn in self.for_rev_list:
                #self.m.setObjective(self.v[f'{self.target_rxn}_for'] - self.v[f'{self.target_rxn}_rev'], GRB.MAXIMIZE)
                    
            #else:
                #self.m.setObjective(self.v[f'{self.target_rxn}_for'], GRB.MAXIMIZE)
    
    
        def return_results(self):
    
            manipulation_target_iter = []
            optimized_biomass_iter = []
            optimized_target_iter = []
    
            # finding multiple manipulation target sets 
            for i in range(self.iter_num):  
                
                print('iterations:', i+1)
                
                self.m.optimize()
                opt_status = self.m.status
                
                if opt_status in (2, 13):
    
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
    
                    print(manipulation_target)
                    
                    for tgt, ste in zip(manipulation_target, manipulation_state):
                        manipulation.append(tgt+ste)
    
                    manipulation_target_iter.append(manipulation)
                    optimized_biomass_iter.append(self.v[f'{self.biomass_rxn}_for'].x)
                    optimized_target_iter.append(self.v[f'{self.target_rxn}_for'].x)
                    
                    # manipulated reactions in previous step
                    for rxn_id in manipulation_target:

                        if '_for' not in rxn_id and '_rev' not in rxn_id:
                            self.m.addConstr(self.yu[f'{rxn_id}_for'] == 1)
                            self.m.addConstr(self.yd[f'{rxn_id}_rev'] == 1)

                        else:
                            self.m.addConstr(self.yu[rxn_id] == 1)
                            self.m.addConstr(self.yd[rxn_id] == 1)
                            

                        if rxn_id[:-4] in self.for_rev_list:
                            self.m.addConstr(self.yk[rxn_id[:-4]] == 1)
                            
                        else:
                            self.m.addConstr(self.yk[rxn_id] == 1)

                    self.m.update()
                        
                   
                else:
                    break
                    
            return manipulation_target_iter, optimized_biomass_iter, optimized_target_iter, opt_status

    
    '''
    Run OptReg and Output results
    '''
    
    optreg = OptReg(model_path = model_path,
                    biomass_lb = biomass_lb,
                    target_rxn = target_rxn,
                    L = max_manipulation,
                    C = regulation_strength,
                    eps = eps, 
                    oxygen_condition = oxygen_condition,
                    fraction_of_optimum = fraction_of_optimum,
                    iter_num = iter_num)
    
    rxn_ub, rxn_lb, vo_u, vo_l = optreg.calculate_cellular_state()
    
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
                    
    df.index = ['biomass', 'objective', 'manipulation targets', 'manipulation target names', 'manipulation target GPRs']
    df.to_csv(f'./{output_path}.csv', encoding = 'cp949')

    df_cell = pd.DataFrame()

    for rxn in model.reactions:
        df_cell[rxn.id] = [rxn_lb[rxn.id], vo_l[rxn.id], vo_u[rxn.id], rxn_ub[rxn.id]]

    df_cell.index = ['V_min', 'Vo_lb', 'Vo_ub', 'V_max']
    
    df_cell = df_cell.transpose()
    df_cell.to_csv(f'./{output_path}_cellular_state.csv', encoding = 'cp949')

    print(df)
    print(opt_status)



