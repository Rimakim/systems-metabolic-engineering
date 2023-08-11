import os
import cobra
from cobra.io import read_sbml_model

import numpy as np
import pandas as pd
from gurobipy import *
from tqdm import tqdm
from utils import argument_parser




if __name__=='__main__':
    
    parser = argument_parser()
    options = parser.parse_args()

    model_dir = options.model_dir
    output_dir = options.output_dir
    biomass_lb = options.biomass_lb
    target_rxn = options.objective_reaction
    oxygen_condition = options.oxygen_condition
    cpu_num = options.cpu_num
    default_bound = options.default_bound
    iter_num = options.iter_num
    max_ko = options.max_knockout

    class OptKnock:

        def __init__(self, f_model, biomass_lb, target_rxn, max_ko, oxygen_condition, iter_num):

            self.model = read_sbml_model(f_model)
            self.target_rxn = target_rxn
            self.biomass_lb = biomass_lb
            self.max_ko = max_ko
            self.iter_num = iter_num
            
            # biomass, target
            biomass_rxn = None

            for rxn in self.model.reactions:
                if rxn.objective_coefficient == 1:
                    biomass_rxn = rxn.id

            self.biomass_rxn = biomass_rxn
            self.model.reactions.get_by_id(self.biomass_rxn).lower_bound = self.biomass_lb

            self.oxygen_condition = oxygen_condition
            # oxygen condition
            if oxygen_condition == 'aerobic':
                pass

            elif oxygen_condition == 'anaerobic':
                self.model.reactions.EX_o2_e.knock_out()

        def model_setting(self, cpu_num, default_bound):
            self.cpu_num = cpu_num
            self.ub_default = default_bound
            self.lb_default = -default_bound

            self.m = Model()
            self.m.reset()
            self.m.params.Threads = cpu_num

            self.m.params.FeasibilityTol = 1e-9
            self.m.params.OptimalityTol = 1e-9
            self.m.params.BarConvTol = 1e-12
            self.m.params.IntFeasTol = 1e-9
            self.m.update()


        def add_constraints(self):

            # S_matrix
            S_matrix = {}

            for rxn in self.model.reactions:

                reactant_list = list(rxn.reactants)
                product_list = list(rxn.products)

                reactant_coff_list = list(rxn.get_coefficients(reactant_list))
                product_coff_list = list(rxn.get_coefficients(product_list))

                met_list = reactant_list + product_list
                coff_list = reactant_coff_list + product_coff_list

                for met, coff in zip(met_list, coff_list):
                    S_matrix[(rxn.id, met.id)] = coff

            pair, coff_value = multidict(S_matrix)
            pair = tuplelist(pair)

            # variables
            rxn_ub = {}
            rxn_lb = {}

            self.v = {}
            vf = {}
            vr = {}

            self.y = {}

            ds = {}
            qu = {}
            ql = {}
            zu = {}
            zl = {}

            for rxn in self.model.reactions:
                rxn_lb[rxn.id] = rxn.lower_bound
                rxn_ub[rxn.id] = rxn.upper_bound

            for rxn in self.model.reactions:
                self.v[rxn.id] = self.m.addVar(lb = rxn_lb[rxn.id], ub = rxn_ub[rxn.id], name = f'{rxn.id} flux')
                vf[rxn.id] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} forward')
                vr[rxn.id] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} reverse')

                self.y[rxn.id] = self.m.addVar(vtype = GRB.BINARY, name = f'{rxn.id} on/off') # on = 0, off = 1

                qu[rxn.id] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} ub dual')
                ql[rxn.id] = self.m.addVar(lb = 0, ub = self.ub_default, name = f'{rxn.id} lb dual')

                zu[rxn.id] = self.m.addVar(lb = self.lb_default, ub = self.ub_default, name = f'{rxn.id} ub z')
                zl[rxn.id] = self.m.addVar(lb = self.lb_default, ub = self.ub_default, name = f'{rxn.id} low z')

            for met in self.model.metabolites:
                ds[met.id] = self.m.addVar(lb = self.lb_default, ub = self.ub_default, name = f'{met.id} stoich')

            self.m.update()

            # primal constraints
            # 1. flux on/off constraints
            for rxn in self.model.reactions:
                rxn_id = rxn.id
                self.m.addConstr(self.v[rxn_id] - vf[rxn_id] + vr[rxn_id] == 0)
                self.m.addConstr(vf[rxn_id] - vr[rxn_id] >= rxn_lb[rxn_id] * (1-self.y[rxn_id]))
                self.m.addConstr(vf[rxn_id] - vr[rxn_id] <= rxn_ub[rxn_id] * (1-self.y[rxn_id]))

            # 2. steady state constraints
            for met in tqdm(self.model.metabolites):
                self.m.addConstr(quicksum((vf[rxn_id] - vr[rxn_id]) * coff_value[rxn_id, met_id] for rxn_id, met_id in pair.select('*', met.id)) == 0, name = f'{met.id} steay state')

            self.m.update()

            # dual constraints
            # 1. dual stoichiometric constraints
            for rxn in tqdm(self.model.reactions):
                if rxn.id == self.biomass_rxn:
                    self.m.addConstr(quicksum(coff_value[rxn_id, met_id] * ds[met_id] for rxn_id, met_id in pair.select(rxn.id, '*'))
                                + qu[rxn.id] - ql[rxn.id] == 1, name = 'biomass ds')
                else:
                    self.m.addConstr(quicksum(coff_value[rxn_id, met_id] * ds[met_id] for rxn_id, met_id in pair.select(rxn.id, '*'))
                                + qu[rxn.id] - ql[rxn.id] == 0, name = f'{rxn.id} ds')

            self.m.update()

            # 2. z variables (z = qu*y)
            for rxn in self.model.reactions:
                rxn_id = rxn.id

                self.m.addConstr(zu[rxn_id] >= 0)
                self.m.addConstr(zu[rxn_id] <= self.ub_default * self.y[rxn_id])
                self.m.addConstr(qu[rxn_id] - self.ub_default * (1-self.y[rxn_id]) <= zu[rxn_id])
                self.m.addConstr(zu[rxn_id] <= qu[rxn_id])

                self.m.addConstr(zl[rxn_id] >= 0)
                self.m.addConstr(zl[rxn_id] <= self.ub_default * self.y[rxn_id])
                self.m.addConstr(ql[rxn_id] - self.ub_default * (1-self.y[rxn_id]) <= zl[rxn_id])
                self.m.addConstr(zl[rxn_id] <= ql[rxn_id])

            self.m.update()

            # max KO
            self.m.addConstr(quicksum(self.y[rxn.id] for rxn in self.model.reactions) <= self.max_ko, name = 'max KO')

            # max primal obj = min dual obj
            self.m.addConstr(quicksum(rxn_ub[rxn.id]*(qu[rxn.id]-zu[rxn.id]) - rxn_lb[rxn.id]*(ql[rxn.id]-zl[rxn.id]) for rxn in self.model.reactions) == vf[self.biomass_rxn] - vr[self.biomass_rxn], 
                             name = 'biomass')
            self.m.update()


            # essential, no gpr associated genes => always on
            essential_list = []
            for rxn in tqdm(self.model.reactions):
                with self.model as m:

                    rxn.knock_out()
                    sol = m.slim_optimize()

                    if sol <= self.biomass_lb:
                        essential_list.append(rxn.id)


            no_gpr_associated_list = [] 
            for rxn in self.model.reactions:
                if rxn.gene_reaction_rule == '':
                    no_gpr_associated_list.append(rxn.id)

            always_on_list = essential_list + no_gpr_associated_list

            for rxn_id in always_on_list:
                self.m.addConstr(self.y[rxn_id] <= 0.1) 


            self.m.setObjective(vf[self.target_rxn] - vr[self.target_rxn], GRB.MAXIMIZE)
            self.m.update()


        def return_result(self):
            
            KO_targets_iter = []
            biomass_val_iter = []
            target_val_iter = []
            
            for i in range(self.iter_num):  
                self.m.optimize()
                
                if self.m.status == 2:

                    v_values = {rxn_id : v_val.x for rxn_id, v_val in self.v.items()}
                    y_values = {rxn_id : y_val.x for rxn_id, y_val in self.y.items()}

                    KO_targets = []
                    for rxn_id, y_val in y_values.items():
                        if y_val > 0.9:
                            KO_targets.append(self.model.reactions.get_by_id(rxn_id))
                    
                    self.m.addConstr(quicksum(self.y[KO_targets[i].id] for i in range(len(KO_targets))) <= len(KO_targets) - 0.1) 

                    KO_targets_iter.append(KO_targets)

                    biomass_val = v_values[self.biomass_rxn]
                    biomass_val_iter.append(biomass_val)

                    target_val = v_values[self.target_rxn]
                    target_val_iter.append(target_val)
               
                else:
                    break
                
            return biomass_val_iter, target_val_iter, KO_targets_iter
        
    
    # run OptKnock
    op = OptKnock(f_model = model_dir, 
                  biomass_lb = biomass_lb, 
                  target_rxn = target_rxn, 
                  max_ko = max_ko,
                  oxygen_condition = oxygen_condition,
                  iter_num = iter_num)

    op.model_setting(cpu_num = cpu_num, default_bound = default_bound)
    op.add_constraints()
    
    biomass_val_iter, target_val_iter, KO_targets_iter = op.return_result()
    
        
    df = pd.DataFrame()
    
    if len(KO_targets_iter) < len(target_val_iter):
        KO_targets_iter = KO_targets_iter + (len(target_val_iter) - len(KO_targets_iter)) * ['-']
    
    for i in range(len(KO_targets_iter)):
        if KO_targets_iter[i] == '-':
            pass
        else:
            df[f'step{i+1}'] = [biomass_val_iter[i], 
                                target_val_iter[i], 
                                [rxn.id for rxn in KO_targets_iter[i]], 
                                [rxn.name for rxn in KO_targets_iter[i]], 
                                [rxn.gene_reaction_rule for rxn in KO_targets_iter[i]]]
        
        
    df.index = ['biomass', 'objective', 'KO target id', 'KO target name', 'KO target GPR']
    df.to_csv(f'./{output_dir}.csv', encoding = 'cp949')
