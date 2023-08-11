import os

import cobra
from cobra.io import read_sbml_model
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from utils import argument_parser




if __name__=='__main__':
    
    parser = argument_parser()
    options = parser.parse_args()

    model_dir = options.model_dir
    output_dir = options.output_dir
    objective = options.biomass_reaction
    target = options.objective_reaction
    num_pt = options.max_step
    oxygen_condition = options.oxygen_condition
    
    
    model = read_sbml_model(model_dir)
    model.objective = objective
    target_rxn = model.reactions.get_by_id(target)
    
    O_exchange = model.reactions.EX_o2_e
    if oxygen_condition == 'aerobic':
        O_exchange.bounds = (-1000, 1000) 

    elif oxygen_condition == 'anaerobic':
        O_exchange.bounds = (0, 1000)
        
        
    if output_dir[-1] == '/':
        output_dir = output_dir[:-1]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    ###########################################################################################

    # flux distribution 구하기
    print('running FBA...')

    with model as m:
        max_prod = m.slim_optimize()    # lys 생산 flux 최댓값

    solution_space = []
    flux_range = np.linspace(0, max_prod, num_pt)

    with model as m:
        for flux in tqdm(flux_range):    
            target_rxn.bounds = (flux, flux)
            sol = m.optimize()

            solution_space.append(sol)

    df = pd.concat([item.fluxes for item in solution_space], axis=1)
    df.columns = range(1, num_pt + 1)

    print('done.')


    # 반응의 방향이 바뀌는 것들 지우기
    print('deleting reactions that Vmax * Vmin < 0...')

    for rxn in tqdm(df.index):
        if min(df.loc[rxn])*max(df.loc[rxn]) < 0:
            df.drop(rxn, inplace = True)

    print('done.')


    # flux가 최대한 일정하게 변화하는 것만 남기기
    # 기울기 = 0 또는 R2 < 0.9 인 reaction 지우기
    print('conducting linear regression analysis...')

    reaction_list = df.index
    Gradient = []
    R2_score = []

    for reaction in tqdm(reaction_list):

        flux_dist = np.array(df.loc[reaction])

        lr = LinearRegression()

        x = flux_range.reshape(-1, 1)
        y = flux_dist.reshape(-1, 1)
        lr.fit(x, y)
        Gradient.append(lr.coef_[0,0])
        R2_score.append(lr.score(x, y))

    df['Gradient'] = Gradient # 기울기: 증가 감소 정도
    df['R2 score'] = R2_score # 결정계수: 얼마나 선형적으로 증가하는지

    idx = df[df['Gradient'] == 0].index
    df.drop(idx, inplace = True)

    idx = df[df['R2 score'] < 0.9].index
    df.drop(idx, inplace = True)
    df.drop(objective, inplace = True)  # objective reaction 제외

    print('done')

    # up/down-regulation target 구하기
    print('extracting target file...')

    df_flux = df[[i for i in range(1, num_pt+1)]]
    up_rxn = []

    for rxn in tqdm(df.index):
        if sum(df_flux.loc[rxn]) * df['Gradient'][rxn] > 0:
           up_rxn.append(rxn)

    up_regulation_data = df.loc[up_rxn]  # Goal of FSEOF
    down_regulation_data = df.drop(up_rxn)


    up_regulation_data.to_csv(f'{output_dir}/FESOF_up_regulation_targets.csv', encoding = 'cp949')
    down_regulation_data.to_csv(f'{output_dir}/FSEOF_down_regulation_targets.csv', encoding = 'cp949')

    print('done.')
