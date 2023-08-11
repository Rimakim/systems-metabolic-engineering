import os

import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

from utils import argument_parser




if __name__=='__main__':
    
    # terminal에서 arg 받아오기
    parser = argument_parser()
    options = parser.parse_args()

    model_dir = options.model_dir
    output_dir = options.output_dir
    objective = options.biomass_reaction
    target = options.objective_reaction
    num_pt = options.max_step
    oxygen_condition = options.oxygen_condition
    
    # option 정의
    model = read_sbml_model(model_dir)      # 균주 모델
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
        
    #########################################################################################3    
    
    with model as m:
        m.objective = target
        max_prod = m.slim_optimize()    # lys 생산 flux 최댓값

    flux_range = np.linspace(0, max_prod, num_pt)
    
    # bias, capacity 구해야 함
    print('running FVA...')
    
    solution_space = []
    inappro_rxn = []
    with model as m:            
        for i, flux in tqdm(enumerate(flux_range)):
        
            target_rxn.bounds = (flux, flux)
            sol = flux_variability_analysis(m, m.reactions[:], fraction_of_optimum=1.0, loopless = False)
        
            for rxn in sol.index:
                if min(sol.loc[rxn])*max(sol.loc[rxn]) < 0:     #min, max 값이 다르면 역반응 일어날 수 있다는 의미이므로 삭제 
                    inappro_rxn.append(rxn)
            
            sol.rename(columns = {'minimum':f'minimum{i+1}'}, inplace = True)
            sol.rename(columns = {'maximum':f'maximum{i+1}'}, inplace = True)
            sol[f'bias {i+1}'] = sol.mean(axis = 1)
            sol[f'capacity {i+1}'] = abs(sol[f'maximum{i+1}'] - sol[f'minimum{i+1}'])

            solution_space.append(sol)
    
    print('done.')
    
    df_flux = pd.concat([item for item in solution_space], axis=1)
    df_flux.drop(inappro_rxn, inplace = True)
    
    
    #각각의 값별로 df를 나눈다.
    df_min_flux = df_flux[[f'minimum{i+1}' for i in range(num_pt)]]
    df_max_flux = df_flux[[f'maximum{i+1}' for i in range(num_pt)]]
    df_flux_bias = df_flux[[f'bias {i+1}' for i in range(num_pt)]]
    df_flux_capacity = df_flux[[f'capacity {i+1}' for i in range(num_pt)]]
    
    
    # df_min, df_max, df_bias 중에서 최대 최소 부호가 다른 반응은 모두 지운다.
    # 반응의 방향이 일정하지 않다는 뜻이므로!
    print('deleting reversible reactions...')
    
    delete_reaction = []

    for idx, rxn in enumerate(df_flux.index):
        if min(df_min_flux.loc[rxn]) * max(df_min_flux.loc[rxn]) < 0:
            delete_reaction.append(rxn)

        elif min(df_max_flux.loc[rxn]) * max(df_max_flux.loc[rxn]) < 0:
            delete_reaction.append(rxn)

        elif min(df_flux_bias.loc[rxn]) * max(df_flux_bias.loc[rxn]) < 0:
            delete_reaction.append(rxn)

    delete_reaction = list(set(delete_reaction))
    
    df_min_flux.drop(delete_reaction, inplace = True)
    df_max_flux.drop(delete_reaction, inplace = True)
    df_flux_bias.drop(delete_reaction, inplace = True)
    df_flux_capacity.drop(delete_reaction, inplace = True)
    
    print('done.')
    
    
    # Linear regression 통해 slope와 R2 score 구한다
    print('running Linear regression analysis...')
    
    flux_slope = []
    R2_score = []
    for rxn in tqdm(df_flux_bias.index):

        lr = LinearRegression()
        x = flux_range.reshape(-1, 1)
        y = np.array(df_flux_bias.loc[rxn]).reshape(-1, 1)

        lr.fit(x, y)
        flux_slope.append(lr.coef_[0,0])
        R2_score.append(lr.score(x, y))

    df_flux_bias['flux slope'] = flux_slope
    df_flux_bias['R2 score'] = R2_score
    
    
    # 꾸준히 증가하는 반응만을 남긴다. (EFSEOF에서와 같은 기준으로 적용)
    idx = df_flux_bias[df_flux_bias['R2 score'] < 0.9].index

    df_min_flux.drop(idx, inplace = True)
    df_max_flux.drop(idx, inplace = True)
    df_flux_bias.drop(idx, inplace = True)
    df_flux_capacity.drop(idx, inplace = True)
    
    print('done.')
    
    # up/down regulation target 구분하기
    # flux가 계속 0인 것은 담지 않는다.
    print('extracting results...')
    
    df_bias_dist = df_flux_bias[[f'bias {i}' for i in range(1, num_pt+1)]]
    df_slope = df_flux_bias['flux slope']

    up_target_rxn = []
    down_target_rxn = []
    for rxn in df_flux_bias.index:
        if sum(df_bias_dist.loc[rxn]) * df_slope.loc[rxn] > 0:
            up_target_rxn.append(rxn)

        elif sum(df_bias_dist.loc[rxn]) * df_slope.loc[rxn] < 0:
            down_target_rxn.append(rxn)

    up_target_min = df_min_flux.loc[up_target_rxn]
    up_target_max = df_max_flux.loc[up_target_rxn]
    up_target_bias = df_flux_bias.loc[up_target_rxn]
    up_target_capacity = df_flux_capacity.loc[up_target_rxn]
    
    # output file
    up_target_min.to_csv(f'{output_dir}/FESOF_up_regulation_min.csv', encoding = 'cp949')
    up_target_max.to_csv(f'{output_dir}/FESOF_up_regulation_max.csv', encoding = 'cp949')
    up_target_bias.to_csv(f'{output_dir}/FESOF_up_regulation_bias.csv', encoding = 'cp949')
    up_target_capacity.to_csv(f'{output_dir}/FESOF_up_regulation_capacity.csv', encoding = 'cp949')
    
    print('done')
