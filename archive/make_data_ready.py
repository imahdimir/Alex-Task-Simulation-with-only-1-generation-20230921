from pathlib import Path

import pandas as pd

sims = Path(
        '/Users/mahdi/Dropbox/1-GitHub/Alex-Task-Simulation-with-only-1-generation-20230921/1/simulations')

def read_r_delta(fp) :
    df = pd.read_csv(sims / fp , sep = '\s')
    return df['r_delta'].iloc[-1]

def read_r_delta_estimate(fp) :
    df = pd.read_csv(sims / fp , sep = '\s' , index_col = 0)
    return df.loc['r']['estimate'] , df.loc['r']['SE']

def main() :
    pass

    ##

    # Read the data
    df = pd.DataFrame(columns = ['r_y' , 'v_indir' , 'r_dir_indir'])

    df['r_y'] = [0 , .25 , .5 , .75] * 4
    df['v_indir'] = [0] * 4 + [.25] * 12
    df['r_dir_indir'] = [0] * 8 + [.5] * 4 + [1] * 4

    ##
    df = df.astype(str)
    df = df.replace({
            '0.0' : '0' ,
            '1.0' : '1'
            })

    ##
    df[['r_delta' , 'r_delta_SE' , 'r_delta_hat' , 'r_delta_hat_SE']] = None

    ##
    for i in range(16) :
        if df.at[i , 'v_indir'] == '0' :
            fp = f'r_y_{df.at[i , "r_y"]}_VCS.txt'
            df.at[i , 'r_delta'] = read_r_delta(fp)

            fp = f'r_y_{df.at[i , "r_y"]}_pgi_v1.am_adj_pars.txt'
            df.at[i , 'r_delta_hat'] , df.at[
                i , 'r_delta_hat_SE'] = read_r_delta_estimate(fp)

        else :
            fp = f'v_indir_{df.at[i , "v_indir"]}_r_dir_indir_{df.at[i , "r_dir_indir"]}_r_y_{df.at[i , "r_y"]}_VCS.txt'
            df.at[i , 'r_delta'] = read_r_delta(fp)

            fp = f'v_indir_{df.at[i , "v_indir"]}_r_dir_indir_{df.at[i , "r_dir_indir"]}_r_y_{df.at[i , "r_y"]}_pgi_v1.am_adj_pars.txt'
            df.at[i , 'r_delta_hat'] , df.at[
                i , 'r_delta_hat_SE'] = read_r_delta_estimate(fp)

    ##
    df = df.fillna('0')

    ##
    df.to_csv('data.csv' , index = False)

    ##
    import matplotlib.pyplot as plt
    import numpy as np

    fig , ax = plt.subplots()

    plt.scatter(df['r_delta'] , df['r_delta_hat'] , c = 'black')

    x = np.linspace(*ax.get_xlim())
    ax.plot(x , x)

    plt.xlabel('r_delta')
    plt.ylabel('r_delta_hat')

    plt.savefig('r_delta_vs_r_delta_hat.png')
    plt.show()

    ##

    ##
