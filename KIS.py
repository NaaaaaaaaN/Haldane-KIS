import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import inv
import matplotlib.pyplot as plt


folder_path = './data'
file_names = os.listdir(folder_path)
file_names = [i.split('.')[0] for i in file_names]
file_names.sort()

descriptions = ['mu1', 'mu2', 'mu3', 'mu4','mu5']
paraname = [['muhat', 'Kx', 'Ki'],['mum', 'Kx', 'Ki'],['a', 'b', 'c'],['mumax', 'alpha', 'xopt'],['gammamax', 'xstar']]
algae = ["Skeletonema costatum", "Isochrysis galbana", "Dunaliella salina}", "Platymonus subcordiformis}", "Chlorococcum sp. FACHB-1556", "Microcystis aeruginosa FACHB-905", "Microcystis wesenbergii FACHB-1112", "Scenedesmus obliquus FACHB-116"]


def mu1(para, x):
    y1 = para[0] * x / (x + para[1] + x**2/para[2])
    return y1

def mu2(para, x):
    y1 = para[0] * para[2] * x / ((x + para[1]) * (x + para[2]))
    return y1

def mu3(para, x):
    y1 = x / (para[0] * x**2 + para[1] * x + para[2])
    return y1

def mu4(para, x):
    y1 = para[0] * x / (x + para[0]/para[1] * (x/para[2] - 1)**2)
    return y1

def mu5(para, x):
    # y1 = para[0] * x * para[1] / (x + para[1])**2
    y1 = 4 * para[0] * x * para[1] / (x + para[1])**2
    return y1

def SSE(para, mu, data):
    if mu in globals():
        SSE = np.sum((globals()[mu](para, data[0]) - data[1])**2)
    return SSE

def dmu1(optpara, x):
# definition of the partial derivative of mu_1 w.r.t. each parameters.
    muhat = optpara[0]
    Kx = optpara[1]
    Ki = optpara[2]
    dmuhat = x / (x + Kx + x**2/Ki )
    dKx = -muhat * x / (x + Kx + x**2/Ki)**2
    dKi = muhat * x**3 / (Ki * x + Ki * Kx + x**2 )**2
    y = np.vstack([dmuhat, dKx, dKi]).T
    return y

def dmu2(optpara, x):
# definition of the partial derivative of mu_2 w.r.t. each parameters.
    mum = optpara[0]
    Kx = optpara[1]
    Ki = optpara[2]
    dmum = Ki * x / (x+Kx) / (x+Ki)
    dKx = -mum * Ki * x / (x+Kx)**2 / (x+Ki)
    dKi = mum * x**2 / (x+Kx) / (x+Ki)**2
    y = np.vstack([dmum, dKx, dKi]).T
    return y

def dmu3(optpara, x):
# definition of the partial derivative of mu_3 w.r.t. each parameters.
    a = optpara[0]
    b = optpara[1]
    c = optpara[2]
    da = -x**3 / (a * x**2 + b * x + c)**2
    db = -x**2 / (a * x**2 + b * x + c)**2
    dc = -x / (a * x**2 + b * x + c)**2
    y = np.vstack([da, db, dc]).T
    return y

def dmu4(optpara, x):
#  definition of the partial derivative of mu_4 w.r.t. each parameters.
    mumax = optpara[0]
    al = optpara[1]
    xopt = optpara[2]
    dmumax = al**2 * xopt**4 * x**2 / (mumax * (xopt-x)**2 + al * xopt**2 * x)**2
    dal = (mumax*xopt)**2 * x * (xopt-x)**2 / ( mumax * (xopt-x)**2 + al * xopt**2 * x)**2
    dxopt = -2 * mumax**2 * al * xopt * x**2 * (xopt-x) / ( mumax * (xopt-x)**2 + al * xopt**2 * x)**2
    y = np.vstack([dmumax, dal, dxopt]).T
    return y

def dmu5(optpara, x):
# definition of the partial derivative of mu_5 w.r.t. each parameters.
    gammax = optpara[0]
    xstar = optpara[1]
    # dgammax = xstar * x / (x+xstar)**2
    # dxstar = gammax * x * (xstar-x) / (x+xstar)**3
    dgammax = 4 * xstar * x / (x+xstar)**2
    dxstar = 4 * gammax * x * (xstar-x) / (x+xstar)**3
    y = np.vstack([dgammax, dxstar]).T
    return y

def computeSensitivity(data, mu, dmu, opt_result, z):
    optpara = opt_result.x
    Q = np.diag(np.array([max(max(data[1])/50, 0.1*i)**(-2) for i in data[1]]))
    status = 'ok'
    if mu in globals() and dmu in globals():
        FIM = globals()[dmu](optpara, data[0]).T @ Q @ globals()[dmu](optpara, data[0])
        d_FIM = np.linalg.det(FIM)
        if d_FIM > 10**(-20):
            V = inv(FIM)
            sd = np.sqrt(np.diagonal(V))
            sd = sd.reshape((len(optpara),1))
            R = V / (sd * sd.T)
            s = np.sqrt(((globals()[mu](optpara, data[0])-data[1]).T @ Q @(globals()[mu](optpara, data[0])-data[1]) ) / ( len(data[0])-len(optpara)))
            sgm = s * np.sqrt(len(data[0]) * np.diagonal(V))
            e = np.sqrt(np.sum(globals()[dmu](optpara, z)**2 * (sgm**2).T, axis=1))
        else:
            e = np.nan
            sgm = np.nan
            R = np.nan
            status = 'singular'
    return e, sgm, R, status

def computeAICc(opt_result, data):
    SSE = opt_result.fun
    optpara = opt_result.x
    AICc = len(data[0]) * np.log(SSE/len(data[0])) + 2*(len(optpara)+1) + 2*(len(optpara)+1)*(len(optpara)+2) / (len(data[0])-len(optpara)-2)
    return AICc

def computeBIC(opt_result, data):
    SSE = opt_result.fun
    optpara = opt_result.x
    BIC = len(data[0]) * np.log(SSE/len(data[0])) + (len(optpara)+1) * np.log(len(data[0]))
    return BIC

def computePEMAC(opt_result, e, data):
    if np.shape(e) != ():
        SSE = opt_result.fun
        optpara = opt_result.x
        PEMAC = 2*len(data[0]) * np.log(np.mean(e)) + len(data[0]) * np.log(SSE/len(data[0])) + 2*(len(optpara)+1)*(len(optpara)+2) / (len(data[0])-len(optpara)-2)
    else:
        PEMAC = np.nan
    return PEMAC

def plot_optimal_experiment(data, mu, dmu, optpara, z):
    h = 3600
    e, sigma, corr, status = computeSensitivity(data, mu, dmu, optpara, z)
    if status != 'singular':
        fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
        plt.plot(data[0], data[1] * h, '+', label='Measurement', markersize=10)
        if len(data[0]) == 12 :
            t = 2.201
        elif len(data[0]) == 13 :
            t = 2.179
        elif len(data[0]) == 22 :
            t = 2.080
        if mu in globals() and dmu in globals():
            plt.plot(z, globals()[mu](optpara.x, z) * h, '-', label='Estimation', linewidth=2)
            plt.plot(z, (globals()[mu](optpara.x, z) - t * e) * h, '--', label='Lower Bound', linewidth=2)
            plt.plot(z, (globals()[mu](optpara.x, z) + t * e) * h, '--', label='Upper Bound', linewidth=2)
        ax.legend()
        plt.show()

def plot_sensitivity(data, mu, dmu, optpara, z, paraname):
    e, sigma, corr, status = computeSensitivity(data, mu, dmu, optpara, z)
    if status != 'singular':
        fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
        if dmu in globals():
            for i in range(len(sigma)) :
                plt.plot(z, (sigma[i] / e * globals()[dmu](optpara.x, z)[:, i])**2, '-', label=f'{paraname[desp][i]}', linewidth=2)
        ax.legend()
        ax.set_xlim([0, 1600])
        ax.set_ylim([-0.05, 1.05])
        plt.show()

xt = np.array(range(-2,1600))


print(f'Experimental data name : {file_names}')
response = input("Please enter the experimental data filename : ")
while response not in file_names:
    print('Experimental data name : {file_names}')
    response = input("Please enter the experimental data filename : ")

if response in file_names:
    data_file = response
    print(f'Descriptions : {descriptions}')
    response = input("Please enter the description name : ")
    while response not in descriptions:
        print(f'Descriptions : {descriptions}')
        response = input("Please enter the description name : ")
    if response in descriptions:
        description = response

        with open(f'{folder_path}/{data_file}.txt', 'r') as file:
            content = file.read()
        lines = content.split('\n')
        lines = [i for i in lines if i !='']
        if data_file[0] == 'A':
            x = [float(i.split('\t')[0]) for i in lines]
            y = [float(i.split('\t')[1])/3600 for i in lines]
        else:
            x = [float(i.split(' ')[0]) for i in lines]
            y = [float(i.split(' ')[1])/3600 for i in lines]
        if y[0] < 0:
            y = np.array(y) - y[0]
        data = [np.array(x), np.array(y)]

        
        import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)

        if description != None:
            desp = int(description.split('u')[-1])
            mu = description
            dmu = 'd'+ description

            if data_file[0] == 'A':
                algae_test = algae[0]
            else:
                algae_test = algae[int(data_file[-1])]
            print(f'\nTitle : Experimental data of {algae_test} ({data_file}) with the description {description}\n')

            if desp < 5 :
                if desp == 1:
                    if data_file[0] == 'A' or data_file[-1] in ['2', '3', '4']:
                        x0 = [0, 1, 1]
                    elif data_file[-1] == '1':
                        x0 = [1, 4000, 50]
                    elif data_file[-1] == '5':
                        x0 = [-1, -1000, -10]
                    elif data_file[-1] == '6':
                        x0 = [1, 2, 0]
                    elif data_file[-1] == '7':
                        x0 = [0, 10, 1]
                elif desp == 2 :
                    x0 = [0, 10, 10]
                elif desp == 3:
                    if data_file[0] == 'A':
                        x0 = [0, 0, 100]
                    elif data_file[-1] == '6':
                        x0 = [0, 0, 0]
                    else:
                        x0 = [0, 0, 1]
                elif desp == 4:
                    if data_file[0] == 'A' or data_file[-1] == '1':
                        x0 = [0, 0, 100]
                    elif data_file[-1] in ['2', '4', '6', '7']:
                        x0 = [0, 0, 10]
                    else:
                        x0 = [0, 0, 100]
                optimal_result = minimize(SSE, x0, args=(mu, data), method='Nelder-Mead')
                desp = desp - 1
                print(f'Optimal parameters of description :\n')
                for i in range(len(paraname[desp])):
                    if i < len(paraname[desp])-1 : 
                        print(f'{paraname[desp][i]}: {optimal_result.x[i]},')
                    else:
                        print(f'{paraname[desp][i]}: {optimal_result.x[i]},\n')
                print(f'SSE : {optimal_result.fun}\n')
            else:
                x0 = [0, 5]
                optimal_result = minimize(SSE, x0, args=(mu, data), method='Nelder-Mead')
                desp = desp - 1
                print(f'Optimal parameters of description :\n'
                    f'{paraname[desp][0]:10}: {optimal_result.x[0]},\n'
                    f'{paraname[desp][1]:10}: {optimal_result.x[1]},\n')
                print(f'SSE : {optimal_result.fun}\n')
                
            e_data, sigma, corr, status = computeSensitivity(data, mu, dmu, optimal_result, data[0])
            AICc = computeAICc(optimal_result, data) 
            BIC = computeBIC(optimal_result, data)
            PEMAC = computePEMAC(optimal_result, e_data, data)

            if status == 'singular':
                print('Fisher information matrix is close to singular.')

            print(f'Correlation matrix of parameter estimation error of description:\n{corr}\n')
            print(f'AICc : {AICc}\n')
            print(f'BIC : {BIC}\n')
            print(f'PEMAC : {PEMAC}\n')

            plot_optimal_experiment(data, mu, dmu, optimal_result, xt)
            plot_sensitivity(data, mu, dmu, optimal_result, xt, paraname)


          

