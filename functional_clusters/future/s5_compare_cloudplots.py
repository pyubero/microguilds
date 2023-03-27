# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 19:23:55 2022

@author: Pablo
"""
from scipy.optimize import least_squares
import numpy as np
from matplotlib import pyplot as plt


def linear(p, t, y=0):
    return p[1] * t  +p[0] - y

def linear_strict(p, t, y=0):
    return p[0] * t  - y

def michaelis(p, t, y=0):
    return ( p[0] * t )/( p[1] + t ) - y

FILENAME_ENRICHMENT = 'data_clade_enrichment_potF.npz'
FILENAME_FILTNODES  = 'data_nodes_sign_filtered.npz'
FILENAME_TREECOMP_1 = 'data_tree_comparison_16S_rplB.npz'
FILENAME_TREECOMP_2 = 'data_tree_comparison_potF_16S.npz'


# Load enrichment data
data = np.load(FILENAME_ENRICHMENT)
F = data['F']
S = data['S']
features = data['features']
ZSCORES = data['ZSCORES'][1:,:]
ZSCORES[ np.isnan(ZSCORES)] = 0
# idx_sign = np.argwhere( np.any(np.abs(ZSCORES[:,6:9])>3, axis=1))[:,0]
# all_sign = idx_sign.copy()

# Load significant nodes filtered by parenthood 
# ... I subtract 1 because we do not want to plot node 0
data = np.load( FILENAME_FILTNODES, allow_pickle=True)
idx_sign = data['sign_nodes']
idx_sign = np.array([ np.array(_,dtype='int')-1 for _ in idx_sign])
all_sign = np.unique( [_ for a in idx_sign for _ in a] )


data = np.load( FILENAME_TREECOMP_2, allow_pickle=True)
n_bichos_potF = data['n_bichos']
depth_lca_16s_potF = data['depth_lca_16s']
depth_lca_gene_potF = data['depth_lca_gene']
leafs_idx_lca_potF = data['leafs_idx_lca']

x2 = depth_lca_16s_potF[1:]
y2 = depth_lca_gene_potF[1:]


data = np.load( FILENAME_TREECOMP_1, allow_pickle=True)
n_bichos_recA = data['n_bichos']
depth_lca_16s_recA = data['depth_lca_16s']
depth_lca_gene_recA = data['depth_lca_gene']
leafs_idx_lca_recA = data['leafs_idx_lca']

x1 = depth_lca_16s_recA[1:]
y1 = depth_lca_gene_recA[1:] 
x_ = np.linspace( 0, 1, 1000)


# Normalize??
y1 = y1/np.nanmax( y1 )
y2 = y2/np.nanmax( y2 )


# from sklearn.linear_model import LinearRegression 
# fit = LinearRegression(fit)

# res_lsq    = least_squares( michaelis, [0.1,0.1], args=(x2, y2))
res_rob1 = least_squares( linear_strict, [0.1,], loss='arctan', f_scale=0.1,
                         args=(x1, y1))
# res_rob2 = least_squares( michaelis, [1.0,1.0], args=(x2, y2))
y_rob = linear_strict( res_rob1.x, x_)
yp_rob= linear_strict( res_rob1.x, x1)

x1 = x1+np.random.randn( *x1.shape)/299
x2 = x2+np.random.randn( *x2.shape)/299



import statsmodels.api as sm
N = len(x1)
p = 2  # plus one because LinearRegression adds an intercept term

X_with_intercept = np.empty(shape=(N, p) )
X_with_intercept[:, 0] = 1
X_with_intercept[:, 1:p] = np.expand_dims(x1,axis=-1)

ols = sm.OLS(y1, X_with_intercept)
ols_result = ols.fit()
params = ols_result.params
ci = ols_result.conf_int(0.01)


def confidence_bands(x,y,xq,alfa=0.95):
    from scipy.stats import t
    x = x1
    y = y1
    alfa = 0.99
    # xq = np.linspace( np.nanmin(x), np.nanmax(x), 100)
    xq =x_
    xmean = np.nanmean(x)
    n = len(y)
    ypred = params[0] + params[1]*x
    tn_2 = t.ppf( (1+alfa)/2, n-2 )
    F1 = np.sum( (y-ypred)**2)/(n-2)
    F2 = 1/n + (xq-xmean)**2/np.sum( (x-xmean)**2 )
    
    DeltaY = tn_2*np.sqrt(F1*F2)
    return DeltaY
DeltaY  = confidence_bands(x1,y1,x_, alfa=0.99)


# Compute sigma of residuals NOT the error of the regression line
ypred = linear(params, x1)
s = 1.96*np.sqrt( np.sum( (y1-ypred)**2)/(len(y1)-1) ) 

# Clustering nodes
# >> OJO << 870 o 869??????
idx_cluster = []
# idx_cluster = -1+ np.array([37, 40, 257, 361, 866, 867, 869]) 
                         # 484, 1145, 1155, 10, 1])
# lbl_cluster = ['CI', 'CIa', 'CIb','CIc','CII','CIIa','CIIb','CIc_p','I1','I2','I3','I4']


# Bivariate
x1,y1
C = np.cov(x1,y1)
eigval, eigvec = np.linalg.eig(C)

x_vec= np.array([1 , 0]) # vector along x-axis
cosrotation = np.dot(x_vec,eigvec[1])/(np.linalg.norm(x_vec)*np.linalg.norm(eigvec[1])); 
rotation =np.pi/2-np.arccos(cosrotation); # rotation angle
xcent = np.nanmean(x1)
ycent = np.nanmean(y1)
alfa = eigvec[1,0]/eigvec[0,0]
y = alfa*(x_-xcent) + ycent



# Figure 1
# Show all significant ones
plt.figure( dpi=300)
plt.plot(x_, y_rob,'r')
plt.plot(x_, y,'g')
plt.plot(x_, linear(params, x_) ,'y')
p1=plt.plot( x1, y1,'.', alpha=0.5)
p2=plt.plot( x2, y2,'.', alpha = 0.5)
plt.plot( x2[all_sign], y2[all_sign],'.', c=p2[0].get_color(), mec='k')
plt.plot( x2[idx_cluster], y2[idx_cluster],'o', ms=5,mfc='none', mec='r', alpha=1)
# for jj, idx in enumerate(idx_cluster):
#     plt.text(x2[idx]+0.04, 
#              y2[idx]+np.random.randn()/60,
#              lbl_cluster[jj],
#              ha='center', va='center',
#              fontsize=8)


# plt.plot( x_, linear(params, x_), c = p1[0].get_color())
# taking the min and max of the ci of the params ---> bad! gotta take into account joint pdf
# taking into account the joint pdf ----> bad, this is not what we want!
# plt.gca().fill_between( x_, linear(params, x_)-DeltaY , linear(params, x_)+DeltaY,
#                         color= p1[0].get_color() , alpha=0.3)
# plt.gca().fill_between( x_, linear(params, x_)-s , linear(params, x_)+s,
#                         color= p1[0].get_color() , alpha=0.3)

plt.ylim( (0, plt.ylim()[1]) )


plt.grid()
plt.xlabel('Node relatedness in 16S')
plt.ylabel('Node relatedness')
plt.legend()


#%% Figure 2
# ## Multiplot with significancy per feature
plt.figure( figsize=(15,20), dpi=200)
for jj in range(len(features)):
    plt.subplot(5,3, jj+1)
    p1=plt.plot( x1, y1,'.', ms=2, alpha=0.5)
    # p2=plt.plot( x2, y2,'.', alpha = 0.5)
    # plt.scatter( x2, y2, s=1.5*ZSCORES[:,jj])
    plt.scatter( x2, y2,
                marker='o',
                c=p2[0].get_color(),
                s= 6*np.abs(ZSCORES[:,jj]) )
    plt.scatter( x2[idx_sign[jj]], y2[idx_sign[jj]],
                marker='o',
                c=p2[0].get_color(),
                s= 6*np.abs(ZSCORES[:,jj][idx_sign[jj]]),
                edgecolors='k',
                zorder=99)
    plt.gca().fill_between( x_, linear(params, x_)-s , linear(params, x_)+s,
                            color= p1[0].get_color() , alpha=0.15,
                            zorder=0)
    
    for idx in idx_sign[jj]:
        plt.text(x2[idx]+0.01, y2[idx]+0.01, "%d" % (idx+1), fontsize=8)
    plt.grid()    
    plt.text(0,1, features[jj])
    if jj==0:
        plt.legend(('rplB','potF'), loc='lower right')
    if jj>13:
        plt.xlabel('Node depth in 16S tree')
    if np.mod(jj,3)==0:
        plt.ylabel('Node depth in gene tree')



