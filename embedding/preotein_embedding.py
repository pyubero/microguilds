# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 09:51:20 2022

@author: logslab
"""



import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from sklearn.manifold import TSNE, Isomap, SpectralEmbedding, LocallyLinearEmbedding, MDS #LocallyLinearEmbedding

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

FILENAME = 'seqs_upper_ref.faln' # "./query3_query_fraction.faln"
MODEL = 'pam70' #['blosum45', 'blosum50', 'blosum62',
                 # 'blosum80', 'blosum90', 'pam250',
                 # 'pam30', 'pam70']
EMBEDDING_MODE = 'Isomap'


# Convert all sequences to upper-case letters
# with open(FILENAME,'r') as f:
#     with open('seqs_upper_qry.faln','w+') as f2:
#         new_label = f.readline().replace('\n','') 
#         new_seq   = f.readline().replace('\n','') 
        
#         while new_seq != '':
#             f2.write( new_label+'\n' )
#             f2.write( new_seq.upper()  +'\n' )
            
#             new_label = f.readline().replace('\n','') 
#             new_seq   = f.readline().replace('\n','') 
        


# Read alignment
# labels = []
# alseqs = []

# with open(FILENAME,'r') as f:
#     new_label = f.readline().replace('\n','') 
#     new_seq   = f.readline().replace('\n','') 
    
#     while new_seq != '':
#         labels.append( new_label )
#         alseqs.append( new_seq   )
        
#         new_label = f.readline().replace('\n','') 
#         new_seq   = f.readline().replace('\n','') 
        

# Read alignment
alignment=AlignIO.read(FILENAME, 'fasta')

# Compute distance matrix...
calculator = DistanceCalculator(MODEL)
dm = calculator.get_distance(alignment)

# ... convert data to labels and a symmetric matrix
labels = dm.names
NSEQS  = len(labels)

distances = np.zeros( (NSEQS, NSEQS ) )



positive = np.zeros( (NSEQS,) )
for jj , label in enumerate(labels):
    if 'alves' in label:
        positive[jj] = 1
    if label in ['31','19','72','khadka159416','khadka150231','khadka163194']:
        positive[jj] = 2
    if label in ['69','khadka186523','khadka145526']:
        positive[jj]=3
    if label in ['3','30','22','59','43','13','41','8','27','63']:
        positive[jj]=4


for rr, row in enumerate( dm.matrix ):
    for cc, value in enumerate(row):
        distances[rr, cc] = value
        distances[cc, rr] = value


# plt.imshow(distances); plt.colorbar( label='Dissimilarity')
# plt.xlabel('Sequence idx')
# plt.ylabel('Sequence idx')
# plt.title( MODEL )


initial_coords = np.zeros( (len(labels), 2) )
t = np.linspace(0,2*np.pi, len(labels) )
initial_coords[:,0] = np.cos(t)
initial_coords[:,1] = np.sin(t)

PERPLEXITY = 10
MAXITER = 10_000
NJOBS = 4

if EMBEDDING_MODE == 'TSNE':
    embedding = TSNE(n_components=2,
                      perplexity = PERPLEXITY,         # : int : 30
                      early_exaggeration = 12, # : float : 12
                      learning_rate='auto', 
                      n_iter = MAXITER,
                      n_iter_without_progress = 200,
                      metric = 'precomputed',
                      init = initial_coords,
                      random_state = 10,
                      verbose=1,
                      n_jobs = NJOBS)
    
elif EMBEDDING_MODE =='Isomap':
    embedding = Isomap( n_neighbors=PERPLEXITY,
                       radius=None,
                       n_components=2,
                       eigen_solver='auto',
                       tol=0, 
                       max_iter=MAXITER,
                       path_method='auto',
                       neighbors_algorithm='auto',
                       n_jobs=NJOBS,
                       metric='precomputed')

elif EMBEDDING_MODE == 'MDS':
    embedding = MDS( n_components=2,
                     metric = True,
                     n_init=10,
                     max_iter=MAXITER,
                     verbose=2,
                     n_jobs=NJOBS,
                     random_state=10,
                     dissimilarity='precomputed' )




X_embedded = embedding.fit_transform( distances )


############################################
############################################

def place_sequence( new_sequence, old_alignment, embedded_coordinates, model_distance, gamma=0.1, step=0.01, tol=1e-3, verbose=True):
    
    # Add new sequence to old alignment...
    old_alignment.append(SeqRecord(new_sequence, id='New_seq'))
    
    # Compute distance matrix...
    dm = DistanceCalculator(model_distance).get_distance(old_alignment)
    true_dist = dm.matrix[-1][:-1]
    
    # Remove sequence from alignment.
    old_alignment = old_alignment[:-1,:]
    
    # Feel free to change the metric for steepest *ascent*
    def metric(v):
        return np.corrcoef(v, true_dist)[0,1]
    
    def compute_distance( p0, points):
        return np.sqrt(np.sum((p0-points)**2, axis=1))

    
    # Steepest ascent through R2
    n_iter = 0                                      #... keeps track of the number of iterations
    p = np.array( [0,0])                            #... set starting point    
    r0 = metric( compute_distance(p, X_embedded) )  #... evaluate metric of starting point

    # Compute gradient
    p_x = p + np.array([ step, 0])
    p_y = p + np.array([ 0, step])
    r_x = metric( compute_distance(p_x, X_embedded))
    r_y = metric( compute_distance(p_y, X_embedded))
    grad = np.array([(r_x - r0), (r_y-r0)])/step
    
    # ... gradient norm is our "error" measure
    grad_norm = np.linalg.norm(grad)

    while grad_norm > tol:
        p = p + gamma*grad
        r0 = metric( compute_distance(p, X_embedded) )  #... evaluate metric of starting point
        
        # Compute gradient
        p_x = p + np.array([ step, 0])
        p_y = p + np.array([ 0, step])
        r_x = metric( compute_distance(p_x, X_embedded))
        r_y = metric( compute_distance(p_y, X_embedded))
        grad = np.array([(r_x - r0), (r_y-r0)])/step
        
        # ... gradient norm is our "error" measure
        grad_norm = np.linalg.norm(grad)
        n_iter += 1
        if verbose==1: 
            print(n_iter, p, r0, grad_norm)
            
        if n_iter>=1000:
            break
        
    if verbose == 1:
        print('')
        
    return p, r0, old_alignment
        

new_seq = '---------------------------------------------------------------------------------------------------MSWIRSSIH--YL----FIVVVAVNSTLLTINAGD---------------------------------SIFYSDWMWTSFVIFNLSQSTMLVVGAIYYLLFT--GVPGTATYYATIMTIYTWVAK-WAWIGGLGYPYDFMIVPVWIPSAM------LLDLAYWATRRNKHAAILIGGSLLGMSMPLFNMINL--LTVHDPLEM---------------AFKYPRPTLPAYLTPI------EPQVGKFYNSPVALA---AGISAVISVPMAALGAKL-NTWTYRWAA--AWS--------KW-D-----'



values, r0, alignment = place_sequence(new_seq, alignment, X_embedded, MODEL, gamma=0.1)
# alignment= alignment[:-1,:]


x_out = values[0]
y_out = values[1]

plt.figure( dpi=600)
plt.scatter(X_embedded[:,0], X_embedded[:,1], c=positive, s=12)
plt.plot(x_out, y_out,'r+')
plt.title("%s - %s" % (MODEL,EMBEDDING_MODE) )
# plt.colorbar()

# Print some tags
# for jj , label in enumerate(labels):
#     if 'alves' in label:
#         pass
#     elif 'khadka' in label:
#         pass
#     else:
#         plt.text( X_embedded[jj,0], X_embedded[jj,1], label, fontsize=4)

# plt.xticks( plt.xticks()[0],  labels='')
# plt.yticks( plt.yticks()[0],  labels='')



############################################
############################################
QUERY_FILENAME = './seqs_upper_qry.faln'
NSEQS = 200

query_labels = []
query_sequences=[]


with open(QUERY_FILENAME,'r') as f :
    for jj in range(NSEQS):
        query_labels.append( f.readline().replace('\n','') )
        query_sequences.append( f.readline().replace('\n','') )

print('Loaded %d query sequences.' % len(query_sequences))





# alignment= alignment[:-1,:]

query_positions = []
query_gof = []
for jj in range(NSEQS):
    print('Localizing sequence [%d] - %s' % (jj, query_labels[jj]) )
    values, r2, alignment = place_sequence( query_sequences[jj] , alignment, X_embedded, MODEL, gamma=0.1)    
    query_positions.append(values)
    query_gof.append(r2)




plt.figure( dpi=600)
plt.scatter(X_embedded[:,0], X_embedded[:,1], c=positive, s=12)
# for pt in query_positions:
#     plt.plot( pt[0], pt[1],'r+')
plt.title("%s - %s" % (MODEL,EMBEDDING_MODE) )

plt.xticks( plt.xticks()[0],  labels='')
plt.yticks( plt.yticks()[0],  labels='')












# # Read alignment
# new_alignment=AlignIO.read(FILENAME, 'fasta')
# new_alignment.append(SeqRecord(new_seq, id='New_seq'))

# # Compute distance matrix...
# new_dm = calculator.get_distance(new_alignment)

# new_true_dist = new_dm.matrix[-1][:-1]


# # Covering all space: really not efficient
# p = np.array( [0,0] )

# dist_embedded = compute_distance(p, X_embedded)


# N = 100
# new_x = np.linspace(-1,1, N)
# new_y = np.linspace(-1,1,N)
# D = np.zeros( (N,N) )
# for jj, _x in enumerate(new_x):
#     for kk, _y in enumerate(new_y):
#         p = [ _x, _y]
#         v = compute_distance(p, X_embedded)
#         D[jj,kk] = np.corrcoef( v, new_true_dist)[0,1]

# idx = np.argwhere(D==D.max())
# row = np.argwhere(D==D.max())[0,:][0]
# col = np.argwhere(D==D.max())[0,:][1]


# x_out = new_x[row]
# y_out = new_y[col]
# print(x_out, y_out, D.max() )




# Steepest ascent
# gamma = 10                              #... learning rate
# step  = 0.01                            #... initial step size
# tol = 1e-3

# p = np.array( [0,0])                    #... set starting point
# v = compute_distance(p, X_embedded)     #... compute distances in embedded space
# r0 = np.corrcoef( v, new_true_dist)[0,1]     #... compute correlation with true distances

# p_x = p + np.array([ step, 0])
# p_y = p + np.array([ 0, step])
# r_x = np.corrcoef(compute_distance(p_x, X_embedded), new_true_dist )[0,1]
# r_y = np.corrcoef(compute_distance(p_y, X_embedded), new_true_dist )[0,1]

# grad = np.array([(r_x - r0), (r_y-r0)])/step
# grad = grad

# grad_norm = np.linalg.norm(grad)
# niter = 0
# while grad_norm > tol:
    
#     p = p + gamma*step*grad
#     v = compute_distance(p, X_embedded)     #... compute distances in embedded space
#     r0 = np.corrcoef( v, new_true_dist)[0,1]     #... compute correlation with true distances
    
#     p_x = p + np.array([ step, 0])
#     p_y = p + np.array([ 0, step])
#     r_x = np.corrcoef(compute_distance(p_x, X_embedded), new_true_dist )[0,1]
#     r_y = np.corrcoef(compute_distance(p_y, X_embedded), new_true_dist )[0,1]
    
#     grad = np.array([(r_x - r0), (r_y-r0)])/step
#     grad = grad#/np.linalg.norm(grad)
    
    # grad_norm = np.linalg.norm(grad)
    # niter += 1
    # print(niter, p, r0, grad_norm)









