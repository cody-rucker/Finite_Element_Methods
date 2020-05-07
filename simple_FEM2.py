"""
Author: ME
DATE: Nov. 29th 2016
"""
#solve -u_xx - u_yy = f on unit square
# u = 0 on the boundary

import numpy as np

h = 1/4  #grid spacing
M = int(1/h) #M+1 is number of nodes in each direction
N = (M-1)**2 #total number of interior nodes
A = np.zeros((N,N))
V = -1*np.ones((M+1,M+1))
for i in range(M-1):
	for j in range(M-1):
		V[j+1,i+1] = i*(M-1) + j

print(V)

lower_nodes = np.zeros((M**2,3))
upper_nodes = np.zeros((M**2,3))
lower_barycenters = np.zeros((M**2,2))
upper_barycenters = np.zeros((M**2,2))

for i in range(M):
	for j in range(M):
		lower_nodes[i*M+j,:] = [V[j,i],V[j,i+1],V[j+1,i+1]]
		upper_nodes[i*M+j,:] = [V[j,i],V[j+1,i+1],V[j+1,i]]
		lower_barycenters[i*M+j,:] = [2*h/3 + i*h,h/3 + j*h]
		upper_barycenters[i*M+j,:] = [h/3+i*h, 2*h/3 + j*h]

Ke_lower = np.array(([0.5, -0.5, 0],[-0.5, 1, -0.5],[0, -0.5, 0.5]))
Ke_upper = np.array(([0.5, 0, -0.5],[0, 0.5, -0.5],[-0.5, -0.5, 1]))
#counter-clockwise orientation for all triangles

s = M**2   #number of lower triangles


def boundary_exile(node_triplet):  #this function takes in a node triplet and deletes any negative ones and then tells me the relevant indices of the element stiffness matrix to use
	indices = np.array([0, 1, 2])
	idx = np.array([[-1]])
       #print(node_triplet[0])
	for j in range(3):
		if node_triplet[j] > -1:
	#		print(node_triplet[j])
			idx = np.append(idx,[[indices[j]]],0)
	sidx = idx.shape
	sidx = int(sidx[0])
	idx = idx[1:sidx]
	nodes = node_triplet[idx]
	number_of_nodes = nodes.shape
	number_of_nodes = int(number_of_nodes[0])
	return nodes, number_of_nodes,idx

#nodes, num_nodes, idx = boundary_exile(np.array([-1, 8, 4]))

for i in range(s):
	nodes,num_nodes,idx = boundary_exile(lower_nodes[i,:])
	for m in range(num_nodes):
		for p in range(num_nodes):
			A[int(nodes[m]),int(nodes[p])] = A[int(nodes[m]),int(nodes[p])] + Ke_lower[idx[m],idx[p]]
	nodes,num_nodes,idx = boundary_exile(upper_nodes[i,:])
	for m in range(num_nodes):
		for p in range(num_nodes):
			A[int(nodes[m]),int(nodes[p])] = A[int(nodes[m]),int(nodes[p])] + Ke_upper[idx[m],idx[p]]

print(A)
#TO DO: assemble load vector b = ∱f ϕ dx
#TO DO: Solve Ax = b for some x
# Develop this to solve -Δu = f(x,y), where f(x,y) = sin(pi*x)*sin(pi*y)
#assembly by loop through elements b = b_1 + b_2 +...



#Below is just me testing out the built in norm function
v1 = np.array([0,1,2])
v2 = np.array([1,2,3])
print(v1)
err_v = np.linalg.norm(v1-v2,2)
print(err_v)
