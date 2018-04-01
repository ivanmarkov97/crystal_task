import pandas as pd
import numpy as np

def fill_dis(dis, matrix):
	for i in range(len(matrix)):
		dis[i] = sum(matrix[i])

def first_crystal_iter(dis, matrix, conn):
	for i in range(len(matrix)):
		conn[0][i] = 2 * matrix[i][0] - dis[i]
	conn[0][0] = None
	max_value = np.nanmax(conn[0, :])
	max_vert = np.argwhere(conn[0, :] == max_value).tolist()
	max_vert = [v[0] for v in max_vert]
	max_vert.insert(0, 0)
	return max_vert

def find_near_place(crystal):
	distances = []
	indexes = []
	for i in range(crystal.shape[0]):
		for j in range(crystal.shape[1]):
			if crystal[i][j] is None:
				distances.append((i*i + j*j)**0.5)
				indexes.append((i, j))
	ind_dis = distances.index(min(distances))
	return indexes[ind_dis]

def place_on_crystal(crystal, max_vert):
	for v in max_vert:
		pos_x, pos_y = find_near_place(crystal)
		crystal[pos_x][pos_y] = v

def add_max_vert(vert_on_crystal, max_vert):
	for v in max_vert:
		vert_on_crystal.append(v)
	vert_on_crystal.sort()

def crystal_iter(dis, matrix, conn, vert_on_crystal, iter_num):
	conn = np.vstack((conn, [0.0 for _ in range(len(matrix))]))
	conn[iter_num][vert_on_crystal] = None
	for i in range(len(matrix)):
		if conn[iter_num][i] != None:
			for v in vert_on_crystal:
				conn[iter_num][i] += matrix[i][v]
			conn[iter_num][i] *= 2
			conn[iter_num][i] -= dis[i]
	max_value = np.nanmax(conn[iter_num, :])
	max_vert = np.argwhere(conn[iter_num, :] == max_value).tolist()
	max_vert = [v[0] for v in max_vert]
	return max_vert, conn

def prepare_part(dis, matrix, conn, crystal, vert_on_crystal):
	fill_dis(dis, matrix)

	max_vert = first_crystal_iter(dis, matrix, conn)
	add_max_vert(vert_on_crystal, max_vert)
	place_on_crystal(crystal, max_vert)

	iter_num = 0
	while len(vert_on_crystal) < len(matrix):
		iter_num += 1
		max_vert, conn = crystal_iter(dis, matrix, conn, vert_on_crystal, iter_num)
		add_max_vert(vert_on_crystal, max_vert)
		print(vert_on_crystal)
		place_on_crystal(crystal, max_vert)

	return crystal, vert_on_crystal


def find_xdist(crystal, v1, v2):
	pos1 = np.argwhere(crystal == v1)[0]
	pos2 = np.argwhere(crystal == v2)[0]
	return abs(pos1[1] - pos2[1])

def find_ydist(crystal, v1, v2):
	pos1 = np.argwhere(crystal == v1)[0]
	pos2 = np.argwhere(crystal == v2)[0]
	return abs(pos1[0] - pos2[0])

def crys_dist(crystal, v1, v2):
	return find_xdist(crystal, v1, v2) + find_ydist(crystal, v1, v2)

def count_Q(crystal, matrix):
	Q = 0
	print(matrix)
	print(crystal)
	for v in range(len(matrix)):
		pos = np.argwhere(crystal == v)[0]
		for i in range(crystal.shape[0]):
			lij = 0
			Sij = 0
			for j in range(crystal.shape[1]):
				#lij = abs(pos[0] - i) + abs(pos[1] - j)
				lij = crys_dist(crystal, v, crystal[i][j])
				Sij = matrix[int(v)][int(crystal[i][j])]
				Q += lij * Sij
	Q /= 2
	print(Q)

def Vi_to_replace(dis, matrix, crystal):
	L = [0.0 for _ in range(len(matrix))]
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			L[i] += matrix[i][j] * crys_dist(crystal, i, j)
		L[i] /= dis[i]
	print("Lv")
	print(L)
	vi = L.index(max(L))
	return vi

def find_center_mass(dis, matrix, crystal, v):
	CM = [0.0, 0.0]
	pos_v = np.argwhere(crystal == v)[0]
	print(pos_v)
	input()
	for i in range(crystal.shape[0]):
		for j in range(crystal.shape[1]):
			if j == pos_v[1]:
				continue
			print("TEST ", v, crystal[i][j], \
						 matrix[v][crystal[i][j]], \
						 find_xdist(crystal, v, crystal[i][j]), \
						 matrix[v][crystal[i][j]] * find_xdist(crystal, v, crystal[i][j])
				 )
			input()
			CM[0] += matrix[v][crystal[i][j]] * find_xdist(crystal, v, crystal[i][j])
			CM[1] += matrix[v][crystal[i][j]] * find_ydist(crystal, v, crystal[i][j])
	CM /= dis[v]
	print(CM)
	input()
	return CM

def find_verts_around_cord(crystal, cord):
	xs, xe = 0, 0
	ys, ye = 0, 0
	if cord[0] < 1:
		xs = 0
		xe = 1
	elif cord[0] > crystal.shape[0]:
		xs = crystal.shape[0] - 1
		xe = crystal.shape[0]
	elif cord[0] - int(cord[0]) == 0:
		xs = int(cord[0]) - 1
		xe = int(cord[0])
	else:
		xs = int(cord[0]) - 1
		xe = int(cord[0]) + 1
	if cord[1] < 1:
		ys = crystal.shape[1] - 1
		ye = crystal.shape[1]
	elif cord[1] > crystal.shape[1]:
		ys = 0
		ye = 1
	elif cord[1] - int(cord[1]) == 0:
		ys = crystal.shape[1] - int(cord[1])
		ye = crystal.shape[1] - int(cord[1]) + 1
	else:
		ys = crystal.shape[1] - (int(cord[1]) + 1)
		ye = crystal.shape[1] - int(cord[1]) + 1
	return crystal[ys:ye, xs:xe]

def optimize_part(dis, matrix, crystal, vert_on_crystal):
	pass

def main():
	df = pd.read_csv("lec.csv", header=None)
	matrix = df.as_matrix()
	dis = np.array([0.0 for _ in range(len(matrix))])
	conn = np.array([[0.0 for _ in range(len(matrix))]])
	crystal = np.array([[None for _ in range(3)] for _ in range(3)])
	vert_on_crystal = []
	crystal, vert_on_crystal = prepare_part(dis, matrix, conn, crystal, vert_on_crystal)
	print(crystal)
	crystal[2][0] = 1
	crystal[1][2] = 5
	count_Q(crystal, matrix)
	vi = Vi_to_replace(dis, matrix, crystal)
	print("Vi", vi)
	#CM = find_center_mass(dis, matrix, crystal, vi)
	print(find_verts_around_cord(crystal, [3.5, 3.5]))
	input()
	#print("CM", CM)
	#input()

	#crystal[0][1] = 3
	#crystal[1][0] = vi
	#count_Q(crystal, matrix)

	print(crystal)
	print(vert_on_crystal)

main()
