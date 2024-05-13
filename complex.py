####################################################################################
#
#                DelaunayComplex
#
####################################################################################
#TODO:
#
#Input: np array of coordinate : [[1.1, 2.2],[2.3, 3.5]] 2D for now
#
#Output: 绘制，边的坐标，边的index
#
#Example:
#
# complex = DelaunayComplex(np.array(data),r) # r 代表complex中的边长最大值
# edge_list_position = complex.edge()#返回边上二点的坐标
# complex.edge_index() #返回边，以点对的形式，eg：[[1,2]] 代表点1和2间存在边
# complex.draw_complex() #绘制
####################################################################################
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

class DelaunayComplex:

    def __init__(self, data, max_r) -> None:
        # lets data X a set finite in R^d
        self.data = data
        #print(self.data)
        self.simplices = Delaunay(self.data).simplices
        self.max_r = max_r

    def mesh(self):
      return self.simplices



    def edge_index(self):
      def distance_btp(p1, p2):
        r0 = pow((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.5)
        return r0

      #simplices = Delaunay(self.data).simplices
      a = self.simplices
      #print(a)
      n = len(self.data.tolist())

      edge_matrix = zeros_array = np.zeros((n, n))

      for triangle in self.simplices:
        edge_matrix[triangle[0], triangle[1]] = 1
        edge_matrix[triangle[0], triangle[2]] = 1
        edge_matrix[triangle[1], triangle[2]] = 1
        edge_matrix[triangle[1], triangle[0]] = 1
        edge_matrix[triangle[2], triangle[0]] = 1
        edge_matrix[triangle[2], triangle[1]] = 1

      #print(edge_matrix)
      edge_list_index = []

      for  i in range(n):
        for j in range(i):
          if edge_matrix[i, j] == 1:
            #print(self.data[i],self.data[j])
            if(distance_btp(self.data[i],self.data[j]) < self.max_r):
              edge_list_index.append([i,j])


      return edge_list_index

    def two_simplex(self):
      edge_matrix= zeros_array = np.zeros((len(self.data),len(self.data)))
      for edge in self.edge_index():
        edge_matrix[edge[0],edge[1]] = 1
        edge_matrix[edge[1],edge[0]] = 1
      def find_triangles(adjacency_matrix):
        num_vertices = len(adjacency_matrix)
        triangles = []
        for i in range(num_vertices):
          for j in range(i + 1, num_vertices):
            for k in range(j + 1, num_vertices):
                # 检查三个顶点是否形成三角形
                if adjacency_matrix[i, j] and adjacency_matrix[j, k] and adjacency_matrix[k, i]:
                    triangle = [i, j, k]
                    triangles.append(triangle)
        return triangles
      two_simplex_list = find_triangles(edge_matrix)
      return two_simplex_list



    def draw_complex(self):
      for edge in self.edge_index():
        plt.plot(self.data[edge, 0], self.data[edge, 1],color= 'blue')
      #plt.triplot(self.data[:, 0], self.data[:, 1], self.simplices)
      plt.scatter(self.data[:, 0], self.data[:, 1], color='r')
      plt.show()
####################################################################################
#
#               HyperCube Cover
#
####################################################################################
#Input: 点坐标和边(with length),
#Output: 生成一组cover，X0,...,Xn
#    cover格式契合BlowupComplex的输入格式，[[point1], [point2], ..., [pointi, pointj],...]，可直接作为BlowupComplex的输入
#
####################################################################################
import random
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt



class HyperCube:

    def __init__(self, point, edge, two_simplex, N, M, r, add_two_simplex) -> None:
      self.point = point #一组带坐标的点,numpy array
      self.edge = edge #edge_index,调用端点坐标时使用 self.point[self.edge[i][0]], self.point[self.edge[i][1]]
      self.N = N #x方向cover个数
      self.M = M #y方向cover个数
      self.r = r #0-2,overlapping degree
      self.edge_matrix= zeros_array = np.zeros((len(self.point),len(self.point)))
      self.bounder_list_x = []
      self.bounder_list_y = []
      self.two_simplex = two_simplex
      self.add_two_simplex = add_two_simplex
      for edge in self.edge:
        self.edge_matrix[edge[0],edge[1]] = 1
        self.edge_matrix[edge[1],edge[0]] = 1

    def distance_btp(p1, p2):
      r0 = pow((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.5)
      return r0

    def cover(self):
      cover = []

      max_x = max(self.point[:,0])
      min_x = min(self.point[:,0])
      max_y = max(self.point[:,1])
      min_y = min(self.point[:,1])
      range_x = (max_x - min_x)/self.N
      range_y = (max_y - min_y)/self.M

      cover_x = []
      #print(type(self.N))
      for i in range(self.N):
        cover_x.append([])
        self.bounder_list_x.append([min_x + i * range_x - range_x*self.r/2,
                                    min_x + (i+1)*range_x + range_x*self.r/2])
        #print(min_x + (i+1)*range_x + range_x*self.r , min_x + i * range_x)
        for j, point in enumerate(self.point):

          if ((min_x + (i+1)*range_x + range_x*self.r/2) >= point[0] >= (min_x + i*range_x-range_x*self.r/2)):
            cover_x[i].append([j])

      cover_y = []
      for i in range(self.M):
        cover_y.append([])
        self.bounder_list_y.append([min_y + i * range_y - range_y*self.r/2,
                                    min_y + (i+1)*range_y + range_y*self.r/2])
        #print(min_y + (i+1)*range_y + range_y*self.r , min_y + i * range_y)
        for j, point in enumerate(self.point):
          if ((min_y + (i+1)*range_y + range_y*self.r/2) >= point[1] >= (min_y + i * range_y - range_y*self.r/2)):
            cover_y[i].append([j])

      cover = [[] for _ in range(self.N * self.M)]
      i = 0


      for x in range(self.N):
        for y in range(self.M):
          for point in cover_x[x]:
            if point in cover_y[y]:
              cover[i].append(point)

          i = i+1


      step_cover = cover
      for j in range(len(step_cover)):

        step_edge_matrix = self.edge_matrix

        for k in range(len(step_cover[j])): #第J个cover 第k个点

          for i in range(len(self.point)):


            #print(i,j)
            if step_edge_matrix[step_cover[j][k][0],i] == 1 and step_edge_matrix[i,step_cover[j][k][0]] == 1:# and [i] in step_cover[j]:
              if step_cover[j][k][0] < i:
                if [step_cover[j][k][0],i] in cover[j]:
                  continue
                else:
                  cover[j].append([step_cover[j][k][0],i])
                  if [i] in step_cover[j]:
                    continue
                  else:
                    cover[j].append([i])

              else:
                if [i, step_cover[j][k][0]] in cover[j]:
                  continue
                else:
                  cover[j].append([i,step_cover[j][k][0]])
                  if [i] in step_cover[j]:
                    continue
                  else:
                    cover[j].append([i])



      triangles_list = self.two_simplex

      step_cover = [[] for _ in range(self.N * self.M)]

      for triangle in triangles_list:
        j = -1
        for small_cover in cover:
          j = j+1

          if [triangle[0],triangle[1]] in small_cover and [triangle[0],triangle[2]] in small_cover and [triangle[1],triangle[2]] in small_cover:
            small_cover.append(triangle)
            continue
          if [triangle[0],triangle[1]] in small_cover or [triangle[0],triangle[2]] in small_cover or [triangle[1],triangle[2]] in small_cover:
            if [triangle[0],triangle[1]] not in small_cover and [triangle[0],triangle[1]] not in step_cover[j]:
              step_cover[j].append([triangle[0],triangle[1]])
            if [triangle[0],triangle[2]] not in small_cover and [triangle[0],triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[0],triangle[2]])
            if [triangle[1],triangle[2]] not in small_cover and [triangle[1],triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[1],triangle[2]])
            if [triangle[0]] not in small_cover and [triangle[0]] not in step_cover[j]:
              step_cover[j].append([triangle[0]])
            if [triangle[1]] not in small_cover and [triangle[1]] not in step_cover[j]:
              step_cover[j].append([triangle[1]])
            if [triangle[2]] not in small_cover and [triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[2]])
            step_cover[j].append(triangle)

      #print(step_cover)
      if (self.add_two_simplex):
        for j in range(len(step_cover)):
          for i in range(len(step_cover[j])):
            cover[j].append(step_cover[j][i])

      return cover

    def draw_cover(self):


        #print(self.bounder_list_x)
        #print(self.bounder_list_y)

        for i in range(self.N):
          for j in range(self.M):
            x1, y1 = self.bounder_list_x[i][0], self.bounder_list_y[j][0]
            x2, y2 = self.bounder_list_x[i][1], self.bounder_list_y[j][1]
            plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 'ro-',alpha = 0.5,
                     color = [ round(i/self.N, 1),round(j/self.M, 1),round(i*j/(self.N*self.M),1)])

        #print(self.edge)
        for edge in self.edge:

          plt.plot(self.point[edge, 0], self.point[edge, 1],color= 'blue')
          #plt.triplot(self.data[:, 0], self.data[:, 1], self.simplices)
        plt.scatter(self.point[:, 0], self.point[:, 1], color='r')
        plt.show()

        return 0



