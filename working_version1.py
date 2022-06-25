from cmath import pi
import math
import sys
import numpy as np
import random
import matplotlib.pyplot as plt

def square_s(a,b):
    return a**2 + b**2

def stack_pop(S):
    x = S[-1]
    S.pop()
    return x

def compute_angle(v1,v2):
    res = math.atan2(v2.x - v1.x, v2.y - v1.y)
    if res < 0:
        res = 360*pi/180 + res
    return res

def incircle_test(p, v1, v2, v3):   #incircle test using vertices as input
    a = v1
    b = v2
    c = v3
    d = p

    mat = [[a.x, a.y, square_s(a.x, a.y),1],
           [b.x, b.y, square_s(b.x, b.y),1],
           [c.x, c.y, square_s(c.x, c.y),1],
           [d.x, d.y, square_s(d.x, d.y),1]]

    return np.linalg.det(mat) > 0   #return True if inside the triangle

def orient_test_vertices(v1,v2,v3):
    mat = [[v1.x, v1.y, 1], [v2.x, v2.y, 1], [v3.x, v3.y, 1]]
    return np.linalg.det(mat) > 0   #return true if left turn

class Vertex:

    def __init__(self,x,y) :
        
        self.x = x
        self.y = y
        self.hlist = []     #list of hedges coming from v
        #to keep track in the randomized construction
        self.prev = None
        self.next = None

    def sort_inc(self):     #sort incident half edges
        self.hlist.sort( key=lambda hedge: hedge.angle)

class Face:

    def __init__(self):
        self.hedge = None
        self.name = None
        self.voronoi = False
        self.center = None

class Hedge:

    def __init__(self,v1,v2):

        self.origin = v1
        self.dest = v2
        self.twin = None
        self.iface = None      
        self.nexthe = None
        self.prevhe = None
        self.len = math.sqrt((v2.x - v1.x)**2 + (v2.y - v1.y)**2)
        self.angle = compute_angle(v1,v2)

class DCEL:
    def __init__(self):     
        self.vertices = []
        self.voronoi_hedges = []
        self.faces_dic = {}
        self.hedges_dic = {}
        self.counter = 0

    def add_face(self,f):
        f.name = self.counter
        self.faces_dic[f.name] = f
        self.counter+=1

    def add_hedge(self,h):
        key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
        if key in self.hedges_dic:
            #find existing one
            hedge = self.hedges_dic[key]
        else:
            #create new one
            hedge = h 
            self.hedges_dic[key] = hedge
            #add half edge to the vertex h edge list
            h.origin.hlist.append(hedge)
            
        #twin setting
        key_t = (h.dest.x,h.dest.y,h.origin.x,h.origin.y)
        #twin exists-> link it
        if key_t in self.hedges_dic:
            twin = self.hedges_dic[key_t]  
        else:
            #create twin and link it
            twin = Hedge(Vertex(h.dest.x,h.dest.y),Vertex(h.origin.x,h.origin.y))
            self.hedges_dic[key_t] = twin
            #add half edge to the vertex h edge list
            h.dest.hlist.append(twin)
        hedge.twin = twin
        twin.twin = hedge
        return hedge,twin

    def remove_hedge(self,h):
        key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
        del self.hedges_dic[key]

    def remove_face(self,f1):
        key = f1.name
        if key not in self.faces_dic:
            print('face not found\n')
            return 0
        else:
            del self.faces_dic[f1.name]


    def read_vertices(self):      #store points from input txt file
        for l in sys.stdin: 
            L = l.strip().split(' ') #put input names of a line into a vector
            x = L[0]
            y = L[1]
            self.vertices.append(Vertex(float(x),float(y)))
        if len(self.vertices) < 3:
            print('not enough points to triangulate\n')
            exit()
        #check convexity
        for v in range(len(self.vertices)-2):
            v1 = self.vertices[v]
            v2 = self.vertices[v+1]
            v3 = self.vertices[v+2]
            if not orient_test_vertices(v1,v2,v3):
                print('input vertices not in convex position\n') 
                exit()
        print('input points in convex position\n')
        return 0

    def flip_test(self, v, h_opposite):
        #find the opposite vertex to compute the inircle test 
        if h_opposite.twin.iface!= None:
            v_opposite = h_opposite.twin.nexthe.dest
            #perform the incircle test
            if incircle_test(v_opposite, v, h_opposite.dest, h_opposite.origin):
                #flip
                self.recursive_flip(v,h_opposite,v_opposite)
            else: return False
        else: return False

    def recursive_flip(self, v, h_opposite,v_opposite):
        #flip 
        if h_opposite.iface!= None:
            self.remove_face(h_opposite.iface)
        if h_opposite.twin.iface!= None:
            self.remove_face(h_opposite.twin.iface)
        self.remove_hedge(h_opposite)
        self.remove_hedge(h_opposite.twin)

        new_h,new_h_t = self.add_hedge(Hedge(v,v_opposite))
        #link hedges
        #v->v_opposite
        h_next,_ = self.add_hedge(Hedge(v_opposite,h_opposite.dest))
        h_prev,_ = self.add_hedge(Hedge(h_opposite.dest,v))
        #set
        new_h.prevhe = h_prev
        new_h.nexthe = h_next

        h_prev.nexthe = new_h
        h_prev.prevhe = h_next

        h_next.prevhe = new_h
        h_next.nexthe = h_prev

        #v_opposite->v
        h_next_t,_ = self.add_hedge(Hedge(v,h_opposite.origin))
        h_prev_t,_ = self.add_hedge(Hedge(h_opposite.origin,v_opposite))
        #set
        new_h_t.prevhe = h_prev_t
        new_h_t.nexthe = h_next_t

        h_prev_t.prevhe = h_next_t
        h_prev_t.nexthe = new_h_t

        h_next_t.prevhe = new_h_t
        h_next_t.nexthe = h_prev_t
        #faces setting
        #face1
        f1 = Face()
        f1.hedge = new_h
        self.add_face(f1)
        new_h.iface = f1
        new_h.prevhe.iface = f1
        new_h.nexthe.iface = f1
              
        #face2 
        f2 = Face()
        f2.hedge = new_h_t
        self.add_face(f2)
        new_h_t.iface = f2
        new_h_t.prevhe.iface = f2
        new_h_t.nexthe.iface = f2
        #recursion call
        h_opposite_1,_ =  self.add_hedge(new_h.nexthe)
        h_opposite_2,_ = self.add_hedge(new_h_t.prevhe)
        self.flip_test(v, h_opposite_1)
        self.flip_test(v, h_opposite_2)

    def new_triangle(self,v):   #add a triangle connecting a input vertex to its conv hull neighbours

        hn,hnt = self.add_hedge(Hedge(v.next,v))
        hp,hpt = self.add_hedge(Hedge(v,v.prev))

        h_opposite,h_opposite_t = self.add_hedge(Hedge(v.prev,v.next))   #internal edge, opposite to v
        # v1_opposite = h_opposite.dest
        # v2_opposite = h_opposite.origin
        #setting prev and next he
        hn.nexthe = hp
        hn.prevhe = h_opposite

        hp.nexthe = h_opposite
        hp.prevhe = hn

        h_opposite.nexthe = hn
        h_opposite.prevhe = hp
        #add the new face
        f = Face()
        f.hedge = h_opposite
        # f.name = len(self.faces)
        self.add_face(f)

        #associate internal face
        hn.iface = f
        hp.iface = f
        h_opposite.iface = f
        self.flip_test(v,h_opposite)

    def incremental_delauny(self):
        #stack 
        S = []
        #connect neighbours vertices
        for p in range(-1,len(self.vertices)-1):
            self.vertices[p].prev = self.vertices[p-1]
            self.vertices[p].next = self.vertices[p+1]
        tmp_vertices = self.vertices.copy()
        #incrementally pluck a vertex and push it to the stack
        while len(tmp_vertices)>2:
            vr = tmp_vertices[random.randint(0, len(tmp_vertices)-1)]
            tmp_vertices.remove(vr)
            #push pr to the stack
            S.append(vr)
            #link the 2 new neighbour vertices
            vr.prev.next = vr.next
            vr.next.prev = vr.prev
        while len(S)>0:     #keep building triangles until the stack is empty
            v = stack_pop(S)
            self.new_triangle(v)

    def compute_center(self,v1,v2,v3):
        
            p1 = (v1.x,v1.y)
            p2 = (v2.x,v2.y)
            p3 = (v3.x,v3.y)
        
            temp = p2[0] * p2[0] + p2[1] * p2[1]
            bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
            cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
            det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

            if abs(det) < 1.0e-6:
                return (None, np.inf)

            # Center of circle
            cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
            cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det

            return cx,cy

    def compute_midpoint(self,v1,v2):
        x_m = (v1.x + v2.x)/2
        y_m = (v1.y + v2.y)/2
        return(Vertex(x_m,y_m)) 

    #main function to compute voronoi diagram
    def voronoi_func(self):
        #compute voronoi centers
        for key in self.hedges_dic:
            h=self.hedges_dic[key]
            if h.iface != None and not h.iface.voronoi:
                h.iface.voronoi = True  #mark as computed
                v1 = h.origin
                v2 = h.dest
                v3 = h.prevhe.origin
                cx,cy = self.compute_center(v1,v2,v3)
                center = Vertex(cx,cy)
                h.iface.center = center 
        #connect neighbours centers
        for key in self.hedges_dic:
            h = self.hedges_dic[key]
            if h.iface == None:
                m = self.compute_midpoint(h.origin,h.dest)
                vh = Hedge(h.twin.iface.center,m)   #connect voronoi center with midpoint of corresp h edge
                self.voronoi_hedges.append(vh)
            elif h.iface.center!= None and h.twin.iface!= None:
                p1 = h.iface.center
                p2 = h.twin.iface.center
                vh = Hedge(p1,p2)
                self.voronoi_hedges.append(vh)

    def plots(self,input='both'):
        #arrays for plotting
        #delaunat h edges
        X = []
        Y = []
        #voronoi h edges
        Xv = []
        Yv = []
        #voronoi center
        X_v = []
        Y_v = []
        #fill the arrays and print DCEL infos
        #half edges
        for key in self.hedges_dic:
            h = self.hedges_dic[key]
            X.append(h.origin.x)
            Y.append(h.origin.y)
            X.append(h.dest.x)
            Y.append(h.dest.y)

        #voronoi half edges
        for h in self.voronoi_hedges:
            Xv.append(h.origin.x)
            Yv.append(h.origin.y)
            Xv.append(h.dest.x)
            Yv.append(h.dest.y)

        #voronoi centers
        for key in self.faces_dic:
            f = self.faces_dic[key]
            if f.center!= None:
                X_v.append(f.center.x)
                Y_v.append(f.center.y)

        if input == 'delaunay':
            for i in range(0,len(X)-2,2):
                plt.plot(X[i:i+2], Y[i:i+2], 'g', label='graph')
        elif input == 'voronoi':
            for i in range(0,len(Xv),2):
                plt.plot(Xv[i:i+2], Yv[i:i+2], 'r', label='graph')
        elif input == 'both':
            for i in range(0,len(X)-2,2):
                plt.plot(X[i:i+2], Y[i:i+2], 'g', label='graph')
            for i in range(0,len(Xv),2):
                plt.plot(Xv[i:i+2], Yv[i:i+2], 'r', label='graph')
            plt.scatter(X_v, Y_v)
        plt.show()

#main execution
def main():
    myDCEL = DCEL()
    myDCEL.read_vertices()
    myDCEL.incremental_delauny()
    myDCEL.voronoi_func()

    print(len(myDCEL.faces_dic),'total faces')
    print(len(myDCEL.hedges_dic),'total half edges')

    #plot voronoi centers
    if len(sys.argv) < 2:
        #plot both
        myDCEL.plots()
    else:
        #plot the one specified by the input
        myDCEL.plots(sys.argv[1])

    #write an output file with the new order of the vertices
    with open('output.txt', 'w') as f:
        for key in myDCEL.hedges_dic:
            h = myDCEL.hedges_dic[key]
            f.write(f"{h.origin.x} {h.origin.y}\n")
            f.write(f"{h.dest.x} {h.dest.y}\n")

if __name__ == "__main__":
    main()


            