from cmath import inf, pi
import math
import sys
import numpy as np
import random
import matplotlib.pyplot as plt

def stack_pop(S):
    x = S[-1]
    S.pop()
    return x
def square_s(a,b):
    return a**2 + b**2

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

def compute_angle(v1,v2):
    res = math.atan2(v2.x - v1.x, v2.y - v1.y)
    if res < 0:
        res = 360*pi/180 + res
    return res 

class Vertex:

    def __init__(self,x,y) :
        
        self.x = x
        self.y = y
        self.hlist = []     #list of hedges coming from v

        #convex hull info
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


#main class: Doubly Connected Edge List           
class DCEL:
    def __init__(self):
        self.f_ext = None      
        self.vertices = []
        self.voronoi_hedges = []
        self.ext_hedges_dic = {}
        self.ext_hedges = []
        self.faces_dic = {}
        self.hedges_dic = {}
        self.counter = 0

    def add_face(self,f):
        f.name = self.counter
        self.faces_dic[f.name] = f
        self.counter+=1

    def add_hedge(self,h):
        key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
        self.hedges_dic[key] = h
    def add_ext_hedge(self,h):
        self.ext_hedges.append(h)
        key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
        self.ext_hedges_dic[key] = h
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


    def find_hedge(self,v):     #find internal half edge opposite to v
        key = (v.next.x,v.next.y,v.prev.x,v.prev.y)
        if key in self.hedges_dic:
            return self.hedges_dic[key]
        print('hedge not found!')
        return 0

    def hedge_already(self,h):     #checks if the hedge is already stored
        key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
        if key in self.hedges_dic:
            return self.hedges_dic[key]
        return False

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

    
    def new_triangle(self,v):   #add a triangle connecting a input vertex to its conv hull neighbours

        hn = Hedge(v,v.next)
        hnt= Hedge(v.next,v)
        hp = Hedge(v.prev,v)
        hpt = Hedge(v,v.prev)

        h_extra = self.find_hedge(v)    #internal edge, opposite to v
        v1_extra = h_extra.dest
        v2_extra = h_extra.origin
        
        #setting twins
        hn.twin = hnt
        hnt.twin = hn
        hp.twin = hpt
        hpt.twin = hp

        #setting prev and next he
        hn.nexthe = h_extra
        hn.prevhe = hp

        hp.nexthe = hn
        hp.prevhe = h_extra

        h_extra.nexthe = hp
        h_extra.prevhe = hn

        #setting incident h edges
        v.hlist.append(hn)
        v.hlist.append(hpt)

        v1_extra.hlist.append(hp)
        v2_extra.hlist.append(hnt)

        #sort the inc h edges
        vs = [v,v1_extra,v2_extra]
        for v in vs:
            v.sort_inc()

        hs = [hn,hnt,hp,hpt]

        #store the hedges
        for h in hs:
            if not self.hedge_already(h):
                self.add_hedge(h)
            else:
                h = self.hedge_already(h)
                
        #add the new face
        f = Face()
        f.hedge = h_extra
        # f.name = len(self.faces)
        self.add_face(f)

        #associate internal face
        hn.iface = f
        hp.iface = f
        h_extra.iface = f

    def add_triangle(self, v1,v2,v3):   #add a triangle from input vertices
        #vertices are given in ccw order

        h1 = Hedge(v1,v2)   #inside he v1->v2
        h2 = Hedge(v2,v1)   #ext twin
        h3 = Hedge(v1,v3)   #ext twin
        h4 = Hedge(v3,v1)   #inside he v3->v1
        h5 = Hedge(v2,v3)   #inside he v2->v3
        h6 = Hedge(v3,v2)   #ext twin
        
        #set next and prev h edges
        h1.nexthe = h5
        h1.prevhe = h4
        h4.nexthe = h1
        h4.prevhe = h5
        h5.nexthe = h4
        h5.prevhe = h1

        #set twins
        h1.twin = h2
        h2.twin = h1
        h3.twin = h4
        h4.twin = h3
        h5.twin = h6
        h6.twin = h5

        #store new half edges
        hs = [h1,h2,h3,h4,h5,h6]
        for h in hs:
            if not self.hedge_already(h):
                self.add_hedge(h)
                v1.hlist.append(h)
            else:
                h = self.hedge_already(h)
                h.origin.hlist.append(h)

        vs = [v1,v2,v3]
        
        for v in vs:
            v.sort_inc()

        #add the new face
        f = Face()
        f.hedge = h5
        self.add_face(f)

        #associate the face with the internal h edges
        h1.iface = f
        h4.iface = f
        h5.iface = f

        return f

    def flip_test(self,v):      #check if the 2 triangle need to be flipped
        h_extra =  self.find_hedge(v)    #internal edge, opposite to v
        v1 = v
        v2 = v.next
        v3 = v.prev

        p = h_extra.twin.nexthe.dest
        self.flip(p,v1,v2,v3)

    def flip(self,p,v1,v2,v3):
        
        if incircle_test(p,v1,v2,v3):
            h_extra = self.hedge_already(Hedge(v2,v3))
            self.remove_face(h_extra.iface)
            self.remove_face(h_extra.twin.iface)
            self.remove_hedge(h_extra)
            self.remove_hedge(h_extra.twin)
            
            h1 = Hedge(p,v1)
            h1t= Hedge(v1,p)
            
            #nexts and prevs
            h_a1= self.hedge_already(Hedge(v1,v2))
            h_a2 = self.hedge_already(Hedge(v2,p))
            h_a3 = self.hedge_already(Hedge(p,v3))
            h_a4 = self.hedge_already(Hedge(v3,v1))

            h1.twin = h1t
            h1t.twin = h1

            h1.nexthe = h_a1
            h_a1.prevhe = h1
            h_a1.nexthe = h_a2

            h1.prevhe = h_a2
            h_a2.nexthe = h1
            h_a2.prevhe = h_a1

            h1t.nexthe = h_a3
            h_a3.prevhe = h1t
            h_a3.nexthe = h_a4

            h1t.prevhe = h_a4
            h_a4.nexthe = h1t
            h_a4.prevhe = h_a3

            #store the new half edges
            for h in[h1,h1t]:
                if not self.hedge_already(h):
                    self.add_hedge(h)
                else:
                    h1 = self.hedge_already(h)

            #faces
            f1 = Face()
            f1.hedge = h1
            self.add_face(f1)

            #associate the face with the internal h edges
            h1.iface = f1
            h_a1.iface = f1
            h_a2.iface = f1

            f2 = Face()
            f2.hedge = h1t
            self.add_face(f2)

            h1t.iface = f2
            h_a3.iface = f2
            h_a4.iface = f2

            #recursevely perform flip test on the other half edges
            if h_a1.nexthe.twin.iface != None and h_a1.nexthe.twin.iface != 'ext_f':
                v1 = v1
                v2 = v2
                v3 = p
                p = h_a1.nexthe.twin.nexthe.dest
                self.flip(p,v1,v2,v3)
                
            if h_a4.prevhe.twin.iface != None and h_a4.prevhe.twin.iface != 'ext_f':
                v1 = v1
                v2 = p
                v3 = v3
                p = h_a4.prevhe.twin.nexthe.dest
                self.flip(p,v1,v2,v3)


    def incremental_delauny(self):
        #stack 
        S = []
        #connect neighbours vertices
        for p in range(-1,len(self.vertices)-1):
            self.vertices[p].prev = self.vertices[p-1]
            self.vertices[p].next = self.vertices[p+1]
        #create external half edges of the hull
        for v in range(-1,len(self.vertices)-1):
            h = Hedge(self.vertices[v+1],self.vertices[v])
            h.iface = self.f_ext
            self.add_ext_hedge(h)
        #connect the ext half hedges
        for h in range(-1,len(self.ext_hedges)-1):
            self.ext_hedges[h].prevhe = self.ext_hedges[h+1]
            self.ext_hedges[h].nexthe = self.ext_hedges[h-1]
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
        #start building triangles once the stack is full
        v1 = stack_pop(S)
        v2 = v1.next
        v3 = v1.prev
        self.add_triangle(v1,v2,v3)
        #add external face
        self.f_ext = Face()
        self.f_ext.hedge = self.ext_hedges[0]
        self.f_ext.name = 'f_ext'
        self.faces_dic[self.f_ext.name] = self.f_ext
        while len(S)>0:     #keep building triangles until the stack is empty
            v = stack_pop(S)
            self.new_triangle(v)
            self.flip_test(v)
            #connect 
        for key in self.hedges_dic:
            h = self.hedges_dic[key]
            if h.iface== None:
                key = (h.origin.x,h.origin.y,h.dest.x,h.dest.y)
                if key in self.ext_hedges_dic:
                    h_temp = self.ext_hedges_dic[(h.origin.x,h.origin.y,h.dest.x,h.dest.y)]  #get the ext hedge
                    self.add_hedge(h_temp)

                    h.twin.twin = h_temp
                    h_temp.twin = h.twin

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
#main execution
myDCEL = DCEL()
myDCEL.read_vertices()
flag_d = False
flag_v = False
if len(sys.argv) < 2:
    #defalut compute Delaunay and Voronoi
    flag_d = True
    flag_v = True
    myDCEL.incremental_delauny()
    myDCEL.voronoi_func()

elif sys.argv[1] == 'delaunay':
    #only compute delaunay
    flag_d = True
    myDCEL.incremental_delauny()

elif sys.argv[1] == 'voronoi':
    #only compute voronoi
    flag_v = True
    myDCEL.incremental_delauny()
    myDCEL.voronoi_func()

#arrays to plot delaunay
X = []
Y = []
#arrays to plot voronoi
Xv = []
Yv = []
print('you should have these voronoi edges',len(myDCEL.voronoi_hedges))
#fill the arrays and print DCEL infos
for h in myDCEL.voronoi_hedges:
    Xv.append(h.origin.x)
    Yv.append(h.origin.y)
  
    Xv.append(h.dest.x)
    Yv.append(h.dest.y)
print("faces and center")
for key in myDCEL.faces_dic:
    if myDCEL.faces_dic[key].center!= None:
        print(myDCEL.faces_dic[key].name,'->',myDCEL.faces_dic[key].center.x,myDCEL.faces_dic[key].center.y)

print('final half edges, in toatal:',len(myDCEL.hedges_dic),'\n')
print('final faces, in toatal:',len(myDCEL.faces_dic),'\n')
for key in myDCEL.hedges_dic:
    h = myDCEL.hedges_dic[key]
    print('hedge:',h.origin.x,h.origin.y,'->',h.dest.x,h.dest.y)
    if h.twin!= None:
        print('twin hedge', h.twin.origin.x,h.twin.origin.y,'->',h.twin.dest.x,h.twin.dest.y)
        if(h.twin.iface):
            print('twin hedge internal face center', h.twin.iface.center.x,h.twin.iface.center.y)
        else:
            print('no twin internal face')
    else:
        print('NO TWIN\n')
    
    X.append(h.origin.x)
    Y.append(h.origin.y)

    X.append(h.dest.x)
    Y.append(h.dest.y)
    print('next hedge:',h.nexthe.origin.x,h.nexthe.origin.y,'->',h.nexthe.dest.x,h.nexthe.dest.y)
    print('prev hedge', h.prevhe.origin.x,h.prevhe.origin.y,'->',h.prevhe.dest.x,h.prevhe.dest.y)
    if h.iface != None:
        print('face:',h.iface.name)
    else:
        print('no face')
    print('______','\n')
print('final vertices\n')            
for v in myDCEL.vertices:
    print(v.x,v.y,'\n')

#PLOTS
if flag_d:
#plot delauny edges
    for i in range(0,len(X)-2,2):
        plt.plot(X[i:i+2], Y[i:i+2], 'g', label='graph')
if flag_v:
    #plot voroni hedges
    for i in range(0,len(Xv),2):
        plt.plot(Xv[i:i+2], Yv[i:i+2], 'r', label='graph')

#write an output file with the new order of the vertices
with open('output.txt', 'w') as f:
    for key in myDCEL.hedges_dic:
        h = myDCEL.hedges_dic[key]
        f.write(f"{h.origin.x} {h.origin.y}\n")
        f.write(f"{h.dest.x} {h.dest.y}\n")

#plot voronoi centers
X_v = []
Y_v = []
for key in myDCEL.faces_dic:
    f = myDCEL.faces_dic[key]
    if f.center!= None:
        X_v.append(f.center.x)
        Y_v.append(f.center.y)
plt.scatter(X_v, Y_v)
plt.show()