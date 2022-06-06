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

def orient_test(p , he):
    mat = [[he.origin.x, he.origin.y, 1], [he.dest.x, he.dest.y, 1], [p.x, p.y, 1]]
    return np.linalg.det(mat) > 0   #return true if left turn

def orient_test_vertices(v1,v2,v3):
    mat = [[v1.x, v1.y, 1], [v2.x, v2.y, 1], [v3.x, v3.y, 1]]
    return np.linalg.det(mat) > 0   #return true if left turn

def orient_test_points(p1,p2,p3):
    mat = [[p1[0], p1[1], 1], [p2[0], p2[1], 1], [p3[0], p3[1], 1]]
    return np.linalg.det(mat) > 0   #return true if left turn


def square_s(a,b):
    return a**2 + b**2

def incircle_test_h(p, he1, he2, he3):  #incircle test using half edges as input
    a = he1.v1
    b = he2.v1
    c = he3.v1
    d = p

    mat = [[a[0], a[1], square_s(a[0], a[1]),1],
           [b[0], b[1], square_s(b[0], b[1]),1],
           [c[0], c[1], square_s(c[0], c[1]),1],
           [d[0], d[1], square_s(d[0], d[1]),1]]

    return np.linalg.det(mat) > 0   #return True if inside the triangle

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
    
def compute_angle(v1,v2):
    res = math.atan2(v2.x - v1.x, v2.y - v1.y)
    if res < 0:
        res = 360*pi/180 + res
    return res          

class Vertex:

    def __init__(self,p) :
        
        self.x = p[0]
        self.y = p[1]
        self.hlist = []     #list of hedges coming from v

        #convex hull info
        self.prev = None
        self.next = None

    def sort_inc(self):     #sort incident half edges
        self.hlist.sort( key=lambda hedge: hedge.angle)

    def dist(self, other_v):
        return np.sqrt((self.x - other_v.x)**2 + (self.y - other_v.y)**2)

class Face:

    def __init__(self):
        self.hedge = None
        self.name = None
        self.voronoi = False
        self.center = None

    def inside(self,p):     #checks whether a pt is inside the face
        h = self.hedge
        if orient_test(p,h):
            while ( h.nexthe != self.hedge):
                h = h.nexthe
                if not orient_test(p,h):    #check for every h edge to have the point on the left
                    return False
            return True
        else:
            return False

#declaration of external face
f_ext = Face()
f_ext.name = 'ext_f'
    
class Hedge:

    def __init__(self,v1,v2, f = f_ext):

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
    def __init__(self, f = f_ext):
        self.f_ext = f
        self.points = []
        self.vertices = []
        self.hedges = []
        self.voronoi_hedges = []
        self.faces = [self.f_ext]
        self.ext_dic = {}

    def check_convexity(self):
        for p in range(len(self.points)-2):
            p1 = self.points[p]
            p2 = self.points[p+1]
            p3 = self.points[p+2]
            if not orient_test_points(p1,p2,p3):
                print('input points not in convex position\n') 
                exit()
        print('input points in convex position\n')
        return 0

    def face_already_all(self,f1):
        for f in self.faces:
            if f.name == f1.name :
                return f
        return False

    def add_face(self,f):
        self.faces.append(f)
        f.name = len(self.faces)-1

    def remove_face(self,f1):
        if not self.face_already_all(f1):
            print('face not found\n')
            return 0
        else:
            f1 = self.face_already_all(f1)
            self.faces.remove(f1)
            #update of all faces names
            for f in range(1,len(self.faces)):
                self.faces[f].name = f


    def intersect(self,v1,v2):  #return true if intersection with any of the existinf half edges (not used)
        for h in self.hedges:
            #check non coincident coordinates
            if ((v1.x != h.origin.x and v1.y != h.origin.y) or (v1.x != h.dest.x and v1.y != v1.y != h.dest.y) or (v2.x != h.origin.x and v2.y != h.origin.y) or (v2.x != h.origin.x and v2.y != h.origin.y)):
                pass
            elif orient_test_vertices(h.origin,h.dest,v1) != orient_test_vertices(h.origin,h.dest,v2) and orient_test_vertices(v1,v2,h.origin) != orient_test_vertices(v1,v2,h.dest):
                return True
        return False


    def find_hedge(self,v):     #find internal half edge opposite to v
        for h in self.hedges:
            if v.next == h.origin and v.prev == h.dest:
                return h
        print('hedge not found!\n')
        return 0

    def vertex_already(self,vi):    #checks if the vertex is already stored
        for v in self.vertices:
            if v.x == vi.x and v.y == vi.y:
                return True
        return False

    def hedge_already(self,h1):     #checks if the hedge is already stored
        for h in self.hedges:
            if h.origin == h1.origin and h.dest == h1.dest:
                return h
        return False

    def read_points(self):      #store points from input txt file
        for l in sys.stdin: 
            L = l.strip().split(' ') #put input names of a line into a vector
            x = L[0]
            y = L[1]
            self.points.append((float(x),float(y)))
        if len(self.points) < 3:
            print('not enough points to triangulate\n')
            exit()
            

    def face_res(self,v):     #returns the face containing the vertex
        for f in self.faces:
            if f.inside(v):
                return f
            else:
                return False

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
                self.hedges.append(h)
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
        if not self.hedge_already(h1):
            self.hedges.append(h1)
            v1.hlist.append(h1)
        else:
            h1 = self.hedge_already(h1)
            v1.hlist.append(h1)
            
        if not self.hedge_already(h2) :
            self.hedges.append(h2)
            v2.hlist.append(h2)
        else:
            h = self.hedge_already(h2)
            v2.hlist.append(h)

        if not self.hedge_already(h3) :
            self.hedges.append(h3)
            v1.hlist.append(h3)
        else:
            h = self.hedge_already(h3)
            v1.hlist.append(h)

        if not self.hedge_already(h4)  :
            self.hedges.append(h4)
            v3.hlist.append(h4)
        else:
            h4 = self.hedge_already(h4)
            v3.hlist.append(h4)

        if not self.hedge_already(h5) :
            self.hedges.append(h5)
            v2.hlist.append(h5)
        else:
            h5 = self.hedge_already(h5)
            v2.hlist.append(h5)        

        if not self.hedge_already(h6) :
            self.hedges.append(h6)
            v3.hlist.append(h6)
        else:
            h = self.hedge_already(h6)
            v3.hlist.append(h)

        vs = [v1,v2,v3]
        
        for v in vs:
            v.sort_inc()

        #add the new face
        f = Face()
        f.hedge = h5
        # f.name = len(self.faces)
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
            self.hedges.remove(h_extra)
            self.hedges.remove(h_extra.twin)
            
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
            if not self.hedge_already(h1):
                self.hedges.append(h1)
            else:
                h1 = self.hedge_already(h1)
                
            if not self.hedge_already(h1t):
                self.hedges.append(h1t)
            else:
                h1t = self.hedge_already(h1t)

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

        #add all vertices to the DCEL
        for p in self.points:
            v = Vertex(p)
            self.vertices.append(v)

        #build vertex neighbours
        self.vertices[0].prev = self.vertices[-1]
        self.vertices[0].next = self.vertices[1]
        self.vertices[-1].next = self.vertices[0]
        self.vertices[-1].prev = self.vertices[-2]

        for p in range(1,len(self.points)-1):
            self.vertices[p].prev = self.vertices[p-1]
            self.vertices[p].next = self.vertices[p+1]

        #create external face half edge
        f_eh = Hedge(self.vertices[1],self.vertices[0])
        self.f_ext.hedge = f_eh
        
        #create external half edges of the hull
        h_0 = Hedge(self.vertices[0],self.vertices[-1])
        h_0.iface = self.f_ext
        self.hedges.append(h_0)

        for v in range(len(self.vertices)-1):
            h = Hedge(self.vertices[v+1],self.vertices[v])
            h.iface = self.f_ext
            self.hedges.append(h)
        
        #connect the half edges
        self.hedges[0].nexthe = self.hedges[-1]
        self.hedges[0].prevhe = self.hedges[1]
        for h in range(1,len(self.hedges)-1):
            self.hedges[h].prevhe = self.hedges[h+1]
            self.hedges[h].nexthe = self.hedges[h-1]
        self.hedges[len(self.hedges)-1].nexthe = self.hedges[-2]
        self.hedges[len(self.hedges)-1].prevhe = self.hedges[0]

        tmp_vertices = self.vertices.copy()

        #incrementally pluck a vertex and push it to the stack
        while len(tmp_vertices)>2:
            pr = tmp_vertices[random.randint(0, len(tmp_vertices)-1)]
            tmp_vertices.remove(pr)
            #push pr to the stack
            S.append(pr)
            #link the 2 new neighbour vertices
            pr.prev.next = pr.next
            pr.next.prev = pr.prev

        #start building triangles once the stack is full
        v1 = stack_pop(S)
        v2 = v1.next
        v3 = v1.prev
        self.add_triangle(v1,v2,v3)
        
        while len(S)>0:     #keep building triangles until the stack is empty
            v = stack_pop(S)
            self.new_triangle(v)
            self.flip_test(v)

        #fix ext hedges connections
        for h in self.hedges:
            if h.iface == f_ext:    #build the  dictionary
                self.ext_dic[(h.origin.x,h.origin.y,h.dest.x,h.dest.y)] = h
        #connect 
        for h in self.hedges:
            if h.twin!= None:
                key = (h.twin.origin.x,h.twin.origin.y,h.twin.dest.x,h.twin.dest.y)
                if key in self.ext_dic:
                    h_temp = self.ext_dic[(h.twin.origin.x,h.twin.origin.y,h.twin.dest.x,h.twin.dest.y)]  #get the ext hedge
                    h.twin = h_temp
                    h_temp.twin = h

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
        return(Vertex((x_m,y_m)))    

    def voronoi_func(self):
        #compute voronoi centers
        for h in self.hedges:
            if h.iface.name != 'ext_f' and not h.iface.voronoi:
                h.iface.voronoi = True  #mark as computed
                v1 = h.origin
                v2 = h.dest
                v3 = h.prevhe.origin
                cx,cy = self.compute_center(v1,v2,v3)
                center = Vertex((cx,cy))
                h.iface.center = center
        
        #connect neighbours centers
        for h in self.hedges:
            if h.iface.name == 'ext_f':
                m = self.compute_midpoint(h.origin,h.dest)
                vh = Hedge(h.twin.iface.center,m)   #connect voronoi center with midpoint of corresp h edge
                self.voronoi_hedges.append(vh)

            else:
                if h.iface.center!= None and h.twin.iface.center!= None:
                    p1 = h.iface.center
                    p2 = h.twin.iface.center
                    vh = Hedge(p1,p2)
                    self.voronoi_hedges.append(vh)


#main execution
myDCEL = DCEL()
myDCEL.read_points()
myDCEL.check_convexity()
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
X= []
Y = []
#arrays to plot voronoi
Xv = []
Yv = []

#fill the arrays and print DCEL infos
for h in myDCEL.voronoi_hedges:
    Xv.append(h.origin.x)
    Yv.append(h.origin.y)
  
    Xv.append(h.dest.x)
    Yv.append(h.dest.y)

print('final half edges, in toatal:',len(myDCEL.hedges),'\n')
for h in myDCEL.hedges:
    print('hedge:',h.origin.x,h.origin.y,'->',h.dest.x,h.dest.y)
    if h.twin!= None:
        print('twin hedge', h.twin.origin.x,h.twin.origin.y,'->',h.twin.dest.x,h.twin.dest.y)
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
    #plot delauny
    for i in range(0,len(X)-2,2):
        plt.plot(X[i:i+2], Y[i:i+2], 'g', label='graph')

if flag_v:
    #plot voroni
    for i in range(0,len(Xv)-2,2):
        plt.plot(Xv[i:i+2], Yv[i:i+2], 'r', label='graph')

#write an output file with the new order of the vertices
with open('output.txt', 'w') as f:
    for h in myDCEL.hedges:
        f.write(f"{h.origin.x} {h.origin.y}\n")
        f.write(f"{h.dest.x} {h.dest.y}\n")

#plot voronoi centers
X_v = []
Y_v = []
for f in myDCEL.faces:
    if f.center!= None:
        X_v.append(f.center.x)
        Y_v.append(f.center.y)
        print('voronoi center at:',f.center.x,f.center.y,'for face ',f.name,'\n')

plt.scatter(X_v, Y_v)
plt.show()





    
        