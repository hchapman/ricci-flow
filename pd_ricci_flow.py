from libpl.pdcode import PlanarDiagram
from ricci_flow import TriangleMesh, ThurstonCPMetric
import numpy as np

class Edge(object):
    def __init__(self, tail, head):
        self.tail = tail
        self.head = head
        self.adj_faces = []
    def __repr__(self):
        return "(%s,%s)"%(self.tail,self.head)
    def __str__(self):
        return repr(self)

class Face(object):
    def __init__(self, edges, signs):
        self.edges = edges
        self.signs = signs
        for e in edges:
            e.adj_faces.append(self)
    def __len__(self):
        return len(self.edges)
    def __repr__(self):
        return repr(self.edges)
    def __str__(self):
        return str(self.edges)

def split_edge(verts, edges, edge, n):
    i,j = edge.tail, edge.head
    N = len(verts)
    new_verts = range(N, N+n)
    verts.extend(new_verts)
    run = [i] + new_verts + [j]

    new_edges = []
    for tail,head in zip(run, run[1:]):
        new_edges.append(Edge(tail,head))
    edges.remove(edge)
    edges.extend(new_edges)

    for face in edge.adj_faces:
        i = face.edges.index(edge)
        if face.signs[i] == 1:
            face.edges = face.edges[:i] + new_edges + face.edges[i+1:]
            face.signs = face.signs[:i] + (1,)*len(new_edges) + face.signs[i+1:]
        else:
            face.edges = face.edges[:i] + list(reversed(new_edges)) + face.edges[i+1:]
            face.signs = face.signs[:i] + (0,)*len(new_edges) + face.signs[i+1:]

    if len(edge.adj_faces) == 1:
        return new_verts
    else:
        return []

def verts_on_face(face):
    verts = []
    for ei,ej in zip(face.edges, face.edges[1:]+face.edges[:1]):
        if ei.tail == ej.tail or ei.tail == ej.head:
            verts.append(ei.tail)
        elif ei.head == ej.tail or ei.head == ej.head:
            verts.append(ei.head)
    return verts

def triangulate(face):
    verts = verts_on_face(face)
    assert len(verts) == len(face)
    if len(face) == 3:
        yield verts
    else:
        e0 = verts[0]
        for e1, e2 in zip(verts[1:], verts[2:]):
            yield (e0, e1, e2)

def pdcode_to_mesh(L, boundary=None):
    N = L.ncross
    verts = range(N)
    edges = [Edge(e.tail, e.head) for e in L.edges]
    faces = [Face([edges[i] for i in face.edges], face.signs) for face in L.faces]

    max_fv = 0
    max_face = None
    for face in faces:
        if len(face) > max_fv:
            max_fv = len(face)
            max_face = face
    faces.remove(max_face)
    for edge in max_face.edges:
        edge.adj_faces.remove(max_face)
    b_verts = verts_on_face(max_face)

    splits = dict()
    for face in faces:
        if len(face) < 3:
            for edge in face.edges:
                splits[edge] = 2

    for edge, nsplit in splits.iteritems():
        b_verts.extend(split_edge(verts, edges, edge, nsplit))

    if boundary is not False:
        assert 1 == len(verts) + len(faces) - len(edges)
    else:
        assert 2 == len(verts) + len(faces) - len(edges)
    assert len(verts) == len(set(verts))
    assert len(edges) == len(set(edges))
    assert len(faces) == len(set(faces))

    m_verts = verts
    m_edges = set()
    m_faces = []

    for edge in edges:
        pass#m_edges.add(frozenset([edge.tail, edge.head]))
    for face in faces:
        print "triangulating face %s"%(face,)
        for tri in triangulate(face):
            print "  triangle: %s"%(tri,)
            i,j,k = tri
            if i != k and j != k and i != j:
                m_faces.append(frozenset(tri))
                #m_edges.add(frozenset([i,k]))
                #m_edges.add(frozenset([i,j]))
                #m_edges.add(frozenset([j,k]))
            else:
                print "    skipping face %s edge %s"%(tri,(i,k))

    print "~~~"
    print "triangulated faces"
    print m_faces

    for face in m_faces:
        i,j,k = face
        m_edges.add(frozenset([i,j]))
        m_edges.add(frozenset([k,j]))
        m_edges.add(frozenset([i,k]))

    m_edges = list(m_edges)
    return TriangleMesh(m_verts, b_verts, m_edges, m_faces), edges

def faces_around_vert(mesh, v_i, first_face, l_ij):
    assert v_i in first_face
    adj_faces = [f_ijk for f_ijk in mesh.faces if v_i in f_ijk]
    print adj_faces

    assert v_i in l_ij
    v_j, = l_ij - set([v_i])

    assert v_j in first_face
    adj_faces.remove(first_face)
    yield first_face
    l_ij = first_face - set([v_j])

    while adj_faces:
        next_face = [face for face in adj_faces if len(l_ij & face) == 2]
        if next_face:
            next_face = next_face[0]
        else:
            print "dead on del"
            break # hit a boundary face
        adj_faces.remove(next_face)
        yield next_face
        v_j, = next_face - l_ij
        l_ij = frozenset([v_i, v_j])

def boundary_edge_on(mesh, v_i):
    return (edge for edge in mesh.edges if
            (v_i in edge and
             iter((edge - set([v_i]))).next() in mesh.b_verts)).next()

def adjacent_faces(faces, face):
    for a_face in faces:
        if len(face & a_face) == 2:
            yield a_face
def adjacent_faces_and_edge(faces, face):
    for a_face in faces:
        e = face & a_face
        if len(e) == 2:
            yield a_face, e

def adj_or_face_and_edge(faces, face):
    i,j,k = face
    for a_face in faces:
        if i in a_face and j in a_face:
            yield a_face, (i,j)
        elif j in a_face and k in a_face:
            yield a_face, (j,k)
        elif k in a_face and i in a_face:
            yield a_face, (k,i)

def orient_faces(faces):
    edges = []
    oriented = []
    to_orient = list(faces)
    adj_queue = set()
    f_0 = to_orient.pop()
    i,j,k = f_0
    edges.extend([(i,j),(j,k),(k,i)])
    oriented.append((i,j,k))
    adj_queue.update(adjacent_faces_and_edge(to_orient, f_0))

    while to_orient:
        # Pop the next face to orient
        F,e = adj_queue.pop()
        if F not in to_orient:
            continue

        v_k, = F-e
        v_i, v_j = e
        if (v_i, v_j) in edges:
            i,j,k = v_j,v_i,v_k
        else:
            i,j,k = v_i,v_j,v_k

        edges.extend([(i,j),(j,k),(k,i)])
        oriented.append((i,j,k))

        to_orient.remove(F)
        adj_queue.update(adjacent_faces_and_edge(to_orient, F))

    return oriented

def u_theta(theta):
    return np.array([np.cos(theta), np.sin(theta)])
def embed_faces(cpm):
    mesh = cpm._mesh
    x = np.zeros((cpm._n, 2))
    phi = dict()
    pi = np.pi

    faces = orient_faces(mesh.faces)
    to_embed = list(faces)
    embed_queue = set()

    F_0 = to_embed.pop()
    i,j,k = F_0
    x[i] = (0,0)
    x[j] = (cpm.length((i,j)),0)
    phi_ik = cpm.angle(set(F_0), i) %(2*pi) # phi_ik

    #print cpm.length((i,k)), cpm.angle(set(F_0),i)
    #print phi_ik
    #print u_theta(phi_ik)
    #print x[i]


    x[k] = cpm.length((i,k))*u_theta(phi_ik)
    phi_jk = (pi - cpm.angle(set(F_0),j)) %(2*pi) # phi_jk

    phi[i,j] = 0
    phi[j,i] = pi
    phi[j,k] = phi_jk
    phi[k,j] = pi+phi_jk
    phi[i,k] = phi_ik
    phi[k,i] = pi+phi_ik
    #print F_0

    #plt.figure()
    #plt.gca().set_aspect('equal')
    #plt.triplot(x[:,0], x[:,1], faces, 'go-')
    #plt.show()

    embed_queue.update(adj_or_face_and_edge(to_embed, F_0))
    while to_embed:
        F,e = embed_queue.pop()

        if F not in to_embed:
            continue
        #print F, e

        i,j = e
        if F[(F.index(i)+1)%3] != j:
            k = F[(F.index(i)+1)%3]
            j,i = e
        else:
            k = F[(F.index(i)+2)%3]
        #print i,j,k

        # We already know x[i], x[j]. Only have to find x[K].
        phi_ik = (phi[i,j] + cpm.angle(set(F), i)) %(2*pi)
        phi_jk = (phi[j,i] - cpm.angle(set(F), j)) %(2*pi)
        if (j,k) not in phi and (i,k) not in phi:
            #print cpm.length((i,k)), cpm.angle(set(F),i)
            #print phi[i,j]/pi, phi_ik/pi
            #print u_theta(phi_ik)
            #print x[i]
            x[k] = x[i] + cpm.length((i,k))*u_theta(phi_ik)
        #print (x[k] - (x[j] + cpm.length((j,k))*u_theta(phi_jk)))

        if (j,k) not in phi:
            phi[j,k] = phi_jk
            phi[k,j] = pi+phi_jk
        else:
            pass
            #print phi[j,k] - phi_jk
            #assert np.isclose(phi[j,k], phi_jk)
        if (i,k) not in phi:
            phi[i,k] = phi_ik
            phi[k,i] = pi+phi_ik
        else:
            pass
            #print phi[i,k] - phi_ik
            #assert np.isclose(phi[i,k], phi_ik)


        to_embed.remove(F)
        embed_queue.update(adj_or_face_and_edge(to_embed, F))
        # plt.figure()
        # plt.gca().set_aspect('equal')
        # plt.triplot(x[:,0], x[:,1], faces, 'go-')
        # plt.show()

    #print x
    return x, faces

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import math

if __name__ == "__main__":
    test_faces = [
        frozenset([0,1,2]),
        frozenset([0,2,3]),
        frozenset([0,3,4]),
        frozenset([0,4,1]),
    ]
    test_verts = range(5)
    test_edges = [
        frozenset([0,1]),
        frozenset([0,2]),
        frozenset([0,3]),
        frozenset([0,4]),
        frozenset([1,2]),
        frozenset([2,3]),
        frozenset([3,4]),
        frozenset([4,1]),
    ]
    test_bverts = range(1,5)

    test_mesh = TriangleMesh(
        test_verts,
        test_bverts,
        test_edges,
        test_faces
    )

    K = np.zeros((5,))
    test_g = ThurstonCPMetric.from_triangle_mesh(test_mesh)
    s = len(test_bverts)
    for v in test_bverts:
        K[v] = (2*np.pi)/s
    # test_flat = test_g.newton(target_K=K, dt=0.05, thresh=1e-8)
    # xy, triangles = embed_faces(test_flat)

    # M = np.zeros((5,5))
    # for k,v in test_flat._l.iteritems():
    #     M[k] = v
    # print M

    # plt.figure()
    # plt.gca().set_aspect('equal')
    # plt.triplot(xy[:,0], xy[:,1], triangles, 'go-')
    # plt.show()

#    import sys
#    sys.exit(0)


    #tref = PlanarDiagram.torus_knot(2,3)
    #tref = PlanarDiagram.unknot(1)
    #with open("6a5and6n1.kt") as f:
    #    tref = PlanarDiagram.read_knot_theory(f)
    tref = PlanarDiagram.db_link(6,2,True)
    mesh, real_edges = pdcode_to_mesh(tref)
    #print mesh.verts
    #print mesh.edges
    #print mesh.faces
    print mesh.chi()

    cpm = ThurstonCPMetric.from_triangle_mesh(mesh)
    print cpm.gb_chi()
    #print cpm.newton().gb_chi()

    dummy_bv = list(mesh.b_verts)
    real_bv = []
    for crossing in tref.crossings:
        i = crossing.index
        if i in dummy_bv:
            dummy_bv.remove(i)
            real_bv.append(i)


    K = np.zeros((cpm._n))
    s = len(mesh.b_verts)
    if len(dummy_bv) >= 2*len(real_bv):
        s = len(dummy_bv)
        r = len(real_bv)
        t = .1#.7#-4/(r+s)
        for v in mesh.b_verts:#real_bv:
            K[v] = -t*np.pi/2
        for v in dummy_bv:
            K[v] = (2*np.pi + t*r*np.pi/2)/s
        # print K
        #sys.exit(0)
    else:
        for v in mesh.b_verts:
            K[v] = (2*np.pi)/s


    assert np.isclose(sum(K),2*np.pi)
    flatm = cpm.newton(target_K=K, dt=1, thresh=1e-2)
    xy, triangles = embed_faces(flatm)

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(xy[:,0], xy[:,1], triangles, 'co:')
    for edge in real_edges:
        plt.plot(xy[[edge.tail, edge.head],0],
                 xy[[edge.tail, edge.head],1], 'k-')
    plt.show()
