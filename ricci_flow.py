import numpy
from scipy import sparse
from numpy.linalg import norm, inv, lstsq, solve
from itertools import combinations, izip, product
import functools
from numpy import dot
from scipy.sparse import linalg

sqrt = numpy.sqrt
cos = numpy.cos
pi = numpy.pi

numpy.seterr(all='raise')

class IdentityDictMap(object):
    def __init__(self, output=1, domain=None):
        self._o = output
        self._domain = domain
    def __getitem__(self, i):
        if self._domain and i not in self._domain:
            raise KeyError()
        return self._o

def partition_face(face):
    i,j,k = face
    yield i, (j,k)
    yield j, (k,i)
    yield k, (i,j)

def triangulate(poly):
    if len(poly) == 3:
        yield poly
    else:
        e0 = poly[0]
        for e1, e2 in izip(poly[1:], poly[2:]):
            yield (e0, e1, e2)

def read_obj_file(f):
    verts = []
    normals = []
    edges = []
    lengths = dict()
    faces = []
    boundaries = []
    for line in f:
        tok = line.split()
        if not tok:
            continue
        elif tok[0] == "#":
            continue
        elif tok[0] == "o":
            print "Reading %s"%tok[1]
        elif tok[0] == "v":
            verts.append(numpy.array([float(x) for x in tok[1:]]))
        elif tok[0] == "vn":
            normals.append(numpy.array([float(x) for x in tok[1:]]))
        elif tok[0] == "f":
            # Throw out normal info for now...
            poly = [int(vstr.split("/")[0])-1 for vstr in tok[1:]]
            for face in triangulate(poly):
                faces.append(frozenset(face))

    maxnorm = max(numpy.linalg.norm(v) for v in verts)
    verts = numpy.array(verts)
    verts = verts/maxnorm/2

    for face in faces:
        for edge in combinations(face, 2):
            elen = norm(verts[edge[0]] - verts[edge[1]])
            edge = frozenset(edge)
            if edge not in edges:
                edges.append(edge)
                lengths[edge] = elen
                boundaries.append(edge)
            else:
                boundaries.remove(edge)

    print "Done... %s boundary edges"%len(boundaries)
    print "%s vertices, %s edges, %s faces"%(len(verts), len(edges), len(faces))
    if boundaries:
        b_verts = sorted(list(reduce(lambda A,B: A.union(B), boundaries)))
    else:
        b_verts = []
    return verts, b_verts, edges, faces, lengths

class TriangleMesh(object):
    def __init__(self, verts, b_verts, edges, faces, bg_geom="euclidean"):
        self.verts = verts
        self.b_verts = b_verts
        self.edges = edges
        self.faces = faces
        self.bg_geom = bg_geom

    def adjacent_edges(self, vert):
        return [edge for edge in self.edges if vert in edge]

    def min_valence(self):
        return min([len([face for face in self.faces if vert in face]) for vert in self.verts])

    def chi(self):
        return len(self.verts) - len(self.edges) + len(self.faces)

class DiscreteRiemannianMetric(object):
    def __init__(self, mesh, length_map):
        self._n = len(mesh.verts)
        self._mesh = mesh

        # Lengths l_{ij}
        self._lmap = sparse.lil_matrix((self._n, self._n))

        # Angles \theta_i^{jk}
        self._theta = dict()

        # Curvatures K_i
        self._K = None

        for (i,j), l in length_map.iteritems():
            self._lmap[i,j] = l
            self._lmap[j,i] = l

        self.update()

    def update(self):
        # Set angles using law of cosines
        for face in self._mesh.faces:
            for i,(j,k) in partition_face(face):
                theta = self.compute_angle(face,i)
                self._theta[(i,j,k)] = theta
                self._theta[(i,k,j)] = theta

        # Set curvatures
        self._K = numpy.zeros((self._n,))
        for vert in self._mesh.verts:
            self._K[vert] = self.curvature(vert)

    def is_ok(self):
        for face in self._mesh.faces:
            edges = [frozenset(e) for e in combinations(face,2)]
            k = len(edges)
            for i in range(k):
                if self.length(edges[(i+2)%k]) + self.length(edges[(i+1)%k]) < self.length(edges[i]):
                    return False
        return True

    def length(self, edge):
        i,j = edge
        return self._lmap[i,j]

    def abc_for_vert(self, face, vert):
        assert(vert in face)
        other_v = list(face - set([vert]))
        edge_a = [vert, other_v[0]]
        edge_b = [vert, other_v[1]]
        edge_c = other_v

        return [self.length(e_i) for e_i in (edge_a,edge_b,edge_c)]

    def angle(self, face, vert):
        j,k = face - set([vert])
        return self._theta[(vert,j,k)]

    def law_of_cosines(self, a,b,c):
        if self._mesh.bg_geom == "euclidean":
            ratio = (a**2 + b**2 - c**2)/(2.0*a*b)

            return numpy.arccos(ratio)

    def compute_angle(self, face, vert):
        a,b,c = self.abc_for_vert(face,vert)
        #print a,b,c
        assert(a != 0 and b != 0 and c != 0)

        try:
            return self.law_of_cosines(a,b,c)
        except FloatingPointError:
            print self.is_ok()
            print "~~~~~~~~~~~~ error"
            print face, vert
            print a,b,c
            print (a**2 + b**2 - c**2)/(2.0*a*b)
            return 0.0001
            raise

    def curvature(self, vert):
        faces = [face for face in self._mesh.faces if vert in face]
        if vert in self._mesh.b_verts:
            return numpy.pi - sum([self.angle(face, vert) for face in faces])
        else:
            return 2*numpy.pi - sum([self.angle(face, vert) for face in faces])

    def total_curvature(self):
        return sum([self.curvature(vert) for vert in self._mesh.verts])

    def face_area(self, face):
        i,j,k = face
        gamma = self._theta[(i,j,k)]
        a,b = self.length((i,j)), self.length((i,k))

        if self._mesh.bg_geom == "euclidean":
            return .5*a*b*numpy.sin(gamma)

    def area(self):
        return sum([self.face_area(face) for face in self._mesh.faces])

    EPSILON_DICT = {
        'euclidean':   0,
        'spherical':   1,
        'hyperbolic': -1,
    }

    def gb_chi(self):
        epsilon = self.EPSILON_DICT[self._mesh.bg_geom]
        return (self.total_curvature() + epsilon*self.area())/(numpy.pi*2)

    def as_cp_metric(self, scheme="inversive"):
        gamma = numpy.zeros((self._n,))
        eta = sparse.lil_matrix((self._n, self._n))

        if scheme == "inversive" and self._mesh.bg_geom == "euclidean":
            for vert in self._mesh.verts:
                gamma[vert] = (1.0/3)*min(self.length(edge) for
                                          edge in self._mesh.adjacent_edges(vert))

            for edge in self._mesh.edges:
                i,j = edge
                struct_c = ((self.length(edge)**2 - gamma[i]**2 - gamma[j]**2)/
                            (2*gamma[i]*gamma[j]))
                eta[i,j] = struct_c
                eta[j,i] = struct_c

            ret = CirclePackingMetric(self._mesh, gamma, eta, IdentityDictMap())

        assert numpy.allclose(self._lmap.todense(), ret._l.todense())
        return ret

class ThurstonCPMetric(DiscreteRiemannianMetric):
    def __init__(self, mesh, radius_map, edge_weights):
        self._n = len(mesh.verts)
        self._mesh = mesh

        self._gamma = radius_map
        self.u = self.conf_factor(radius_map)
        self._l = dict()

        self._theta = dict()
        self._phi = edge_weights
        self.update()

    @classmethod
    def from_triangle_mesh(cls, mesh):
        """Create a new Thurston's CP Metric using a triangle mesh
        without metric (i.e. length) data"""
        n = len(mesh.verts)
        gamma = numpy.array([1 for v in mesh.verts])
        phi = sparse.dok_matrix((n,n))
        for edge in mesh.edges:
            i,j = edge
            phi[i,j] = 0
            phi[j,i] = 0

        return cls(mesh, gamma, phi)

    @classmethod
    def from_riemannian_metric(cls, g):
        pre_gamma = [[] for _ in g._mesh.verts]
        mesh = g._mesh
        n = len(mesh.verts)

        CP_RAD_SCHEME = 0

        if CP_RAD_SCHEME == 0:
            for face in mesh.faces:
                for i, opp_edge in partition_face(face):
                    j,k = opp_edge
                    pre_gamma[i].append(
                        .5*(g.length((k,i)) + g.length((i,j)) - g.length((j,k))))
            gamma = numpy.array([(1.0/len(g_ijk))*sum(g_ijk) for g_ijk in pre_gamma])

        elif CP_RAD_SCHEME == 1:
            gamma = numpy.array(
                [(2.5/3.0)*min(g.length(edge) for edge in mesh.adjacent_edges(vert)) for
                 vert in mesh.verts])

        print gamma

        # For Thurston's CP scheme, all adjacent circles
        # should be at least tangent (otherwise we're hitting)
        # an impossible triangle with the law of cosines
        for i,j in mesh.edges:
            assert gamma[i]+gamma[j] >= g.length((i,j))

        phi = sparse.dok_matrix((n, n))
        for i,j in mesh.edges:
            if mesh.bg_geom == "euclidean":
                g_i,g_j,l_ij = gamma[i], gamma[j], g.length((i,j))
                eta = .5*(l_ij**2 - g_i**2 - g_j**2)/(g_i*g_j)
                #eta = min(1, eta)
                phi_ij = numpy.arccos(eta)
                #phi_ij = min(pi/2, phi_ij)
                phi[i,j] = phi_ij
                phi[j,i] = phi_ij

        print phi.todense()
        ret = ThurstonCPMetric(mesh, gamma, phi)

        # This new metric should approximate the old
        #assert numpy.allclose(g._lmap.todense(), ret._l.todense())

        return ret

    def _constant_curvature_goal(self):
        return (2*self._mesh.chi()*numpy.pi)/self._n

    def gradient_descent(self, target_K=None, epsilon=0.05, thresh=0.01):
        if target_K is None:
            # Constant curvature
            target_K = self._constant_curvature_goal()

        g = ThurstonCPMetric(self._mesh, self._gamma, self._phi)

        K = g._K
        deltaK = target_K - K
        while numpy.max(deltaK) > thresh:
            g.u = g.u + epsilon*(deltaK)
            g.u = g.u - sum(g.u)/g._n
            g.update()
            K = g._K
            deltaK = target_K - K
            print numpy.max(deltaK)

        return g

    def newton(self, target_K=None, dt=0.05, thresh=1e-4):
        if target_K is None:
            # Constant curvature
            target_K = self._constant_curvature_goal()

        g = ThurstonCPMetric(self._mesh, self._gamma, self._phi)

        K = g._K
        DeltaK = target_K - K
        while numpy.max(numpy.abs(DeltaK)) > thresh:
            H = self.hessian()
            deltau = sparse.linalg.lsqr(H, DeltaK)[0]
            g.u -= dt*deltau
            g.update()
            K = g._K
            DeltaK = target_K - K
            print numpy.max(numpy.abs(DeltaK))

        return g

    def conf_factor(self, gamma):
        return numpy.log(gamma)

    def update(self):
        self._gamma = numpy.exp(self.u)

        for edge in self._mesh.edges:
            i,j = edge
            l = self.compute_length(edge)
            self._l[i,j] = l
            self._l[j,i] = l

        super(ThurstonCPMetric, self).update()

    def length(self, edge):
        i,j = edge
        return self._l[i,j]

    def compute_length(self, edge):
        i,j = list(edge)
        g_i, g_j = self._gamma[[i,j]]
        if self._mesh.bg_geom == "euclidean":
            return sqrt(2*g_i*g_j*cos(self._phi[i,j]) + g_i**2 + g_j**2)

    def _s(self, x):
        if self._mesh.bg_geom == "euclidean":
            return x

    def _tau2(self, l_jk, g_j, g_k):
        return .5*(l_jk**2 + g_j**2 - g_k**2)

    def _Theta(self, face):
        i,j,k = face
        theta = functools.partial(self.angle, face)
        cos = numpy.cos
        return numpy.array(
            ((-1,            cos(theta(k)), cos(theta(j))),
             (cos(theta(k)), -1,            cos(theta(i))),
             (cos(theta(j)), cos(theta(i)), -1           )
         ))

    def hessian(self):
        n = len(self._mesh.verts)
        H = dict()#sparse.dok_matrix((n,n))
        t = self._tau2
        for face in self._mesh.faces:
            i,j,k = face
            l_k, l_i, l_j = self._l[i,j], self._l[j,k], self._l[k,i]
            g_i, g_j, g_k = self._gamma[[i,j,k]]
            th_i, th_j, th_k = (
                self.angle(face, i),
                self.angle(face, j),
                self.angle(face, k))

            A = self.face_area(face)
            L = numpy.diag((l_i, l_j, l_k))
            D = numpy.array(
                ((0,              t(l_i,g_j,g_k), t(l_i,g_k,g_j)),
                 (t(l_j,g_i,g_k), 0,              t(l_j,g_k,g_i)),
                 (t(l_k,g_i,g_j), t(l_k,g_j,g_i), 0            )))
            Theta = numpy.cos(numpy.array(
                ((pi,   th_k, th_j),
                 (th_k, pi,   th_i),
                 (th_j, th_i, pi))))

            Tijk = -.5/A * (L.dot(Theta).dot(inv(L)).dot(D))
            for a,row in izip((i,j,k), Tijk):
                for b,dtheta in izip((i,j,k), row):
                    if (a,b) in H:
                        H[a,b] += dtheta
                    else:
                        H[a,b] = dtheta
        Hm = sparse.dok_matrix((n,n))
        for (du_i,dtheta_j), val in H.iteritems():
            Hm[du_i, dtheta_j] = val
        return Hm.tocsr()

class CirclePackingMetric(DiscreteRiemannianMetric):
    def __init__(self, mesh, radius_map, struct_coeff, scheme_coeff):
        self._n = len(mesh.verts)
        self._mesh = mesh
        self._eta = struct_coeff
        self._eps = scheme_coeff
        self.u = self.conf_factor(radius_map)
        self._gamma = radius_map
        self._theta = dict()
        self._l = sparse.dok_matrix((self._n, self._n))
        self.update()

    def _s(self, x):
        if self._mesh.bg_geom == "euclidean":
            return x

    def _tau(self, i,j,k):
        if self._mesh.bg_geom == "euclidean":
            return .5*(self._l[j,k]**2 +
                       self._eps[j] * (self._gamma[j]**2) +
                       self._eps[k] * (self._gamma[k]**2))

    def _L(self, i,j,k):
        return numpy.diag((self._l[j,k],
                           self._l[i,k],
                           self._l[i,j]))

    def _Theta(self, face):
        i,j,k = face
        theta = functools.partial(self.angle, face)
        cos = numpy.cos
        return numpy.array(
            ((-1,            cos(theta(k)), cos(theta(j))),
             (cos(theta(k)), -1,            cos(theta(i))),
             (cos(theta(j)), cos(theta(i)), -1           )
         ))

    def _D(self, i,j,k):
        return numpy.array(
            ((0,                self._tau(i,j,k), self._tau(i,k,j)),
             (self._tau(j,i,k), 0,                self._tau(j,k,i)),
             (self._tau(k,i,j), self._tau(k,j,i), 0               )
         ))

    def conf_factor(self, gamma):
        return numpy.log(gamma)

    def update_with_u(self, u):
        self.u = u
        self._gamma = numpy.exp(u)
        self.update()

    def update(self):
        for edge in self._mesh.edges:
            i,j = edge
            l = self.compute_length(edge)
            self._l[i,j] = l
            self._l[j,i] = l

        super(CirclePackingMetric, self).update()

    def curvature_array(self):
        K = numpy.zeros(len(self._mesh.verts))
        for v_i in self._mesh.verts:
            K[v_i] = self.curvature(v_i)
        return K

    def hessian(self):
        n = len(self._mesh.verts)
        H = sparse.dok_matrix((n,n))
        for face in self._mesh.faces:
            i,j,k = face
            A = self.face_area(face)
            L = self._L(i,j,k)
            D = self._D(i,j,k)
            Theta = self._Theta(face)

            Tijk = -.5/A * (L.dot(Theta).dot(inv(L)).dot(D))
            for a,row in izip((i,j,k), Tijk):
                for b,dtheta in izip((i,j,k), row):
                    H[a,b] += dtheta
        return H

    def length(self, edge):
        i,j = edge
        return self._l[i,j]

    def compute_length(self, edge):
        i,j = list(edge)
        if self._mesh.bg_geom == "euclidean":
            return numpy.sqrt(2 * self._eta[i,j] * numpy.exp(self.u[i] + self.u[j]) +
                              self._eps[i] * numpy.exp(2 * self.u[i]) +
                              self._eps[j] * numpy.exp(2 * self.u[j]))

def unified_ricci_flow(mesh, g, Kbar, thresh=.001, dt=.05):
    print g.is_ok()
    cpm = g.as_cp_metric()
    print cpm.is_ok()
    print cpm.gb_chi()
    Kbar = numpy.array(Kbar)
    DeltaK = cpm.curvature_array() - Kbar

    newton = True

    while DeltaK.max() > thresh:
        #print cpm.angle(frozenset([27,28,22]), 27)
        #print cpm.length([27,28]), cpm.length([22,28]), cpm.length([27,22])
        #print cpm.face_area(frozenset([27,28,22]))
        #print cpm.is_ok()

        if newton:
            H = cpm.hessian()
            #H = H.todense()

            # Newton's method
            deltau = sparse.linalg.lsqr(H, DeltaK)[0]
            cpm.update_with_u( cpm.u + dt*deltau )
            #print [cpm.curvature(i) for i in mesh.verts],

        else:
            # Gradient Descent
            cpm.update_with_u( cpm.u + dt*DeltaK )

        DeltaK = cpm.curvature_array() - Kbar
        print "Max in \DeltaK: %s"%DeltaK.max()
        #raw_input("continue... ")
    return cpm


# verts = numpy.array([0,1,2,3])
# edges = numpy.array([frozenset([0,1]),frozenset([0,2]),frozenset([0,3]),
#          frozenset([1,2]),frozenset([1,3]),frozenset([2,3])])
# faces = numpy.array([frozenset(verts[[0,2,3]]),
#                      frozenset(verts[[0,1,3]]),
#                      frozenset(verts[[0,2,1]]),
#                      frozenset(verts[[1,2,3]]),])

# mesh = TriangleMesh(verts, edges, faces)
# g = DiscreteRiemannianMetric(
#     mesh,
#     {e:1 for e in edges}
# )

# print g.length(edges[2])
# print g.angle(faces[0], 0)*180/numpy.pi
# print g.curvature(0)
# print g.face_area(faces[0])
# print g.area()
# print g.total_curvature()
# print g.gb_chi()

def run_thcpm(objname="chair1"):
    with open("%s.obj"%objname) as fp:
        v,bv,e,f,l = read_obj_file(fp)

    mesh = TriangleMesh(range(len(v)),bv,e,f)
    g = DiscreteRiemannianMetric(mesh, l)

    cp = ThurstonCPMetric.from_triangle_mesh(mesh)
    cp.update()

    unif_cp = cp.newton()

def run(objname="bunny_1k"):
    fnames = [
        "tetra.obj",
        "cube.obj",
        "chair1.obj",
        "torus.obj",
        "teapot-low.obj",
        "wt_teapot.obj"]
    fnames = ["%s.obj"%objname]

    for fname in fnames:
        with open(fname) as f:
            verts, b_verts, edges, faces, lengths = read_obj_file(f)

        mesh = TriangleMesh(range(len(verts)),b_verts,edges,faces)
        g = DiscreteRiemannianMetric(mesh, lengths)

        print "Min vertex valence: %s" % mesh.min_valence()
        print "Mesh chi: %s, G.B. chi: %s" % (mesh.chi(),g.gb_chi())

    uniform_K = (2*mesh.chi()*numpy.pi)/len(verts)
    K = [uniform_K]*len(verts)

    res = unified_ricci_flow(mesh,g,K)


if __name__ == "__main__":
    run_thcpm()
