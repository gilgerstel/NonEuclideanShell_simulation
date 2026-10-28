import os
import numpy as np
from meshpy.triangle import MeshInfo, build


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import inspect, re
import sympy as sp
from sympy.printing.octave import octave_code


u, v = sp.symbols('u v')
_sympy_namespace = {name: obj for name, obj in sp.__dict__.items()}

def make_flag_func(umin, umax, vmin, vmax, tolu=0.05, tolv=0):
    def flag(u,v):
        if u < umin + tolu or u > umax - tolu or \
           v < vmin + tolv or v > vmax - tolv:
            return -1
        else:
            return -2
    return flag

import os
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from meshpy.triangle import MeshInfo, build

def create_open_annulus_input_files(
    sim_name,
    output_dir,
    r_in,
    r_out,
    vmin = -3*np.pi/4,
    vmax =  3*np.pi/4,
    max_area = 0.001,
    n_outer  = 160,   # boundary samples (outer arc)
    n_inner  =  80,   # boundary samples (inner arc)
    n_cut    =   2,   # samples along each radial cut (>=2; 2 = just endpoints)
    band_dist = 0.05  # Euclidean band thickness around each cut ray segment
):
    """
    Open annulus sector:
        r_in <= r <= r_out, vmin <= theta <= vmax

    Writes:
      {sim_name}_p-distmesh.dat   (u=r, v=theta)
      {sim_name}_t-distmesh.dat   (0-based triangles)
      {sim_name}_Vertices         4×N: [index; u; v; flag]
      {sim_name}_Faces            10×M (neighbor format)

    Flags:
      -1 : within Euclidean distance <= band_dist of the *ray segments*
           at v=vmin and v=vmax (projection t in [r_in, r_out])
      -2 : otherwise
    """
    os.makedirs(output_dir, exist_ok=True)

    # ---- 1) Boundary polygon in (x,y), CCW, no duplicate consecutive points ----
    th_outer = np.linspace(vmin, vmax, n_outer, endpoint=True)
    outer = np.c_[r_out*np.cos(th_outer), r_out*np.sin(th_outer)]

    cut_max_r = np.linspace(r_out, r_in, max(2, n_cut), endpoint=True)
    radial_vmax = np.c_[cut_max_r*np.cos(vmax), cut_max_r*np.sin(vmax)][1:]  # drop dup

    th_inner = np.linspace(vmax, vmin, n_inner, endpoint=True)
    inner = np.c_[r_in*np.cos(th_inner), r_in*np.sin(th_inner)][1:]          # drop dup

    cut_min_r = np.linspace(r_in, r_out, max(2, n_cut), endpoint=True)
    radial_vmin = np.c_[cut_min_r*np.cos(vmin), cut_min_r*np.sin(vmin)]
    radial_vmin = radial_vmin[1:-1]  # drop both ends (dup with inner end & outer start)

    pts = np.vstack([outer, radial_vmax, inner, radial_vmin])
    Np = len(pts)
    facets = [(i, (i+1) % Np) for i in range(Np)]

    # ---- 2) Triangulate (uniform-ish area via max_area) ----
    info = MeshInfo()
    info.set_points(pts.tolist())
    info.set_facets(facets)
    mesh = build(info, refinement_func=lambda verts, area: area > max_area)

    P_xy = np.array(mesh.points, dtype=float)   # (N,2)
    T_raw = np.array(mesh.elements, dtype=int)  # (M,3)

    # ---- 3) Polar coords limited to [vmin, vmax] ----
    x, y = P_xy[:,0], P_xy[:,1]
    u = np.hypot(x, y)
    v = np.arctan2(y, x)
    v[v < vmin] += 2*np.pi
    v[v > vmax] -= 2*np.pi
    v = np.clip(v, vmin, vmax)

    # ---- 4) Save raw point/triangle files ----
    with open(os.path.join(output_dir, f"{sim_name}_p-distmesh.dat"), 'w') as f:
        for (ru, th) in zip(u, v):
            f.write(f"{ru:.7f}\t{th:.7f}\n")

    with open(os.path.join(output_dir, f"{sim_name}_t-distmesh.dat"), 'w') as f:
        for (a, b, c) in T_raw:
            f.write(f"{a+1}\t{b+1}\t{c+1}\n") 

    # ---- 5) Flags: band around *ray segments* (not infinite lines) ----
    # Ray direction & normal at angle a:
    def ray_dir(a): return np.array([np.cos(a), np.sin(a)])
    def ray_nrm(a): return np.array([-np.sin(a), np.cos(a)])  # perpendicular

    dirs = [ray_dir(vmin), ray_dir(vmax)]
    nrms = [ray_nrm(vmin), ray_nrm(vmax)]

    flags = np.full(len(P_xy), -2, dtype=int)
    for d, n in zip(dirs, nrms):
        # perpendicular distance to the line through origin at angle a
        d_perp = np.abs(n[0]*x + n[1]*y)
        # projection coordinate along the ray direction
        t = d[0]*x + d[1]*y
        # keep only points whose closest point on the line lies on the *segment* [r_in, r_out]
        near_segment = (d_perp <= band_dist) & (t >= r_in) & (t <= r_out)
        flags[near_segment] = -1

    # ---- 6) Vertices file ----
    N = len(u)
    V = np.zeros((4, N), dtype=float)
    V[0, :] = np.arange(N)
    V[1, :] = u
    V[2, :] = v
    V[3, :] = flags

    with open(os.path.join(output_dir, f"{sim_name}_Vertices"), 'w') as f:
        f.write(f"{N}\n")
        for j in range(N):
            f.write(f"{int(V[0,j])}\t{V[1,j]:.7f}\t{V[2,j]:.7f}\t{int(V[3,j])}\n")

    # ---- 7) Faces (10×M) with neighbors ----
    M = len(T_raw)
    edge_map = {}
    for tri_idx, tri in enumerate(T_raw):
        for k in range(3):
            a, b = tri[(k+1)%3], tri[(k+2)%3]
            edge_map.setdefault(tuple(sorted((a,b))), []).append(tri_idx)

    neigh_tri  = -np.ones((M,3), dtype=int)
    neigh_vert = -np.ones((M,3), dtype=int)
    for i, tri in enumerate(T_raw):
        for k in range(3):
            a, b = tri[(k+1)%3], tri[(k+2)%3]
            adj = edge_map.get(tuple(sorted((a,b))), [])
            other = next((j for j in adj if j != i), None)
            if other is not None:
                neigh_tri[i, k] = other
                opp = [vv for vv in T_raw[other] if vv not in (a, b)][0]
                neigh_vert[i, k] = opp

    T = np.vstack([np.arange(M), T_raw.T, neigh_vert.T, neigh_tri.T]).astype(int)
    with open(os.path.join(output_dir, f"{sim_name}_Faces"), 'w') as f:
        f.write(f"{M}\n")
        for j in range(M):
            f.write("\t".join(str(x) for x in T[:, j]) + "\n")

        # ---- 8) Plot triangle mesh with visible triangles + fixed in red ----
    x = P_xy[:, 0]
    y = P_xy[:, 1]

    fig, ax = plt.subplots(figsize=(6, 6))
    # light gray mesh lines so triangles are visible
    ax.triplot(x, y, T_raw, color='0.6', lw=0.4)

    # red dots for fixed vertices
    fixed_idx = np.where(flags == -1)[0]
    if fixed_idx.size:
        ax.scatter(x[fixed_idx], y[fixed_idx], s=20, c='red', marker='o',
                   zorder=3, label='fixed')

    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{sim_name} open annulus sector (triangles; fixed in red)')
    if fixed_idx.size:
        ax.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


def create_input_files(sim_name,
                       output_dir,
                       domain_params,
                       max_area=0.001,
                       tols=(0,0)):
    """
    Creates:
      {sim_name}_p-distmesh.dat   raw distmesh points
      {sim_name}_t-distmesh.dat   raw distmesh triangles
      {sim_name}_Vertices         MATLAB-style 4×N array
      {sim_name}_Faces            MATLAB-style 10×M array
    """
    os.makedirs(output_dir, exist_ok=True)
    # 1) build boundary
    
    umin,umax,vmin,vmax = domain_params
    pts    = [(umin,vmin),(umax,vmin),(umax,vmax),(umin,vmax)]
    facets = [(i,(i+1)%4) for i in range(4)]

    # 2) refine
    info = MeshInfo()
    info.set_points(pts)
    info.set_facets(facets)
    mesh = build(info, refinement_func=lambda verts, area: area>max_area)
    P_raw = np.array(mesh.points)      # (N,2)
    T_raw = np.array(mesh.elements)    # (M,3)

    # 3) save distmesh files
    pfile = os.path.join(output_dir, f"{sim_name}_p-distmesh.dat")
    tfile = os.path.join(output_dir, f"{sim_name}_t-distmesh.dat")
    with open(pfile,'w') as f:
        # f.write(f"{len(P_raw)}\n")
        for (u,v) in P_raw:
            f.write(f"{u:.7f}\t{v:.7f}\n")
    with open(tfile,'w') as f:
        for (v0,v1,v2) in T_raw:
            f.write(f"{v0+1}\t{v1+1}\t{v2+1}\n")

    N = len(P_raw)
    V = np.zeros((4, N), dtype=float)
    V[0, :] = np.arange(N)
    V[1:3, :] = P_raw.T
    # apply your flag_func to each (u,v)
    tolu, tolv = tols
    flag_func = make_flag_func(umin, umax, vmin, vmax, tolu=tolu, tolv=tolv)

    V[3, :] = [flag_func(u, v) for u, v in P_raw]

    # write out the Vertices file
    vfile = os.path.join(output_dir, f"{sim_name}_Vertices")
    with open(vfile, 'w') as f:
        f.write(f"{N}\n")
        for j in range(N):
            # MATLAB fprintf('%5.0f\t%10.7f\t%10.7f\t%5.0f\n', P(:,j))
            f.write(f"{int(V[0,j])}\t{V[1,j]:.7f}\t{V[2,j]:.7f}\t{int(V[3,j])}\n")
    # 4) build MATLAB-style P,T
    M = len(T_raw)

    # 4a) build an edge→triangles map
    edge_map = {}
    for tri_idx, tri in enumerate(T_raw):
        for k in range(3):
            v1 = tri[(k+1)%3]
            v2 = tri[(k+2)%3]
            key = tuple(sorted((v1,v2)))
            edge_map.setdefault(key, []).append(tri_idx)

    # 4b) initialize neighbor arrays
    neigh_tri  = -np.ones((M,3), dtype=int)
    neigh_vert = -np.ones((M,3), dtype=int)

    # 4c) fill them
    for i, tri in enumerate(T_raw):
        for k in range(3):
            # edge opposite vertex k in triangle i
            a, b = tri[(k+1)%3], tri[(k+2)%3]
            adj = edge_map.get(tuple(sorted((a,b))), [])
            # find the “other” triangle in adj
            other = next((j for j in adj if j != i), None)
            if other is not None:
                neigh_tri[i, k] = other
                # the vertex in 'other' that is *not* a or b
                opp = [v for v in T_raw[other] if v not in (a,b)][0]
                neigh_vert[i, k] = opp

    # 4d) stack into 10×M array
    #   row 0: triangle index
    #   rows 1–3: vertices
    #   rows 4–6: opposite-vertex indices from each neighbor
    #   rows 7–9: neighbor-triangle indices
    T = np.vstack([
        np.arange(M),           # row 0
        T_raw.T,                # rows 1–3
        neigh_vert.T,           # rows 4–6
        neigh_tri.T             # rows 7–9
    ]).astype(int)

    # 5) write out the Faces file
    ffile = os.path.join(output_dir, f"{sim_name}_Faces")
    with open(ffile, 'w') as f:
        f.write(f"{M}\n")
        for j in range(M):
            f.write("\t".join(str(x) for x in T[:,j]) + "\n")
    
      # 6) **Plot the triangle mesh** in the UV‐plane (inline, no saving)
    fig, ax = plt.subplots()
    ax.triplot(P_raw[:,0], P_raw[:,1], T_raw, color='black', lw=0.5)

    # highlight clamped points (flag = -1)
    clamped_idx = np.where(V[3,:] == -1)[0]
    ax.scatter(P_raw[clamped_idx,0], P_raw[clamped_idx,1],
               c='red', s=15, label='clamped')

    ax.set_xlabel('u')
    ax.set_ylabel('v')
    ax.set_title(f'{sim_name} UV-mesh')
    ax.legend(loc='upper right', markerscale=1.5)
    plt.tight_layout()
    plt.show()



def make_in_file(params, output_dir, sim_name):
    lines = []
    lines.append(params['vertices_file'])
    lines.append(params['faces_file'])
    lines.append(str(params['loops']))
    lines.append(str(params['save_every']))

    # flatten abar, bbar and pos0, casting to strings
   # Reference first, second, and third fundamental forms
    lines.extend(str(x) for x in params['abar'])
    lines.extend(str(x) for x in params['bbar'])
    lines.extend(str(x) for x in params['cbar'])

    lines.append(str(params['c_prefactor']))

    lines.extend([
        str(params['K1111']),
        str(params['K1112']),
        str(params['K1122']),
        str(params['K1212']),
        str(params['K1222']),
        str(params['K2222']),
    ])

    lines.extend([
        str(params['thickness']),
        str(params['E']),
        str(params['nu'])
    ])
    lines.extend(str(x) for x in params['pos0'])
    lines.extend([str(params['eta']), str(params['lambdaG']), str(params['muG']) ])

    # γ‐table
    for k in (0,1):
        for i in (0,1):
            for j in (0,1):
                lines.append(str(params['gamma'][k][i][j]))

    # final three
    lines.extend([
        str(params['thickness_adjust']),
        str(params['metric_adjust']),
        str(params['connection_adjust']),
        str(int(params['restart']))
    ])

    in_file = os.path.join(output_dir, f'run_input_{sim_name}.txt')
    with open(in_file,'w') as f:
        f.write("\n".join(lines))
    return in_file


def plot_uv_surface(x_func, y_func, z_func, umin, umax, vmin, vmax, nu=50, nv=50):
    u = np.linspace(umin, umax, nu)
    v = np.linspace(vmin, vmax, nv)
    U, V = np.meshgrid(u, v)
    X = x_func(U, V)
    Y = y_func(U, V)
    Z = z_func(U, V)

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    # compute the span of each axis
    # rx = X.max() - X.min()
    # ry = Y.max() - Y.min()
    # rz = Z.max() - Z.min()

    # set the box aspect to the actual data-ratios
    # ax.set_box_aspect((rx, ry, rz))
    # ax.set_box_aspect((1, 1, 1))

    
    surf = ax.plot_surface(X, Y, Z, 
                           rstride=1, cstride=1, edgecolor='none')
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    fig.colorbar(surf, ax=ax, shrink=0.5, label='Z value')
    plt.tight_layout()
    return ax

import numpy as np
import os
import matplotlib.pyplot as plt
from meshpy.triangle import MeshInfo, build



import numpy as np
import os
import matplotlib.pyplot as plt
from meshpy.triangle import MeshInfo, build

import numpy as np
import os
import matplotlib.pyplot as plt
from meshpy.triangle import MeshInfo, build

def create_annulus_input_files(
    sim_name,
    output_dir,
    r_in,
    r_out,
    max_area=0.001,
    n_outer=120,
    n_inner=60,
    tols=(0, 0)  # (in_flag, out_flag)
):
    """
    Writes annulus mesh files with columns: index, u=r, v=theta, flag,
    and Faces in the 10-column format of your example.
    """
    os.makedirs(output_dir, exist_ok=True)

    # -- Boundary points in (x, y) --
    theta_outer = np.linspace(0, 2*np.pi, n_outer, endpoint=False)
    theta_inner = np.linspace(0, 2*np.pi, n_inner, endpoint=False)
    outer = np.c_[r_out * np.cos(theta_outer), r_out * np.sin(theta_outer)]
    inner = np.c_[r_in * np.cos(theta_inner[::-1]), r_in * np.sin(theta_inner[::-1])]
    pts = np.vstack([outer, inner])

    facets_outer = [(i, (i+1)%n_outer) for i in range(n_outer)]
    facets_inner = [(n_outer + i, n_outer + (i+1)%n_inner) for i in range(n_inner)]
    facets = facets_outer + facets_inner
    hole_marker = [[0.0, 0.0]]

    info = MeshInfo()
    info.set_points(pts.tolist())
    info.set_facets(facets)
    info.set_holes(hole_marker)

    mesh = build(info, refinement_func=lambda verts, area: area > max_area)
    P_cart = np.array(mesh.points)
    T_raw = np.array(mesh.elements)

    # Convert all mesh points to (u, v) = (r, theta)
    u = np.sqrt(P_cart[:, 0]**2 + P_cart[:, 1]**2)
    v = np.arctan2(P_cart[:, 1], P_cart[:, 0])
    v = np.mod(v, 2*np.pi)  # wrap theta to [0, 2pi)

    # Save raw point and triangle files (still in r/theta)
    pfile = os.path.join(output_dir, f"{sim_name}_p-distmesh.dat")
    tfile = os.path.join(output_dir, f"{sim_name}_t-distmesh.dat")
    with open(pfile, 'w') as f:
        for (rval, tval) in zip(u, v):
            f.write(f"{rval:.7f}\t{tval:.7f}\n")
    with open(tfile, 'w') as f:
        for (v0, v1, v2) in T_raw:
            f.write(f"{v0}\t{v1}\t{v2}\n")

    # ---- Vertices file ----
    N = len(u)
    V = np.zeros((4, N), dtype=float)
    V[0, :] = np.arange(N)
    V[1, :] = u
    V[2, :] = v

    in_flag, out_flag = tols
    def assign_flag(r, theta):
        if in_flag == 0 and out_flag == 0:
            return -2
        if abs(r - r_in) < 1e-8 and in_flag != 0:
            return in_flag
        elif abs(r - r_out) < 1e-8 and out_flag != 0:
            return out_flag
        else:
            return -2

    V[3, :] = [assign_flag(rval, tval) for rval, tval in zip(u, v)]

    vfile = os.path.join(output_dir, f"{sim_name}_Vertices")
    with open(vfile, 'w') as f:
        f.write(f"{N}\n")
        for j in range(N):
            f.write(f"{int(V[0, j])}\t{V[1, j]:.7f}\t{V[2, j]:.7f}\t{int(V[3, j])}\n")

    # ---- Faces file (10 columns) ----
    M = len(T_raw)
    edge_map = {}
    for tri_idx, tri in enumerate(T_raw):
        for k in range(3):
            v1 = tri[(k+1)%3]
            v2 = tri[(k+2)%3]
            key = tuple(sorted((v1, v2)))
            edge_map.setdefault(key, []).append(tri_idx)

    neigh_tri = -np.ones((M, 3), dtype=int)
    neigh_vert = -np.ones((M, 3), dtype=int)
    for i, tri in enumerate(T_raw):
        for k in range(3):
            a, b = tri[(k+1)%3], tri[(k+2)%3]
            adj = edge_map.get(tuple(sorted((a, b))), [])
            other = next((j for j in adj if j != i), None)
            if other is not None:
                neigh_tri[i, k] = other
                opp = [v for v in T_raw[other] if v not in (a, b)][0]
                neigh_vert[i, k] = opp

    # 10-row face matrix (like your files)
    T = np.vstack([
        np.arange(M),        # 0: face index
        T_raw.T,             # 1–3: triangle vertices
        neigh_vert.T,        # 4–6: opposite-vertex indices
        neigh_tri.T          # 7–9: neighbor-triangle indices
    ]).astype(int)

    ffile = os.path.join(output_dir, f"{sim_name}_Faces")
    with open(ffile, 'w') as f:
        f.write(f"{M}\n")
        for j in range(M):
            f.write("\t".join(str(x) for x in T[:, j]) + "\n")

    # ---- Plot mesh in (r, theta) coordinates (optional) ----
    x = u * np.cos(v)
    y = u * np.sin(v)

    fig, ax = plt.subplots()
    ax.triplot(x, y, T_raw, color='black', lw=0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    ax.set_title(f'{sim_name} annulus mesh (x, y)')
    plt.tight_layout()
    plt.show()



# def func_to_matlab_str(func):
#     """
#     Turn a lambda like `lambda u,v: np.sin(u)*np.cos(v)` or `lambda u,v: u**2`
#     into a MATLAB/Octave‐style string: 'sin(u)*cos(v)' or 'u^2', etc.
#     """
#     # grab the RHS of the lambda
#     src = inspect.getsource(func).strip()
#     m = re.match(r'.*lambda\s*[\w, ]*:\s*(.*)', src)
#     if not m:
#         raise ValueError("Expected a one-liner lambda(u,v): …")
#     expr = m.group(1)

#     # 3) strip off the "np." (or "numpy.") prefix everywhere
#     expr = re.sub(r'\b(?:np|numpy)\.([A-Za-z_]\w*)', r'\1', expr)

#     # 4) sympify with *all* Sympy names available
#     sym = sp.sympify(expr, locals={**_sympy_namespace, 'u': u, 'v': v})

#     # 5) print out Octave/MATLAB code
#     return octave_code(sym)
def func_to_matlab_str(func):
    """
    Turn a lambda like `lambda u,v: np.sin(u)*np.cos(v)` or `lambda u,v: u**2`
    into a MATLAB‐style string: 'sin(u)*cos(v)' or 'u^2', etc.
    """
    # 1) grab the RHS of the lambda
    src = inspect.getsource(func).strip()
    m   = re.match(r'.*lambda\s*[\w, ]*:\s*(.*)', src)
    if not m:
        raise ValueError("Expected a one-liner lambda(u,v): …")
    expr = m.group(1)

    # 2) strip numpy prefixes
    expr = re.sub(r'\b(?:np|numpy)\.([A-Za-z_]\w*)', r'\1', expr)

    # 3) sympify with u,v in scope
    sym = sp.sympify(expr, locals={**_sympy_namespace, 'u': u, 'v': v})

    # 4) use Python‐printer, then tweak **→^ for MATLAB
    code = str(sym)
    code = code.replace('**', '^')
    code = code.replace(" ", "")
        # 5) fix scientific notation: e.g. 1.0e-5 → 1.0*10^(-5)
    code = re.sub(r'(\d+(\.\d+)?)[eE]([-+]?\d+)', r'\1*10^(\3)', code)
    return code


def calc_Gamma(a, u, v):
    """
    give a 2x2 matrix a with sympy functions
    returns the Levi Chivita connectoin
    """
    inva  = a.inv()
    da_du = a.diff(u)
    da_dv = a.diff(v)

    Gamma = [[[sp.simplify(sp.Rational(1, 2) * sum(
                    inva[k, l] * (
                    (da_du[l, j] if i == 0 else da_dv[l, j]) +
                    (da_du[l, i] if j == 0 else da_dv[l, i]) -
                    (da_du[i, j] if l == 0 else da_dv[i, j])
                    )
                for l in range(2)))
            for j in range(2)]
            for i in range(2)]
            for k in range(2)]

    Gamma_str = [[[octave_code(Gamma[k][i][j])
                for j in range(2)]
                for i in range(2)]
                for k in range(2)]
    return Gamma_str