#%%

import numpy as np
import matplotlib.pyplot as plt

#%%

def vertices(xs, ys, xc, yc, z_in, z_out):

    point_labels = ['s0b', 's1b', 's2b', 's3b',
                    'r0b', 'r1b', 'r2b', 'r3b',
                    's0t', 's1t', 's2t', 's3t',
                    'r0t', 'r1t', 'r2t', 'r3t']
    
    cont0 = len(xs)
    lines_1 = [f'({xs[i]} {ys[i]} {z_in}) //point {i}\n' for i in range(cont0)]

    cont1 = len(xc)

    lines_2 = [f'({xc[i]} {yc[i]} {z_in}) //point {i + cont0}\n' for i in range(cont1)]

    lines_3 = [f'({xs[i]} {ys[i]} {z_out}) //point {i + cont0 + cont1}\n' for i in range(cont0)]

    lines_4 = [f'({xc[i]} {yc[i]} {z_out}) //point {i + cont0 + cont1 + cont0}\n' for i in range(cont1)]

    lines = ['vertices\n','(\n\n', '//square(z_in)\n'] + lines_1 + ['\n', '//circle(z_in)\n'] + lines_2 + ['\n', '//square(z_out)\n'] + lines_3 + ['\n', '//circle(z_out)\n'] + lines_4 + ['\n);\n\n']
    
    labeled_points = {point_labels[i]:i for i in range(2*(cont0 + cont1))}

    return lines, labeled_points

def evert(cxs, cys, cxc, cyc, z_in, z_out):

    #inner circle curvature
    cont0 = len(cxc)
    
    lines_1 = []
    lines_2 = []
    for i in range(cxc.shape[0]):
        line1 = ''
        line2 = ''
        for j in range(cxc.shape[1]):
            line1 += f"({cxc[i,j]} {cyc[i,j]} {z_in}) "
            line2 += f"({cxc[i,j]} {cyc[i,j]} {z_out}) "
        if cxc.shape[1] > 1:
            line1 = f"({line1})"
            line2 = f"({line2})"
        lines_1.append(line1)
        lines_2.append(line2)

    #inner square curvature
    cont1 = len(cxs)
    lines_3 = [[f'({cxs[i]} {cys[i]} {z_in})' for i in range(cont1)]]
    lines_4 = [[f'({cxs[i]} {cys[i]} {z_out})' for i in range(cont1)]]

    lines = [lines_1] + [lines_2] + lines_3 + lines_4

    return lines

def hex2D(p1, p2, p3, p4, labeled_points):
    line = f"hex ({labeled_points[p1+'b']} {labeled_points[p2+'b']} {labeled_points[p3+'b']} {labeled_points[p4+'b']} {labeled_points[p1+'t']} {labeled_points[p2+'t']} {labeled_points[p3+'t']} {labeled_points[p4+'t']})"
    return line

def blocks(labeled_points, ni, ns, nz, innerCircle_grading='(1 1 1)'):

    bk_lines = ['blocks\n', '(\n\n', '//block 0\n']
    bk_lines.append(hex2D('s1', 's0', 's3', 's2', labeled_points))
    bk_lines.append(f" square ({ns} {ns} {nz})")
    bk_lines.append(f" simpleGrading (1 1 1)\n\n") 
    bk_lines.append('//block1\n')

    bk_lines.append(hex2D('s0', 'r0', 'r3', 's3', labeled_points))
    bk_lines.append(f" innerCircle ({ni} {ns} {nz})")
    bk_lines.append(f" simpleGrading {innerCircle_grading}\n\n") 
    bk_lines.append('//block2\n')

    bk_lines.append(hex2D('s3', 'r3', 'r2', 's2', labeled_points))
    bk_lines.append(f" innerCircle ({ni} {ns} {nz})")
    bk_lines.append(f" simpleGrading {innerCircle_grading}\n\n") 
    bk_lines.append('//block3\n')

    bk_lines.append(hex2D('s2', 'r2', 'r1', 's1', labeled_points))
    bk_lines.append(f" innerCircle ({ni} {ns} {nz})")
    bk_lines.append(f" simpleGrading {innerCircle_grading}\n\n") 
    bk_lines.append('//block4\n')

    bk_lines.append(hex2D('s1', 'r1', 'r0', 's0', labeled_points))
    bk_lines.append(f" innerCircle ({ni} {ns} {nz})")
    bk_lines.append(f" simpleGrading {innerCircle_grading}\n\n"),
    bk_lines.append(');\n\n')

    return bk_lines

def edges(cxs, cys, cxc, cyc, z_in, z_out, labeled_points, arc):

    evert_lines = evert(cxs, cys, cxc, cyc, z_in, z_out)
    ed_lines = ['edges\n','(\n\n','//Circle edges\n']

    ed_lines.append(f"{arc} {labeled_points['r3b']} {labeled_points['r0b']} {evert_lines[0][0]}\n")
    ed_lines.append(f"{arc} {labeled_points['r0b']} {labeled_points['r1b']} {evert_lines[0][1]}\n")
    ed_lines.append(f"{arc} {labeled_points['r1b']} {labeled_points['r2b']} {evert_lines[0][2]}\n")
    ed_lines.append(f"{arc} {labeled_points['r2b']} {labeled_points['r3b']} {evert_lines[0][3]}\n")


    ed_lines.append(f"{arc} {labeled_points['r3t']} {labeled_points['r0t']} {evert_lines[1][0]}\n")
    ed_lines.append(f"{arc} {labeled_points['r0t']} {labeled_points['r1t']} {evert_lines[1][1]}\n")
    ed_lines.append(f"{arc} {labeled_points['r1t']} {labeled_points['r2t']} {evert_lines[1][2]}\n")
    ed_lines.append(f"{arc} {labeled_points['r2t']} {labeled_points['r3t']} {evert_lines[1][3]}\n")
    
    ed_lines.append('\n//Square edges\n')

    ed_lines.append(f"arc {labeled_points['s3b']} {labeled_points['s0b']} {evert_lines[2][0]}\n")
    ed_lines.append(f"arc {labeled_points['s0b']} {labeled_points['s1b']} {evert_lines[2][1]}\n")
    ed_lines.append(f"arc {labeled_points['s1b']} {labeled_points['s2b']} {evert_lines[2][2]}\n")
    ed_lines.append(f"arc {labeled_points['s2b']} {labeled_points['s3b']} {evert_lines[2][3]}\n")


    ed_lines.append(f"arc {labeled_points['s3t']} {labeled_points['s0t']} {evert_lines[3][0]}\n")
    ed_lines.append(f"arc {labeled_points['s0t']} {labeled_points['s1t']} {evert_lines[3][1]}\n")
    ed_lines.append(f"arc {labeled_points['s1t']} {labeled_points['s2t']} {evert_lines[3][2]}\n")
    ed_lines.append(f"arc {labeled_points['s2t']} {labeled_points['s3t']} {evert_lines[3][3]}\n")
    
    ed_lines.append('\n')
    ed_lines.append(');\n\n')
    
    return ed_lines

def btquad(p1, p2, labeled_points):
    return f"({labeled_points[p1+'b']} {labeled_points[p2+'b']} {labeled_points[p2+'t']} {labeled_points[p1+'t']})"

def bottomQuad(p1, p2, p3, p4, labeled_points):
    return f"({labeled_points[p1+'b']} {labeled_points[p2+'b']} {labeled_points[p3+'b']} {labeled_points[p4+'b']})"

def topQuad(p1, p2, p3, p4, labeled_points):
    return f"({labeled_points[p1+'t']} {labeled_points[p4+'t']} {labeled_points[p3+'t']} {labeled_points[p2+'t']})"

def boundaries(labeled_points):

    
    bdry_lines = ['boundary\n',
                  '(\n',
                  '    walls\n',
                  '    {\n',
                  '        type wall;\n',
                  '        faces\n',
                  '        (\n']
    
    bdry_lines.append(f"        {btquad('r0', 'r3', labeled_points)}\n")
    bdry_lines.append(f"        {btquad('r1', 'r0', labeled_points)}\n")
    bdry_lines.append(f"        {btquad('r2', 'r1', labeled_points)}\n")
    bdry_lines.append(f"        {btquad('r3', 'r2', labeled_points)}\n")

    bdry_lines = bdry_lines + ['        );\n\n', 
                               '    }\n',
                               '    inlet\n',
                               '    {\n',
                               '        type cyclic;\n',
                               '        neighbourPatch outlet;\n',
                               '        faces\n',
                               '        (\n']
    
    bdry_lines.append(f"        {bottomQuad('s3', 's0', 's1', 's2', labeled_points)}\n")
    bdry_lines.append(f"        {bottomQuad('s3', 'r3', 'r0', 's0', labeled_points)}\n")
    bdry_lines.append(f"        {bottomQuad('s2', 'r2', 'r3', 's3', labeled_points)}\n")
    bdry_lines.append(f"        {bottomQuad('s1', 'r1', 'r2', 's2', labeled_points)}\n")
    bdry_lines.append(f"        {bottomQuad('s0', 'r0', 'r1', 's1', labeled_points)}\n")

    bdry_lines = bdry_lines + ['        );\n\n', 
                               '    }\n',
                               '    outlet\n',
                               '    {\n',
                               '        type cyclic;\n',
                               '        neighbourPatch inlet;\n',
                               '        faces\n',
                               '        (\n']
    
    bdry_lines.append(f"        {topQuad('s3', 's0', 's1', 's2', labeled_points)}\n")
    bdry_lines.append(f"        {topQuad('s3', 'r3', 'r0', 's0', labeled_points)}\n")
    bdry_lines.append(f"        {topQuad('s2', 'r2', 'r3', 's3', labeled_points)}\n")
    bdry_lines.append(f"        {topQuad('s1', 'r1', 'r2', 's2', labeled_points)}\n")
    bdry_lines.append(f"        {topQuad('s0', 'r0', 'r1', 's1', labeled_points)}\n")

    bdry_lines = bdry_lines + ['        );\n\n',
                               '    }\n',
                               ');\n',
                               '\n',
                               'mergePatchPairs\n',
                               '(\n',
                               ');\n\n']
    
    return bdry_lines

def of_header(scale):

    header = ['/*--------------------------------*- C++ -*----------------------------------*\\\n',
    '| =========                 |                                                 |\n',
    '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
    '|  \\    /   O peration     | Version:  2306                                  |\n',
    '|   \\  /    A nd           | Website:  www.openfoam.com                      |\n',
    '|    \\/     M anipulation  |                                                 |\n',
    '\*---------------------------------------------------------------------------*/\n',
    '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n',
    'FoamFile\n',
    '{\n',
    '    version         2;\n',
    '    format          ascii;\n',
    '    class           dictionary;\n',
    '    object          blockMeshDict;\n',
    '}\n',
    '\n',
    f'scale {scale};\n\n']

    return header


#%%
save_path = './'

#circle radius
r = 1

# number of cells at the side of the inner square
ns = 60

#number of cells between square and circle
ni = 40

#number of cells in height direction (flow direction)
nz = 1

# the parameter ig controls the elongation of cells in between walls and the
# inner square. 
# for example using 0.5 will cause the last cell before the square to be twice as big
# as the one closest to the wall.
ig = 0.5

#inner square half side
s = 0.4

#curvature factor
# this is a constant that multiplies the square half length 
# i.e. the height of the central point in the inner square curvature will be cf*(square half height)
cf = 1.2

#height of cylinder
h = 0.1
#inlet z
zin = 0.0
#outlet z
zout = h + zin

#circle points
xm, ym = (r*np.cos(np.pi/4), r*np.sin(np.pi/4))

#inner square side curvature
sc = cf*s

#inner square x and y positions
xs = [s, -1.0*s, -1.0*s, s]
ys = [-1.0*s, -1.0*s, s, s]

#inner square x and y position of middle curvature
cxs = [sc, 0.0, -1.0*sc, 0.0]
cys = [0, -1.0*sc, 0.0, sc]

#number of points per half quadrant to interpolate to
N = 1
#right sided points, bottom sided points
rs = np.array([r,0]).reshape((1,2))
bs = np.array([0,-r]).reshape((1,2))

#circle y and z positions
xc = [xm, -xm, -xm, xm]
yc = [-ym, -ym, ym, ym]

# this is the specification of the curve blockMesh will use to join the inner square points at each side
# default is the arc of a circle passing by the three points at each square side
ARC = 'arc'

#istarts in [r,0] point and moves in counter clockwise
cxc = np.array([[r, 0, -r, 0]]).T
cyc = np.array([[0, -r, 0, r]]).T

circ_x = [r*np.cos(th) for th in np.linspace(0, 2*np.pi, 100)]
circ_y = [r*np.sin(th) for th in np.linspace(0, 2*np.pi, 100)]

#verification plots
plt.figure()
plt.scatter(xc, yc, label='circle')
plt.scatter(xs, ys, label='square')
plt.scatter(cxc, cyc, label='circle curvature', s=20)
plt.scatter(cxs, cys, label='square curvature')
plt.legend()
plt.gca().set_aspect('equal')


# generating the blockMesh file
verts, labeled_points = vertices(xs, ys, xc, yc, z_in=zin, z_out=zout)
blockMesh = []
blockMesh += of_header(scale=1.0)
blockMesh += verts
blockMesh += blocks(labeled_points, ni=ni, ns=ns, nz=nz, innerCircle_grading=f'({ig} 1 1)')
blockMesh += edges(cxs, cys, cxc, cyc, zin, zout, labeled_points, arc = ARC)
blockMesh += boundaries(labeled_points)


print(f'Saving blockMeshDict to {save_path}')
with open(f'{save_path}/blockMeshDict','w+') as file:
    [file.write(line) for line in blockMesh]

#%%