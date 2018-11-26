import bpy
from math import (
        atan, acos, cos,
        sin, tan, pi,
        radians, sqrt
        )
from mathutils import (
        geometry, Vector,
        )

#Initial assigments

def cart2pol(v):
    v=Vector((v[0],v[1]))
    theta = atan(v.y/v.x)
    rho = v.length
    return theta,rho

def pol2cart(theta,rho, z):
    x = rho*cos(theta)
    y = rho*sin(theta)
    return Vector((x, y, z))

def calc_max_beta(a):
    beta = 0
    for I in range(9):
        inc = 10**-I
        while abs(a) > tan(beta)-beta:
            beta += inc
        beta -= inc
    return beta

#Start calculations

#add_tooth(t, d, radius, Ad, De, base, p_angle, crown=0.0, res=64)
#Calculates the shape of one tooth, result in polar coordinates
#t      Tooth width
#d      Z axis value
def add_tooth(t, radius, Ra, Rd, base, p_angle, res):
    p_angle= radians(p_angle)
    rinv = radius*cos(p_angle)
    
    k = -t/4-(tan(p_angle)-p_angle)
    
    max_beta = calc_max_beta(k)
    beta = acos(rinv/Ra)
    if beta > max_beta:
        beta = max_beta
    
    A = [k+I*tan(beta)/res for I in range(res+1)]

    verts = []
    verts_polar = []
    jj = -1
    for ii in A:
        x = rinv*(cos(ii)+(ii-k)*sin(ii))
        y = rinv*(sin(ii)-(ii-k)*cos(ii))
        verts.append((x,y))
        
    verts_polar = [(cart2pol(verts[I])) for I in range(len(verts))]

    verts_polar += [(-verts_polar[I][0],verts_polar[I][1])for I in reversed(range(len(verts_polar)))]
    return verts_polar
    
def add_gear(teethNum, Dp, Ad, De, base, p_angle, t_res, r_res, 
             width=1, skew=0, conangle=0, crown=0.0):
    # Initial calculations
    t=2*pi/teethNum
    radius = teethNum/Dp/2
    Ra = radius + Ad
    Rd = radius - De
    Rb = Rd - base
    
    # Generate vertex for single tooth (in polar coordinates)
    verts_pol = add_tooth(t, radius, Ra, Rd, base, p_angle, t_res)
    verts_pol = [(verts_pol[0][0], Rd)] + verts_pol + [(verts_pol[-1][0], Rd)]
    
    # Store the number of vertex per tooth
    tooth_vert_cnt = len(verts_pol)
    
    # Generate vertex for all teeth
    theta = (verts_pol[-1][0]-verts_pol[0][0])/(r_res+1)
    beta = (t-(verts_pol[-1][0]-verts_pol[0][0]))/(r_res+1)
    ring_verts_pol=[]
    verts=[]
    for k in [t*I for I in range(teethNum)]:
        ver_pol = [(verts_pol[I][0]+k, verts_pol[I][1]) for I in range(len(verts_pol))]
        ring_verts_pol.extend([(ver_pol[0][0]+theta*I, Rd) for I in range(1,r_res+1)])
        ring_verts_pol.extend([(ver_pol[-1][0]+beta*I, Rd) for I in range(1,r_res+1)])
        
        verts.extend(pol2cart(ver_pol[I][0], ver_pol[I][1], width) for I in range(len(ver_pol)))
    
    # Generate inner vertex
    ring_verts=[]
    ring_verts.extend(
                [pol2cart(ring_verts_pol[I][0], Rd, width) for I in range(len(ring_verts_pol))]
                )
    ring_verts.extend(
                reversed(
                [pol2cart(ring_verts_pol[I][0], Rd, -width) for I in range(len(ring_verts_pol))]
                ))
    # Store the number of vertex of the deddendum circle
    Rd_circ_cnt = len(ring_verts)
    theta = t/2/(r_res+1)
    alpha = -theta*((r_res+1)/2)
    ring_verts.extend(
                [pol2cart(alpha + theta*I, Rb, width)for I in range((r_res+1)*teethNum*2)]
                )
    ring_verts.extend(
                reversed(
                [pol2cart(alpha + theta*I, Rb, -width)for I in range((r_res+1)*teethNum*2)]
                ))
    # Store the number of vertices of the base circle
    Rb_circ_cnt = len(ring_verts) - Rd_circ_cnt
                
    # Store the number of vertex of all teeth (one side only)
    teeth_vert_cnt = len(verts)
    
    verts.extend(reversed([(verts[I][0], verts[I][1], -verts[I][2])for I in range(len(verts))]))
    verts.extend(ring_verts)

# Create faces
  # Create bottom and upper faces 
    # for each tooth 
    faces=[]
    for I in range(teethNum*2):
        if I < teethNum:
            faces.append(list(range(I*tooth_vert_cnt,(I+1)*tooth_vert_cnt)) + 
                        [teeth_vert_cnt*2 + r_res*2*I + J for J in reversed(range(r_res))])
        else:
            faces.append(list(range(I*tooth_vert_cnt,(I+1)*tooth_vert_cnt)) + 
                        [teeth_vert_cnt*2 + r_res*2*I+J for J in reversed(range(r_res, 2*r_res))])
    
    # for th ring
    state = 0
    J = 0
    K = 0
    strt=teeth_vert_cnt*2
    for I in range(teethNum*(r_res+1)*4):
        if I == teethNum*(r_res+1)*2-1:
            if r_res != 0:
                J += r_res
                faces.append([strt+Rd_circ_cnt/2-1, 0, strt+Rd_circ_cnt, strt+Rd_circ_cnt+I])
                state = 0
            else:
                faces.append([teeth_vert_cnt-1, 0, strt, strt+Rd_circ_cnt+I])
                K+=1
        elif I == teethNum*(r_res+1)*4-1:
            if r_res != 0:
                faces.append([strt-1, strt+Rd_circ_cnt/2, strt+Rd_circ_cnt+Rb_circ_cnt/2, len(verts)-1])
            else:
                faces.append([strt-1, teeth_vert_cnt, strt+Rd_circ_cnt+Rb_circ_cnt/2, len(verts)-1 ])
        elif r_res == 0:
            faces.append([tooth_vert_cnt*((K+1)//2)-(K%2), tooth_vert_cnt*((K+2)//2)-((K-1)%2), strt+Rd_circ_cnt+I+1, strt+Rd_circ_cnt+I])
            K += 1
        elif state == 0:
            faces.append([tooth_vert_cnt*((K+1)//2)-(K%2), strt+J, strt+Rd_circ_cnt+I+1, strt+Rd_circ_cnt+I])
            K += 1
        elif state == r_res:
            J += r_res
            faces.append([strt+J-1, tooth_vert_cnt*((K+1)//2)-(K%2), strt+Rd_circ_cnt+I+1, strt+Rd_circ_cnt+I])
            state = -1
            
        else:
            faces.append([strt+J+state-1, strt+J+state, strt+Rd_circ_cnt+I+1, strt+Rd_circ_cnt+I])
            
        state += 1
    
  # Create side faces
    # Outer faces
    for I in range(teethNum):
        for J in range(tooth_vert_cnt-1):
            K= I * tooth_vert_cnt +J
            faces.append([strt-K-1, strt-K-2, K+1, K])
        state = 0
        for J in range(r_res+1):
            if I == teethNum-1 and J == r_res:
                if r_res != 0:
                    faces.append([strt-1, 0, strt-1+(I+1)*2*r_res, strt+Rd_circ_cnt-(I+1)*2*r_res])
                else:
                    faces.append([strt-1, 0, teeth_vert_cnt-1, teeth_vert_cnt])
            elif state == 0:
                if r_res != 0:
                    faces.append([tooth_vert_cnt*(I+1)-1, strt-tooth_vert_cnt*(I+1), strt+Rd_circ_cnt-(I*2+1)*r_res-1, strt+(I*2+1)*r_res])
                else:
                    faces.append([tooth_vert_cnt*(I+1)-1, strt-tooth_vert_cnt*(I+1), strt-tooth_vert_cnt*(I+1)-1, tooth_vert_cnt*(I+1)])
            elif state == r_res:
                faces.append([strt-tooth_vert_cnt*(I+1)-1, tooth_vert_cnt*(I+1), strt-1+(I+1)*2*r_res, strt+Rd_circ_cnt-(I+1)*2*r_res])
            else:
                faces.append([strt+state+(I*2+1)*r_res, strt+state-1+(I*2+1)*r_res, strt+Rd_circ_cnt-(I*2+1)*r_res-state, strt+Rd_circ_cnt-(I*2+1)*r_res-state-1])
            state += 1
            
    #Create inner faces
    strt = strt+Rd_circ_cnt
    for I in range(int(Rb_circ_cnt/2)):
        a=1
        #faces.append([strt+state+(I*2+1)*points, 0,1,2])
    
    return verts,faces

    
#Call function
verts, faces = add_gear(12, 6, 0.1, 0.1, 0.2, 20, 8, 2)

edges = []
name="New Object Mesh"
me = bpy.data.meshes.new(name)
me.from_pydata(verts, edges, faces)

ob = bpy.data.objects.new(name, me)
ob.location = (0,0,0)
ob.show_name = False

# Link object to scene and make active
scn = bpy.context.scene
scn.objects.link(ob)