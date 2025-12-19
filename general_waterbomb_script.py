# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

 # I imported these
import numpy as np
import regionToolset  
from abaqus import session
import re
import os
import subprocess

#regular imports                               
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior



for t in range(1,6):

    step_path = fr"C:\University of Tokyo ESEP\origami\cad\rectangles\step\waterbomb_base_60_{t}.STEP"

    # Extract the filename
    filename = os.path.basename(step_path)

    # Use regex to extract both theta0 and n
    match = re.search(r'waterbomb_base_(\d+)_(\d+)', filename)
    if match:
        theta0 = int(match.group(1))
        n = int(match.group(2))
        print(f"Theta0 extracted: {theta0}")
        print(f"n extracted: {n}")
    else:
        raise ValueError("Theta0 and n could not be extracted from filename.")

    step = mdb.openStep(step_path
        , 
        scaleFromFile=OFF)

    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-1', 
        geometryFile=step, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-2', 
        geometryFile=step, bodyNum=2, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-3', 
        geometryFile=step, bodyNum=3, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-4', 
        geometryFile=step, bodyNum=4, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-5', 
        geometryFile=step, bodyNum=5, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-6', 
        geometryFile=step, bodyNum=6, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-7', 
        geometryFile=step, bodyNum=7, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    mdb.models['Model-1'].PartFromGeometryFile(name=f'waterbomb_base_{theta0}-8', 
        geometryFile=step, bodyNum=8, combine=False, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)



    p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-8']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.111464, 
        farPlane=0.203825, width=0.0725191, height=0.0447918, cameraPosition=(
        -0.0609749, 0.0781987, 0.148167), cameraUpVector=(-0.630002, 0.170846, 
        -0.757568), cameraTarget=(-0.023126, -0.0152819, 0.0258567))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=0.111922, 
        farPlane=0.19948, width=0.0728169, height=0.0449757, cameraPosition=(
        -0.147833, 0.0793698, 0.0487164), cameraUpVector=(0.679878, 0.472986, 
        -0.560402), cameraTarget=(-0.0226395, -0.0152885, 0.0264138))





    # User-defined input values

    # GEOMETRY
    theta_deg = theta0
    theta = np.radians(theta_deg)
    length = 0.1
    breadth = length/n
    diag = np.sqrt(length**2 + breadth**2)



    # STEP and BC
    # STEP
    if n == 1:
        num_steps = 110
    elif n == 2:
        num_steps = 70
    elif n == 3:
        num_steps = 40
    elif n == 4:
        num_steps = 32
    else:
        num_steps = 25
    
    start_step_size = 0.0001
    max_step_size = 0.0001
    min_step_size = 0.0000001
    #BC
    vertical_disp = 5


    # MATERIAL
    E = 93000000000
    v = 0.35



    # Values from MATLAB-based conversion
    alpha = np.arctan(breadth / length)                              

    theta_dash = np.arccos((breadth/length)*np.cos(theta))    # theta for longer valley

    a_coef = (1/np.tan(theta))**2 + (1/np.tan(theta_dash))**2 + 1
    b_coef = 2*( np.cos(theta)*np.sin(alpha)/(np.sin(theta))**2 + np.cos(alpha)*np.cos(theta_dash)/(np.sin(theta_dash))**2)
    c_coef = (np.sin(alpha)/np.sin(theta))**2 + (np.cos(alpha)/np.sin(theta_dash))**2 - 1

    print(f"breadth:{breadth}")
    print(f"alpha:{np.rad2deg(alpha)}")
    print(f"theta:{np.rad2deg(theta)}")
    print(f"theta_dash:{np.rad2deg(theta_dash)}")
    print(f"a_coef:{a_coef}")
    print(f"b_coef:{b_coef}")
    print(f"c_coef:{c_coef}")

    # z = (b_coef + np.sqrt(b_coef**2 - 4 * a_coef * c_coef)) / (2 * a_coef)
    # y = -(np.cos(alpha) - z * np.sin(theta)) / np.cos(theta)

    y = (-b_coef + np.sqrt(b_coef**2 - 4 * a_coef * c_coef)) / (2 * a_coef)
    x = -(y*np.cos(theta) + np.sin(alpha))/np.sin(theta)
    z = (np.cos(alpha) + y*np.cos(theta_dash))/np.sin(theta_dash)
    # mag_dir_vec = 0.055902
    # x = -0.024482/mag_dir_vec
    # y = -0.007596/mag_dir_vec
    # z = 0.049678/mag_dir_vec


    print(f"({x},{y},{z})")


    # 1. Direction vectors for crease edges 1 to 8
    dir_vectors = {
        1: np.array([-np.sin(theta), -np.cos(theta), 0]),
        2: np.array([x, y, z]),
        3: np.array([0, -np.cos(theta_dash), np.sin(theta_dash)]),
        4: np.array([-x, y, z]),
        5: np.array([np.sin(theta), -np.cos(theta), 0]),
        6: np.array([-x, y, -z]),
        7: np.array([0, -np.cos(theta_dash), -np.sin(theta_dash)]),
        8: np.array([x, y, -z])
    }

    # 2. Midpoints of square sides (before folding) — a, b, c, d
    valley_midpoints = {
        'a': (breadth / 2) * dir_vectors[1],
        'b': (length / 2) * dir_vectors[3],
        'c': (breadth / 2) * dir_vectors[5],
        'd': (length / 2) * dir_vectors[7]
    }

    # 3. Midpoints of crease edges 1 to 8
    crease_midpoints = {
        1: (breadth / 4) * dir_vectors[1],
        2: (diag / 4) * dir_vectors[2],
        3: (length / 4) * dir_vectors[3],
        4: (diag / 4) * dir_vectors[4],
        5: (breadth / 4) * dir_vectors[5],
        6: (diag / 4) * dir_vectors[6],
        7: (length / 4) * dir_vectors[7],
        8: (diag / 4) * dir_vectors[8],
        9: (breadth / 4) * dir_vectors[1],
    }

    # 4. Corner coordinates (mountain end points) — p, q, r, s
    corner_coords = {
        'p': (diag / 2) * dir_vectors[2],
        'q': (diag / 2) * dir_vectors[4],
        'r': (diag / 2) * dir_vectors[6],
        's': (diag / 2) * dir_vectors[8]
    }

    # 5. Direction vectors 9 to 16 based on valleys and mountain corners
    def unit_vec(vec):
        return vec / np.linalg.norm(vec)

    dir_vectors.update({
        9:  unit_vec(corner_coords['p'] - valley_midpoints['a']),
        10: unit_vec(corner_coords['p'] - valley_midpoints['b']),
        11: unit_vec(corner_coords['q'] - valley_midpoints['b']),
        12: unit_vec(corner_coords['q'] - valley_midpoints['c']),
        13: unit_vec(corner_coords['r'] - valley_midpoints['c']),
        14: unit_vec(corner_coords['r'] - valley_midpoints['d']),
        15: unit_vec(corner_coords['s'] - valley_midpoints['d']),
        16: unit_vec(corner_coords['s'] - valley_midpoints['a'])
    })

    # 6. Midpoints of outer edges (non-crease) — 9 to 16
    outer_midpoints = {
        9:  (valley_midpoints['a'] + corner_coords['p']) / 2,
        10: (valley_midpoints['b'] + corner_coords['p']) / 2,
        11: (valley_midpoints['b'] + corner_coords['q']) / 2,
        12: (valley_midpoints['c'] + corner_coords['q']) / 2,
        13: (valley_midpoints['c'] + corner_coords['r']) / 2,
        14: (valley_midpoints['d'] + corner_coords['r']) / 2,
        15: (valley_midpoints['d'] + corner_coords['s']) / 2,
        16: (valley_midpoints['a'] + corner_coords['s']) / 2
    }



    ## PARTITIONING EDGES

    print(f"valley_mid_pt = {valley_midpoints['a']}")

    for i in range(8):
        print(f"i = {i}")
        p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-{i+1}']

        # Coordinates for the 3 edges of each triangle
        coord_1 = crease_midpoints[i+1]
        coord_2 = crease_midpoints[i+2]
        coord_3 = outer_midpoints[i+9]
        coords = [coord_1, coord_2, coord_3]

        # Set naming rules
        if i % 2 == 0:
            set_names = ['valley_crease', 'mountain_crease', 'outer_edge']
        else:
            set_names = ['mountain_crease', 'valley_crease', 'outer_edge']

        for j in range(3):
            print(f"j = {j}")
            #coord = coords[j]
            coord = np.round(coords[j], 7)
            set_name = set_names[j]

            # Partition the edge at 0.1 and 0.88
            #edge = p.edges.findAt((coord,))[0]
            edge = p.edges.getClosest((coord,))[0][0]
            p.PartitionEdgeByParam(edges=(edge,), parameter=0.1)
            #edge = p.edges.findAt((coord,))[0]
            edge = p.edges.getClosest((coord,))[0][0]
            p.PartitionEdgeByParam(edges=(edge,), parameter=0.88)

            # Use getClosest first to retrieve the actual point on the edge
            edge = p.edges.getClosest((coord,))[0][0]
            true_coord = edge.pointOn[0]  # Abaqus-trusted coordinate

            # Now use this exact point for findAt
            edge = p.edges.findAt((true_coord,))
            p.Set(edges=edge, name=f"{set_name}")



    # MAKING SETS FOR VALLEY END VERTICE (to be used for BC)

    for i in range(8):
        print(i)
        p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-{i+1}']

        # Choose the coordinate based on index
        if i in [0, 7]:
            key = 'a'
        elif i in [1, 2]:
            key = 'b'
        elif i in [3, 4]:
            key = 'c'
        elif i in [5, 6]:
            key = 'd'

        approx_coord = valley_midpoints[key]

        # Use getClosest to find the actual vertex and extract the trusted coordinate
        closest_vertex = p.vertices.getClosest((approx_coord,))[0][0]
        print(f"approx_coord = {approx_coord}")
        true_coord = closest_vertex.pointOn[0]
        print(f"true_coord = {true_coord}")
        # Now use the trusted coordinate with findAt
        vertex = p.vertices.findAt((true_coord,))
        p.Set(vertices=vertex, name="valley_vertex")





    # MAKING SETS MIDDLE TOP VERTEX ( to be used for BC to push down )

    coord = (0,0,0)

    for i in range(8):
        p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-{i+1}']

        vertex = p.vertices.findAt(((coord),)) 
        p.Set(vertices=vertex, name="centre_top_vertex")



    ## MATERIAL AND SECTION MAKING

    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Elastic(table=((E, v), ))

    mdb.models['Model-1'].HomogeneousShellSection(name='Section-1', 
    preIntegrate=OFF, material='Material-1', thicknessType=UNIFORM, 
    thickness=0.003, thicknessField='', nodalThicknessField='', 
    idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
    thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
    integrationRule=SIMPSON, numIntPts=5)






    # SECTION ASSIGNMENT

    for i in range(8):
        p = mdb.models['Model-1'].parts[f'waterbomb_base_{theta0}-{i+1}']

        # Select all faces in the part (this works because your part is a shell with one face)
        faces = p.faces[:]

        # Create a set from these faces
        region = p.Set(faces=faces, name='Set-1')

        # Assign section to this region
        p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)



    ## STEP

    mdb.models['Model-1'].StaticRiksStep(
        name='Step-1',
        previous='Initial',
        nlgeom=ON,  # Enable geometric nonlinearity
        initialArcInc = start_step_size,
        minArcInc = min_step_size,
        maxArcInc = max_step_size,
        maxNumInc = num_steps
    )



    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')



    ## MESH + INSTANCE CREATION

    element_size = 0.004

    total_elements = 0
    total_nodes = 0

    # Reference to model and root assembly
    model = mdb.models['Model-1']
    a = model.rootAssembly

    # Step 1: Mesh all parts before creating instances
    for i in range(8):
        part_name = f'waterbomb_base_{theta0}-{i+1}'
        p = model.parts[part_name]
        
        # Seeding and meshing
        p.seedPart(size=element_size, deviationFactor=0.1, minSizeFactor=0.1)
        p.generateMesh()
        
        # Count elements and nodes
        num_elements = len(p.elements)
        num_nodes = len(p.nodes)
        
        total_elements += num_elements
        total_nodes += num_nodes

        print(f"Part {part_name}: Elements = {num_elements}, Nodes = {num_nodes}")

    print(f"\n Total Elements = {total_elements}")
    print(f"Total Nodes = {total_nodes}")

    # Step 2: Now create independent instances (must come after meshing)
    a.DatumCsysByDefault(CARTESIAN)

    for i in range(8):
        part_name = f'waterbomb_base_{theta0}-{i+1}'
        inst_name = f'{part_name}-1'
        p = model.parts[part_name]
        a.Instance(name=inst_name, part=p, dependent=ON)



    ## INTERACTION

    # Reference point

    for i in range(8):
        part_name = f'waterbomb_base_{theta0}-{i+1}'
        p = mdb.models['Model-1'].parts[part_name]
        a = mdb.models['Model-1'].rootAssembly
        inst = a.instances[f'{part_name}-1']

        # Only get edge sets if needed
        if i % 2 == 0:
            edge_types = ['valley_crease', 'mountain_crease', 'outer_edge']
        else:
            edge_types = ['outer_edge']

        for j, edge_type in enumerate(edge_types):
            edge_set = p.sets[edge_type]
            edge = edge_set.edges[0]  # Assumes single edge per set

            v1, v2 = edge.getVertices()
            # coord1 = p.vertices[v1].pointOn[0]
            # coord2 = p.vertices[v2].pointOn[0]

            coord1 = inst.vertices[v1].pointOn[0]
            coord2 = inst.vertices[v2].pointOn[0]


            rp1 = a.ReferencePoint(point=coord1)
            rp2 = a.ReferencePoint(point=coord2)

            if edge_type == 'valley_crease':
                if i in [4]:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_{i+1}_upper')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_{i+1}_lower')
                else:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_{i+1}_lower')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_{i+1}_upper')
            elif edge_type == 'mountain_crease':
                if i in [4]:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_{i+2}_lower')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_{i+2}_upper')
                else:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_{i+2}_upper')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_{i+2}_lower')
            elif edge_type == 'outer_edge':
                if i in [0, 1, 2, 6]:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_outer_edge_{2*i+2}')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_outer_edge_{2*i+1}')
                else:
                    a.Set(referencePoints=(a.referencePoints[rp1.id],), name=f'rp_outer_edge_{2*i+1}')
                    a.Set(referencePoints=(a.referencePoints[rp2.id],), name=f'rp_outer_edge_{2*i+2}')



    # MAKE WIRE

    a = mdb.models['Model-1'].rootAssembly

    for i in range(1, 9):
        # Access the reference points by their set names
        rp_upper_set = a.sets[f'rp_{i}_upper']
        rp_lower_set = a.sets[f'rp_{i}_lower']

        # Get the actual reference point objects
        rp_upper = rp_upper_set.referencePoints[0]
        rp_lower = rp_lower_set.referencePoints[0]

        # Create the wire between these two reference points
        wire = a.WirePolyLine(points=((rp_upper, rp_lower),), mergeType=IMPRINT, meshable=OFF)

    print("Wire Made")

    # Creating sets of wire

    for i in range(1, 9):
        # Use the known midpoint coordinate of the crease (edge i)
        midpoint_coord = crease_midpoints[i]

        # Find the edge created by WirePolyLine using the midpoint
        edge = a.edges.findAt((midpoint_coord,))

        # Create a set for the edge
        a.Set(edges=(edge,), name=f'Wire_{i}_{theta0}_{n}')

    print("Sets of wire made")

    # MAKING CONNECTOR SECTION

    spring_stiffness = 1

    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    mdb.models['Model-1'].ConnectorSection(name=f'ConnSect_{theta0}_{n}', assembledType=HINGE)
    elastic_0 = connectorBehavior.ConnectorElasticity(components=(4, ), table=((
        spring_stiffness, ), ))
    mdb.models['Model-1'].sections[f'ConnSect_{theta0}_{n}'].setValues(behaviorOptions =(
        elastic_0, ) )
    mdb.models['Model-1'].sections[f'ConnSect_{theta0}_{n}'].behaviorOptions[0].ConnectorOptions(
        )

    print("connector section made")

    # ASSIGN SECTION + ORIENTATION FOR 8 WIRES
    a = mdb.models['Model-1'].rootAssembly
    model = mdb.models['Model-1']

    # Dictionary to store datum IDs by wire name
    datum_cs_dict = {}

    for i in range(1, 9):
        print(i)

        # Step 1: Create Datum CSYS
        try:
            origin_rp = a.sets[f'rp_{i}_lower'].referencePoints[0]
            point1_rp = a.sets[f'rp_{i}_upper'].referencePoints[0]

            if i == 8:
                point2_rp = a.sets['rp_1_lower'].referencePoints[0]
            else:
                point2_rp = a.sets[f'rp_{i+1}_lower'].referencePoints[0]

            datum = a.DatumCsysByThreePoints(
                origin=origin_rp,
                point1=point1_rp,
                point2=point2_rp,
                name=f'Datum_csys_w{i}_{theta0}_{n}',  # GUI label only
                coordSysType=CARTESIAN
            )
            datum_id = datum.id
            datum_cs_dict[f'Datum_csys_w{i}_{theta0}_{n}'] = datum_id  #  Save the ID for future reference

        except KeyError as e:
            print(f"Reference point missing for wire {i}: {e}")
            continue

        # Step 2: Get wire set
        wire_set_name = f'Wire_{i}_{theta0}_{n}'
        if wire_set_name not in a.sets:
            print(f"Set {wire_set_name} not found.")
            continue
        region = a.sets[wire_set_name]

        # Step 3: Assign Section
        csa = a.SectionAssignment(sectionName=f'ConnSect_{theta0}_{n}', region=region)

        # Step 4: Assign Connector Orientation
        model.rootAssembly.ConnectorOrientation(
            region=csa.getSet(),
            localCsys1=model.rootAssembly.datums[datum_id]
        )

    # The dictionary datum_cs_dict can now be reused in coupling constraints or elsewhere.

    print("section_assigned")


    # COUPLING CONSTRAINTS USING LOCAL DATUM CSYS


    a = mdb.models['Model-1'].rootAssembly
    model = mdb.models['Model-1']

    for i in range(8):
        edge_type = 'valley_crease' if i % 2 == 0 else 'mountain_crease'

        # Reference points
        try:
            rp_upper = a.sets[f'rp_{i+1}_upper'].referencePoints[0]
            rp_lower = a.sets[f'rp_{i+1}_lower'].referencePoints[0]
        except KeyError as e:
            print(f"Missing reference point for i={i+1}: {e}")
            continue

        # Datum CSYS to use
        datum_name = f'Datum_csys_w{i+1}_{theta0}_{n}'
        if datum_name not in datum_cs_dict:
            print(f"Missing datum CSYS: {datum_name}")
            continue
        datum_id = datum_cs_dict[datum_name]

        # Wraparound case for i=0
        if i == 0:
            inst1 = a.instances[f'waterbomb_base_{theta0}-8-1']
            inst2 = a.instances[f'waterbomb_base_{theta0}-1-1']
        else:
            inst1 = a.instances[f'waterbomb_base_{theta0}-{i}-1']
            inst2 = a.instances[f'waterbomb_base_{theta0}-{i+1}-1']

        try:
            # Region from sets
            region_ctrl_upper = regionToolset.Region(referencePoints=(rp_upper,))
            region_ctrl_lower = regionToolset.Region(referencePoints=(rp_lower,))

            region_surf_1 = inst1.sets[edge_type]
            region_surf_2 = inst2.sets[edge_type]

            # Coupling 1 (upper)
            model.Coupling(
                name=f'Coupling_{i}_upper',
                controlPoint=region_ctrl_upper,
                surface=region_surf_1,
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
                rotationalCouplingType=ROTATIONAL_STRUCTURAL,
                weightingMethod=UNIFORM,
                localCsys=model.rootAssembly.datums[datum_id],
                u1=ON, u2=ON, u3=ON,
                ur1=ON, ur2=ON, ur3=ON
            )

            # Coupling 2 (lower)
            model.Coupling(
                name=f'Coupling_{i}_lower',
                controlPoint=region_ctrl_lower,
                surface=region_surf_2,
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
                rotationalCouplingType=ROTATIONAL_STRUCTURAL,
                weightingMethod=UNIFORM,
                localCsys=model.rootAssembly.datums[datum_id],
                u1=ON, u2=ON, u3=ON,
                ur1=ON, ur2=ON, ur3=ON
            )

        except KeyError as e:
            print(f"Missing set for coupling at i={i}: {e}")




    ## LOAD and BC

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)

    a = mdb.models['Model-1'].rootAssembly

    for i in range(8):


        region1 = a.instances[f'waterbomb_base_{theta0}-{i+1}-1'].sets['valley_vertex']
        region2 = a.instances[f'waterbomb_base_{theta0}-{i+1}-1'].sets['centre_top_vertex']

        mdb.models['Model-1'].DisplacementBC(
            name=f'Tip-{i+1}',
            createStepName='Step-1',
            region=region2,
            u1=0.0, u2=-vertical_disp, u3=0.0,
            ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=UNSET,
            fixed=OFF,
            distributionType=UNIFORM,
            fieldName='',
            localCsys=None)

        if i in [1,2,5,6]:

            mdb.models['Model-1'].DisplacementBC(
                name=f'BC_y-{i+1}',
                createStepName='Step-1',
                region=region1,
                u1=0.0, u2=0.0, u3=UNSET,
                ur1=UNSET, ur2=UNSET, ur3=UNSET,
                amplitude=UNSET,
                fixed=OFF,
                distributionType=UNIFORM,
                fieldName='',
                localCsys=None)
            
        else:

            mdb.models['Model-1'].DisplacementBC(
                name=f'BC_x-{i+1}',
                createStepName='Step-1',
                region=region1,
                u1=UNSET, u2=0.0, u3=0.0,
                ur1=UNSET, ur2=UNSET, ur3=UNSET,
                amplitude=UNSET,
                fixed=OFF,
                distributionType=UNIFORM,
                fieldName='',
                localCsys=None)




    ## HISTORY OUTPUT

    # TIP DISPLACEMENT

    regionDef = a.instances[f'waterbomb_base_{theta0}-{i+1}-1'].sets['centre_top_vertex']

    mdb.models['Model-1'].HistoryOutputRequest(name='Vertical_disp_of_tip', 
        createStepName='Step-1', variables=('U2', ), region=regionDef, 
        sectionPoints=DEFAULT, rebar=EXCLUDE)


    # REACTION FORCE

    for i in range(8):

        regionDef = a.instances[f'waterbomb_base_{theta0}-{i+1}-1'].sets['centre_top_vertex']

        mdb.models['Model-1'].HistoryOutputRequest(name=f'RF_{i+1}', 
            createStepName='Step-1', variables=('RF2', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)




    ## JOB

    mdb.Job(name=f'Waterbomb_base_Job_{theta0}_{n}', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB)

    # SUBMIT JOB

    job_name = f'Waterbomb_base_Job_{theta0}_{n}'
    mdb.jobs[job_name].submit(consistencyChecking=ON)
    mdb.jobs[job_name].waitForCompletion()




    ## RESULTS AND POST PROCESSING



    # # Path to your ODB
    # # Now open ODB and run post-processing
    # odb_path = fr"C:\University of Tokyo ESEP\Single_cell_simulations\Waterbomb_base_Job_{theta0}_{n}.odb"
    # odb = session.openOdb(name=odb_path)
    # session.viewports['Viewport: 1'].setValues(displayedObject=odb)

    # # Store extracted data
    # all_curves = []

    # # Step 1: Loop over all 8 crests
    # for i in range(1,9):
    #     # Prefix to match
    #     expected_prefix = f'Node WATERBOMB_BASE_{theta0}_{n}-{i}-1'
    #     # Node WATERBOMB_BASE_60-1-1.1
    #     #expected_prefix = 'Reaction force: RF2 PI: WATERBOMB_BASE_60-1-1 Node 1 in NSET CENTRE_TOP_VERTEX'

    #     #print(odb.steps['Step-1'].historyRegions.keys())

    #     # Step 2: Find matching key
    #     matching_key = None
    #     for key in odb.steps['Step-1'].historyRegions.keys():
    #         if expected_prefix in key:
    #             matching_key = key
    #             break

    #     if matching_key is None:
    #         print(f"  No match found for crest {i+1}")
    #         continue

    #     # Step 3: Get RF2 data from that history region
    #     region = odb.steps['Step-1'].historyRegions[matching_key]
    #     if 'RF2' not in region.historyOutputs:
    #         print(f"  RF2 not found in {matching_key}")
    #         continue

    #     rf2_data = region.historyOutputs['RF2'].data
    #     xy_data = session.XYData(data=rf2_data, name=f'RF2_{i}')
        
    #     # Step 4: Store for plotting
    #     curve = session.Curve(xyData=xy_data)
    #     all_curves.append(curve)

    # ## OTHER 3 RESULTS

    # # Open the ODB
    # odb_path = fr"C:\University of Tokyo ESEP\Single_cell_simulations\Waterbomb_base_Job_{theta0}_{n}.odb"
    # odb = session.openOdb(name=odb_path)

    # # Extract strain energy for whole model
    # se_data = session.XYDataFromHistory(
    #     name='strain_energy',
    #     odb=odb,
    #     outputVariableName='Strain energy: ALLSE for Whole Model',
    #     __linkedVpName__='Viewport: 1'
    # )

    # # Extract external work for whole model
    # wk_data = session.XYDataFromHistory(
    #     name='work_done_external',
    #     odb=odb,
    #     outputVariableName='External work: ALLWK for Whole Model',
    #     __linkedVpName__='Viewport: 1'
    # )

    # # Extract vertical displacement (U2) of the center top vertex
    # # Loop to find the matching key (like you did for RF2)
    # disp_key = None
    # expected_disp_prefix = f'Node WATERBOMB_BASE_{theta0}_{n}-1-1'
    # for key in odb.steps['Step-1'].historyRegions.keys():
    #     if expected_disp_prefix in key:
    #         disp_key = key
    #         break

    # if disp_key is not None and 'U2' in odb.steps['Step-1'].historyRegions[disp_key].historyOutputs:
    #     u2_data_raw = odb.steps['Step-1'].historyRegions[disp_key].historyOutputs['U2'].data
    #     u2_data = session.XYData(data=u2_data_raw, name='u2_disp')
    # else:
    #     print("Could not find U2 displacement data")
    #     u2_data = None

    # # Combine all curves (RFs already in all_curves)
    # if u2_data:
    #     all_curves.extend([
    #         session.Curve(xyData=se_data),
    #         session.Curve(xyData=wk_data),
    #         session.Curve(xyData=u2_data)
    #     ])
    # else:
    #     all_curves.extend([
    #         session.Curve(xyData=se_data),
    #         session.Curve(xyData=wk_data)
    #     ])

    # ## SAVE RESULTS
    # # Ensure XY Report output is in separate tables

    # report_path = f'abaqus_{theta0}_{n}.rpt'
    # if os.path.exists(report_path):
    #     with open(report_path, 'w') as f:
    #         pass  # truncate the file

    # session.xyReportOptions.setValues(layout=SEPARATE_TABLES)

    # # Collect all required result names
    # xy_names = [f'RF2_{i}' for i in range(1, 9)]  # RF2_Crest_1 to RF2_Crest_8
    # xy_names += ['strain_energy', 'work_done_external', 'u2_disp']

    # # Fetch xyData objects from session using names
    # xy_data_objects = []
    # for name in xy_names:
    #     if name in session.xyDataObjects:
    #         xy_data_objects.append(session.xyDataObjects[name])
    #     else:
    #         print(f"Warning: {name} not found in xyDataObjects")

    # # Write to report
    # session.writeXYReport(fileName=f'abaqus_{theta0}_{n}.rpt', xyData=xy_data_objects)


    # Path to your ODB
    odb_path = fr"C:\University of Tokyo ESEP\Single_cell_simulations\Waterbomb_base_Job_{theta0}_{n}.odb"
    odb = session.openOdb(name=odb_path)
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)

    # Store extracted data
    all_curves = []

    # Step 1: Loop over all 8 crests to extract RF2
    for i in range(1, 9):
        expected_prefix = f'Node WATERBOMB_BASE_{theta0}-{i}-1'

        matching_key = None
        for key in odb.steps['Step-1'].historyRegions.keys():
            if expected_prefix in key:
                matching_key = key
                break

        if matching_key is None:
            print(f"  No match found for crest {i}")
            continue
        
        print(matching_key)
        region = odb.steps['Step-1'].historyRegions[matching_key]
        if 'RF2' not in region.historyOutputs:
            print(f"  RF2 not found in {matching_key}")
            continue

        rf2_data = region.historyOutputs['RF2'].data
        xy_data = session.XYData(data=rf2_data, name=f'RF2_{i}')
        #session.xyDataObjects[f'RF2_{i}'] = xy_data  # Register for report
        all_curves.append(session.Curve(xyData=xy_data))

    # Step 2: Strain energy and external work
    se_data = session.XYDataFromHistory(
        name='strain_energy',
        odb=odb,
        outputVariableName='Strain energy: ALLSE for Whole Model',
        __linkedVpName__='Viewport: 1'
    )
    wk_data = session.XYDataFromHistory(
        name='work_done_external',
        odb=odb,
        outputVariableName='External work: ALLWK for Whole Model',
        __linkedVpName__='Viewport: 1'
    )
    all_curves.extend([
        session.Curve(xyData=se_data),
        session.Curve(xyData=wk_data)
    ])

    # Step 3: Vertical displacement (U2)
    disp_key = None
    expected_disp_prefix = f'Node WATERBOMB_BASE_{theta0}-8-1'

    print("\n🔍 Available history region keys:")
    for key in odb.steps['Step-1'].historyRegions.keys():
        print(f"  {key}")
        if expected_disp_prefix in key:
            disp_key = key  # Save it if it matches
            print(f"✅ Found potential match: {disp_key}")

    if disp_key is not None:
        region = odb.steps['Step-1'].historyRegions[disp_key]
        print(f"\n📦 Outputs available in selected region:\n  {list(region.historyOutputs.keys())}")

        if 'U2' in region.historyOutputs:
            u2_data_raw = region.historyOutputs['U2'].data
            u2_data = session.XYData(data=u2_data_raw, name='u2_disp')
            all_curves.append(session.Curve(xyData=u2_data))
        else:
            print("⚠️  'U2' not found in selected region outputs.")
            u2_data = None
    else:
        print("❌ No matching key found for center top vertex U2.")
        u2_data = None


    session.xyReportOptions.setValues(layout=SEPARATE_TABLES)

    # Prepare list of xyDataObjects to write
    xy_names = [f'RF2_{i}' for i in range(1, 9)]
    xy_names += ['strain_energy', 'work_done_external']
    if u2_data is not None:
        xy_names.append('u2_disp')

    xy_data_objects = []
    for name in xy_names:
        if name in session.xyDataObjects:
            xy_data_objects.append(session.xyDataObjects[name])
        else:
            print(f" Warning: {name} not found in xyDataObjects")

    # Final report
    report_path = f'abaqus_{theta0}_{n}.rpt'
    session.writeXYReport(fileName=report_path, xyData=xy_data_objects)
    print(f"Report written to {report_path}")


    
