
from base_io import *
from polygons import *
from schism_structure import *
from schism_source import *
from schism_input import *
from schism_mesh import *
from gr3 import *
import yaml
#import osgeo.ogr, osgeo.osr
from pyproj import Proj, transform
from scipy.interpolate import interp1d
import datetime
import numpy as np
import re
from math import *
import copy
import os
import sys
class SchismSetup(object):
    """ A class to manage SCHISM input data

    """
    def __init__(self):
        """ Constructor
        """
        self._input = SchismInput()
        self._verbose = 1

    @property
    def input(self):
        return self._input

    @property
    def mesh(self):
        """ Mesh

            :getter: Return the mesh.
        """
        return self._input.mesh

    @mesh.setter
    def mesh(self, value):
        """ Mesh setter
        """
        self._input.mesh = value

    def load(self, dir):
        """ Read SCHISM input

            Parameters
            ----------
            dir: String
                An directory where the SCHISM inputs reside
        """
        # Assuming we are dealing SCHISM format
        # Read the grid
        gr3_reader = Gr3IO()
        fname = os.path.join(dir, "hgrid.gr3")
        if not os.path.exists(fname):
            raise Exception("A grid file is not found!")
        self._input.mesh = gr3_reader.read(fname)
        # Read hydraulics if exists
        fname = os.path.join(dir, "hydraulics.in")
        if os.path.exists(fname):
            structure_reader = SchismStructureIO(self._input)
            structure_reader.read(fname)

        # Read source if exists
        fname = os.path.join(dir, "source_sink.in")
        if os.path.exists(fname):
            source_reader = SchismSourceIO(self._input)
            source_reader.read(fname)

    def get_closest_node_i_from_new_mesh(self, mesh, node_i):
        old_node = self._input.mesh.nodes[node_i,]
        return mesh.find_closest_nodes(old_node[:2], count = 1, boundary = 1)

    def adopt_new_mesh(self, fname):
        """ Adopt a new grid.
                fname = the file name of a new hgrid file
        """
        gr3_reader = Gr3IO()
        if not os.path.exists(fname):
            raise Exception("The given mesh file is not found")
        new_mesh = gr3_reader.read(fname, 1)

        mesh = self._input.mesh
        ## Boundary build
        # First delete boundary information of the new mesh
        new_mesh.clear_boundary()
        if mesh.n_boundaries == 0:
            raise Exception("No boundary is in the original mesh.")
        # Open boundary
        for b in mesh.boundaries:
            # First Node
            node_1 = self.get_closest_node_i_from_new_mesh(new_mesh, b.nodes[0])
            # End Node
            node_2 = self.get_closest_node_i_from_new_mesh(new_mesh, b.nodes[-1])
            # Get a path between the two
            new_nodes = new_mesh.shortest_path(node_1, node_2, boundary_only = 1)
            new_mesh.add_boundary(new_nodes, b.btype)

        ## Structure adoption
        new_structures = []
        for struct in self._input.structures:
            # Ref node
            new_ref_pair = tuple(self.get_closest_node_i_from_new_mesh(
                new_mesh, node_i) for node_i in struct.reference_pair)

            # Node pairs.  This one is slightly trickier
            original_up = [node_i[0] for node_i in struct.node_pairs]
            original_down = [node_i[1] for node_i in struct.node_pairs]
            # Step 1: Create a new path of the upstream side
            node_prev = None
            new_up = None
            for node in original_up:
                if node_prev is None:
                    node_prev = node
                else:
                    path = list(new_mesh.shortest_path(node_prev, node))
                    if new_up is None:
                        new_up = path
                    else:
                        if new_up[-1] == path[0]:
                            new_up.extend(path[1:])
                        else:
                            raise Exception('Path segment does not match.')
                    node_prev = node
            # Step 2: Create a new downstream node pair
            new_down_inter = []
            for i in range(len(original_up)):
                path = list(new_mesh.shortest_path(original_up[i],
                                                   original_down[i]))
                new_down_inter.append(path[1])
            # Step 3: Build middle nodes of the downstream side
            node_prev = None
            new_down = None
            for node in new_down_inter:
                if node_prev is None:
                    node_prev = node
                else:
                    path = list(new_mesh.shortest_path(node_prev, node))
                    if new_down is None:
                        new_down = path
                    else:
                        if new_down[-1] == path[0]:
                            new_down.extend(path[1:])
                        else:
                            raise Exception('Path segment does not match.')
                    node_prev = node
            # Step 4: Check the pairs
            new_pairs = []
            if not len(new_up) == len(new_down):
                raise Exception('Number of pairs does not match')
            for i in range(len(new_up)):
                path = list(new_mesh.shortest_path(new_up[i], new_down[i]))
                if not len(path) == 2:
                    raise Exception('Having trouble with ')
                new_pairs.append((new_up[i], new_down[i]))


            # Create a new structure
            struct_new = copy.deepcopy(struct)
            struct_new.reference_pair = new_ref_pair
            struct_new.node_pairs = new_pairs
            new_structures.append(struct_new)

        ## Source/sink adoption
        # TODO: This is tricky... I'll think more about this

        # Set the new mesh
        self._input.mesh = new_mesh
        # Set new structures
        self._structures = new_structures

    def write(self, dir, grid_fname = "hgrid.gr3",
              struct_fname = "hydraulics.in",
              source_fname = "source_sink.in"):
        """ Write SCHISM output files into the given directory.
        """

        # Grid output
        gr3_writer = Gr3IO()
        if not os.path.exists(dir):
            os.makedirs(dir)
        fname = os.path.join(dir, grid_fname)
        gr3_writer.write(self._input.mesh, fname)
        # Structure output
        struct_writer = SchismStructureIO(self._input)
        fname = os.path.join(dir, struct_fname)
        struct_writer.write(fname)
        # Source sink output
        source_writer = SchismSourceIO(self._input)
        fname = os.path.join(dir, source_fname)
        source_writer.write(fname)


    def write_hgrid(self, fname, attr = None, boundary = True):
        """ Write a hgrid file only.
            fname = the output file name
            attr = attribute array in numpy format. If None, original data
            will be kept.
            boundary = If True, boundary info will be appended at the end of
            the file. Otherwise, not appended.
        """
        gr3_writer = Gr3IO()
        gr3_writer.write(self._input.mesh, fname, attr, boundary)


    def write_hgrid_ll(self, fname, input_epsg=26910, output_epsg=4269):
        """ Write a hgrid.ll, lat-long mesh file, of the current mesh.

            Parameters
            ----------
            fname: str
                the output file name
            input_epsg: int, optional
                input EPSG. default value is 26310, NAD83/UTM10N.
            output_epsg: int, optional
                output EPSG. default value is 4269, NAD83
        """
        
	#inSpatialRef = osgeo.osr.SpatialReference()
        #inSpatialRef.ImportFromEPSG(input_epsg)
        #outSpatialRef = osgeo.osr.SpatialReference()
        #outSpatialRef.ImportFromEPSG(output_epsg)
        inProj = Proj(init='epsg:'+str(input_epsg))
        outProj = Proj(init='epsg:'+str(output_epsg))
	#coordTransform = osgeo.osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

        new_mesh = SchismMesh()
        new_mesh._nodes = np.copy(self.mesh.nodes)
        new_mesh._elems = np.copy(self.mesh.elems)

        #point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        for i, node in enumerate(self.mesh.nodes):
            #point.AddPoint(node[0], node[1])
            #point.Transform(coordTransform)
            new_mesh.nodes[i,0],new_mesh.nodes[i,1]=transform(inProj,outProj,node[0],node[1])
	    #new_mesh.nodes[i,0] = point.GetX()
            #new_mesh.nodes[i,1] = point.GetY()
        gr_writer = Gr3IO()
        gr_writer.write(new_mesh, fname)


    def write_fluxregions(self, lines, out_fname = 'fluxflag.prop'):
        """ Create and write flux_regions.gr3 with the given lines
            lines = list of pairs of (X, Y) of line segments
            out_fname = output file name
        """
        down_paths, up_paths = self._create_fluxlines(lines)
        mesh = self._input.mesh
        flags = np.empty((mesh.n_elems(),1))
        flags.fill(-1)
        for i, (down_path, up_path) in enumerate(zip(down_paths, up_paths)):
            for j in range(len(down_path) - 1):
                edge_i = mesh._find_edge([down_path[j], down_path[j+1]])
                element_1, element_2 = mesh._edges[edge_i][3:5]
                down = False
                for node_i in mesh.elems[element_1]:
                    if not (node_i in down_path or node_i in up_path):
                        down = True
                        break
                if down:
                    flags[element_1] = i
                    flags[element_2] = i + 1
                else:
                    flags[element_2] = i
                    flags[element_1] = i + 1

        index = np.arange(flags.shape[0]).reshape((flags.shape[0],-1))
        index += 1
        elementflags = np.concatenate((index, flags), axis=1)
        np.savetxt(out_fname, elementflags, fmt='%d')

    def _create_fluxlines(self, lines):
        """ This creates a flux line.
            lines = line segments
            return = region arrays of nodes
        """
        mesh = self._input.mesh
        nodes = mesh.nodes
        flags = np.empty((mesh.n_nodes()))
        flags.fill(-1)
        up_paths = []
        down_paths = []
        count = 0
        for name, line in lines.iteritems():
            line = map(float, line.split())
            down_path, up_path = mesh.find_two_neighboring_paths(line)
            for node_i in up_path:
                if not flags[node_i] == -1:
                    raise Exception("Flux lines overlap")
                flags[node_i] = count
            for node_i in down_path:
                if not flags[node_i] == -1:
                    raise Exception("Flux lines overlap")
                flags[node_i] = count + 1
            if len(down_path) < 1 or len(up_path) < 1:
                print()
            up_paths.append(down_path)
            down_paths.append(up_path)
            count += 1

        return down_paths, up_paths

    def _area_of_poligon(x):
        """ Based on shoelace formula
        """
        sum = 0.
        for i in range(len(x)):
            sum += x[i][0] * x[i+1][1] - x[i+1][0] * x[i][1]

        return fabs(sum) * .5

    def _clip(subjectPolygon, clipPolygon):
        """  Sutherland-Hodgman polygon clipping algorithm
             from Rosetta code
        """
        def inside(p):
            return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])

        def computeIntersection():
            dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ]
            dp = [ s[0] - e[0], s[1] - e[1] ]
            n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
            n2 = s[0] * e[1] - s[1] * e[0]
            n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
            return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]

        outputList = subjectPolygon
        cp1 = clipPolygon[-1]

        for clipVertex in clipPolygon:
            cp2 = clipVertex
            inputList = outputList
            outputList = []
            s = inputList[-1]

            for subjectVertex in inputList:
                e = subjectVertex
                if inside(e):
                    if not inside(s):
                        outputList.append(computeIntersection())
                    outputList.append(e)
                elif inside(s):
                    outputList.append(computeIntersection())
                s = e
            cp1 = cp2
        return(outputList)

    def creart_sources_from_user_input(self):
        """ Create a source from a user input
        """
        # Read a user input
        source_inputs = self._read_source_user_input()
        source_th_inputs = self._read_source_th_user_input()
        out_fname = 'source_sink.in'
        self._write_source_file(source_inputs, out_fname)
        out_fname = 'vsource.th'
        self._write_volume_source_history_file(out_fname)

    def _read_source_user_input(self, fname = 'source_sink.user'):
        """ Read in a user input
        """
        f = open(fname, 'r')
        inputs = []
        for l in f:
            l = l.strip()
            if l[0] == '#':
                continue
            tokens = l.split()
            if len(tokens) < 4:
                raise Exception('the file is corrupted.')
            name = tokens[0]
            coord = (float(tokens[1]), float(tokens[2]))
            if tokens[3].upper() == "SOURCE":
                type = 1
            elif tokens[3].upper() == "SINK":
                type = -1
            else:
                raise Exception('the file is corrupted.')
            inputs.append((name, coord, type))
        f.close()
        return inputs

    def _read_source_th_user_input(self, fname = 'source_sink_th.user'):
        """ Read in a user input for time history of source/sink
        """
        f = open(fname, 'r')

        f.close()

    def _write_source_file(self, inputs, fname = 'source_sink.in'):
        """ Write a source/sink file with the give user input data
        """
        f = open(fname, 'w')
        # Count the number of source/sink
        n_sources = sum(item[2] == 1 for item in inputs)
        n_sinks = sum(item[2] == -1 for item in inputs)
        # Sources
        buf = "%d   ! total # of elements with sources\n" % n_sources
        f.write(buf)
        for item in inputs:
            if item[2] == 1:
                element_i = self.mesh.find_elem(item[1])
                buf = "%d  ! %s\n" % (element_i, item[0])
                f.write(buf)
        # Sink
        buf = "%d   ! total # of elements with sinks\n" % n_sinks
        f.write(buf)
        for item in inputs:
            if item[2] == -1:
                element_i = self.mesh.find_elem(item[1])
                buf = "%d  ! %s\n" % (element_i, item[0])
                f.write(buf)

        f.flush()
        f.close()

    def _write_volume_source_history_file(self, fname = 'msource.th'):
        """ Write a time history file of sources/sinks.
        """
        pass

    def create_structures(self, structure_data):
        """ Create structures
        """
        self._input.clear_structures()
        nudging = float(structure_data['nudging'])
        self._input._nudging = nudging
        structures = structure_data['structures']
        for name, structure in structures.iteritems():
            struct = SchismStructure()
            struct.name = name
            struct.coords = map(float, structure['line'].split())
            struct.type = structure['type'].lower()
            struct.properties = structure['configuration']
            #struct.use_timeseries = int(structure['use time series'])
            # Find node pairs
            up_path, down_path = \
                     self.mesh.find_two_neighboring_paths(struct.coords)
            struct.node_pairs = zip(up_path, down_path)
            # Reference pair
            k = 'reference'
            ref = structure[k].lower() if k in structure.keys() else 'self'
            if ref == 'self':
                struct.reference_pair = self._generate_reference_pair(up_path,
                                                                   down_path)
            else:
                struct.reference_pair = [None, None]
            self._input.add_structure(struct)

        for struct in self._input.structures:
            if struct.reference_pair[0] is None:
                ref = structures[struct.name]['reference']
                for s in self._input.structures:
                    if s.name == ref:
                        struct.reference_pair = s.reference_pair
                        break

    def write_structures(self, fname = "hydraulics.in"):
        """ Write a SCHISM structure file
            fname = output file name
        """
        struct_writer = SchismStructureIO(self._input)
        struct_writer.write(fname)

    def _generate_reference_pair(self, up_path, down_path):
        """ Generate a new referece pair from the current node pairs.
            For now, it picks the neighboring nodes around the
            middle node pair.
            node_pairs = the list of node pairs
            return = the new reference pair
        """
        ref_up = self._find_reference_node(up_path, down_path)
        ref_down = self._find_reference_node(down_path, up_path)
        return (ref_up, ref_down)

    def _find_reference_node(self, path1, path2):
        # TODO: The safety of this code needs to check further
        mesh = self.mesh
        # path1
        center_node_i = path1[len(path1)/2]
        neighbors = mesh.get_neighbor_nodes(center_node_i)
        candidates = []
        for n in neighbors:
            if not (n in path1 or n in path2):
                candidates.append(n)
        if len(candidates) < 1:
            raise Exception("No reference node founde")
        elif len(candidates) == 1:
            return candidates[0]
        else:
            # If there are multiple candidates, a most perpendicular one
            # to the line constructed by the two end points of the path1
            # is chosen.
            tangent = mesh.nodes[path1[-1]][:2] - mesh.nodes[path1[0]][:2]
            tangent /= np.linalg.norm(tangent)
            min = 1.
            hit = -1
            for c in candidates:
                vec = mesh.nodes[c][:2] - mesh.nodes[center_node_i][:2]
                vec /= np.linalg.norm(vec)
                dot = fabs(np.dot(vec, tangent))
                if dot < min:
                    min = dot
                    hit = c
            return hit

    def _parse_attribute(self, expr):
        """ Parse expression that can be understood by the tool
        """
        expr = re.sub(r"\bx\b", "node[0]", expr)
        expr = re.sub(r"\by\b", "node[1]", expr)
        expr = re.sub(r"\bz\b", "node[2]", expr)
        return expr

    def _partition_nodes_with_polygons(self, polygons, default):
        """ Partition the grid with the given polygons.
            Each node (not element) will be assigned with an integer ID
            which is the index of the polygons.
            If some polygons overlap, the latter one will trump the former one.
            The area that are not covered by any of the given polygons
            will have a default negative one attribute.
            polygons = a list of Polygon instances
            return = numpy array with the attributes
        """
        
        mesh = self.mesh
        if isinstance(default, str) and default.lower() == 'none':
            default = None
        if default is None:
            # Use depth
            attr = np.copy(mesh.nodes[:, 2])
        else:
            # Fill default values
            attr = np.empty(mesh.n_nodes())
            attr.fill(default)
        for name, polygon in polygons.iteritems():
            vertices = map(float, polygon['vertices'].split())
            if len(vertices) % 2 != 0:
                raise ValueError('The number of coordinates in vertices are wrong.')
            poly_type = polygon['type'].lower() if 'type' in polygon.keys() else "none"
            attribute = polygon['attribute']
            poly = Polygon(vertices, name, attribute, poly_type)
            if isinstance(attribute, str):
                is_eqn = True
                expr = self._parse_attribute(attribute)
            else:
                expr = None
                is_eqn = False
            box = poly.box()
	
            nodes = mesh.find_nodes_in_box(box)
            empty = True
       	     
        for node_i in nodes:
                node = mesh.nodes[node_i]
                flag = poly.check_point_inside_polygon(node[:2])
                if is_eqn:
                    try:
                        value = eval(expr)
                    except:
                        print ("Eqn:", attribute)
                        raise ValueError("The polygon equation does not seem to be well-formed..")
                else:
                    value = attribute
                if flag:
                    empty = False
                    if poly.type == "none":
                        attr[node_i] = value
                    elif poly.type == "min":
                        if attr[node_i] < value:
                            attr[node_i] = value
                    elif poly.type == "max":
                        if attr[node_i] > value:
                            attr[node_i] = value
                    else:
                        raise Exception('Not supported polygon type')


        n_missed = sum(1 for i, a in enumerate(attr) if a == default)
        if n_missed > 0:
            print( "Warning: there are %d nodes that do not belong " \
            "to any polygon..." % n_missed)
            if default is not None:
                print( "Default value of %.f is used for them." % default)
        return attr


    def create_node_partitioning(self, gr3_fname, polygon_data):
        """ Create a gr3 file with node partitioning using
            polygons in the polygon file.
            gr3_fname = output gr3 file name
            polygons = polygons
        """
        option_name = 'default'
        if option_name in polygon_data.keys():
            default = polygon_data[option_name]
        else:
            default = None
        polygons = polygon_data['polygons']
        attr = self._partition_nodes_with_polygons(polygons, default)
        if gr3_fname != 'hgrid.gr3':
            self.write_hgrid(gr3_fname, attr, False)
        else:
            self.write_hgrid(gr3_fname, attr)

    def create_prop_partitioning(self, prop_fname, polygon_data):
        """ Create a prop file with element partitioning using
            polygons in the polygon file.
            prop_fname = output prop file name
            polygons = polygons
        """
        option_name = 'default'
        default = polygon_data[option_name]
        polygons = polygon_data['polygons']
        attr = self._partition_nodes_with_polygons(polygons, default)
        mesh = self.mesh
        elementflags = np.empty((mesh.n_elems(), 1))
        elementflags.fill(0.)
        for i, element in enumerate(mesh.elems):
            flags = attr[element]
            if np.amax(flags) == 1.:
                elementflags[i] = 1.
        index = np.arange(mesh.n_elems()).reshape((mesh.n_elems(), 1))
        index += 1
        np.savetxt(prop_fname, np.concatenate((index, elementflags), axis=1), fmt="%d")


    def modify_depth(self, polygon_data):
        """ Modify depth with the polygon information
        """
        polygons = polygon_data['polygons']
        attr = self._partition_nodes_with_polygons(polygons, None)
        self.mesh.nodes[:,2] = attr


    def create_source_sink_in(self, source_sinks, out_fname = "source_sink.in"):
        """ Create source_sink.in from source/sink location information
            in_fname = input file name
        """
        # TODO: I need to use a common source/sink I/O routines.
        key = 'sources'
        sources = source_sinks[key] if key in source_sinks.keys() else dict()
        key = 'sinks'
        sinks = source_sinks[key] if key in source_sinks.keys() else dict()

        fout = open(out_fname, 'w')
        buf = "%d ! total # of elems with sources\n" % len(sources.keys())
        fout.write(buf)
        for name, coord in sources.iteritems():
            coord = map(float, coord.split())
            element_i = self.mesh.find_elem(coord)
            if element_i is None:
                element_i = self.mesh.find_closest_elems(coord)
                buf = "%d ! %s, nudged\n" % (element_i + 1, name)
                fout.write(buf)
            else:
                buf = "%d ! %s \n" % (element_i + 1, name)
                fout.write(buf)

        buf = "\n%d ! total # of elems with sinks\n" % len(sinks.keys())
        fout.write(buf)
        for name, coord in sinks.iteritems():
            coord = map(float, coord.split())
            element_i = self.mesh.find_elem(coord)
            if element_i is None:
                element_i = self.mesh.find_closest_elems(coord)
                buf = "%d ! %s, nudged\n" % (element_i + 1, name)
                fout.write(buf)
            else:
                buf = "%d ! %s \n" % (element_i + 1, name)
                fout.write(buf)

        fout.flush()
        fout.close()


    def reorder_open_boundaries(self, order):
        """ Reorder open boundaries with the given order
            order = list of open boundary names
            TODO: How to put the name of the boundaries is not settled.
        """
        open_boundaries = list(boundary for boundary in self.mesh.boundaries
                           if boundary.btype == OPEN_BOUNDARY)
        names = []
        for boundary in open_boundaries:
            # extract names from the comments
            p1 = boundary.comment.find("\"")
            p2 = boundary.comment.find("\"", p1 + 1)
            name = boundary.comment[p1 + 1:p2]
            names.append(name)

        new_order = [names.index(i) for i in order]
        self.mesh.boundaries[:len(open_boundaries)] = \
                  [open_boundaries[i] for i in new_order]

    def _interpolate(self, times, data, dt_out, out_fname):
        """ Interpolate tide data.
            The code is copied and modified from interpolate_elev.py.
            times = an array of time stamps
            data = an array of tide data
            dt_out = delta t for output
            out_fname = output file name

            TODO: The code may need to be cleaned up.
        """
        ntimes = times.shape[0]
        stime = times[0]
        etime = times[-1]
        dt_in = times[1] - times[0]

        new_times = np.arange(stime, etime, dt_out)
        nout = new_times.shape[0]

        g = open(out_fname, 'w')

        # There is a limit on the size of the array for 1d
        # spline interpolation (interp1d).
        # Testing on PC indicates it will work up to 2500.
        # The larger the array, the slower the performance.
        # We use 1000 in this script.
        # To minimize edge effect, at 20% on each side.
        max_size_spline = 1000
        add_size = int(1000 * 0.2)

        dt_ratio = int(dt_in / dt_out)
        loop_step = max_size_spline * dt_ratio
        if nout % loop_step != 0:
            iloop = nout/loop_step+1
        else:
            iloop = nout/loop_step



        for j in range(0,iloop):


            if j == 0:
                iin_start = 0
            else:
                iin_start = j * max_size_spline - add_size

            if j == iloop - 1:         #for the last run of the loop
                iin_end = ntimes - 1    #set the end to the last input element
            else:
                iin_end = (j+1) * max_size_spline + add_size
                if iin_end > (ntimes-1):
                    iin_end = ntimes-1

            iout_start = j * loop_step
            iout_end = iout_start + loop_step
            # to avoid "ValueError:
            # A value in x_new is above the interpolation range."
            # reduce output end time
            if iout_end * dt_out > (iin_end-1) * dt_in:
                iout_end = (iin_end-1) * dt_ratio

            spline = interp1d(times[iin_start:iin_end],
                              data[iin_start:iin_end],kind='cubic')
            est_times = new_times[iout_start:iout_end]
            # print est_times[0],est_times[-1],times[iin_start],times[iin_end]
            new_data=spline(est_times)

            for i in range(len(new_data)):
                g.write("%s    %s\n" % (est_times[i] + dt_out, new_data[i]))

        g.close()

    def _adjust_ocean_bc(self):
        """
        """
        # Read ocean map
        fname = 'ocean.gr3'
        gr3io = Gr3IO()
        ocean_mesh = gr3io.read(fname)

        # Read ocean tide harmonics
        fname = 'ocean.ap'
        f = open(fname, 'r')
        f.readline()
        l = f.readline().strip()
        n_nodes = int(l)
        if n_nodes != ocean_mesh.n_nodes():
            raise Exception('Node numbers in the mesh and tide harmonic files'
                            ' are not identical.')
        l = f.readline().strip()
        n_harmonics = int(l)
        harmonics = {}
        for i in range(n_harmonics):
            l = f.readline().strip()
            name = l
            harmonic = np.empty([n_nodes, 2])
            for j in range(n_nodes):
                l = f.readline().strip()
                tokens = l.split()
                harmonic[j, ] = map(float, tokens)
            harmonics[name] = harmonic

        # Calculate weights
        # Get ocean boundary: Assume the first one is the ocean
        ocean_boundary = self.mesh.boundaries[0]
        for node_i in ocean_boundary.nodes:
            element_i = ocean_mesh.find_elem(self.mesh.nodes[node_i][:2])
            if element_i is None:
                element_i = ocean_mesh.find_elem_with_tolerance(\
                    self.mesh.nodes[node_i][:2], 0.1)

    def interpolate_tide(self, time_start, time_end, dt, in_fname, out_fname):
        """ Interpolate an observed tide time series with the given dt.
            time_start = start time in Python datetime
            time_end = end time in Python datetime
            dt = delta t for the interpolation output
            in_fname = input file name of tide records
            out_fname = output file name
        """
        ## Read the tide input
        f = open(in_fname, 'r')
        # ignore first three header lines
        for i in range(3):
            f.readline()
        timestamps = []
        values = []
        hit = 0
        for l in f:
            # Parse it
            tokens = l.split(',')
            time_str = tokens[0] + "," + tokens[1]
            t = datetime.datetime.strptime(time_str, "%m/%d/%Y,%H:%M")
            if t < time_start:
                pass
            elif t <= time_end:
                delta =  t - time_start
                delta_in_sec = delta.days  * 86400 + delta.seconds
                timestamps.append(delta_in_sec)
                values.append(float(tokens[2]))
            else:
                break
        f.close()

        times = np.array(timestamps)
        data = np.array(values)

        # Interpolate
        self._interpolate(times, data, dt, out_fname)

    def create_open_boundaries(self, segments):
        """ Create open boundaries with the given two end points file
            fname = the end points file
        """
        self.mesh.clear_boundary()
        boundary_only = True
        for name, segment in segments.iteritems():
            p = map(float, segment.split())
            if len(p) < 4:
                raise Exception('Open boundary file seems corrupted.')
            start = self.mesh.find_closest_nodes(p[0:2], 1, boundary_only)
            end = self.mesh.find_closest_nodes(p[2:4], 1, boundary_only)
            path = self.mesh.shortest_path(start, end, boundary_only)
            comment = "! Open boundary \"%s\"" % name
            self.mesh.add_boundary(path, OPEN_BOUNDARY, comment)

    def _read_flow_xsects(self, fname):
        """ Read flow transects information from a file
            fname = file name to read in
            return = array of line segmentns
        """
        f = open(fname, 'r')
        lines = []
        for l in f:
            l = l.strip()
            if len(l) == 0:
                continue
            if l[0] == '#':             # Comment line
                continue
            tokens = l.split(',')
            if len(tokens) < 8:
                raise Exception()
            vals = map(float, tokens[4:8])
            lines.append((vals[0:2],vals[2:4]))
        f.close()
        return lines

    def create_flux_regions(self, fluxlines, out_fname):
        """ Create fluxflag.gr3
            in_fname = file name for flow transect information
            out_fname = output file name
        """
        # lines = self._read_flow_xsects(in_fname)
        self.write_fluxregions(fluxlines, out_fname)

    def trim_to_left_of_mesh(self, line_segments):
        """ Trim mesh using line_segments.
            This function trims the mesh on the left sides
            of the line segments. The left side here means left when you look
            at the second end point of a line segment from the first one.
            An actual path to trimming is a nodal path that is on the
            right side of the line segment.
            To manage torus like mesh topology, the argument takes an array
            of line segments. The user need to provides a set of line segments
            to make sure the left side of the line segments does not
            cover the whole mesh.
            To trim multiple parts of a mesh, use this function
            multiple times.
            This function clears up any boundary definitions associated with
            the original grid. It is user's responsibility to create them
            again.

            line_segments = array of line segments defined by two end points
        """

        paths = []
        for i, l in enumerate(line_segments):
            p = self.mesh.find_two_neighboring_paths(l)
            paths.append(p[0])
        self.mesh.trim_to_left(paths)
        self.mesh.clear_boundary()



def load_schism(dir):
    """ A load function
    """
    s = SchismSetup()
    s.load(dir)
    return s


def load_gr3(fname):
    """ Read a mesh in GR3 format, and return a SCHISM input instance.

        Parameters
        ----------
        fname: str

        Returns
        -------
        SCHISM input instance
    """
    s = SchismSetup()
    # Read the grid
    gr3io = Gr3IO()
    gr3io.verbose = 2
    if not os.path.exists(fname):
        raise Exception("A grid file is not found!")
    s.mesh = gr3io.read(fname)
    return s