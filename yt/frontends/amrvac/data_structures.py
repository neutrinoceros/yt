"""
AMRVAC data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import stat
import weakref
import struct

import numpy as np

from yt.data_objects.grid_patch import \
   AMRGridPatch
from yt.geometry.grid_geometry_handler import \
   GridIndex
from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.static_output import \
   Dataset
from yt.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    boltzmann_constant_cgs, \
    mu_0

from .fields import AMRVACFieldInfo
from .datfile_utils import get_header, get_tree_info



class AMRVACGrid(AMRGridPatch):
    """A class to populate AMRVACHierarchy.grids, setting parent/children relations """
    _id_offset = 0

    def __init__(self, id, index, level):
        #<level> should use yt's convention (start from 0)
        super(AMRVACGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "AMRVACGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        start_index = (self.LeftEdge - self.ds.domain_left_edge)/self.dds
        self.start_index = np.rint(start_index).astype('int64').ravel()
        return self.start_index


class AMRVACHierarchy(GridIndex):
    grid = AMRVACGrid
    def __init__(self, ds, dataset_type="amrvac"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # the index file *is* the datfile
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64

        super(AMRVACHierarchy, self).__init__(ds, dataset_type)

    def _detect_output_fields(self):
        """Parse field names from datfile header, which is stored in self.dataset.parameters"""
        # required method
        self.field_list = [(self.dataset_type, f) for f in self.dataset.parameters["w_names"]]

    def _count_grids(self):
        """Set self.num_grids from datfile header."""
        # required method
        self.num_grids = self.dataset.parameters['nleafs']

    def _parse_index(self):
        """Populate self.grid_* attributes from tree info from datfile header."""
        # required method
        with open(self.index_filename, "rb") as istream:
            vaclevels, morton_indices, block_offsets = get_tree_info(istream)
            assert len(vaclevels) == len(morton_indices) == len(block_offsets) == self.num_grids

        self.block_offsets = block_offsets
        # YT uses 0-based grid indexing, lowest level = 0 (AMRVAC uses 1 for lowest level)
        ytlevels = np.array(vaclevels, dtype="int32") - 1
        self.grid_levels.flat[:] = ytlevels
        self.max_level = np.max(ytlevels)
        assert self.max_level == self.dataset.parameters["levmax"] - 1

        # some aliases for left/right edges computation in the coming loop
        domain_width = self.dataset.parameters["xmax"] - self.dataset.parameters["xmin"]
        block_nx = self.dataset.parameters["block_nx"]
        xmin = self.dataset.parameters["xmin"]
        dx0 = domain_width / self.dataset.parameters["domain_nx"] # dx at coarsest grid level (YT level 0)
        dim = self.dataset.dimensionality

        self.grids = np.empty(self.num_grids, dtype='object')
        for igrid, (ytlevel, morton_index) in enumerate(zip(ytlevels, morton_indices)):
            dx = dx0 / self.dataset.refine_by**ytlevel
            left_edge = xmin + (morton_index-1) * block_nx * dx

            # edges and dimensions are filled in a dimensionality-agnostic way
            self.grid_left_edge[igrid, :dim] = left_edge
            self.grid_right_edge[igrid, :dim] = left_edge + block_nx * dx
            self.grid_dimensions[igrid, :dim] = block_nx
            self.grids[igrid] = self.grid(igrid, self, ytlevels[igrid])

    def _populate_grid_objects(self):
        # required method
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()



class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(self, filename, dataset_type='amrvac',
                units_override=None, unit_system="cgs",
                geometry_override=None):
        # note: geometry_override is specific to this frontend
        self._geometry_override = geometry_override
        super(AMRVACDataset, self).__init__(filename, dataset_type,
                                            units_override=units_override, unit_system=unit_system)
        self.fluid_types += ('amrvac',)
        # refinement factor between a grid and its subgrid
        self.refine_by = 2

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """At load time, check whether data is recognized as AMRVAC formatted."""
        # required class method
        validation = False
        if args[0].endswith(".dat"):
            try:
                with open(args[0], mode="rb") as istream:
                    fmt = "=i"
                    [datfile_version] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
                    if 3 <= datfile_version < 6:
                        fmt = "=ii"
                        offset_tree, offset_blocks = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
                        istream.seek(0,2)
                        file_size = istream.tell()
                        validation = offset_tree < file_size and offset_blocks < file_size
            except:
                pass
        return validation

    def parse_geometry(self, geometry_string):
        """Transform a string such as "Polar_2D" or "Cartesian_1.75D" to yt's standard equivalent (i.e. respectively "polar", "cartesian")."""
        # frontend specific method
        geom = geometry_string.split("_")[0].lower()
        if geom not in ("cartesian", "polar", "cylindrical", "spherical"):
            raise ValueError
        return geom

    def _parse_parameter_file(self):
        # required method
        self.unique_identifier = int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # populate self.parameters with header data
        with open(self.parameter_filename, 'rb') as istream:
            self.parameters.update(get_header(istream))

        self.current_time = self.parameters['time']
        self.dimensionality = self.parameters['ndim']

        # force 3D for this definition
        dd = np.ones(3, dtype="int64")
        dd[:self.dimensionality] = self.parameters['domain_nx']
        self.domain_dimensions = dd

        # the following parameters may not be present in the datfile,
        # dependending on format version
        if self.parameters["datfile_version"] < 5:
            mylog.warning("This data format does not contain geometry or periodicity info")
        if self.parameters.get("staggered", False):
            mylog.warning("'staggered' flag was found, but is currently ignored (unsupported)")

        # parse geometry
        # by order of decreasing priority, we use
        # - geometry_override
        # - "geometry" parameter from datfile
        # - if all fails, default to "cartesian"
        geom_candidates = {"param": None, "override": None}
        amrvac_geom = self.parameters.get("geometry", None)
        if amrvac_geom is None:
            mylog.warning("Could not find a 'geometry' parameter in source file.")
        else:
            geom_candidates.update({"param": self.parse_geometry(amrvac_geom)})

        if self._geometry_override is not None:
            try:
                geom_candidates.update({"override": self.parse_geometry(self._geometry_override)})
            except ValueError:
                mylog.error("Unknown value for geometry_override (will be ignored).")
    
        if geom_candidates["override"] is not None:
            mylog.warning("Using override geometry, this may lead to surprising results for inappropriate values.")
            self.geometry = geom_candidates["override"]
        elif geom_candidates["param"] is not None:
            mylog.info("Using parameter geometry")
            self.geometry = geom_candidates["param"]
        else:
            mylog.warning("No geometry parameter supplied or found, defaulting to cartesian.")
            self.geometry = "cartesian"

        # parse peridiocity
        per = self.parameters.get("periodic", np.array([False, False, False]))
        missing_dim = 3 - len(per)
        self.periodicity = np.append(per, [False]*missing_dim)

        self.gamma = self.parameters.get("gamma", 5.0/3.0)

        # parse domain edges
        dle = np.zeros(3)
        dre = np.ones(3)
        dle[:self.dimensionality] = self.parameters['xmin']
        dre[:self.dimensionality] = self.parameters['xmax']
        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # defaulting to non-cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_matter = 0.0
        self.omega_lambda = 0.0
        self.hubble_constant = 0.0

    # units stuff ===============================================================================
    def _set_code_unit_attributes(self):
        # required method
        # note: this method is never defined in the parent abstract class Dataset
        # but it is called in Dataset.set_code_units(), which is part of Dataset.__init__()
        # so it must be defined here.

        # note2: this gets called later than Dataset._override_code_units()
        # This is the reason why it uses setdefaultattr: it will only fill the gaps lefts by the "override",
        # instead of overriding them again

        # First check if overrides have been supplied, if that's the case use those instead.
        # Assume cgs values and let YT handle conversion if supplied in an 'mks' unit system.
        if self.units_override:
            unit_length = self.units_override.get('length_unit', (1, 'cm'))
            unit_numberdensity = self.units_override.get('numberdensity_unit', (1, 'cm**-3'))
            unit_velocity = self.units_override.get('velocity_unit', (0, 'cm*s**-1'))
            unit_temperature = self.units_override.get('temperature_unit', (1, 'K'))
            mylog.info('Overriding numberdensity_unit: {:1.0e} {}.'.format(unit_numberdensity[0], unit_numberdensity[1]))
        else:
            # Next, check for unit normalisations in the .dat file. Do not use .get() with default value added
            # to explicitly print a warning to the console.
            try:
                unit_length = (self.parameters['unit_length'], 'cm')
                unit_numberdensity = (self.parameters['unit_numberdensity'], 'cm**-3')
                unit_temperature = (self.parameters['unit_temperature'], 'K')
                unit_velocity = (self.parameters['unit_velocity'], 'cm*s**-1')
            # If nothing present in datfile, use default values and print warning
            except KeyError:
                mylog.warning('No unit normalisations found in datfile and no overrides specified. Using defaults.')
                unit_length = (1, 'cm')
                unit_numberdensity = (1, 'cm**-3')
                unit_temperature = (1, 'K')
                unit_velocity = (0, 'cm*s**-1')

        # Convert to YT quantities
        unit_length = self.quan(unit_length[0], unit_length[1])
        unit_numberdensity = self.quan(unit_numberdensity[0], unit_numberdensity[1])
        unit_temperature = self.quan(unit_temperature[0], unit_temperature[1])
        unit_velocity = self.quan(unit_velocity[0], unit_velocity[1])

        # constants needed for conversion to AMRVAC normalisations
        He_abundance = 0.1  # hardcoded in AMRVAC
        mp = mass_hydrogen_cgs
        kb = boltzmann_constant_cgs

        unit_density = (1.0 + 4.0*He_abundance) * mp * unit_numberdensity
        if unit_velocity == 0:
            unit_pressure = ((2.0 + 3.0*He_abundance) * unit_numberdensity * kb * unit_temperature).to('dyn*cm**-2')
            unit_velocity = (np.sqrt(unit_pressure / unit_density)).to('cm*s**-1')
        else:
            unit_pressure = (unit_density * unit_velocity**2).to('dyn*cm**-2')
            unit_temperature = (unit_pressure / ((2.0 + 3.0*He_abundance) * unit_numberdensity * kb)).to('K')
        unit_time = unit_length / unit_velocity
        unit_mass = unit_density * unit_length**3
        unit_magneticfield = (np.sqrt(4*np.pi * unit_pressure)).to('gauss')

        setdefaultattr(self, "length_unit", unit_length)
        setdefaultattr(self, "mass_unit", unit_mass)
        setdefaultattr(self, "time_unit", unit_time)

        setdefaultattr(self, "velocity_unit", unit_velocity)
        setdefaultattr(self, "density_unit", unit_density)
        setdefaultattr(self, "numberdensity_unit", unit_numberdensity)

        setdefaultattr(self, "temperature_unit", unit_temperature)
        setdefaultattr(self, "pressure_unit", unit_pressure)
        setdefaultattr(self, "magnetic_unit", unit_magneticfield)
