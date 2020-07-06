"""
AMRVAC-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.geometry.selection_routines import \
    GridSelector
from .datfile_utils import get_single_block_field_data


class AMRVACIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = 'amrvac'

    def __init__(self, ds):
        BaseIOHandler.__init__(self, ds)
        self.ds = ds
        self.datfile = ds.parameter_filename
        header = self.ds.parameters
        self.block_shape = np.append(header["block_nx"], header["nw"])

    def _read_data(self, grid, field):
        """Retrieve field data from a grid.

        Parameters
        ----------
        grid : yt.frontends.amrvac.data_structures.AMRVACGrid
            The grid from which data is to be read.
        field : str
            A field name.

        Returns
        -------
        data : np.ndarray
            A 3D array of float64 type representing grid data.

        """
        ileaf = grid.id
        offset = grid._index.block_offsets[ileaf]
        field_idx = self.ds.parameters['w_names'].index(field)
        with open(self.datfile, "rb") as istream:
            data = get_single_block_field_data(istream, offset, self.block_shape, field_idx)

        # Always convert data to 3D, as grid.ActiveDimensions is always 3D
        while len(data.shape) < 3:
            data = data[..., np.newaxis]
        return data

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError