from itertools import groupby

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


# http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
def grid_sequences(grids):
    g_iter = sorted(grids, key=lambda g: g.id)
    for _, g in groupby(enumerate(g_iter), lambda i_x1: i_x1[0] - i_x1[1].id):
        seq = [v[1] for v in g]
        yield seq


ii = [0, 1, 0, 1, 0, 1, 0, 1]
jj = [0, 0, 1, 1, 0, 0, 1, 1]
kk = [0, 0, 0, 0, 1, 1, 1, 1]


class IOHandlerParthenon(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "parthenon"

    def __init__(self, ds):
        super().__init__(ds)
        self._handle = ds._handle

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        f = self._handle
        rv = {}
        for field in fields:
            # Always use *native* 64-bit float.
            rv[field] = np.empty(size, dtype="=f8")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )
        last_dname = None
        for field in fields:
            ftype, fname = field
            dname, fdi = self.ds._field_map[fname]
            if dname != last_dname:
                ds = f[f"/{dname}"]
            ind = 0
            for chunk in chunks:
                for gs in grid_sequences(chunk.objs):
                    start = gs[0].id - gs[0]._id_offset
                    end = gs[-1].id - gs[-1]._id_offset + 1
                    data = ds[start:end, fdi, :, :, :].transpose()
                    for i, g in enumerate(gs):
                        ind += g.select(selector, data[..., i], rv[field], ind)
            last_dname = dname
        return rv

    def _read_chunk_data(self, chunk, fields):
        f = self._handle
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        if len(fields) == 0:
            return rv
        for field in fields:
            ftype, fname = field
            dname, fdi = self.ds._field_map[fname]
            ds = f[f"/{dname}"]
            for gs in grid_sequences(chunk.objs):
                start = gs[0].id - gs[0]._id_offset
                end = gs[-1].id - gs[-1]._id_offset + 1
                data = ds[start:end, fdi, :, :, :].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray(data[..., i], "=f8")
        return rv
