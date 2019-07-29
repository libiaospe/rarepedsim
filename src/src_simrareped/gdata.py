# $File: gdata.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)

from tarfile import TarFile, TarInfo, DEFAULT_FORMAT, ENCODING, \
     ReadError, CompressionError
import os
# import spower.hickle as serialization
import cPickle as serialization
import numpy as np
from scipy.sparse import csr_matrix

class GData(dict):
    def __init__(self, data, name = None):
        if type(data) is str:
            data = serialization.load(open(data))
        elif type(data) is dict:
            pass
        else:
            # is file stream
            data = serialization.load(data)
        self.update(data)
        # This is the key name for "marker"
        self.__marker = name

    def sink(self, filename):
        # serialization.dump(dict(self), filename, 'w')
        serialization.dump(dict(self), open(filename, 'wb'), protocol = -1)
        
    def compress(self):
        def convert(item):
            # implementation of 2D list to sparse matrix
            try:
                return csr_matrix(np.array(item, dtype=np.int8))
            except:
                return item
        self[self.__marker] = convert(self[self.__marker])

    def decompress(self, tolist = True):
        def convert(item):
            # implementation of sparse matrix to numpy 2D array or list
            try:
                # return item.todense().view(np.ndarray)
                if tolist:
                    return item.todense().tolist()
                else:
                    return item.todense()
            except:
                return item
        self[self.__marker] = convert(self[self.__marker])


class GFile(TarFile):
    def __init__(self, name=None, mode='r', fileobj=None,
                 format=DEFAULT_FORMAT, tarinfo=TarInfo,
                 dereference=False, ignore_zeros=False,
                 encoding=ENCODING, errors=None,
                 pax_headers=None, debug=0, errorlevel=0):
        TarFile.__init__(self, name, mode, fileobj,
                 format, tarinfo,
                 dereference, ignore_zeros,
                 encoding, errors,
                 pax_headers, debug, errorlevel)

    @classmethod
    def open(cls, name=None, fileobj=None, bufsize=10240, **kwargs):
        '''Disable the write interface'''
        for comptype in cls.OPEN_METH:
            func = getattr(cls, cls.OPEN_METH[comptype])
            if fileobj is not None:
                saved_pos = fileobj.tell()
            try:
                return func(name, "r", fileobj, **kwargs)
            except (ReadError, CompressionError), e:
                if fileobj is not None:
                    fileobj.seek(saved_pos)
                continue
            raise ReadError("file could not be opened successfully.")

    def getdata(self, name):
        return GData(self.extractfile(name), name)

if __name__ == '__main__':
    import sys
    with GFile.open(sys.argv[1]) as g:
        names = g.getnames()
        sys.stderr.write('\n'.join(names))
        data = g.getdata(names[0])
        print data.keys()
        print data
        data.decompress()
        print data
        data.sink('myfile')
