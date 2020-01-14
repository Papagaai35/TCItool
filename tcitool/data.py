import xarray as xr

class DataStore(object):
    def __init__(self,file_or_xarray=None,**kwargs):
        """Inits this DataStore

        Args:
            filename_or_xarray: A file path describing the location of the file
            to be loaded (if string), or a xarray.Dataset containing the data.
            **kwargs: will be passed to xarray.open_dataset
        """
        self._ds = None
        if file_or_xarray is not None:
            self.load(file_or_xarray,**kwargs)

    def ds():
        doc = "The xarray.Dataset containing the data"
        def fget(self):
            if self._ds is None:
                raise ValueError("ds has not been set yet.")
            return self._ds
        def fset(self, value):
            if not isinstance(value, xr.Dataset):
                raise TypeError("You may only set ds with a xarray.Dataset")
            self._ds = value
        return locals()
    ds = property(**ds())
    def shape():
        doc = "The shape of the data in the xarray.Dataset"
        def fget(self):
            params = list(self.ds.keys())
            shp = self.ds[params[0]].shape
            for param in params:
                if self.ds[param].shape != shp:
                    raise ValueError("Different parameters have different "
                        "shapes")
            return shp
        return locals()
    shape = property(**shape())
    def dims():
        doc = "The dimensions of the data in the xarray.Dataset"
        def fget(self):
            params = list(self.ds.keys())
            dims = self.ds[params[0]].dims
            for param in params:
                if self.ds[param].dims != dims:
                    raise ValueError("Different parameters have different "
                        "dimensions")
            return dims
        return locals()
    dims = property(**dims())
    def chunks():
        doc = "The chunk size of the data in the xarray.Dataset"
        def fget(self):
            params = list(self.ds.keys())
            chunks = sum(self.ds[params[0]].chunks, ())
            for param in params:
                if sum(self.ds[param].chunks, ()) != chunks:
                    raise ValueError("Different parameters have different "
                        "chunk sizes")
            return chunks
        return locals()
    chunks = property(**chunks())

    def __getitem__(self,*args,**kwargs):
        """Alias for xarray.Dataset.__getitem__"""
        return self.ds.__getitem__(*args,**kwargs)

    def __setitem__(self,*args,**kwargs):
        """Alias for xarray.Dataset.__setitem__"""
        return self.ds.__setitem__(*args,**kwargs)

    def __delitem__(self,*args,**kwargs):
        """Alias for xarray.Dataset.__delitem__"""
        return self.ds.__delitem__(*args,**kwargs)

    def __contains__(self,*args,**kwargs):
        """Alias for xarray.Dataset.__contains__"""
        return self.ds.__contains__(*args,**kwargs)

    def get(self,key,default):
        return self.ds[key] if key in self else default

    def load(self,file_or_xarray,**kwargs):
        """Loads data form a xarray or file

        Args:
            filename_or_xarray: A file path describing the location of the file
            to be loaded (if string), or a xarray.Dataset containing the data.
            **kwargs: will be passed to xarray.open_dataset
        """
        self.ds = ( file_or_xarray
                    if isinstance(file_or_xarray, xr.Dataset)
                    else xr.open_dataset(file_or_xarray,**kwargs))
        self.transpose_default()

    def save(self,filepath,**kwargs):
        kwargs.update({'path':filepath})
        self.ds.to_netcdf(**kwargs)

    def buffer(self,filepath,save_kwargs=None,load_kwargs=None):
        if save_kwargs is None:
            save_kwargs = {}
        if not isinstance(save_kwargs,dict):
            raise TypeError("'save_kwargs' must be a dictionary")
        if load_kwargs is None:
            load_kwargs = {}
        if not isinstance(load_kwargs,dict):
            raise TypeError("'load_kwargs' must be a dictionary")

        self.save(filepath,**save_kwargs)
        self.load(filepath,**save_kwargs)

    def merge(self,datastore_or_xarray):
        """Merge with an other DataStore or xarray.Dataset"""
        if isinstance(datastore_or_xarray,DataStore):
            self.ds = xr.merge([self.ds,datastore_or_xarray.ds])
        elif isinstance(datastore_or_xarray,xr.Dataset):
            self.ds = xr.merge([self.ds,datastore_or_xarray])
        else:
            raise TypeError("datastore_or_xarray must be a DataStore or "
                            "xarray.Dataset")
        self.transpose_default()

    def copy_empty(self):
        return self.ds.coords.to_dataset()

    def has_keys(self,*args):
        return all(map(lambda param: param in self.ds,args))

    def get_coord_var(self,coord):
        coordxarray_1d = xr.DataArray(self.ds.coords[coord].values,
            dims=[coord],attrs=self.ds.coords[coord].attrs)
        _, coordxarray_md = xr.broadcast(self.ds,coordxarray_1d)
        if self.get_chunk_size():
            coordxarray_md = coordxarray_md.chunk(self.get_chunk_size())
        return coordxarray_md

    def get_chunk_size(self):
        return {coord: chunks[0] for coord, chunks in self.ds.chunks.items()}

    def refresh_chunks(self):
        if self.get_chunk_size():
            self.ds = self.ds.chunk(self.get_chunk_size())


    def persist(self):
        self._ds.persist()

    def transpose_default(self,preferd_order=None):
        preferd_order = (['time','latitude','longitude']
                         if preferd_order is None
                         else preferd_order)
        if (len(self.ds.coords)==len(preferd_order) and
                all(coord in self.ds.coords for coord in preferd_order)):
            self.ds = self.ds.transpose(*preferd_order)

    def table_repr(self):
        table_shape = [2,4,5,9,5]
        header_line = ['Nr','Name','Units', 'Long name', 'Shape']
        table = []
        keys = list(self.ds.keys())
        for keyi in range(len(keys)):
            key = keys[keyi]
            tableline = (str(keyi),key,
                self.ds[key].attrs.get('units','?'),
                self.ds[key].attrs.get('long_name','?'),
                str(self.ds[key].shape))
            for i in range(4):
                table_shape[i] = max(table_shape[i],len(tableline[i]))
            table.append(tableline)
        fmt_string = "{:%d} | {:<%d} | {:<%d} | {:<%d} | {:<%d}"%tuple(table_shape)
        print(fmt_string.format(*header_line))
        print('-'*(sum(table_shape)+12))
        for l in table:
            print(fmt_string.format(*l))
