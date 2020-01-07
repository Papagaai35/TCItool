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

    def copy_empty(self):
        return self.ds.coords.to_dataset()

    def has_keys(self,*args):
        return all(map(lambda param: param in self.ds,args))
