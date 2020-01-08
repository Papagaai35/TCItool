import tcitool

class GeneratorRegistry(object):
    def __init__(self,tool):
        self.tool = tool
        self.generators = {}
        tcitool.CommonMeteoGenerators.register_generators(self)
        tcitool.SolarGenerators.register_generators(self)

    def register(self,func,provides=None,requires=None,options=None):
        """Registers a callable, that may be used as a generator in this
        GeneratorRegistry

        Args:
            func: callable
            provides: a (list of) string(s), describing the data-keys generated
                by this generator-method
            requires: a (list of) string(s), describing the data-keys needed to
                be able to generate the provides, by this generator-method
            options: a (list of) string(s) describing the options that need to
                be set in the tool, to be able to calculate the provides
        """
        kwargs = {'provides':provides,'requires':requires,'options':options}
        for kwargn in kwargs.keys():
            kwargv = kwargs[kwargn]
            if kwargv is None or (isinstance(kwargv, list) and len(kwargv)==0):
                kwargv = []
            if isinstance(kwargv, str):
                kwargv = [kwargv]
            if not all(map(lambda elem: isinstance(elem,str),provides)):
                raise TypeError("'%s' should be a list of stings"%kwargn)
            kwargs[kwargn] = kwargv

        if not callable(func):
            raise TypeError("'func' should be callable")
        kwargs['func'] = func

        for provides_param in kwargs['provides']:
            if provides_param not in self.generators:
                self.generators[provides_param] = []
            self.generators[provides_param].append(kwargs)

    def find_and_run(self,param):
        """Finds a relevant generator, to satisfy the data need, and runs it.

        Finds a generator, that generates the parameter, where all requirements
        are satisfied (requires and options). If multiple generators exsist, the
        frist generator (in the order of registering) will be used.

        Args:
            param: string describing the data parameter needed.

        Returns:
            True if succeeded
            False if the nessesary `requires` or `options` were not pressent
            None if no generator-method could be found
        """
        if param in self.generators:
            for gen in self.generators[param]:
                if (self.tool.data.has_keys(*gen['requires']) and
                    self.tool.has_options(*gen['options'])):
                    gen['func'](self.tool)
                    return True
            return False
        return None
