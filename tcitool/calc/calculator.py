import sys
import numpy as np

class Calculator(object):
    """A calculator calculates a thermal comfort index.

    A calculator calculates a thermal comfort index, and may use Generators, and
    Functions to do that calculation. Calculators need to be registerd in the
    main-class (Tool), and need to implement this function.

    Attributes:
        tool: reference to the tool used
        data: the datastore of this calculator
        name: name of the class
        export_params: a dict[(str,str)] of parameters that will be exported to
            the DataStore of tool (tool.data). Keys indicate local names, values
            names after export.
    """
    def __init__(self,tool):
        self.tool = tool
        self.data = tool.data.copy_empty()
        self.name = self.__class__.__name__
        self.export_params = {}

    def require_data(self,*args):
        self.tool.require_data(*args,operation_name=self.name)

    def export(self):
        params = list(self.export_params.keys())
        self.tool.data.merge(
            self.data[params].rename_vars(**self.export_params)
        )

    def run(self):
        self.preface()
        self.main()
        self.postface()

    def preface(self):
        pass
    def main(self):
        raise NotImplementedError()
    def postface(self):
        self.export()

class OptimizationCalculator(Calculator):
    def main(self):
        self.optimize()
