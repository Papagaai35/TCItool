class MissingDataError(Exception):
    """Error generated when certain parameters are nessesary for a calculation,
    but one/more parameters are not available"""
    pass

class UnknownCalculatorWarning(UserWarning):
    """Warning generated when a calculation is requested, but the calculator is
    not found"""
    pass
