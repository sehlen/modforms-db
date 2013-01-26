from schema import *
from sage.all import loads, dumps, sage_eval

class ModularSymbols_ambient(ModularSymbols_ambient_DB):
    def sage_object():
        return NotImplementedError()

class ModularSymbols_oldspace_factor(ModularSymbols_oldspace_factor_DB):
    def sage_object():
        return NotImplementedError()
    
class ModularSymbols_newspace_factor(ModularSymbols_newspace_factor_DB):
    def sage_object():
        return NotImplementedError()

class Coefficient(Coefficient_DB):
    def sage_object():
        return NotImplementedError()

class NumberField(NumberField_DB):
    def sage_object():
        return NotImplementedError()

class ModularSymbols_base_field(ModularSymbols_base_field_DB):
    def sage_object():
        return NotImplementedError()

class CoefficientField(CoefficientField_DB):
    def sage_object():
        return NotImplementedError()

class AlgebraicNumber(AlgebraicNumber_DB):
    @property
    def vector(self):
        if self._value is not None:
            v = loads(str(self._value))
            v = sage_eval(v)
            return v

    def _set_value_from_sage(self,x):
        v = x.list()
        v = dumps(str(v))
        self._value = v
