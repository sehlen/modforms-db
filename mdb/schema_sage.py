from schema import *
from sage.all import sage_eval
from conversions import *

class SageObject_DB(object):
    def sage_object(self):
        r"""
            Returns a Sage object corresponding to self.
        """
        return NotImplementedError()

    def from_sage_object(self, so):
        r"""
            Sets the properties of self from a corresponding sage object.
        """
        return NotImplementedError()

class ModularSymbols_ambient(ModularSymbols_ambient_DB, SageObject_DB):
    def sage_object(self):
        d=self.as_dict()
        return dict_to_ambient_sage(d)

    def as_dict(self):
        d=dict()
        d['basis'] = self.get_basis()
        d['manin'] = self.get_manin()
        d['rels'] = self.get_rels()
        d['mod2term'] = self.get_mod2term()
        d['space'] = (self.level, self.weight, self.character)
        return d

    def from_sage(self, M):
        d=sage_ambient_to_dict(M)
        self.level = int(M.level())
        self.weight = int(M.weight())
        self.character = int(dirichlet_character_conrey_galois_orbit_rep(M.character()))
        self.set_basis(d['basis'])
        self.set_manin(d['manin'])
        self.get_rels(d['rels'])
        self.get_mod2term(d['mod2term'])
        self.dimension_modular_forms = int(M.dimension())
        self.dimension_new_cusp_forms = int(M.cuspidal_subspace().new_subspace().dimension())
        self.dimension_cusp_forms = int(M.cuspidal_subspace().dimension())

class ModularSymbols_oldspace_factor(ModularSymbols_oldspace_factor_DB, SageObject_DB):
    def sage_object():
        return NotImplementedError()
    
class ModularSymbols_newspace_factor(ModularSymbols_newspace_factor_DB, SageObject_DB):
    def sage_object(self):
        d['B'] = self.get_B()
        d['Bd'] = self.get_Bd()
        d['v'] = self.get_v()
        d['nz'] = self.get_nz()
        d['ambient'] = self.ambient.as_dict()
        return dict_to_factor_sage(d)
    
    def from_sage(self, M, names='a'):
        d=sage_factor_to_dict(M)
        self.set_B(d['B'])
        self.set_Bd(d['Bd'])
        self.get_v(d['v'])
        self.set_nz(d['nz'])
        self.dimension = int(M.dimension())
        #self.has_cm = has_cm(M)
        # now we set the coefficient field
        extension_field = M.eigenvalue(1,name=names).parent()
        if extension_field != M.base_ring(): # .degree() != 1 and rings.is_NumberField(extension_field):
            assert extension_field.base_field() == M.base_ring()
            minpoly = extension_field.relative_polynomial()
            degree = M.base_ring().degree()*minpoly.degree()
        else:
            minpoly = extension_field.defining_polynomial()
            degree = extension_field.degree()
        self.coefficient_field = NumberField()
        self.coefficient_field.minimal_polynomial = str(minpoly)
        self.coefficient_field.degree = degree

        

class Coefficient(Coefficient_DB, SageObject_DB):
    def sage_object():
        return NotImplementedError()

class NumberField(NumberField_DB, SageObject_DB):
    def sage_object(self):
        minpoly = self.minimal_polynomial
        

class ModularSymbols_base_field(ModularSymbols_base_field_DB, SageObject_DB):
    def sage_object():
        return NotImplementedError()

class CoefficientField(CoefficientField_DB, SageObject_DB):
    def sage_object():
        return NotImplementedError()

class AlgebraicNumber(AlgebraicNumber_DB, SageObject_DB):
    
    def as_vector(self):
        return sage_eval(self.get_value())
    
    def _set_value_from_sage(self,x):
        v = x.list()
        self.set_value(v)
