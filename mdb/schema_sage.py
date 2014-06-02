from schema import ModularSymbols_ambient_DB_class,Coefficient_DB_class,NumberField_DB_class,ModularSymbols_base_field_DB_class,CoefficientField_DB_class,AlgebraicNumber_DB_class,ModularSymbols_oldspace_factor_DB_class,ModularSymbols_newspace_factor_DB_class
from sage.all import sage_eval
from conversions import *

#Q = NumberField_DB_class()
#Q.minimal_polynomial = 'x'
#Q.degree = int(1)

class SageObject_DB_class(object):
    def sage_object(self):
        r"""
            Returns a Sage object corresponding to self.
        """
        return NotImplementedError()

    def from_sage(self, so):
        r"""
            Sets the properties of self from a corresponding sage object.
        """
        return NotImplementedError()

class ModularSymbols_ambient(ModularSymbols_ambient_class, SageObject_DB_class):
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
        self.character = int(dirichlet_character_to_int(M.character(), convention='Conrey'))
        self.set_basis(d['basis'])
        self.set_manin(d['manin'])
        self.set_rels(d['rels'])
        self.set_mod2term(d['mod2term'])
        self.dimension_modular_forms = int(M.dimension())
        self.dimension_new_cusp_forms = int(M.cuspidal_subspace().new_subspace().dimension())
        self.dimension_cusp_forms = int(M.cuspidal_subspace().dimension())
        self.base_field = ModularSymbols_base_field()
        self.base_field.minimal_polynomial = str(M.base_ring().defining_polynomial())
        self.base_field.degree = int(M.base_ring().degree())
        for N in M.cuspidal_subspace().new_subspace().decomposition():
            NN = ModularSymbols_newspace_factor(dimension = N.dimension())
            self.newspace_factors.append(NN)
            NN.from_sage(N)

class ModularSymbols_oldspace_factor(ModularSymbols_oldspace_factor_class, SageObject_DB_class):
    def sage_object():
        return NotImplementedError()
    
class ModularSymbols_newspace_factor(ModularSymbols_newspace_factor_class, SageObject_DB_class):
    def sage_object(self):
        d['B'] = self.get_B()
        d['Bd'] = self.get_Bd()
        d['v'] = self.get_v()
        d['nz'] = self.get_nz()
        d['ambient'] = self.ambient.as_dict()
        return dict_to_factor_sage(d)
    
    def from_sage(self, M, names='a'):
        d=factor_to_dict_sage(M)
        self.set_B(d['B'])
        self.set_Bd(d['Bd'])
        self.set_v(d['v'])
        self.set_nz(d['nz'])
        self.dimension = int(M.dimension())
        #self.has_cm = has_cm(M)
        # now we set the coefficient field
        extension_field = M.eigenvalue(1,name=names).parent()
        if extension_field != M.base_ring(): # .degree() != 1 and rings.is_NumberField(extension_field):
            assert extension_field.base_field() == M.base_ring()
            minpoly = extension_field.relative_polynomial()
            degree = int(minpoly.degree())
        else:
            minpoly = extension_field.defining_polynomial()
            degree = extension_field.degree()
        self.coefficient_field = CoefficientField()
        self.coefficient_field.minimal_polynomial = str(minpoly)
        self.coefficient_field.base_field = self.ambient.base_field
        self.coefficient_field.degree = degree

        if hasattr(M,'_HeckeModule_free_module__eigenvalues'):
            for n,c in M._HeckeModule_free_module__eigenvalues.iteritems():
                print n,c
                value = str(c[c.keys()[0]].list())
                print value
                cc = Coefficient(index=int(n))
                a = AlgebraicNumber()
                a.number_field = self.coefficient_field
                a.set_value(value)
                cc.value = a
                self.coefficients.append(cc)

class Coefficient(Coefficient_class, SageObject_class):
    def sage_object():
        return NotImplementedError()

class NumberField(NumberField_class, SageObject_DB_class):
    def sage_object(self):
        minpoly = self.minimal_polynomial
        

class ModularSymbols_base_field(ModularSymbols_base_field_class, SageObject_DB_class):
    def sage_object():
        return NotImplementedError()

class CoefficientField(CoefficientField_class, SageObject_DB_class):
    def sage_object():
        return NotImplementedError()

class AlgebraicNumber(AlgebraicNumber_class, SageObject_class):
    
    def as_vector(self):
        return sage_eval(self.get_value())
    
    def _set_value_from_sage(self,x):
        v = x.list()
        self.set_value(v)
