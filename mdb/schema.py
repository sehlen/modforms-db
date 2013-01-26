import inspect
import os
from elixir import *
from sqlalchemy.ext.associationproxy import AssociationProxy
import bz2 as comp

prefix='schema.'
f = inspect.getabsfile(inspect.currentframe())
datadir = "/".join(os.path.dirname(f).rsplit("/")[0:-1]+["data"])
#metadata.bind = "sqlite:///data/modularforms.sqlite"
metadata.bind = "sqlite:///{0}/modularforms.sqlite".format(datadir)
metadata.bind.echo = True

class ModularSymbols_ambient_DB(Entity):
    # primary key
    level = Field(Integer, primary_key=True)
    weight = Field(Integer, primary_key=True)
    # the character is an integer following the Conrey naming scheme
    character = Field(Integer, default=0, primary_key=True)
    
    # data to reconstruct the space
    # TODO: add documentation!
    basis = Field(Binary)
    manin = Field(Binary)
    rels = Field(Binary)
    mod2term = Field(Binary)

    # the field of values of the character
    has_one('base_field', of_kind='{0}ModularSymbols_base_field_DB'.format(prefix))
    
    # dimensions for convenience
    dimension_modular_forms = Field(Integer, default=-1)
    dimension_cusp_forms = Field(Integer, default=-1)
    dimension_new_cusp_forms = Field(Integer, default=-1)

    # decomposition data
    newform_orbits = OneToMany('ModularSymbols_newspace_factor_DB')
    oldspace_factors = OneToMany('ModularSymbols_oldspace_factor_DB', inverse='ambient')
    oldspaces = AssociationProxy('oldspace_factors', 'factor',
                                        creator= lambda (factor,multiplicity):
                                        ModularSymbols_oldspace_factor(factor=factor,multiplicity=multiplicity))

    def __repr__(self):
        return 'Modular smybols ambient space of level {0}, weight {1}, character {2} and dimension {3}'.format(
            self.level, self.weight, self.character, self.dimension_modular_forms)

class ModularSymbols_oldspace_factor_DB(Entity):
    multiplicity = Field(Integer)
    ambient = ManyToOne('{0}ModularSymbols_ambient_DB'.format(prefix))
    factor = ManyToOne('{0}ModularSymbols_ambient_DB'.format(prefix))
    
    def __repr__(self):
        return 'Oldspace factor of level {0}, weight {1}, character {2}, dimension {3} with multiplicity {4}'\
               + ' of Modular forms ambient space of level {5}, weight {6} and character {7}' .format(
            self.factor.level, self.factor.weight, self.factor.character, self.factor.multiplicity,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)

    
class ModularSymbols_newspace_factor_DB(Entity):
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB'.format(prefix))
    # data to rectonstruct the ModularSymbols space
    _B = Field(Binary, colname = 'B')
    _Bd = Field(Binary, colname = 'Bd')
    _v = Field(Binary, colname = 'v')
    _nz = Field(Binary, colname = 'nz')
    #
    @property
    def B(self):
        if self._B is not None:
            return comp.decompress(str(self._B))
        
    def set_B(self, B):
        if B is not None:
            self._B = comp.compress(str(B))

    @property
    def Bd(self):
        if self._Bd is not None:
            return comp.decompress(str(self._Bd))
        
    def set_Bd(self, Bd):
        if Bd is not None:
            self._Bd = comp.compress(str(Bd))

    @property
    def v(self):
        if self._v is not None:
            return comp.decompress(str(self._v))
        
    def set_v(self, v):
        if v is not None:
            self._v = comp.compress(str(v))

    @property
    def nz(self):
        if self._nz is not None:
            return comp.decompress(str(self._nz))
        
    def set_nz(self, nz):
        if nz is not None:
            self._nz = comp.compress(str(nz))
    #
    dimension = Field(Integer)
    has_cm = Field(Boolean) # has complex multiplication?
    has_one('coefficient_field', of_kind='{0}CoefficientField_DB'.format(prefix))
    has_many('coefficients', of_kind='{0}Coefficient_DB'.format(prefix))

    def __repr__(self):
        return 'Newspace factor of dimension {0} of Modular forms ambient space of level {1}, weight {2}, character {3} and dimension {4}'.format(
            self.dimension,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)

class Coefficient_DB(Entity):
    belongs_to('newform', of_kind='{0}ModularSymbols_newspace_factor_DB'.format(prefix))
    index = Field(Integer)
    has_one('value', of_kind='{0}AlgebraicNumber_DB'.format(prefix))

class NumberField_DB(Entity):
    has_many('extensions', of_kind='{0}NumberField_DB'.format(prefix))
    belongs_to('base_field', of_kind='{0}NumberField_DB'.format(prefix))
    minimal_polynomial = Field(String, primary_key=True) # relative to the base field
    degree = Field(Integer) # relative to the base field
    is_cyclotomic = Field(Boolean)

    def __repr__(self):
        return 'Number field with minimal polynomial {0} over its base field.'.format(self.minimal_polynomial)

class ModularSymbols_base_field_DB(NumberField_DB):
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB'.format(prefix))

class CoefficientField_DB(NumberField_DB):
    belongs_to('newspace', of_kind='{0}ModularSymbols_newspace_factor_DB'.format(prefix))

class AlgebraicNumber_DB(Entity):
    belongs_to('number_field', of_kind='{0}NumberField_DB'.format(prefix))
    belongs_to('coefficient', of_kind='{0}Coefficient_DB'.format(prefix))

    # _value is the coefficient vector in terms of a power basis
    # the 
    # it is stored with bz2 compression
    _value = Field(Binary, colname='value')

    @property
    def value(self, v=None):
        if v is not None:
            self._value = comp.compress(str(v))
        else:
            if self._value is not None:
                return comp.decompress(str(self._value))
    
    def __repr__(self):                    
        return 'Algebraic Number {0}, element of Number Field with defining polynomial {1} over its base field.'.format(
            self.value, self.number_field.minimal_polynomial)
