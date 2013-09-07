import inspect
import os
import bz2 as comp
from elixir import has_one,OneToMany,ManyToOne,belongs_to,has_many
from sqlalchemy.ext.associationproxy import AssociationProxy
from accessor import MixinAccessors
from sqlalchemy.ext.hybrid import hybrid_property
from mdb import db
# Types from db
Text = db.Text
String = db.String
Boolean = db.Boolean
Binary = db.Binary
# Note: We keep the namespace for db.Integer to avoid confusion with Sage's Integer.

compress = True
prefix = ''


class ModularSymbols_generic_DB_Mixin(object): #MixinAccessors): #,ModularSymbols_ambient_DB_class_base):
    r"""
    Ambient modular symbols space.
    """
#    __metaclass__ = Accessors
#    __tablename__ = 'ModularSymbols_ambient_DB_class'
    __table_args__ = {'useexisting': True}
    # primary keys
    level = db.Column(db.Integer, primary_key=True)
    weight = db.Column(db.Integer, primary_key=True)
    # the character is an integer following the Conrey naming scheme
    character = db.Column(db.Integer, default=int(0), primary_key=True)
    sign = db.Column(db.Integer, default=int(1))

    # the field of values of the character
    has_one('base_field', of_kind='ModularSymbols_base_field_DB_class')

    # dimensions for convenience
    dimension_modular_forms = db.Column(db.Integer, default=-1)
    dimension_cusp_forms = db.Column(db.Integer, default=-1)
    dimension_new_cusp_forms = db.Column(db.Integer, default=-1)

    # decomposition data
    is_ambient = db.Column(Boolean,default=False)
    is_newspace_factor = db.Column(Boolean,default=False)
    is_oldspace_factor = db.Column(Boolean,default=False)

    def __repr__(self):
        return 'Modular smybols generic space of level {0}, weight {1}, character {2} and dimension {3}'.format(
            self.level, self.weight, self.character, self.dimension_modular_forms)



#class ModularSymbols_ambient_DB_class(ModularSymbols_generic_DB_Mixin,db.Model):
class ModularSymbols_ambient_DB_class(ModularSymbols_generic_DB_Mixin,db.Model):
    r"""
    Ambient modular symbols space.
    """
    #__metaclass__ = Accessors
    __tablename__ = 'ModularSymbols_ambient_DB_class'
    __table_args__ = {'useexisting': True}
    compress = True
    # data to reconstruct the space
    # TODO: add documentation!
    __bfields = ['_basis', '_manin', '_rels', '_mod2term']
    if compress:
        _basis = db.Column(Binary, name = 'basis')
        _manin = db.Column(Binary, name = 'manin')
        _rels = db.Column(Binary,  name = 'rels')
        _mod2term = db.Column(Binary, name = 'mod2term')
        _COMPRESSED = __bfields
    else:
        _basis = db.Column(Text, name = 'basis')
        _manin = db.Column(Text, name = 'manin')
        _rels = db.Column(Text, name = 'rels')
        _mod2term = db.Column(Text, name = 'mod2term')
        _READ_WRITE = __bfields

    #__mapper_args__ = {
    #    'polymorphic_on' : _basis
    #    #('basis','manin,rels,mod2term)
   # }
    # the field of values of the character
    #    has_one('base_field', of_kind='{0}ModularSymbols_base_field_DB_class'.format(prefix))
    is_ambient = True
    # decomposition data
    newspace_factors = OneToMany('ModularSymbols_newspace_factor_DB_class')
    oldspace_factors = OneToMany('ModularSymbols_oldspace_factor_DB_class', inverse='ambient')
    oldspaces = AssociationProxy('oldspace_factors', 'factor',
                                 creator= lambda (factor,multiplicity):
                                     ModularSymbols_oldspace_factor(factor=factor,multiplicity=multiplicity))
    @hybrid_property
    def basis(self):
        return self._basis
    @basis.setter
    def basis(self,basis):
        if compress:
            self._basis = comp.compress(str(basis))
        else:
            self._basis = str(basis)
    @hybrid_property
    def manin(self):
        return self._manin
    @manin.setter
    def manin(self,manin):
        if compress:
            self._manin = comp.compress(str(manin))
        else:
            self._manin = str(manin)
    @hybrid_property
    def rels(self):
        return self._rels
    @rels.setter
    def rels(self,rels):
        if compress:
            self._rels = comp.compress(str(rels))
        else:
            self._rels = str(rels)
    @hybrid_property
    def mod2term(self):
        return self._mod2term
    @mod2term.setter
    def mod2term(self,mod2term):
        if compress:
            self._mod2term = comp.compress(str(mod2term))
        else:
            self._mod2term = str(mod2term)     

            
    def __repr__(self):
        return 'Modular smybols ambient space of level {0}, weight {1}, character {2} and dimension {3}'.format(
            self.level, self.weight, self.character, self.dimension_modular_forms)



class ModularSymbols_oldspace_factor_DB_class(ModularSymbols_generic_DB_Mixin):
    r"""
    An oldspace factor is a newspace `factor` that has an
    inclusion map into `ambient`. The number of different inclusion
    is the `multiplicity`
    """
    multiplicity = db.Column(db.Integer)
    ambient = ManyToOne('ModularSymbols_ambient_DB_class')
    factor = ManyToOne('ModularSymbols_ambient_DB_class')
    is_oldspace_factor = True
    #ambient = ManyToOne('{0}ModularSymbols_ambient_DB_class'.format(prefix))
    #factor = ManyToOne('{0}ModularSymbols_ambient_DB_class'.format(prefix))


    def __repr__(self):
        return 'Oldspace factor of level {0}, weight {1}, character {2}, dimension {3} with multiplicity {4}'\
            + ' of Modular forms ambient space of level {5}, weight {6} and character {7}' .format(
            self.factor.level, self.factor.weight, self.factor.character, self.factor.multiplicity,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)


class ModularSymbols_newspace_factor_DB_class(ModularSymbols_generic_DB_Mixin):
    r"""
    A single Galois orbit contained in `ambient`.
    """
    #__metaclass__ = Accessors

    is_newspace_factor = True

    #belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB_class'.format(prefix))
    belongs_to('ambient', of_kind='ModularSymbols_ambient_DB_class')
    # data to rectonstruct the ModularSymbols space

    __bfields = ['_B', '_Bd', '_v', '_nz']
    if compress:
        _B = db.Column(Binary, name = 'B')
        _Bd = db.Column(Binary, name = 'Bd')
        _v = db.Column(Binary, name = 'v')
        _nz = db.Column(Binary, name = 'nz')
        _COMPRESSED = __bfields
    else:
        _B = db.Column(Text, name = 'B')
        _Bd = db.Column(Text, name = 'Bd')
        _v = db.Column(Text, name = 'v')
        _nz = db.Column(Text, name = 'nz')
        _READ_WRITE = __bfields

    dimension = db.Column(db.Integer)
    has_cm = db.Column(Boolean) # has complex multiplication?
    #has_one('coefficient_field', of_kind='{0}Coefficientdb.Column_DB_class'.format(prefix))
    #has_many('coefficients', of_kind='{0}Coefficient_DB_class'.format(prefix))
    has_one('coefficient_field', of_kind='CoefficientField_DB_class')
    has_many('coefficients', of_kind='Coefficient_DB_class')

    def __repr__(self):
        return 'Newspace factor of dimension {0} of Modular forms ambient space of level {1}, weight {2}, character {3} and dimension {4}'.format(
            self.dimension,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)

class Coefficient_DB_class(db.Model):
    r"""
    A coefficient of a newform (factor).
    """
    __table_args__ = {'useexisting': True}
    belongs_to('newform', of_kind='ModularSymbols_newspace_factor_DB_class',primary_key=True)
    index = db.Column(db.Integer,primary_key=True)
    has_one('value', of_kind='AlgebraicNumber_DB_class',primary_key=True)


class NumberField_DB_class(db.Model):
    r"""
    A relative number field, represented by the minimal polynomial
    of a generator over its base field.
    """
    __table_args__ = {'useexisting': True}
    has_many('extensions', of_kind='{0}NumberField_DB_class'.format(prefix))
    belongs_to('base_field', of_kind='{0}NumberField_DB_class'.format(prefix))
    minimal_polynomial = db.Column(String, primary_key=True) # relative to the base field
    degree = db.Column(db.Integer) # relative to the base field
    is_cyclotomic = db.Column(Boolean)

    def __repr__(self):
        return 'Number field with minimal polynomial {0} over its base field.'.format(self.minimal_polynomial)

class ModularSymbols_base_field_DB_class(NumberField_DB_class):
    r"""
    The base field of a Moular symbols space.
    This will always be a cyclotomic field, determined by the field of values
    of the Dirichlet character.
    """
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB_class'.format(prefix))

class CoefficientField_DB_class(NumberField_DB_class):
    r"""
    The field of definition of a newform.
    """
    belongs_to('newspace', of_kind='{0}ModularSymbols_newspace_factor_DB_class'.format(prefix),primary_key=True)


class AlgebraicNumber_DB_class(db.Model,MixinAccessors):
    r"""
    An algebraic number is represented by a vector of coefficients
    in terms of a power basis of the number_field it is contained in.
    """
    #__metaclass__ = Accessors
    __tablename__ = 'AlgebraicNumber_DB_class'
    __table_args__ = {'useexisting': True}
    id = db.Column("id", db.Integer, primary_key=True)

    belongs_to('number_field', of_kind='{0}NumberField'.format(prefix),primary_key=True)
    belongs_to('coefficient', of_kind='{0}Coefficient_DB_class'.format(prefix),primary_key=True)

    # _value is the coefficient vector in terms of a power basis
    # of the number field (with specified minimal polynomial)
    __bfields = ['_value']
    if compress:
        _value = db.Column(Binary, name='value')
        _COMPRESSED = __bfields
    else:
        _value = db.Column(Text, name='value')
        _READ_WRITE = __bfields

    def __repr__(self):
        return 'Algebraic Number {0}, element of Number Fields with defining polynomial {1} over its base field.'.format(
            self.value, self.number_field.minimal_polynomial)
