import inspect
import os
from elixir import *
from sqlalchemy.ext.associationproxy import AssociationProxy
from accessor import Accessors
from sqlalchemy.ext.declarative import declarative_base

prefix='schema.'
f = inspect.getabsfile(inspect.currentframe())
datadir = "/".join(os.path.dirname(f).rsplit("/")[0:-1]+["data"])
print datadir
# Commented out:
#metadata.bind = "sqlite:///{0}/modularforms.sqlite".format(datadir)
#metadata.bind.echo = True


compress = True

from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy,DeclarativeMeta

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///{0}/modularforms.sqlite".format(datadir)
db = SQLAlchemy(app)

class MyMeta(Accessors,DeclarativeMeta):
    pass

class ModularSymbols_ambient_DB(db.Model):
    r"""
        Ambient modular symbols space.
    """
    __metaclass__ = MyMeta #Accessors
    
    # primary key
    level = Field(Integer, primary_key=True)
    weight = Field(Integer, primary_key=True)
    # the character is an integer following the Conrey naming scheme
    character = Field(Integer, default=0, primary_key=True)
    sign = Field(Integer, default=1)
    
    # data to reconstruct the space
    # TODO: add documentation!
    __bfields = ['_basis', '_manin', '_rels', '_mod2term']
    if compress:
        _basis = Field(Binary, colname = 'basis')
        _manin = Field(Binary, colname = 'manin')
        _rels = Field(Binary, colname = 'rels')
        _mod2term = Field(Binary, colname = 'mod2term')
        _COMPRESSED = __bfields
    else:
        _basis = Field(Text, colname = 'basis')
        _manin = Field(Text, colname = 'manin')
        _rels = Field(Text, colname = 'rels')
        _mod2term = Field(Text, colname = 'mod2term')
        _READ_WRITE = __bfields

    # the field of values of the character
    has_one('base_field', of_kind='{0}ModularSymbols_base_field_DB'.format(prefix))
    
    # dimensions for convenience
    dimension_modular_forms = Field(Integer, default=-1)
    dimension_cusp_forms = Field(Integer, default=-1)
    dimension_new_cusp_forms = Field(Integer, default=-1)

    # decomposition data
    newspace_factors = OneToMany('ModularSymbols_newspace_factor_DB')
    oldspace_factors = OneToMany('ModularSymbols_oldspace_factor_DB', inverse='ambient')
    oldspaces = AssociationProxy('oldspace_factors', 'factor',
                                        creator= lambda (factor,multiplicity):
                                        ModularSymbols_oldspace_factor(factor=factor,multiplicity=multiplicity))

    def __repr__(self):
        return 'Modular smybols ambient space of level {0}, weight {1}, character {2} and dimension {3}'.format(
            self.level, self.weight, self.character, self.dimension_modular_forms)

class ModularSymbols_oldspace_factor_DB(db.Model):
    r"""
        An oldspace factor is a newspace `factor` that has an
        inclusion map into `ambient`. The number of different inclusion
        is the `multiplicity`
    """
    multiplicity = Field(Integer)
    ambient = ManyToOne('{0}ModularSymbols_ambient_DB'.format(prefix))
    factor = ManyToOne('{0}ModularSymbols_ambient_DB'.format(prefix))
    
    def __repr__(self):
        return 'Oldspace factor of level {0}, weight {1}, character {2}, dimension {3} with multiplicity {4}'\
               + ' of Modular forms ambient space of level {5}, weight {6} and character {7}' .format(
            self.factor.level, self.factor.weight, self.factor.character, self.factor.multiplicity,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)

    
class ModularSymbols_newspace_factor_DB(db.Model):
    r"""
        A single Galois orbit contained in `ambient`.
    """
    __metaclass__ = Accessors
    

    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB'.format(prefix))
    # data to rectonstruct the ModularSymbols space

    __bfields = ['_B', '_Bd', '_v', '_nz']
    if compress:
        _B = Field(Binary, colname = 'B')
        _Bd = Field(Binary, colname = 'Bd')
        _v = Field(Binary, colname = 'v')
        _nz = Field(Binary, colname = 'nz')
        _COMPRESSED = __bfields
    else:
        _B = Field(Text, colname = 'B')
        _Bd = Field(Text, colname = 'Bd')
        _v = Field(Text, colname = 'v')
        _nz = Field(Text, colname = 'nz')
        _READ_WRITE = __bfields

    dimension = Field(Integer)
    has_cm = Field(Boolean) # has complex multiplication?
    has_one('coefficient_field', of_kind='{0}CoefficientField_DB'.format(prefix))
    has_many('coefficients', of_kind='{0}Coefficient_DB'.format(prefix))

    def __repr__(self):
        return 'Newspace factor of dimension {0} of Modular forms ambient space of level {1}, weight {2}, character {3} and dimension {4}'.format(
            self.dimension,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)

class Coefficient_DB(db.Model):
    r"""
        A coefficient of a newform (factor).
    """
    belongs_to('newform', of_kind='{0}ModularSymbols_newspace_factor_DB'.format(prefix))
    index = Field(Integer)
    has_one('value', of_kind='{0}AlgebraicNumber_DB'.format(prefix))

class NumberField_DB(db.Model):
    r"""
        A relative number field, represented by the minimal polynomial
        of a generator over its base field.
    """
    has_many('extensions', of_kind='{0}NumberField_DB'.format(prefix))
    belongs_to('base_field', of_kind='{0}NumberField_DB'.format(prefix))
    minimal_polynomial = Field(String, primary_key=True) # relative to the base field
    degree = Field(Integer) # relative to the base field
    is_cyclotomic = Field(Boolean)

    def __repr__(self):
        return 'Number field with minimal polynomial {0} over its base field.'.format(self.minimal_polynomial)

class ModularSymbols_base_field_DB(NumberField_DB):
    r"""
        The base field of a Moular symbols space.
        This will always be a cyclotomic field, determined by the field of values
        of the Dirichlet character.
    """
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_DB'.format(prefix))

class CoefficientField_DB(NumberField_DB):
    r"""
        The field of definition of a newform.
    """
    belongs_to('newspace', of_kind='{0}ModularSymbols_newspace_factor_DB'.format(prefix))

class AlgebraicNumber_DB(db.Model):
    r"""
        An algebraic number is represented by a vector of coefficients
        in terms of a power basis of the number_field it is contained in.
    """
    __metaclass__ = Accessors
    
    belongs_to('number_field', of_kind='{0}NumberField_DB'.format(prefix))
    belongs_to('coefficient', of_kind='{0}Coefficient_DB'.format(prefix))

    # _value is the coefficient vector in terms of a power basis
    # of the number field (with specified minimal polynomial)
    __bfields = ['_value']
    if compress:
        _value = Field(Binary, colname='value')
        _COMPRESSED = __bfields
    else:
        _value = Field(Text, colname='value')
        _READ_WRITE = __bfields
    
    def __repr__(self):                    
        return 'Algebraic Number {0}, element of Number Field with defining polynomial {1} over its base field.'.format(
            self.value, self.number_field.minimal_polynomial)
