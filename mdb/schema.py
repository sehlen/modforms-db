import inspect
import os
from elixir import *
from sqlalchemy.ext.associationproxy import AssociationProxy
prefix='schema.'
f = inspect.getabsfile(inspect.currentframe())
datadir = "/".join(os.path.dirname(f).rsplit("/")[0:-1]+["data"])
#metadata.bind = "sqlite:///data/modularforms.sqlite"
metadata.bind = "sqlite:///{0}/modularforms.sqlite".format(datadir)
metadata.bind.echo = True

class ModularSymbols_ambient(Entity):
    # primary key
    level = Field(Integer, primary_key=True)
    weight = Field(Integer, primary_key=True)
    character = Field(Integer, default=0, primary_key=True)
    
    # data to reconstruct the space
    # TODO: add documentation!
    basis = Field(Text)
    manin = Field(Text)
    rels = Field(Text)
    mod2term = Field(Text)

    # the field of values of the character
    has_one('base_field', of_kind='{0}ModularSymbols_base_field'.format(prefix))
    
    # dimensions for convenience
    dimension_modular_forms = Field(Integer, default=-1)
    dimension_cusp_forms = Field(Integer, default=-1)
    dimension_new_cusp_forms = Field(Integer, default=-1)

    # decomposition data
    newform_orbits = OneToMany('ModularSymbols_newspace_factor')
    oldspace_factors = OneToMany('ModularSymbols_oldspace_factor', inverse='ambient')
    oldspaces = AssociationProxy('oldspace_factors', 'factor',
                                        creator= lambda (factor,multiplicity):
                                        ModularSymbols_oldspace_factor(factor=factor,multiplicity=multiplicity))

class ModularSymbols_oldspace_factor(Entity):
    multiplicity = Field(Integer)
    ambient = ManyToOne('{0}ModularSymbols_ambient'.format(prefix))
    factor = ManyToOne('{0}ModularSymbols_ambient'.format(prefix))
    
class ModularSymbols_newspace_factor(Entity):
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient'.format(prefix))
    # data to rectonstruct the ModularSymbols space
    B = Field(Text)
    Bd = Field(Text)
    v = Field(Text)
    nz = Field(Text)
    #
    dimension = Field(Integer)
    has_cm = Field(Boolean) # has complex multiplication?
    has_one('coefficient_field', of_kind='{0}CoefficientField'.format(prefix))
    has_many('coefficients', of_kind='{0}Coefficient'.format(prefix))

class Coefficient(Entity):
    belongs_to('newform', of_kind='{0}ModularSymbols_newspace_factor'.format(prefix))
    index = Field(Integer)
    has_one('value', of_kind='{0}AlgebraicNumber'.format(prefix))

class NumberField(Entity):
    has_many('extensions', of_kind='{0}NumberField'.format(prefix))
    belongs_to('base_field', of_kind='{0}NumberField'.format(prefix))
    minimal_polynomial = Field(Text) # relative to the base field
    degree = Field(Integer) # relative to the base field
    is_cyclotomic = Field(Boolean)

class ModularSymbols_base_field(NumberField):
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient'.format(prefix))

class CoefficientField(NumberField):
    belongs_to('newspace', of_kind='{0}ModularSymbols_newspace_factor'.format(prefix))

class AlgebraicNumber(Entity):
    belongs_to('number_field', of_kind='{0}NumberField'.format(prefix))
    belongs_to('coefficient', of_kind='{0}Coefficient'.format(prefix))
    value = Field(Text) # this is the coefficient vector in terms of a power basis
    
