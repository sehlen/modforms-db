from elixir import *
from sqlalchemy.ext.associationproxy import AssociationProxy
prefix='schema.'

metadata.bind = "sqlite:///modularforms.sqlite"
metadata.bind.echo = True

class ModularSymbols_ambient(Entity):
    level = Field(Integer, primary_key=True)
    weight = Field(Integer, primary_key=True)
    character = Field(Integer, default=0, primary_key=True)
    dimension_modular_forms = Field(Integer, default=-1)
    dimension_cusp_forms = Field(Integer, default=-1)
    dimension_new_cusp_forms = Field(Integer, default=-1)
    number_of_new_orbits = Field(Integer, default=-1)
    oldspaces_associations = ManyToMany('OldspaceAssociation')
#    newspaces_associations = OneToMany('OldspaceAssociation', inverse='ambient')
    oldspace_factors = AssociationProxy('oldspaces_associations', 'factor',
                                        creator= lambda factor: OldspaceAssociation(factor=factor,ambient=self))
#    constituent_of = AssociationProxy('oldspaces_associations', 'ambient',
#                                        creator= lambda ambient: OldspaceAssociation(ambient=ambient))
#    has_many('oldspace_factors', through='{0}Modularsymbols_oldspace_factor'.format(prefix),
#             via='{0}ModularSymbols_ambient'.format(prefix))
    has_one('base_field', of_kind='{0}ModularSymbols_base_field'.format(prefix))

class OldspaceAssociation(Entity):
    multiplicity = Field(Integer)
    ambient = ManyToOne('{0}ModularSymbols_ambient'.format(prefix))
    factor = ManyToOne('{0}ModularSymbols_ambient'.format(prefix))

#class ModularSymbols_oldspace_factor(Entity):
#    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient'.format(prefix))
#    has_one('factor', of_kind='{0}ModularSymbols_ambient'.format(prefix))
#    multiplicity = Field(Integer)
    
class ModularSymbols_newspace_factor(Entity):
    belongs_to('ambient', of_kind='{0}ModularSymbols_ambient'.format(prefix))
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
    
