import inspect
import os
import bz2 as comp
#from elixir import has_one,ManyToOne,belongs_to,has_many
from sqlalchemy.ext.associationproxy import AssociationProxy
from sqlalchemy.ext.hybrid import hybrid_property
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy,DeclarativeMeta

from mdb import db
# Schemas for number fields.
from nf_schema import NumberField_class,get_number_field,AlgebraicNumber_class,NumberField_DB
from conversions import extract_args_kwds

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

compress = True
prefix = ''



class ModularSymbols_generic(db.Model):
    r"""
    Ambient modular symbols space.
    """
    __tablename__ = 'ModularSymbols_generic'
    __table_args__ = {'useexisting': True}
    # primary keys
    id = db.Column(db.Integer,primary_key=True)
    level = db.Column(db.Integer)
    weight = db.Column(db.Integer)
    # the character is an integer following the Conrey naming scheme
    character = db.Column(db.Integer)
    db.UniqueConstraint('level','weight','character', name='uix_1')

    sign = db.Column(db.Integer, default=int(1))

    # the field of values of the character
    #has_one('base_field', of_kind='ModularSymbols_base_field_DB_class')
    
    base_field = db.relationship('ModularSymbols_base_field_class')    
    base_field_id = db.Column(db.Integer,db.ForeignKey('NumberField_class.id'))   
    # dimensions for convenience
    dimension_modular_forms = db.Column(db.Integer, default=-1)
    dimension_cusp_forms = db.Column(db.Integer, default=-1)
    dimension_new_cusp_forms = db.Column(db.Integer, default=-1)

    # decomposition data
    is_ambient = db.Column(db.Boolean,default=False)
    is_newspace_factor = db.Column(db.Boolean,default=False)
    is_oldspace_factor = db.Column(db.Boolean,default=False)

    def __repr__(self):
        return 'Modular smybols generic space of level {0}, weight {1}, character {2} and dimension {3}'.format(
            self.level, self.weight, self.character, self.dimension_modular_forms)

#class ModularSymbols_ambient_DB_class(ModularSymbols_generic_DB_Mixin,db.Model):
class ModularSymbols_ambient_class(ModularSymbols_generic): #,db.Model):
    r"""
    Ambient modular symbols space.
    """
    #__metaclass__ = Accessors
    __tablename__ = 'ModularSymbols_ambient_class'
    __table_args__ = {'extend_existing': True}
    compress = True
    id = db.Column(db.Integer, db.ForeignKey('ModularSymbols_generic.id'), primary_key=True) 
    # data to reconstruct the space
    # TODO: add documentation!
    __bfields = ['_basis', '_manin', '_rels', '_mod2term']
    if compress:
        _basis = db.Column(db.Binary, name = 'basis')
        _manin = db.Column(db.Binary, name = 'manin')
        _rels = db.Column(db.Binary,  name = 'rels')
        _mod2term = db.Column(db.Binary, name = 'mod2term')
        _COMPRESSED = __bfields
    else:
        _basis = db.Column(db.Text, name = 'basis')
        _manin = db.Column(db.Text, name = 'manin')
        _rels = db.Column(db.Text, name = 'rels')
        _mod2term = db.Column(db.Text, name = 'mod2term')
        _READ_WRITE = __bfields

    is_ambient = True
    # decomposition data
    newspace_factors = db.relationship('ModularSymbols_newspace_factor_class',
                                        backref='ambient',
                                        primaryjoin="ModularSymbols_newspace_factor_class.ambient_id==ModularSymbols_ambient_class.id")
    
#    dummies = db.relationship('Dummy')
    _oldspace_factors = db.relationship('ModularSymbols_oldspace_factor_class', backref='ambient',
                                        primaryjoin="ModularSymbols_oldspace_factor_class.ambient_id==ModularSymbols_ambient_class.id")
    oldspace_factors = AssociationProxy('_oldspace_factors', 'factor',
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




class ModularSymbols_newspace_factor_class(ModularSymbols_generic):
    r"""
    A single Galois orbit contained in `ambient`.
    """

    __tablename__ = 'ModularSymbols_newspace_factor_class'
    __table_args__ = {'useexisting': True}
    is_newspace_factor = True
    id = db.Column(db.Integer, db.ForeignKey('ModularSymbols_generic.id'), primary_key=True)
    ambient_id = db.Column(db.Integer,db.ForeignKey('ModularSymbols_ambient_class.id'))
    # data to reconstruct the ModularSymbols space

    __bfields = ['_B', '_Bd', '_v', '_nz']
    if compress:
        _B = db.Column(db.Binary, name = 'B')
        _Bd = db.Column(db.Binary, name = 'Bd')
        _v = db.Column(db.Binary, name = 'v')
        _nz = db.Column(db.Binary, name = 'nz')
        _COMPRESSED = __bfields
    else:
        _B = db.Column(db.Text, name = 'B')
        _Bd = db.Column(db.Text, name = 'Bd')
        _v = db.Column(db.Text, name = 'v')
        _nz = db.Column(db.Text, name = 'nz')
        _READ_WRITE = __bfields

    dimension = db.Column(db.Integer)
    has_cm = db.Column(db.Boolean) # has complex multiplication?
    coefficient_field = db.relationship('CoefficientField_class',backref='newform')
    
    coefficients_list = db.relationship('Coefficient_class',backref='newform')    
    coefficients = AssociationProxy('coefficients_list', 'value',
                                    creator= lambda (index,value):
                                        Coefficient_DB(index=index,value=value))

    
    @hybrid_property
    def B(self):
        return self._B
    @B.setter
    def B(self,B):
        if compress:
            self._B = comp.compress(str(B))
        else:
            self._B = str(B)
    @hybrid_property
    def Bd(self):
        return self._Bd
    @Bd.setter
    def Bd(self,Bd):
        if compress:
            self._Bd = comp.compress(str(Bd))
        else:
            self._Bd = str(Bd)
    @hybrid_property
    def v(self):
        return self._v
    @v.setter
    def v(self,v):
        if compress:
            self._v = comp.compress(str(v))
        else:
            self._v = str(v)
    @hybrid_property
    def nz(self):
        return self._nz
    @nz.setter
    def nz(self,nz):
        if compress:
            self._nz = comp.compress(str(nz))
        else:
            self._nz = str(nz)

            
    def __repr__(self):
        return 'Newspace factor of dimension {0} of Modular forms ambient space of level {1}, weight {2}, character {3} and dimension {4}'.format(
            self.dimension,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)


    
class ModularSymbols_oldspace_factor_class(ModularSymbols_generic):
    r"""
    An oldspace factor is a newspace `factor` that has an
    inclusion map into `ambient`. The number of different inclusion
    is the `multiplicity`
    """
    id = db.Column(db.Integer, db.ForeignKey('ModularSymbols_generic.id'), primary_key=True)
    #ambient = db.relationship('ModularSymbols_ambient_class')
    ambient_id = db.Column(db.Integer,db.ForeignKey('ModularSymbols_ambient_class.id'))
    
    factor_id = db.Column(db.Integer,db.ForeignKey('ModularSymbols_newspace_factor_class.id'))
    factor = db.relationship('ModularSymbols_newspace_factor_class',
                             primaryjoin="ModularSymbols_oldspace_factor_class.factor_id==ModularSymbols_newspace_factor_class.id")
                           
    multiplicity = db.Column(db.Integer)
    is_oldspace_factor = True

    def __repr__(self):
        return 'Oldspace factor of level {0}, weight {1}, character {2}, dimension {3} with multiplicity {4}'\
            + ' of Modular forms ambient space of level {5}, weight {6} and character {7}' .format(
            self.factor.level, self.factor.weight, self.factor.character, self.factor.multiplicity,
            self.ambient.level, self.ambient.weight, self.ambient.character, self.ambient.dimension_modular_forms)


    
class Coefficient_class(db.Model):
    r"""
    A coefficient of a newform (factor).
    """
    __tablename__ = 'Coefficient_class'
    __table_args__ = {'useexisting': True}
    id = db.Column(db.Integer,primary_key=True)
    #belongs_to('newform', of_kind='ModularSymbols_newspace_factor_DB_class',primary_key=True)
    newform_id = db.Column(db.Integer,db.ForeignKey('ModularSymbols_newspace_factor_class.id'))
    index = db.Column(db.Integer)
    value = db.relationship('AlgebraicNumber_class',backref='coefficient')
    value_id = db.Column(db.Integer,db.ForeignKey('AlgebraicNumber_class.id'))
#has_one('value', of_kind='AlgebraicNumber_DB_class',primary_key=True)

    
def Coefficient_DB(*args,**kwds):
    print "args=",args
    print "kwds=",kwds
    kwds['value']=['index','value']
    n,val = extract_args_kwds(*args,**kwds)
    q = Coefficient_class.query.filter_by(index=n,value=val)
    if q.count()==1:
        return q.first()
    elif q.count()>1:
        raise ValueError,"Coefficient not uniquely described!"
    return Coefficient_class(index=n,value=val)
        

class ModularSymbols_base_field_class(NumberField_class):
    r"""
    The base field of a Moudlar symbols space.
    This will always be a cyclotomic field, determined by the field of values
    of the Dirichlet character.
    """
    #belongs_to('ambient', of_kind='{0}ModularSymbols_ambient_class'.format(prefix))
    id = db.Column(db.Integer,primary_key=True)
    number_field = db.relationship('NumberField_class')
    number_field_id = db.Column(db.Integer,db.ForeignKey('NumberField_class.id'))
    
class CoefficientField_class(NumberField_class):
    r"""
    The field of definition of a newform.
    """
    id = db.Column(db.Integer, primary_key=True)
    newform_id = db.Column(db.Integer, db.ForeignKey('ModularSymbols_newspace_factor_class.id'),unique=True)
    number_field = db.relationship('NumberField_class')
    number_field_id = db.Column(db.Integer,db.ForeignKey('NumberField_class.id'))    
    #belongs_to('newspace', of_kind='{0}ModularSymbols_newspace_factor_DB_class'.format(prefix),primary_key=True)


def ModularSymbols_base_field(*args,**kwds): #ModularSymbols_base_field_DB_class, SageObject_DB_class):
    return get_number_field(ModularSymbols_base_field_class,*args,**kwds)

def CoefficientField(*args,**kwds):
    return get_number_field(CoefficientField_class,*args,**kwds)
