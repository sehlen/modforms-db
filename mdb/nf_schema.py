import inspect
import os
import sage
import bz2 as comp
#from accessor import MixinAccessors
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy,DeclarativeMeta
from flask import ext
#from flask.ext.sqlalchemy.sqlalchemy.ext.associationproxy import AssociationProxy
AssociationProxy = ext.sqlalchemy.sqlalchemy.ext.associationproxy.AssociationProxy
hybrid_property = ext.sqlalchemy.sqlalchemy.ext.hybrid.hybrid_property 
#from flask.ext.sqlalchemy.declarative import declarative_base
from db import db

Text = db.Text
String = db.String
Boolean = db.Boolean
Binary = db.Binary
# Note: We keep the namespace for db.Integer to avoid confusion with Sage's Integer.

compress = True
prefix = ''

#
extensions = db.Table('extensions',
    db.Column('extension_id', db.Integer, db.ForeignKey('NumberField_DB.id'),primary_key=True),
    db.Column('base_field_id', db.Integer, db.ForeignKey('NumberField_DB.id'),primary_key=True),
                      info={'useexisting':True})



class NumberField_DB_class(db.Model):
    r"""
    A relative number field, represented by the minimal polynomial
    of a generator over its base field.
    """
    __tablename__ = 'NumberField_DB'
    __table_args__ = {'useexisting': True}
    id = db.Column(db.Integer,primary_key=True)
    extensions = db.relationship('NumberField_DB_class',
                                 secondary = extensions,
                                 backref = db.backref('base_field',lazy='joined'),  #remote_side=[id,base_field_id]) #,
                                 primaryjoin=id==extensions.c.extension_id,
                                 secondaryjoin=id==extensions.c.base_field_id)
    extensions_list = AssociationProxy('extensions_list','relative_polynomial',
                                    creator= lambda x:
                                        NumberField_DB(x))
    name = db.Column(String)
    relative_polynomial = db.Column(String,unique=True) #, primary_key=True)  # relative to the base field
    power_basis = db.Column(String) 
    degree = db.Column(db.Integer) # relative to the base field
    is_cyclotomic = db.Column(Boolean)
    algebraic_numbers = db.relationship('AlgebraicNumber_DB_class',backref='number_field',
                                       cascade = "all, delete-orphan")
    
    def __repr__(self):
        return 'Number field with minimal polynomial {0} over its base field.'.format(self.relative_polynomial)


class AlgebraicNumber_DB_class(db.Model):
    r"""
    An algebraic number is represented by a vector of coefficients
    in terms of a power basis of the number_field it is contained in.
    """
    __tablename__ = 'AlgebraicNumber_DB_class'
    __table_args__ = {'useexisting': True}
    id = db.Column("id", db.Integer, primary_key=True)
    value = db.Column(String)
    coefficient_id = db.Column(db.Integer, db.ForeignKey('Coefficient_DB_class.id'))
    number_field_id = db.Column(db.Integer, db.ForeignKey('NumberField_DB.id'))
    #algebraic_number = db.relationship('AlgebraicNumber_DB_class',backref='algebraic_number')
    # _value is the coefficient vector in terms of a power basis
    # of the number field (with specified minimal polynomial)
    #vector = db.relationship('vector_components')
    #def value(self):
    #    for x in self.vector:
            
        
    # __bfields = ['_value']
    # if compress:
    #     _value = db.Column(Binary, name='value')
    #     _COMPRESSED = __bfields
    # else:
    #     _value = db.Column(Text, name='value')
    #     _READ_WRITE = __bfields
    # @hybrid_property
    # def value(self):
    #     return self._value
    # @value.setter
    # def value(self,value):
    #     if compress:
    #         self._value = comp.compress(str(value))
    #     else:
    #         self._value = str(value)
            
    def __repr__(self):
        return 'Algebraic Number {0}, element of Number Fields with defining polynomial {1} over its base field.'.format(
            self.value, self.number_field.relative_polynomial)

#def vector_components(db.Model):
#    id = db.Column(db.Integer,primary_key=1)
#    index = db.Column(Integer)
#    value = db.Column(Integer)
    
## Constructor functions

def NumberField_DB(F,**kwds): #NumberField_DB_class, SageObject_DB_class):
    r"""
    Return a number field.
    """
    return get_number_field(NumberField_DB_class,F,**kwds)

def nf_dict_from_obj(obj,**kwds):
    try:
        print "obj=",obj
        res = dict(obj)
    except:
        res = {}
        try:
            res['degree'] = obj.degree()
        except (TypeError,AttributeError):
            try:
                res['degree'] = int(obj.degree)
            except (AttributeError,TypeError):
                pass
        try:
            res['relative_polynomial'] = obj.relative_polynomial()
            res['is_cyclotomic'] = obj.relative_polynomial().is_cyclotomic()
            res['power_basis'] = obj.power_basis()
        except (AttributeError, TypeError):
            try:
                res['relative_polynomial'] = obj.relative_polynomial
                res['is_cyclotomic'] = obj.is_cyclotomic
                res['power_basis'] = obj.power_basis
            except AttributeError:
                pass
    return res
    
def get_number_field(cls,data={},**kwds):
    r"""
    Returns a number field of class cls if in the database, otherwise creates and inserts it.
    """
    if not issubclass(cls,NumberField_DB_class):
        raise ValueError,"Expected a subclass of NumberField_DB_class! Got: {0}".format(cls)
    d = nf_dict_from_obj(data)
    print "d=",d
    d.update(kwds)
    degree = d.get('degree')
    rel_polynomial = d.get('relative_polynomial')    
    is_cyclotomic = d.get('is_cyclotomic',False)
    power_b = d.get('power_basis')
    if degree == 1:
        rel_polynomial = 'x'
    if degree == None or rel_polynomial == None:
        raise ValueError,"Could not create number field from {0} with keywords: {1}".format(data,kwds)
    power_b = str(power_b)
    rel_polynomial = str(rel_polynomial)
    rel_polynomial = rel_polynomial.replace(" ","") # normalize the polynomial
    print "\n minpol={0} \n".format(rel_polynomial)
    if rel_polynomial<>'':
        s = cls.query.filter_by(relative_polynomial=rel_polynomial)
        n = s.count()
        print "n=",n
        if s.count()>0:
            return s.first()
    res = cls(degree=int(degree),is_cyclotomic=is_cyclotomic,
              relative_polynomial=rel_polynomial,
              power_basis = power_b)
    if res not in db.session:
        db.session.add(res)
    if kwds.get('do_commit',1)==1:
        db.session.commit()
    print "con=",db.session.connection()
    
    if degree > 1: # We add base field if it doesn't exist and we have data.
        if isinstance(data,sage.rings.number_field.number_field_base.NumberField):
            base_field = data.base_field()
        else:
            base_field = data.get('base_field')
        base_field = get_number_field(cls,base_field)
        if res not in base_field.extensions:
            base_field.extensions.append(res)
            if kwds.get('do_commit',1)==1:
                db.session.commit()
    return res
    
    

def AlgebraicNumber_DB(data={},**kwds): 
    r"""
    Construct an algebraic number.
    """
    if isinstance(data,dict):
        nf = data['number_field']
        value = data['value']
    else:
        try:
            nf = data.parent()
            value = data.list()
        except (AttributeError,TypeError):
            try:
                nf = data.number_field
                value = data.value
            except AttributeError:
                raise ValueError,"Could not construct algebraic number from input!"
    NF = NumberField_DB(nf)
    value = str(value)
    print "NF=",NF,
    print "id=",NF.id
    print "value=",value
    q = AlgebraicNumber_DB_class.query.filter_by(number_field_id=NF.id,value=value)
    if q.count()>0:
        return q.first()
    else:
        newnumber = AlgebraicNumber_DB_class(value=value)
        db.session.add(newnumber)
        NF.algebraic_numbers.append(newnumber)
    if kwds.get('do_commit',1)==1:
        db.session.commit()
    return newnumber

    #     def as_vector(self):
#         return sage_eval(self.get_value())
    
#     def _set_value_from_sage(self,x):
#         v = x.list()
#         self.set_value(v)
    
