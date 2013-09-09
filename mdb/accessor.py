import bz2 as comp
#from elixir import EntityMeta
from db import db

MetaBase = type(db.Model)


def _addMethod(fldName, clsName, verb, methodMaker, dict):
    """Make a get or set method and add it to dict."""
    compiledName = _getCompiledName(fldName, clsName)
    methodName = _getMethodName(fldName, verb)
    dict[methodName] = methodMaker(compiledName)
    
def _getCompiledName(fldName, clsName):
    """Return mangled fldName if necessary, else no change."""
    # If fldName starts with 2 underscores and does *not* end with 2 underscores...
    if fldName[:2] == '__' and fldName[-2:] != '__':
        return "_%s%s" % (clsName, fldName)
    else:
        return fldName

def _getMethodName(fldName, verb):
    """'_salary', 'get'  => 'get_salary'"""
    s = fldName.lstrip('_') # Remove leading underscores
    return verb + '_' + s

def _makeGetterCompressed(compiledName):
    """Return a method that gets compiledName's value."""
    return lambda self: comp.decompress(str(self.__dict__[compiledName]))

def _makeSetterCompressed(compiledName):
    """Return a method that sets compiledName's value."""    
    return lambda self, value: setattr(self, compiledName, comp.compress(str(value)))

def _makeGetter(compiledName):
    """Return a method that gets compiledName's value."""
    return lambda self: self.__dict__[compiledName]

def _makeSetter(compiledName):
    """Return a method that sets compiledName's value."""    
    return lambda self, value: setattr(self, compiledName, value)

#class Accessors(MetaBase):
class MixinAccessors(object):
    """Adds accessor methods to a class."""
    def __new__(self):
        cls = self.__class__ # , clsName, bases, dict):
        dict = self.__dict__
        print "dict=",dict
        clsName = self.__name__
        for fldName in dict.get('_READ',[]) + dict.get('_READ_WRITE',[]):
            _addMethod(fldName, clsName, 'get', _makeGetter, dict)
        for fldName in dict.get('_WRITE', []) + dict.get('_READ_WRITE', []):
            _addMethod(fldName, clsName, 'set', _makeSetter, dict)
        for fldName in dict.get('_COMPRESSED',[]):
            #print 'adding get and set for ', fldName
            _addMethod(fldName, clsName, 'set', _makeSetterCompressed, dict)
            _addMethod(fldName, clsName, 'get', _makeGetterCompressed, dict)
#        return EntityMeta.__new__(cls, clsName, bases, dict)
#        return MetaBase.__new__(cls, clsName, bases, dict)
