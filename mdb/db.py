r"""
File defining the database to use for holomorphic modular forms.
"""
import inspect,os
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
prefix='schema.'
f = inspect.getabsfile(inspect.currentframe())
datadir = "/".join(os.path.dirname(f).rsplit("/")[0:-1]+["data"])
print datadir

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///{0}/modularforms.sqlite".format(datadir)
db = SQLAlchemy(app)
db.session.bind.echo=True
