import os
import sys
import pdb
import logging
try:
    import simplejson
except ImportError:
    import json as simplejson

try:
    import sqlite as sqlite3
except ImportError:
    import sqlite3

class Admin(object):
    def __init__(self, jobname):
        self.job = jobname
        self.initlog()
        self.initdb()

    def initlog(self):
        logfile = "%s.txt" % self.job
##        if os.path.isfile(logfile):
##            exitwitherror()
        self.log = open(logfile, "a")

    def initdb(self):
        createtable = False
        db = "%s.db" % self.job
        if not os.path.isfile(db):
            createtable = True
        self.conn = sqlite3.connect(db)
        self.db = self.conn.cursor()
        if createtable:
            self.db.execute("""create table polymer
(smiles text primary key, gen integer, sequence text, json text)""")
            self.conn.commit()

    def storedata(self, key, gen, sequence="", json=""):
        self.db.execute("""insert into polymer
values ('%s', %d, '%s', '%s')""" % (key, gen, sequence, json))
        self.conn.commit()

    def getdata(self, key):
        self.db.execute("""select * from polymer where
smiles='%s'""" % key)
        fetch = self.db.fetchone()
        if fetch:
            return fetch[1:]
        else:
            return None

    def getalldata(self):
        self.db.execute("select * from polymer")
        fetch = self.db.fetchall()
        return fetch

    def deletedata(self, key):
        self.db.execute("""delete from polymer where
smiles='%s'""" % key)
        self.conn.commit()


def testdb(admin):
    pdb.set_trace()
    sys.exit(0)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        sys.exit(0)

    jobname = sys.argv[1]
    admin = Admin(jobname)
    testdb(admin)
