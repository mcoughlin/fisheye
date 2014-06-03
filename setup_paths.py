
import os,sys,glob

dir = "/home/coughlin/.local/lib/python2.6/site-packages/"
eggs = glob.glob(os.path.join(dir,"*egg*"))

pythonpath = ""
for egg in eggs:
    pythonpath = "%s:%s"%(pythonpath,egg)

print pythonpath

