#!@PYTHON@

#
# TODO: CLEAN THIS UP, merge packagepython.py and release.py 
#



import fnmatch
import sys
import os
import re
import string
import getopt
import pipes


_debug = 0

_prune = ['(*)']

def find(pattern, dir = os.curdir):
        list = []
        names = os.listdir(dir)
        names.sort()
        for name in names:
                if name in (os.curdir, os.pardir):
                        continue
                fullname = os.path.join(dir, name)
                if fnmatch.fnmatch(name, pattern):
                        list.append(fullname)
                if os.path.isdir(fullname) and not os.path.islink(fullname):
                        for p in _prune:
                                if fnmatch.fnmatch(name, p):
                                        if _debug: print "skip", `fullname`
                                        break
                        else:
                                if _debug: print "descend into", `fullname`
                                list = list + find(pattern, fullname)
        return list


topdir = ''
def gulp_file(f):
	try:
		i = open(f)
		i.seek (0, 2)
		n = i.tell ()
		i.seek (0,0)
	except:
		sys.stderr.write ("can't open file: %s\n" % f)
		return ''
	s = i.read (n)
	if len (s) <= 0:
		sys.stderr.write ("gulped emty file: %s\n" % f)
	i.close ()
	return s


def mailaddress():
	try:
		return os.environ['MAILADDRESS']
	except KeyError:
		return '(address unknown)'



class Flags:
	def __init__ (self):
		self.to_version = 0
		self.from_version = 0
		self.package = 0

flags = Flags ()

def help ():
	sys.stdout.write (
		'Generate a patch to go to current version\n'
		'  -f, --from=FROM      old is FROM\n'
		'  -h, --help           print this help\n'
		'      --outdir=DIR     generate in DIR\n'
		'  -o, --output=NAME    write patch to NAME\n'
		'  -p, --package=DIR    specify package\n'
		'  -r, --release        diff against latest release\n'  
		'  -t, --to=TO          to version TO\n'  
		'  -F, --dir-from=FROM  diff from directory FROM\n'  
		'  -T, --dir-to=TO      diff to directory TO\n'  
		)

def cleanup ():
	global from_diff, to_diff, original_dir
	os.chdir ('/tmp/package-diff')
	sys.stderr.write ('Cleaning ... ')
	os.system ('rm -fr %s %s' % (from_diff, to_diff))
	sys.stderr.write ('\n')
	os.chdir (original_dir)

def untar (fn):
	# os.system ('pwd');
	try:
		open (fn)
	except:
		sys.stderr.write ("Can't find tarball: %s\n" % fn)
		cleanup ()
		sys.exit (1)
	sys.stderr.write ("Untarring: %s\n" % fn)
	os.system ('gzip --quiet -dc ' + fn + '| tar xf - ')
	sys.stderr.flush ()

def remove_automatic (dirnames):
	files = []

	for d in dirnames:
		try:
			for p in pats:
				files = files + find (p, d)
		except:
			sys.stderr.write ("Can't find dir: %s\n" % d)
			cleanup ()
			sys.exit (1)

	dirs = map (lambda d: find ('out*', d), dirnames)
	dirs = reduce (lambda x,y:  x + y, dirs)
	
	#print dirs

	for d in dirs:
		if os.path.isdir (d):
			files = files + find ('*', d)
		
	for f in files:
		try:
			os.remove (f)
		except:
			sys.stderr.write ("can't remove: `" + f + "'\n'")

def dirname (v):
	return flags.package.name + '-' + version_tuple_to_str (v)

def tarball(v):
	return dirname (v)  + '.tar.gz'

def released_tarball(v):
	return flags.package.release_dir + tarball (v)


def remove_configure (dir):
	os.chdir (dir)

	# ugh
	os.system ('rm -f *.ly')

	# should do 'make distclean ?'
	os.system ('rm -rf debian/out GNUmakefile config.cache config.h config.hh config.log config.make config.status configure stepmake/GNUmakefile stepmake/config.hh stepmake/config.log stepmake/config.status stepmake/config.make')

	# ugh: symlinks
	os.system ('rm -f stepmake/stepmake/stepmake stepmake/stepmake/bin')


# ugh, how to get rid of .tex files generated by lily?
pats = ['configure', '*.diff', '*.diff.gz', '*.pyc',  '*.txt']

header = """Generated by %s,
From = %s, To = %s

usage 

    cd %s-source-dir; patch -E -p1 < %s

Patches do not contain automatically generated files 
or (urg) empty directories, 
i.e., you should run

	./autogen.sh
	./configure [your options.]

	
"""


def makediff (fromdir, todir, patch_name):
	remove_automatic ([fromdir, todir])
	
	# ugh
	remove_configure (fromdir)
	remove_configure (todir)
	os.chdir (todir)
	
	fromname = fromdir
	toname = todir
	if os.path.dirname (fromname) == os.path.dirname (toname):
		fromname = os.path.basename (fromname)
		toname = os.path.basename (toname)
		fromdir = '../' + fromname

	f = open (patch_name, 'wb')
	f.write (header %
		(mailaddress (),   fromname, toname, 
		 flags.package.name, os.path.basename (patch_name)))

	f.close ()
			
	sys.stderr.write ('diffing to %s... ' % patch_name)
	os.system ('pwd')
	print ('diff -purN %s . >> %s' % (fromdir, patch_name))
	os.system ('diff -purN %s . >> %s' % (fromdir, patch_name))
	os.system ('gzip --quiet -9f %s' % patch_name)
   

os.environ['GZIP'] = '-q'
print 'argv: ' + string.join (sys.argv[1:])
(options, files) = getopt.getopt (sys.argv[1:], 
	'hF:f:o:p:rT:t:', ['conf=', 'from=', 'dir-from=', 'dir-to=', 'help', 'outdir=', 'output=', 'package=', 'release', 'to='])

patch_name = ''
conf = ''
outdir = ''
from_src = ''
to_src = ''
release=0
from_version=0
to_version=0
for opt in options:
	o = opt[0]
	a = opt[1]
	if o == '--from' or o == '-f':
		 from_version = a
	elif o == '--to' or o == '-t':
		 to_version = a
	elif o == '--dir-from' or o == '-F':
		 from_src = a;
	elif o == '--dir-to' or o == '-T':
		 to_src = a;
	elif o == '--help' or o == '-h':
		help ()
		sys.exit (0)
	elif o == '--outdir':
		outdir = a
	elif o == '--conf':
		conf = a
	elif o == '--output' or o == '-o':
		patch_name = a
	elif o == '-p' or o == '--package':
		topdir = a
	elif o == '--release' or o == '-r':
		release=1
	else:
		raise getopt.error

sys.path.append (topdir + '/stepmake/bin')
from packagepython import *
package = Package (topdir)
flags.package = package
packager = Packager ()

if from_src:
	from_package = Package (from_src)
	flags.from_version = from_package.version
if from_version:
	flags.from_version = version_str_to_tuple (from_version)
	from_src = ''

if to_src:
	to_package = Package (to_src)
	flags.to_version = to_package.version
if to_version:
	flags.to_version = version_str_to_tuple (to_version)
	to_src = ''

if not flags.to_version:
	flags.to_version = package.version

if not flags.from_version:
	flags.from_version = prev_version (flags.to_version)

# urg
if release:
	flags.from_version = (flags.from_version[0], 
		flags.from_version[1], flags.from_version[2], '');

import tempfile
original_dir = os.getcwd ();

os.system ('rm -rf /tmp/package-diff') 
try:
	os.mkdir ('/tmp/package-diff')
except:
	pass

from_diff = dirname (flags.from_version)
to_diff =  dirname (flags.to_version)
if to_diff == from_diff:
	if from_src:
	      from_diff = from_diff + '.src'
	elif to_src:
	      to_diff = to_diff + '.src'
	else:
	      sys.stderr.write (patch_name + ': nothing to do: to == from = ' + from_diff + '\n')
	      sys.exit (1)

def compat_abspath (path):
	return os.path.normpath (os.path.join (os.getcwd (), path))

if conf and not outdir:
	outdir = 'out-' + conf

if not patch_name:
	patch_name = os.path.join (outdir, '%s-%s-%s.diff' % (package.name,
							      version_tuple_to_str (flags.from_version),
							      version_tuple_to_str (flags.to_version)))

	patch_name = compat_abspath (patch_name)

from_diff = '/tmp/package-diff/' + from_diff
to_diff =  '/tmp/package-diff/' + to_diff

if not from_src:
	os.chdir ('/tmp/package-diff')
	untar (released_tarball (flags.from_version))
	os.chdir (original_dir)
else:
	sys.stderr.write ('copying ' + from_src + ' to ' + from_diff + '\n')
	# os.system ('cp -pr %s %s' % (srcdir, from_diff))
	os.system ('mkdir -p %s '% (from_diff))
	os.chdir (from_src)
	os.system ('tar cf - --exclude out --exclude out-www . \
		| tar -xf - -C %s' % from_diff)


if not to_src:
	os.chdir ('/tmp/package-diff')
	untar (released_tarball (flags.to_version))
	os.chdir (original_dir)
else:
	sys.stderr.write ('copying ' + to_src + ' to ' + to_diff + '\n')
	os.system ('mkdir -p %s '% (to_diff))
	# os.system ('cp -pr %s %s'  (to_src, to_diff))%
	os.chdir (to_src)
	os.system ('tar -cf - --exclude out --exclude out-www . \
		. | tar -xf - -C %s ' % to_diff)

os.chdir (to_diff)
makediff (from_diff, to_diff, patch_name) 

cleanup ()

