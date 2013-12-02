#!/usr/bin/env python3

"""
Parameter file format

Parameters for tasks and programs follow a hierarchical namespace, a tree 
similar to a directory structure, but with segments separated by the period 
character: "."
E.g., the parameters of the task "foo" live under
tasks.foo
Any parameters specified in the path to a node are considered when looking 
for a parameter in a node; parameters late in the path take precedence. 
Hence with:

threads = 2
tasks.sieve.threads = 1

the sieving task will use the value 1 for the threads parameter, while other 
tasks use the value 2 (unless they specify a different value in their subtree).

Tasks run programs, and those have their own node in the parameter tree.
For example, with
threads = 2
tasks.sieve.las.threads = 1
the theads=1 parameter would apply only to the las program, but not to any 
other programs run during the sieving tasks, if there were any. The name of 
the node of a program is equal to the name of the binary executable.
"""

import os
import re
import abc
import logging
import cadologger

logger = logging.Logger("Parameters")

parse_array = []

class BoolParam(object):
    def __init__(self, value):
        if value is True or isinstance(value, str) and \
                value.lower() in ["yes", "true", "on", "1"]:
            self.value = True
        elif value is False or isinstance(value, str) and \
                value.lower() in ["no", "false", "off", "0"]:
            self.value = False
        else:
            raise ValueError("Could not parse %s as truth value" % value)
    
    def __bool__(self):
        return self.value

    def __str__(self):
        return str(self.value)

class Parameters(object):
    """ Class that stores parameters for cadofactor in hierarchical dictionaries
    """
    # TODO: Add a self.used dictionary with the same structure and keys as
    # self.data, but booleans as the values which indicate whether a parameter
    # has ever been accessed with myparams(). Add a function that tests that
    # all parameters have been accessed, so that a warning can be printed
    # about parameters in the parameter file that are not used by anything,
    # which might indicate a misspelling, etc.
    
    old_key_map =  {
        "alambda": "tasks.sieve.alambda",
        "alim": "alim",
        "bindir": "tasks.execpath",
        "bwc_interleaving": "tasks.linalg.bwc.interleaving", 
        "bwc_interval": "tasks.linalg.bwc.interval",
        "bwc_mm_impl": None, # This parameter seems to be gone
        "bwc_mn": "tasks.linalg.bwc.mn",
        "bwc_shuffled_product": "tasks.linalg.bwc.shuffled_product",
        "bwmt": "tasks.linalg.bwc.threads",
        "bwstrat": "tasks.filter.merge.forbw",
        "checkrange": None,
        "coverNmax": "tasks.merge.coverNmax",
        "degree": "tasks.polyselect.degree",
        "delay": None,
        "dup_rm": "tasks.filter.duplicates2.rm",
        "excesspurge": "excesspurge",
        "filterlastrels": None, # FIXME: implement this behaviour
        "firstcheck": "tasks.sieve.rels_wanted",
        "hosts": None, # FIXME: how to implement this
        "I": "tasks.sieve.I",
        "keep": "tasks.filter.merge.keep",
        "keeppurge": "tasks.filter.purge.keep",
        "keeprelfiles": None,
        "linalg" : None, # Should we allow different linalg packages?
        "lpba": "lpba",
        "lpbr": "lpbr",
        "machines": None, # FIXME: how to handle this?
        "maxlevel": "tasks.filter.maxlevel",
        "mfba": "tasks.sieve.mfba",
        "mfbr": "tasks.sieve.mfbr",
        "mpi": None, # FIXME: Implement this
        "n": "N",
        "name": "name",
        "nchar": "tasks.linalg.characters.nchar",
        "nkermax": None, # FIXME: implement this
        "nslices_log": "tasks.filter.nslices_log",
        "nthchar": "tasks.linalg.characters.threads",
        "poly_max_threads": "tasks.polyselect.threads", 
        "parallel": None, # We (currently) always use client/server
        "polsel_admax": "tasks.polyselect.admax",
        "polsel_admin": "tasks.polyselect.admin",
        "polsel_adrange": "tasks.polyselect.adrange",
        "polsel_delay": None,
        "polsel_incr": "tasks.polyselect.incr",
        "polsel_maxnorm": "tasks.polyselect.maxnorm",
        "polsel_nice": "slaves.niceness", # polsel_nice and sievenice overwrite each other
        "polsel_nq": "tasks.polyselect.nq",
        "polsel_P": "tasks.polyselect.P",
        "qmin": "tasks.sieve.qmin",
        "qrange": "tasks.sieve.qrange",
        "ratio": "tasks.filter.ratio",
        "ratq": "tasks.sieve.ratq",
        "rlambda": "tasks.sieve.rlambda",
        "rlim": "rlim",
        "scriptpath": "slaves.scriptpath",
        "serveraddress": "server.address",
        "sieve_max_threads": "tasks.sieve.threads",
        "sievenice": "slaves.niceness", # polsel_nice and sievenice overwrite each other
        "slaves": "slaves.hostnames",
        "skip": "tasks.purge.skip",
        "wdir": "tasks.workdir",
        "expected_factorization": None
    }
    
    key_types = {
        "admin": int,
        "admax": int,
        "adrange": int,
        "qmin": int,
        "qmax": int,
        "qrange": int,
        "maxwu": int,
        "alim": int,
        "rlim": int,
        "rels_wanted": int,
        "nslices_log": int,
        "lognrels": int,
        "skip": int,
        "N": int,
        "nrclients": int,
        "threads": int,
        "port": int,
        "verbose": BoolParam,
        "rm": BoolParam,
        "quite": BoolParam,
        "sizeonly": BoolParam,
        "nopowers": BoolParam,
        "ratq": BoolParam,
        "bzip": BoolParam,
        "raw": BoolParam,
        "complete": BoolParam,
        "wipeout": BoolParam,
        "dryrun": BoolParam,
        "keepoldresult": BoolParam,
        "compression": BoolParam,
        "ssl": BoolParam,
        "threaded": BoolParam,
    }

    def __init__(self, *args, **kwargs):
        self.data = dict(*args, **kwargs)
        self._have_read_defaults = False
    
    def myparams(self, keys, path):
        ''' From the hierarchical dictionary params, generate a flat 
        dictionary with those parameters which are listed in keys and that 
        are found along path. 
        path is specified as a string with path segments separated '.'
        
        >>> d = {'a':1,'b':2,'c':3,'foo':{'a':3},'bar':{'a':4,'baz':{'a':5}}}
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'foo')
        {'a': 3, 'b': 2}
        
        >>> Parameters(d).myparams(keys=('a', 'b'), path = 'bar.baz')
        {'a': 5, 'b': 2}
        '''
        # path can be an array of partial paths, i.e., each entry can contain
        # one or more path segments separated by '.'. First join
        # them all, then split them again
        if isinstance(path, str):
            joinpath = path
        else:
            joinpath = '.'.join(path)
        splitpath = joinpath.split('.')
        
        source = self.data
        result = self._extract_by_keys(source, keys)
        for segment in splitpath:
            if not segment in source:
                break
            source = source[segment]
            result.update(self._extract_by_keys(source, keys))
        return result
    
    @staticmethod
    def _extract_by_keys(source, keys):
        return {key:source[key] for key in keys 
                if key in source and not isinstance(source[key], dict)}
    
    def __iter__(self):
        return self._recurse_iter(self.data, [])
    
    @staticmethod
    def _recurse_iter(source, path):
        for key in source:
            if isinstance(source[key], dict):
                for y in Parameters._recurse_iter(source[key], path + [key]):
                    yield y
            else:
                yield (path, key, source[key])
    
    def _get_subdict(self, path):
        source = self.data
        for d in path:
            if not d in source:
                return None
            assert isinstance(source[d], dict)
            source = source[d]
        return source
    
    def find(self, path, regex):
        source = self._get_subdict(path)
        if not source:
            source = {}
        pattern = re.compile(regex)
        result = [[l[0], l[1]] for l in self._recurse_iter(source, path)
            if pattern.search(l[1])]
        return result
    
    def _insertkey(self, path, value):
        ''' path is a path with segments delimited by '.' or an 
        array of pieces of the path, 
        value is inserted in the hierarchical parameters dictionary at the 
        location specified by path
        Keys overwrite previously existing keys, but a conflict between a key
        and a sub-dictionary causes a KeyError (similar to open('filepath','w')
        when filepath exists as a subdirectory).
        
        >>> p=Parameters({'a':1, 'b':2, 'foo':{'a':3}, 'bar':{'a':4}})
        >>> p._insertkey('c', 3)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.c', 5)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.baz.c', 6)
        >>> p.data == {'a': 1, 'b': 2, 'c': 3, 'bar': {'a': 4, 'c': 5, 'baz': {'c': 6}}, 'foo': {'a': 3}}
        True
        >>> p._insertkey('bar.baz.c.d', 6)
        Traceback (most recent call last):
        KeyError: 'Subdirectory c already exists as key'
        >>> p._insertkey('bar', 6)
        Traceback (most recent call last):
        KeyError: 'Key bar already exists as subdictionary'
        '''
        
        if isinstance(path, str):
            joinpath = path
        else:
            joinpath = '.'.join(path)
        splitpath = joinpath.split('.')
        key = splitpath.pop()
        dest = self.data
        for segment in splitpath:
            if segment in dest.keys():
                if not isinstance(dest[segment], dict):
                    raise KeyError('Subdirectory %s already exists as key' 
                                   % segment)
                dest = dest[segment]
            else:
                dest[segment] = {}
                dest = dest[segment]
        if key in dest.keys() and isinstance(dest[key], dict):
            raise KeyError('Key %s already exists as subdictionary' % key)
        dest[key] = value
    
    @staticmethod
    def subst_env_var(fqn, value):
        """ Substitute strings like '${HOME}' in value by the corresponding
        shell environment variables
        """
        while True:
            match = re.search(r"^(.*)\$\{(.*)\}(.*)$", value)
            if not match:
                break
            (prefix, varname, postfix) = match.groups()
            if not varname in os.environ:
                raise KeyError('Shell environment variable ${%s} referenced '
                               'in key %s is not defined (maybe not exported?)'
                               % (varname, fqn))
            value = prefix + os.environ[varname] + postfix
        return value
    
    def _subst_reference(self, path, key, value):
        """ Substitute strings like '$(somekey)' in a value by the value of
        "somekey" found along the current path, e.g.,
        foo.bar.k = $(m)
        foo.m = 5
        k = $(m)
        m = 3
        results in foo.bar.k = 5 and k = 3
        """
        while isinstance(value, str):
            match = re.search(r"^(.*)\$\((.*)\)(.*)$", value)
            if not match:
                break
            (prefix, varname, postfix) = match.groups()
            if key == varname:
                raise KeyError("Self-referential substitution $(%s) in key %s")
            result = self.myparams([varname], path)
            if not result:
                raise KeyError('Key $(%s) referenced in key %s is not defined'
                               % (varname, '.'.join(path + [key])))
            value = prefix + result[varname] + postfix
        return value
    
    def _subst_references(self, dic, path):
        for key in dic:
            if isinstance(dic[key], dict):
                self._subst_references(dic[key], path + [key])
            else:
                dic[key] = self._subst_reference(path, key, dic[key])
    
    def translate_old_key(self, key):
        """ If key is in the translation table, translate it.
        """
        # If allow_new is True, we allow new-style parameters to occur, too,
        # simply by not translating anything that does not occur in the
        # translation table.
        allow_new = True
        if allow_new:
            return self.old_key_map.get(key, key)
        else:
            return self.old_key_map[key]

    def _convert_one_type(self, path, key, orig_value):
        """ If a particular data type is registered in for this parameter,
        convert its value to that type, otherwise return the value unchanged.
        If the dataype is int, and the characters 'e' or '.' occur in the value,
        we convert to float first to convert scientific notation.
        """
        # If this value was converted before, don't try to convert again, as
        # the code below assumes it is a str
        if not isinstance(orig_value, str):
            return orig_value
        param = ".".join(path + [key])
        # print ("Trying to convert param=%s, value=%s" % (param, orig_value))
        datatype = self.key_types.get(key, None)
        value = orig_value
        if datatype is None:
            return value
        if datatype is int and ("e" in value or "." in value):
            value = float(value)
            if float(int(value)) != value:
                raise ValueError("Value %s for parameter %s cannot be "
                                 "converted to int without loss"
                                 % (orig_value, param))
        try:
            value = datatype(value)
            # print ("Converted param=%s, value=%r to %r, type %s" %
            #       (param, orig_value, value, datatype.__name__))
        except ValueError as err:
            ValueError("Cannot convert value %s for parameter %s to type %s"
                       % (value, param, datatype.__name__))
        return value

    def _convert_types(self, dic, path):
        for key in dic:
            if isinstance(dic[key], dict):
                self._convert_types(dic[key], path + [key])
            else:
                dic[key] = self._convert_one_type(path, key, dic[key])

    def parseline(self, line, old_format):
        (line2, comment) = line.split('#', 1) if '#' in line else (line, None)
        line2 = line2.strip()
        if not line2:
            return (None, None, comment, None)
        if not '=' in line2:
            raise Exception('Invalid line, missing "=": %s' % line)
        # Which one is worse?
        # (key, value) = re.match(r'(\S+)\s*=\s*(\S+)', line).groups()
        (key, value) = (s.strip() for s in line2.split('=', 1))
        oldkey = None
        if old_format:
            oldkey = key
            key = self.translate_old_key(key)
            if key is None:
                value = None
        return (key, value, comment, oldkey)

    def readparams(self, infile, old_format = False):
        """ 
        Read configuration file lines from infile, which must be an iterable.
        An open file handle, or an array of strings, work.
        
        >>> p = Parameters()
        >>> p.readparams(DEFAULTS_OLD, True)
        >>> p.data["tasks"]["sieve"]["rels_wanted"]
        1
        >>> p.myparams(["degree", "incr"], "tasks.polyselect") == \
        {'incr': '60', 'degree': '5'}
        True
        """
        for line in infile:
            line = line.strip('\n')
            (key ,value, comment, oldkey) = self.parseline(line, old_format)
            # print ("%s = %s # %s", (key ,value, comment))
            parse_array.append((key ,value, comment, oldkey))
            if key is None:
                continue
            value = self.subst_env_var(key, value)
            self._insertkey(key, value)
        self._subst_references(self.data, [])
        self._convert_types(self.data, [])

    def read_old_defaults(self):
        """ Read the DEFAULTS_OLD parameter table to set default values as they
        were set by the Perl script, as some of the old parameter files may
        assume those defaults being effective.
        """
        logger.debug("Reading old default parameters table")
        self.readparams(DEFAULTS_OLD, True)

    def readfile(self, filename, old_format = False):
        """ Read parameters from a file """
        logger.debug("Reading parameter file %s", filename)
        with open(filename, "r") as handle:
            self.readparams(handle, old_format)

    def __str_internal__(self):
        ''' Returns all entries of the dictionary dic as key=sep strings
        in an array
        '''
        return ("%s = %s" % (".".join(path + [key]), value) for
                (path, key, value) in self)
    
    def __str__(self):
        r = self.__str_internal__()
        return "\n".join(r)


class UseParameters(metaclass=abc.ABCMeta):
    """ Mix-in class for objects that take parameters from the parameter file
    """
    @abc.abstractproperty
    def param_nodename(self):
        """ The name of this object's node in the parameter tree """
        pass
    
    class MyParameters():
        """ Class that encapsules info on this node's parameters
        
        It stores a reference to the parameter dictionaries and info on this
        node's path and node name in the parameter tree.
        It can be initialised from a MyParameters or a Parameters instance;
        in the latter case the path prefix defaults to empty. If path_prefix
        is specified, it is always used as the path prefix.
        """
        def __init__(self, parent, name, path_prefix = None):
            self.name = name
            if isinstance(parent, UseParameters.MyParameters):
                self.parameters = parent.parameters
                self.path_prefix = parent.get_param_path()
            else:
                # Assume parent is a Parameters instance
                self.parameters = parent
                self.path_prefix = []
            if not path_prefix is None:
                self.path_prefix = path_prefix
        
        def get_param_path(self):
            if self.name is None:
                return self.path_prefix[:]
            else:
                return self.path_prefix + [self.name]
        
        def get_parameters(self):
            return self.parameters
        
        def myparams(self, keys, extrapath = None):
            path = self.get_param_path()
            if not extrapath is None:
                if isinstance(extrapath, str):
                    path += [extrapath]
                else:
                    path += extrapath
            return self.parameters.myparams(keys, path)
        
        def __str__(self):
            return "MyParameters: name = %s, prefix = %s" % (
                self.name, self.path_prefix)
        
    
    def __init__(self, *args, parameters, path_prefix = None,
                 **kwargs):
        self.parameters = UseParameters.MyParameters(
            parameters, self.param_nodename, path_prefix)
        super().__init__(*args, **kwargs)


DEFAULTS_OLD = (
    # global
    #'wdir         = undef',
    #'bindir      = undef',
    #'name         = undef',
    #'machines     = undef',
    #'n                = undef',
    'parallel     = 0',

    # polyselect using Kleinjung's algorithm
    'degree         = 5',
    'polsel_nq      = 1000',
    'polsel_incr    = 60',
    # 'polsel_admin   = 0', # 0 is default anyway
    #'polsel_admax   = undef',
    'polsel_adrange = 1e7',
    'polsel_delay   = 120',
    #'polsel_P       = undef',
    'polsel_maxnorm = 1e9',
    'polsel_nice    = 10',

    # sieve
    'rlim         = 8000000',
    'alim         = 8000000',
    'lpbr         = 29',
    'lpba         = 29',
    'mfbr         = 58',
    'mfba         = 58',
    'rlambda      = 2.3',
    'alambda      = 2.3',
    'I            = 13',
    'qmin         = 12000000',
    'qrange       = 1000000',
    'checkrange   = 1',
    'firstcheck   = 1',

    'delay        = 120',
    'sievenice    = 19',
    'keeprelfiles = 0',
    'sieve_max_threads = 2',
    # 'poly_max_threads = 1', # 1 is the default anyway
    # 'ratq	 = 0', # 0 is the default anyway

    # filtering
    # 'skip         = -1', # should be about bwc_mn - 32
    # 'keep         = -1', # should be 128 + skip
    'keeppurge    = 208', # should be 160 + #ideals <= FINAL_BOUND (cf purge.c)
    'maxlevel     = 15',
    'ratio        = 1.5',
    'bwstrat      = 3',
    'coverNmax    = 100',
    # 'nslices_log  = 1', # 1 is the default anyway
    'filterlastrels = 1',
    # 'dup_rm = $on_mingw', # Give -rm parameter to dup2 when on MinGW

    # linalg
    'linalg       = bwc',
    'bwmt         = 2',
    'mpi          = 0',
    'hosts	 = ""',
    'bwc_interval = 1000',
    'bwc_mm_impl = bucket',
    'bwc_interleaving = 0',
    # bwc_mn should be 64 or 128
    'bwc_mn       = 64',
    # shuffled product is expected to be better in most cases', at least
    # when we use MPI. Since it is the preferred communication algorithm
    # for large runs', we prefer to force its use also for mid-range
    # examples.
    'bwc_shuffled_product = 1',

    # characters
    'nkermax      = 30',
    'nchar        = 50',
    'nthchar      = 2',

    # holy grail
    #'expected_factorization = undef',

    # logfile
    #'logfile = undef',
)

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        import doctest
        doctest.testmod()
    elif len(sys.argv) == 2:
        argsold = True
        paramfile = sys.argv[1]
        parameters = Parameters()
        if argsold:
            parameters.read_old_defaults()
        parameters.readfile(paramfile, old_format = argsold)
        
        # Keep only the last occurence of each key. First reverse
        parse_array.reverse()
        # Now keep only the first occurence of each key
        keys_seen = set()
        filtered = []
        for (key, value, comment, oldkey) in parse_array:
            # We keep entries without an oldkey, as those correspond to comment
            # or blank lines in the input file, which we try to preserve.
            # This is also why we use two arrays instead of an OrderedDict.
            # Lines with an oldkey that maps to None are discarded.
            keep = oldkey is None or not (key is None or key in keys_seen)
            if False:
                print("%skeeping %s=%s #%s /%s" %
                      ("" if keep else "not ", key, value, comment, oldkey) )
            keys_seen.add(key)
            if keep:
                filtered.append((key, value, comment, oldkey))
        
        # Revert again to original order
        filtered.reverse()
        
        for (key ,value, comment, oldkey) in filtered:
            output = []
            if not key is None:
                output.append("%s = %s" % (key, value))
            if not comment is None:
                output.append("#%s" % comment)
            print("\t\t".join(output))
