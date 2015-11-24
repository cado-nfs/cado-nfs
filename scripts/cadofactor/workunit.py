class _ScalarKey(object):
    """ Keys that must occur at most once """
    (takes_value, takes_multiple, takes_checksum) = True, False, False
class _SignalKey(object):
    """ Keys that must occur at most once and which take no value """
    (takes_value, takes_multiple, takes_checksum) = False, False, False
class _ListKey(object):
    """ Keys that can occur zero or more times """
    (takes_value, takes_multiple, takes_checksum) = True, True, False
class _ChecksummedKey(object):    
    """ Keys that can be accompanied by a CHECKSUM line. These keys can occur 
    zero or more times, like those in LIST_KEYS
    """
    (takes_value, takes_multiple, takes_checksum) = True, True, True

class Workunit(object):
    # Keys and their type as a tuple in their preferred ordering
    KEYS = (("WORKUNIT", _ScalarKey), ("TERMINATE", _SignalKey), 
        ("FILE", _ChecksummedKey),
        ("EXECFILE", _ChecksummedKey),
        ("SUGGEST_EXECFILE", _ScalarKey),
        ("COMMAND", _ListKey), ("RESULT", _ListKey), ("DELETE", _ListKey), 
        ("CHECKSUM", _ListKey))
    # The type for CHECKSUM does not really matter so long as it's not 
    # _SignalKey, as we don't add the value of CHECKSUMs to the dict directly
    
    def __init__(self, text):
        """ Init a workunit from the text of a WU file """
        KEYS = dict(self.__class__.KEYS)
        wu = {}
        checksum_key = None

        for line in text.splitlines():
            # Drop leading/trailing whitespace, incl. CR/LF. Split first word 
            # and rest of line
            s = line.strip().split(" ", 1)
            key = s[0]

            if not key in KEYS:
                raise Exception("Error: key " + key + " not recognized")
            if len(s) == 2:
                # Drop leading whitespace from value
                value = s[1].lstrip()
                if not KEYS[key].takes_value:
                    raise Exception("Key " + key + " with value " + str(value))
            else:
                value = None
                if KEYS[key].takes_value:
                    raise Exception("Key " + key + " without value")

            # Handle CHECKSUM keys separately
            if key == "CHECKSUM":
                if not checksum_key:
                    raise Exception("Extraneous " + key)
                # Store this checksum along with its FILE or EXECFILE entry
                assert wu[checksum_key][-1][1] == None
                wu[checksum_key][-1][1] = value
                checksum_key = None
                continue
            # if this isn't a CHECKSUM line, forget what the last line was
            checksum_key = None

            if not KEYS[key].takes_multiple:
                if key in wu:
                    raise Exception("Key " + key + " redefined")
                wu[key] = value
            else:
                if KEYS[key].takes_checksum:
                    # Append the filename, plus None as the checksum
                    value = [value, None]
                    checksum_key = key
                if not key in wu:
                    wu[key] = []
                wu[key].append(value)
        # END: for line in text.splitlines()
        self.wudata = wu
        # The WORKUNIT key must be given
        if not "WORKUNIT" in wu:
            raise Exception("Workunit has no WORKUNIT line")
    
    def  __str__(self):
        """ Produce text for a WU, as could be stored in a WU file """
        result = ""
        for (key, keytype) in self.__class__.KEYS:
            if not key in self.wudata:
                continue
            if not keytype.takes_value:
                result = result + key + "\n"
            elif not keytype.takes_multiple:
                result = result + key + " " + self.wudata[key] + "\n"
            elif keytype.takes_checksum:
                for (name, checksum) in self.wudata[key]:
                    result = result + key + " " + name + "\n"
                    if not checksum == None:
                        result = result + "CHECKSUM " + checksum + "\n"
            else:
                for value in self.wudata[key]:
                    result = result + key + " " + value + "\n"
        return result
    
    def get_id(self):
        """ Get the Workunit ID """
        return self.wudata["WORKUNIT"]

    def set_id(self, wuid):
        """ Set the Workunit ID, overriding the old value """
        self.wudata["WORKUNIT"] = wuid

    def get(self, *args, **kwargs):
        """ Delegates to the wudata dictionary """
        return self.wudata.get(*args, **kwargs)

def wu_test():
    """ Dummy function to test workunit parser 
    
    >>> Workunit("")
    Traceback (most recent call last):
    Exception: Workunit has no WORKUNIT line
    
    >>> Workunit("FOO")
    Traceback (most recent call last):
    Exception: Error: key FOO not recognized
    
    >>> str(Workunit("WORKUNIT a\\n"))
    'WORKUNIT a\\n'
    
    >>> Workunit("WORKUNIT\\n")
    Traceback (most recent call last):
    Exception: Key WORKUNIT without value
    
    >>> Workunit("FILE\\n")
    Traceback (most recent call last):
    Exception: Key FILE without value
    
    >>> Workunit("EXECFILE\\n")
    Traceback (most recent call last):
    Exception: Key EXECFILE without value
    
    >>> Workunit("RESULT\\n")
    Traceback (most recent call last):
    Exception: Key RESULT without value
    
    >>> Workunit("WORKUNIT a\\nWORKUNIT b\\n")
    Traceback (most recent call last):
    Exception: Key WORKUNIT redefined

    >>> str(Workunit("WORKUNIT a\\nFILE foo\\n"))
    'WORKUNIT a\\nFILE foo\\n'
    
    >>> str(Workunit("WORKUNIT a\\nCHECKSUM 0x1\\n"))
    Traceback (most recent call last):
    Exception: Extraneous CHECKSUM
    
    >>> str(Workunit("WORKUNIT a\\nFILE foo\\nCHECKSUM 0x1\\n"))
    'WORKUNIT a\\nFILE foo\\nCHECKSUM 0x1\\n'
    
    >>> str(Workunit("WORKUNIT a\\nEXECFILE foo\\nCHECKSUM 0x1\\n"))
    'WORKUNIT a\\nEXECFILE foo\\nCHECKSUM 0x1\\n'
    
    >>> str(Workunit("WORKUNIT a\\nFILE foo\\nCHECKSUM 0x1\\nCHECKSUM 0x1\\n"))
    Traceback (most recent call last):
    Exception: Extraneous CHECKSUM
    
    >>> str(Workunit("WORKUNIT a\\nFILE foo\\nEXECFILE efoo\\nCHECKSUM 0x1\\nFILE bar\\nCOMMAND ${DLDIR}foo bar >> ${WORKDIR}baz\\nCOMMAND ${DLDIR}another >> ${WORKDIR}result\\nRESULT baz\\nRESULT result\\n"))
    'WORKUNIT a\\nFILE foo\\nFILE bar\\nEXECFILE efoo\\nCHECKSUM 0x1\\nCOMMAND ${DLDIR}foo bar >> ${WORKDIR}baz\\nCOMMAND ${DLDIR}another >> ${WORKDIR}result\\nRESULT baz\\nRESULT result\\n'

    >>> str(Workunit("WORKUNIT a\\nTERMINATE\\n"))
    'WORKUNIT a\\nTERMINATE\\n'

    >>> str(Workunit("WORKUNIT a\\nTERMINATE x\\n"))
    Traceback (most recent call last):
    Exception: Key TERMINATE with value x

    """
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
