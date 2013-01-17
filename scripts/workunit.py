class Workunit:
    # Keys that must occur only once
    SCALAR_KEYS = ("WORKUNIT",)
    # Keys that can occur multiple times
    LIST_KEYS = ("COMMAND", "RESULT")
    # Keys that can be accompanied by a CHECKSUM line. These keys can occur 
    # multiple times, like those in LIST_KEYS
    CHECKSUM_KEYS = ("FILE", "EXECFILE")
    
    def __init__(self, text):
        """ Init a workunit from the text of a WU file """
        wu = {}
        for key in self.LIST_KEYS + self.CHECKSUM_KEYS:
            wu[key] = [] # Init to empty list
        checksum_key = None
        for line in text.splitlines():
            # Split first word and rest of line
            (key, value) = line.split(" ", 1)
            # Drop leading/trailing whitespace, incl. CR/LF
            value = value.lstrip().rstrip()
            # Handle CHECKSUM keys separately
            if key == "CHECKSUM":
                if not checksum_key:
                    raise Exception("Extraneous " + key)
                # Store this checksum along with its FILE or EXECFILE entry
                wu[checksum_key][-1][1] = value
                checksum_key = None
                continue
            # if this isn't a CHECKSUM line, forget what the last line was
            checksum_key = None
            if key in self.SCALAR_KEYS:
                if key in wu:
                    raise Exception("Key " + key + " redefined")
                wu[key] = value
            elif key in self.LIST_KEYS:
                wu[key].append(value)
            elif key in self.CHECKSUM_KEYS:
                # Append the filename, plus None as the checksum
                wu[key].append([value, None])
                checksum_key = key
            else:
                raise Exception("Error: key " + key + " not recognized")
        self.data = wu
    
    def  __str__(self):
        """ Produce text for a WU, as could be stored in a WU file """
        str = ""
        for key in self.SCALAR_KEYS:
            if key in self.data:
                str = str + key + " " + self.data[key] + "\n"
        for key in self.CHECKSUM_KEYS:
            if key in self.data:
                for (name, checksum) in self.data[key]:
                    str = str + key + " " + name + "\n"
                    if not checksum == None:
                        str = str + "CHECKSUM " + checksum + "\n"
        for key in self.LIST_KEYS:
            if key in self.data:
                for value in self.data[key]:
                    str = str + key + " " + value + "\n"
        return str
    
    def get_id(self):
        return self.data["WORKUNIT"]
    
    