import logging

from logging import DEBUG, INFO, WARNING, ERROR, FATAL, CRITICAL
MESSAGE = DEBUG - 1
TRANSACTION = MESSAGE - 1
COMMAND = DEBUG + 1
# TRANSACTION = INFO
logging.addLevelName(MESSAGE, "Message")
logging.addLevelName(TRANSACTION, "Transaction")
logging.addLevelName(COMMAND, "Command")


class ANSI(object):
    """ Class defining some ANSI control sequences, for example for 
    changing text colour """
    CSI = '\x1b[' # ANSI Control Sequence Introducer. Not the TV show
    SGR = 'm' # Set Graphics Rendition code
    NORMAL = CSI + '0' + SGR
    BLACK = CSI + '30' + SGR 
    GREY = CSI + '30;1'
    RED = CSI + '31' + SGR 
    BRIGHTRED = CSI + '31;1' + SGR
    GREEN = CSI + '32' + SGR
    BRIGHTGREEN = CSI + '32;1' + SGR
    YELLOW = CSI + '33' + SGR
    BRIGHTYELLOW = CSI + '33;1' + SGR
    BLUE = CSI + '34' + SGR
    BRIGHTBLUE = CSI + '34;1' + SGR
    MAGENTA = CSI + '35' + SGR
    BRIGHTMAGENTA = CSI + '35;1' + SGR
    CYAN = CSI + '36' + SGR
    BRIGHTCYAN = CSI + '36;1' + SGR
    # There is another grey with code "37" - white without intensity
    # Not sure if it is any different from "30;1" aka "bright black"
    WHITE = CSI + '37;1' + SGR

class ScreenFormatter(logging.Formatter):
    """ Class for formatting logger records for screen output, optionally
    with colorized logger level name (like cadofct.pl used to). """
    colours = {
        MESSAGE: ANSI.GREEN,
        TRANSACTION: ANSI.CYAN,
        INFO : ANSI.BRIGHTGREEN,
        WARNING : ANSI.BRIGHTYELLOW,
        ERROR : ANSI.BRIGHTRED,
        FATAL : ANSI.BRIGHTRED,
        COMMAND : ANSI.BRIGHTBLUE
    }
    
    # Format string that switches to a different colour (with ANSI code 
    # specified in the 'colour' key of the log record) for the log level name, 
    # then back to default text rendition (ANSI code in 'nocolour')
    colourformatstr = \
        '%(padding)s%(colour)s%(levelnametitle)s%(nocolour)s:%(name)s: %(message)s'
    # Format string that does not use colour changes
    nocolourformatstr = \
        '%(padding)s%(levelnametitle)s:%(name)s: %(message)s'
    
    def __init__(self, colour = True):
        self.useColour = colour
        if colour:
            super().__init__(fmt=self.colourformatstr)
        else:
            super().__init__(fmt=self.nocolourformatstr)
    
    def format(self, record):
        # Add attributes to record that our format string expects
        if self.useColour:
            record.colour = self.colours.get(record.levelno, ANSI.NORMAL)
            record.nocolour = ANSI.NORMAL
        if record.levelno == COMMAND:
            record.levelnametitle = 'Running Command'
        else:
            record.levelnametitle = record.levelname.title()
        if hasattr(record, "indent"):
            record.padding = " " * record.indent
        else:
            record.padding = ""
        return super().format(record)


class FileFormatter(logging.Formatter):
    """ Class for formatting a log record for writing to a log file. No colours 
    here, but we add the process ID and a time stamp """
    formatstr = \
       'PID%(process)s %(asctime)s %(levelnametitle)s:%(name)s: %(message)s' 
    
    def format(self, record):
        record.levelnametitle = record.levelname.title()
        return super().format(record)
    
    def __init__(self):
        super().__init__(fmt=self.formatstr)


class CmdFileFormatter(logging.Formatter):
    """ Class for formatting a log record for writing to a command log file. 
        No colours here, but we add the process ID, the spanwed process id, 
        and a time stamp """
    formatstr = \
       '# PPID%(process)s PID%(childpid)s %(asctime)s\n%(message)s' 
    
    def __init__(self):
        super().__init__(fmt=self.formatstr)


class ScreenHandler(logging.StreamHandler):
    def __init__(self, lvl = logging.INFO, colour = True, **kwargs):
        super().__init__(**kwargs)
        self.setLevel(lvl)
        self.setFormatter(ScreenFormatter(colour = colour))


class FileHandler(logging.FileHandler):
    def __init__(self, filename, lvl = logging.DEBUG, **kwargs):
        super().__init__(filename, **kwargs)
        self.setLevel(lvl)
        self.setFormatter(FileFormatter())


class CmdFileFilter(logging.Filter):
    def filter(self, record):
        return record.levelno == COMMAND


class CmdFileHandler(logging.FileHandler):
    def __init__(self, filename , **kwargs):
        super().__init__(filename, **kwargs)
        cmdfilter = CmdFileFilter()
        self.addFilter(cmdfilter)
        self.setFormatter(CmdFileFormatter())


class MyLogger(logging.Logger):
    def cmd(self, msg, pid, *args, extra=None, **kwargs):
        """ Log a message with a level of cadologger.COMMAND """
        if not extra is None:
            extra = extra.copy()
        else:
            extra = {}
        extra["childpid"] = pid
        
        self.log(COMMAND, msg, *args, extra=extra, **kwargs)
    
    def message(self, msg, *args, **kwargs):
        """ Log a message with a level of cadologger.MESSAGE """
        self.log(MESSAGE, msg, *args, **kwargs)
    
    def transaction(self, msg, *args, **kwargs):
        """ Log a message with a level of cadologger.TRANSACTION """
        self.log(TRANSACTION, msg, *args, **kwargs)
    

logging.setLoggerClass(MyLogger)
root = logging.getLogger()
root.setLevel(logging.NOTSET)
del(root)

def init_test_logger():
    logger = logging.getLogger()
    logger.addHandler(ScreenHandler(lvl = logging.NOTSET))
    logger.addHandler(FileHandler(filename = "log", lvl = logging.NOTSET))

if __name__ == '__main__':
    init_test_logger()

    logger = logging.getLogger()
    logger.info("An Info Center!")
    logger.warning("Beware")
    logger.error("All hope abandon", extra={"indent" : 4})
