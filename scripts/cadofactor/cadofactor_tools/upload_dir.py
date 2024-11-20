import os
import logging

logger = logging.getLogger(__name__)


class UploadDirProvider(object):
    """
    Generate name of desired upload subdirectory

    If nrsubdir == 0, then the base upload directory is used.
    If nrsubdir > 0 (incl. when nrsubdir == 1), subdirectories numbered
    uploaddir_base + "/0/", uploaddir_base + "/1/" etc are used, where
    variable i determines the subdirectory.
    """
    def __init__(self, uploaddir_base, nrsubdir):
        self.uploaddir_base = uploaddir_base
        self.nrsubdir = nrsubdir

        if not os.path.isdir(self.uploaddir_base):
            logger.debug("Creating upload directory %s" % self.uploaddir_base)
            os.mkdir(self.uploaddir_base)

        for i in range(nrsubdir):
            uploaddir_i = self(i)
            if not os.path.isdir(uploaddir_i):
                logger.debug("Creating upload subdirectory %s" % uploaddir_i)
                os.mkdir(uploaddir_i)

    def __call__(self, i):
        if self.nrsubdir == 0:
            return self.uploaddir_base
        else:
            # Empty final segment to make the path end in a directory separator
            return os.path.join(self.uploaddir_base,
                                "%d" % (i % self.nrsubdir))
