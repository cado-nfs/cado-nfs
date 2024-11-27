import json
import hashlib

r"""
Workunit:
    type: object
    properties:
        id:
            type: string
            description: name of the workunit, used internally in the
                         server database
        commands:
            type: array
            items:
                type: string
                description: command that the client is expected to run.
                             Substitutions are performed based on the files.
        timeout:
            type: integer
            description: how long before the server will assume that the
                         client has failed and is gone, from which point
                         on this workunit will be reassigned.
        files:
            type: object
            properties:
                filename:
                    type: string
                checksum:
                    type: string
                algorithm:
                    type: string
                    enum:
                        - sha1
                        - sha256
                        - sha3_256
                    description: checksumming algorithm that was used to
                                 compute the checksum (as in hashlib)
                upload:
                    type: boolean
                    description: whether this should be uploaded by the
                                 client.  Exclusive with checksum
                download:
                    type: boolean
                    description: whether this should be downloaded by the
                                 client.  Normally comes with checksum
                suggest_path:
                    type: string
                    description: suggestion of the subpath where this file
                                 may be found in the cado build tree. Used
                                 a priori for executables.
            required:
              - filename
"""


class Workunit(dict):
    def __init__(self, *args, **kwargs):
        if args and type(args[0]) is str:
            assert len(args) == 1
            assert not kwargs
            super().__init__(json.loads(args[0]))
        else:
            super().__init__(*args, **kwargs)
        self._check_schema()

    def _check_schema(self):
        # check schema
        if 'id' not in self:
            raise Exception("Workunit has no id")

        for key, value in self.items():
            if key not in ['id', 'files', 'commands', 'timeout']:
                raise Exception(f"key {key} not recognized")
            if key == 'timeout':
                if type(value) is int:
                    pass
                else:
                    try:
                        tt = float(value)  # noqa: F841
                    except ValueError:
                        raise Exception("timeout must be numerical")

        for key, value in self.get('files', {}).items():
            for field, info in value.items():
                if field not in ['filename',
                                 'checksum',
                                 'upload',
                                 'download',
                                 'algorithm',
                                 'suggest_path']:
                    raise Exception(f"field {field} not recognized")
                if field == 'algorithm':
                    if getattr(hashlib, info, None) is None:
                        raise Exception(f"hash algorithm {info} unknown")
            if value.get('upload', False) and 'checksum' in value:
                raise Exception("fields 'upload' and 'checksum' are exclusive")
            if 'checksum' in value and 'algorithm' not in value:
                raise Exception("field 'checksum' requires 'algorithm'")
            if 'algorithm' in value and 'checksum' not in value:
                raise Exception("field 'algorithm' requires 'checksum'")
            if 'filename' not in value:
                raise Exception("field 'filename' is required")
            if 'suggest_path' in value and value.get('output', False):
                raise Exception("field 'suggest_path'"
                                " is meaningless for output files")

    def __str__(self):
        return json.dumps(self)

    def get_id(self):
        return self['id']

    def set_id(self, wuid):
        self['id'] = wuid


def wu_test():
    """ Dummy function to test workunit parser

    >>> Workunit()
    Traceback (most recent call last):
    Exception: Workunit has no id

    >>> Workunit(id='a1', foo='a')
    Traceback (most recent call last):
    Exception: key foo not recognized
    """
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
