#!/usr/bin/env python3

import json
import sys
import re
import os

if __name__ == '__main__':
    for filename in sys.argv[1:]:
        sz = os.stat(filename).st_size
        print(f">>> Processing {filename} ({sz} bytes)", file=sys.stderr)
        with open(filename, "r") as f:
            j = json.loads(f.read())

        all_files = [ f["file"] for f in j['files'] ]

        remove_directories = [
                r"utils/embedded",
                r"gf2x/",
                r"linalg/bwc/flint-fft",
                r"linalg/bwc/mpfq"
                ]
        matcher = r"^(" + r"|".join(remove_directories) + r")"
        removed_files = [ f for f in all_files if re.match(matcher, f) ]

        jfiles = { x["file"]:x for x in j["files"] }

        for f in removed_files:
            print(f">>> removing {f}", file=sys.stderr)
            del jfiles[f]

        j["files"] = list(jfiles.values())

        with open(filename, "w") as f:
            print(json.dumps(j), file=f)

        sz = os.stat(filename).st_size
        print(f">>> Written {filename} ({sz} bytes)", file=sys.stderr)
