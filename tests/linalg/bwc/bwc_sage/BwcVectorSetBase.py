import os
import re


class BwcVectorSetBase(object):
    def __init__(self, objectbase, params, dirname):
        self.__objectbase = objectbase
        self.__what = objectbase._what
        self.__whats = objectbase._whats
        self.__pattern = objectbase.pattern

        self.params = params
        self.Vs = []
        self.dirname = dirname
        self.vfiles = []

        doing = f"Scanning for {self.__whats} in {dirname}"
        print(doing)
        for x in os.listdir(dirname):
            if re.match(self.__pattern, x):
                fn = os.path.join(dirname, x)
                self.Vs.append(self.__objectbase(params, fn))

        self.Vs = sorted(self.Vs, key=lambda x: tuple(x))
        print(f"{doing} ... {len(self.Vs)} {self.__whats} found")

    def read(self):
        for v in self.Vs:
            v.read()

    def __iter__(self):
        return iter(self.Vs)

    def __getitem__(self, i):
        return self.Vs[i]

    def __len__(self):
        return len(self.Vs)
