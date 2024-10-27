import re


# A memory image of the reconstructlog.dlog file.
class LogBase(object):
    def __init__(self, general):
        self.general = general
        self.known = {}
        self.badideals = []
        self.SMs = [ [], [] ]
        self.fullcolumn = None
        try:
            print("--- Reading %s to find which are the known logs ---"
                  % general.log())

            def process(line):
                index, p, side, r, *value = line.split()
                if p == "bad" and side == "ideals":
                    self.badideals.append(int(r))
                elif p == "SM":
                    self.SMs[int(side)].append(int(value[0]))
                else:
                    # for rational side, we actually don't have the root.
                    # We force it to be on side 0 (legacy, again...)
                    if r == "rat":
                        assert int(side) == 0
                        r = -1
                    else:
                        r = int(r, 16)
                    self.known[(int(p, 16), r, int(side))] = int(value[0])

            with open(general.log(), 'r') as file:
                line = file.readline()
                m = re.match(r"^(\w+) added column (\d+)$", line)
                if m:
                    self.fullcolumn = int(m.groups()[1])
                else:
                    self.fullcolumn = None
                    process(line)
                for i, line in enumerate(file):
                    if i % 1000000 == 0:
                        print("Reading line %d" % i)
                    process(line)
            print("Found",
                  f"{len(self.badideals)} bad ideals,",
                  f"{len(self.known)} known logs,",
                  f"and {len(self.SMs[0])},{len(self.SMs[1])} SMs",
                  "in general.log()")
        except Exception:
            raise ValueError("Error while reading %s" % general.log())

    def has(self, p, r, side):
        return (p, r, side) in self.known

    def get_log(self, p, r, side):
        if self.general.has_rational_side() and side == 0:
            r = -1
        if (p, r, side) in self.known:
            return self.known[(p, r, side)]
        else:
            return None

    def add_log(self, p, r, side, log):
        if self.general.has_rational_side() and side == 0:
            r = -1
        self.known[(p, r, side)] = log

    def bad_ideal(self, i):
        return self.badideals[i]

    def allSM(self, i):
        SM = self.SMs[0] + self.SMs[1]
        return SM[i]

    def nSM(self, side):
        return len(self.SMs[side])

    def SM(self, side, i):
        return self.SMs[side][i]

    def full_column(self):
        return self.fullcolumn
