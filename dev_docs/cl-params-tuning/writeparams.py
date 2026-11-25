import json, sys
from mkparams import emit
C = json.load(open("choices.json"))
def fmt(x):
    if x >= 365 * 3: return f"{x/365:.0f} core-years"
    if x >= 100: return f"{x:.0f} core-days"
    if x >= 3: return f"{x:.0f} core-days"
    if x >= 0.2: return f"{x*24:.0f} core-hours"
    return f"{x*24:.1f} core-hours"
STATUS = {80: "Checked with a complete computation.", 90: "Checked with a complete computation.",
          100: "Sieving and filtering statistics measured, not run to completion.",
          110: "Sieving statistics measured; filtering simulated; not run.",
          120: "Sieving statistics measured; filtering simulated; not run.",
          130: "Extrapolated (sieving speed measured, yield and filtering modelled); not run.",
          140: "Extrapolated (sieving speed measured, yield and filtering modelled); not run.",
          150: "Extrapolated (sieving speed measured, yield and filtering modelled); not run.",
          160: "Extrapolated (sieving speed measured, yield and filtering modelled); not run."}
QR = {80: 100, 90: 300, 100: 250, 110: 200, 120: 70, 130: 60, 140: 50, 150: 30, 160: 20}
for c in C:
    c = dict(c)
    if c['d'] == 90 and len(sys.argv) > 2:
        continue
    c["margin"] = 1.0
    c["qrange"] = QR[c["d"]]
    c["notes"] = (
        "Parameters chosen (November 2026) to minimise sieving + linear algebra\n"
        "together: unlike for factoring, the linear algebra cost grows like the\n"
        "cube of the matrix size (the determinants have about 1.9 bits per row, so\n"
        "the number of moduli grows with the matrix).\n"
        f"{STATUS[c['d']]}\n"
        "\n"
        "Expected for a typical discriminant, in core time of a 2.1 GHz Xeon\n"
        f"(Emerald Rapids): about {fmt(c['total'])} (sieving {fmt(c['sieve'])},\n"
        f"linear algebra {fmt(c['la'])}), {c['nsq']:.2g} special-q's, a matrix with about\n"
        f"{c['N']} rows, {c['mem']} GB per sieving process.\n"
        "The sieving time depends on which small primes split: 30 % more or less\n"
        "is common between discriminants of the same size.\n")
    fn, info = emit(c, sys.argv[1])
    print(fn, info)
