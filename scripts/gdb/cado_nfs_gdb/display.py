display_limit = 40


def truncate_output(r):
    # not completely satisfactory
    if display_limit == 0:
        return r
    ell = len(r)
    if ell >= display_limit:
        r = r[0:display_limit-15] \
            + "...<%d>..." % (ell - (display_limit - 12)) + r[ell - 7:ell]
    return r


def update_display_limit(ell):
    global display_limit
    display_limit = ell
