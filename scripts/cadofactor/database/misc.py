
def join3(L, pre=None, post=None, sep=", "):
    """
    If any parameter is None, it is interpreted as the empty string
    >>> join3 ( ('a'), pre="+", post="-", sep=", ")
    '+a-'
    >>> join3 ( ('a', 'b'), pre="+", post="-", sep=", ")
    '+a-, +b-'
    >>> join3 ( ('a', 'b'))
    'a, b'
    >>> join3 ( ('a', 'b', 'c'), pre="+", post="-", sep=", ")
    '+a-, +b-, +c-'
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = ""
    return sep.join([pre + k + post for k in L])


def dict_join3(d, sep=None, op=None, pre=None, post=None):
    """
    If any parameter is None, it is interpreted as the empty string
    >>> t = dict_join3({"a": "1", "b": "2"},
    ...                sep=",", op="=", pre="-", post="+")
    >>> t == '-a=1+,-b=2+' or t == '-b=2+,-a=1+'
    True
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = ""
    if op is None:
        op = ""
    return sep.join([pre + op.join(k) + post for k in d.items()])
