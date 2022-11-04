

# wait4 is nonstandard, and bsd-specific. Not all platforms have it,
# unfortunately, and we _are_ using it in cado_popen. So let's do some
# work to detect it.
search_for_function(wait4 HAVE_WAIT4)
