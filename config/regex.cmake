#
# We do require POSIX, and regex.h is provided by posix. Hence we should
# not have to care. Now, regex.h is arguably a seldom used header. It
# seems that some installs manage to botch it, and that does not prevent
# release: the intel compiler suite, when python is also enabled,
# exposes a failing regex.h file that takes priority over all others.
CHECK_INCLUDE_FILES (regex.h HAVE_REGEX_H)
