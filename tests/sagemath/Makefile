
check:
	if ! [ -d /tmp/venv ] ; then python3 -m venv /tmp/venv ; fi
	if ! [ -x /tmp/venv/bin/flake8 ] ; then /tmp/venv/bin/pip3 install flake8 ; fi
	/tmp/venv/bin/flake8 --ignore=E201,E202,W503,W504,E731,E741 cado_sage

# E201 / E202 enforce the absence of whitespace in lists. It's quite
# often a minus in terms of readability in my opinion
#
# W503 / W504 just look up
# https://stackoverflow.com/questions/57074300/what-is-the-recommended-way-to-break-a-long-if-statement-w504-line-break-after
# for an absurd discussion on what the best practice is. I don't care.
#
# E731: lambdas are not evil
#
# E741: I definitely want to use I.


# For the moment most files do not have doctests.
doctests:
	PYTHONPATH=. /usr/local/bin/sage -t cado_sage
