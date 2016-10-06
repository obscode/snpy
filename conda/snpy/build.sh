# the main event
LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"
$PYTHON setup.py install
